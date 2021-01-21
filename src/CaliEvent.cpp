#include <iostream>
#include <fstream>
#include "CaliEvent.h"
#include "detconf.h"

/// Declaration of calibration parameters
static double cell_bigd_sma[900];
static double cell_sm_ad[900];

static double bigpmt_spe[900];
static double bigpmt_da[900];
static double bigpmt_tcalib[900];

void CaliEvent::GetWCDA_QT_Cali()
{
  TFile *file = new TFile("./califile/qcal_r20201110.root");
  printf("    -->input wcda adc file!!\n","");
  TTree *tree = (TTree *)file->Get("tfit");
  int igcell;
  int imode;
  double ratio;
  tree->SetBranchAddress("igcell",&igcell);
  tree->SetBranchAddress("imode",&imode);
  tree->SetBranchAddress("ratio",&ratio);

  int n = tree->GetEntries();
  for(int i=0;i<n;i++)
  {
    tree->GetEntry(i);
    if(igcell>=900) continue;
    if(imode==0)  bigpmt_da[igcell] = 1./ratio;
    if(imode==1)  cell_bigd_sma[igcell] = 1./ratio;
    if(imode==2)  cell_sm_ad[igcell] = 1./ratio;
    if(imode==-1) bigpmt_spe[igcell] = ratio;
  }
  file->Close();
  /*
  for(int i=0;i<900;i++) {
    printf("%d %lf %lf %lf %lf \n",i,bigpmt_da[i] ,cell_bigd_sma[i],cell_sm_ad[i],bigpmt_spe[i]);
  }
  */

  ifstream tfile;
  tfile.open("./califile/tcal_P1LED_2020121000_1.txt");
  printf("    -->input wcda tdc file!!\n","");
  char buf[500];
  double tmp6;
  while (tfile.getline(buf,200)) {
    sscanf(buf,"%d %lf\n",&igcell, &tmp6);
    bigpmt_tcalib[igcell]=tmp6;
  }
  tfile.close();
}

void CaliEvent::WCDA_Cali( SaveEvent *wcdaevent, LHCALEvent *lhcalevent)
{
  double peBA_same_igFlag[4000]={0.};
  double peBD_same_igFlag[4000]={0.};

  TClonesArray *Hitlist;
  TClonesArray *Hitlist_plus;
  SaveWCDAHit *Hit;
  SavePLUSHit *HitPlus;

  Hitlist = wcdaevent->hits;
  Hitlist_plus = wcdaevent->hitsplus;

  int nHits = wcdaevent->nhit;
  int nHits_plus = wcdaevent->nhitplus;

  int ifee, ich, iconf;
  int igcell;
  double big_charge_anode;
  double big_charge_dynode;
  for(int j=0; j<nHits; j++)
  {
    Hit = (SaveWCDAHit *)((*Hitlist)[j]);

    //////// get big igcell /////////
    ifee = Hit->Fee();
    ich = Hit->Ch();

    igcell = findig( ifee, ich );
    if(igcell<0||igcell>=900) { continue; }  // only 1# pool calibrated.

    /////////////// Charge calibration  ///////////////////
    big_charge_anode = Hit->AnodeCharge();
    big_charge_dynode = Hit->DynodeCharge();
    peBA_same_igFlag[igcell] = big_charge_anode / bigpmt_spe[igcell];
    peBD_same_igFlag[igcell] = big_charge_dynode * bigpmt_da[igcell] / bigpmt_spe[igcell];

    double q0 = peBA_same_igFlag[igcell];
    double q1 = peBD_same_igFlag[igcell];
    double q = q0;
    int s_igcell = igcell; //

    double q_sm=-10;
    if( big_charge_anode>3500 && q1>0. ) { q=q1;}

    if( big_charge_dynode>=2900)
    {
      //cout<<"###########big=="<<q1;
      int fee=0,db=0,pmt=0;
      //////// get small igcell /////////
      for(int i=0; i<nHits_plus;i++)
      {
        HitPlus = (SavePLUSHit *)((*Hitlist_plus)[i]);

        fee=HitPlus->Fee(), db=HitPlus->Db(), pmt=HitPlus->Pmt();
        if(fee<0||db<0||pmt<0) continue;

	s_igcell = findplusig( fee, db, pmt);
        if(s_igcell<0||s_igcell>=900) { continue; }

        if(s_igcell==igcell)
        {
          if((HitPlus->AnodePeak()-HitPlus->AnodePed())>0&&(HitPlus->AnodePeak()-HitPlus->AnodePed())<3000){
            q_sm =  bigpmt_da[i]*cell_bigd_sma[igcell]*(HitPlus->AnodePeak()-HitPlus->AnodePed())/bigpmt_spe[igcell];
          }
          else{
            q_sm =  bigpmt_da[i]*cell_bigd_sma[igcell]*cell_sm_ad[igcell]* (HitPlus->DynodePeak()-HitPlus->DynodePed())/ bigpmt_spe[igcell];
          }
          break;
        }
      }
      //cout<<"###########sm=="<<q_sm<<endl;
      if(q_sm>0) q=q_sm;
    }
    if (q<=0){ q = 1.e-20;}

    /////////////// Time cali for big pmt ///////////////////
    double tp = Hit->DeltaTime() / 1000.0 -bigpmt_tcalib[igcell];
    //The curve correction is completed in the reconstruction program (ZZK)

    lhcalevent->cellig.push_back(igcell);
    lhcalevent->cellpe.push_back(q);
    lhcalevent->cellt.push_back(tp+2000);
  }

  lhcalevent->iwcdaevt = (long long)(wcdaevent->sec+wcdaevent->subsec)*1.e9;
  lhcalevent->nhitplus = wcdaevent->nhitplus;
  lhcalevent->wcda_trigerflag = wcdaevent->flag;
  lhcalevent->wcda_phi = wcdaevent->phi;
  lhcalevent->wcda_theta = wcdaevent->theta;

}
