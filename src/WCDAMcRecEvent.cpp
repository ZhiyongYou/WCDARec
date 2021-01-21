#include "WCDAMcRecEvent.h"
#include "WCDAMcEvent.h"
//double WFCTAX = -180.;// m
//double WFCTAY = 86;//m
double x_wfcta_inwcda=-161.801; //m
double y_wfcta_inwcda=-121.487; //m
double WCDAANGLE=(TMath::Pi()/2-29.442*TMath::DegToRad());
double WCDAANGLE_ = 29.442;
double sinwcdaangle = sin(WCDAANGLE);
double coswcdaangle = cos(WCDAANGLE);
ClassImp(WCDAMcRecEvent);

WCDAMcRecEvent::WCDAMcRecEvent()
{
}

WCDAMcRecEvent::~WCDAMcRecEvent()  {
}
void WCDAMcRecEvent::SetTelCenter(double wfctax,double wfctay)
{
   WFCTAX = wfctax;
   WFCTAY = wfctay;
}

void WCDAMcRecEvent::WCDA2Corsika(double azimuth,double corex, double corey,double *corsikaazi,double *corsikax,double *corsikay)
{
   double WCDA1CenterX = -76; //m 
   double WCDA1CenterY = -56; //m
   corex += WCDA1CenterX;
   corey += WCDA1CenterY; 
   corex = corex - x_wfcta_inwcda ;
   corey = corey - y_wfcta_inwcda ;
   double x = corex*coswcdaangle+corey*sinwcdaangle;
   double y = -corex*sinwcdaangle+corey*coswcdaangle;
   *corsikax = x+WFCTAX;
   *corsikay = y+WFCTAY;
//printf("%f %f\n",WFCTAX,WFCTAY);
   double azi = (azimuth+WCDAANGLE_ - 360*round((azimuth+WCDAANGLE_)/360)) - 90;
   if(azi<0) azi += 360;
   if(azi>360) azi -= 360; 
   *corsikaazi = azi;

}
void WCDAMcRecEvent::Corsika2WCDA(double azimuth,double corex, double corey,double *wcdaazi,double *wcdax,double *wcday)
{
	double WCDA1CenterX = -76; //m 
	double WCDA1CenterY = -56; //m

	double x_coretowfcta_incorsika = corex - WFCTAX;
	double y_coretowfcta_incorsika = corey - WFCTAY;

	double x_coretowfcta_inwcda = x_coretowfcta_incorsika*coswcdaangle - y_coretowfcta_incorsika*sinwcdaangle;
	double y_coretowfcta_inwcda = x_coretowfcta_incorsika*sinwcdaangle + y_coretowfcta_incorsika*coswcdaangle;
	double x_core_inwcda = x_coretowfcta_inwcda + x_wfcta_inwcda;
	double y_core_inwcda = y_coretowfcta_inwcda + y_wfcta_inwcda;

	*wcdaazi = -1000;
	*wcdax = x_core_inwcda - WCDA1CenterX;
	*wcday = y_core_inwcda - WCDA1CenterY;
}

void WCDAMcRecEvent::GetRecResult(double evrec_[70])
{
  double corsikaazi,corsikax,corsikay;
  for(int i=0; i<70; i++){
      evrec[i] = evrec_[i];
      
  }
  //printf("%f %f\n",evrec_[4],evrec_[5]);
  //printf("%f %f\n",evrec[4],evrec[5]);
  WCDA2Corsika(evrec[57],evrec[10], evrec[11],&corsikaazi,&corsikax,&corsikay);
  //Corsika2WCDA(evrec[57],evrec[10], evrec[11],&corsikaazi,&corsikax,&corsikay);
  azi_zzk_incorsika = corsikaazi;
  recx_zzk_incorsika = corsikax;
  recy_zzk_incorsika = corsikay;
}
void WCDAMcRecEvent::SetEvent(WCDAMcEvent *wcdaevent)
{
	Corsika2WCDA(wcdaevent->azic, wcdaevent->mcxc, wcdaevent->mcyc, &azic_inwcda, &mcxc_inwcda, &mcyc_inwcda);
	Corsika2WCDA(wcdaevent->azic, wcdaevent->mthx, wcdaevent->mthy, &azic_inwcda, &mthx_inwcda, &mthy_inwcda);
  fileid = wcdaevent->fileid;
  runid = wcdaevent->runid;
  eventid = wcdaevent->eventid;
  coreset= wcdaevent->coreset;
  e0= wcdaevent->e0;
  mczen= wcdaevent->mczen;
  mcazi= wcdaevent->mcazi;
  mcxc= wcdaevent->mcxc;
  mcyc= wcdaevent->mcyc;
  mthx= wcdaevent->mthx;
  mthy= wcdaevent->mthy;
  fitflag= wcdaevent->fitflag;
  recflag= wcdaevent->recflag;
  nhit= wcdaevent->nhit;
  npmtfire= wcdaevent->npmtfire;
  nfitc= wcdaevent->nfitc;
  ndetc= wcdaevent->ndetc;
  zenc= wcdaevent->zenc;
  azic= wcdaevent->azic;
  recxc= wcdaevent->recxc;
  recyc= wcdaevent->recyc;
  recxc1= wcdaevent->recxc1;
  recyc1= wcdaevent->recyc1;
  npea= wcdaevent->npea;
  npec= wcdaevent->npec;
  dcore= wcdaevent->dcore;
  dcore1= wcdaevent->dcore1;
  omega= wcdaevent->omega;
  compactness= wcdaevent->compactness;
  npmtshower= wcdaevent->npmtshower;
  nnpmtnoise= wcdaevent->nnpmtnoise;
  mynfitc= wcdaevent->mynfitc;
  nbig= wcdaevent->nbig;


}
void WCDAMcRecEvent::SetEvent(LHCALEvent *lhacalevent)
{

}

