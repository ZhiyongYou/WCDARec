/*  2016-04-18: reconstruction code for LHAASO fastMC
 *  2017-11-08: modified for G4KM2A geant4 simulation
 *  If you find any bug please send email to
 *  chensz@ihep.ac.cn
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "WCDAMcRecEvent.h"
#include "LHCALEvent.h"
#include "SaveEvent.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "WCDARec.h"
#include "CaliEvent.h"

const double DetectorAreaMD=36; //m^2
const double DetectorAreaED=1; //m^2
double WFCTAX = -180.;// m
double WFCTAY = 86;//m

void Coordinate_Corsika2WCDA(std::vector<double>& cellx, std::vector<double>& celly);

int main(int argc, char *argv[])
{
	if(argc<3)
	{
		printf("Usage: %s outfile in_file ... \n", argv[0]);
		return 1;
	}

	// WCDA Events  
	SaveEvent *wcdaevent = new SaveEvent();
	LHCALEvent *lhcalevent = new LHCALEvent();
	//output file and tree  rec event 
	WCDAMcRecEvent *wcdarecevent = new WCDAMcRecEvent();

	char Name1[300]="root://eos01.ihep.ac.cn/";
	char outfilename[500];
	strcpy(outfilename,Name1);
	strcat(outfilename,argv[1]);

	TFile *newfile = TFile::Open(outfilename, "recreate");
	newfile->SetCompressionSettings(101);

	WCDARec* wcdareconstruct = new WCDARec();
	//wcdareconstruct->SetWaterEff(argv[5]);
	wcdarecevent->SetTelCenter(WFCTAX,WFCTAY); 

	Long64_t sec;           //second of this event
	Double_t subsec;        //Subsecond of this event
	vector<int>* wcda_ig = new vector<int>();
	vector<double>* wcda_t = new vector<double>();
	vector<double>* wcda_pe = new vector<double>();
	vector<double>* wcda_x = new vector<double>();
	vector<double>* wcda_y = new vector<double>();
	TTree *EventRec = new TTree("Rec", "LHAASO Reconstructed Data");
	EventRec->Branch("WCDARecEvent", &wcdarecevent);
	EventRec->Branch("sec", &sec, "sec/L");
	EventRec->Branch("subsec", &subsec, "subsec/D");
	EventRec->Branch("nhitplus", &lhcalevent->nhitplus);
	//EventRec->Branch("iwcdaevt",&lhcalevent->iwcdaevt);
	EventRec->Branch("wcda_trigerflag",&lhcalevent->wcda_trigerflag);
	EventRec->Branch("wcda_phi",&lhcalevent->wcda_phi);
	EventRec->Branch("wcda_theta",&lhcalevent->wcda_theta);

//	EventRec->Branch("wcda_ig","vector<int>",&wcda_ig);
//	EventRec->Branch("wcda_t","vector<double>",&wcda_t);
//	EventRec->Branch("wcda_pe","vector<double>",&wcda_pe);
//	EventRec->Branch("wcda_x","vector<double>",&wcda_x);
//	EventRec->Branch("wcda_y","vector<double>",&wcda_y);
	EventRec->SetAutoSave(1000000);

	double evrec[70];  // reconstruction results by ZengZk method 
	//import the files and reconstructed them
	for(int nFile=2;nFile<argc;nFile++)
	{
		char Name2[300];
		strcpy(Name2,Name1);
		strcat(Name2,argv[nFile]);
		TFile *hfile= TFile::Open(Name2);
		if(!hfile) {   printf("%s does not exist\n",argv[nFile]);   continue;   }
		if(hfile->IsZombie()||hfile->GetEND()<50) {   printf("%s file error!!\n",argv[nFile]); hfile->Close();    continue;   }

		TTree *EventTree = (TTree *)hfile->Get("wcda");
		if(EventTree==nullptr) {  printf("%s is null file\n",argv[nFile]); hfile->Close();    continue;   }
		std::cout << argv[nFile] << " read succeed" << std::endl;
		EventTree->SetBranchAddress("wcda", &wcdaevent);
		int nentries = Int_t(EventTree->GetEntriesFast());
		printf("Total event %d\n",nentries);

		CaliEvent *calievent = new CaliEvent();
		calievent->GetWCDA_QT_Cali();

		for(int i=0; i<nentries; i++)
		{
			lhcalevent->Init();
			EventTree->GetEntry(i);
			calievent->WCDA_Cali( wcdaevent, lhcalevent );
			sec = wcdaevent->sec;
			subsec = wcdaevent->subsec;

			double maxpe = -1000;
			for(int ii=0;ii<lhcalevent->cellig.size();ii++)
			{
				maxpe = maxpe < lhcalevent->cellpe.at(ii) ? lhcalevent->cellpe.at(ii) : maxpe;
			}
			if(maxpe<300)
				continue;
//			wcda_ig->clear();
//			wcda_t->clear();
//			wcda_pe->clear();
//			wcda_x->clear();
//			wcda_y->clear();
			// ****reconstruction of WCDA events ** //
			// *** Get Reconstruction results by ZengZk ***//
			wcdareconstruct->Init();
			wcdareconstruct->SetWCDAEvent(lhcalevent->cellig, lhcalevent->cellpe, lhcalevent->cellt);
			//wcdareconstruct->SetWCDAEvent(*(lhcalevent->mcig),*(lhcalevent->mcq),*(lhcalevent->mct),*(lhcalevent->mcx),*(lhcalevent->mcy),*(lhcalevent->mcz));
			wcdareconstruct->AfterPulseClean();
			wcdareconstruct->TimeResidualClean();
			wcdareconstruct->CoreReconstruction();
			wcdareconstruct->DirectionReconstruction();
			wcdareconstruct->RecEnergy();
			wcdareconstruct->GetRecResult(evrec);
//			wcdareconstruct->GetHits(*wcda_x, *wcda_y, *wcda_pe, *wcda_ig, *wcda_t);

			wcdarecevent->wcdaMaxPe = wcdareconstruct->GetMaxPe();
			wcdarecevent->GetRecResult(evrec);

			//printf("%f %f\n",wcdarecevent->evrec[4],wcdarecevent->evrec[5]);
			// ***reconstruction of WFCTA events ** //            
			EventRec->Fill(); 
			if(0==i%1000)
				std::cerr << "iEntry: " << i << " / " << nentries << std::endl;
			if(i%100==0||i==nentries-1) printf("---------->Reconstructed event %d\n",i);
		}
		hfile->Close();
		delete calievent;
	}

	delete wcdaevent;
	delete lhcalevent;
	delete wcda_ig;
	delete wcda_t;
	delete wcda_pe;
	delete wcda_x;
	delete wcda_y;


	delete wcdarecevent;
	delete wcdareconstruct;

	newfile->Write();
	newfile->Close();
	return 0;
}

void Coordinate_Corsika2WCDA(std::vector<double>& cellx, std::vector<double>& celly)
{
	//double WCDA1CenterX = -76; //m 
	//double WCDA1CenterY = -56; //m
	WCDAMcRecEvent wfctamcrec;
	wfctamcrec.SetTelCenter(WFCTAX,WFCTAY);
	double azimuth;
	double wcdaazi = 0;
	double cell_x, cell_y;
	double cell_x_wcda, cell_y_wcda;
	for(int i=0;i<cellx.size();i++)
	{
		cell_x = cellx.at(i);
		cell_y = celly.at(i);
		wfctamcrec.Corsika2WCDA(azimuth, cell_x, cell_y, &wcdaazi, &cell_x_wcda, &cell_y_wcda);
		//cell_x_wcda += WCDA1CenterX;
		//cell_y_wcda += WCDA1CenterY;
		cellx.at(i) = cell_x_wcda;
		celly.at(i) = cell_y_wcda;
	}
}

