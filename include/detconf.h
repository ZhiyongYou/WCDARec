#include	<iostream>
#include  <vector>
#include	<stdlib.h>
#include	<fstream>
#include	<math.h>

#include	"TGraph.h"
#include	"TStyle.h"
#include	"TString.h"
#include	"TH1F.h"
#include	"TH2F.h"
#include	"TF1.h"
#include	"TF2.h"
#include	"TCanvas.h"
#include	"TMath.h"
#include	"TTree.h"
#include	"TFile.h"
#include 	"TSpectrum.h"
#include	"TLine.h"
#include	"wcdaconf.h"

#define	pi	3.1415926
#define	cspeed 0.2997
#define rad	57.3
using namespace std;
//	vars

//	wcdaconf 
//	FEE channel x[m] y[m] scc g d igcell  z[m] dx[m] dy[m] dz[m] //iconf PMT IP

static	TString	fname;
static	TString	mode;
static	const	TString	single	=	"single";
static	const	TString	event		=	"event";
static	const	int	mask				=	-1;

static	int	ipool;
static	int	igcell;
static	int	sccgd;
static	int	s,cc,g,d;
static	int	id;
static	UInt_t  second0, coarse_time0, low_th_fine_time0, high_th_fine_time0;
static	UInt_t  second, coarse_time, low_th_fine_time, high_th_fine_time;
static	Float_t	x,y,z;
static	int			maxbin;
static	TF1*		func;
static	TF1*		fq;
static	TF1*		ft;


//	vectors
static	vector<Float_t>	vq[nwcdaconf];
static	vector<Float_t>	vtl[nwcdaconf];
static	vector<Float_t>	vth[nwcdaconf];
static	vector<Float_t>	vx[nwcdaconf];
static	vector<Float_t>	vy[nwcdaconf];

//	arrays
static	Int_t		vuse[nwcdaconf] = {0};

//	calibration info
//	qcal
static	TString	qcalfile;
static	Float_t	spe[nwcdaconf]	= {0.};
static	Float_t	k[nwcdaconf]		=	{0.};
static	Float_t	b[nwcdaconf]		=	{0.};
int loadqcal(TString qcalfile)	{
  TTree* t = new TTree();
  t->ReadFile(qcalfile,"iigcell/I:gain/D:ad/D");
  int     iigcell;
  double  gain,ad;
  t->SetBranchAddress("iigcell",&iigcell);
  t->SetBranchAddress("gain",&gain);
  t->SetBranchAddress("ad",&ad);
  for (int ix=0;ix<t->GetEntries();ix++)  {
    t->GetEntry(ix);
	  spe[iigcell] = gain; k[iigcell] = ad; 
  }
	return 0;
}

//	tcal
static	TString	tcalfile;

//	for drawing & fitting
static	TH2F*	h2rate;
static	TH2F*	h2spe;
static	TH2F*	h2dar;
static	TH2F*	h2mueff;

//	charge	calibration : already done by Zhiguo.Yao's program
static	int	used	=	-1;

//	event	data,	offline / hardware cailibration
static	Float_t	tmean;
static	Float_t	tmeanl,tmeanh;
static	Float_t	toffset;
static	Float_t	tresidual;
static	Float_t	qtpar[3]	=	{0.,0.,0.};	//	fitting function : T(q) = [0]*pow(q,[1]) + [2];
static	Float_t	tinl[nwcdaconf]	=	{-1000.};
static	Float_t	tcrl[nwcdaconf]	=	{-1000.};
static	Float_t	tinh[nwcdaconf]	=	{-1000.};
static	Float_t	tcrh[nwcdaconf]	=	{-1000.};
static	Float_t	qin[nwcdaconf]	=	{-1000.};
static	Float_t	qcr[nwcdaconf]	=	{-1000.};
static	Float_t	q[nwcdaconf]		=	{-1000.};
static	Float_t	t[nwcdaconf]		=	{-1000.};
static	Float_t	tl[nwcdaconf]		=	{-1000.};
static	Float_t	th[nwcdaconf]		=	{-1000.};
static	Float_t	itl;
static	Float_t	ith;
static	Float_t qq,tt; // tmp varibles



static	TF1*	fqt;
static	TF2*	fplane;
static	Float_t	theta,phi;
static	Float_t	evt_time;

// static	TH1F*	hq;
static	TH1F*	hq;
static	TH1F*	htl;
static	TH1F*	hth;
static	TH1F*	htlh;
static	TH1F*	htl0;
static	TH1F*	hth0;
static	TH1F*	htlall;
static	TH1F*	hthall;

static	TH1F*	hnhit;
static	int	nhit_mean,nhit_min,nhit_max,nhit_rms;
static	TH2F*	h2ledq;
static	TH2F*	h2toff;

static  TString _2dxtit   = "X [m]";
static  TString _2dytit   = "Y [m]";
static  TString _qhtit    = "HG ADC";
static  TString _qltit    = "LG ADC";
static  TString _ytit     = "Entries";
static  TString _qtit     = "NPE";
static  TString _ttit     = "Time [ns]";
static  TString _rtit     = "Rate [kHz]";
static  TString _ctit     = "igcell";

//	funcs
double	mint(double v[4])	{
	double min = 1000;
  double e = 0;
  double e1,e2,e3,e4;
  double x;

  if (fabs(v[3]-v[0])<1.) {
    x = (v[0]+v[1]+v[2]+v[3])/4.;
    e1 = v[0] - x;
    e2 = v[1] - x;
    e3 = v[2] - x;
    e4 = v[3] - x;
    e = sqrt(e1*e1+e2*e2+e3*e3+e4*e4)/2.;
  }

  if (fabs(v[3]-v[0])>1.&&(fabs(v[3]-v[1])<1.||fabs(v[2]-v[0])<1.)) {
    for (int i=0;i<4;i++) {
      for (int j=i+2;j<4;j++) {
        if (min>fabs(v[j]-v[i]))  {
          min = fabs(v[j]-v[i]);
            x = (v[j]+v[i]+v[j-1])/3.;
          e1 = v[j]-x;
          e2 = v[j-1]-x;
          e3 = v[i]-x;
          e = sqrt(e1*e1+e2*e2+e3*e3)/sqrt(3.);
        }
      }
    }
  }
  if (fabs(v[0]-v[2])>1.||fabs(v[1]-v[3])>1.) {
    for (int i=0;i<4;i++) {
      for (int j=i+1;j<4;j++) {
        if (min>fabs(v[j]-v[i]))  {
          min = fabs(v[j]-v[i]);
          x = (v[j]+v[i])/2.;
          e1 = v[j]-x;
          e2 = v[i]-x;
          e = sqrt(e1*e1+e2*e2)/sqrt(2.);
        }
      }
    }
  }
  return  e;
}

//	qt correction
/*
double	qtc(int icc,double qi,double q0)	{
  double  p0,p1,p2;
  p0 = qtpar[icc][0];
  p1 = qtpar[icc][1];
  p2 = qtpar[icc][2];
  return
    p0*TMath::Exp(p1*qi)*TMath::Power(qi,p2) - p0*TMath::Exp(p1*q0)*TMath::Power(q0,p2);
}
*/

int findid(UInt_t ifee,UInt_t ich)  {
  id = 0;
  int i=0;
  int j=0;
  igcell = 0;
	UInt_t xfee = 0;
	if (ich > 9) return -1;
  //  wcda1
  if (ifee < 1000)  {
      id = (ifee-1)*9+ich-1; 
			igcell = wcdaconf[id][7];
			x	=	wcdaconf[id][2];
			y	=	wcdaconf[id][3];
			z	=	wcdaconf[id][8];
     // cout << ifee << " | " << ich << " = " << id << " " << igcell << endl;
      return id;
    }
  //  wcda2
  if (ifee >= 1000 && ifee < 2000)  {
      id = (ifee%1000-1)*9+900+ich-1; 
			igcell = wcdaconf[id][7];
     // cout << ifee << " | " << ich << " = " << id << " " << igcell << endl;
			x	=	wcdaconf[id][2];
			y	=	wcdaconf[id][3];
			z	=	wcdaconf[id][8];
      return id;
  }
  //  wcda3
  if (ifee >= 2001) {
    xfee  = ifee%2000; 
		i = xfee/16; 
		j = xfee%16;
    if (j==0) return -1;
    if (j<14&&j!=0)   { 
			id = i*132+(j-1)*9+ich-1+1800; 
			igcell = wcdaconf[id][7]; 
			x	=	wcdaconf[id][2];
			y	=	wcdaconf[id][3];
			z	=	wcdaconf[id][8];
		}
    if (j==14)  {
      if (ich<=3) { 
				id = i*132+(j-1)*9+ich-1+1800; 
				igcell = wcdaconf[id][7];  
				x	=	wcdaconf[id][2];
				y	=	wcdaconf[id][3];
				z	=	wcdaconf[id][8];
			}
      if (ich>3)  { 
				id = i*132+(j-1)*9+ich-4+1800; 
				igcell = wcdaconf[id][7];  
				x	=	wcdaconf[id][2];
				y	=	wcdaconf[id][3];
				z	=	wcdaconf[id][8];
			}
    }
    if (j>14)   { 
				id = i*132+(j-1)*9+ich-4+1800; 
				igcell = wcdaconf[id][7];    
				x	=	wcdaconf[id][2];
				y	=	wcdaconf[id][3];
				z	=	wcdaconf[id][8];
		}
  	//	cout << ifee << " | " << ich << " = " << id << " " << igcell << endl;
  	return id;
  }
}

int findig(int sccgd)	{
	int s;
	int cc;
	int g;
	int d;

	s		= sccgd/10000;
	cc	= sccgd%10000/100;
	g		= sccgd%100/10;
	d		= sccgd%10;
	int ig = 0;
	ig	=	d+g*9+cc*36;
	return ig;
}

int findig(int ifee,int ich)  {
  id = 0;
  int i=0;
  int j=0;
  igcell = 0;
  UInt_t xfee = 0;
  if (ich > 9) return -1;
  //  wcda1
  if (ifee < 1000)  {
    id = (ifee-1)*9+ich-1;
    igcell = wcdaconf[id][7];
    return igcell;
  }
  //  wcda2
  if (ifee >= 1000 && ifee < 2000)  {
    id = (ifee%1000-1)*9+900+ich-1;
    igcell = wcdaconf[id][7];
    return igcell;
  } 
  //  wcda3
  if (ifee >= 2001) {
    xfee  = ifee%2000;
    i = xfee/16;
    j = xfee%16;
    if (j==0) return -1;

    if (j<14&&j!=0)   {
      id = i*132+(j-1)*9+ich-1+1800;
      igcell = wcdaconf[id][7];
    }
    if (j==14)  {
      if (ich<=3) {
        id = i*132+(j-1)*9+ich-1+1800;
        igcell = wcdaconf[id][7];
      }
      if (ich>3)  {
        id = i*132+(j-1)*9+ich-4+1800;
        igcell = wcdaconf[id][7];
      }
    }
    if (j>14)   {
      id = i*132+(j-1)*9+ich-4+1800;
      igcell = wcdaconf[id][7];
    }

    return igcell;
  }
}

int findplusig(int fee,int db, int pmt)	{
	int id = fee*1000+db*100+pmt;
	if (fee<=25)	{
			for (int i=0;i<900;i++)	{
				if (plusconf[i][0] == id)	{
					sccgd	=	plusconf[i][2];
					igcell = sccgd%10000/100*36+sccgd%100/10*9+sccgd%10;
					return igcell;
					break;
				}
			}
	}
	if (fee>1000&&fee<2000)	{
			for (int i=900;i<1800;i++)	{
				if (plusconf[i][0] == id)	{
					sccgd	=	plusconf[i][2];
					igcell = sccgd%10000/100*36+sccgd%100/10*9+sccgd%10 + 900;
					return igcell;
					break;
				}
			}
	}
	if (fee>=2000)	{
			for (int i=1800;i<3120;i++)	{
				if (plusconf[i][0] == id)	{
					sccgd	=	plusconf[i][2];
					igcell = sccgd%10000/100*36+sccgd%100/10*9+sccgd%10 + 1800;
					return igcell;
					break;
				}
			}
	}
}

int initled	()	{
	tinl[nwcdaconf] = {-1000.};
	tcrl[nwcdaconf] = {-1000.};
	tinh[nwcdaconf] = {-1000.};
	tcrh[nwcdaconf] = {-1000.};
	qin[nwcdaconf]	= {-1000.};
	qcr[nwcdaconf]	= {-1000.};
	return 1;
}

float qpeakfit(TH1F* h,Int_t i)	{ // h is hq, i is igcell [nwcdaconf]
	if (h->GetEntries() <= 100) {vuse[i] = 0;return -1000.;}
	if (h->GetEntries() >= 100)	{
		h->GetXaxis()->SetRangeUser(5,100);
		maxbin	= h->GetMaximumBin();
		double mean		=	h->GetBinCenter(maxbin);
		double p1,p2,p3;
		//	double meana		=	h->GetMean();
		//	double m;
		//	if (fabs(mean-meana) < 100 ) m = mean;
		//	if (fabs(mean-meana) > 100 ) m = meana;
		//	fq = new TF1("fq","expo+gaus(2)");
		//	fq->SetParameters(5.,-0.01,1000,m,m*0.1);
		//	fq = new TF1("fq","gaus");
		h->Fit("gaus","Q","",mean-10,mean+10);
		fq	=	h->GetFunction("gaus");
		p1	=	fq->GetParameter(0);
		p2	=	fq->GetParameter(1);
		p3	=	fq->GetParameter(2);
		h->Fit("gaus","Q","",p2-p3,p2+p3);
		fq	=	h->GetFunction("gaus");
		p1	=	fq->GetParameter(0);
		p2	=	fq->GetParameter(1);
		p3	=	fq->GetParameter(2);
		h->Fit("gaus","Q","",p2-p3,p2+p3);
		fq	=	h->GetFunction("gaus");
		p1	=	fq->GetParameter(0);
		p2	=	fq->GetParameter(1);
		p3	=	fq->GetParameter(2);
		h->Fit("gaus","Q","",p2-p3,p2+p3);
		fq	=	h->GetFunction("gaus");
		p1	=	fq->GetParameter(0);
		p2	=	fq->GetParameter(1);
		p3	=	fq->GetParameter(2);

		return float(p2);
	}
}


float tpeakfit(TH1F* h,Int_t i)	{ // h is hq, i is igcell [nwcdaconf]
	if (h->GetEntries() <= 1000) {return -1000.;}
	int n;
	if (h->GetEntries() >= 1000)	{
		ft = new TF1();
		double peak,sigma;
		maxbin=	h->GetMaximumBin();
	  peak 	= h->GetBinCenter(maxbin);
	  n 	= h->GetBinContent(maxbin);
		if (n<200) {return -1000.;}
		h->Fit("gaus","Q","",peak-2.5,peak+2.5);
		ft		=	h->GetFunction("gaus");
		peak 	= ft->GetParameter(1);
		sigma	=	ft->GetParameter(2);	
		h->Fit("gaus","Q","",peak-sigma,peak+sigma);
		ft		=	h->GetFunction("gaus");
		peak 	= ft->GetParameter(1);
		sigma	=	ft->GetParameter(2);	
		h->Fit("gaus","Q","",peak-sigma,peak+sigma);
		ft		=	h->GetFunction("gaus");
		peak 	= ft->GetParameter(1);
		vuse[i] = 1;
		return float(peak);
	}
}

/*
//  find the peakx of led nhit
static  double* _peak_x;
static  double* _peak_y;
double peakx(TH1D* h) {
  h->GetXaxis()->SetRangeUser(50,2000);
  double    rms   = h->GetRMS();
  TSpectrum* sp   = new TSpectrum();
  const int npeak = sp->Search(h,0.1*rms,"",0.01);
  int v[npeak];
  _peak_x         = sp->GetPositionX();
  _peak_y         = sp->GetPositionY();
  TMath::Sort(npeak,_peak_y,v,kTRUE);
  return  _peak_x[v[0]];
}
*/
