#ifndef __WCDAMCRECEVENT_HH__
#define __WCDAMCRECEVENT_HH__
#include "TObject.h"
#include <vector>
#include <map>
#include "WCDAMcEvent.h"
#include "LHCALEvent.h"
#include "TMath.h"
using namespace std;

class WCDAMcRecEvent : public TObject {
private:
  double WFCTAX;//!
  double WFCTAY;//!
public:  
  int fileid;
  int runid;
  int eventid;
  int coreset;
  double e0;
  double mczen;
  double mcazi;
  double mcxc;
  double mcyc;
  double mcxc_inwcda;
  double mcyc_inwcda;
  double mthx;
  double mthy;
  double mthx_inwcda;
  double mthy_inwcda;
  int fitflag;
  int recflag;
  int nhit;
  int npmtfire;
  int nfitc;
  int ndetc;
  double zenc;
  double azic;
  double azic_inwcda;
  double recxc;
  double recyc;
  double recxc1;
  double recyc1;
  int npea;
  int npec;
  double dcore;
  double dcore1;
  double omega;
  double compactness;
  int npmtshower;
  int nnpmtnoise;
  int mynfitc;
  int nbig;
  double azi_zzk_incorsika;
  double recx_zzk_incorsika;
  double recy_zzk_incorsika;
 
  double evrec[70];
  double wcdaMaxPe; 
  
  void SetTelCenter(double wfctax, double wfctay);
  WCDAMcRecEvent();
  ~WCDAMcRecEvent();
  void SetEvent(WCDAMcEvent *wcdaevent);
  void SetEvent(LHCALEvent *lhacalevent);
  void GetRecResult(double evrec[70]);

  void WCDA2Corsika(double azimuth,double corex, double corey,double *corsikaazi,double *corsikax,double *corsikay);
  void Corsika2WCDA(double azimuth,double corex, double corey,double *wcdaazi,double *wcdax,double *wcday);

  ClassDef(WCDAMcRecEvent,1);
};

#endif

