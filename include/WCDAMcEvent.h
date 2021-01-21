#ifndef __WCDAMCEVENT_HH__
#define __WCDAMCEVENT_HH__
#include "TObject.h"
#include <vector>
#include <map>
using namespace std;

class WCDAMcEvent : public TObject {
public:  
  int fileid;
  int runid;
  int eventid;
  int coreset;
  double e0;
  double mczen;
  double mcazi;
 // double mcazi1;
  double mcxc;
  double mcyc;
  double mthx;
  double mthy;
  int fitflag;
  int recflag;
  int nhit;
  int npmtfire;
  int nfitc;
  int ndetc;
  double zenc;
  double azic;
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
  double CXPE[10];
  
  vector<double>* mcq = new vector<double>();
  vector<double>* mct = new vector<double>();
  vector<int>* mcig = new vector<int>();
  vector<double>* mcx = new vector<double>(); 
  vector<double>* mcy = new vector<double>(); 
  vector<double>* mcz = new vector<double>(); 
  

  WCDAMcEvent();
  ~WCDAMcEvent();
  void SetEvent();
  ClassDef(WCDAMcEvent,1);
};

#endif

