#ifndef __LHCALEVENT_HH__
#define __LHCALEVENT_HH__
#include "TObject.h"

class LHCALEvent : public TObject {

public:
//wcda
  std::vector<int> cellig;
  std::vector<double> cellpe;
  std::vector<double> cellt;

  int nhitplus;
  unsigned long int iwcdaevt;
  int wcda_trigerflag;
  float wcda_phi;
  float wcda_theta;

  LHCALEvent();
  ~LHCALEvent();
  void Init();

  ClassDef(LHCALEvent,1);
};

#endif

