#include "LHCALEvent.h"

ClassImp(LHCALEvent);

LHCALEvent::LHCALEvent()
{
}

LHCALEvent::~LHCALEvent()
{
}

void LHCALEvent::Init()
{
  nhitplus = -1;
  iwcdaevt = 0;
  wcda_trigerflag = -9;
  wcda_phi = -10.;
  wcda_theta = -10.;
  cellig.clear();
  cellpe.clear();
  cellt.clear();
}
