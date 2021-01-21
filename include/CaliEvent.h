#ifndef CALIEVENT_H
#define CALIEVENT_H

#include "TObject.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"

#include "SaveEvent.h"
#include "LHCALEvent.h"

using namespace std;

class CaliEvent : public TObject {

  public:
    void GetWCDA_QT_Cali();
    void WCDA_Cali( SaveEvent *wcdaevent, LHCALEvent *lhcalevent);
};
#endif // CALIEVENT_H
