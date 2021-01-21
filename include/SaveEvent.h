#ifndef _SaveEvent
#define _SaveEvent 1


#include "TObject.h"
#include "TClonesArray.h"
#include "iostream"

//define it, for compatible reason
#define TREE_WITH_LONG_TIME

class SaveEvent : public TObject {
public:
  SaveEvent();
  virtual ~SaveEvent();

  Int_t    ievent;        //Event number of this event
  Int_t    flag;          //flag of this event
  Float_t  theta;         //Zenith of this event
  Float_t  phi;           //Azimuth of this event
  Long64_t sec;           //second of this event
  Double_t subsec;        //Subsecond of this event

  Int_t nhitpluso;        //Original number of PLUS hits
  Int_t nhitplus;         //Number of PLUS hits
  TClonesArray *hitsplus; //->array of all the PLUS hits

  Int_t nhito;            //Original number of WCDA hits
  Int_t nhit;             //Number of WCDA hits
  TClonesArray *hits;     //->array of all the WCDA hits

  // To add a hit ...
  void AddHit(
    Int_t fee,
    Int_t ch,
    Int_t anode_charge,
    Int_t dynode_charge,
    Int_t charge_over_range_id,
    Int_t delta_time,
    Int_t delta_time_hl
  );

  void AddHitPlus(
    Int_t fee,
    Int_t flag,
    Int_t db,
    Int_t pmt,
    Int_t anode_peak,
    Int_t anode_ped,
    Int_t anode_time,
    Int_t dynode_peak,
    Int_t dynode_ped,
    Int_t dynode_time,  
    Int_t delta_time
  );

  void Init() {
    nhitpluso = 0;
    hitsplus->Clear();
    nhitplus = 0;

    nhito = 0;
    hits->Clear();
    nhit = 0;

    ievent = 0;
    flag = 0;
    theta = 0;
    phi = 0;
    sec = 0;
    subsec = 0;
  };

  ClassDef(SaveEvent,1) 
};


class SaveWCDAHit : public TObject {
private:
  Int_t fee;
  Int_t ch;
  Int_t anode_charge;
  Int_t dynode_charge;
  Int_t charge_over_range_id;
  //in picosecond, delta time with the event time (time of a median hit)
  Int_t delta_time;
  Int_t delta_time_hl;

public:
  // constructor and distructor
  SaveWCDAHit() {};
  SaveWCDAHit(
    Int_t fee,
    Int_t ch,
    Int_t anode_charge,
    Int_t dynode_charge,
    Int_t charge_over_range_id,
    Int_t delta_time,
    Int_t delta_time_hl
  );
 
  virtual ~SaveWCDAHit() {};

  void update(
    Int_t fee,
    Int_t ch,
    Int_t anode_charge,
    Int_t dynode_charge,
    Int_t charge_over_range_id,
    Int_t delta_time,
    Int_t delta_time_hl
  );

public:
  void SetFee(Int_t val) { this->fee = val; }
  Int_t Fee() { return this->fee; }

  void SetCh(Int_t val) { this->ch = val; }
  Int_t Ch() { return this->ch; }

  void SetAnodeCharge(Int_t val) { this->anode_charge = val; }
  Int_t AnodeCharge() { return this->anode_charge; }

  void SetDynodeCharge(Int_t val) { this->dynode_charge = val; }
  Int_t DynodeCharge() { return this->dynode_charge; }

  void SetChargeOverRangeId(Int_t val) { this->charge_over_range_id = val; }
  Int_t ChargeOverRangeId() { return this->charge_over_range_id; }

  void SetDeltaTime(Int_t val) { this->delta_time = val; }
  Int_t DeltaTime() { return this->delta_time; }

  void SetDeltaTimeHL(Int_t val) { this->delta_time_hl = val; }
  Int_t DeltaTimeHL() { return this->delta_time_hl; }

  ClassDef(SaveWCDAHit,1) 
};


class SavePLUSHit : public TObject {
private:
  Int_t fee;
  Int_t flag;
  Int_t db;
  Int_t pmt;

  Int_t anode_peak;
  Int_t anode_ped;
  Int_t anode_time;
  Int_t dynode_peak;
  Int_t dynode_ped;
  Int_t dynode_time;

  //in picosecond, delta time with the event time (time of a median hit)
  Int_t delta_time;

public:
  // constructor and distructor
  SavePLUSHit() {};
  SavePLUSHit(
    Int_t fee,
    Int_t flag,
    Int_t db,
    Int_t pmt,
    Int_t anode_peak,
    Int_t anode_ped,
    Int_t anode_time,
    Int_t dynode_peak,
    Int_t dynode_ped,
    Int_t dynode_time,
    Int_t delta_time
  );

  virtual ~SavePLUSHit() {};

  void update(
    Int_t fee,
    Int_t flag,
    Int_t db,
    Int_t pmt,
    Int_t anode_peak,
    Int_t anode_ped,
    Int_t anode_time,
    Int_t dynode_peak,
    Int_t dynode_ped,
    Int_t dynode_time,
    Int_t delta_time
  );
 
  void SetFee(Int_t fee) { this->fee = fee; }
  Int_t Fee() { return fee; }

  void SetFlag(Int_t flag) { this->flag = flag; }
  Int_t Flag() { return flag; }

  void SetDb(Int_t db) { this->db = db; }
  Int_t Db() { return db; }

  void SetPmt(Int_t pmt) { this->pmt = pmt; } 
  Int_t Pmt() { return pmt; }

  void SetDeltaTime(Int_t val) { this->delta_time = val; }
  Int_t DeltaTime() { return delta_time; }

  void SetAnodePeak(Int_t anode_peak) { this->anode_peak = anode_peak; }
  Int_t AnodePeak() { return anode_peak; }

  void SetAnodePed(Int_t anode_ped) { this->anode_ped = anode_ped; } 
  Int_t AnodePed() { return anode_ped; } 

  void SetAnodeTime(UInt_t anode_time) { this->anode_time = anode_time; }
  UInt_t AnodeTime() { return anode_time; }

  void SetDynodePeak(UInt_t dynode_peak) { this->dynode_peak = dynode_peak; }
  UInt_t DynodePeak() { return dynode_peak; }

  void SetDynodePed(UInt_t dynode_ped) { this->dynode_ped = dynode_ped; }
  UInt_t DynodePed() { return dynode_ped; }

  void SetDynodeTime(UInt_t dynode_time) { this->dynode_time = dynode_time; }
  UInt_t DynodeTime() { return dynode_time; }

  ClassDef(SavePLUSHit,1)
};
#endif
