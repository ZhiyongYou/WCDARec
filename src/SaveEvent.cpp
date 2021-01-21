#include <math.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include "SaveEvent.h"
#include <iostream>

ClassImp(SaveEvent)
ClassImp(SaveWCDAHit)
ClassImp(SavePLUSHit)

SaveEvent::SaveEvent() {
  nhit = 0;
  nhito = 0;
  hits = new TClonesArray("SaveWCDAHit",8000);
  hits->Clear();

  nhitplus = 0;
  nhitpluso = 0;
  hitsplus = new TClonesArray("SavePLUSHit",8000);
  hitsplus->Clear();
}


SaveEvent::~SaveEvent() {
  hits->Clear(); 
  delete  hits;

  hitsplus->Clear(); 
  delete  hitsplus;
}


void SaveEvent::AddHit(
    Int_t fee,
    Int_t ch,
    Int_t anode_charge,
    Int_t dynode_charge,
    Int_t charge_over_range_id,
    Int_t delta_time,
    Int_t delta_time_hl) {
  new((*hits)[nhit++]) SaveWCDAHit(
    fee,
    ch,
    anode_charge,
    dynode_charge,
    charge_over_range_id,
    delta_time,
    delta_time_hl
  );
}


void SaveEvent::AddHitPlus(
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
    Int_t delta_time) {
  new((*hitsplus)[nhitplus++]) SavePLUSHit(
    fee,
    flag,
    db,
    pmt,
    anode_peak,
    anode_ped,
    anode_time,
    dynode_peak,
    dynode_ped,
    dynode_time,
    delta_time
  );
}


SaveWCDAHit::SaveWCDAHit(
    Int_t fee,
    Int_t ch,
    Int_t anode_charge,
    Int_t dynode_charge,
    Int_t charge_over_range_id,
    Int_t delta_time,
    Int_t delta_time_hl) {
  update(
    fee,
    ch,
    anode_charge,
    dynode_charge,
    charge_over_range_id,
    delta_time,
    delta_time_hl
  );
}


SavePLUSHit::SavePLUSHit(
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
    Int_t delta_time) {
  update(
    fee,
    flag,
    db,
    pmt,
    anode_peak,
    anode_ped,
    anode_time,
    dynode_peak,
    dynode_ped,
    dynode_time,
    delta_time
  );
}


void SaveWCDAHit::update(
    Int_t fee,
    Int_t ch,
    Int_t anode_charge,
    Int_t dynode_charge,
    Int_t charge_over_range_id,
    Int_t delta_time,
    Int_t delta_time_hl) {
  this->fee = fee;
  this->ch = ch;
  this->anode_charge = anode_charge;
  this->dynode_charge = dynode_charge;
  this->charge_over_range_id = charge_over_range_id;
  this->delta_time = delta_time;
  this->delta_time_hl = delta_time_hl;
}


void SavePLUSHit::update(
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
    Int_t delta_time) {
  this->fee = fee;
  this->flag = flag;
  this->db = db;
  this->pmt = pmt;
  this->delta_time = delta_time;
  this->anode_peak = anode_peak;
  this->anode_ped = anode_ped;
  this->anode_time = anode_time;
  this->dynode_peak = dynode_peak;
  this->dynode_ped = dynode_ped;
  this->dynode_time = dynode_time;
}
