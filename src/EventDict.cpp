// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIsrcdIEventDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "WCDAMcEvent.h"
#include "WCDAMcRecEvent.h"
#include "LHCALEvent.h"
#include "SaveEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_WCDAMcEvent(void *p = 0);
   static void *newArray_WCDAMcEvent(Long_t size, void *p);
   static void delete_WCDAMcEvent(void *p);
   static void deleteArray_WCDAMcEvent(void *p);
   static void destruct_WCDAMcEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCDAMcEvent*)
   {
      ::WCDAMcEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCDAMcEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCDAMcEvent", ::WCDAMcEvent::Class_Version(), "WCDAMcEvent.h", 8,
                  typeid(::WCDAMcEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCDAMcEvent::Dictionary, isa_proxy, 4,
                  sizeof(::WCDAMcEvent) );
      instance.SetNew(&new_WCDAMcEvent);
      instance.SetNewArray(&newArray_WCDAMcEvent);
      instance.SetDelete(&delete_WCDAMcEvent);
      instance.SetDeleteArray(&deleteArray_WCDAMcEvent);
      instance.SetDestructor(&destruct_WCDAMcEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCDAMcEvent*)
   {
      return GenerateInitInstanceLocal((::WCDAMcEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WCDAMcEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_LHCALEvent(void *p = 0);
   static void *newArray_LHCALEvent(Long_t size, void *p);
   static void delete_LHCALEvent(void *p);
   static void deleteArray_LHCALEvent(void *p);
   static void destruct_LHCALEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LHCALEvent*)
   {
      ::LHCALEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LHCALEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LHCALEvent", ::LHCALEvent::Class_Version(), "LHCALEvent.h", 5,
                  typeid(::LHCALEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LHCALEvent::Dictionary, isa_proxy, 4,
                  sizeof(::LHCALEvent) );
      instance.SetNew(&new_LHCALEvent);
      instance.SetNewArray(&newArray_LHCALEvent);
      instance.SetDelete(&delete_LHCALEvent);
      instance.SetDeleteArray(&deleteArray_LHCALEvent);
      instance.SetDestructor(&destruct_LHCALEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LHCALEvent*)
   {
      return GenerateInitInstanceLocal((::LHCALEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LHCALEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCDAMcRecEvent(void *p = 0);
   static void *newArray_WCDAMcRecEvent(Long_t size, void *p);
   static void delete_WCDAMcRecEvent(void *p);
   static void deleteArray_WCDAMcRecEvent(void *p);
   static void destruct_WCDAMcRecEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCDAMcRecEvent*)
   {
      ::WCDAMcRecEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCDAMcRecEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCDAMcRecEvent", ::WCDAMcRecEvent::Class_Version(), "WCDAMcRecEvent.h", 11,
                  typeid(::WCDAMcRecEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCDAMcRecEvent::Dictionary, isa_proxy, 4,
                  sizeof(::WCDAMcRecEvent) );
      instance.SetNew(&new_WCDAMcRecEvent);
      instance.SetNewArray(&newArray_WCDAMcRecEvent);
      instance.SetDelete(&delete_WCDAMcRecEvent);
      instance.SetDeleteArray(&deleteArray_WCDAMcRecEvent);
      instance.SetDestructor(&destruct_WCDAMcRecEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCDAMcRecEvent*)
   {
      return GenerateInitInstanceLocal((::WCDAMcRecEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WCDAMcRecEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SaveEvent(void *p = 0);
   static void *newArray_SaveEvent(Long_t size, void *p);
   static void delete_SaveEvent(void *p);
   static void deleteArray_SaveEvent(void *p);
   static void destruct_SaveEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SaveEvent*)
   {
      ::SaveEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SaveEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SaveEvent", ::SaveEvent::Class_Version(), "SaveEvent.h", 12,
                  typeid(::SaveEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SaveEvent::Dictionary, isa_proxy, 4,
                  sizeof(::SaveEvent) );
      instance.SetNew(&new_SaveEvent);
      instance.SetNewArray(&newArray_SaveEvent);
      instance.SetDelete(&delete_SaveEvent);
      instance.SetDeleteArray(&deleteArray_SaveEvent);
      instance.SetDestructor(&destruct_SaveEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SaveEvent*)
   {
      return GenerateInitInstanceLocal((::SaveEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SaveEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SaveWCDAHit(void *p = 0);
   static void *newArray_SaveWCDAHit(Long_t size, void *p);
   static void delete_SaveWCDAHit(void *p);
   static void deleteArray_SaveWCDAHit(void *p);
   static void destruct_SaveWCDAHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SaveWCDAHit*)
   {
      ::SaveWCDAHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SaveWCDAHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SaveWCDAHit", ::SaveWCDAHit::Class_Version(), "SaveEvent.h", 78,
                  typeid(::SaveWCDAHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SaveWCDAHit::Dictionary, isa_proxy, 4,
                  sizeof(::SaveWCDAHit) );
      instance.SetNew(&new_SaveWCDAHit);
      instance.SetNewArray(&newArray_SaveWCDAHit);
      instance.SetDelete(&delete_SaveWCDAHit);
      instance.SetDeleteArray(&deleteArray_SaveWCDAHit);
      instance.SetDestructor(&destruct_SaveWCDAHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SaveWCDAHit*)
   {
      return GenerateInitInstanceLocal((::SaveWCDAHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SaveWCDAHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SavePLUSHit(void *p = 0);
   static void *newArray_SavePLUSHit(Long_t size, void *p);
   static void delete_SavePLUSHit(void *p);
   static void deleteArray_SavePLUSHit(void *p);
   static void destruct_SavePLUSHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SavePLUSHit*)
   {
      ::SavePLUSHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SavePLUSHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SavePLUSHit", ::SavePLUSHit::Class_Version(), "SaveEvent.h", 140,
                  typeid(::SavePLUSHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SavePLUSHit::Dictionary, isa_proxy, 4,
                  sizeof(::SavePLUSHit) );
      instance.SetNew(&new_SavePLUSHit);
      instance.SetNewArray(&newArray_SavePLUSHit);
      instance.SetDelete(&delete_SavePLUSHit);
      instance.SetDeleteArray(&deleteArray_SavePLUSHit);
      instance.SetDestructor(&destruct_SavePLUSHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SavePLUSHit*)
   {
      return GenerateInitInstanceLocal((::SavePLUSHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SavePLUSHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr WCDAMcEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCDAMcEvent::Class_Name()
{
   return "WCDAMcEvent";
}

//______________________________________________________________________________
const char *WCDAMcEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCDAMcEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCDAMcEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCDAMcEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LHCALEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LHCALEvent::Class_Name()
{
   return "LHCALEvent";
}

//______________________________________________________________________________
const char *LHCALEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHCALEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LHCALEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LHCALEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LHCALEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHCALEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LHCALEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LHCALEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCDAMcRecEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCDAMcRecEvent::Class_Name()
{
   return "WCDAMcRecEvent";
}

//______________________________________________________________________________
const char *WCDAMcRecEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcRecEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCDAMcRecEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcRecEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCDAMcRecEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcRecEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCDAMcRecEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCDAMcRecEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SaveEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SaveEvent::Class_Name()
{
   return "SaveEvent";
}

//______________________________________________________________________________
const char *SaveEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SaveEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SaveEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SaveEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SaveEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SaveEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SaveEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SaveEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SaveWCDAHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SaveWCDAHit::Class_Name()
{
   return "SaveWCDAHit";
}

//______________________________________________________________________________
const char *SaveWCDAHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SaveWCDAHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SaveWCDAHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SaveWCDAHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SaveWCDAHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SaveWCDAHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SaveWCDAHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SaveWCDAHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SavePLUSHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SavePLUSHit::Class_Name()
{
   return "SavePLUSHit";
}

//______________________________________________________________________________
const char *SavePLUSHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SavePLUSHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SavePLUSHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SavePLUSHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SavePLUSHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SavePLUSHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SavePLUSHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SavePLUSHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void WCDAMcEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCDAMcEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCDAMcEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCDAMcEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCDAMcEvent(void *p) {
      return  p ? new(p) ::WCDAMcEvent : new ::WCDAMcEvent;
   }
   static void *newArray_WCDAMcEvent(Long_t nElements, void *p) {
      return p ? new(p) ::WCDAMcEvent[nElements] : new ::WCDAMcEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCDAMcEvent(void *p) {
      delete ((::WCDAMcEvent*)p);
   }
   static void deleteArray_WCDAMcEvent(void *p) {
      delete [] ((::WCDAMcEvent*)p);
   }
   static void destruct_WCDAMcEvent(void *p) {
      typedef ::WCDAMcEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCDAMcEvent

//______________________________________________________________________________
void LHCALEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class LHCALEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LHCALEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(LHCALEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LHCALEvent(void *p) {
      return  p ? new(p) ::LHCALEvent : new ::LHCALEvent;
   }
   static void *newArray_LHCALEvent(Long_t nElements, void *p) {
      return p ? new(p) ::LHCALEvent[nElements] : new ::LHCALEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_LHCALEvent(void *p) {
      delete ((::LHCALEvent*)p);
   }
   static void deleteArray_LHCALEvent(void *p) {
      delete [] ((::LHCALEvent*)p);
   }
   static void destruct_LHCALEvent(void *p) {
      typedef ::LHCALEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LHCALEvent

//______________________________________________________________________________
void WCDAMcRecEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCDAMcRecEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCDAMcRecEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCDAMcRecEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCDAMcRecEvent(void *p) {
      return  p ? new(p) ::WCDAMcRecEvent : new ::WCDAMcRecEvent;
   }
   static void *newArray_WCDAMcRecEvent(Long_t nElements, void *p) {
      return p ? new(p) ::WCDAMcRecEvent[nElements] : new ::WCDAMcRecEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCDAMcRecEvent(void *p) {
      delete ((::WCDAMcRecEvent*)p);
   }
   static void deleteArray_WCDAMcRecEvent(void *p) {
      delete [] ((::WCDAMcRecEvent*)p);
   }
   static void destruct_WCDAMcRecEvent(void *p) {
      typedef ::WCDAMcRecEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCDAMcRecEvent

//______________________________________________________________________________
void SaveEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class SaveEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SaveEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(SaveEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SaveEvent(void *p) {
      return  p ? new(p) ::SaveEvent : new ::SaveEvent;
   }
   static void *newArray_SaveEvent(Long_t nElements, void *p) {
      return p ? new(p) ::SaveEvent[nElements] : new ::SaveEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_SaveEvent(void *p) {
      delete ((::SaveEvent*)p);
   }
   static void deleteArray_SaveEvent(void *p) {
      delete [] ((::SaveEvent*)p);
   }
   static void destruct_SaveEvent(void *p) {
      typedef ::SaveEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SaveEvent

//______________________________________________________________________________
void SaveWCDAHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class SaveWCDAHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SaveWCDAHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(SaveWCDAHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SaveWCDAHit(void *p) {
      return  p ? new(p) ::SaveWCDAHit : new ::SaveWCDAHit;
   }
   static void *newArray_SaveWCDAHit(Long_t nElements, void *p) {
      return p ? new(p) ::SaveWCDAHit[nElements] : new ::SaveWCDAHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_SaveWCDAHit(void *p) {
      delete ((::SaveWCDAHit*)p);
   }
   static void deleteArray_SaveWCDAHit(void *p) {
      delete [] ((::SaveWCDAHit*)p);
   }
   static void destruct_SaveWCDAHit(void *p) {
      typedef ::SaveWCDAHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SaveWCDAHit

//______________________________________________________________________________
void SavePLUSHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class SavePLUSHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SavePLUSHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(SavePLUSHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SavePLUSHit(void *p) {
      return  p ? new(p) ::SavePLUSHit : new ::SavePLUSHit;
   }
   static void *newArray_SavePLUSHit(Long_t nElements, void *p) {
      return p ? new(p) ::SavePLUSHit[nElements] : new ::SavePLUSHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_SavePLUSHit(void *p) {
      delete ((::SavePLUSHit*)p);
   }
   static void deleteArray_SavePLUSHit(void *p) {
      delete [] ((::SavePLUSHit*)p);
   }
   static void destruct_SavePLUSHit(void *p) {
      typedef ::SavePLUSHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SavePLUSHit

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 216,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 216,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_EventDict_Impl() {
    static const char* headers[] = {
"WCDAMcEvent.h",
"WCDAMcRecEvent.h",
"LHCALEvent.h",
"SaveEvent.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/lhaaso.ihep.ac.cn/anysw/slc5_ia64_gcc73/external/root/6.14.00/include",
"/workfs/ybj/youzhiyong/WCDA/WCDARec/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "EventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$WCDAMcEvent.h")))  WCDAMcEvent;
class __attribute__((annotate("$clingAutoload$LHCALEvent.h")))  __attribute__((annotate("$clingAutoload$WCDAMcRecEvent.h")))  LHCALEvent;
class __attribute__((annotate("$clingAutoload$WCDAMcRecEvent.h")))  WCDAMcRecEvent;
class __attribute__((annotate("$clingAutoload$SaveEvent.h")))  SaveEvent;
class __attribute__((annotate("$clingAutoload$SaveEvent.h")))  SaveWCDAHit;
class __attribute__((annotate("$clingAutoload$SaveEvent.h")))  SavePLUSHit;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "EventDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "WCDAMcEvent.h"
#include "WCDAMcRecEvent.h"
#include "LHCALEvent.h"
#include "SaveEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LHCALEvent", payloadCode, "@",
"SaveEvent", payloadCode, "@",
"SavePLUSHit", payloadCode, "@",
"SaveWCDAHit", payloadCode, "@",
"StrDup", payloadCode, "@",
"WCDAMcEvent", payloadCode, "@",
"WCDAMcRecEvent", payloadCode, "@",
"operator+", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("EventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_EventDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_EventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_EventDict() {
  TriggerDictionaryInitialization_EventDict_Impl();
}
