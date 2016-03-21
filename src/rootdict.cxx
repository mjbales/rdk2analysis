// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME rootdict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
#include "RDK2CutSet.h"
#include "RDK2Set.h"
#include "RDK2MCAnalyzer.h"
#include "RDK2ExpAnalyzer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RDK2CutSet(void *p = 0);
   static void *newArray_RDK2CutSet(Long_t size, void *p);
   static void delete_RDK2CutSet(void *p);
   static void deleteArray_RDK2CutSet(void *p);
   static void destruct_RDK2CutSet(void *p);
   static void streamer_RDK2CutSet(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RDK2CutSet*)
   {
      ::RDK2CutSet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RDK2CutSet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RDK2CutSet", ::RDK2CutSet::Class_Version(), "RDK2CutSet.h", 10,
                  typeid(::RDK2CutSet), DefineBehavior(ptr, ptr),
                  &::RDK2CutSet::Dictionary, isa_proxy, 16,
                  sizeof(::RDK2CutSet) );
      instance.SetNew(&new_RDK2CutSet);
      instance.SetNewArray(&newArray_RDK2CutSet);
      instance.SetDelete(&delete_RDK2CutSet);
      instance.SetDeleteArray(&deleteArray_RDK2CutSet);
      instance.SetDestructor(&destruct_RDK2CutSet);
      instance.SetStreamerFunc(&streamer_RDK2CutSet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RDK2CutSet*)
   {
      return GenerateInitInstanceLocal((::RDK2CutSet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RDK2CutSet*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RDK2Set(void *p = 0);
   static void *newArray_RDK2Set(Long_t size, void *p);
   static void delete_RDK2Set(void *p);
   static void deleteArray_RDK2Set(void *p);
   static void destruct_RDK2Set(void *p);
   static void streamer_RDK2Set(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RDK2Set*)
   {
      ::RDK2Set *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RDK2Set >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RDK2Set", ::RDK2Set::Class_Version(), "RDK2Set.h", 41,
                  typeid(::RDK2Set), DefineBehavior(ptr, ptr),
                  &::RDK2Set::Dictionary, isa_proxy, 16,
                  sizeof(::RDK2Set) );
      instance.SetNew(&new_RDK2Set);
      instance.SetNewArray(&newArray_RDK2Set);
      instance.SetDelete(&delete_RDK2Set);
      instance.SetDeleteArray(&deleteArray_RDK2Set);
      instance.SetDestructor(&destruct_RDK2Set);
      instance.SetStreamerFunc(&streamer_RDK2Set);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RDK2Set*)
   {
      return GenerateInitInstanceLocal((::RDK2Set*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RDK2Set*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RDK2MCAnalysis(void *p = 0);
   static void *newArray_RDK2MCAnalysis(Long_t size, void *p);
   static void delete_RDK2MCAnalysis(void *p);
   static void deleteArray_RDK2MCAnalysis(void *p);
   static void destruct_RDK2MCAnalysis(void *p);
   static void streamer_RDK2MCAnalysis(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RDK2MCAnalysis*)
   {
      ::RDK2MCAnalysis *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RDK2MCAnalysis >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RDK2MCAnalysis", ::RDK2MCAnalysis::Class_Version(), "RDK2MCAnalyzer.h", 17,
                  typeid(::RDK2MCAnalysis), DefineBehavior(ptr, ptr),
                  &::RDK2MCAnalysis::Dictionary, isa_proxy, 16,
                  sizeof(::RDK2MCAnalysis) );
      instance.SetNew(&new_RDK2MCAnalysis);
      instance.SetNewArray(&newArray_RDK2MCAnalysis);
      instance.SetDelete(&delete_RDK2MCAnalysis);
      instance.SetDeleteArray(&deleteArray_RDK2MCAnalysis);
      instance.SetDestructor(&destruct_RDK2MCAnalysis);
      instance.SetStreamerFunc(&streamer_RDK2MCAnalysis);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RDK2MCAnalysis*)
   {
      return GenerateInitInstanceLocal((::RDK2MCAnalysis*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RDK2MCAnalysis*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RDK2ExpAnalysis(void *p = 0);
   static void *newArray_RDK2ExpAnalysis(Long_t size, void *p);
   static void delete_RDK2ExpAnalysis(void *p);
   static void deleteArray_RDK2ExpAnalysis(void *p);
   static void destruct_RDK2ExpAnalysis(void *p);
   static void streamer_RDK2ExpAnalysis(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RDK2ExpAnalysis*)
   {
      ::RDK2ExpAnalysis *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RDK2ExpAnalysis >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RDK2ExpAnalysis", ::RDK2ExpAnalysis::Class_Version(), "RDK2ExpAnalyzer.h", 20,
                  typeid(::RDK2ExpAnalysis), DefineBehavior(ptr, ptr),
                  &::RDK2ExpAnalysis::Dictionary, isa_proxy, 16,
                  sizeof(::RDK2ExpAnalysis) );
      instance.SetNew(&new_RDK2ExpAnalysis);
      instance.SetNewArray(&newArray_RDK2ExpAnalysis);
      instance.SetDelete(&delete_RDK2ExpAnalysis);
      instance.SetDeleteArray(&deleteArray_RDK2ExpAnalysis);
      instance.SetDestructor(&destruct_RDK2ExpAnalysis);
      instance.SetStreamerFunc(&streamer_RDK2ExpAnalysis);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RDK2ExpAnalysis*)
   {
      return GenerateInitInstanceLocal((::RDK2ExpAnalysis*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RDK2ExpAnalysis*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RDK2CutSet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RDK2CutSet::Class_Name()
{
   return "RDK2CutSet";
}

//______________________________________________________________________________
const char *RDK2CutSet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2CutSet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RDK2CutSet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2CutSet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RDK2CutSet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2CutSet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RDK2CutSet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2CutSet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RDK2Set::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RDK2Set::Class_Name()
{
   return "RDK2Set";
}

//______________________________________________________________________________
const char *RDK2Set::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2Set*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RDK2Set::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2Set*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RDK2Set::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2Set*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RDK2Set::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2Set*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RDK2MCAnalysis::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RDK2MCAnalysis::Class_Name()
{
   return "RDK2MCAnalysis";
}

//______________________________________________________________________________
const char *RDK2MCAnalysis::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2MCAnalysis*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RDK2MCAnalysis::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2MCAnalysis*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RDK2MCAnalysis::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2MCAnalysis*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RDK2MCAnalysis::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2MCAnalysis*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RDK2ExpAnalysis::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RDK2ExpAnalysis::Class_Name()
{
   return "RDK2ExpAnalysis";
}

//______________________________________________________________________________
const char *RDK2ExpAnalysis::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2ExpAnalysis*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RDK2ExpAnalysis::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RDK2ExpAnalysis*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RDK2ExpAnalysis::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2ExpAnalysis*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RDK2ExpAnalysis::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RDK2ExpAnalysis*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RDK2CutSet::Streamer(TBuffer &R__b)
{
   // Stream an object of class RDK2CutSet.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> eELow;
      R__b >> eEHigh;
      R__b >> pELow;
      R__b >> pEHigh;
      R__b >> pTLow;
      R__b >> pTHigh;
      R__b >> gELow;
      R__b >> gEHigh;
      R__b >> bGELow;
      R__b >> bGEHigh;
      R__b.CheckByteCount(R__s, R__c, RDK2CutSet::IsA());
   } else {
      R__c = R__b.WriteVersion(RDK2CutSet::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << eELow;
      R__b << eEHigh;
      R__b << pELow;
      R__b << pEHigh;
      R__b << pTLow;
      R__b << pTHigh;
      R__b << gELow;
      R__b << gEHigh;
      R__b << bGELow;
      R__b << bGEHigh;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RDK2CutSet(void *p) {
      return  p ? new(p) ::RDK2CutSet : new ::RDK2CutSet;
   }
   static void *newArray_RDK2CutSet(Long_t nElements, void *p) {
      return p ? new(p) ::RDK2CutSet[nElements] : new ::RDK2CutSet[nElements];
   }
   // Wrapper around operator delete
   static void delete_RDK2CutSet(void *p) {
      delete ((::RDK2CutSet*)p);
   }
   static void deleteArray_RDK2CutSet(void *p) {
      delete [] ((::RDK2CutSet*)p);
   }
   static void destruct_RDK2CutSet(void *p) {
      typedef ::RDK2CutSet current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RDK2CutSet(TBuffer &buf, void *obj) {
      ((::RDK2CutSet*)obj)->::RDK2CutSet::Streamer(buf);
   }
} // end of namespace ROOT for class ::RDK2CutSet

//______________________________________________________________________________
void RDK2Set::Streamer(TBuffer &R__b)
{
   // Stream an object of class RDK2Set.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> numFiles;
      eventsDateString.Streamer(R__b);
      eventsTypeString.Streamer(R__b);
      pResultsID.Streamer(R__b);
      eResultsID.Streamer(R__b);
      gResultsID.Streamer(R__b);
      savedEPCut.Streamer(R__b);
      R__b >> filesLoaded;
      R__b >> isHmg;
      R__b >> isRad;
      R__b >> isEventGenerator;
      R__b >> epRate;
      R__b.CheckByteCount(R__s, R__c, RDK2Set::IsA());
   } else {
      R__c = R__b.WriteVersion(RDK2Set::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << numFiles;
      eventsDateString.Streamer(R__b);
      eventsTypeString.Streamer(R__b);
      pResultsID.Streamer(R__b);
      eResultsID.Streamer(R__b);
      gResultsID.Streamer(R__b);
      savedEPCut.Streamer(R__b);
      R__b << filesLoaded;
      R__b << isHmg;
      R__b << isRad;
      R__b << isEventGenerator;
      R__b << epRate;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RDK2Set(void *p) {
      return  p ? new(p) ::RDK2Set : new ::RDK2Set;
   }
   static void *newArray_RDK2Set(Long_t nElements, void *p) {
      return p ? new(p) ::RDK2Set[nElements] : new ::RDK2Set[nElements];
   }
   // Wrapper around operator delete
   static void delete_RDK2Set(void *p) {
      delete ((::RDK2Set*)p);
   }
   static void deleteArray_RDK2Set(void *p) {
      delete [] ((::RDK2Set*)p);
   }
   static void destruct_RDK2Set(void *p) {
      typedef ::RDK2Set current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RDK2Set(TBuffer &buf, void *obj) {
      ((::RDK2Set*)obj)->::RDK2Set::Streamer(buf);
   }
} // end of namespace ROOT for class ::RDK2Set

//______________________________________________________________________________
void RDK2MCAnalysis::Streamer(TBuffer &R__b)
{
   // Stream an object of class RDK2MCAnalysis.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      threeBodySet.Streamer(R__b);
      fourBodySet.Streamer(R__b);
      theCutSet.Streamer(R__b);
      R__b >> epPer3Decay;
      R__b >> epPer3DecayError;
      R__b >> epgPer4Decay;
      R__b >> epgPer4DecayError;
      R__b >> epbgPer4Decay;
      R__b >> epbgPer4DecayError;
      R__b >> epgPerEP;
      R__b >> epgPerEPError;
      R__b >> epbgPerEP;
      R__b >> epbgPerEPError;
      R__b >> branchingRatio;
      R__b >> ep_eEHist;
      R__b >> ep_pEHist;
      R__b >> ep_pTHist;
      R__b >> epg_eEHist;
      R__b >> epg_pEHist;
      R__b >> epg_pTHist;
      R__b >> epbg_eEHist;
      R__b >> epbg_pEHist;
      R__b >> epbg_pTHist;
      R__b >> epg_gEAvgHist;
      R__b >> epg_gEAvgVarHist;
      int R__i;
      for (R__i = 0; R__i < 12; R__i++)
         R__b >> epg_gEDetHists[R__i];
      for (R__i = 0; R__i < 12; R__i++)
         R__b >> epg_gEDetVarHists[R__i];
      R__b >> epbg_bGEAvgHist;
      R__b >> epbg_bGEAvgVarHist;
      for (R__i = 0; R__i < 3; R__i++)
         R__b >> epbg_bGEDetHists[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b >> epbg_bGEDetVarHists[R__i];
      R__b.CheckByteCount(R__s, R__c, RDK2MCAnalysis::IsA());
   } else {
      R__c = R__b.WriteVersion(RDK2MCAnalysis::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      threeBodySet.Streamer(R__b);
      fourBodySet.Streamer(R__b);
      theCutSet.Streamer(R__b);
      R__b << epPer3Decay;
      R__b << epPer3DecayError;
      R__b << epgPer4Decay;
      R__b << epgPer4DecayError;
      R__b << epbgPer4Decay;
      R__b << epbgPer4DecayError;
      R__b << epgPerEP;
      R__b << epgPerEPError;
      R__b << epbgPerEP;
      R__b << epbgPerEPError;
      R__b << branchingRatio;
      R__b << ep_eEHist;
      R__b << ep_pEHist;
      R__b << ep_pTHist;
      R__b << epg_eEHist;
      R__b << epg_pEHist;
      R__b << epg_pTHist;
      R__b << epbg_eEHist;
      R__b << epbg_pEHist;
      R__b << epbg_pTHist;
      R__b << epg_gEAvgHist;
      R__b << epg_gEAvgVarHist;
      int R__i;
      for (R__i = 0; R__i < 12; R__i++)
         R__b << epg_gEDetHists[R__i];
      for (R__i = 0; R__i < 12; R__i++)
         R__b << epg_gEDetVarHists[R__i];
      R__b << epbg_bGEAvgHist;
      R__b << epbg_bGEAvgVarHist;
      for (R__i = 0; R__i < 3; R__i++)
         R__b << epbg_bGEDetHists[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b << epbg_bGEDetVarHists[R__i];
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RDK2MCAnalysis(void *p) {
      return  p ? new(p) ::RDK2MCAnalysis : new ::RDK2MCAnalysis;
   }
   static void *newArray_RDK2MCAnalysis(Long_t nElements, void *p) {
      return p ? new(p) ::RDK2MCAnalysis[nElements] : new ::RDK2MCAnalysis[nElements];
   }
   // Wrapper around operator delete
   static void delete_RDK2MCAnalysis(void *p) {
      delete ((::RDK2MCAnalysis*)p);
   }
   static void deleteArray_RDK2MCAnalysis(void *p) {
      delete [] ((::RDK2MCAnalysis*)p);
   }
   static void destruct_RDK2MCAnalysis(void *p) {
      typedef ::RDK2MCAnalysis current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RDK2MCAnalysis(TBuffer &buf, void *obj) {
      ((::RDK2MCAnalysis*)obj)->::RDK2MCAnalysis::Streamer(buf);
   }
} // end of namespace ROOT for class ::RDK2MCAnalysis

//______________________________________________________________________________
void RDK2ExpAnalysis::Streamer(TBuffer &R__b)
{
   // Stream an object of class RDK2ExpAnalysis.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      threeBodySet.Streamer(R__b);
      fourBodySet.Streamer(R__b);
      theCutSet.Streamer(R__b);
      R__b >> epCounts;
      R__b >> epgCounts;
      R__b >> epgCountsError;
      R__b >> epbgCounts;
      R__b >> epbgCountsError;
      R__b >> epgPerEP;
      R__b >> epgPerEPError;
      R__b >> epbgPerEP;
      R__b >> epbgPerEPError;
      R__b.ReadStaticArray((int*)epCountsForLiveBGO);
      R__b.ReadStaticArray((int*)epCountsForLiveBAPD);
      R__b >> averageEPCountsForLiveBGO;
      R__b >> averageEPCountsForLiveBAPD;
      R__b >> numLiveBGO;
      R__b >> numLiveBAPD;
      R__b >> ep_eEHist;
      R__b >> ep_pEHist;
      R__b >> ep_pTHist;
      R__b >> epg_eEHist;
      R__b >> epg_pEHist;
      R__b >> epg_pTHist;
      R__b >> epbg_eEHist;
      R__b >> epbg_pEHist;
      R__b >> epbg_pTHist;
      R__b >> epg_gEAvgHist;
      R__b >> epg_gEAvgVarHist;
      int R__i;
      for (R__i = 0; R__i < 12; R__i++)
         R__b >> epg_gEDetHists[R__i];
      for (R__i = 0; R__i < 12; R__i++)
         R__b >> epg_gEDetVarHists[R__i];
      R__b >> epbg_bGEAvgHist;
      R__b >> epbg_bGEAvgVarHist;
      for (R__i = 0; R__i < 3; R__i++)
         R__b >> epbg_bGEDetHists[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b >> epbg_bGEDetVarHists[R__i];
      PIDString.Streamer(R__b);
      groupString.Streamer(R__b);
      R__b.ReadStaticArray((double*)bgoWindows);
      R__b.ReadStaticArray((double*)bapdWindows);
      R__b >> chainLoaded;
      R__b >> expChain;
      R__b.CheckByteCount(R__s, R__c, RDK2ExpAnalysis::IsA());
   } else {
      R__c = R__b.WriteVersion(RDK2ExpAnalysis::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      threeBodySet.Streamer(R__b);
      fourBodySet.Streamer(R__b);
      theCutSet.Streamer(R__b);
      R__b << epCounts;
      R__b << epgCounts;
      R__b << epgCountsError;
      R__b << epbgCounts;
      R__b << epbgCountsError;
      R__b << epgPerEP;
      R__b << epgPerEPError;
      R__b << epbgPerEP;
      R__b << epbgPerEPError;
      R__b.WriteArray(epCountsForLiveBGO, 12);
      R__b.WriteArray(epCountsForLiveBAPD, 3);
      R__b << averageEPCountsForLiveBGO;
      R__b << averageEPCountsForLiveBAPD;
      R__b << numLiveBGO;
      R__b << numLiveBAPD;
      R__b << ep_eEHist;
      R__b << ep_pEHist;
      R__b << ep_pTHist;
      R__b << epg_eEHist;
      R__b << epg_pEHist;
      R__b << epg_pTHist;
      R__b << epbg_eEHist;
      R__b << epbg_pEHist;
      R__b << epbg_pTHist;
      R__b << epg_gEAvgHist;
      R__b << epg_gEAvgVarHist;
      int R__i;
      for (R__i = 0; R__i < 12; R__i++)
         R__b << epg_gEDetHists[R__i];
      for (R__i = 0; R__i < 12; R__i++)
         R__b << epg_gEDetVarHists[R__i];
      R__b << epbg_bGEAvgHist;
      R__b << epbg_bGEAvgVarHist;
      for (R__i = 0; R__i < 3; R__i++)
         R__b << epbg_bGEDetHists[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b << epbg_bGEDetVarHists[R__i];
      PIDString.Streamer(R__b);
      groupString.Streamer(R__b);
      R__b.WriteArray(bgoWindows, 6);
      R__b.WriteArray(bapdWindows, 6);
      R__b << chainLoaded;
      R__b << expChain;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RDK2ExpAnalysis(void *p) {
      return  p ? new(p) ::RDK2ExpAnalysis : new ::RDK2ExpAnalysis;
   }
   static void *newArray_RDK2ExpAnalysis(Long_t nElements, void *p) {
      return p ? new(p) ::RDK2ExpAnalysis[nElements] : new ::RDK2ExpAnalysis[nElements];
   }
   // Wrapper around operator delete
   static void delete_RDK2ExpAnalysis(void *p) {
      delete ((::RDK2ExpAnalysis*)p);
   }
   static void deleteArray_RDK2ExpAnalysis(void *p) {
      delete [] ((::RDK2ExpAnalysis*)p);
   }
   static void destruct_RDK2ExpAnalysis(void *p) {
      typedef ::RDK2ExpAnalysis current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RDK2ExpAnalysis(TBuffer &buf, void *obj) {
      ((::RDK2ExpAnalysis*)obj)->::RDK2ExpAnalysis::Streamer(buf);
   }
} // end of namespace ROOT for class ::RDK2ExpAnalysis

namespace {
  void TriggerDictionaryInitialization_rootdict_Impl() {
    static const char* headers[] = {
"RDK2CutSet.h",
"RDK2Set.h",
"RDK2MCAnalyzer.h",
"RDK2ExpAnalyzer.h",
0
    };
    static const char* includePaths[] = {
"/home/mjbales/work/software/root/root-6.04.14/include",
"/home/mjbales/work/code/rdk2analysis/src/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Class to analyze RDK)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Class to analyze RDK)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RDK2CutSet.h")))  RDK2CutSet;
class __attribute__((annotate(R"ATTRDUMP(Class to analyze a set of RDK events, three body or four body)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Class to analyze a set of RDK events, three body or four body)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RDK2Set.h")))  RDK2Set;
class __attribute__((annotate(R"ATTRDUMP(Class to analyze RDK)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Class to analyze RDK)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RDK2MCAnalyzer.h")))  RDK2MCAnalysis;
class __attribute__((annotate(R"ATTRDUMP(Class to analyze RDK experimental data)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Class to analyze RDK experimental data)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RDK2ExpAnalyzer.h")))  RDK2ExpAnalysis;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "RDK2CutSet.h"
#include "RDK2Set.h"
#include "RDK2MCAnalyzer.h"
#include "RDK2ExpAnalyzer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RDK2CutSet", payloadCode, "@",
"RDK2ExpAnalysis", payloadCode, "@",
"RDK2MCAnalysis", payloadCode, "@",
"RDK2Set", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("rootdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_rootdict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_rootdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_rootdict() {
  TriggerDictionaryInitialization_rootdict_Impl();
}
