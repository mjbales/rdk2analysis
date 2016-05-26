#ifndef INCLUDE_RDK2IO_H_
#define INCLUDE_RDK2IO_H_

#include "TString.h"
#include "TChain.h"
#include "RDK2Constants.h"

TString getNextEventFileName(TString inpFileName);                       //Based on Changbo's event naming pattern, increments the set number
TString getNextResultFileName(TString inpFileName);                      //Based on Matt's result file naming pattern, increments the set number
TChain* makeResultChain(TString firstFileName, TString treeName,unsigned int numEventFiles);
TChain* makeEventChain(TString firstFileName,TString treeName,unsigned int numEventFiles, bool isEventGenerator=false);
TChain* makeLilBChain(TString eventString,TString treeName,unsigned int numEventFiles);
TChain* makeExpChain(TString PIDString,const int numFiles,const int* seriesIDArray);
void makeGeantMacroFilesForClusterEventGenerator(int startSet, int numSets, TString outResultTypeString, ParType partype, TString eventString, TString addLine = "");
void makeGeantMacroFilesForClusterEventGeneratorMakeEventsOnly(int startSet,int numSets, TString eventString,TString addLine="");


#endif /* INCLUDE_RDK2IO_H_ */
