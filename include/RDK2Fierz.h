#ifndef INCLUDE_RDK2FIERZ_H_
#define INCLUDE_RDK2FIERZ_H_

#include "TString.h"

#include "MRHistDim.h"

#include "RDK2CutSet.h"

void createHomogenousEventSettings(int numFiles,TString eventSettingsString);
void createHomogenousEventsUsingSettings(int numFiles, int numEvents, TString eventSettingsString);
void createLilBTreesFromHmgEventFiles(int numFiles, TString eventSettingsString, int numLilB, double startLilB, double incrementLilB);
void makeDWCutHists(int numDWCut, TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet);
void plotLilBExpFitToMC(TString expHistFilePath,TString mcAnalIDString, HistDim histDim,int numDWCut,double startLilB, double incrementLilB, double beginEnergy, double endEnergy);

#endif /* INCLUDE_RDK2FIERZ_H_ */
