#ifndef INCLUDE_RDK2FIERZ_H_
#define INCLUDE_RDK2FIERZ_H_

#include "TString.h"

#include "MRHistDim.h"

#include "RDK2CutSet.h"
#include "RDK2MCAnalyzer.h"

void createHomogenousEventSettings(int numFiles,TString eventSettingsString);
void createHomogenousEventsUsingSettings(int numFiles, int numEvents, TString eventSettingsString);
void createLilBTreesFromHmgEventFiles(int numFiles, TString eventSettingsString, int numLilB, double startLilB, double incrementLilB,int littleBSet);
void makeEE0Plots(int numDWCut, TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet);
void makeDWCutHists(int numDWCut, TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet);
void plotLilBExpFitToMC(TString plotTitleString,TString expHistFilePath, TString mcAnalIDString, HistDim histDim, int numDWCut, double startLilB, double incrementLilB, double beginEnergy, double endEnergy);
void analyzeElectronBackscatteringSBDOnly(TString resultsFilePath);
void analyzeAllEBackscatteringSBDOnly();
void makeBSPlotsSBDOnly();
void makeBSPlotsSBDOnly2();
void makeIncidentEventsFromResults(TString resultsFilePath,TString eventFilePath);
void makeAllIncidentEventsFromResults(int numFiles, TString RIDString);
void makeSecondaryMacros(int numFiles, TString RIDInString, TString RIDOutString,TString libString);
void makeDWProbHists(TString inpName, int num3, TString p3RID, TString e3RID, RDK2CutSet inpCutSet,int numLilB, double startLilB, double incrementLilB);
double fitLittleB(TString nameString, TH1* expHist, TH1** mcHists,  int numLilB, double startLilB, double endLilB, double beginEnergy, double endEnergy,bool saveHistsAndPlots);
void makeBootStrapLittleB(TString expHistFilePath, HistDim histDim, RDK2MCAnalyzer& theAnalyzer,  int numLilB, double startLilB, double endLilB, double beginEnergy, double endEnergy);
#endif /* INCLUDE_RDK2FIERZ_H_ */
