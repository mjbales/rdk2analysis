#ifndef RDK2ANALYSISPLOTTING_H_INCLUDED
#define RDK2ANALYSISPLOTTING_H_INCLUDED

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TPad.h"
#include "TPaveText.h"

#include "RDK2MCAnalyzer.h"
#include "RDK2ExpAnalyzer.h"
#include "MRMisc.h"

#include <vector>

const int GRAPHICS_NUM_TYPES=12;
enum PlotVarType{PLOTVAR_PT,PLOTVAR_EE,PLOTVAR_PE,PLOTVAR_GE,PLOTVAR_GEVAR};
enum PlotType{PLOT_COMP,PLOT_RESID,PLOT_NORMRESID,PLOT_DETS};
const int colorList[GRAPHICS_NUM_TYPES]= { kBlue, kRed, kGreen, kCyan,kBlack,kMagenta, kPink, kOrange, kViolet,kAzure,kYellow,kGray};
const int markerList[GRAPHICS_NUM_TYPES]={20,22,23,25,24,32,33,3,5,28,26,30};

class RDK2AnalysisPlotter
{
    public:
    RDK2AnalysisPlotter();
    RDK2AnalysisPlotter(TString inpTitleString,int inpNumExp, RDK2ExpAnalyzer** inpExp,int inpNumMC, RDK2MCAnalyzer** inpMC);
    RDK2AnalysisPlotter(TString inpTitleString, TString expIDList,TString mcIDList);
    ~RDK2AnalysisPlotter();
    void clearPlotHists();
    void makeAll();
    void makeEPPlot(CoDet detType);
    void makeEPGPlot(CoDet detType);
    TH1* getExpHist(int i,CoDet detType, PlotVarType plotVarType,int det=0);
    TH1* getMCHist(int i,CoDet detType, PlotVarType plotVarType,int det=0);
//    void drawPlot(CoDet detType, PlotVarType plotVarType,PlotType plotType);


    protected:
    TString titleString;
    int numExp;
    int numMC;
    std::vector<TH1*> plotHists;

    RDK2ExpAnalyzer** exp;
    RDK2MCAnalyzer** mc;
    RDK2CutSet extremeCutSet;
    bool deleteAllAnalysis;

    TPaveText* makeTitleBox(CoDet detType);
    void drawHist(int numExp, TH1** expHists,int numMC, TH1** mcHists);
    void drawPlot(CoDet detType, PlotVarType plotVarType,PlotType plotType);
    void drawPlotComp(CoDet detType, PlotVarType plotVarType);
    void drawPlotResid(CoDet detType, PlotVarType plotVarType,bool normalize);
    void drawPlotDets(CoDet detType, PlotVarType plotVarType);

    void setDisplayQuantitiesForHist(int i,double factor,TH1* inpHist, CoDet detType,PlotType plotType, PlotVarType plotVarType);
    int getTitleBoxLines(CoDet detType);
//    void determineExtremeCutSet();
};

void createExpAnalysisAndPlot(TString inpName, TString inpTitle, TString inpPIDString,TString inpGroupString, RDK2CutSet inpCutSet,const double inpBGOWindows[6],const double inpBAPDWindows[6]);
void createMCAnalysisAndPlot(TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4,TString p4RID, TString e4RID,TString g4RID, RDK2CutSet inpCutSet);
void createAndPlotPESections(TString nameString, RDK2Set* inpThreeBodySet, int numSections, double begin, double end,int pOption=0, bool useExp=false);
#endif // RDK2ANALYSISPLOTTING_H_INCLUDED
