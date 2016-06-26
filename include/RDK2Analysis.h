#ifndef RDK2ANALYSIS_H_INCLUDED
#define RDK2ANALYSIS_H_INCLUDED

#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TF1.h"
#include "TChain.h"
#include "TGraph.h"
#include "TRandom3.h"

#include "RDK2Set.h"
#include "RDK2Constants.h"




const double KevinGausConstx1[12]={1.0117587419,1.2969635549,1.0582115564,1.0866376782,1.0219480886,1.4117148459,1.0931917292,0.8982337281,1.1857472628,1.1184266421,1.2077854238,0.7276814981};
const double KevinGausConstx2[12]={-211.2134917937,-64.9081602237,-126.6368362482,-47.7886212268,-19.0498820874,-340.2122513963,-76.0295788977,-64.9416983908,70.1920045111,-51.6530784506,-16.934915825,36.8215632563};
const double KevinGausConstx3[12]={0.9465278357,0.9379644556,0.9451269649,0.9396714702,0.9407854443,0.9297717251,0.9402597255,0.9448462323,0.9389171743,0.9379893338,0.9357536715,0.9455026732};
const double KevinGausConstx4[12]={33.7345338594,39.7115920595,34.5275142523,35.8415364719,35.6507843239,44.2876380116,37.7001367425,35.4813420564,39.79751993,40.2861988094,41.025304126,34.3575402677};
const double KevinGausCal_Fac[12]={1.0176856608,1.0193338075,1.0181775132,1.0222475976,1.0213198384,1.0223011845,1.0193857171,1.0173191223,1.0182445502,1.0186146803,1.0200374549,1.0179941233};

const double KevinGammaConstx1[12]={1.0015113229,1.2569448129,1.002297337,1.069226213,0.9889569799,1.36390385,0.8613405652,1.1565009698,1.0566707343,1.1476167503,0.7215809135};
const double KevinGammaConstx2[12]={-205.0752503292,-38.3214217695,-86.3300869419,-31.1777992687,9.6754677238,-305.9085976251,-38.7312967521,94.7997257994,-6.3567448688,25.3275788728,40.6734049323};
const double KevinGammaConstx3[12]={0.9501680633,0.9413172521,0.9485032541,0.9434849403,0.944389067,0.9337603817,0.9481730684,0.9418624337,0.9412491603,0.9390578475,0.9484160588};
const double KevinGammaConstx4[12]={31.51360719,37.7891437158,32.5324530242,33.5993431603,33.5447453048,41.8834617358,33.4848795456,38.1690302086,38.3895651946,39.1195148492,32.629600562};
const double KevinGammaCal_Fac[12]={1.0165913258,1.0165913258,1.0165913258,1.0165913258,1.0165913258,1.0165913258,1.0165913258,1.0165913258,1.0165913258,1.0165913258,1.0165913258};

const int MattKhodyukParNumRanges=8;
const double MattKhodyukParLow[MattKhodyukParNumRanges]={0,13,17,18,90,91,95,100};
const double MattKhodyukParHigh[MattKhodyukParNumRanges]={13,17,18,90,91,95,100,800};
const double MattKhodyukParP0[MattKhodyukParNumRanges]={0.354676,4.62433,0.304524,0.453489,0.754255,0.41428,0.412629,0.78667};
const double MattKhodyukParP1[MattKhodyukParNumRanges]={0.020557,-0.26468,0.0134874,0.0166349,0.00602818,0.00341621,0.00321131,0.00113219};
const double MattKhodyukParP2[MattKhodyukParNumRanges]={0.00103943,-0.0168161,0.000546591,-0.000152853,4.1124e-05,2.55743e-05,2.26598e-05,-1.55463e-06};
const double MattKhodyukParP3[MattKhodyukParNumRanges]={3.948e-05,7.87402e-05,1.8438e-05,-1.04802e-06,1.72795e-07,1.55279e-07,1.29642e-07,-8.30498e-10};
const double MattKhodyukParP4[MattKhodyukParNumRanges]={2.61034e-07,7.28475e-05,3.25152e-07,1.15096e-08,-1.20206e-09,3.87577e-10,3.08751e-10,1.71459e-12};
const double MattKhodyukParP5[MattKhodyukParNumRanges]={-1.50273e-07,4.13368e-06,-2.27699e-08,1.6232e-10,-4.76543e-11,-9.53956e-12,-6.90997e-12,2.31777e-15};
const double MattKhodyukParP6[MattKhodyukParNumRanges]={-2.02813e-08,-2.86927e-07,-3.6416e-09,-1.4141e-12,-9.06227e-13,-2.4887e-13,-1.697e-13,-2.81793e-18};

enum EDepModel {NO_EDEP_MODEL,KEVIN_GAUS_EDEP,MATT_EXP_EDEP,KEVIN_GAUS_EDEP_NOEP,KE_FIRST_MODEL,KE_FIRST_NOEP_MODEL,KHODYUK_EDEP_MODEL,MATT_EDGE_EDEP_MODEL,
    KHODYUK_EDEP_MODEL_DOUBLE,KHODYUK_EDEP_MODEL_HALF,KEVIN_GAMMA_EDEP,KHODYUK_EDEP_MODEL_WGAMMA,KHODYUK_PARAM_MODEL_WGAUSS,KHODYUK_PARAM_MODEL_WGAMMA,
    KHODYUK_PARAM_MODEL_NORES,KHODYUK_PARAM_MODEL_DBLGAUSS,KHODYUK_PARAM_MODEL_DBLRES, EDEP_WRESONLY,ALLEXP_PARAM_MODEL,ALLEXP_PARAM_MODEL_HIGH,ALLEXP_PARAM_MODEL_LOW};


double getMLS(TH1I* unscaledHist, TH1I* scaledHist,int begin, int end, int loopNum,double& bestScale);
double getChiSquared(TH1I* unscaledHist, TH1I* scaledHist,int begin, int end,double scale);
void getResiduals(TH1I* dataHist, TH1I* scaledMCHist, double* res,int numBins);
void getStandardizedResiduals(TH1I* dataHist, TH1I* scaledMCHist, double* res,int numBins);

void calcResultsStats(TString firstFileName,int numFiles,bool homogeneous,bool fourBody);
void combineGammaEnergyHists();
TCanvas* makeRadRangeComparison();
void makeEG0Combined();
void createEGOpeningAngleGraph(int type);
TCanvas* plotManyHistogramsTogether(int numHists,TH1** theHists, TLegend* theLegend,TString inpTitle="",TString canvasName="c2");
void makeHistWithGeantandMRK(TString mrkFileName,TString geantProtonFileName,TString geantElectronFileName,TString littlebFileName, TString eventFileName,unsigned int numFiles);
void makeMRKAndGeantHistograms(TString axis,ParType parType);
void checkGeantRKAccuracyToStandard(TString standardFileName,TString testFileName,int numEvents,double accuracyReq);
void makeRadHighHist();
TCanvas* plotBNLSpectra500(int num);
TCanvas* plotBNLSpectra1500(int num);
TCanvas* plotBNLSpectra750(int num);
TCanvas* plotBNLSpectra1400(int num);
TCanvas* plotBNLSpectra600(int num);
TCanvas* plotBNLSpectra1066(int num);
TCanvas* plotBNLSpectra400(int num);
TCanvas* plotee0VsExp();
TCanvas* ploteeFVsExp();
//void makeSIO2Tree(TString gammaFileName);
//void makeSIO2Trees(TString firstFileName,int numFiles);
void makeGammaEnergyHistsBGO(TString protonString, TString electronString, TString gammaString, int numEPGFiles,TString gammaLevelString,double numEPExpEvents,double epPer3BodyEvent,TString tofNameString="1.4us",const int numBGOs=12,EDepModel eDepModel=KEVIN_GAUS_EDEP);
void makeBAPDCEGraphs();
void makeBAPDProgression(TString radString,bool logGraph=false);
void makeCombinedBGOGraphs(TString radhighstring, bool logGraph);
void compareCalibrationMCs(TString isotopeString,TString resultsFileName,bool useLog,bool dontRecalcMC);
void checkFalseTriple(TString detType,const int numDets);
void combineBGORanges(TString radlowString,TString radmidString,TString radhighString,TString analysisNameString,int multiplicity,int numBGOSUsed=12);
void compareBGOSingleCrystal(TString resultsFileName,bool useLog);
void compareProtonSpectrum(TString geantResultsFile, TString inpStrimString,bool useLog);
void compareTestEm5(TString particle,TString material,int ver1,int ver2,TString histName,bool useLog);
void compareAllHistsTestEm5(TString particle,TString material,int ver1,int ver2);
TH1D* makeEPProtonTOF(int numFiles, TString protonString, TString electronString,TString additionalCuts="",int firstFile=0, bool useOldGeantData=false);
TH1D* makeEPProtonTOFOldMRK(int numFiles,int firstFile=0,TString datecode="100421");
void makeEPGToFComparison(TString expToFHistFileName, TString midString,TString highString,TString tofString,int multiplicity=1);
TH1D* makeEPHist(int numFiles, TString protonString, TString electronString,TString drawString,TString cutsString,int numBins,double bottomHist,double topHist);

double applyEDepModel(double e_dep,int detectorNumber, TRandom3* rootRanGen,EDepModel eDepModel = KEVIN_GAUS_EDEP,TGraph* responseGraph = nullptr, int modelNumber=3); //Presumes det number is 1-12
void makeMultiplicityFile(TString firstGammaFile, int numEPGFiles); //includes Kevin's Energy Resolution
void analyzeBGORDK2(TString expToFHistFileName,TString expBGOSpectraHistFileName,TString tofNameString,int multiplicity);
void makeSingleBGOSpectrum(TString filePath,TString histName, bool useKevinResolution=false);
void makeBGOSim3DHists(int numEPGFiles,TString protonMidString, TString electronMidString, TString gammaMidString,TString protonHighString, TString electronHighString, TString gammaHighString,double epPer3BodyEvent,EDepModel eDepModel=KEVIN_GAUS_EDEP);
void makeBGOSim3DHists(RDK2Set* hmgSet,RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet);
void makeMCTextHistogramsFrom3D(TString analysisIDString,TString processingRunString, TString mcRIDString, AnalysisCuts theCuts, int multiplicity);
void makeBGOSingleEnergyGraphs();
void examinePLO();
void examinePeaksFromRDKPositions(TString fileName);
void plotQuadraticVsExpResponse();
double getBGOHighToLowRatio(TH1* theHist,double& ratioError); // presumes 10 - 800 or 1000 keV range with 5 keV bins  Low is 10-200 keV high is 200-800 keV
void getBGOHighToLowRatioList(int numHists, TString* histFileList);
TH1* combineAndScaleRadHists(TString newHistName,int numFiles, TH1* radLowHist, TH1* radMidHist, TH1* radHighHist, double epPer3BodyEvent);
TH1* combineAndScaleRadSetsToHist(RDK2Set* hmgSet,RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet,TString drawString,TCut epCut,TCut gCut,HistDim inpDimX,HistDim inpDimY,HistDim inpDimZ);
TH1* combineAndScaleRadSetsToHist(RDK2Set* hmgSet,RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet,TString drawString,TCut epCut,TCut gCut,HistDim inpDimX,HistDim inpDimY);
TH1* combineAndScaleRadSetsToHist(RDK2Set* hmgSet,RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet,TString drawString,TCut epCut,TCut gCut,HistDim inpDimX);
TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString drawString,TCut epCut,TCut gCut,HistDim inpDimX);
TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString drawString,TCut epCut,TCut gCut,HistDim inpDimX,HistDim inpDimY);
TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString drawString,TCut epCut,TCut gCut,HistDim inpDimX,HistDim inpDimY,HistDim inpDimZ);
TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString drawString,TCut epCut,TCut gCut,TH1* histPrototype);
//TGraph* getKhodyukBGOResponseGraph(); //No longer use
TGraph* getKhodyukBGOResponseGraph2();
TGraph* getKhodyukBGOResponseGraphOfficial();
TGraph* getMattBGOResponseGraph();
void plotKhodyukBGOResponseGraph();
TH1D* makeBAPDCombinedHist(TString histName,RDK2Set* hmgSet,RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet,HistDim inpDim,TCut inpEPCut,bool doKEFirst=true);
TH1D* createHistWBGOResolution(TH1D* inpHist, int detectorNumber,TRandom3* rootRanGen,EDepModel eDepModel); //Currently only applying EDepModel
void makeAllMCHistograms(TString dateString,RDK2Set* hmgSet,RDK2Set* lowSet,RDK2Set* midSet,RDK2Set* highSet);
void makeAllMCHistograms(TString dateString,RDK2Set* threeBodySet,RDK2Set* fourBodySet);
void makeAllStandardComparisonGraphs(TString mcDateString,RDK2Set* hmgSet,RDK2Set* lowSet,RDK2Set* midSet,RDK2Set* highSet,TString expDateString,TString expPIDString,TString expGroupIDString);
void makeAllStandardComparisonGraphs(TString mcDateString,RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString expDateString,TString expPIDString,TString expGroupIDString);
TGraph* makeKevinResolutionGraph(int detectorNumber);
double evalMattKhodyukParam(double energyIn);
double evalMattAllExpParamModel(double energyIn, EDepModel inpEDepModel);//Model 1
double evalMattAllExpParamModel2(double energyIn, EDepModel inpEDepModel);//Model 2
double evalMattAllExpParamModel3(double energyIn, EDepModel eDepModel);//Model 3
double evalMattAllExpParamModel4(double energyIn, EDepModel eDepModel);//Model 4-testing only <16.5
double evalMattAllExpParamModel5(double energyIn, EDepModel eDepModel);//Model 5-testing only <13.4
void makeIndividualDetectorHistograms(TString mcAnalysisDateString,RDK2Set* threeBodySet, RDK2Set* fourBodySet);
void createHistogramsForMCCalibraitonSourceData(int numFiles,TString runString);
void createDerivedResultsFile(int startingEvent, TTree* inpResultsTree,TString filePathOut, int numberEvents=1000000);
void createDerivedResultsFilesForCal(int numFiles,TString runString);
TH1D* createCalMCHistWithRes(int detectorNum, TString mcRunString, TString isotopeString,bool prebakeKhodyuk, double isotopeRatio);
int getMCGaussPeaksToFit(TString seriesName, bool prebakeKhodyuk, TF1** gaussPeakOut);
double deconvolutePrebakeKhodyuk(int numMCFiles, TString mcRunString, int BGONum, TF1* gaussPeakInOut);
Double_t fitFuncForMCHist(Double_t* v, Double_t* par);
TF1* fitMCtoExp(TH1* expHist, TH1* mcHist,int binMin,int binMax);
double getRandom();
TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString quantityToDraw,TCut epCut,TCut inpGCut,HistDim inpDimX);
TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString quantityToDraw,TCut epCut,TCut inpGCut,HistDim inpDimX,HistDim inpDimY);
TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString quantityToDraw,TCut epCut,TCut inpGCut,HistDim inpDimX,HistDim inpDimY,HistDim inpDimZ);
TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet,RDK2Set* fourBodySet,TString quantityToDraw,TCut epCut,TCut inpGCut,TH1* histToFill);
TH1D* createEG0HitPercentHist(RDK2Set*threeBodySet,RDK2Set* fourBodySet,TCut epCut,const HistDim inpDimX);
TH1D* expSpecToTheorySpec(TH1D* expHist, RDK2Set* threeBodySet, RDK2Set* fourBodySet,TCut epCut,TString histPath2D, const HistDim inpDimX,const HistDim inpDimY);
TH1D* expSpecToTheorySpecWError(TH1D* expHist, RDK2Set* threeBodySet, RDK2Set* fourBodySet,TCut epCut,TString histPath2D);
void createPlotOfRIDSequence(int numFiles,int beginningRID,TString dateString="");
double BAPDCollectionEfficiencyModel(double inpEnergy, double depth);
void makeCombinedBGOAndBAPDPlot(TString expBGOPath, double bgoCorrection, double epBGOCount, TString expBAPDPath, double bapdCorrection, double epBAPDCount, int mcAIDNum,bool logx=true, bool logy=true, bool absScale=false, bool varWidth=true, bool usePDF=true);
void makeEGTimingPlot(TString histPath,TString bhistPath);
void makeZ0Hist();
void makeOtherParameterPlot(TString expID, double correction, int mcAIDNum, TString coTypeStr, TString plotType, TString imageType="pdf", bool halfPlotNoResids=false);
void makeOtherParameterPlots(TString bgoExpID, double bgoCorrection, TString apdExpID, double apdCorrection, int mcID, bool halfPlotNoResids=false);
TH1D* getExpEEBSHist();
void outputParameterizedLOModelToText();
void outputBGODataToFile();
#endif // RDK2ANALYSIS_H_INCLUDED
