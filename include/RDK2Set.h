#ifndef RDK2SET_H
#define RDK2SET_H

#include "TNamed.h"
#include "TString.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"

#include "MRHistDim.h"

#include "RDK2Constants.h"
#include "RDK2CutSet.h"

const HistDim BLANK_DIM={0,0,0};
const HistDim BAPD_DIM={80,0.,20};
const HistDim BAPD_DIM2={400,0.,20};
const HistDim BAPD_DIM3={300,0.,30};
const HistDim BAPD_DIM4={140,0.,14};
const HistDim BGO_DIM={800,0.,800};
const HistDim BGO_DIM2={8000,0.,800};
const HistDim EDEPE_DIM={800,0.,800};
const HistDim EDEPE_DIM2={8000,0.,800};
const HistDim TOFP_DIM={625,-20e-9,625*40e-9-20e-9};
const HistDim TOFP_DIM2={621,-20e-9,621*40e-9-20e-9};
const HistDim TOFP_DIM3={750,-20e-9,750*40e-9-20e-9};
const HistDim TOFP_DIM4={620,-20e-9,620*40e-9-20e-9};
const HistDim TOFP_DIM5={625,0,620*40e-9};
const HistDim TOFPCCHAN_DIM={300,-.5,300-.5};
const HistDim BGOCAL_DIM={10000,0,1000};
const HistDim BGOEXPCAL_DIM={2000,15,30285};
const HistDim EDEPP_DIM={500,0,50};
const HistDim EDEPP_DIM2={8000,0,800};
const HistDim EGTIME_DIM={1501,-500,1000};
const HistDim EGTIME_DIM2={1800,-800.5,999.5};
const double BGO_VARHIST_LOWEDGES[21]={10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,800};
const double BAPD_VARHIST_LOWEDGES[21]={10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,800}; //This is not correct

class RDK2Set : public TNamed
{
public:
	RDK2Set();
	RDK2Set(TString inpName, TString inpTitle,int inp_numFiles,TString inp_eventsDateString,TString inp_eventsTypeString,TString inp_pResultsID,
			TString inp_eResultsID,TString inp_gResultsID);
	virtual ~RDK2Set();

	const RDK2Set& operator=(const RDK2Set& inp);

	inline int GetNumFiles() { return numFiles; }
	inline int GetNumEvents(){return eventsChain->GetEntries();}
	inline void SetNumFiles(int val) { numFiles = val; }
	inline TString GetEventDateString() { return eventsDateString; }
	inline void SetEventDateString(TString val) { eventsDateString = val; }
	inline TString GetEventTypeString() { return eventsTypeString; }
	inline void SetEventTypeString(TString val) { eventsTypeString = val; }
	inline TString GetPResultsID() { return pResultsID; }
	inline void SetPResultsID(TString val) { pResultsID = val; }
	inline TString GetEResultsID() { return eResultsID; }
	inline void SetEResultsID(TString val) { eResultsID = val; }
	inline TString GetGResultsID() { return gResultsID; }
	inline void SetGResultsID(TString val) { gResultsID = val; }
	inline void SetIsHmg(bool inp){ isHmg = inp; }
	inline void SetIsRad(bool inp){ isRad = inp; }
	inline void SetIsEventGenerator(bool inp){ isEventGenerator = inp; }

	TH1* makeHist(TString histName,TString drawString,TCut inpCut,HistDim inpHistDimX); //creates a histogram for the set based on the drawString and cuts
	TH1* makeHist(TString histName,TString drawString,TCut inpCut,HistDim inpHistDimX, HistDim inpHistDimY); //creates a histogram for the set based on the drawString and cuts
	TH1* makeHist(TString histName,TString drawString,TCut inpCut,HistDim inpHistDimX, HistDim inpHistDimY,HistDim inpHistDimZ); //creates a histogram for the set based on the drawString and cuts
	TH1* makeHist(TString histName,TString drawString,TCut inpCut,TH1* histPrototype);
	TH2D* makeHist2D(TString histName,TString drawString,TCut inpCut,HistDim inpHistDim1,HistDim inpHistDim2);
	TH2D* makeHist2D(TString histName,TString drawString,TCut inpCut,int numX, double* varArrayX, int numY, double* varArrayY);

	void makeDerivedResultsFiles(ParType parType);

	double calcEPRate(TCut epCut);
	double calcEPGRate(TCut epCut, TCut gCut);
	double calcEPGRateFor12IndividualBGOs(TCut epCut);
	double calcEPGRateFor3IndividualBAPDs(TCut epCut, double gammaEnergyCutLow = 0.5, double gammaEnergyCutHigh = 20);

	TH1* makeEEHistDecayWidthProb(TString nameString, RDK2CutSet inpCutSet,double littleb,double normConst, HistDim histDim);

	double makeAllStandardHists(TString nameString, double& error, RDK2CutSet inpCutSet, CoDet detType, double inpEPPer3Decay = 0, double branchingRatio = 0);

	void setEventList(TCut inpCut);

protected:
	void loadEventsChain();
	void loadLittleBChain();
	void loadPResultsChain();
	void loadEResultsChain();
	void loadGResultsChain();
	void loadChainsAndFriend();
	void loadDervChains();

	TChain* eventsChain;  //! EventChain

	TChain* pResultsChain; //! Proton Chain
	TChain* eResultsChain; //! Electron Chain
	TChain* gResultsChain; //! Gamma Chain
	TChain* littleBChain; //! little b Chain
	TChain* pDervChain; //! Proton Dervived Chain
	TChain* eDervChain; //! Electron Derived Chain
	TChain* gDervChain; //! Gamma Chain

	int numFiles;  //Number of files used for each chain
	TString eventsDateString; //The date TString on the events
	TString eventsTypeString; //What type of events are loaded
	TString pResultsID; //Proton results ID
	TString eResultsID;  //Electron results ID
	TString gResultsID;  //Gamma results ID

	TCut savedEPCut; //stored to double check whether EP rate should be repeated

	bool filesLoaded;  //are the files loaded into the chains
	bool isHmg; //whether it is homegenous ep events
	bool isRad; //true if 4B, false if 3B (hmg or not)
	bool isEventGenerator;  //Was this created using the event generator
	double epRate;  //What is the ep rate
public:
	ClassDef(RDK2Set,1);
	//Class to analyze a set of RDK events, three body or four body
};

TH1D* makeTH1DFromDim(TString name, TString title, HistDim inpHistDim);
TH2D* makeTH2DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim);
TH3D* makeTH3DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim, HistDim inpZHistDim);
HistDim getDimFromHist(TH1* inpHist);
int getVarArrayFromHist(TH1* inpHist, double** outArray);
#endif // RDK2SET_H
