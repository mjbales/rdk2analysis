#ifndef RDK2MCANALYSIS_H
#define RDK2MCANALYSIS_H

#define MC_NUMBGO 12
#define MC_NUMBAPD 3

//#define MC_NUMBGO 1
//#define MC_NUMBAPD 1
#include "TString.h"
#include "TH1.h"

#include "RDK2CutSet.h"
#include "RDK2Set.h"
#include "RDK2Constants.h"
#include "MRHistDim.h"

class RDK2MCAnalyzer : public TNamed
{
public:

	///Constructors/Destructors
	RDK2MCAnalyzer();
	RDK2MCAnalyzer(TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet, bool hmgEvents=false);
	RDK2MCAnalyzer(TString inpName, TString inpTitle, RDK2Set inp3Set, RDK2Set inp4Set, RDK2CutSet inpCutSet, bool hmgEvents=false);
	virtual ~RDK2MCAnalyzer();

	///Setters
	inline void SetCutSet(RDK2CutSet inp){theCutSet=inp;}
	inline void SetBranchingRatio(double inp){branchingRatio=inp;}
	inline void SetAddToAllCut(TCut inp){addToAllCut=inp;}

	///Getters
	inline RDK2CutSet GetCutSet(){return theCutSet;}
	inline double GetBranchingRatio(){return branchingRatio;}
	inline TH1* GetEP_eEHist(){return ep_eEHist;}
	inline TH1* GetEP_pEHist(){return ep_pEHist;}
	inline TH1* GetEP_pTHist(){return ep_pTHist;}
	inline TH1* GetEPG_eEHist(){return epg_eEHist;}
	inline TH1* GetEPG_pEHist(){return epg_pEHist;}
	inline TH1* GetEPG_pTHist(){return epg_pTHist;}
	inline TH1* GetEPBG_eEHist(){return epbg_eEHist;}
	inline TH1* GetEPBG_pEHist(){return epbg_pEHist;}
	inline TH1* GetEPBG_pTHist(){return epbg_pTHist;}
	inline TH1* GetEPG_gEAvgHist(){return epg_gEAvgHist;}
	inline TH1* GetEPG_gEAvgVarHist(){return epg_gEAvgVarHist;}
	inline TH1* GetEPG_gEDetHists(int inp){return epg_gEDetHists[inp-1];}
	inline TH1* GetEPG_gEDetVarHists(int inp){return epg_gEDetVarHists[inp-1];}
	inline TH1* GetEPBG_bGEAvgHist(){return epbg_bGEAvgHist;}
	inline TH1* GetEPBG_bGEAvgVarHist(){return epbg_bGEAvgVarHist;}
	inline TH1* GetEPBG_bGEDetHists(int inp){return epbg_bGEDetHists[inp-1];}
	inline TH1* GetEPBG_bGEDetVarHists(int inp){return epbg_bGEDetVarHists[inp-1];}
	inline TString GetP3ID(){return threeBodySet.GetPResultsID();}
	inline TString GetE3ID(){return threeBodySet.GetEResultsID();}
	inline TString GetP4ID(){return fourBodySet.GetPResultsID();}
	inline TString GetE4ID(){return fourBodySet.GetEResultsID();}
	inline TString GetG4ID(){return fourBodySet.GetGResultsID();}
	inline double GetEPPer3Decay(double& error){error=epPer3DecayError; return epPer3Decay;}
	inline double GetEPGPer4Decay(double& error){error=epgPer4DecayError; return epgPer4Decay;}
	inline double GetEPGPerEP(double& error){error=epgPerEPError; return epgPerEP;}
	inline double GetEPBGPer4Decay(double& error){error=epbgPer4DecayError; return epbgPer4Decay;}
	inline double GetEPBGPerEP(double& error){error=epbgPerEPError; return epbgPerEP;}

	TH1* MakeHist(TString histName, EveType inpEveType, TString drawString, TCut inpCut, TH1* histPrototype);
	TH1* MakeHist(TString histName, EveType inpEveType, TString drawString, TCut inpCut, HistDim inpHistDimX, HistDim inpHistDimY = BLANK_DIM, HistDim inpHistDimZ = BLANK_DIM);

	TH1* MakeRatioHist(TString histName, TString drawString, TCut inpCut, TH1* histPrototype);
	TH1* MakeRatioHist(TString histName, TString drawString, TCut inpCut, HistDim inpHistDimX, HistDim inpHistDimY = BLANK_DIM, HistDim inpHistDimZ = BLANK_DIM);

	TH1* MakeRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, TH1* histPrototype);  //Make Rad Hist using cutset
	TH1* MakeRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, HistDim inpHistDimX, HistDim inpHistDimY = BLANK_DIM, HistDim inpHistDimZ = BLANK_DIM); //Make Rad Hist using cutset

	TH1* MakeFalseRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, TH1* histPrototype);  //Make Rad Hist using cutset
	TH1* MakeFalseRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, HistDim inpHistDimX, HistDim inpHistDimY = BLANK_DIM, HistDim inpHistDimZ = BLANK_DIM); //Make Rad Hist using cutset

	TH1* MakeRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, TH1* histPrototype);  //Make RadAverage Hist using cutset
	TH1* MakeRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, HistDim inpHistDimX, HistDim inpHistDimY = BLANK_DIM, HistDim inpHistDimZ = BLANK_DIM); //Make Rad Average Hist using cutset

	TH1* MakeFalseRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, TH1* histPrototype);  //Make RadAverage Hist using cutset
	TH1* MakeFalseRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, HistDim inpHistDimX, HistDim inpHistDimY = BLANK_DIM, HistDim inpHistDimZ = BLANK_DIM); //Make Rad Average Hist using cutset

	TH1* MakeStandardHist(TString histName, CoDet detType, TString drawString, HistDim inpHistDimX);

	TH1* MakeMultiplicityHist();

	void MakeAllStandardHists();

	TH1* makeEEHistDecayWidthProb(TString nameString, double littleb, double normConstant,HistDim histDim);

	double CalcEPer3Decay(double& error);
	double CalcPPer3Decay(double& error);
	double CalcEPer4Decay(double& error);
	double CalcPPer4Decay(double& error);
	double CalcGPer4Decay(CoDet detType, double& error);
	double CalcEPPer4Decay(double& error);
	double CalcEPPer3Decay(double& error); //Calculates epPerDecay using theCutSet
	double CalcEPGPer4Decay(CoDet detType, double& error); //Calculates epgPerDecay using theCutSet
	double CalcEPGPerEP(CoDet detType, double& error); //Calculates epgPerEP using theCutSet
	double CalcEPFalseGPerEP(CoDet detType, double& error); //Calculates detected EPGs from electrons from three body events
	void CalcRates(); //Calculates EPG and EP rates
	void PrintRates(); //Prints rates

protected:
	///Variables
	RDK2Set threeBodySet; //Three body MC RDK2Set
	RDK2Set fourBodySet;  //Four body MC RDK2Set

	RDK2CutSet theCutSet;  //E, P, and G cut set

	double epPer3Decay; //EP Detection Rate per decay
	double epPer3DecayError; //EP Detection Rate per decay error

	double epgPer4Decay; //EPG Detection Rate per decay
	double epgPer4DecayError; //EPG Detection Rate per decay error

	double epbgPer4Decay; //EPG Detection Rate per decay
	double epbgPer4DecayError; //EPG Detection Rate per decay error

	double epgPerEP; //EPG/EP detection ratio
	double epgPerEPError; //EPG/EP detection ratio error

	double epbgPerEP; //EPG/EP detection ratio
	double epbgPerEPError; //EPG/EP detection ratio error

	double branchingRatio; //Branching Ratio for fourbody decays

	int currentDWCut;  //Index of little b tree

	TCut addToAllCut;

	///Hists
	TH1* ep_eEHist;
	TH1* ep_pEHist;
	TH1* ep_pTHist;
	TH1* epg_eEHist;
	TH1* epg_pEHist;
	TH1* epg_pTHist;
	TH1* epbg_eEHist;
	TH1* epbg_pEHist;
	TH1* epbg_pTHist;
	TH1* epg_gEAvgHist;
	TH1* epg_gEAvgVarHist;
	TH1* epg_gEDetHists[MC_NUMBGO];
	TH1* epg_gEDetVarHists[MC_NUMBGO];
	TH1* epbg_bGEAvgHist;
	TH1* epbg_bGEAvgVarHist;
	TH1* epbg_bGEDetHists[MC_NUMBAPD];
	TH1* epbg_bGEDetVarHists[MC_NUMBAPD];

	///Methods
	void ResetHists();
	void NULLHists();

public:
	ClassDef(RDK2MCAnalyzer,1);//Class to analyze RDK
};

RDK2MCAnalyzer* getRDK2MCAnalysis(TString inpString);
#endif // RDK2MCANALYSIS_H
