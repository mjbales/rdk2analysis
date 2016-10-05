#include "RDK2MCAnalyzer.h"

#include <iostream>
#include <time.h>

#include "TString.h"
#include "TH1.h"
#include "TNamed.h"

#include "MRMisc.h"
#include "MRText.h"
#include "MRIO.h"

#include "RDK2Constants.h"
#include "RDK2Set.h"
#include "RDK2Vars.h"
#include "RDK2CutSet.h"

using namespace std;

ClassImp(RDK2MCAnalyzer);

RDK2MCAnalyzer::RDK2MCAnalyzer()
{
	SetName("RDK2MCAnalysis");
	epPer3Decay = 0;
	epPer3DecayError = 0;
	epgPer4Decay = 0;
	epgPer4DecayError = 0;
	epgPerEP = 0;
	epgPerEPError = 0;
	epbgPer4Decay = 0;
	epbgPer4DecayError = 0;
	epbgPerEP = 0;
	epbgPerEPError = 0;
	branchingRatio = 0.010277;
	currentDWCut = 0;
	NULLHists();
}



RDK2MCAnalyzer::RDK2MCAnalyzer(TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet, bool hmgEvents)
: RDK2MCAnalyzer{}
{
	SetName(inpName);
	SetTitle(inpTitle);

	threeBodySet.SetNumFiles(num3);
	threeBodySet.SetPResultsID(p3RID);
	threeBodySet.SetEResultsID(e3RID);
	if(hmgEvents)
	{
		threeBodySet.SetEventTypeString("Hmg_3B");
		threeBodySet.SetEventDateString("160323");
		threeBodySet.SetIsHmg(true);
	}
	else
	{
		threeBodySet.SetEventTypeString("STD_3B");
		threeBodySet.SetEventDateString("130429");
		threeBodySet.SetIsHmg(false);

	}
	threeBodySet.SetIsRad(false);

	fourBodySet.SetNumFiles(num4);
	fourBodySet.SetPResultsID(p4RID);
	fourBodySet.SetEResultsID(e4RID);
	fourBodySet.SetGResultsID(g4RID);
	fourBodySet.SetEventTypeString("STD_4B_100eV_Endpoint");
	fourBodySet.SetEventDateString("130429");
	fourBodySet.SetIsRad(true);

	theCutSet = inpCutSet;
}

RDK2MCAnalyzer::RDK2MCAnalyzer(TString inpName, TString inpTitle, RDK2Set inp3Set, RDK2Set inp4Set, RDK2CutSet inpCutSet, bool hmgEvents)
: RDK2MCAnalyzer{}
{
	threeBodySet = inp3Set;
	fourBodySet = inp4Set;
	theCutSet = inpCutSet;
	SetName(inpName);
	SetTitle(inpTitle);
}

RDK2MCAnalyzer::~RDK2MCAnalyzer()
{
	ResetHists();
}

void RDK2MCAnalyzer::NULLHists()
{
	ep_eEHist = nullptr;
	ep_pEHist = nullptr;
	ep_pTHist = nullptr;
	epg_eEHist = nullptr;
	epg_pEHist = nullptr;
	epg_pTHist = nullptr;
	epbg_eEHist = nullptr;
	epbg_pEHist = nullptr;
	epbg_pTHist = nullptr;
	epg_gEAvgHist = nullptr;
	epg_gEAvgVarHist = nullptr;
	for (int i = 0; i < MC_NUMBGO; i++)
	{
		epg_gEDetHists[i] = nullptr;
		epg_gEDetVarHists[i] = nullptr;
	}

	for (int i = 0; i < MC_NUMBAPD; i++)
	{
		epbg_bGEDetHists[i] = nullptr;
		epbg_bGEDetVarHists[i] = nullptr;
	}
	epbg_bGEAvgHist = nullptr;
	epbg_bGEAvgVarHist = nullptr;
}

void RDK2MCAnalyzer::ResetHists()
{
	delete ep_eEHist;
	delete ep_pEHist;
	delete ep_pTHist;
	delete epg_eEHist;
	delete epg_pEHist;
	delete epg_pTHist;
	delete epbg_eEHist;
	delete epbg_pEHist;
	delete epbg_pTHist;
	delete epg_gEAvgHist;
	delete epg_gEAvgVarHist;
	for (int i = 0; i < MC_NUMBGO; i++)
	{
		delete epg_gEDetHists[i];
		delete epg_gEDetVarHists[i];
	}
	for (int i = 0; i < MC_NUMBAPD; i++)
	{
		delete epbg_bGEDetHists[i];
		delete epbg_bGEDetVarHists[i];
	}

	delete epbg_bGEAvgHist;
	delete epbg_bGEAvgVarHist;
	NULLHists();
}

TH1* RDK2MCAnalyzer::MakeHist(TString histName, EveType inpEveType, TString drawString, TCut inpCut, HistDim inpHistDimX, HistDim inpHistDimY, HistDim inpHistDimZ)
{
	gROOT->cd();
	TH1* histPrototype;
	TH1* outHist;

	if(inpHistDimZ.numBins > 0)
	{
		histPrototype = makeTH3DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY, inpHistDimZ);
	}
	else if(inpHistDimY.numBins > 0)
	{
		histPrototype = makeTH2DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY);
	}
	else
	{
		histPrototype = makeTH1DFromDim("prototype" + histName, "", inpHistDimX);
	}

	outHist = MakeHist(histName, inpEveType, drawString, inpCut, histPrototype);
	delete histPrototype;
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeHist(TString histName, EveType inpEveType, TString drawString, TCut inpCut, TH1* histPrototype)
{
	TH1* outHist = nullptr;
	TCut finalCut=inpCut+addToAllCut;
	if(inpEveType == THREEBODY)
	{
		outHist = (TH1*) threeBodySet.makeHist(histName, drawString, finalCut, histPrototype);
	}
	else
	{
		outHist = (TH1*) fourBodySet.makeHist(histName, drawString, finalCut, histPrototype);
	}
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeRatioHist(TString histName, TString drawString, TCut inpCut, HistDim inpHistDimX, HistDim inpHistDimY, HistDim inpHistDimZ)
{
	gROOT->cd();
	TH1* histPrototype;
	TH1* outHist;

	if(inpHistDimZ.numBins > 0)
	{
		histPrototype = makeTH3DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY, inpHistDimZ);
	}
	else if(inpHistDimY.numBins > 0)
	{
		histPrototype = makeTH2DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY);
	}
	else
	{
		histPrototype = makeTH1DFromDim("prototype" + histName, "", inpHistDimX);
	}

	outHist = MakeRatioHist(histName, drawString, inpCut, histPrototype);
	delete histPrototype;
	return outHist;

}

TH1* RDK2MCAnalyzer::MakeRatioHist(TString histName, TString drawString, TCut inpCut, TH1* histPrototype)
{
	TH1* outHist = nullptr;
	outHist = MakeHist(histName, FOURBODY, drawString, inpCut, histPrototype);
	double error;
	if(epPer3Decay == 0)
	{
		CalcEPPer3Decay(error);
	}
	double num4Events = fourBodySet.GetNumEvents();
	outHist->Scale(branchingRatio / (epPer3Decay * num4Events));

	return outHist;
}

TH1* RDK2MCAnalyzer::MakeRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, TH1* histPrototype)
{
	TH1* outHist = nullptr;
	TCut additionalCut;

	if(detType == DET_EPG || detType == DET_G)
	{
		if(drawString == "")
		{
			if(gOption == 0)
			{
				drawString += "gDChn.BGO" + int2str(det) + "EDepBlurLOAllExpParam";
			}
			else if(gOption == 1) //multiplicity 2
			{
				drawString += "gDChn.BGO" + int2str(det) + "EDepBlurLOAllExpParam";
				additionalCut = "gDChn.BGOmultiplicityLOAllExpParam==2";
			}
			else if(gOption == 2) //Nonprop high
			{
				drawString += "gDChn.BGO" + int2str(det) + "EDepBlurLOAllExpParamHigh";
			}
			else if(gOption == 3) //Nonprop Low
			{
				drawString += "gDChn.BGO" + int2str(det) + "EDepBlurLOAllExpParamLow";
			}
			else
			{
				cout << "Error G Option not valid!" << endl;
			}
		}
		TCut theCuts;
		if(detType == DET_G)
		{
			theCuts = theCutSet.CreateMCGCut(gOption, det, detType) && additionalCut;
		}
		else
		{
			theCuts = theCutSet.CreateMCEPGCut(1, 0, gOption, det, detType) && additionalCut;
		}
		outHist = MakeRatioHist(histName, drawString, theCuts, histPrototype);

	}
	else if(detType == DET_EPBG || detType == DET_BG)
	{
		if(drawString == "")
		{
			if(gOption == 0)
			{
				drawString += "gChn.BAPD" + int2str(det) + "EDepTotal";
			}
			else if(gOption == 1) //multiplicity 2
			{
				drawString += "gChn.BAPD" + int2str(det) + "EDepTotal";
				additionalCut = "gDChn.BAPDmultiplicity==2";
			}
			else
			{
				cout << "Error G Option not valid!" << endl;
			}
		}

		TCut theCuts;
		if(detType == DET_BG)
		{
			theCuts = theCutSet.CreateMCGCut(gOption, det, detType) && additionalCut;
		}
		else
		{
			theCuts = theCutSet.CreateMCEPGCut(1, 0, gOption, det, detType) && additionalCut;
		}
		outHist = MakeRatioHist(histName, drawString, theCuts, histPrototype);
	}

	return outHist;
}

TH1* RDK2MCAnalyzer::MakeRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, HistDim inpHistDimX, HistDim inpHistDimY, HistDim inpHistDimZ)
{
	gROOT->cd();
	TH1* histPrototype;
	TH1* outHist;

	if(inpHistDimZ.numBins > 0)
	{
		histPrototype = makeTH3DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY, inpHistDimZ);
	}
	else if(inpHistDimY.numBins > 0)
	{
		histPrototype = makeTH2DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY);
	}
	else
	{
		histPrototype = makeTH1DFromDim("prototype" + histName, "", inpHistDimX);
	}

	outHist = MakeRadHist(histName, detType, drawString, gOption, det, histPrototype);
	delete histPrototype;
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeFalseRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, TH1* histPrototype)
{
	TH1* outHist = nullptr;

	if(detType == DET_EPG)
	{
		if(drawString == "")
		{
			if(gOption == 0)
			{
				drawString += "eDChn.BGO" + int2str(det) + "EDepBlurLOAllExpParam";
			}
			else
			{
				cout << "Error G Option not valid!" << endl;
			}
		}
		outHist = MakeHist(histName, THREEBODY, drawString, theCutSet.CreateMCEPFalseGCut(1, 0, gOption, det, detType), histPrototype);

	}
	else if(detType == DET_EPBG)
	{
		if(drawString == "")
		{
			if(gOption == 0)
			{
				drawString += "eChn.BAPD" + int2str(det) + "EDepTotal";
			}
			else
			{
				cout << "Error G Option not valid!" << endl;
			}
		}

		outHist = MakeHist(histName, THREEBODY, drawString, theCutSet.CreateMCEPFalseGCut(1, 0, gOption, det, detType), histPrototype);
	}

	return outHist;
}

TH1* RDK2MCAnalyzer::MakeFalseRadHist(TString histName, CoDet detType, TString drawString, int gOption, int det, HistDim inpHistDimX, HistDim inpHistDimY, HistDim inpHistDimZ)
{
	gROOT->cd();
	TH1* histPrototype;
	TH1* outHist;

	if(inpHistDimZ.numBins > 0)
	{
		histPrototype = makeTH3DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY, inpHistDimZ);
	}
	else if(inpHistDimY.numBins > 0)
	{
		histPrototype = makeTH2DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY);
	}
	else
	{
		histPrototype = makeTH1DFromDim("prototype" + histName, "", inpHistDimX);
	}

	outHist = MakeFalseRadHist(histName, detType, drawString, gOption, det, histPrototype);
	delete histPrototype;
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, TH1* histPrototype)
{
	int numDets = 0;
	if(detType == DET_EPG || detType == DET_G)
	{
		numDets = MC_NUMBGO;
	}
	else if(detType == DET_EPBG || detType == DET_BG)
	{
		numDets = MC_NUMBAPD;
	}
	TH1* outHist = (TH1*) histPrototype->Clone(histName);
	outHist->Sumw2();
	for (int i = 0; i < numDets; i++)
	{
		TH1* tempHist = MakeRadHist("radIndivTemp", detType, drawString, gOption, i + 1, histPrototype);
		outHist->Add(tempHist, 1. / numDets);
		delete tempHist;
	}
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, HistDim inpHistDimX, HistDim inpHistDimY, HistDim inpHistDimZ)
{
	gROOT->cd();
	TH1* histPrototype;
	TH1* outHist;

	if(inpHistDimZ.numBins > 0)
	{
		histPrototype = makeTH3DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY, inpHistDimZ);
	}
	else if(inpHistDimY.numBins > 0)
	{
		histPrototype = makeTH2DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY);
	}
	else
	{
		histPrototype = makeTH1DFromDim("prototype" + histName, "", inpHistDimX);
	}

	outHist = MakeRadAvgHist(histName, detType, drawString, gOption, histPrototype);
	delete histPrototype;
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeFalseRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, TH1* histPrototype)
{
	int numDets = 0;
	if(detType == DET_EPG)
	{
		numDets = MC_NUMBGO;
	}
	else if(detType == DET_EPBG)
	{
		numDets = MC_NUMBAPD;
	}
	TH1* outHist = (TH1*) histPrototype->Clone(histName);
	outHist->Sumw2();
	for (int i = 0; i < numDets; i++)
	{
		TH1* tempHist = MakeFalseRadHist("radIndivTemp", detType, drawString, gOption, i + 1, histPrototype);
		outHist->Add(tempHist, 1. / numDets);
		delete tempHist;
	}
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeFalseRadAvgHist(TString histName, CoDet detType, TString drawString, int gOption, HistDim inpHistDimX, HistDim inpHistDimY, HistDim inpHistDimZ)
{
	gROOT->cd();
	TH1* histPrototype;
	TH1* outHist;

	if(inpHistDimZ.numBins > 0)
	{
		histPrototype = makeTH3DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY, inpHistDimZ);
	}
	else if(inpHistDimY.numBins > 0)
	{
		histPrototype = makeTH2DFromDim("prototype" + histName, "", inpHistDimX, inpHistDimY);
	}
	else
	{
		histPrototype = makeTH1DFromDim("prototype" + histName, "", inpHistDimX);
	}

	outHist = MakeFalseRadAvgHist(histName, detType, drawString, gOption, histPrototype);
	delete histPrototype;
	return outHist;
}

TH1* RDK2MCAnalyzer::MakeStandardHist(TString histName, CoDet detType, TString drawString, HistDim inpHistDimX)
{
	TH1* outHist;
	if(detType == DET_EPBG || detType == DET_EPG)
	{
		outHist = MakeRadAvgHist(histName, detType, drawString, 0, inpHistDimX);
	}
	else
	{
		outHist = MakeHist(histName, THREEBODY, drawString, theCutSet.CreateMCEPCut(1, 0), inpHistDimX);

	}
	return outHist;
}

//void RDK2MCAnalysis::makeAllStandardHists()
//{
//    cout << "Making all standard histograms for: " << GetName() << endl;
//    time_t startTime,endTime;
//    time (&startTime);  //Start Time for Run
//    epPer3Decay=threeBodySet.makeAllStandardHists(GetName(),epPer3DecayError,theCutSet,DET_EP);
//    epgPer4Decay=fourBodySet.makeAllStandardHists(GetName(),epgPer4DecayError,theCutSet,DET_EPG,epPer3Decay,branchingRatio);
//    epbgPer4Decay=fourBodySet.makeAllStandardHists(GetName(),epbgPer4DecayError,theCutSet,DET_EPBG,epPer3Decay,branchingRatio);
//    double tempError;
//    CalcEPGPerEP(DET_EPG,tempError);
//    CalcEPGPerEP(DET_EPBG,tempError);
//    PrintRates();
//    time (&endTime);  //End Time for Run
//
//    cout << "New Method: " << endTime-startTime << endl;
//    typeAnythingToContinue("");
//}

TH1* RDK2MCAnalyzer::MakeMultiplicityHist()
{
	HistDim multdim = { 13, -.5, 12.5 };
	TH1* outHist = MakeRatioHist("multHist", "gDChn.BGOmultiplicityLOAllExpParam", theCutSet.CreateMCEPCut(), multdim);
	return outHist;
}


TH1* RDK2MCAnalyzer::makeEEHistDecayWidthProb(TString nameString,double littleb,double normConst,HistDim histDim)
{

	TH1* eEHist=threeBodySet.makeEEHistDecayWidthProb(nameString, theCutSet,littleb, normConst,histDim);

	return eEHist;

}
void RDK2MCAnalyzer::MakeAllStandardHists()
{
	cout << "Making all standard histograms for: " << GetName() << endl;
	ResetHists();
	time_t startTime, endTime;
	time(&startTime);  //Start Time for Run

//    CalcRates();  //Calculate EP and EPG/EPBG rates

	TString HistsDir = HISTS_DIR;
	TString baseHistName = TString("Hist_MC_") + GetName() + "_";
	TString endHistName = ".txt";

	TH1** eEHists[3] = { &ep_eEHist, &epg_eEHist, &epbg_eEHist };
	TH1** pEHists[3] = { &ep_pEHist, &epg_pEHist, &epbg_pEHist };
	TH1** pTHists[3] = { &ep_pTHist, &epg_pTHist, &epbg_pTHist };

	CoDet detType;
	TString coString, yAxisString;
	int numDets = 0;
	int numTypesToRun=1;
	threeBodySet.setEventList(theCutSet.CreateMCEPCut(1, 0));
	if(fourBodySet.GetNumFiles() > 0)
	{
		numTypesToRun=3;
		fourBodySet.setEventList(theCutSet.CreateMCEPCut(1, 0));
	}




	for (int i = 0; i < numTypesToRun; i++)
	{
		if(i == 0)
		{
			detType = DET_EP;
			coString = "EP";
			yAxisString = "Counts";

		}
		else if(i == 1)
		{

			detType = DET_EPG;
			coString = "EPG";
			numDets = MC_NUMBGO;
			yAxisString = "#frac{epg}{ep}";
		}
		else
		{
			detType = DET_EPBG;
			coString = "EPBG";
			numDets = MC_NUMBAPD;
			yAxisString = "#frac{epbg}{ep}";
		}

		*(pTHists[i]) = MakeStandardHist(coString + "_pTHist", detType, "pChn.SBDTimeFirst-eChn.SBDTimeFirst", TOFP_DIM5);
		(*(pTHists[i]))->SetTitle("P-E Time Difference " + coString + ";P-E Time Difference (s);" + yAxisString);
		convertTH1ToTXT(*(pTHists[i]), HistsDir + baseHistName + coString + "_ToFP" + endHistName);
		if(i == 0)
		{
			epPer3Decay = (*(pTHists[i]))->IntegralAndError(1, (*(pTHists[i]))->GetNbinsX(), epPer3DecayError);
			epPer3Decay /= (double) threeBodySet.GetNumEvents();
			epPer3DecayError /= (double) threeBodySet.GetNumEvents();
		}

		*(pEHists[i]) = MakeStandardHist(coString + "_pEHist", detType, "pDChn.SBDEDepBlur", EDEPP_DIM2);  //modified to blur 140430
		(*(pEHists[i]))->SetTitle("Proton Energy " + coString + ";SBD Energy Channel (782=782 keV);" + yAxisString);
		convertTH1ToTXT(*(pEHists[i]), HistsDir + baseHistName + coString + "_EDepP" + endHistName);

		*(eEHists[i]) = MakeStandardHist(coString + "_eEHist", detType, "eDChn.SBDEDepBlur", EDEPE_DIM2);
		(*(eEHists[i]))->SetTitle("Electron Energy " + coString + ";SBD Energy Channel (782=782 keV);" + yAxisString);
		convertTH1ToTXT(*(eEHists[i]), HistsDir + baseHistName + coString + "_EDepE" + endHistName);

		if(i == 1)
		{
			epg_gEAvgHist = makeTH1DFromDim("epg_gEAvgHist", "epg_gEAvgHist", BGO_DIM2);
			epg_gEAvgHist->SetTitle("Average BGO Gamma Energy " + coString + ";BGO Energy Channel (511=511 keV);" + yAxisString);
			epg_gEAvgVarHist = new TH1D("epg_gEAvgVarHist", "epg_gEAvgHist", 20, BGO_VARHIST_LOWEDGES);
			epg_gEAvgVarHist->SetTitle("Average BGO Gamma Energy " + coString + ";BGO Energy Channel (511=511 keV);" + yAxisString);
			for (int i = 0; i < numDets; i++)
			{
				epg_gEDetHists[i] = MakeRadHist("epg_gEDet" + int2str(i + 1) + "Hist", detType, "", 0, i + 1, BGO_DIM2);
				epg_gEDetHists[i]->SetTitle("BGO Det " + int2str(i + 1) + " Gamma Energy " + coString + ";BGO Energy Channel (511=511 keV);" + yAxisString);
				epg_gEDetVarHists[i] = MakeRadHist("epg_gEDet" + int2str(i + 1) + "VarHist", detType, "", 0, i + 1, epg_gEAvgVarHist);
				epg_gEDetVarHists[i]->SetTitle("BGO Det " + int2str(i + 1) + " Gamma Energy " + coString + ";BGO Energy Channel (511=511 keV);" + yAxisString);
				convertTH1ToTXT(epg_gEDetHists[i], HistsDir + baseHistName + coString + "_EDepGDet" + int2str(i + 1) + endHistName);
				convertTH1ToTXT(epg_gEDetVarHists[i], HistsDir + TString("VarHist_MC_") + GetName() + "_" + coString + "_EDepGDet" + int2str(i + 1) + endHistName);
				epg_gEAvgHist->Add(epg_gEDetHists[i], 1. / numDets);
				epg_gEAvgVarHist->Add(epg_gEDetVarHists[i], 1. / numDets);
			}
			convertTH1ToTXT(epg_gEAvgHist, HistsDir + baseHistName + coString + "_EDepGAvg" + endHistName);
			convertTH1ToTXT(epg_gEAvgVarHist, HistsDir + TString("VarHist_MC_") + GetName() + "_" + coString + "_EDepGAvg" + endHistName);

			epgPer4Decay = epg_gEAvgHist->IntegralAndError(1, epg_gEAvgHist->GetNbinsX(), epgPer4DecayError);
			epgPer4Decay *= epPer3Decay / branchingRatio;
			epgPer4DecayError *= epPer3Decay / branchingRatio;
			double temp;
			CalcEPGPerEP(DET_EPG, temp);

		}
		else if(i == 2)
		{
			epbg_bGEAvgHist = makeTH1DFromDim("epbg_bGEAvgHist", "epbg_bGEAvgHist", BAPD_DIM3);
			epbg_bGEAvgHist->SetTitle("Average BAPD Gamma Energy " + coString + ";BAPD Energy Channel (5.9=5.9 keV);" + yAxisString);
			epbg_bGEAvgVarHist = new TH1D("epbg_bGEAvgVarHist", "epbg_bGEAvgHist", 20, BGO_VARHIST_LOWEDGES); //Incorrect low edges
			epbg_bGEAvgVarHist->SetTitle("Average BAPD Gamma Energy " + coString + ";BAPD Energy Channel (5.9=5.9 keV);" + yAxisString);
			double weights[3] = { 0.8540694268, 1.0729652866, 1.0729652866 }; //Email from Jeff 151122

			for (int i = 0; i < numDets; i++)
			{
				epbg_bGEDetHists[i] = MakeRadHist("epbg_bGEDet" + int2str(i + 1) + "Hist", detType, "", 0, i + 1, BAPD_DIM3);
				epbg_bGEDetHists[i]->SetTitle("BAPD Det " + int2str(i + 1) + " Gamma Energy " + coString + ";BAPD Energy Channel (5.9=5.9 keV);" + yAxisString);
				epbg_bGEDetVarHists[i] = MakeRadHist("epbg_bGEDet" + int2str(i + 1) + "VarHist", detType, "", 0, i + 1, epbg_bGEAvgVarHist);
				epbg_bGEDetVarHists[i]->SetTitle("BAPD Det " + int2str(i + 1) + " Gamma Energy " + coString + ";BAPD Energy Channel (5.9=5.9 keV);" + yAxisString);
				convertTH1ToTXT(epbg_bGEDetHists[i], HistsDir + baseHistName + coString + "_EDepBGDet" + int2str(i + 1) + endHistName);
				convertTH1ToTXT(epbg_bGEDetVarHists[i], HistsDir + TString("VarHist_MC_") + GetName() + "_" + coString + "_EDepBGDet" + int2str(i + 1) + endHistName);
				epbg_bGEAvgHist->Add(epbg_bGEDetHists[i], weights[i] / numDets);
				epbg_bGEAvgVarHist->Add(epbg_bGEDetVarHists[i], weights[i] / numDets);
			}
			convertTH1ToTXT(epbg_bGEAvgHist, HistsDir + baseHistName + coString + "_EDepBGAvg" + endHistName);
			convertTH1ToTXT(epbg_bGEAvgVarHist, HistsDir + TString("VarHist_MC_") + GetName() + "_" + coString + "_EDepBGAvg" + endHistName);

			epbgPer4Decay = epbg_bGEAvgHist->IntegralAndError(1, epbg_bGEAvgHist->GetNbinsX(), epbgPer4DecayError);
			epbgPer4Decay *= epPer3Decay / branchingRatio;
			epbgPer4DecayError *= epPer3Decay / branchingRatio;
			double temp;
			CalcEPGPerEP(DET_EPBG, temp);
		}

	}

	PrintRates();  //Print the rates at the end
	time(&endTime);  //End Time for Run

}

double RDK2MCAnalyzer::CalcPPer3Decay(double& error)
{
	HistDim simpleDim = { 3, -10, 800 };
	TH1* pRateHist = MakeHist("pRateHist", THREEBODY, "pChn.SBDEDepTotal", theCutSet.CreateMCPCut(1), simpleDim);
	double pPer3DecayError;
	double pPer3Decay = pRateHist->IntegralAndError(1, 3, pPer3DecayError);
	delete pRateHist;
	pPer3Decay /= (double) threeBodySet.GetNumEvents();
	pPer3DecayError /= (double) threeBodySet.GetNumEvents();
	error = pPer3DecayError;
	return pPer3Decay;
}

double RDK2MCAnalyzer::CalcEPer3Decay(double& error)
{
	HistDim simpleDim = { 3, -10, 800 };
	TH1* eRateHist = MakeHist("eRateHist", THREEBODY, "eChn.SBDEDepTotal", theCutSet.CreateMCECut(0), simpleDim);
	double ePer3DecayError;
	double ePer3Decay = eRateHist->IntegralAndError(1, 3, ePer3DecayError);
	delete eRateHist;
	ePer3Decay /= (double) threeBodySet.GetNumEvents();
	ePer3DecayError /= (double) threeBodySet.GetNumEvents();
	error = ePer3DecayError;
	return ePer3Decay;
}

double RDK2MCAnalyzer::CalcPPer4Decay(double& error)
{
	HistDim simpleDim = { 3, -10, 800 };
	TH1* pRateHist = MakeHist("pRateHist", FOURBODY, "pChn.SBDEDepTotal", theCutSet.CreateMCPCut(1), simpleDim);
	double pPer4DecayError;
	double pPer4Decay = pRateHist->IntegralAndError(1, 3, pPer4DecayError);
	delete pRateHist;
	pPer4Decay /= (double) fourBodySet.GetNumEvents();
	pPer4DecayError /= (double) fourBodySet.GetNumEvents();
	error = pPer4DecayError;
	return pPer4Decay;
}

double RDK2MCAnalyzer::CalcEPer4Decay(double& error)
{
	HistDim simpleDim = { 3, -10, 800 };
	TH1* eRateHist = MakeHist("eRateHist", FOURBODY, "eChn.SBDEDepTotal", theCutSet.CreateMCECut(0), simpleDim);
	double ePer4DecayError;
	double ePer4Decay = eRateHist->IntegralAndError(1, 3, ePer4DecayError);
	delete eRateHist;
	ePer4Decay /= (double) fourBodySet.GetNumEvents();
	ePer4DecayError /= (double) fourBodySet.GetNumEvents();
	error = ePer4DecayError;
	return ePer4Decay;
}

double RDK2MCAnalyzer::CalcEPPer3Decay(double& error)
{
	HistDim simpleDim = { 3, -10, 800 };
	TH1* epRateHist = MakeHist("epRateHist", THREEBODY, "pChn.SBDEDepTotal", theCutSet.CreateMCEPCut(1, 0), simpleDim);
	epPer3Decay = epRateHist->IntegralAndError(1, 3, epPer3DecayError);
	delete epRateHist;
	epPer3Decay /= (double) threeBodySet.GetNumEvents();
	epPer3DecayError /= (double) threeBodySet.GetNumEvents();
	error = epPer3DecayError;
	return epPer3Decay;
}

double RDK2MCAnalyzer::CalcEPPer4Decay(double& error)
{
	HistDim simpleDim = { 3, -10, 800 };
	TH1* epRateHist = MakeHist("epRateHist", FOURBODY, "pChn.SBDEDepTotal", theCutSet.CreateMCEPCut(1, 0), simpleDim);
	double epPer4Decay = epRateHist->IntegralAndError(1, 3, epPer3DecayError);
	delete epRateHist;
	epPer4Decay /= (double) fourBodySet.GetNumEvents();
	error /= (double) fourBodySet.GetNumEvents();
	return epPer4Decay;
}

double RDK2MCAnalyzer::CalcGPer4Decay(CoDet detType, double& error)
{
	double answer = 0;
	error = 0;
	HistDim simpleDim = { 3, -10, 800 };
	int numDets = 0;
	TString drawString;
	if(detType == DET_G)
	{
		numDets = 12;
	}
	else if(detType == DET_BG)
	{
		numDets = 3;
	}

	for (int i = 0; i < numDets; i++)
	{
		if(detType == DET_G)
		{

			drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam";
		}
		else if(detType == DET_BG)
		{
			drawString = "gChn.BAPD" + int2str(i + 1) + "EDepTotal";
		}
		TCut theCut = theCutSet.CreateMCGCut(0, i + 1, detType);

		TH1* gRateHist = MakeHist("gRateHist", FOURBODY, drawString, theCut, simpleDim);
		double detError;
		answer += gRateHist->IntegralAndError(1, 3, detError);
		error += detError * detError;
		delete gRateHist;
	}
	answer /= (double) fourBodySet.GetNumEvents();
	error = sqrt(error) / (double) fourBodySet.GetNumEvents();
	;
	return answer;
}

double RDK2MCAnalyzer::CalcEPGPer4Decay(CoDet detType, double& error)
{
	double answer = 0;
	error = 0;
	HistDim simpleDim = { 3, -10, 800 };
	int numDets = 0;
	TString drawString;
	if(detType == DET_G || detType == DET_EPG)
	{
		numDets = 12;
	}
	else if(detType == DET_BG || detType == DET_EPBG)
	{
		numDets = 3;
	}

	for (int i = 0; i < numDets; i++)
	{
		if(detType == DET_G || detType == DET_EPG)
		{

			drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam";
		}
		else if(detType == DET_BG || detType == DET_EPBG)
		{
			drawString = "gChn.BAPD" + int2str(i + 1) + "EDepTotal";
		}
		TCut theCut = theCutSet.CreateMCGCut(0, i + 1, detType) && theCutSet.CreateMCEPCut(1, 0);

		TH1* gRateHist = MakeHist("gRateHist", FOURBODY, drawString, theCut, simpleDim);
		double detError;
		answer += gRateHist->IntegralAndError(1, 3, detError);
		error += detError * detError;
		delete gRateHist;
	}
	answer /= (double) fourBodySet.GetNumEvents();
	error = sqrt(error) / (double) fourBodySet.GetNumEvents();
	return answer;
}

double RDK2MCAnalyzer::CalcEPGPerEP(CoDet detType, double& error)
{
	double tempError;
	double answer;
	if(epPer3Decay == 0)
	{
		CalcEPPer3Decay(tempError);
	}

	if((epgPer4Decay == 0 && detType == DET_EPG) || (epbgPer4Decay == 0 && detType == DET_EPBG))
	{
		CalcEPGPer4Decay(detType, tempError);
	}

	if(detType == DET_EPG)
	{
		answer = epgPerEP = (branchingRatio * epgPer4Decay) / epPer3Decay;
		error = epgPerEPError = branchingRatio * getBasicErrorProgation(epgPer4Decay, epPer3Decay, epgPer4DecayError, epPer3DecayError, DIVISION);
	}
	else
	{
		answer = epbgPerEP = (branchingRatio * epbgPer4Decay) / epPer3Decay;
		error = epbgPerEPError = branchingRatio * getBasicErrorProgation(epbgPer4Decay, epPer3Decay, epbgPer4DecayError, epPer3DecayError, DIVISION);
	}

	return answer;
}

double RDK2MCAnalyzer::CalcEPFalseGPerEP(CoDet detType, double& error)
{
	double answer;
	HistDim simpleDim = { 3, -10, 800 };
	TH1* epgRateHist = MakeFalseRadAvgHist("epgRateHist", detType, "", 0, simpleDim);
	answer = epgRateHist->IntegralAndError(1, 3, error);
	delete epgRateHist;
	answer /= threeBodySet.GetNumEvents() * epPer3Decay;
	error /= threeBodySet.GetNumEvents() * epPer3Decay;
	return answer;
}

void RDK2MCAnalyzer::CalcRates()
{
	double tempError;
	CalcEPPer3Decay(tempError);
	CalcEPGPer4Decay(DET_EPG, tempError);
	CalcEPGPerEP(DET_EPG, tempError);
	CalcEPGPer4Decay(DET_EPBG, tempError);
	CalcEPGPerEP(DET_EPBG, tempError);
}

void RDK2MCAnalyzer::PrintRates()
{
	cout << GetName() << " Rates:" << endl;
	cout << "EP Per Three Body Decay: " << epPer3Decay << " +/- " << epPer3DecayError << endl;
	if(fourBodySet.GetNumFiles() > 0)
	{
		cout << "EPG Per Four Body Decay: " << epgPer4Decay << " +/- " << epgPer4DecayError << endl;
		cout << "EPG Per EP: " << epgPerEP << " +/- " << epgPerEPError << endl;
		cout << "EPBG Per Four Body Decay: " << epbgPer4Decay << " +/- " << epbgPer4DecayError << endl;
		cout << "EPBG Per EP: " << epbgPerEP << " +/- " << epbgPerEPError << endl;
	}
}

RDK2MCAnalyzer* getRDK2MCAnalysis(TString inpString)
{
	TString filePath = ANALYSIS_DIR;
	filePath += inpString + ".root";
	TFile f(filePath);
	RDK2MCAnalyzer* theAnal = (RDK2MCAnalyzer*) f.Get(inpString);
	gROOT->cd();
	theAnal = (RDK2MCAnalyzer*) theAnal->Clone(inpString + "temp");
	theAnal->SetName(inpString);
	f.Close();
	return theAnal;
}
