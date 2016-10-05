#include "RDK2Set.h"

#include <iostream>

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TCut.h"
#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TEventList.h"
#include "TVector3.h"

#include "MRIO.h"
#include "MRText.h"
#include "MRMisc.h"

#include "RDK2Vars.h"
#include "RDK2Constants.h"
#include "RDK2CutSet.h"
#include "RDK2IO.h"
#include "RDK2Analysis.h"
#include "RDK2Events.h"

using namespace std;

RDK2Set::RDK2Set()
{
	//ctor
	SetName("RDK2Set");
	filesLoaded = false;
	isHmg = false;
	isRad = false;
	isEventGenerator = true;
	eventsChain = nullptr;
	pResultsChain = nullptr;
	eResultsChain = nullptr;
	gResultsChain = nullptr;
	littleBChain = nullptr;
	pDervChain = nullptr;
	eDervChain = nullptr;
	gDervChain = nullptr;
	numFiles = 0;
	epRate = 0;
	eventsDateString = "";
	eventsTypeString = "";
	savedEPCut = "";

}

RDK2Set::RDK2Set(TString inpName, TString inpTitle, int inp_numFiles, TString inp_eventsDateString, TString inp_eventsTypeString, TString inp_pResultsID, TString inp_eResultsID, TString inp_gResultsID)
{
	SetName(inpName);
	SetTitle(inpTitle);
	numFiles = inp_numFiles;
	eventsDateString = inp_eventsDateString;
	eventsTypeString = inp_eventsTypeString;
	pResultsID = inp_pResultsID;
	eResultsID = inp_eResultsID;
	gResultsID = inp_gResultsID;

	if(eventsTypeString == "STD_3B")
	{
		isRad = false;
		isHmg = false;
		isEventGenerator = true;
	}
	else if(eventsTypeString == "STD_4B_100eV_Endpoint" || "STD_4B_1keV_Endpoint" || "Deconv_4B")
	{
		isRad = true;
		isHmg = false;
		isEventGenerator = true;
	}
	else if(eventsTypeString == "radlow" || eventsTypeString == "radmid" || eventsTypeString == "radhigh")
	{
		isRad = true;
		isHmg = false;
		isEventGenerator = false;
	}
	else if(eventsTypeString.Contains("Hmg"))
	{
		isRad = false;
		isHmg = true;
		isEventGenerator = true;
	}

	eventsChain = nullptr;
	pResultsChain = nullptr;
	eResultsChain = nullptr;
	gResultsChain = nullptr;
	littleBChain = nullptr;
	pDervChain = nullptr;
	eDervChain = nullptr;
	gDervChain = nullptr;
	filesLoaded = false;
	epRate = 0;
	savedEPCut = "";

}

RDK2Set::~RDK2Set()
{
	if(eventsChain == nullptr) delete eventsChain;

	if(pResultsChain == nullptr) delete pResultsChain;

	if(eResultsChain == nullptr) delete eResultsChain;

	if(gResultsChain == nullptr) delete gResultsChain;

	if(littleBChain == nullptr) delete littleBChain;
}

ClassImp(RDK2Set);

const RDK2Set& RDK2Set::operator=(const RDK2Set& inp)
{
	TString newName = inp.GetName();
	newName += "New";
	SetName(newName);
	filesLoaded = false;
	isHmg = inp.isHmg;
	isRad = inp.isRad;
	isEventGenerator = inp.isEventGenerator;
	eventsChain = nullptr;
	pResultsChain = nullptr;
	eResultsChain = nullptr;
	gResultsChain = nullptr;
	littleBChain = nullptr;
	pDervChain = nullptr;
	eDervChain = nullptr;
	gDervChain = nullptr;
	numFiles = inp.numFiles;
	epRate = 0;
	eventsDateString = "";
	eventsTypeString = "";
	savedEPCut = "";
	return *this;
}

void RDK2Set::loadEventsChain()
{

	int numBodies;
	TString filePath = EVENTS_DIR;

	if(isRad)
	{
		numBodies = 4;
	}
	else
	{
		numBodies = 3;
	}

	if(isEventGenerator)
	{

		filePath = eventsTypeString + "_" + eventsDateString;
		if(eventsChain == nullptr)
		{
			eventsChain = makeEventChain(filePath, "t", numFiles, true);
		}
	}
	else
	{
		filePath += "MJB-Evn-" + eventsTypeString + "-n1.B" + int2str(numBodies) + "-0-20" + eventsDateString + ".root";

		if(eventsChain == nullptr)
		{
			eventsChain = makeEventChain(filePath, "t", numFiles);
		}
	}

}

void RDK2Set::loadPResultsChain()
{

	TString filePath = RESULTS_DIR;

	if(isEventGenerator)
	{
		filePath += "Results_" + pResultsID + "_P_" + eventsTypeString + "_" + eventsDateString + "_Set0.root";
	}
	else
	{
		filePath += "Results" + eventsDateString + "_" + eventsTypeString + "_P_Set0_" + pResultsID + ".root";
	}

	if(pResultsChain == nullptr)
	{
		pResultsChain = makeResultChain(filePath, "g", numFiles);
	}
}

void RDK2Set::loadEResultsChain()
{
	TString filePath = RESULTS_DIR;

	if(isEventGenerator)
	{
		filePath += "Results_" + eResultsID + "_E_" + eventsTypeString + "_" + eventsDateString + "_Set0.root";
	}
	else
	{
		filePath += "Results" + eventsDateString + "_" + eventsTypeString + "_E_Set0_" + eResultsID + ".root";
	}

	if(eResultsChain == nullptr)
	{
		eResultsChain = makeResultChain(filePath, "g", numFiles);
	}
}

void RDK2Set::loadGResultsChain()
{
	TString filePath = RESULTS_DIR;

	if(isEventGenerator)
	{
		filePath += "Results_" + gResultsID + "_G_" + eventsTypeString + "_" + eventsDateString + "_Set0.root";
	}
	else
	{
		filePath += "Results" + eventsDateString + "_" + eventsTypeString + "_G_Set0_" + gResultsID + ".root";
	}

	if(gResultsChain == nullptr)
	{
		gResultsChain = makeResultChain(filePath, "g", numFiles);
	}
}

void RDK2Set::loadLittleBChain()
{

	TString eventString = eventsTypeString + "_" + eventsDateString;
	if(littleBChain == nullptr)
	{
		littleBChain = makeLilBChain(eventString, DEFAULT_LIL_B_TREENAME, numFiles);
	}

}

void RDK2Set::loadDervChains()
{
	TString resultsID, letter;
	ParType parType;

	for (int i = 0; i < 3; i++)
	{

		if(i == 0)
		{
			letter = "P";
			resultsID = pResultsID;
			parType = PROTON;
		}
		else if(i == 1)
		{
			parType = ELECTRON;
			letter = "E";
			resultsID = eResultsID;
		}
		else if(i == 2)
		{
			parType = PHOTON;
			letter = "G";
			resultsID = gResultsID;
			if(!isRad) //No Gamma if hmg
				break;
		}

		//If first file exists we're going to assume they all exist
		TString derivedPath = RESULTS_DIR;
		if(isEventGenerator)
		{
			derivedPath += "Results_" + resultsID + "_" + letter + "_" + eventsTypeString + "_" + eventsDateString + "_Derv_Set0.root";
		}
		else
		{
			derivedPath += "Results" + eventsDateString + "_" + eventsTypeString + "_" + letter + "_Set0_" + resultsID + "_Derv.root";
		}
		if(!FileExists(derivedPath.Data())) makeDerivedResultsFiles(parType);
		if(i == 0)
		{
			pDervChain = makeResultChain(derivedPath.Data(), DEFAULT_RESULTSDERV_TREENAME, numFiles);
		}
		else if(i == 1)
		{
			eDervChain = makeResultChain(derivedPath.Data(), DEFAULT_RESULTSDERV_TREENAME, numFiles);
		}
		else if(i == 2)
		{
			gDervChain = makeResultChain(derivedPath.Data(), DEFAULT_RESULTSDERV_TREENAME, numFiles);
		}
	}

}

void RDK2Set::loadChainsAndFriend()
{
	loadEventsChain();

	if(pResultsID != "")
	{
		loadPResultsChain();
		eventsChain->AddFriend(pResultsChain, "pChn");
	}

	if(eResultsID != "")
	{
		loadEResultsChain();
		eventsChain->AddFriend(eResultsChain, "eChn");
	}

	if(gResultsID[0] == 'R' || eResultsID[0] == 'R' || pResultsID[0] == 'R')  //Only make Derv for recent run ids
	{
		loadDervChains();
		eventsChain->AddFriend(pDervChain, "pDChn");
		eventsChain->AddFriend(eDervChain, "eDChn");
	}

	if(isHmg)
	{
		loadLittleBChain();
		eventsChain->AddFriend(littleBChain, "bChn");
	}

	if(isRad && gResultsID != "")
	{
		loadGResultsChain();
		eventsChain->AddFriend(gResultsChain, "gChn");
		if(gResultsID[0] == 'R')  //Only make Derv for recent run ids
			eventsChain->AddFriend(gDervChain, "gDChn");
	}
	filesLoaded = true;

}

TH1* RDK2Set::makeHist(TString histName, TString drawString, TCut inpCut, HistDim inpHistDimX, HistDim inpHistDimY)
{
	HistDim inpHistDimZ = { 0, 0, 0 };
	return makeHist(histName, drawString, inpCut, inpHistDimX, inpHistDimY, inpHistDimZ);
}

TH1* RDK2Set::makeHist(TString histName, TString drawString, TCut inpCut, HistDim inpHistDimX)
{
	HistDim inpHistDimY = { 0, 0, 0 };
	HistDim inpHistDimZ = { 0, 0, 0 };
	return makeHist(histName, drawString, inpCut, inpHistDimX, inpHistDimY, inpHistDimZ);
}

TH1* RDK2Set::makeHist(TString histName, TString drawString, TCut inpCut, HistDim inpHistDimX, HistDim inpHistDimY, HistDim inpHistDimZ)
{
	gROOT->cd();
	TH1* histPrototype;
	TH1* outHist;

	if(inpHistDimZ.numBins > 0)
		histPrototype = makeTH3DFromDim("prototype2", "", inpHistDimX, inpHistDimY, inpHistDimZ);
	else if(inpHistDimY.numBins > 0)
		histPrototype = makeTH2DFromDim("prototype2", "", inpHistDimX, inpHistDimY);
	else
		histPrototype = makeTH1DFromDim("prototype2", "", inpHistDimX);

	outHist = makeHist(histName, drawString, inpCut, histPrototype);
	delete histPrototype;
	return outHist;

}

TH1* RDK2Set::makeHist(TString histName, TString drawString, TCut inpCut, TH1* histPrototype)
{
	gROOT->cd();
	TH1* outHist = (TH1*) histPrototype->Clone(histName);
	outHist->Reset();
	if(!filesLoaded)
	{
		loadChainsAndFriend();
	}
	TString fullDrawString = drawString + " >> " + histName;

	cout << "Drawing: " << fullDrawString << endl;
	cout << "using cuts: " << inpCut.GetTitle() << endl;
	cout << "Events drawn: " << eventsChain->Draw(fullDrawString, inpCut,"goff") << endl;
	cout << "Hist Mean: " << outHist->GetMean() << endl;
	outHist->Sumw2();
	return outHist;
}

TH2D* RDK2Set::makeHist2D(TString histName, TString drawString, TCut inpCut, HistDim inpHistDim1, HistDim inpHistDim2)
{
	if(!filesLoaded)
	{
		loadChainsAndFriend();
	}

	TH2D* outHist = new TH2D(histName, histName, inpHistDim1.numBins, inpHistDim1.binLowEdge, inpHistDim1.binHighEdge, inpHistDim2.numBins, inpHistDim2.binLowEdge, inpHistDim2.binHighEdge);

	cout << "Drawing: " << drawString << endl;
	cout << "using cuts: " << inpCut.GetTitle() << endl;
	eventsChain->Project(histName, drawString, inpCut,"goff");
	outHist->Sumw2();
	return outHist;
}

TH2D* RDK2Set::makeHist2D(TString histName, TString drawString, TCut inpCut, int numX, double* varArrayX, int numY, double* varArrayY)
{
	if(!filesLoaded)
	{
		loadChainsAndFriend();
	}

	TH2D* outHist = new TH2D(histName, histName, numX, varArrayX, numY, varArrayY);

	cout << "Drawing: " << drawString << endl;
	cout << "using cuts: " << inpCut.GetTitle() << endl;
	eventsChain->Project(histName, drawString, inpCut);
	outHist->Sumw2();
	return outHist;
}

void RDK2Set::makeDerivedResultsFiles(ParType parType)
{

	TChain* theChain = nullptr;

	TString letter;

	TString resultsID;

	int numEventsPerFile = 0;

	if(parType == PROTON)
	{

		loadPResultsChain();
		theChain = pResultsChain;
		letter = "P";
		resultsID = pResultsID;
		numEventsPerFile = pResultsChain->GetEntries() / GetNumFiles();
	}
	else if(parType == ELECTRON)
	{

		loadEResultsChain();
		theChain = eResultsChain;
		letter = "E";
		resultsID = eResultsID;
		numEventsPerFile = eResultsChain->GetEntries() / GetNumFiles();
	}
	else if(parType == PHOTON)
	{
		if(isHmg) return;

		loadGResultsChain();
		theChain = gResultsChain;
		letter = "G";
		resultsID = gResultsID;
		numEventsPerFile = gResultsChain->GetEntries() / GetNumFiles();

	}

	cout << "Making Derived Files for " << resultsID << " " << letter << endl;

	for (int j = 0; j < numFiles; j++)
	{
		TString derivedPath = RESULTS_DIR;
		if(isEventGenerator)
		{
			derivedPath += "Results_" + resultsID + "_" + letter + "_" + eventsTypeString + "_" + eventsDateString + "_Derv_Set" + int2str(j) + ".root";
		}
		else
		{
			derivedPath += "Results" + eventsDateString + "_" + eventsTypeString + "_" + letter + "_Set" + int2str(j) + "_" + resultsID + "_Derv.root";
		}

		createDerivedResultsFile(j * numEventsPerFile, theChain, derivedPath.Data());
	}

}

double RDK2Set::calcEPRate(TCut epCut)
{
	if(savedEPCut == epCut) return epRate;

	HistDim inpHistDim = { 100, 0, 1000 };

	TString drawString;
	if(eResultsID[0] != 'R')
		drawString = "eChn.SBDEDepTotal";
	else
		drawString = "eDChn.SBDEDepBlur";

	TCut totalCut;
	if(isHmg)
	{
		totalCut = epCut && STD_bZeroCut;
	}
	else
	{
		totalCut = epCut;
	}

	TH1D* electronEnergyHist = (TH1D*) makeHist("eeHist", drawString, totalCut, inpHistDim);

	epRate = electronEnergyHist->Integral();

	delete electronEnergyHist;

	if(isHmg)
	{
		HistDim libHistDim = { 3, -1, 1 };
		TH1D* lilbHist = (TH1D*) makeHist("lilbHist", "bChn.dwcutb20", STD_bZeroCut, libHistDim);
		epRate /= lilbHist->Integral();
		delete lilbHist;
	}
	else
	{
		epRate /= eResultsChain->GetEntries();
	}

	cout << "EP rate for " << GetName() << " is " << epRate << endl;

	savedEPCut = epCut;

	return epRate;
}

double RDK2Set::calcEPGRate(TCut epCut, TCut gCut)
{

	HistDim inpHistDim = { 100, 0, 1000 };

	TString drawString;
	if(eResultsID[0] != 'R')
		drawString = "eChn.SBDEDepTotal";
	else
		drawString = "eDChn.SBDEDepBlur";

	TCut totalCut = epCut && gCut;

	TH1D* electronEnergyHist = (TH1D*) makeHist("eeHist", drawString, totalCut, inpHistDim);

	double epgRate = electronEnergyHist->Integral();

	delete electronEnergyHist;

	epgRate /= eResultsChain->GetEntries();

	cout << "EPG rate for " << GetName() << " is " << epgRate << endl;

	return epgRate;
}

double RDK2Set::calcEPGRateFor12IndividualBGOs(TCut epCut)
{
	double epgRate = 0;
	for (int i = 0; i < 12; i++)
	{
		TString gCut = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam" + "> 10";
		epgRate += calcEPGRate(epCut, TCut(gCut));
	}
	cout << "Total EPG Rate for 12 detectors: " << epgRate << endl;
	return epgRate;
}

double RDK2Set::calcEPGRateFor3IndividualBAPDs(TCut epCut, double gammaEnergyCutLow, double gammaEnergyCutHigh)
{
	double epgRate = 0;
	for (int i = 0; i < 3; i++)
	{
		TString gCut = "gChn.BAPD" + int2str(i + 1) + "EDepTotal" + "> " + d2str(gammaEnergyCutLow) + " && gChn.BAPD" + int2str(i + 1) + "EDepTotal" + "< " + d2str(gammaEnergyCutHigh);
		epgRate += calcEPGRate(epCut, TCut(gCut));
	}
	cout << "Total EPG Rate for 3 BAPD detectors: " << epgRate << endl;
	return epgRate;
}

TH1* RDK2Set::makeEEHistDecayWidthProb(TString nameString, RDK2CutSet inpCutSet,double littleb, double normConst, HistDim histDim)
{
	TString detType =DET_EP;
	if(!filesLoaded)
		loadChainsAndFriend();

	TString HistsDir = HISTS_DIR;
	TString baseHistName = TString("MC_Hist_HmgProb_") + GetName();
	TString endHistName = ".txt";

	double pE, eE, pT;
	//Begin by using SBDEDepTotal for P energy as we once were
	pDervChain->SetBranchAddress("SBDEDepBlur", &pE);
	eDervChain->SetBranchAddress("SBDEDepBlur", &eE);
	pResultsChain->SetBranchAddress("SBDTimeFirst", &pT);

	double ee0,ep0,mxe0,mye0,mze0,mxp0,myp0,mzp0; // Original event variables
	eventsChain->SetBranchAddress("ee0",&ee0);
	eventsChain->SetBranchAddress("ep0",&ep0);
	eventsChain->SetBranchAddress("mxe0",&mxe0);
	eventsChain->SetBranchAddress("mye0",&mye0);
	eventsChain->SetBranchAddress("mze0",&mze0);
	eventsChain->SetBranchAddress("mxp0",&mxp0);
	eventsChain->SetBranchAddress("myp0",&myp0);
	eventsChain->SetBranchAddress("mzp0",&mzp0);


	TString coString = "EP";

	TH1* eEHist = makeTH1DFromDim(nameString, HistsDir + baseHistName + coString + "_EDepE" + endHistName, histDim);
	vector<double> binProb(histDim.numBins+2,0.);


	TVector3 eVec,pVec,nVec;
	double cos_en,prob;
	///Fill Hists
	int numEvents = eventsChain->GetEntries();
	int foundBin=-1;
	for (int i = 0; i < numEvents; i++)
	{
		eventsChain->GetEntry(i);

		if(inpCutSet.CheckEE(eE)) //Passes e cut
		{

			if(inpCutSet.CheckPT(pT) && inpCutSet.CheckPE(pE))
			{

				//Calc prob
				eVec.SetXYZ(mxe0,mye0,mze0);
				pVec.SetXYZ(mxp0,myp0,mzp0);
				nVec=-(eVec+pVec);
				cos_en = eVec.Dot(nVec)/(eVec.Mag()*nVec.Mag());
				prob = getJTWProb(ee0,ep0,cos_en,littleb)*normConst;

				foundBin=eEHist->FindBin(eE);
				if(foundBin >= 0 && foundBin < histDim.numBins+2)
				{
					binProb[foundBin] = binProb[foundBin] + prob;
					eEHist->Fill(eE);
				}

			}

		}

	}

	//Scale histogram by probabilies
	eEHist->Sumw2();
	for(int i=0;i<binProb.size();i++){
		eEHist->SetBinContent(i,eEHist->GetBinContent(i)*binProb[i]);
		eEHist->SetBinError(i,eEHist->GetBinError(i)*binProb[i]);
	}
	scaleHistTo(eEHist,1.);


	return eEHist;

}

double RDK2Set::makeAllStandardHists(TString nameString, double& error, RDK2CutSet inpCutSet, CoDet detType, double inpEPPer3Decay, double branchingRatio)
{
	loadChainsAndFriend();
	int numDets = 0;
	TString HistsDir = HISTS_DIR;
	TString baseHistName = TString("MC_Hist_") + GetName();
	TString endHistName = ".txt";

	double pE, eE, rE[12], pT;
	//Begin by using SBDEDepTotal for P energy as we once were
	pResultsChain->SetBranchAddress("SBDEDepTotal", &pE);
	eDervChain->SetBranchAddress("SBDEDepBlur", &eE);
	pResultsChain->SetBranchAddress("SBDTimeFirst", &pT);

	TString coString;
	TH1* radAvgHist = nullptr;
	TH1* radAvgVarHist = nullptr;
	TH1* radDetHists[12] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
	TH1* radDetVarHists[12] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
	if(detType == DET_EP)
	{
		coString = "EP";
	}
	else if(detType == DET_EPG)
	{
		//        numDets=12;
		numDets = 1;
		coString = "EPG";
		radAvgHist = makeTH1DFromDim("gEAvgHist", HistsDir + baseHistName + coString + "_EDepGAvg" + endHistName, BGO_DIM);
		radAvgVarHist = new TH1D("gEAvgVarHist", HistsDir + TString("MC_VarHist_") + GetName() + coString + "_EDepGAvg" + endHistName, 20, BGO_VARHIST_LOWEDGES);
		for (int i = 0; i < numDets; i++)
		{
			radDetHists[i] = makeTH1DFromDim("gEDet" + int2str(i + 1) + "Hist", HistsDir + baseHistName + coString + "_EDepGDet" + int2str(i + 1) + endHistName, BGO_DIM);
			radDetVarHists[i] = new TH1D("gEDet" + int2str(i + 1) + "VarHist", HistsDir + TString("MC_VarHist_") + GetName() + coString + "_EDepGDet" + int2str(i + 1) + endHistName, 20, BGO_VARHIST_LOWEDGES);
			gDervChain->SetBranchAddress("BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam", &(rE[i]));
		}
	}
	else if(detType == DET_EPBG)
	{
		//        numDets=3;
		numDets = 1;
		coString = "EPBG";
		radAvgHist = makeTH1DFromDim("bGEAvgHist", HistsDir + baseHistName + coString + "_EDepBGAvg" + endHistName, BAPD_DIM);
		radAvgVarHist = new TH1D("bGEAvgVarHist", HistsDir + TString("MC_VarHist_") + GetName() + coString + "_EDepBGAvg" + endHistName, 20, BGO_VARHIST_LOWEDGES);  //Need to fix this dimension later
		for (int i = 0; i < numDets; i++)
		{
			radDetHists[i] = makeTH1DFromDim("bGEDet" + int2str(i + 1) + "Hist", HistsDir + baseHistName + coString + "_EDepBGDet" + int2str(i + 1) + endHistName, BAPD_DIM);
			radDetVarHists[i] = new TH1D("bGEDet" + int2str(i + 1) + "VarHist", HistsDir + TString("MC_VarHist_") + GetName() + coString + "_EDepBGDet" + int2str(i + 1) + endHistName, 20, BGO_VARHIST_LOWEDGES);  //Need to fix this dimension later
			gResultsChain->SetBranchAddress("BAPD" + int2str(i + 1) + "EDepTotal", &(rE[i]));
		}
	}

	TH1* eEHist = makeTH1DFromDim("eEHist", HistsDir + baseHistName + coString + "_EDepE" + endHistName, EDEPE_DIM);
	TH1* pEHist = makeTH1DFromDim("pEHist", HistsDir + baseHistName + coString + "_EDepP" + endHistName, EDEPP_DIM);
	TH1* pTHist = makeTH1DFromDim("pTHist", HistsDir + baseHistName + coString + "_ToFP" + endHistName, TOFP_DIM2);

	///Fill Hists
	int numEvents = pResultsChain->GetEntries();
	int numHits = 0;
	for (int i = 0; i < numEvents; i++)
	{
		eDervChain->GetEntry(i);
		if(inpCutSet.CheckEE(eE)) //Passes e cut
		{
			pResultsChain->GetEntry(i);
			if(inpCutSet.CheckPT(pT) && inpCutSet.CheckPE(pE))
			{
				if(detType == DET_EP)
				{
					numHits++;
					eEHist->Fill(eE);
					pEHist->Fill(pE);
					pTHist->Fill(pT);
				}
				else if(detType == DET_EPG)
				{
					gDervChain->GetEntry(i);
					for (int j = 0; j < numDets; j++)
					{
						if(inpCutSet.CheckGE(rE[j]))
						{
							numHits++;
							eEHist->Fill(eE);
							pEHist->Fill(pE);
							pTHist->Fill(pT);
							radAvgHist->Fill(rE[j]);
							radAvgVarHist->Fill(rE[j]);
							radDetHists[j]->Fill(rE[j]);
							radDetVarHists[j]->Fill(rE[j]);
						}
					}

				}
				else if(detType == DET_EPBG)
				{
					gResultsChain->GetEntry(i);
					for (int j = 0; j < numDets; j++)
					{
						if(inpCutSet.CheckBGE(rE[j]))
						{
							numHits++;
							eEHist->Fill(eE);
							pEHist->Fill(pE);
							pTHist->Fill(pT);
							radAvgHist->Fill(rE[j]);
							radAvgVarHist->Fill(rE[j]);
							radDetHists[j]->Fill(rE[j]);
							radDetVarHists[j]->Fill(rE[j]);
						}
					}

				}
			}

		}

	}

	///Make hists ratio hist as necessary
	if(detType == DET_EPG || detType == DET_EPBG)
	{
		double scaleFactor = branchingRatio / (inpEPPer3Decay * numEvents);
		eEHist->Sumw2();
		eEHist->Scale(scaleFactor / numDets);
		pEHist->Sumw2();
		pEHist->Scale(scaleFactor / numDets);
		pTHist->Sumw2();
		pTHist->Scale(scaleFactor / numDets);
		radAvgHist->Sumw2();
		radAvgHist->Scale(scaleFactor / numDets);
		radAvgVarHist->Sumw2();
		radAvgVarHist->Scale(scaleFactor / numDets);
		convertTH1ToTXT(radAvgHist, radAvgHist->GetTitle());
		convertTH1ToTXT(radAvgVarHist, radAvgVarHist->GetTitle());
		for (int i = 0; i < numDets; i++)
		{
			radDetHists[i]->Sumw2();
			radDetHists[i]->Scale(scaleFactor);
			convertTH1ToTXT(radDetHists[i], radDetHists[i]->GetTitle());
			radDetVarHists[i]->Sumw2();
			radDetVarHists[i]->Scale(scaleFactor);
			convertTH1ToTXT(radDetVarHists[i], radDetVarHists[i]->GetTitle());
			delete radDetHists[i];
			delete radDetVarHists[i];
		}
		delete radAvgHist;
		delete radAvgVarHist;

	}

	convertTH1ToTXT(eEHist, eEHist->GetTitle());
	convertTH1ToTXT(pEHist, pEHist->GetTitle());
	convertTH1ToTXT(pTHist, pTHist->GetTitle());

	delete eEHist;
	delete pEHist;
	delete pTHist;

	double answer = numHits / (double) numEvents;
	error = answer / sqrt(numHits);
	return answer;
}

void RDK2Set::setEventList(TCut inpCut)
{
	if(!filesLoaded)
	{
		loadChainsAndFriend();
	}
	eventsChain->SetEventList(0); //reset
	eventsChain->Draw(">>theEventList", inpCut);
	TEventList* theEventList;
	gDirectory->GetObject("theEventList", theEventList);
	eventsChain->SetEventList(theEventList);

}

TH1D* makeTH1DFromDim(TString name, TString title, HistDim inpHistDim)
{
	return new TH1D(name, title, inpHistDim.numBins, inpHistDim.binLowEdge, inpHistDim.binHighEdge);
}

TH2D* makeTH2DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim)
{
	return new TH2D(name, title, inpXHistDim.numBins, inpXHistDim.binLowEdge, inpXHistDim.binHighEdge, inpYHistDim.numBins, inpYHistDim.binLowEdge, inpYHistDim.binHighEdge);
}

TH3D* makeTH3DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim, HistDim inpZHistDim)
{
	return new TH3D(name, title, inpXHistDim.numBins, inpXHistDim.binLowEdge, inpXHistDim.binHighEdge, inpYHistDim.numBins, inpYHistDim.binLowEdge, inpYHistDim.binHighEdge, inpZHistDim.numBins, inpZHistDim.binLowEdge, inpZHistDim.binHighEdge);
}

HistDim getDimFromHist(TH1* inpHist)
{
	HistDim outDim;
	outDim.numBins = inpHist->GetNbinsX();
	outDim.binLowEdge = inpHist->GetBinLowEdge(1);
	outDim.binHighEdge = inpHist->GetBinWidth(outDim.numBins) + inpHist->GetBinLowEdge(outDim.numBins);
	return outDim;
}

int getVarArrayFromHist(TH1* inpHist, double** outArray)
{
	const int numBins = inpHist->GetNbinsX();
	*outArray = new double[numBins + 1];

	for (int i = 0; i < numBins; i++)
	{
		(*outArray)[i] = inpHist->GetBinLowEdge(i + 1);

	}
	(*outArray)[numBins] = inpHist->GetBinLowEdge(numBins) + inpHist->GetBinWidth(numBins);

	return numBins;
}

