#include "RDK2Fierz.h"

#include <iostream>
#include <cmath>

#include "TString.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TH1.h"
#include "TCanvas.h"

#include "MRIO.h"
#include "MRMisc.h"
#include "MRGraphics.h"

#include "RDK2Constants.h"
#include "RDK2Events.h"
#include "RDK2IO.h"
#include "RDK2MCAnalyzer.h"

void createHomogenousEventSettings(int numRandomSeeds, TString eventSettingsString)
{
	TString filePath = MRKGEANTDATA_DIR + TString("3BRandomNumbers.txt");
	ifstream randomNumbersFile(filePath);
	MRKEvents theEvents;

	TString ranFilePath = TString(MRKGEANTDATA_DIR) + "EventSettings_Hmg_3B_" + eventSettingsString + ".root";
	TFile ranGenFile(ranFilePath, "recreate");
	ranGenFile.cd();
	TVectorD fourBody(1);
	TVectorD homogeneousEvents(1);
	TVectorD littleb(1);
	TVectorD normalizationConstant(1);
	TVectorD minimumGammaEnergy(1);
	TVectorD maximumGammaEnergy(1);

	fourBody[0] = 0;
	homogeneousEvents[0] = 1;
	littleb[0] = 0;
	normalizationConstant[0] = 0;

	minimumGammaEnergy[0] = 0;
	maximumGammaEnergy[0] = 0;

	fourBody.Write("fourBody", TObject::kOverwrite);
	homogeneousEvents.Write("homogeneousEvents", TObject::kOverwrite);
	littleb.Write("littleb", TObject::kOverwrite);
	normalizationConstant.Write("normalizationConstant", TObject::kOverwrite);
	minimumGammaEnergy.Write("minimumGammaEnergy", TObject::kOverwrite);
	maximumGammaEnergy.Write("maximumGammaEnergy", TObject::kOverwrite);
	for (int i = 0; i < numRandomSeeds; i++)
	{
		int randomSeed;
		randomNumbersFile >> randomSeed;
		TRandom3* ranGen = new TRandom3(randomSeed);
		TString ranGenName = "ranGen" + int2str(i);
		cout << ranGenName << endl;
		ranGen->SetName(ranGenName);
		ranGen->Write("", TObject::kOverwrite);
		ranGenFile.Write("", TObject::kOverwrite);

	}
	ranGenFile.Close();
	randomNumbersFile.close();
}

void createHomogenousEventsUsingSettings(int numFiles, int numEvents, TString eventSettingsString)
{
	TString fluxFilePath = TString(MRKGEANTDATA_DIR) + "Fluxmap_RDK2_Matt1.txt";
	TString eventSettingsPath = TString(MRKGEANTDATA_DIR) + "EventSettings_Hmg_3B_" + eventSettingsString + ".root";
	for (int i = 0; i < numFiles; i++)
	{
		MRKEvents theEvents;
		TString eventFilePath = TString(MRKGEANTDATA_DIR) + "Events_Hmg_3B_" + eventSettingsString + "_Set" + TString::Itoa(i, 10) + ".root";
		theEvents.loadEventSettingsAndMakeFile(eventSettingsPath, i, eventFilePath, fluxFilePath, numEvents);
	}
}

void createLilBTreesFromHmgEventFiles(int numFiles, TString eventSettingsString, int numLilB, double startLilB, double incrementLilB)
{
	for (int i = 0; i < numFiles; i++)
	{
		TString eventFilePath = TString(EVENTS_DIR) + "Events_" + eventSettingsString + "_Set" + TString::Itoa(i, 10) + ".root";
		if(!FileExists(eventFilePath))
		{
			cout << eventFilePath << " does not exist!" << endl;
			return;
		}
		MRKEvents theEvents;
		theEvents.loadEvents(eventFilePath, DEFAULT_EVENTS_TREENAME);
		theEvents.setNormConstant(6.0e-12*0.7);
		theEvents.makeDecayWidthCutMultiple(numLilB, startLilB, incrementLilB);
	}
}

void makeDWCutHists(int numDWCut, TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet)
{
	RDK2MCAnalyzer theMCAnalyzer(inpName, inpTitle, num3, p3RID, e3RID, num4, p4RID, e4RID, g4RID, inpCutSet, true);
	for (int i = 0; i < numDWCut; i++)
	{
		theMCAnalyzer.SetName(inpName + "_dwcut" + int2str(i));
		theMCAnalyzer.SetAddToAllCut(TCut("bChn.dwcutb" + int2str(i) + TString(" == 0")));
		theMCAnalyzer.MakeAllStandardHists();
	}
}

void plotLilBExpFitToMC(TString expHistFilePath, TString mcAnalIDString, HistDim histDim, int numDWCut, double startLilB, double incrementLilB, double beginEnergy, double endEnergy)
{
	int numRebin = 100;
	if(!FileExists(expHistFilePath))
	{
		cout << "Could not find: " << expHistFilePath << endl;
		return;
	}
	TH1* expHist = getTabSeperatedHist(expHistFilePath, histDim);
	expHist->SetName("expHist");
	expHist->SetTitle("Exp.;Electron Pulse Height (~keV);Experimental counts");
	expHist->Rebin(numRebin);
	expHist->GetXaxis()->SetRangeUser(beginEnergy, endEnergy);

	double res[histDim.numBins] = { };

	double littleBValues[histDim.numBins] = { };
	double chiSquaredTestResults[histDim.numBins] = { };
	double kolmogorovTestResults[histDim.numBins] = { };

	for (int i = 0; i < numDWCut; i++)
	{
		littleBValues[i] = startLilB + i * incrementLilB;

		int beginBin = lround((beginEnergy - histDim.binLowEdge) / (histDim.getBinWidth() * numRebin)) + 1;
		int endBin = lround((endEnergy - histDim.binLowEdge) / (histDim.getBinWidth() * numRebin)) + 1;

		TString mcHistFilePath = TString(HISTS_DIR) + "Hist_MC_" + mcAnalIDString + "_dwcut" + int2str(i) + "_EP_EDepE.txt";
		if(!FileExists(mcHistFilePath))
		{
			cout << "Could not find: " << mcHistFilePath << endl;
			return;
		}
		TH1* mcHist = getTabSeperatedHist(mcHistFilePath, histDim);
		mcHist->SetName("mcHist");
		mcHist->SetTitle("Hmg MC (scaled);Electron Pulse Height (~keV);Experimental counts");
		mcHist->Rebin(numRebin);
		scaleHistogramsTogether(expHist, mcHist, beginBin, endBin);

		mcHistFilePath = TString(HISTS_DIR) + "Hist_MC_AID87_EP_EDepE.txt";
		TH1* mcOldHist = getTabSeperatedHist(mcHistFilePath, histDim);
		mcOldHist->SetName("mcOldHist");
		mcOldHist->SetTitle("PRL MC (scaled);Electron Pulse Height (~keV);Experimental counts");
		mcOldHist->Rebin(numRebin);
		scaleHistogramsTogether(expHist, mcOldHist, beginBin, endBin);

		mcOldHist->GetXaxis()->SetRangeUser(beginEnergy, endEnergy);
		mcHist->GetXaxis()->SetRangeUser(beginEnergy, endEnergy);

		//Plot histograms together
		TH1* expHists[1] = { expHist };
		TH1* mcHists[2] = { mcOldHist, mcHist };
		TString imagePath = TString(GRAPHS_DIR) + "Graph_Exp_MC" + mcAnalIDString + "_dwcut" + int2str(i) + "_EP_EDepE.png";
		plotExpVersusMCToImage(1, expHists, 2, mcHists, "Comparison with little b = " + d2str(littleBValues[i]), imagePath);

		int ndf, igood;
		double rootChiSquaredResult;
		chiSquaredTestResults[i] = expHist->Chi2TestX(mcHist, rootChiSquaredResult, ndf, igood, "WW", res);
		cout << "P Result: " << chiSquaredTestResults[i] << endl;

		kolmogorovTestResults[i] = expHist->KolmogorovTest(mcHist);

		delete mcHist;
	}

	TGraph chiSquaredGraph(numDWCut, littleBValues, chiSquaredTestResults);
	chiSquaredGraph.SetName("chiSquaredGraph");

	TGraph kolmogorovGraph(numDWCut, littleBValues, kolmogorovTestResults);
	kolmogorovGraph.SetName("kolmogorovGraph");

	TCanvas theCanvas;

	chiSquaredGraph.Draw("ALP");

	TString imagePath = TString(GRAPHS_DIR) + "Graph_ChiSquaredTest_" + mcAnalIDString + ".png";
	theCanvas.SaveAs(imagePath);
	delete expHist;
}
