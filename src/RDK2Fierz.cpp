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
	RDK2Events theEvents;

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
		RDK2Events theEvents;
		TString eventFilePath = TString(MRKGEANTDATA_DIR) + "Events_Hmg_3B_" + eventSettingsString + "_Set" + TString::Itoa(i, 10) + ".root";
		theEvents.loadEventSettingsAndMakeFile(eventSettingsPath, i, eventFilePath, fluxFilePath, numEvents);
	}
}

void createLilBTreesFromHmgEventFiles(int numFiles, TString eventSettingsString, int numLilB, double startLilB, double incrementLilB, int littleBSet)
{
	for (int i = 0; i < numFiles; i++)
	{
		TString eventFilePath = TString(EVENTS_DIR) + "Events_" + eventSettingsString + "_Set" + TString::Itoa(i, 10) + ".root";
		if(!FileExists(eventFilePath))
		{
			cout << eventFilePath << " does not exist!" << endl;
			return;
		}
		RDK2Events theEvents;
		theEvents.loadEvents(eventFilePath, DEFAULT_EVENTS_TREENAME);
		theEvents.setNormConstant(6.0e-12 * 0.7);
		theEvents.makeDecayWidthCutMultiple(numLilB, startLilB, incrementLilB, littleBSet);
	}
}

void makeEE0Plots(int numDWCut, TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet)
{
	RDK2MCAnalyzer theMCAnalyzer(inpName, inpTitle, num3, p3RID, e3RID, num4, p4RID, e4RID, g4RID, inpCutSet, true);
	TH1* theHists[numDWCut];
	for (int i = 0; i < numDWCut; i++)
	{
		theMCAnalyzer.SetName(inpName + "_dwcut" + int2str(i));
		theMCAnalyzer.SetAddToAllCut(TCut("bChn.dwcutb" + int2str(i) + TString(" == 0")));
		HistDim theDim = { 100, 0, 800 };
		theHists[i] = theMCAnalyzer.MakeHist("ee0Hist", THREEBODY, "ee0", "", theDim);
		cout << "Num Events: " << theHists[i]->Integral() << endl;
	}

	for (int i = 0; i < numDWCut; i++)
	{
		double littleBValue = -.1 + i * .01;
		TH1* drawnHists[2] = { theHists[10], theHists[i] };
		TString imagePath = TString(GRAPHS_DIR) + "Graph_Exp_MC_" + inpName + "_dwcut" + int2str(i) + "_EE0.png";
		plotExpVersusMCToImage(0, nullptr, 2, drawnHists, "Comparison with little b = " + d2str(littleBValue), imagePath);
	}
}

void makeDWCutHists(int numDWCut, TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4, TString p4RID, TString e4RID, TString g4RID, RDK2CutSet inpCutSet)
{
	RDK2MCAnalyzer theMCAnalyzer(inpName, inpTitle, num3, p3RID, e3RID, num4, p4RID, e4RID, g4RID, inpCutSet, true);

	for (int i = 0; i < numDWCut; i++)
	{
		cout << "working on:" << i << endl;
		theMCAnalyzer.SetAddToAllCut(TCut("bChn.dwcutb" + int2str(i) + TString(" == 0")));
		theMCAnalyzer.SetName(inpName + "_dwcut" + int2str(i));
//		theMCAnalyzer.MakeAllStandardHists();
		TH1* eEHist = theMCAnalyzer.MakeStandardHist("ep_eEHist", DET_EP, "eDChn.SBDEDepBlur", EDEPE_DIM2);
		eEHist->SetTitle("Electron Energy EP;SBD Energy Channel (782=782 keV);Counts");

		convertTH1ToTXT(eEHist, HISTS_DIR + TString("Hist_MC_") + theMCAnalyzer.GetName() + "_EP_EDepE.txt");
		delete eEHist;
	}
}

void plotLilBExpFitToMC(TString plotTitleString,TString expHistFilePath, TString mcAnalIDString, HistDim histDim, int numDWCut, double startLilB, double incrementLilB, double beginEnergy, double endEnergy)
{
	int numRebin = 100;
	if(!FileExists(expHistFilePath))
	{
		cout << "Could not find: " << expHistFilePath << endl;
		return;
	}
	TH1* expHist = getTabSeperatedHist(expHistFilePath, histDim);
	expHist->SetName("expHist");

	expHist->Rebin(numRebin);
	expHist->GetXaxis()->SetRangeUser(beginEnergy, endEnergy);

//	Set bin error to constant percentage
	for (int i = 0; i < expHist->GetNbinsX(); i++)
	{
		expHist->SetBinError(i + 1, expHist->GetBinContent(11) * 0.01);
	}

	double res[histDim.numBins];
	fill_n(res, histDim.numBins, 0);

	double littleBValues[histDim.numBins];
	fill_n(littleBValues, histDim.numBins, 0);
	double chiSquaredTestResults[histDim.numBins];
	fill_n(chiSquaredTestResults, histDim.numBins, 0);

//	TH1* mcHists[numDWCut]={};

	for (int i = 0; i < numDWCut; i++)
	{
		littleBValues[i] = startLilB + i * incrementLilB;

		int beginBin = lround((beginEnergy - histDim.binLowEdge) / (histDim.getBinWidth() * numRebin)) + 1.;
		int endBin = lround((endEnergy - histDim.binLowEdge) / (histDim.getBinWidth() * numRebin)) + 1.;

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
		expHist->SetTitle("Exp.;Electron Pulse Height (~keV);Experimental counts");
		TString imagePath = TString(GRAPHS_DIR) + "Graph_Exp_MC" + mcAnalIDString + "_dwcut" + int2str(i) + "_EP_EDepE.png";
		plotExpVersusMCToImage(1, expHists, 2, mcHists, "Comparison with little b = " + d2str(littleBValues[i]), imagePath);

		int ndf, igood;
		double rootChiSquaredResult;
		expHist->Chi2TestX(mcHist, rootChiSquaredResult, ndf, igood, "WW", res);
		chiSquaredTestResults[i] = rootChiSquaredResult;
		cout << "ChiSquared Result: " << chiSquaredTestResults[i] << endl;

		delete mcHist;
		delete mcOldHist;
	}

	TGraph chiSquaredGraph(numDWCut, littleBValues, chiSquaredTestResults);
	chiSquaredGraph.SetName("chiSquaredGraph");

	TCanvas theCanvas;
	gPad->SetLogy(1);

	chiSquaredGraph.SetTitle(plotTitleString+";b;#tilde{#chi}^{2}");

	chiSquaredGraph.Draw("ALP");

	TString imagePath = TString(GRAPHS_DIR) + "Graph_ChiSquaredTest_" + mcAnalIDString + ".png";
	theCanvas.SaveAs(imagePath);
	delete expHist;
}

void analyzeElectronBackscatteringSBDOnly(TString resultsFilePath)
{
	//Read in results file and tree
	TString filePath = resultsFilePath;
	cout << "Analyzing electron spectrum from: " << filePath << endl;

	TFile resultsFile(filePath, "READ");

	TTree* resultsTree = (TTree*) resultsFile.Get("g");

	TH1D* eHist = makeTH1DFromDim("eHist", "eHist", EDEPE_DIM);
	resultsTree->Draw("SBDEDepTotal >> eHist", "SBDEDepTotal > 0");

	TString histFilePath = TString(HISTS_DIR) + "Hist_E_" + fileNameNoExtFromFullPath(resultsFilePath) + ".txt";

	convertTH1ToTXT(eHist, histFilePath);

	resultsFile.Close();

}

void analyzeAllEBackscatteringSBDOnly()
{
	vector<TString> angleStrings = { "90", "75", "60", "45", "30", "15" };
	vector<TString> verStrings = { "4.9.6.p02", "4.9.6.p03", "4.10.02" };
	vector<TString> libStrings = { "GS","Penelope","LivStd", "WVI","SingleScatter"};
	vector<TString> energyStrings = { "100", "200", "400","600","800"};

	for (auto angleStr : angleStrings)
	{
		for (auto verStr : verStrings)
		{
			for (auto energyStr : energyStrings)
			{
				for (auto libStr : libStrings)
				{
					TString filePath = TString(MRKGEANTDATA_DIR) + "sbdonly/SBDOnly_E_" + verStr + "_"+libStr+"_" + energyStr + "keV_" + angleStr + "deg.root";
					if(FileExists(filePath) && !TFile(filePath).IsZombie())
					{
						analyzeElectronBackscatteringSBDOnly(filePath);
					}
				}
			}
		}
	}
}

void makeBSPlotsSBDOnly()
{
//	vector<TString> angleStrings = { "90", "75", "60", "45", "30", "15" };
	vector<TString> angleStrings = { "90", "45" };
//	vector<TString> verStrings = { "4.9.6.p03", "4.9.6.p03", "4.10.02" };
	vector<TString> verStrings = { "4.10.02" };
	vector<TString> libStrings = { "GS","Penelope","LivStd", "WVI","SingleScatter"};
	vector<TString> energyStrings = { "100", "200", "400","600" };
	const size_t numLibs = libStrings.size();
	for (TString angleStr : angleStrings)
	{
		for (TString energyStr : energyStrings)
		{
			for (TString verStr : verStrings)
			{
				TH1* theHists[numLibs];
				for (size_t i = 0; i < numLibs; i++)
				{
					TString filePath = TString(HISTS_DIR) + "Hist_E_SBDOnly_E_" + verStr + "_"+libStrings[i]+"_" + energyStr + "keV_" + angleStr + "deg.txt";
					double ratio=0;
					if(FileExists(filePath))
					{
						theHists[i] = getTabSeperatedHist(filePath, EDEPE_DIM);
						ratio=theHists[i]->Integral(101,801)*100. /1000000.;
					}
					else
					{
						theHists[i]=makeTH1DFromDim("aHist"+TString::Itoa(i,10),"aHist",EDEPE_DIM);
					}
					theHists[i]->SetTitle(libStrings[i]+": "+d2str(ratio,1,1)+"% >100 keV;Energy Absorbed (keV);Num. Events");
					theHists[i]->Rebin(10);
					theHists[i]->GetXaxis()->SetRangeUser(0,energyStr.Atof());
				}
				TString graphTitleString = TString::Format("Angle from Surface: %s deg, Energy: %s keV", angleStr.Data(), energyStr.Data());
				TString graphFilePath = TString(GRAPHS_DIR) + "Hist_E_SBDOnly_v4.10.02_" + energyStr + "keV_" + angleStr + "deg.pdf";
				plotExpVersusMCToImage(0, NULL, numLibs, theHists, graphTitleString, graphFilePath,true);
				for (size_t i = 0; i < numLibs; i++)
				{
					delete theHists[i];
				}
			}

		}
	}

}

void makeBSPlotsSBDOnly2()
{
	vector<TString> angleStrings = { "90", "75", "60", "45", "30", "15" };
	vector<TString> verStrings = { "4.9.6.p02", "4.9.6.p03", "4.10.02" };
	vector<TString> libStrings = { "LivStd"};
	vector<TString> energyStrings = { "100", "200", "400", "600","800"};
	const size_t numVer = verStrings.size();
	for (TString angleStr : angleStrings)
	{
		for (TString energyStr : energyStrings)
		{
			for (TString libStr : libStrings)
			{
				TH1* theHists[numVer];
				for (size_t i = 0; i < verStrings.size(); i++)
				{
					TString filePath = TString(HISTS_DIR) + "Hist_E_SBDOnly_E_" + verStrings[i] + "_"+libStr+"_" + energyStr + "keV_" + angleStr + "deg.txt";
					theHists[i] = getTabSeperatedHist(filePath, EDEPE_DIM);
					double ratio=theHists[i]->Integral(101,801)*100. /1000000.;
					theHists[i]->SetTitle(verStrings[i]+": "+d2str(ratio,1,1)+"% >100 keV;Energy Absorbed (keV);Num. Events");
					theHists[i]->Rebin(10);
					theHists[i]->GetXaxis()->SetRangeUser(0,energyStr.Atof());
				}
				TString graphTitleString = TString::Format("Angle from Surface: %s deg, Energy: %s keV", angleStr.Data(), energyStr.Data());
				TString graphFilePath = TString(GRAPHS_DIR) + "Hist_E_SBDOnly_"+libStr+"_" + energyStr + "keV_" + angleStr + "deg.pdf";
				plotExpVersusMCToImage(0, NULL, numVer, theHists, graphTitleString, graphFilePath,true);
				for (size_t i = 0; i < numVer; i++)
				{
					delete theHists[i];
				}
			}

		}
	}

}

void makeIncidentEventsFromResults(TString resultsFilePath,TString eventFilePath)
{
	if(!FileExists(resultsFilePath))
	{
		cout << "File: " << resultsFilePath << " does not exist." << endl;
		return;
	}
	TFile resultsFile(resultsFilePath,"READ");

	TTree* resultsTree=(TTree*) resultsFile.Get(DEFAULT_RESULTS_TREENAME);

	double SBDKEFirst,SBDMomXFirst,SBDMomYFirst,SBDMomZFirst;

	resultsTree->SetBranchAddress("SBDKEFirst",&SBDKEFirst);
	resultsTree->SetBranchAddress("SBDMomXFirst",&SBDMomXFirst);
	resultsTree->SetBranchAddress("SBDMomYFirst",&SBDMomYFirst);
	resultsTree->SetBranchAddress("SBDMomZFirst",&SBDMomZFirst);

	TFile eventFile(eventFilePath,"RECREATE");

	if(eventFile.IsZombie()) {
		cout << "Error finding event file: " << eventFilePath << endl;
		return;

	}
	eventFile.cd();

	TTree* eventTree=new TTree(DEFAULT_EVENTS_TREENAME,"Event Tree from results file");
	int n;
	double x0=0.;
	double y0=0.;
	double z0=0.0002;
	double ee0=0.;
	double ep0=0.;
	double eg0=0.;
	double mxe0=0.;
	double mye0=0.;
	double mze0=0.;
	double mxp0=0.;
	double myp0=0.;
	double mzp0=0.;
	double mxg0=0.;
	double myg0=0.;
	double mzg0=0.;

	eventTree->Branch("n",&n,"n/I");
	eventTree->Branch("x0",&x0,"x0/D");
	eventTree->Branch("y0",&y0,"y0/D");
	eventTree->Branch("z0",&z0,"z0/D");

	eventTree->Branch("ee0",&ee0,"ee0/D");
	eventTree->Branch("ep0",&ep0,"ep0/D");
	eventTree->Branch("eg0",&eg0,"eg0/D");

	eventTree->Branch("mxe0",&mxe0,"mxe0/D");
	eventTree->Branch("mye0",&mye0,"mye0/D");
	eventTree->Branch("mze0",&mze0,"mze0/D");

	eventTree->Branch("mxp0",&mxp0,"mxp0/D");
	eventTree->Branch("myp0",&myp0,"myp0/D");
	eventTree->Branch("mzp0",&mzp0,"mzp0/D");

	eventTree->Branch("mxg0",&mxg0,"mxg0/D");
	eventTree->Branch("myg0",&myg0,"myg0/D");
	eventTree->Branch("mzg0",&mzg0,"mzg0/D");

	double norm=1;
	for(int i=0;i<resultsTree->GetEntries();++i)
	{
		resultsTree->GetEntry(i);
		n=i;
		if(SBDKEFirst>100)
		{
			ee0=SBDKEFirst;
			norm=sqrt(SBDMomXFirst*SBDMomXFirst+SBDMomYFirst*SBDMomYFirst+SBDMomZFirst*SBDMomZFirst);
			mxe0=SBDMomXFirst*COS_BEND/norm-SBDMomZFirst*SIN_BEND/norm;
			mye0=SBDMomYFirst/norm;
			mze0=SBDMomZFirst*COS_BEND/norm+SBDMomXFirst*SIN_BEND/norm;
		}
		else
		{
			ee0=0;
			mxe0=0;
			mye0=0;
			mze0=0;
		}
		eventTree->Fill();
	}
	eventTree->Write("",TObject::kOverwrite);
	eventFile.Close();
	resultsFile.Close();
}

void makeAllIncidentEventsFromResults(int numFiles, TString RIDString)
{
	for(int i=0;i<numFiles;i++)
	{
		TString resultsFilePath=TString(RESULTS_DIR)+"Results_"+RIDString+"_E_Hmg_3B_160323_Set"+TString::Itoa(i,10)+".root";
		TString eventsFilePath=TString("/media/mjbales/mjbexternal/rdk/secondaryevents/")+"Events_Secondary"+RIDString+"_Hmg_3B_160323_Set"+TString::Itoa(i,10)+".root";
		cout << "Making " << eventsFilePath << " .... ";
		makeIncidentEventsFromResults(resultsFilePath,eventsFilePath);
		cout << " done!" << endl;
	}
}

void makeSecondaryMacros(int numFiles, TString RIDInString, TString RIDOutString,TString libString)
{
	ofstream macroFile("/home/mjbales/work/rdk/mrkgeantrun/sbdonly/secondary/"+RIDOutString+".mac");
	macroFile << "/run/verbose 0" << endl;
	macroFile << "/event/verbose 0" << endl;
	macroFile << "/tracking/verbose 0" << endl;
	macroFile << "#" << endl;
	if(libString=="GS")
	{
		macroFile << "/MRK/phys/addPhysics emstandardGS" << endl;
	}
	else if(libString=="SingleScatter")
	{
		macroFile << "/MRK/phys/addPhysics emstandardSS" << endl;
	}
	else if(libString=="WVI")
	{
		macroFile << "/MRK/phys/addPhysics emstandardWVI" << endl;
	}
	else if(libString=="Penelope")
	{
		macroFile << "/MRK/phys/addPhysics empenelope" << endl;
	}
	macroFile << "/MRK/field/addNamedFields None" << endl;
	macroFile << "/MRK/det/setmodel SBDOnly" << endl;
	macroFile << "/MRK/det/setSiDeadLayerLengthCmd 20 nm" << endl;
	macroFile << "/process/eLoss/StepFunction .2 .1 mm" << endl;
	macroFile << "/process/eLoss/minKinEnergy 10 eV" << endl;
	macroFile << "#/process/em/lowestMuHadEnergy 0.1 eV" << endl;
	macroFile << "/process/msc/RangeFactor .01" << endl;
	macroFile << "/process/msc/Skin 3" << endl;
	macroFile << "/MRK/phys/setCuts 1 nm" << endl;
	macroFile << "/cuts/setMaxCutEnergy 99 eV" << endl;
	macroFile << "/MRK/file/genMode 1" << endl;
	macroFile << "/MRK/field/update" << endl;
	macroFile << "/run/initialize" << endl;
	macroFile << "/gun/particle e-" << endl;
	macroFile << "#" << endl;
	macroFile << "######################################################" << endl;
	macroFile << "#";
	for(int i=0;i<100;i++)
	{
		macroFile << endl << "/MRK/file/inpfile /media/mjbales/mjbexternal/rdk/secondaryevents/Events_Secondary"+RIDInString+"_Hmg_3B_160323_Set"+TString::Itoa(i,10)+".root" << endl;
		macroFile << "/MRK/file/rootout Results_"+RIDOutString+"_E_Hmg_3B_160323_Set"+TString::Itoa(i,10)+".root" << endl;
		macroFile << "/run/beamOn 1000000";
	}
	macroFile.close();
}


double fitLittleB(TString nameString, TH1* expHist, TH1** mcHists,  int numLilB, double startLilB, double endLilB, double beginEnergy, double endEnergy,bool saveHistsAndPlots)
{


	int beginBin = lround((beginEnergy - expHist->GetBinLowEdge(1)) / (expHist->GetBinWidth(1) )) + 1.;
	int endBin = lround((endEnergy - expHist->GetBinLowEdge(1)) / (expHist->GetBinWidth(1) )) + 1.;

	int numBins=expHist->GetNbinsX();
	double res[numBins];
	fill_n(res, numBins, 0);

	double littleBValues[numBins];
	fill_n(littleBValues, numBins, 0);
	double chiSquaredTestResults[numBins];
	fill_n(chiSquaredTestResults, numBins, 0);

	double incrementLilB = (startLilB-endLilB)/numLilB;
	TString mcHistFilePath;
	for (int i = 0; i < numLilB; i++)
	{
		littleBValues[i] = startLilB + i * incrementLilB;

		mcHists[i]->SetTitle("Hmg MC (scaled);Electron Pulse Height (~keV);Experimental counts");

		scaleHistogramsTogether(expHist, mcHists[i], beginBin, endBin);
		mcHists[i]->GetXaxis()->SetRangeUser(beginEnergy, endEnergy);


		//Plot histograms together
		if(saveHistsAndPlots)
		{
			TH1* expHists[1] = { expHist };
			TH1* mcHistsToPlot[1] = { mcHists[i] };
			expHist->SetTitle("Exp.;Electron Pulse Height (~keV);Experimental counts");
			TString imagePath = TString(GRAPHS_DIR) + "Graph_Exp_MC" + nameString + "_dwcut" + int2str(i) + "_EP_EDepE.png";
			plotExpVersusMCToImage(1, expHists, 1, mcHistsToPlot, "Comparison with little b = " + d2str(littleBValues[i]), imagePath);
		}

		int ndf, igood;
		double rootChiSquaredResult;
		expHist->Chi2TestX(mcHists[i], rootChiSquaredResult, ndf, igood, "WW", res);
		chiSquaredTestResults[i] = rootChiSquaredResult;
		cout << "ChiSquared Result: " << chiSquaredTestResults[i] << endl;

	}


	TGraph chiSquaredGraph(numLilB, littleBValues, chiSquaredTestResults);
	chiSquaredGraph.SetName("chiSquaredGraph");
	TF1 quadratic("quadratic", "[2] * x * x + [1] * x + [0]");

	if(saveHistsAndPlots){
		TCanvas theCanvas;
		gPad->SetLogy(1);
		chiSquaredGraph.SetTitle(nameString+TString(";b;#tilde{#chi}^{2}"));
		chiSquaredGraph.Draw("ALP");
		chiSquaredGraph.Fit(&quadratic);
		TString imagePath = TString(GRAPHS_DIR) + "Graph_ChiSquaredTest_" + nameString + ".png";
		theCanvas.SaveAs(imagePath);
	}
	else{
		chiSquaredGraph.Fit(&quadratic,"QN");
	}

	double bestLittleb=-0.5*quadratic.GetParameter(2)/quadratic.GetParameter(1);

	return bestLittleb;
}

void makeBootStrapLittleB(TString expHistFilePath, HistDim histDim, RDK2MCAnalyzer& theAnalyzer,  int numLilB, double startLilB, double endLilB, double beginEnergy, double endEnergy)
{
	RDK2Events theEvents;

	int numRebin = 100;
	if(!FileExists(expHistFilePath))
	{
		cout << "Could not find: " << expHistFilePath << endl;
		return;
	}

	TH1* expHist = getTabSeperatedHist(expHistFilePath, histDim);
	expHist->SetName("expHist");
	expHist->Rebin(numRebin);
	expHist->GetXaxis()->SetRangeUser(beginEnergy, endEnergy);

	TH1* mcHists[numLilB];
	for(int i=0;i<numLilB;++i)
	{

		double littleB=startLilB+((endLilB-startLilB)/(double) (numLilB-1))*i;
//		double normConst=theEvents.calcNormConstant(THREEBODY,littleB,0);
		double normConst=5.94e-12; //Should be small enough
		TString mcHistPath=TString(HISTS_DIR)+"Hist_MC_"+theAnalyzer.GetName()+"_dwcut"+TString::Itoa(i,10)+"_EP_EDepEProb.txt";
		if(FileExists(mcHistPath))
		{
			cout << "Loading prob hist for little b: "<< littleB << endl;
			mcHists[i]=getTabSeperatedHist(mcHistPath,histDim);
		}
		else
		{
			cout << "Making prob hist for little b: "<< littleB  << endl;
			mcHists[i]=theAnalyzer.makeEEHistDecayWidthProb(theAnalyzer.GetName(),littleB,normConst,histDim);
			convertTH1ToTXT(mcHists[i],mcHistPath);
		}
		mcHists[i]->Rebin(numRebin);
		mcHists[i]->SetName("mcHist"+TString::Itoa(i,10));
		mcHists[i]->GetXaxis()->SetRangeUser(beginEnergy, endEnergy);
	}

	double bestLittleB0=fitLittleB(theAnalyzer.GetName(), expHist, mcHists,  numLilB, startLilB, endLilB, beginEnergy, endEnergy,true);
	cout << "Best b: "<< bestLittleB0 << endl;

//	int numberOfBootstraps=1000;
//	int numLittleBHistBins=100;
//	TH1* littleBHist=new TH1D("littleBHist","Bootstrap best fit little b on experimental data;b;Realizations",numLittleBHistBins,startLilB,endLilB);

//	for(int i=0;i<numberOfBootstraps;++i)
//	{
//		cout << "Beginning bootstrap " << i+1 << endl;
//		TH1* bootStrapExpHist=makeTH1DFromDim("bootStrapHist","bootStrapHist",histDim);
//		bootStrapExpHist->Rebin(numRebin);
//		bootStrapExpHist->FillRandom(expHist,round(expHist->Integral()));
//		double bestLittleB=fitLittleB(theAnalyzer.GetName(), expHist, mcHists,  numLilB, startLilB, endLilB, beginEnergy, endEnergy,false);
//		littleBHist->Fill(bestLittleB);
//		delete bootStrapExpHist;
//	}
//	TString imagePath = TString(GRAPHS_DIR) + "Graph_Bootstrap_BestLittleB_" + theAnalyzer.GetName() + ".png";
//	TH1* littleBHists[1]={littleBHist};
//	plotExpVersusMCToImage(0, nullptr, 1, littleBHists, "Bootstrap on experimental data, best fit little b's", imagePath);
}
