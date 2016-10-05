#include "RDK2Analysis.h"

#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphQQ.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TString.h"
#include <TLatex.h>
#include <TPaveText.h>
#include <TArrow.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#include "MRText.h"
#include "MRIO.h"
#include "MRMisc.h"
#include "MRGraphics.h"
#include "MRPhys.h"

#include "RDK2Set.h"
#include "RDK2MCAnalyzer.h"
#include "RDK2Constants.h"
#include "RDK2Vars.h"
#include "RDK2IO.h"
#include "cVector.h"

using namespace std;

//user must ensure that begin and end are reasonable
double getMLS(TH1I* unscaledHist, TH1I* toBeScaledHist, int begin, int end, int loopNum, double& bestScale)
{
	gROOT->cd();
	double middleScale = unscaledHist->GetBinContent(begin) / toBeScaledHist->GetBinContent(begin);
	double bestChiSquared, currentScale, currentChiSquared;
	int ndf = end - begin;
	bestChiSquared = REALLY_BIG_DBL;

	for (int j = 0; j < loopNum; j++)
	{
		currentScale = middleScale * (.5 + .5 * (((double) j) / (((double) loopNum) * .5)));
		currentChiSquared = getChiSquared(unscaledHist, toBeScaledHist, begin, end, currentScale);
		currentChiSquared = currentChiSquared / ((double) ndf);
		if(currentChiSquared < bestChiSquared || j == 0)
		{
			bestChiSquared = currentChiSquared;
			bestScale = currentScale;
		}

	}
	return bestChiSquared;
}

double getChiSquared(TH1I* unscaledHist, TH1I* scaledHist, int begin, int end, double scale)
{
	double chiSquaredValue = 0;
	for (int i = 0; i < end - begin + 1; i++)
		chiSquaredValue += pow(unscaledHist->GetBinContent(i + begin) - scale * scaledHist->GetBinContent(i + begin), 2) / (pow(unscaledHist->GetBinError(i + begin), 2) + pow(scale * scaledHist->GetBinError(i + begin), 2));
	return chiSquaredValue;
}

void getResiduals(TH1I* dataHist, TH1I* scaledMCHist, double* res, int numBins)
{
	for (int i = 1; i <= numBins; i++)
		res[i - 1] = (scaledMCHist->GetBinContent(i) - dataHist->GetBinContent(i)) / scaledMCHist->GetMeanError();
}

void getStandardizedResiduals(TH1I* dataHist, TH1I* scaledMCHist, double* res, int numBins)
{
	double mean = 0, sd = 0;
	getResiduals(dataHist, scaledMCHist, res, numBins);
	for (int i = 0; i < numBins; i++)
		mean += res[i];
	mean = mean / numBins;
	for (int i = 0; i < numBins; i++)
		sd += pow(res[i] - mean, 2);
	sd = sqrt(sd / numBins);
	for (int i = 0; i < numBins; i++)
		res[i] = res[i] / sd;
}

void calcResultsStats(TString firstFileName, int numFiles, bool homogeneous, bool fourBody)
{

	int e, ep, epb, epg, b, g, eb, eg;
	int numEvents;
	int hitcodepCount[MAXHITCODE + 1], hitcodeeCount[MAXHITCODE + 1], hitcodegCount[MAXHITCODE + 1];

	for (int j = 0; j <= MAXHITCODE; j++)
	{
		hitcodeeCount[j] = hitcodegCount[j] = hitcodepCount[j] = 0;
	}
	e = ep = epb = epg = b = g = eb = eg = numEvents = 0;

	TChain* resultChain = makeResultChain(firstFileName, "r", numFiles);
	Double_t tofp;
	Int_t hitcodep, hitcodee, hitcodeg;

	resultChain->SetBranchAddress("tofp", &tofp);
	resultChain->SetBranchAddress("hitcodep", &hitcodep);
	resultChain->SetBranchAddress("hitcodee", &hitcodee);
	if(fourBody) resultChain->SetBranchAddress("hitcodeg", &hitcodeg);

	TChain *dervChain;
	Int_t dwcut;
	Double_t eeBS;

	TString firstDervFileName = addBeforeExtension(firstFileName, "_Derv");
	dervChain = makeResultChain(firstDervFileName, "rd", numFiles);
	dervChain->SetBranchAddress("eeBS", &eeBS);
	if(homogeneous)
	{
		dervChain->SetBranchAddress("dwcutb10", &dwcut);
	}

	int numEntries = resultChain->GetEntries();
	for (int i = 0; i < numEntries; i++)
	{

		resultChain->GetEntry(i);
		dervChain->GetEntry(i);

		if(!homogeneous || dwcut == 0)
		{
			if(eeBS > 25) hitcodeeCount[hitcodee] += 1;
			if(fourBody) hitcodegCount[hitcodeg] += 1;
			numEvents += 1;

			if(hitcodep == 1)
			{
				if(tofp < 1.4e-6)
					hitcodepCount[8] += 1;
				else if(tofp > 25e-6)
					hitcodepCount[7] += 1;
				else
					hitcodepCount[1] += 1;
			}
			else
				hitcodepCount[hitcodep] += 1;

			//Bare
			if(fourBody && (hitcodeg >= 15 && hitcodeg <= 17))
			{
				b += 1;
				if(hitcodee == 1 && eeBS > 25)
				{
					eb += 1;
					if(hitcodep == 1 && tofp > 1.4e-6 && tofp < 25e-6)
					{
						epb += 1;
					}
				}
			}
			//BGO
			if(fourBody && (hitcodeg >= 18 && hitcodeg <= 25))
			{
				g += 1;
				if(hitcodee == 1 && eeBS > 25)
				{
					eg += 1;
					if(hitcodep == 1 && tofp > 1.4e-6 && tofp < 25e-6)
					{
						epg += 1;
					}
				}
			}
			//Normal
			if(!fourBody)
			{
				if(hitcodee == 1 && eeBS > 25)
				{
					e += 1;
					if(hitcodep == 1 && tofp > 1.4e-6 && tofp < 25e-6)
					{
						ep += 1;
					}
				}
			}

		}

	}

	double ep2e, epCoHitRatio, epgCoHitRatio, epbCoHitRatio;
	ep2e = epCoHitRatio = epgCoHitRatio = epbCoHitRatio = 0.0;
	if(!fourBody)
	{
		ep2e = ((double) ep) / ((double) e);
		epCoHitRatio = ((double) ep) / ((double) numEvents);
	}
	else
	{
		epgCoHitRatio = ((double) epg) / ((double) numEvents);
		epbCoHitRatio = ((double) epb) / ((double) numEvents);
	}

	TString statFileName = firstFileName + "_Stat.txt";
	ofstream statFile;
	statFile.open(statFileName);

	statFile << "Statitistics file for results beginning with file: " << statFileName << endl;
	statFile << "Num files: " << numFiles << endl;
	statFile << "Four Body?: " << fourBody << endl;
	statFile << "Homogeneous Events?: " << homogeneous << endl;
	statFile << "Number of events (after dwcut): " << numEvents << endl;
	if(fourBody)
	{
		statFile << "b: " << b << endl;
		statFile << "eb: " << eb << endl;
		statFile << "epb: " << epb << endl;
		statFile << "epb CoHitRatio: " << epbCoHitRatio << endl;
		statFile << "g: " << g << endl;
		statFile << "eg: " << eg << endl;
		statFile << "epg: " << epg << endl;
		statFile << "epg CoHitRatio: " << epgCoHitRatio << endl;
	}
	else
	{
		statFile << "e: " << e << endl;
		statFile << "ep: " << ep << endl;
		statFile << "ep/e: " << ep2e << endl;
		statFile << "ep CoHitRatio: " << epCoHitRatio << endl;
	}
	statFile << endl << "Hitcode Ratio:" << endl;
	statFile << "Hitcode\tProton\tElectron";
	if(fourBody) statFile << "\tGamma";
	statFile << endl;

	for (int k = 0; k <= MAXHITCODE; k++)
	{
		statFile << k << "\t" << ((double) hitcodepCount[k]) / ((double) numEvents) << "\t" << ((double) hitcodeeCount[k]) / ((double) numEvents);
		if(fourBody) statFile << "\t" << ((double) hitcodegCount[k]) / ((double) numEvents);
		statFile << endl;
	}
	statFile << endl;

	statFile.close();

	delete resultChain;
	delete dervChain;

}

void combineGammaEnergyHists()
{
	gROOT->cd();

	TFile* mrk_lowFile = new TFile("mrk_ege_hists0.root", "read");
	TFile* mrk_highFile = new TFile("mrk_ege_hists1.root", "read");
	TH1D* mrk_bareLow_hist = (TH1D*) mrk_lowFile->Get("mrk_bareLow_hist");
	TH1D* mrk_bareHigh_hist = (TH1D*) mrk_highFile->Get("mrk_bareHigh_hist");

	TFile mrk_combinedBAPDHistFile("mrk_combinedBAPDHistFile.root", "recreate");
	TH1D* mrkHist = new TH1D("mrkHist", "BGamma Energies by Bare APD", 100, 0, 50);
	mrkHist->Sumw2();

	cout << "Low/high ratio: " << mrk_bareLow_hist->Integral() / mrk_bareHigh_hist->Integral() << endl;

	*mrkHist = *mrk_bareLow_hist + *mrk_bareHigh_hist;
	mrkHist->SetName("mrkHist");
	mrkHist->SetTitle("BGamma Energies by Bare APD");

	mrkHist->GetXaxis()->SetRangeUser(0, 20);
	mrkHist->Write();
	mrkHist->SetLineColor(kOrange);

	TFile* geant_lowFile = new TFile("geant_ege_hists0.root", "read");
	TFile* geant_highFile = new TFile("geant_ege_hists1.root", "read");
	TH1D* geant_bareLow_hist = (TH1D*) geant_lowFile->Get("geant_bareLow_hist");
	TH1D* geant_bareHigh_hist = (TH1D*) geant_highFile->Get("geant_bareHigh_hist");

	TFile geant_combinedBAPDHistFile("geant_combinedBAPDHistFile.root", "recreate");
	TH1D* geantHist = new TH1D("geantHist", "BGamma Energies by Bare APD", 100, 0, 50);
	geantHist->Sumw2();

	cout << "Low/high ratio: " << geant_bareLow_hist->Integral() / geant_bareHigh_hist->Integral() << endl;

	*geantHist = *geant_bareLow_hist + *geant_bareHigh_hist;
	geantHist->SetName("geantHist");
	geantHist->SetTitle("BGamma Energies by Bare APD");

	geantHist->GetXaxis()->SetRangeUser(0, 20);
	geantHist->Write();

	convertTH1ToTXT(mrkHist, "BAPD_MRK_Hist.txt");
	convertTH1ToTXT(geantHist, "BAPD_GEANT_Hist.txt");

	TCanvas* c1 = new TCanvas("Gamma Plot", "Gamma Plot", 10, 10, 800, 600);
	c1->SetTickx(1);
	c1->SetTicky(1);
	c1->SetGridx(1);
	c1->SetGridy(1);
	mrkHist->Draw();
	geantHist->Draw("same");

	gROOT->cd();

	TH1D* jeffHist = getTabSeperatedHist("BAPD_EXP_DATA_JEFF1.txt", 60, 0., 30.);
	TH1D* benHist = getTabSeperatedHist("BAPD_EXP_DATA_BEN1.txt", 94, 0., 47.);
	TH1D* tomHist = getTabSeperatedHist("BAPD_CALC_TOM1.txt", 9, 0.5, 5.);

	jeffHist->SetLineColor(kBlue);
	benHist->SetLineColor(kRed);
	tomHist->SetLineColor(kGreen);

	jeffHist->SetMarkerColor(kBlue);
	benHist->SetMarkerColor(kRed);
	tomHist->SetMarkerColor(kGreen);
	jeffHist->Draw("same");
	benHist->Draw("same");
	tomHist->Draw("same");
	TLegend *l = new TLegend(0.5525025, 0.672956, 0.8421001, 0.8449686, nullptr, "brNDC");
	l->AddEntry(geantHist->GetName(), "Geant MC Hist", "L");
	l->AddEntry(mrkHist->GetName(), "MRK MC Hist", "L");
	l->AddEntry(tomHist->GetName(), "Tom's Calc", "L");
	l->AddEntry(benHist->GetName(), "Ben's Exp Data", "L");
	l->AddEntry(jeffHist->GetName(), "Jeff's Exp Data", "L");
	l->SetFillColor(kWhite);
	l->Draw();
	c1->Update();
	gStyle->SetOptStat(0);
	c1->SaveAs("BareAPD_Combined_Graph2.png");
	c1->Close();
	mrk_combinedBAPDHistFile.Close();
	mrk_lowFile->Close();
	mrk_highFile->Close();
	delete mrk_lowFile;
	delete mrk_highFile;

	geant_combinedBAPDHistFile.Close();
	geant_lowFile->Close();
	geant_highFile->Close();
	delete geant_lowFile;
	delete geant_highFile;
}

TCanvas* makeRadRangeComparison()
{

	TFile f0("MJB-Evn-TEST-n1.B3-00-20100513.root", "READ");
	TFile f1("MJB-Evn-TEST-n1.B3-01-20100513.root", "READ");
	TFile f2("MJB-Evn-TEST-n1.B3-02-20100513.root", "READ");
	TFile f3("MJB-Evn-TEST-n1.B3-03-20100513.root", "READ");
	TFile f4("MJB-Evn-TEST-n1.B3-04-20100513.root", "READ");
	TFile f5("MJB-Evn-TEST-n1.B3-05-20100513.root", "READ");

	TTree* t0 = (TTree*) f0.Get("t");
	TTree* t1 = (TTree*) f1.Get("t");
	TTree* t2 = (TTree*) f2.Get("t");
	TTree* t3 = (TTree*) f3.Get("t");
	TTree* t4 = (TTree*) f4.Get("t");
	TTree* t5 = (TTree*) f5.Get("t");

	gROOT->cd();
	double bottomRange = 0;
	double topRange = 780;
	TH1I* h0 = new TH1I("h0", "3 Body Decay", 100, bottomRange, topRange);
	TH1I* h1 = new TH1I("h1", "4 Body Decay - 0.1 to 2 keV", 100, bottomRange, topRange);
	TH1I* h2 = new TH1I("h2", "4 Body Decay - 15 to 340 keV", 100, bottomRange, topRange);
	TH1I* h3 = new TH1I("h3", "4 Body Decay - 300 to 340 keV", 100, bottomRange, topRange);
	TH1I* h4 = new TH1I("h4", "4 Body Decay - 2 to 780 keV", 100, bottomRange, topRange);
	TH1I* h5 = new TH1I("h5", "4 Body Decay - 0.01 to 2 keV", 100, bottomRange, topRange);

	h1->SetLineColor(kRed);
	h2->SetLineColor(kGreen);
	h3->SetLineColor(kBlue);
	h4->SetLineColor(kOrange);
	h5->SetLineColor(kViolet);
	h3->GetXaxis()->SetTitle("Electron Kinetic Energy (keV)");
	h3->GetYaxis()->SetTitle("Events per 1 million");
	h3->GetYaxis()->SetTitleOffset(1.5);
	h2->GetXaxis()->SetTitle("Electron Kinetic Energy (keV)");
	h2->GetYaxis()->SetTitle("Events per 1 million");
	h2->GetYaxis()->SetTitleOffset(1.5);

	t0->Draw("ee0 >> h0");
	t1->Draw("ee0 >> h1");
	t2->Draw("ee0 >> h2");
	t3->Draw("ee0 >> h3");
	t4->Draw("ee0 >> h4");
	t5->Draw("ee0 >> h5");

	f0.Close();
	f1.Close();
	f2.Close();
	f3.Close();
	f4.Close();
	f5.Close();

	h0->SetStats(0);
	h1->SetStats(0);
	h2->SetStats(0);
	h3->SetStats(0);
	h4->SetStats(0);
	h5->SetStats(0);

	TCanvas *c1_n2 = new TCanvas("c1_n2", "Electron KE Comparison", 14, 49, 800, 600);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	c1_n2->Range(-122.4971, -3631.076, 837.6914, 31897.61);
	c1_n2->SetBorderSize(2);
	c1_n2->SetLeftMargin(0.1275761);
	c1_n2->SetRightMargin(0.04);
	c1_n2->SetTopMargin(0.04);
	c1_n2->SetFrameFillColor(0);

	h3->Draw();
	h1->Draw("same");
	h2->Draw("same");
	h0->Draw("same");
	h4->Draw("same");
	h5->Draw("same");

	TLegend *l = new TLegend(0.5525025, 0.672956, 0.9421001, 0.9449686, nullptr, "brNDC");
	l->AddEntry(h0->GetName(), h0->GetTitle(), "L");
	l->AddEntry(h1->GetName(), h1->GetTitle(), "L");
	l->AddEntry(h2->GetName(), h2->GetTitle(), "L");
	l->AddEntry(h3->GetName(), h3->GetTitle(), "L");
	l->AddEntry(h4->GetName(), h4->GetTitle(), "L");
	l->AddEntry(h5->GetName(), h5->GetTitle(), "L");

	l->Draw();
	c1_n2->Update();

	return c1_n2;

}

void makeEG0Combined()
{
	gROOT->cd();
	TString fileName;
	fileName = EVENTS_DIR;
	fileName += "MJB-Evn-radmid-n1.B4-0-20120208.root";
	TChain* eventMidChain = makeEventChain(fileName, DEFAULT_EVENTS_TREENAME, 20);
	fileName = EVENTS_DIR;
	fileName += "MJB-Evn-radlow-n1.B4-0-20120208.root";
	TChain* eventLowChain = makeEventChain(fileName, DEFAULT_EVENTS_TREENAME, 20);
	fileName = EVENTS_DIR;
	fileName += "MJB-Evn-radhigh-n1.B4-0-20120208.root";
	TChain* eventHighChain = makeEventChain(fileName, DEFAULT_EVENTS_TREENAME, 20);
	TH1D* lowHist = makeTH1DFromDim("lowHist", "Gamma Energies (Lower Energy Events)", BGO_DIM);
	TH1D* midHist = makeTH1DFromDim("midHist", "Gamma Energies (Middle Energy Events)", BGO_DIM);
	TH1D* highHist = makeTH1DFromDim("highHist", "Gamma Energies (Middle Energy Events)", BGO_DIM);

	cout << "Creating mid histogram" << endl;
	eventMidChain->Draw("eg0 >> midHist");
	cout << "Creating low histogram" << endl;
	eventLowChain->Draw("eg0 >> lowHist");

	cout << "Creating high histogram" << endl;
	eventHighChain->Draw("eg0 >> highHist");

	double lowBranchingRatio = 0.00429224;
	cout << "Low Branching Ratio: " << lowBranchingRatio << endl;

	double midBranchingRatio = 0.00318548;
	cout << "Mid Branching Ratio: " << midBranchingRatio << endl << endl;

	double highBranchingRatio = 0.00252823;
	cout << "High Branching Ratio: " << highBranchingRatio << endl << endl;

	midHist->Sumw2();
	lowHist->Sumw2();
	highHist->Sumw2();

	TH1D* eg0Hist = makeTH1DFromDim("eg0Hist", "Gamma Energies", BGO_DIM);

	midHist->Scale(midBranchingRatio / eventMidChain->GetEntries());
	lowHist->Scale(lowBranchingRatio / eventLowChain->GetEntries());
	highHist->Scale(highBranchingRatio / eventHighChain->GetEntries());

	delete eventHighChain;
	delete eventLowChain;
	delete eventMidChain;

	*eg0Hist = *midHist + *lowHist + *highHist;

	fileName = HISTS_DIR;
	fileName += "eg0Hist20120208.txt";
	convertTH1ToTXT(eg0Hist, fileName);

	TCanvas* theCanvas = new TCanvas();
	lowHist->Draw();
	midHist->Draw("same");
	highHist->Draw("same");
	theCanvas->Update();
	gROOT->cd();
	TString filePath = GRAPHS_DIR;
	filePath += "ee0_radlowToHighNew.png";
	theCanvas->SaveAs(filePath);

	theCanvas->Close();
	delete lowHist;
	delete midHist;
	delete highHist;
	delete eg0Hist;
}

void createEGOpeningAngleGraph(int type)
{
	int numFiles = 4;
	int numEvents = numFiles * 1000000;
	TString typeName;
	TChain* eventChain;
	if(type == 0)
	{
		typeName = "RadLow";
		eventChain = makeEventChain("MJB-Evn-radlow-n1.B4-00-20100517.root", DEFAULT_EVENTS_TREENAME, numFiles);
	}
	else
	{
		typeName = "RadMid";
		eventChain = makeEventChain("MJB-Evn-radmid-n1.B4-00-20100924.root", DEFAULT_EVENTS_TREENAME, numFiles);
	}

	double mxg0, myg0, mzg0, mxe0, mye0, mze0;
	eventChain->SetBranchAddress("mxg0", &mxg0);
	eventChain->SetBranchAddress("myg0", &myg0);
	eventChain->SetBranchAddress("mzg0", &mzg0);
	eventChain->SetBranchAddress("mxe0", &mxe0);
	eventChain->SetBranchAddress("mye0", &mye0);
	eventChain->SetBranchAddress("mze0", &mze0);

	TH1D angleHist("angleHist", "Opening angle between gamma and electron", 100, 0, PI);

	cVector3D mg, me;
	double sinAngle, angle;
	for (int i = 0; i < numEvents; i++)
	{
		eventChain->GetEntry(i);
		mg.setVal(mxg0, myg0, mzg0);
		me.setVal(mxe0, mye0, mze0);

		sinAngle = mg.dotProd(me);
		angle = acos(sinAngle);
		angleHist.Fill(angle);
	}

	convertTH1ToTXT(&angleHist, "egOpeningAngleHist.txt");
	TCanvas* c1 = new TCanvas("Gamma Plot", "Gamma Plot", 10, 10, 800, 600);
	c1->SetTickx(1);
	c1->SetTicky(1);
	c1->SetGridx(1);
	c1->SetGridy(1);
	angleHist.Draw();

	c1->SaveAs("egAngleHist.png");

	c1->Close();

	delete c1;
	delete eventChain;
}

TCanvas* plotManyHistogramsTogether(int numHists, TH1** theHists, TLegend* theLegend, TString inpTitle, TString canvasName)
{
	gROOT->cd();
	int colorList[12] = { 1, 920, 632, 416, 600, 616, 432, 800, 820, 840, 860, 880 };
	gStyle->SetOptStat(0);
	TCanvas* c2 = new TCanvas(canvasName, "plotManyHistCanvas", 10, 10, 1200, 1024);
	c2->SetFillColor(kGray);
	c2->cd();
	c2->SetLeftMargin(.1);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);

	TString plotTitle;
	if(inpTitle == "")
		plotTitle = theHists[0]->GetTitle();
	else
		plotTitle = inpTitle;

	if(numHists < 1) return nullptr;

	theHists[0]->SetLineColor(colorList[0]);

//    double maximum=theHists[0]->GetMaximum();
	//double maximum=10000;
	theHists[0]->Draw("e1");
	theHists[0]->GetYaxis()->SetTitleOffset(1.8);
	theHists[0]->SetLineColor(kRed);

	if(theLegend != nullptr)
	{
		theLegend->AddEntry(theHists[0], theHists[0]->GetName(), "PL");
		theLegend->SetBorderSize(0);
	}

	for (int i = 1; i < numHists; i++)
	{

		//theScale=theHists[0]->Integral(450,800)/theHists[i]->Integral(450,800);
		//theScale=theHists[0]->Integral(21,301)/theHists[i]->Integral(21,301);
		//theHists[i]->Scale(theScale);
		theHists[i]->SetLineColor(colorList[(i + 2) % 12]);
		if(theLegend != nullptr) theLegend->AddEntry(theHists[i], theHists[i]->GetName(), "PL");
		theHists[i]->Draw("e1same");
		//if(theHists[i]->GetMaximum() > maximum)
		//    maximum=theHists[i]->GetMaximum();
		//theHists[i]->Rebin(16);

	}
	// theHists[0]->Rebin(16);
	// theHists[0]->GetYaxis()->SetRangeUser(0.,maximum*1.2);
	if(theLegend != nullptr)
	{
		theLegend->SetFillColor(kWhite);
		theLegend->Draw();
	}
	//TLatex* latexTitle = new TLatex();
	//latexTitle->SetTextSize(.06);
	//double yCoord=maximum*1.05;
	//double xCoord=theHists[0]->GetBinLowEdge(theHists[0]->GetXaxis()->GetFirst());
	// latexTitle->DrawLatex(.1,yCoord,plotTitle); // coordinates are absolute here

	c2->Update();

	return c2;
}

void makeHistWithGeantandMRK(TString mrkFileName, TString geantProtonFileName, TString geantElectronFileName, TString littlebFileName, TString eventFileName, unsigned int numFiles)
{
	TString resultsDir = RESULTS_DIR;
	TString eventsDir = EVENTS_DIR;
	TString littlebDir = LITTLEB_DIR;

	TChain* mrkChain = makeResultChain(resultsDir + mrkFileName, "r", numFiles);
	if(mrkChain->IsZombie())
	{
		cout << "Error in opening " << mrkFileName << endl;
		return;
	}

	TChain* geantProtonChain = makeResultChain(resultsDir + geantProtonFileName, "SBDFluxTree", numFiles);
	if(geantProtonChain->IsZombie())
	{
		cout << "Error in opening " << geantProtonFileName << endl;
		return;
	}

	TChain* geantElectronChain = makeResultChain(resultsDir + geantElectronFileName, "SBDFluxTree", numFiles);
	if(geantElectronChain->IsZombie())
	{
		cout << "Error in opening " << geantElectronFileName << endl;
		return;
	}

	TChain* littlebChain = makeEventChain(littlebDir + littlebFileName, "lilbTree", numFiles);
	if(littlebChain->IsZombie())
	{
		cout << "Error in opening " << littlebFileName << endl;
		return;
	}

	TChain* eventChain = makeEventChain(eventsDir + eventFileName, "t", numFiles);
	if(eventChain->IsZombie())
	{
		cout << "Error in opening " << eventFileName << endl;
		return;
	}

	gROOT->cd();
	double protonToFBinMin = 2.02e-6, electronToFBinMin = 0, protonToFBinMax = 24.98e-6, electronToFBinMax = 12e-9;
	double protonKEBinMin = 7, electronKEBinMin = 50, protonKEBinMax = 31, electronKEBinMax = 800;
	double protonKE0BinMin = 0, electronKE0BinMin = 0, protonKE0BinMax = 0.8, electronKE0BinMax = 800;
	int numBinsToF = 41, numBinsKE = 100, numBinsKE0 = 100;

	TString saveFileName = HISTS_DIR;
	saveFileName += "mrkAndGeant3DHists120904.root";
	TFile saveFile(saveFileName, "recreate");
	saveFile.cd();
	TH3I* mrkP3DHist = new TH3I("mrkP3DHist", "MRK Monte Carlo for Proton", numBinsToF, protonToFBinMin, protonToFBinMax, numBinsKE, protonKEBinMin, protonKEBinMax, numBinsKE0, protonKE0BinMin, protonKE0BinMax);
	TH3I* mrkE3DHist = new TH3I("mrkE3DHist", "MRK Monte Carlo for Electron", numBinsToF, electronToFBinMin, electronToFBinMax, numBinsKE, electronKEBinMin, electronKEBinMax, numBinsKE0, electronKE0BinMin, electronKE0BinMax);
	TH3I* geantP3DHist = new TH3I("geantP3DHist", "Geant4 Monte Carlo for Proton", numBinsToF, protonToFBinMin, protonToFBinMax, numBinsKE, protonKEBinMin, protonKEBinMax, numBinsKE0, protonKE0BinMin, protonKE0BinMax);
	TH3I* geantE3DHist = new TH3I("geantE3DHist", "Geant4 Monte Carlo for Electron", numBinsToF, electronToFBinMin, electronToFBinMax, numBinsKE, electronKEBinMin, electronKEBinMax, numBinsKE0, electronKE0BinMin, electronKE0BinMax);

	TString xAxisString = "Time of Flight (s)";
	TString yAxisString = "Final Kinetic Energy (keV)";
	TString zAxisString = "Initial Kinetic Energy (keV)";

	mrkP3DHist->GetXaxis()->SetTitle(xAxisString);
	mrkE3DHist->GetXaxis()->SetTitle(xAxisString);
	geantP3DHist->GetXaxis()->SetTitle(xAxisString);
	geantE3DHist->GetXaxis()->SetTitle(xAxisString);

	mrkP3DHist->GetYaxis()->SetTitle(yAxisString);
	mrkE3DHist->GetYaxis()->SetTitle(yAxisString);
	geantP3DHist->GetYaxis()->SetTitle(yAxisString);
	geantE3DHist->GetYaxis()->SetTitle(yAxisString);

	mrkP3DHist->GetZaxis()->SetTitle(zAxisString);
	mrkE3DHist->GetZaxis()->SetTitle(zAxisString);
	geantP3DHist->GetZaxis()->SetTitle(zAxisString);
	geantE3DHist->GetZaxis()->SetTitle(zAxisString);

	double vxe, vye, vze;
	double vxp, vyp, vzp;
	double xpf;
	double xef;
	double eKE, pKE;
	double tofp, tofe;
	double SBDTimeElectron, SBDKEElectron;
	double SBDTimeProton, SBDKEProton;
	double SBDPosXProton, SBDPosXElectron;
	int geantProtonFluxEvent, geantElectronFluxEvent;
	int hitcodep, hitcodee;
	int dwcutb;
	double z0, ee0, ep0;
	char SBDHitCodeProton, SBDParProton;
	char SBDHitCodeElectron, SBDParElectron;

	mrkChain->SetBranchAddress("hitcodep", &hitcodep);
	mrkChain->SetBranchAddress("hitcodee", &hitcodee);
	mrkChain->SetBranchAddress("tofp", &tofp);
	mrkChain->SetBranchAddress("tofe", &tofe);
	mrkChain->SetBranchAddress("vxef", &vxe);
	mrkChain->SetBranchAddress("vyef", &vye);
	mrkChain->SetBranchAddress("vzef", &vze);
	mrkChain->SetBranchAddress("vxpf", &vxp);
	mrkChain->SetBranchAddress("vypf", &vyp);
	mrkChain->SetBranchAddress("vzpf", &vzp);
	mrkChain->SetBranchAddress("xpf", &xpf);
	mrkChain->SetBranchAddress("xef", &xef);

	geantProtonChain->SetBranchAddress("Time", &SBDTimeProton);
	geantProtonChain->SetBranchAddress("KE", &SBDKEProton);
	geantProtonChain->SetBranchAddress("Event", &geantProtonFluxEvent);
	geantProtonChain->SetBranchAddress("PosX", &SBDPosXProton);
	geantProtonChain->SetBranchAddress("HitCode", &SBDHitCodeProton);
	geantProtonChain->SetBranchAddress("Particle", &SBDParProton);

	geantElectronChain->SetBranchAddress("Time", &SBDTimeElectron);
	geantElectronChain->SetBranchAddress("KE", &SBDKEElectron);
	geantElectronChain->SetBranchAddress("Event", &geantElectronFluxEvent);
	geantElectronChain->SetBranchAddress("PosX", &SBDPosXElectron);
	geantElectronChain->SetBranchAddress("HitCode", &SBDHitCodeElectron);
	geantElectronChain->SetBranchAddress("Particle", &SBDParElectron);

	littlebChain->SetBranchAddress("dwcutb20", &dwcutb);

	eventChain->SetBranchAddress("z0", &z0);
	eventChain->SetBranchAddress("ee0", &ee0);
	eventChain->SetBranchAddress("ep0", &ep0);

	int lastFluxCheckedProton = 0;
	int lastFluxCheckedElectron = 0;
	int currentFile = -1;
	int targetEvent;
	bool firstProtonFluxEventFound, firstElectronFluxEventFound;
	int previousEvent;

	for (int i = 0; i < mrkChain->GetEntries(); i++)
	{
		//Ensure that the FluxTree (since it is not for every event, carries over correctly when file changes)
		if(i % 1000000 == 0)
		{
			currentFile++;

			previousEvent = 0;
			for (int j = lastFluxCheckedProton; j < geantProtonChain->GetEntries() && currentFile != 0; j++)
			{
				geantProtonChain->GetEntry(j);
				if(previousEvent >= geantProtonFluxEvent)
				{
					lastFluxCheckedProton = j;
					break;
				}
				previousEvent = geantProtonFluxEvent;

			}

			previousEvent = 0;
			for (int j = lastFluxCheckedElectron; j < geantElectronChain->GetEntries() && currentFile != 0; j++)
			{
				geantElectronChain->GetEntry(j);
				if(previousEvent >= geantElectronFluxEvent)
				{
					lastFluxCheckedElectron = j;
					break;
				}
				previousEvent = geantElectronFluxEvent;
			}
		}

		//First ensure that the event is not a cut event (we are using little b =0 for now
		littlebChain->GetEntry(i);
		if(dwcutb == 0)
		{
			mrkChain->GetEntry(i);
			eventChain->GetEntry(i);

			//If both the proton and electron hit, fill out the MRK Histograms
			if(hitcodee == 1 && hitcodep == 1)
			{

				pKE = convertVelocityToKE(vxp, vyp, vzp, PMASSE);
				eKE = convertVelocityToKE(vxe, vye, vze, EMASSE);

				if(pKE > 7 && pKE < 31 && eKE > 50 && eKE < 800)
				{
					mrkP3DHist->Fill(tofp, pKE, ep0);
					mrkE3DHist->Fill(tofe, eKE, ee0);
				}
			}

			//I want to look at what Geant does in isolation so use hit information from MRK for now
			if(hitcodee == 1)
			{
				firstProtonFluxEventFound = false;
				//proton
				for (int j = lastFluxCheckedProton; j < geantProtonChain->GetEntries(); j++)
				{

					geantProtonChain->GetEntry(j);
					targetEvent = (i - 1000000 * currentFile);
					if(geantProtonFluxEvent == targetEvent && !firstProtonFluxEventFound && SBDHitCodeProton == 1 && SBDParProton == 1)
					{
						geantP3DHist->Fill(SBDTimeProton, SBDKEProton, ep0);
						firstProtonFluxEventFound = true;
					}
					else if(geantProtonFluxEvent > i - 1000000 * currentFile)
					{
						lastFluxCheckedProton = j;
						break;
					}

				}
			}
			if(hitcodep == 1)
			{
				firstElectronFluxEventFound = false;
				//electron
				for (int j = lastFluxCheckedElectron; j < geantElectronChain->GetEntries(); j++)
				{
					geantElectronChain->GetEntry(j);
					targetEvent = (i - 1000000 * currentFile);
					if(geantElectronFluxEvent == targetEvent && !firstElectronFluxEventFound && SBDHitCodeElectron == 1 && SBDParElectron == 2)
					{
						geantE3DHist->Fill(SBDTimeElectron, SBDKEElectron, ee0);
						firstElectronFluxEventFound = true;
					}
					else if(geantElectronFluxEvent > i - 1000000 * currentFile)
					{
						lastFluxCheckedElectron = j;
						break;
					}

				}
			}

		}

	}

	saveFile.cd();
	mrkE3DHist->Write();
	mrkP3DHist->Write();
	geantE3DHist->Write();
	geantP3DHist->Write();

	saveFile.Close();

	delete mrkChain;
	delete eventChain;
	delete geantProtonChain;
	delete geantElectronChain;
	delete littlebChain;
}

void makeMRKAndGeantHistograms(TString axis, ParType parType)
{
	TString parString;
	if(parType == ELECTRON)
	{
		parString = "E";
	}
	else
	{
		parString = "P";
	}
	TString histFileName = HISTS_DIR;
	histFileName += "mrkAndGeant3DHists120904.root";
	TFile histFile(histFileName, "read");

	TString mrkHistString = "mrk" + parString + "3DHist";
	TString geantHistString = "geant" + parString + "3DHist";

	TH3I* mrk3DHist = (TH3I*) histFile.Get(mrkHistString);
	TH3I* geant3DHist = (TH3I*) histFile.Get(geantHistString);

	histFileName = HISTS_DIR;
	histFileName += "old/mrkAndGeant3DHistsNew.root";
	TFile histFileNew(histFileName, "read");
	TH3I* geant3DHistNew = (TH3I*) histFileNew.Get(geantHistString);

	TLegend* l = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	TCanvas* c2;

	TH1* theHists[3];

	theHists[2] = geant3DHistNew->Project3D(axis);
	TString newName = theHists[2]->GetName();
	newName += "_NEW";
	theHists[2]->SetName(newName);
	theHists[0] = mrk3DHist->Project3D(axis);
	TString HistDir = HISTS_DIR;

	theHists[1] = geant3DHist->Project3D(axis);

	convertTH1ToTXT(theHists[0], HistDir + "MRK_Geant4_Comparison_ToFHist_MRK.txt");
	convertTH1ToTXT(theHists[1], HistDir + "MRK_Geant4_Comparison_ToFHist_GEANT4OLD.txt");
	convertTH1ToTXT(theHists[2], HistDir + "MRK_Geant4_Comparison_ToFHist_GEANT4NEW.txt");

	c2 = plotManyHistogramsTogether(3, theHists, l);

	TString axisVar;
	if(axis == "x")
	{
		axisVar = "tof";
	}
	else if(axis == "y")
	{
		axisVar = "ke";
	}

	TString pictureName = GRAPHS_DIR;
	pictureName += axisVar + parString + "ComplexComparison1107022.png";

	c2->SaveAs(pictureName);
	c2->Close();

	delete c2;
	delete l;

	histFile.Close();

}

void checkGeantRKAccuracyToStandard(TString standardFileName, TString testFileName, int numEvents, double accuracyReq)
{
	TFile standardFile(standardFileName, "read");
	TFile testFile(testFileName, "read");

	TTree* standardFluxTree = (TTree*) standardFile.Get("FluxTree");
	TTree* testFluxTree = (TTree*) testFile.Get("SBDFluxTree");

	double standardTime, standardKE;
	int standardEvent;
	standardFluxTree->SetBranchAddress("SBDTime", &standardTime);
	standardFluxTree->SetBranchAddress("SBDKE", &standardKE);
	standardFluxTree->SetBranchAddress("SBDEvn", &standardEvent);

	double testTime, testKE;
	int testEvent;
	testFluxTree->SetBranchAddress("Time", &testTime);
	testFluxTree->SetBranchAddress("KE", &testKE);
	testFluxTree->SetBranchAddress("Event", &testEvent);

	int standardFluxEvent = 0, testFluxEvent = 0;
	int agreementCount = 0, hitButNotAccurateCount = 0;
	int standardHitOnlyCount = 0, testHitOnlyCount = 0;
	bool standardFluxTreeEnd = false, testFluxTreeEnd = false;

	for (int i = 0; i < numEvents; i++)
	{
		standardFluxTree->GetEntry(standardFluxEvent);
		for (; standardEvent < i;) //get it to or past the approriate real (not Flux) event
		{
			standardFluxEvent++;
			if(standardFluxEvent >= standardFluxTree->GetEntries())
			{
				standardFluxTreeEnd = true;
				break;
			}
			standardFluxTree->GetEntry(standardFluxEvent);
		}
		testFluxTree->GetEntry(testFluxEvent);
		for (; testEvent < i;) //get it to or past the appropriate real event
		{
			testFluxEvent++;
			if(testFluxEvent >= testFluxTree->GetEntries())
			{
				testFluxTreeEnd = true;
				break;
			}
			testFluxTree->GetEntry(testFluxEvent);
		}

		if(standardEvent == i && !standardFluxTreeEnd)
		{
			if(testEvent == i) //Event match
			{
				if(abs(standardTime / testTime - 1.) < accuracyReq)
				{
					if(abs(standardKE / testKE - 1.) < accuracyReq)
						agreementCount++;
					else
						hitButNotAccurateCount++;
				}
				else
					hitButNotAccurateCount++;
			}
			else
				standardHitOnlyCount++;

		}
		else
		{
			if(testEvent != i && !testFluxTreeEnd)
				agreementCount++;
			else
				testHitOnlyCount++;
		}
		if(standardFluxTreeEnd && testFluxTreeEnd) agreementCount++;

	}
	double agreementPercentage = agreementCount * 100. / numEvents;
	double agreementPercentageWOAcuracyReq = (agreementCount + hitButNotAccurateCount) * 100. / numEvents;
	double standardHitOnlyPercentage = standardHitOnlyCount * 100. / numEvents;
	double testHitOnlyPercentage = testHitOnlyCount * 100. / numEvents;
	cout << "Agreement @ " << accuracyReq * 100 << "% accuracy on ToF and KE is " << agreementPercentage << "%." << endl;
	cout << "Agreement with no accuracy requirement " << agreementPercentageWOAcuracyReq << "%." << endl;
	cout << "Standard Hit Only % " << standardHitOnlyPercentage << "%." << endl;
	cout << "Test Hit Only % " << testHitOnlyPercentage << "%." << endl;
}

void makeRadHighHist()
{
	int numFiles = 1;
	TString geantRunString = "geantliv110414";
	TChain* resultsPChain;
	TChain* resultsEChain;
	TChain* resultsGChain;
	TString firstPFile = "Results110414_radhigh_P_Set0_" + geantRunString + ".root";
	TString firstEFile = "Results110414_radhigh_E_Set0_" + geantRunString + ".root";
	TString firstGFile = "Results110414_radhigh_G_Set0_" + geantRunString + ".root";
	resultsPChain = makeResultChain(firstPFile, "g", numFiles);
	resultsEChain = makeResultChain(firstEFile, "g", numFiles);
	resultsGChain = makeResultChain(firstGFile, "g", numFiles);
	resultsEChain->AddFriend(resultsPChain, "gp");
	resultsEChain->AddFriend(resultsGChain, "gg");

	TString drawString = "g.SBDEDepTotal >> ";

	TString finalDrawString;
	TCut ekeCut, pCut, gCut;
	ekeCut = "g.SBDEDepTotal > 25";
	pCut = "gp.SBDEDepTotal > 9.3 && gp.SBDEDepTotal < 30.9 && gp.SBDTimeFirst > 1.4e3 && gp.SBDTimeFirst < 25e3"; //SBDTime First is in ns at the moment
	gCut = "gg.BGO1EDepTotal > 10 || gg.BGO2EDepTotal > 10 || gg.BGO3EDepTotal > 10 || gg.BGO4EDepTotal > 10 || gg.BGO5EDepTotal > 10 || gg.BGO6EDepTotal > 10 || gg.BGO7EDepTotal > 10 || gg.BGO8EDepTotal > 10 || gg.BGO9EDepTotal > 10 || gg.BGO10EDepTotal > 10 || gg.BGO11EDepTotal > 10 || gg.BGO12EDepTotal > 10";
	//pCut="gp.SBDEDepTotal > 0 && gp.SBDTimeFirst > 1.4e3 && gp.SBDTimeFirst < 25e3"; //SBDTime First is in ns at the moment

	TString mcHistString = "EKE_EPG_BGO";
	finalDrawString = drawString + mcHistString;
	TFile* saveFile = new TFile("EKE_EPG_BGO_110414.root", "recreate");
	saveFile->cd();
	TH1D* mcHist = new TH1D(mcHistString, "Electron Kinetic Energy in EPG coincidence (BGOS)", 50, 0, 1000);

	cout << "drawString: " << finalDrawString << endl;
	TCut theCut;
	theCut = ekeCut && pCut && gCut;
	cout << "cutString: " << theCut.GetTitle() << endl;
	resultsEChain->Draw(finalDrawString, theCut);

	saveFile->cd();
	mcHist->Write();
	TH1D* expHist = getTabSeperatedHist("JeffEPGElectronSpectrum.txt", 250, 0., 800.);
	expHist->Write();

	saveFile->Close();

	delete resultsEChain;
	delete resultsPChain;
	delete resultsGChain;

}

TCanvas* plotBNLSpectra500(int num)
{
	const int numHists = 18;
	TH1* theHists[numHists];

	TFile* f = new TFile("/home/mjbales/school/rdk/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root", "read");

	//theHists[0]=(TH1I*) f->Get("apr11_10_20");
	//theHists[1]=(TH1I*) f->Get("apr11_9_2");
	theHists[0] = (TH1I*) f->Get("apr11_9_3");
	theHists[1] = (TH1I*) f->Get("apr11_9_4");
	theHists[2] = (TH1I*) f->Get("apr11_9_6");
	theHists[3] = (TH1I*) f->Get("apr11_9_7");
	//theHists[6]=(TH1I*) f->Get("apr11_10_13");
	//theHists[7]=(TH1I*) f->Get("apr11_10_19");
	//theHists[8]=(TH1I*) f->Get("apr11_10_21");
	//theHists[9]=(TH1I*) f->Get("apr11_7_1");
	//theHists[10]=(TH1I*) f->Get("apr11_9_1");
	// theHists[11]=(TH1I*) f->Get("apr11_10_18");
	//theHists[12]=(TH1I*) f->Get("apr11_10_17");
	// theHists[13]=(TH1I*) f->Get("apr11_10_16");

	theHists[4] = (TH1I*) f->Get("apr11_9_42");
	theHists[5] = (TH1I*) f->Get("apr11_9_43");
	// theHists[16]=(TH1I*) f->Get("apr11_10_14");
	//theHists[17]=(TH1I*) f->Get("apr11_10_15");

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	TCanvas* c2 = plotManyHistogramsTogether(num, theHists, theLegend);

	theHists[0]->SetTitle("BNL 500 eV Runs");
	theHists[0]->GetXaxis()->SetTitle("Channel #");
	theHists[0]->GetYaxis()->SetTitle("Counts (AU)");
	theHists[0]->GetXaxis()->SetRangeUser(0., 100.);
	theHists[0]->GetYaxis()->SetRangeUser(0, 10000);

	return c2;

}

TCanvas* plotBNLSpectra1500(int num)
{
	const int numHists = 20;
	TH1* theHists[numHists];

	TFile* f = new TFile("/home/mjbales/school/rdk/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root", "read");

	/*
	 theHists[0]=(TH1I*) f->Get("apr11_7_16");
	 theHists[1]=(TH1I*) f->Get("apr11_7_17");
	 theHists[2]=(TH1I*) f->Get("apr11_7_18");
	 theHists[3]=(TH1I*) f->Get("apr11_7_20");
	 theHists[4]=(TH1I*) f->Get("apr11_7_21");
	 theHists[5]=(TH1I*) f->Get("apr11_7_22");
	 theHists[6]=(TH1I*) f->Get("apr11_7_23");
	 theHists[7]=(TH1I*) f->Get("apr11_7_24");
	 theHists[8]=(TH1I*) f->Get("apr11_7_25");
	 */

	theHists[0] = (TH1I*) f->Get("apr11_8_4");
	theHists[1] = (TH1I*) f->Get("apr11_8_5");
	theHists[2] = (TH1I*) f->Get("apr11_8_6");
	theHists[3] = (TH1I*) f->Get("apr11_8_7");
	theHists[4] = (TH1I*) f->Get("apr11_8_8");
	theHists[5] = (TH1I*) f->Get("apr11_8_9");
	theHists[6] = (TH1I*) f->Get("apr11_8_11");
	theHists[7] = (TH1I*) f->Get("apr11_8_12");
	theHists[8] = (TH1I*) f->Get("apr11_8_13");
	theHists[9] = (TH1I*) f->Get("apr11_8_14");
	theHists[10] = (TH1I*) f->Get("apr11_6_1");

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(num, theHists, theLegend);

	theHists[0]->SetTitle("BNL 1500 keV Runs");
	theHists[0]->GetXaxis()->SetTitle("Channel #");
	theHists[0]->GetYaxis()->SetTitle("Counts (AU)");
	theHists[0]->GetXaxis()->SetRangeUser(0, 200);
	theHists[0]->GetYaxis()->SetRangeUser(0, 8000);

	c2->Draw();

	return c2;

}

TCanvas* plotBNLSpectra750(int num)
{
	const int numHists = 21;
	TH1* theHists[numHists];

	TFile* f = new TFile("/home/mjbales/school/rdk/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root", "read");

	theHists[0] = (TH1I*) f->Get("apr11_7_3");
	theHists[1] = (TH1I*) f->Get("apr11_7_5");
	theHists[2] = (TH1I*) f->Get("apr11_7_7");
	theHists[3] = (TH1I*) f->Get("apr11_7_8");
	theHists[4] = (TH1I*) f->Get("apr11_7_9");
	theHists[5] = (TH1I*) f->Get("apr11_7_10");
	theHists[6] = (TH1I*) f->Get("apr11_7_11");
	theHists[7] = (TH1I*) f->Get("apr11_7_12");
	theHists[8] = (TH1I*) f->Get("apr11_7_13");
	theHists[9] = (TH1I*) f->Get("apr11_7_14");
	theHists[10] = (TH1I*) f->Get("apr11_8_21");
	theHists[11] = (TH1I*) f->Get("apr11_8_22");
	theHists[12] = (TH1I*) f->Get("apr11_9_12");
	theHists[13] = (TH1I*) f->Get("apr11_9_13");
	theHists[14] = (TH1I*) f->Get("apr11_10_36");
	theHists[15] = (TH1I*) f->Get("apr11_10_37");
	theHists[16] = (TH1I*) f->Get("apr11_10_38");
	theHists[17] = (TH1I*) f->Get("apr11_10_39");
	theHists[18] = (TH1I*) f->Get("apr11_10_40");
	theHists[19] = (TH1I*) f->Get("apr11_10_41");
	theHists[20] = (TH1I*) f->Get("apr11_7_6");

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(num, theHists, theLegend);

	theHists[0]->SetTitle("BNL 750 eV Runs");
	theHists[0]->GetXaxis()->SetTitle("Channel #");
	theHists[0]->GetYaxis()->SetTitle("Counts (AU)");
	theHists[0]->GetXaxis()->SetRangeUser(0, 100);
	theHists[0]->GetYaxis()->SetRangeUser(0, 8000);

	c2->Draw();

	return c2;

}

TCanvas* plotBNLSpectra1400(int num)
{
	const int numHists = 13;
	TH1* theHists[numHists];

	TFile* f = new TFile("/home/mjbales/school/rdk/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root", "read");

	theHists[0] = (TH1I*) f->Get("apr11_9_19");
	theHists[1] = (TH1I*) f->Get("apr11_10_22");
	theHists[2] = (TH1I*) f->Get("apr11_10_23");
	theHists[3] = (TH1I*) f->Get("apr11_10_24");
	theHists[4] = (TH1I*) f->Get("apr11_10_26");
	theHists[5] = (TH1I*) f->Get("apr11_10_27");
	theHists[6] = (TH1I*) f->Get("apr11_10_28");
	theHists[7] = (TH1I*) f->Get("apr11_10_29");
	theHists[8] = (TH1I*) f->Get("apr11_10_30");
	theHists[9] = (TH1I*) f->Get("apr11_10_31");
	theHists[10] = (TH1I*) f->Get("apr11_10_32");
	theHists[11] = (TH1I*) f->Get("apr11_10_34");
	theHists[12] = (TH1I*) f->Get("apr11_10_35");

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(num, theHists, theLegend);

	theHists[0]->SetTitle("BNL 1400 eV Runs");
	theHists[0]->GetXaxis()->SetTitle("Channel #");
	theHists[0]->GetYaxis()->SetTitle("Counts (AU)");
	theHists[0]->GetXaxis()->SetRangeUser(0, 200);
	theHists[0]->GetYaxis()->SetRangeUser(0, 4000);

	c2->Draw();

	return c2;

}

TCanvas* plotBNLSpectra600(int num)
{
	const int numHists = 15;
	TH1* theHists[numHists];

	TFile* f = new TFile("/home/mjbales/school/rdk/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root", "read");

	theHists[0] = (TH1I*) f->Get("apr11_9_8");
	theHists[1] = (TH1I*) f->Get("apr11_9_9");
	theHists[2] = (TH1I*) f->Get("apr11_9_10");
	theHists[3] = (TH1I*) f->Get("apr11_9_35");
	theHists[4] = (TH1I*) f->Get("apr11_9_36");
	theHists[5] = (TH1I*) f->Get("apr11_9_37");
	theHists[6] = (TH1I*) f->Get("apr11_9_38");
	theHists[7] = (TH1I*) f->Get("apr11_9_39");
	theHists[8] = (TH1I*) f->Get("apr11_9_40");
	theHists[9] = (TH1I*) f->Get("apr11_9_41");
	theHists[10] = (TH1I*) f->Get("apr11_5_18");
	theHists[11] = (TH1I*) f->Get("apr11_5_8");
	theHists[12] = (TH1I*) f->Get("apr11_5_18");
	theHists[13] = (TH1I*) f->Get("apr11_6_13");
	theHists[14] = (TH1I*) f->Get("apr11_6_16");

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(num, theHists, theLegend);

	theHists[0]->SetTitle("BNL 600 eV Runs");
	theHists[0]->GetXaxis()->SetTitle("Channel #");
	theHists[0]->GetYaxis()->SetTitle("Counts (AU)");
	theHists[0]->GetXaxis()->SetRangeUser(0, 100);
	theHists[0]->GetYaxis()->SetRangeUser(0, 5000);

	c2->Draw();

	return c2;

}

TCanvas* plotBNLSpectra1066(int num)
{
	const int numHists = 6;
	TH1* theHists[numHists];

	TFile* f = new TFile("/home/mjbales/school/rdk/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root", "read");

	theHists[0] = (TH1I*) f->Get("apr11_7_26");
	theHists[1] = (TH1I*) f->Get("apr11_7_28");
	theHists[2] = (TH1I*) f->Get("apr11_8_16");
	theHists[3] = (TH1I*) f->Get("apr11_8_17");
	theHists[4] = (TH1I*) f->Get("apr11_8_19");
	theHists[5] = (TH1I*) f->Get("apr11_8_20");

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(num, theHists, theLegend);

	theHists[0]->SetTitle("BNL 1066 eV Runs");
	theHists[0]->GetXaxis()->SetTitle("Channel #");
	theHists[0]->GetYaxis()->SetTitle("Counts (AU)");
	theHists[0]->GetXaxis()->SetRangeUser(0, 200);
	theHists[0]->GetYaxis()->SetRangeUser(0, 16000);

	c2->Draw();

	return c2;

}

TCanvas* plotBNLSpectra400(int num)
{
	const int numHists = 13;
	TH1* theHists[numHists];

	TFile* f = new TFile("/home/mjbales/school/rdk/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root", "read");

	theHists[0] = (TH1I*) f->Get("apr11_5_10");
	theHists[1] = (TH1I*) f->Get("apr11_5_21");
	theHists[2] = (TH1I*) f->Get("apr11_9_44");
	theHists[3] = (TH1I*) f->Get("apr11_9_45");
	theHists[4] = (TH1I*) f->Get("apr11_10_1");
	theHists[5] = (TH1I*) f->Get("apr11_10_2");
	theHists[6] = (TH1I*) f->Get("apr11_10_3");
	theHists[7] = (TH1I*) f->Get("apr11_10_5");
	theHists[8] = (TH1I*) f->Get("apr11_10_6");
	theHists[9] = (TH1I*) f->Get("apr11_10_7");
	theHists[10] = (TH1I*) f->Get("apr11_10_8");
	theHists[11] = (TH1I*) f->Get("apr11_10_9");
	theHists[12] = (TH1I*) f->Get("apr11_10_10");

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(num, theHists, theLegend);

	theHists[0]->SetTitle("BNL 400 eV Runs");
	theHists[0]->GetXaxis()->SetTitle("Channel #");
	theHists[0]->GetYaxis()->SetTitle("Counts (AU)");
	theHists[0]->GetXaxis()->SetRangeUser(0, 60);
	theHists[0]->GetYaxis()->SetRangeUser(0, 10000);

	c2->Draw();

	return c2;

}

TCanvas* plotee0VsExp()
{
	const int numHists = 2;
	TH1* theHists[numHists];

	TFile* ee0File = new TFile("apsFile.root", "read");
	theHists[0] = (TH1D*) ee0File->Get("bZeroHist");
	theHists[1] = getExpEEBSHist();
	theHists[0]->Scale(1. / theHists[0]->Integral());
	theHists[1]->Scale(1. / theHists[1]->Integral());

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");
	theHists[1]->SetTitle("Detected Electron Kinetic Energy");
	theHists[0]->SetTitle("Initial Electron Kinetic Energy");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);

	theHists[0]->GetXaxis()->SetTitle("Electron Kinetic Energy (keV)");
	theHists[0]->GetYaxis()->SetTitle("Arbitary units");
	theHists[0]->GetXaxis()->SetRangeUser(0, 790);

	c2->Draw();

	return c2;

}

TCanvas* ploteeFVsExp()
{
	const int numHists = 2;
	TH1* theHists[numHists];

	TFile* eeFFile = new TFile("littleb_analysis_rob3.root", "read");
	theHists[0] = (TH1D*) eeFFile->Get("MC_EEBS_EP20");
	theHists[1] = getExpEEBSHist();
	theHists[0]->Scale(theHists[1]->Integral(201, 800) / theHists[0]->Integral(201, 800));
	for (int i = 0; i <= 1000; i++)
	{
		if(theHists[1]->GetBinContent(i) <= 1) theHists[1]->SetBinContent(i, 1);
		if(theHists[0]->GetBinContent(i) <= 1) theHists[0]->SetBinContent(i, 1);
	}

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");
	theHists[1]->SetTitle("Detected Electron Kinetic Energy");
	theHists[0]->SetTitle("Simulated Detector Response");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);

	theHists[0]->GetXaxis()->SetTitle("Electron Kinetic Energy (keV)");
	theHists[0]->GetYaxis()->SetTitle("Arbitary Units");
	theHists[0]->GetXaxis()->SetRangeUser(100, 790);
	theHists[1]->SetLineColor(kBlue);
	theHists[0]->SetLineColor(kGreen);

	c2->Draw();

	return c2;

}

TH1D* getExpEEBSHist()
{
	int data;
	float error;
	string s;
	ifstream itxFile;
	TH1D* theHist = new TH1D("Exp_EEBS_EP", "Exp_EEBS_EP", 1000, 0, 1000);
	theHist->Sumw2();
	TString fileName = HISTS_DIR;
	fileName += "eEnergy_v4d22++.itx";
	itxFile.open(fileName);

	for (int j = 0; j < 3; j++)
		getline(itxFile, s, '\r'); //Skip to experimental data

	int entries = 0;
	for (int i = 1; i <= 1000; i++)
	{
		getline(itxFile, s, '\r'); //Data line
		sscanf(s.data(), "\t%i\t%f", &data, &error);
		if(data >= 0)
		{
			theHist->SetBinContent(i, data);
			theHist->SetBinError(i, error);
			entries += data;
		}
	}
	theHist->SetEntries(entries);

	itxFile.close();

	return theHist;
}

//Creates files with a root tree that contains a boolean that represents whether or not a gamma event should be rejected due to the SIO2 layer
//void makeSIO2Tree(TString gammaFileName)
//{
//
//    TH1D* sio2Hist=getSiO2Hist();
//    TRandom3* ranGen = new TRandom3(0);
//    int numMultipleEntering=0;
//
//    TFile* inputFile=new TFile(gammaFileName,"read");
//    TTree* BAPDFluxTree1=(TTree*) inputFile->Get("BAPDFluxTree1");
//    TTree* BAPDFluxTree2=(TTree*) inputFile->Get("BAPDFluxTree2");
//    TTree* BAPDFluxTree3=(TTree*) inputFile->Get("BAPDFluxTree3");
//
//    if(BAPDFluxTree1 == NULL || BAPDFluxTree2 == NULL || BAPDFluxTree3 == NULL)
//    {
//        typeAnythingToContinue("ERROR CHAIN NOT FOUND!");
//        return;
//    }
//
//    int currentBAPD1=0;
//    int currentBAPD2=0;
//    int currentBAPD3=0;
//    double MomX,MomY,MomZ;
//    double KE;
//    int Event;
//    char HitCode;
//    BAPDFluxTree1->SetBranchAddress("MomX",&MomX);
//    BAPDFluxTree2->SetBranchAddress("MomX",&MomX);
//    BAPDFluxTree3->SetBranchAddress("MomX",&MomX);
//
//    BAPDFluxTree1->SetBranchAddress("MomY",&MomY);
//    BAPDFluxTree2->SetBranchAddress("MomY",&MomY);
//    BAPDFluxTree3->SetBranchAddress("MomY",&MomY);
//
//    BAPDFluxTree1->SetBranchAddress("MomZ",&MomZ);
//    BAPDFluxTree2->SetBranchAddress("MomZ",&MomZ);
//    BAPDFluxTree3->SetBranchAddress("MomZ",&MomZ);
//
//    BAPDFluxTree1->SetBranchAddress("KE",&KE);
//    BAPDFluxTree2->SetBranchAddress("KE",&KE);
//    BAPDFluxTree3->SetBranchAddress("KE",&KE);
//
//    BAPDFluxTree1->SetBranchAddress("HitCode",&HitCode);
//    BAPDFluxTree2->SetBranchAddress("HitCode",&HitCode);
//    BAPDFluxTree3->SetBranchAddress("HitCode",&HitCode);
//
//    BAPDFluxTree1->SetBranchAddress("Event",&Event);
//    BAPDFluxTree2->SetBranchAddress("Event",&Event);
//    BAPDFluxTree3->SetBranchAddress("Event",&Event);
//    TFile* outputFile;
//    outputFile=new TFile(addBeforeExtension(gammaFileName,"_SIO2"),"recreate");
//    TTree* outputTree;
//    outputTree=new TTree("SIO2Tree", "SIO2 Layer Absorbtion Tree    FALSE=ABSORBED");
//    bool SIO2Trans=true;
//    bool hitDet;
//    outputTree->Branch("SIO2Trans",&SIO2Trans,"SIO2Trans/O");
//
//    double thicknessMultiplier;
//    double transCoeff;
//    for(int i=0;i<1000000;i++)
//    {
//        SIO2Trans=true;
//        hitDet=false;
//
//
//        //BAPDFluxTree1
//        while(currentBAPD1 < BAPDFluxTree1->GetEntries())
//        {
//            BAPDFluxTree1->GetEntry(currentBAPD1);
//            if(Event > i)
//                break;
//
//            if(Event== i && HitCode == 1 && KE > .1 && KE < 30)
//            {
//                if(hitDet)
//                    numMultipleEntering++;
//                else
//                    hitDet=true;
//                thicknessMultiplier=sqrt(MomX*MomX+MomY*MomY+MomZ*MomZ)/abs(MomZ); //Due to angle of incidence
//                //Data sheet is based on 100 nm.  To achieve 35 nm, I use 35./100.
//                transCoeff=pow(sio2Hist->Interpolate(KE),(thicknessMultiplier*35./100.));
//
//                if(ranGen->Rndm() > transCoeff)
//                {
//
//                    SIO2Trans=false;
//
//                }
//            }
//            currentBAPD1++;
//        }
//
//        //BAPDFluxTree2
//        while(currentBAPD2 < BAPDFluxTree2->GetEntries())
//        {
//
//            BAPDFluxTree2->GetEntry(currentBAPD2);
//            if(Event > i)
//                break;
//
//            if(Event== i && HitCode == 1 && KE > .1 && KE < 30)
//            {
//                if(hitDet)
//                    numMultipleEntering++;
//                else
//                    hitDet=true;
//                thicknessMultiplier=sqrt(MomX*MomX+MomY*MomY+MomZ*MomZ)/abs(MomZ); //Due to angle of incidence
//                //Data sheet is based on 100 nm.  To achieve 35 nm, I use 35./100.
//                double transCoeff=pow(sio2Hist->Interpolate(KE),(thicknessMultiplier*35./100.));
//                if(ranGen->Rndm() > transCoeff)
//                {
//                    SIO2Trans=false;
//                }
//            }
//            currentBAPD2++;
//        }
//
//        //BAPDFluxTree3
//        while(currentBAPD3 < BAPDFluxTree3->GetEntries())
//        {
//            BAPDFluxTree3->GetEntry(currentBAPD3);
//            if(Event > i)
//                break;
//
//            if(Event== i && HitCode == 1 && KE > .1 && KE < 30)
//            {
//                if(hitDet)
//                    numMultipleEntering++;
//                else
//                    hitDet=true;
//                thicknessMultiplier=sqrt(MomX*MomX+MomY*MomY+MomZ*MomZ)/abs(MomZ); //Due to angle of incidence
//                //Data sheet is based on 100 nm.  To achieve 35 nm, I use 35./100.
//                double transCoeff=pow(sio2Hist->Interpolate(KE),(thicknessMultiplier*35./100.));
//                if(ranGen->Rndm() > transCoeff)
//                {
//                    SIO2Trans=false;
//
//                }
//            }
//            currentBAPD3++;
//        }
//
//        outputTree->Fill();
//        if(numMultipleEntering >1)
//            cout << numMultipleEntering << " entered the BAPDS" << endl;
//        numMultipleEntering=0;
//    }
//
//    inputFile->Close();
//    delete inputFile;
//    outputFile->Write();
//    outputFile->Close();
//    delete outputFile;
//    delete sio2Hist;
//
//}

//void makeSIO2Trees(TString firstFileName,int numFiles)
//{
//    TString currentFileName=firstFileName;
//    for(int i=0;i<numFiles;i++)
//    {
//
//        makeSIO2Tree(currentFileName);
//        currentFileName=getNextResultFileName(currentFileName);
//    }
//}

void makeBAPDCEGraphs()
{
	TH1* theHists[6];
	int energies[6] = { 1500, 900, 700, 600, 500, 350 };
	TString inpFileName, histName;
	TFile* inpFile[6];
	TString histFileName;

	for (int i = 0; i < 6; i++)
	{

		histName = "BAPD" + int2str(energies[i]) + "eV";
		theHists[i] = new TH1D(histName, histName, 110, 0., 2200);
		inpFileName = RESULTS_DIR;
		inpFileName += "BAPDCEModel_" + int2str(energies[i]) + "keV.root";
		if(!FileExists(inpFileName))
		{
			cout << inpFileName << " does not exist.  Ending" << endl;
			return;
		}
		inpFile[i] = new TFile(inpFileName, "read");

		TTree* theTree = (TTree*) inpFile[i]->Get("g");
		double BAPD3EDepTotal;
		theTree->SetBranchAddress("BAPD3EDepTotal", &BAPD3EDepTotal);

		for (int j = 0; j < theTree->GetEntries(); j++)
		{
			theTree->GetEntry(j);
			theHists[i]->Fill(BAPD3EDepTotal * 1000);
		}
		histFileName = HISTS_DIR;
		histFileName += "BAPDCEModel_" + int2str(energies[i]) + "eV.txt";
		convertTH1ToTXT(theHists[i], histFileName);
	}

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	// TCanvas* c2 = plotManyHistogramsTogether(numHists,theHists, theLegend);
	theHists[0]->SetTitle("BAPD Collection Efficiency Model;Energy Detected (keV);Arbitrary Counts");
	TCanvas* c2 = plotManyHistogramsTogether(6, theHists, theLegend);

	TString imageFileName = GRAPHS_DIR;
	imageFileName += "BAPDCEModel_110825.png";
	c2->SaveAs(imageFileName);

	for (int i = 0; i < 6; i++)
		inpFile[i]->Close();

}

void makeBAPDProgression(TString radString, bool logGraph)
{
	const int numHists = 3;
	TH1* theHists[numHists];
	const int numMCHists = 4;
	TH1* mcHists[numMCHists];
	TString HistsDir = HISTS_DIR;
	TString name;
//    name+="hist_geantliv0608_Flux.png.txt";
//    theHists[0]=getCSVHist(name,200,0,20.,'\n');
//    name=HISTS_DIR;
//    name+="hist_BareAPDw35nmSiO2.txt";
//    theHists[1]=getCSVHist(name,200,0,20.,'\n');
//    name=HISTS_DIR;

	theHists[0] = getTabSeperatedHist(HistsDir + "Exp_BAPD_Ben_ViaTom_120714.txt", 78, 0, 19.5);
	theHists[0]->Rebin(2);
	theHists[0]->SetTitle("BAPD Spectrum Evolution;Energy Detected (eV);EPG/EP Rate");
	setHistNameTitle(theHists[0], "Exp_Ben");
	theHists[0]->Scale(1. / 2.396e7);

	theHists[1] = getTabSeperatedHist(HistsDir + "Exp_BAPD_Jeff_ViaTom_120714.txt", 120, 0, 30);
	theHists[1]->Rebin(2);
	setHistNameTitle(theHists[1], "Exp_Jeff");
	theHists[1]->Scale(1. / 2.396e7);

	theHists[2] = getTabSeperatedHist(HistsDir + "EXP_Matt_JeffM1_AID1EPG_EDepGBAPD.txt", 100, 0, 50);
	setHistNameTitle(theHists[2], "Exp_Matt");

	mcHists[0] = getCSVHist(HistsDir + "hist_" + radString + "_BAPD.txt", 200, 0, 20, '\n');
	mcHists[0]->Rebin(5);
	setHistNameTitle(mcHists[0], "MC_Old");
	mcHists[0]->Scale(1. / 2.396e7);

	mcHists[1] = getTabSeperatedHist(HistsDir + "hist_RID93_BAPD.txt", 200, 0, 20);
	mcHists[1]->Rebin(5);
	setHistNameTitle(mcHists[1], "MC_New");

	mcHists[2] = getCSVHist(HistsDir + "hist_geantliv0608_Flux.png.txt", 200, 0, 20, '\n');
	mcHists[2]->Rebin(5);
	setHistNameTitle(mcHists[2], "Flux_Old");
	mcHists[2]->Scale(1. / 2.396e7);

	mcHists[3] = getTabSeperatedHist(HistsDir + "hist120909_bapdRID93.txt", 200, 0, 20);
	mcHists[3]->Rebin(5);
	setHistNameTitle(mcHists[3], "Flux_New");

	theHists[0]->GetXaxis()->SetRangeUser(0, 20);
	//theHists[0]->GetYaxis()->SetRangeUser(1,400);

	plotExpVersusMCToImage(numHists, theHists, numMCHists, mcHists, "BAPD Comparison", "BAPDEvolution_120909.png");

}

void makeCombinedBGOGraphs(TString radhighString, bool logGraph)
{
	TString f1Name = HISTS_DIR;
	f1Name += "hists_gammaEnergy_BGO_" + radhighString + "Combined.root";
	TFile f1(f1Name, "read");
	//TString f2Name=HISTS_DIR;
	// f2Name+="hists_gammaEnergy_BGO_geantliv110718wBGOPos.root";
	// TFile f2(f2Name,"read");

	const int numHists = 1;
	TH1* theHists[numHists];

	theHists[0] = (TH1D*) f1.Get("hist_BGO_M1_Combined");
	theHists[0]->SetName("Geant4");
	// theHists[1]=(TH1*) f2.Get("hist_BGO_M1");
	// theHists[1]->SetName("Geant4wBGOPos");
	theHists[0]->GetXaxis()->SetRangeUser(10, 799.99);
	theHists[0]->SetTitle("BGO Spectrum;Photon Energy Absorbed (keV);Counts per 1e6 EP Events");
	theHists[0]->GetYaxis()->SetRangeUser(10e-2, 144);

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);
	TString imageFileName = GRAPHS_DIR;
	imageFileName += "BGOMC_" + radhighString;
	if(logGraph)
	{
		c2->SetLogy();
		imageFileName += "_Log";
	}

	imageFileName += ".png";
	c2->SaveAs(imageFileName);

	TString txtFileName = HISTS_DIR;
	txtFileName += "hist_BGO_M1_wPosEff_" + radhighString + ".txt";
	convertTH1ToTXT(theHists[0], txtFileName);

}

void combineBGORanges(TString radlowString, TString radmidString, TString radhighString, TString analysisNameString, int multiplicity, int numBGOsUsed)
{
	TString genHistName;
	if(multiplicity == -1)
		genHistName = "hist_BGO_Sum";
	else if(multiplicity == -2)
		genHistName = "hist_BGO_KEFirst";
	else if(multiplicity == -3)
		genHistName = "hist_BGO_eg0";
	else if(multiplicity == -4)
		genHistName = "hist_BGO_ToF";
	else if(multiplicity == -5)
		genHistName = "hist_BGO_EKE";
	else if(multiplicity == -6)
		genHistName = "hist_BGO_Z0";
	else
		genHistName = "hist_BGO_M" + int2str(multiplicity);

	TString f1Name = HISTS_DIR;
	f1Name += "hists_gammaEnergy_BGO_radhigh_" + radhighString + "_" + analysisNameString + ".root";
	TFile f1(f1Name, "read");

	TString f2Name = HISTS_DIR;
	f2Name += "hists_gammaEnergy_BGO_radmid_" + radmidString + "_" + analysisNameString + ".root";
	;
	TFile f2(f2Name, "read");

	TFile* f3 = nullptr;
	TH1D* radLowHist;

	TH1D* radHighHist = (TH1D*) f1.Get(genHistName);
	TString highName = genHistName + "_radhigh";
	radHighHist->SetName(highName);

	int numBins = radHighHist->GetNbinsX();
	double lowExtent = getHistExtentLow(radHighHist);
	double highExtent = getHistExtentHigh(radHighHist);

	TH1D* radMidHist = (TH1D*) f2.Get(genHistName);
	TString midName = genHistName + "_radmid";
	radMidHist->SetName(midName);

	TString newFileName = HISTS_DIR;
	newFileName += genHistName + "_" + radhighString + "_" + analysisNameString + "_Combined.root";

	TFile newFile(newFileName, "recreate");
	newFile.cd();

	TString combinedName = genHistName + "_Combined";

	newFile.cd();
	TH1D* newHist = new TH1D(combinedName, "", numBins, lowExtent, highExtent);

	if(radlowString != "")
	{
		TString f3Name = HISTS_DIR;
		f3Name += "hists_gammaEnergy_BGO_radlow_" + radlowString + "_" + analysisNameString + ".root";
		;
		f3 = new TFile(f3Name, "read");
		radLowHist = (TH1D*) f3->Get(genHistName);
		TString lowName = genHistName + "_radlow";
		radLowHist->SetName(lowName);
		*newHist = *radHighHist + *radMidHist + *radLowHist;
	}
	else
	{
		*newHist = *radHighHist + *radMidHist;
	}
	newHist->SetName(combinedName);

	newFile.cd();
	newHist->Write("", TObject::kOverwrite);
	newFile.Write("", TObject::kOverwrite);
	TString txtFileName = HISTS_DIR;
	txtFileName += combinedName + "_" + radhighString + "_" + analysisNameString + ".txt";

	convertTH1ToTXT(newHist, txtFileName);
	newFile.Close();
	TString histName, newHistName;
	if(multiplicity == 1)
	{
		for (int i = 0; i < numBGOsUsed; i++)
		{
			histName = "hist_BGO_D" + int2str(i + 1);
			radHighHist = (TH1D*) f1.Get(histName);
			newHistName = histName + "_radhigh";
			radHighHist->SetName(newHistName);
			radMidHist = (TH1D*) f2.Get(histName);
			newHistName = histName + "_radmid";
			radMidHist->SetName(newHistName);

			newFile.cd();
			newHist = new TH1D("hist_BGO_M1_Combined", "", numBins, lowExtent, highExtent);
			if(radlowString != "")
			{
				radLowHist = (TH1D*) f3->Get(histName);
				newHistName = histName + "_radlow";
				radLowHist->SetName(newHistName);
				*newHist = *radHighHist + *radMidHist + *radLowHist;
			}
			else
			{
				*newHist = *radHighHist + *radMidHist;
			}
			newHistName = histName + "_Combined";
			newHist->SetName(newHistName);

			txtFileName = HISTS_DIR;
			txtFileName += histName + "_" + radhighString + "_" + analysisNameString + "_Combined.txt";

			convertTH1ToTXT(newHist, txtFileName);

		}
	}

	f1.Close();
	f2.Close();
	if(f3 != nullptr)
	{
		f3->Close();
		f3 = nullptr;
	}

}

void checkFalseTriple(TString detType, const int numDets)
{
	TString fileName;

	TString gChainString;
	int numFiles = 200;

	TH1D* mattHist = new TH1D("MattGeant4", "", 200, 0, 1000);

	fileName = RESULTS_DIR;
	fileName += "Results100421_hmg_P_Set0_g4empen11910.root";
	TChain* resultsChainP = makeResultChain(fileName, "g", numFiles);
	fileName = RESULTS_DIR;
	fileName += "Results100421_hmg_E_Set0_g4empen11910.root";
	TChain* resultsChainE = makeResultChain(fileName, "g", numFiles);

	TString lilbString = LITTLEB_DIR;
	lilbString += "MJB-Evn-hmg-n1.B3-0-20100421_lilb.root";
	TChain* lilbChain = makeEventChain(lilbString, "lilbTree", numFiles);

	double SBDEDepE, SBDEDepP, SBDTimeFirstP;
	double GDetEDepP[numDets];
	double GDetEDepE[numDets];
	int lilbcut;
	resultsChainE->SetBranchAddress("SBDEDepTotal", &SBDEDepE);
	resultsChainP->SetBranchAddress("SBDEDepTotal", &SBDEDepP);
	resultsChainP->SetBranchAddress("SBDTimeFirst", &SBDTimeFirstP);
	lilbChain->SetBranchAddress("dwcutb20", &lilbcut);

	int falseTripleDetCount[numDets];
	int detNum[numDets];

	for (int i = 0; i < numDets; i++)
	{
		gChainString = detType + int2str(i + 1) + "EDepTotal";
		resultsChainP->SetBranchAddress(gChainString, &(GDetEDepP[i]));
		resultsChainE->SetBranchAddress(gChainString, &(GDetEDepE[i]));
		falseTripleDetCount[i] = 0;
		detNum[i] = i + 1;
	}

	for (int j = 0; j < resultsChainP->GetEntries(); j++)
	{
		if(j % 1000000 == 0) cout << "File " << j / 1000000 << " started!" << endl;
		lilbChain->GetEntry(j);
		if(lilbcut == 0)
		{
			resultsChainP->GetEntry(j);
			resultsChainE->GetEntry(j);
			if(SBDEDepE < 25 || SBDEDepP < 9.3 || SBDEDepP > 30.9 || SBDTimeFirstP < 1.4e-6 || SBDTimeFirstP > 25e-6)
			{
			}
			else
			{
				for (int i = 0; i < numDets; i++)
				{
					mattHist->Fill(GDetEDepP[i]);
					mattHist->Fill(GDetEDepE[i]);
					if(GDetEDepE[i] + GDetEDepP[i] > 9.5) falseTripleDetCount[i]++;
				}
			}
		}
	}

	const int numHists = 1;
	TH1* theHists[numHists];

	theHists[0] = mattHist;
	//mattHist->Rebin(10);
	TString titleString = "False Triple " + detType + " spectrum (Three Body Events b=0);Light Output (keV?);Counts for 200 million decays";
	mattHist->SetTitle(titleString);
	mattHist->GetYaxis()->SetRange(2, 200);

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);
	TString imageFileName = GRAPHS_DIR;
	imageFileName += detType + "FalseTriple_110921";

	//c2->SetLogy();
	//imageFileName+="_Log";

	imageFileName += ".png";
	c2->SaveAs(imageFileName);

	TString histName = HISTS_DIR;
	histName += detType + "FalseTriple_110921.txt";

	convertTH1ToTXT(mattHist, histName);

	TCanvas* c3 = new TCanvas("c2", "plotManyHistCanvas", 10, 10, 1200, 1024);
	TGraph detCountGraph(numDets, detNum, falseTripleDetCount);
	detCountGraph.SetTitle("False Triple Coincidence Location of Events");
	detCountGraph.GetXaxis()->SetTitle("Detector Number");
	detCountGraph.GetYaxis()->SetTitle("Number of Events (arbitrary units)");
	detCountGraph.Draw("AB");
	imageFileName = GRAPHS_DIR;
	imageFileName += detType + "FalseTriple_110921DetLoc";
	imageFileName += ".root";
	c3->SaveAs(imageFileName);

	delete resultsChainE;
	delete resultsChainP;

}

void compareBGOSingleCrystal(TString resultsFileName, bool useLog)
{
	TString fileName = HISTS_DIR;
	fileName += "Cs137v1Rob.root";
	if(!FileExists(fileName)) return;
	TFile robFile(fileName, "read");

	TTree* robTree = (TTree*) robFile.Get("ptracTree");
	int numRobEntries = robTree->GetEntries();
	TH1D* robHist = new TH1D("RobMCNP", "", 10000, 0, 1000);
	robTree->Draw("EnergyDetector1*1000 >> RobMCNP");

	fileName = RESULTS_DIR;
	fileName += resultsFileName;
	TChain* resultsGammaChain = makeResultChain(fileName, "g", 1);

	double BGOEDepG[1];

	TString gChainString;
	for (int i = 0; i < 1; i++)
	{
		gChainString = "BGO" + int2str(i + 1) + "EDepTotal";
		resultsGammaChain->SetBranchAddress(gChainString, &(BGOEDepG[i]));
	}

	TH1D* mattHist = new TH1D("MattGeant4", "", 10000, 0, 1000);

	int numEntries = resultsGammaChain->GetEntries();
	for (int j = 0; j < numEntries; j++)
	{
		resultsGammaChain->GetEntry(j);
		for (int i = 0; i < 1; i++)
		{
			mattHist->Fill(BGOEDepG[i]);
		}
	}

	const int numHists = 2;
	TH1* theHists[numHists];

	theHists[0] = robHist;
	robHist->Rebin(10);
	//robHist->GetYaxis()->SetRangeUser(1,2000000);
	theHists[1] = mattHist;
	mattHist->Rebin(10);
	mattHist->Sumw2();
	double scaleFactor = (double) numRobEntries / (double) numEntries;
	mattHist->Scale(scaleFactor); //Rob uses 10 million events
	robHist->SetTitle("Calibration Setup Single BGO 662 keV;Light Detected (keV?);Counts (AU)");

	//double ratio = robHist->Integral(641,1000)/mattHist->Integral(641,1000);

	double robPeak = robHist->Integral(641, 1000);
	double mattPeak = mattHist->Integral(641, 1000);
	double residual = (robPeak - mattPeak) / robPeak;
	double residualError = residual * sqrt((1 / (robPeak + mattPeak)) + ((1 / robPeak)));
	cout << "rob peak %: " << robPeak / (double) numRobEntries << endl;
	cout << "matt peak % : " << mattPeak / (double) numRobEntries << endl;
	cout << "difference:" << residual << " +//- " << residualError << endl;

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);
	TString imageFileName = GRAPHS_DIR;
	imageFileName += resultsFileName;
	if(useLog)
	{
		c2->SetLogy();
		imageFileName += "_Log";
	}

	imageFileName += ".png";
	c2->SaveAs(imageFileName);

	TString histName = HISTS_DIR;
	histName += "Geant4CalSetup110922SingleBGOLiv.txt";

	convertTH1ToTXT(mattHist, histName);

	delete resultsGammaChain;
}

void compareProtonSpectrum(TString geantResultsFile, TString inpStrimString, bool useLog)
{
	const int numHists = 2;
	TH1* theHists[numHists];
	TH1D srimHist("SRIM", "SRIM", 260, 0., 26);
	theHists[1] = &srimHist;

	TString strimFileName = HISTS_DIR;
	strimFileName += inpStrimString;
	ifstream srimFile;
	srimFile.open(strimFileName);

	TString dataString;

	for (int i = 0; i < 250; i++)
	{
		getline(srimFile, dataString);
		srimHist.SetBinContent(i + 1, str2double(dataString));

	}
	TString dirTestEM5string = "/home/mjbales/school/geant/geant4/work/examples/extended/electromagnetic/TestEm5/sbdbench/";
	TString fullPathString;

	fullPathString = dirTestEM5string + geantResultsFile;
	TFile file1(fullPathString, "read");

	theHists[0] = (TH1*) file1.Get("10");
	theHists[0]->SetName("Geant4.9.4p02");

	srimHist.Sumw2();
	srimHist.Rebin(10);
	srimHist.Scale(1e-5);

	double total = 0, average = 0;
	for (int i = 1; i <= theHists[0]->GetNbinsX(); i++)
	{
		total += theHists[0]->GetBinContent(i);
		average += theHists[0]->GetBinContent(i) * theHists[0]->GetBinCenter(i);
	}
	cout << "Geant Total: " << total << "  Avg: " << average << " kev" << endl;

	total = 0;
	average = 0;
	for (int i = 1; i <= theHists[1]->GetNbinsX(); i++)
	{
		total += theHists[1]->GetBinContent(i);
		average += theHists[1]->GetBinContent(i) * theHists[1]->GetBinCenter(i);
	}
	cout << "SRIM Total: " << total << "  Avg: " << average << " kev" << endl;

	TLegend* theLegend = new TLegend(0.45, 0.70, 0.65, 0.85, nullptr, "brNDC");

	theHists[0]->SetTitle("SRIM & Geant Comparison - 25 keV Protons on 20nm Au(emOpt2lib);Energy Transmitted (keV);Events Detected per Events Incident");
	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);
	TString imageFileName = GRAPHS_DIR;
	imageFileName += geantResultsFile + "_Comparison.png";
	if(useLog)
	{
		c2->SetLogy();
		imageFileName += "_Log";
	}

	imageFileName += ".png";
	c2->SaveAs(imageFileName);

}

void compareTestEm5(TString particle, TString material, int ver1, int ver2, TString histName, bool useLog)
{
	const int numHists = 2;
	TString dirTestEM5string = "/home/mjbales/school/geant/geant4/work/examples/extended/electromagnetic/TestEm5/sbdbench/";
	TString fullPathString, file;

	fullPathString = dirTestEM5string + particle + "_" + material + "_0deg_test" + int2str(ver1) + ".root";
	TFile file1(fullPathString, "read");
	fullPathString = dirTestEM5string + particle + "_" + material + "_0deg_test" + int2str(ver2) + ".root";
	TFile file2(fullPathString, "read");

	TH1* theHists[numHists];

	theHists[0] = (TH1*) file1.Get(histName);
	theHists[1] = (TH1*) file2.Get(histName);

	TString outHistName = particle + "_" + material + "_0deg_test" + int2str(ver1);
	theHists[0]->SetName(outHistName);
	outHistName = particle + "_" + material + "_0deg_test" + int2str(ver2);
	theHists[1]->SetName(outHistName);

	theHists[0]->Rebin(2);
	theHists[1]->Rebin(2);
	if(useLog) theHists[0]->GetYaxis()->SetRangeUser(1.e-5, 1.);

	TLegend* theLegend = new TLegend(0.25, 0.70, 0.45, 0.85, nullptr, "brNDC");

	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);

	TString imageFileName = GRAPHS_DIR;
	imageFileName += "TestEmComparison/TestEm5Comp_" + particle + "_" + material + "_0deg_test" + int2str(ver1) + "_" + int2str(ver2) + "_" + histName;
	if(useLog)
	{
		c2->SetLogy();
		imageFileName += "_Log";

	}
	imageFileName += ".png";

	c2->SaveAs(imageFileName);

	c2->Close();

}

void compareAllHistsTestEm5(TString particle, TString material, int ver1, int ver2)
{
	compareTestEm5(particle, material, ver1, ver2, "1", true);
	compareTestEm5(particle, material, ver1, ver2, "10", true);
	compareTestEm5(particle, material, ver1, ver2, "20", false);
	compareTestEm5(particle, material, ver1, ver2, "30", false);
	compareTestEm5(particle, material, ver1, ver2, "40", false);
}

TH1D* makeEPProtonTOF(int numFiles, TString protonString, TString electronString, TString additionalCuts, int firstFile, bool useOldGeantData)
{

	TChain* resultsEChain;
	TChain* resultsPChain;
	TChain* littleBChain;

	TString eventString = "100421";

	TString fileNameE = RESULTS_DIR;
	fileNameE += "Results" + eventString + "_hmg_E_Set" + int2str(firstFile) + "_" + electronString + ".root";
	TString fileNameP = RESULTS_DIR;
	fileNameP += "Results" + eventString + "_hmg_P_Set" + int2str(firstFile) + "_" + protonString + ".root";
	TString fileNamelilb = LITTLEB_DIR;
	fileNamelilb += "MJB-Evn-hmg-n1.B3-" + int2str(firstFile) + "-20100421_lilb.root";

	resultsEChain = makeResultChain(fileNameE, "g", numFiles);
	resultsPChain = makeResultChain(fileNameP, "g", numFiles);
	littleBChain = makeEventChain(fileNamelilb, "lilbTree", numFiles);
	resultsPChain->AddFriend(resultsEChain, "ge");
	resultsPChain->AddFriend("lilbTree");

	TH1D* mcToFHist = new TH1D("mcEPToF", "mcEPToF", 41, 2.02e-6, 24.98e-6); //Binning to match exp

	//standard
	TString cutString, drawString;

	if(useOldGeantData)
	{
		drawString = "g.SBDTimeFirst*1e-9 >> mcEPToF";
		cutString = "ge.SBDEDepTotal > " + d2str(ELECTRON_MIN_ENERGY) + " && g.SBDEDepTotal > " + d2str(PROTON_MIN_ENERGY) + " && g.SBDEDepTotal < " + d2str(PROTON_MAX_ENERGY) + " && g.SBDTimeFirst > " + d2str(PROTON_ARRIVAL_MIN * 1.e9) + " && g.SBDTimeFirst < " + d2str(PROTON_ARRIVAL_MAX * 1.e9) + " && lilbTree.dwcutb20 == 0";
	}
	else
	{
		drawString = "g.SBDTimeFirst >> mcEPToF";
		cutString = "ge.SBDEDepTotal > " + d2str(ELECTRON_MIN_ENERGY) + " && g.SBDEDepTotal > " + d2str(PROTON_MIN_ENERGY) + " && g.SBDEDepTotal < " + d2str(PROTON_MAX_ENERGY) + " && g.SBDTimeFirst > " + d2str(PROTON_ARRIVAL_MIN) + " && g.SBDTimeFirst < " + d2str(PROTON_ARRIVAL_MAX) + " && lilbTree.dwcutb20 == 0";
	}

	//nonstandard
	//TString cutString = "ge.SBDEDepTotal > 400 && ge.SBDEDepTotal < 800 && g.SBDEDepTotal > "+ d2str(PROTON_MIN_ENERGY)+" && g.SBDEDepTotal < "+ d2str(PROTON_MAX_ENERGY)+" && g.SBDTimeFirst > "+d2str(PROTON_ARRIVAL_MIN*1.e9)+" && g.SBDTimeFirst < "+ d2str(PROTON_ARRIVAL_MAX*1.e9)+" && lilbTree.dwcutb20 == 0";

	if(additionalCuts != "") cutString += " && " + additionalCuts;

	cout << "Making ToF Hist with cutstring: " << cutString << endl;

	resultsPChain->Draw(drawString, cutString);

	delete resultsPChain;
	delete resultsEChain;
	delete littleBChain;

	mcToFHist->Sumw2();

	return mcToFHist;

}

TH1D* makeEPProtonTOFOldMRK(int numFiles, int firstFile, TString dateCode)
{

	TChain* resultsChain;

	TString fileName = RESULTS_DIR;
	fileName += "Results" + dateCode + "_hmg_EP_Set0.root";

	resultsChain = makeResultChain(fileName, "r", numFiles);

	TString fileNamelilb = LITTLEB_DIR;
	fileNamelilb += "MJB-Evn-hmg-n1.B3-" + int2str(firstFile) + "-20100421_lilb.root";
	TChain* littleBChain = makeEventChain(fileNamelilb, "lilbTree", numFiles);

	double vxpf, vypf, vzpf, vxef, vyef, vzef, tofp;
	int hitcodep, hitcodee, dwcutb20;

	//Branches for Protons

	resultsChain->SetBranchAddress("vxpf", &vxpf);
	resultsChain->SetBranchAddress("vypf", &vypf);
	resultsChain->SetBranchAddress("vzpf", &vzpf);
	resultsChain->SetBranchAddress("tofp", &tofp);
	resultsChain->SetBranchAddress("hitcodep", &hitcodep);

	//Branches for Electrons
	resultsChain->SetBranchAddress("vxef", &vxef);
	resultsChain->SetBranchAddress("vyef", &vyef);
	resultsChain->SetBranchAddress("vzef", &vzef);
	resultsChain->SetBranchAddress("hitcodee", &hitcodee);
	littleBChain->SetBranchAddress("dwcutb20", &dwcutb20);

	double eef, epf;

	TH1D* mcToFHist = new TH1D("mcEPToF", "mcEPToF", 41, 2.02e-6, 24.98e-6); //Binning to match exp

	for (int i = 0; i < resultsChain->GetEntries(); i++)
	{

		littleBChain->GetEntry(i);

		//standard
		if(dwcutb20 == 0)
		{
			resultsChain->GetEntry(i);
			if(hitcodee == HIT_DETECTOR && hitcodep == HIT_DETECTOR)
			{
				eef = convertVelocityToKE(vxef, vyef, vzef, EMASSE);
				epf = convertVelocityToKE(vxpf, vypf, vzpf, PMASSE);
				if(eef > 50 && eef < 500 && epf > 7 && epf < 31) mcToFHist->Fill(tofp);
			}
		}

	}

	delete resultsChain;
	delete littleBChain;

	mcToFHist->Sumw2();
	return mcToFHist;

}

void makeEPGToFComparison(TString expToFHistFileName, TString midString, TString highString, TString tofNameString, int multiplicity)
{
	const int numHists = 2;
	TH1* theHists[numHists];

	TString expEPGHistName = HISTS_DIR;
	expEPGHistName += expToFHistFileName;
	if(!FileExists(expEPGHistName))
	{
		typeAnythingToContinue(expEPGHistName + " does not exist!");
	}
	theHists[0] = getTabSeperatedHist(expEPGHistName, 1000, 0, 1000 * 40e-9);
	theHists[0]->SetTitle("Proton Time of Flight;Time of Flight (s);Counts for Experiment Run Time");
	theHists[0]->SetName("RoughExpEPGHist");
	theHists[0]->Sumw2();

	TString mcEPGHistName = HISTS_DIR;
	mcEPGHistName += "hist_tof_BGO_M" + int2str(multiplicity) + "_radmid_" + midString + "_" + tofNameString + ".txt";
	if(!FileExists(mcEPGHistName))
	{
		typeAnythingToContinue(mcEPGHistName + " does not exist!");
	}
	TH1D* midHist = getTabSeperatedHist(mcEPGHistName, 1000, 0, 1000 * 40e-9);
	mcEPGHistName = HISTS_DIR;
	mcEPGHistName += "hist_tof_BGO_M" + int2str(multiplicity) + "_radhigh_" + highString + "_" + tofNameString + ".txt";
	if(!FileExists(mcEPGHistName))
	{
		typeAnythingToContinue(mcEPGHistName + " does not exist!");
	}
	TH1D* highHist = getTabSeperatedHist(mcEPGHistName, 1000, 0, 1000 * 40e-9);

	theHists[1] = new TH1D("mcEPHistwCut", "mcEPHistwCut", 1000, 0, 1000 * 40e-9);
	*(TH1D*) theHists[1] = *midHist + *highHist;
	theHists[1]->SetName("mcEPGHist");
	theHists[1]->Sumw2();

	for (int i = 0; i < numHists; i++)
	{

		theHists[i]->Rebin(5);

	}

	theHists[0]->GetXaxis()->SetRangeUser(0., 25e-6);
	double topYRange;
	if(tofNameString == "2us")
		topYRange = 4000;
	else if(tofNameString == "4us")
		topYRange = 3000;
	else
		topYRange = 7000;
	topYRange = topYRange * pow(.1, multiplicity - 1);  //Rough Graph scaling
	theHists[0]->GetYaxis()->SetRangeUser(0, topYRange);

	TLegend* theLegend = new TLegend(0.6620603, 0.708042, 0.8505025, 0.8391608, nullptr, "brNDC");

	TCanvas* c2 = plotManyHistogramsTogether(numHists, theHists, theLegend);

	TString imageName = GRAPHS_DIR;
	imageName += "ToFEPGComparison_" + highString + "_" + tofNameString + ".png";
	c2->SaveAs(imageName);

	TString textName = HISTS_DIR;
	textName += "hist_EPG_TOF_M" + int2str(multiplicity) + "_Combined_" + highString + "_" + tofNameString + ".txt";
	convertTH1ToTXT(theHists[1], textName);

}

TH1D* makeEPHist(int numFiles, TString protonString, TString electronString, TString gammaString, TString drawString, TString cutsString, int numBins, double bottomHist, double topHist)
{
	/*
	 TChain* resultsEChain;
	 TChain* resultsPChain;
	 TChain* littleBChain;
	 TChain* eventsChain;

	 TString eventString="100421";


	 TString fileNameE=RESULTS_DIR;
	 fileNameE+="Results"+eventString+"_hmg_E_Set0_"+electronString+".root";
	 TString fileNameP=RESULTS_DIR;
	 fileNameP+="Results"+eventString+"_hmg_P_Set0_"+protonString+".root";
	 TString fileNameEvents=EVENTS_DIR;
	 fileNameEvents+="MJB-Evn-hmg-n1.B3-0-20"+eventString+".root";
	 TString fileNamelilb=LITTLEB_DIR;
	 fileNamelilb+="MJB-Evn-hmg-n1.B3-0-20"+eventString+"_lilb.root";



	 resultsEChain = makeResultChain(fileNameE,"g",numFiles);
	 resultsPChain = makeResultChain(fileNameP,"g",numFiles);
	 eventsChain = makeEventChain(fileNameEvents,"t",numFiles);
	 littleBChain = makeEventChain(fileNamelilb,"lilbTree",numFiles);
	 eventsChain->AddFriend(resultsEChain,"ge");
	 eventsChain->AddFriend("lilbTree");
	 eventsChain->AddFriend(resultsPChain,"gp");*/

	TChain* resultsEChain;
	TChain* resultsPChain;
	TChain* littleBChain;

	TString eventString = "100421";

	TString fileNameE = RESULTS_DIR;
	fileNameE += "Results" + eventString + "_hmg_E_Set0_" + electronString + ".root";
	TString fileNameP = RESULTS_DIR;
	fileNameP += "Results" + eventString + "_hmg_P_Set0_" + protonString + ".root";
	TString fileNamelilb = LITTLEB_DIR;
	fileNamelilb += "MJB-Evn-hmg-n1.B3-0-20100421_lilb.root";

	resultsEChain = makeResultChain(fileNameE, "g", numFiles);
	resultsPChain = makeResultChain(fileNameP, "g", numFiles);
	littleBChain = makeEventChain(fileNamelilb, "lilbTree", numFiles);
	resultsPChain->AddFriend(resultsEChain, "ge");
	resultsPChain->AddFriend("lilbTree");

	TString fileNameEvents = EVENTS_DIR;
	fileNameEvents += "MJB-Evn-hmg-n1.B3-0-20" + eventString + ".root";
	TChain* eventsChain = makeEventChain(fileNameEvents, "t", numFiles);
	resultsPChain->AddFriend("t");

	TString histName = "mcEPHist";
	TH1D* mcToFHist = new TH1D(histName, histName, numBins, bottomHist, topHist); //Binning to match exp

	TString completeDrawString = drawString + " >> " + histName;

	cout << "Complete Draw String: " << completeDrawString << endl;
	cout << "Complete Cuts String: " << cutsString << endl;
	resultsPChain->Draw(completeDrawString, cutsString);

	delete resultsEChain;
	delete resultsPChain;
	delete littleBChain;
	delete eventsChain;

	return mcToFHist;

}

//Normalized to 511 keV
//If rootRanGen=NULL it outputs the mean
double applyEDepModel(double e_dep, int detectorNumber, TRandom3* rootRanGen, EDepModel eDepModel, TGraph* responseGraph, int modelNumber) //Presumes det number is 1-12
{

	if(detectorNumber < 1 || detectorNumber > 12)
	{
		cout << "Error in applyEDepModel, Outside Detector Range" << endl;
		return 0;
	}

	double e_ref = 662;

	double Ak = 0;
	if(eDepModel == KHODYUK_PARAM_MODEL_WGAUSS || eDepModel == KHODYUK_PARAM_MODEL_WGAMMA || eDepModel == KHODYUK_PARAM_MODEL_NORES || eDepModel == KHODYUK_PARAM_MODEL_DBLGAUSS)
	{
		Ak = e_dep * evalMattKhodyukParam(e_dep);
	}
	else if(eDepModel == ALLEXP_PARAM_MODEL || eDepModel == ALLEXP_PARAM_MODEL_HIGH || eDepModel == ALLEXP_PARAM_MODEL_LOW)
	{
		if(modelNumber == 1)
		{
			Ak = e_dep * evalMattAllExpParamModel(e_dep, eDepModel);
		}
		else if(modelNumber == 2)
		{
			Ak = e_dep * evalMattAllExpParamModel2(e_dep, eDepModel);
		}
		else if(modelNumber == 3)
		{
			Ak = e_dep * evalMattAllExpParamModel3(e_dep, eDepModel);
		}
		else if(modelNumber == 4)
		{
			Ak = e_dep * evalMattAllExpParamModel4(e_dep, eDepModel);
		}
		else if(modelNumber == 5)
		{
			Ak = e_dep * evalMattAllExpParamModel5(e_dep, eDepModel);
		}
	}
	else if(eDepModel == MATT_EXP_EDEP)
	{
		//Fit to data from MOSZYNSKI et al.: INTRINSIC ENERGY RESOLUTION AND LIGHT YIELD NONPROPORTIONALITY OF BGO
		//IEEE TRANSACTIONS ON NUCLEAR SCIENCE, VOL. 51, NO. 3, JUNE 2004
		const double A1 = -.44704;
		const double t1 = 17.79;
		const double A2 = -.115223;
		const double t2 = 322;
		const double y0 = 1.0235;
		Ak = e_dep * (A1 * exp(-e_dep / t1) + A2 * exp(-e_dep / t2) + y0);
	}
	else if(eDepModel == KHODYUK_EDEP_MODEL || eDepModel == KHODYUK_EDEP_MODEL_DOUBLE || eDepModel == KHODYUK_EDEP_MODEL_HALF || eDepModel == KHODYUK_EDEP_MODEL_WGAMMA)
	{
		if(e_dep > 100)
		{
			const double A1 = -.285656766;
			const double t1 = -.008672992;
			const double y0 = 1.000916925;

			Ak = e_dep * (A1 * exp(e_dep * t1) + y0);
		}
		else
		{
			Ak = e_dep * responseGraph->Eval(e_dep);
		}

	}
	else if(eDepModel == MATT_EDGE_EDEP_MODEL)
	{
		if(e_dep > 100)
		{
			const double A1 = -0.3904193963;
			const double t1 = -.008672992;
			const double y0 = 1.000916925;

			Ak = e_dep * (A1 * exp(e_dep * t1) + y0);
		}
		else
		{
			Ak = e_dep * responseGraph->Eval(e_dep);
		}

	}
	else if(eDepModel == KEVIN_GAMMA_EDEP)
	{
		Ak = KevinGammaConstx3[detectorNumber - 1] * e_dep + KevinGammaConstx4[detectorNumber - 1] * pow(e_dep / e_ref, 2);

		Ak = KevinGammaCal_Fac[detectorNumber - 1] * Ak;
	}
	else if(eDepModel == KEVIN_GAUS_EDEP)
	{
		Ak = KevinGausConstx3[detectorNumber - 1] * e_dep + KevinGausConstx4[detectorNumber - 1] * pow(e_dep / e_ref, 2);

		Ak = KevinGausCal_Fac[detectorNumber - 1] * Ak;
	}
	else
	{
		Ak = e_dep;
	}

	if(rootRanGen == nullptr) return Ak;

	double FWHM = 2 * sqrt(2 * log(2));
	double sige = 0.1 / FWHM * sqrt(e_ref / e_dep) * e_dep;

	double vare, output;

	if(eDepModel == KEVIN_GAMMA_EDEP || eDepModel == KHODYUK_EDEP_MODEL_WGAMMA || eDepModel == KHODYUK_PARAM_MODEL_WGAMMA)
	{
		vare = KevinGammaConstx1[detectorNumber - 1] * sige * sige + KevinGammaConstx2[detectorNumber - 1] * pow(e_dep / e_ref, 2);
		vare = max(vare, 1e-14);
		sige = sqrt(vare);
		sige = KevinGammaCal_Fac[detectorNumber - 1] * sige;

		double rangeLow = Ak - sige * 5;
		double rangeHigh = Ak + sige * 5;
		if(rangeLow < 0)
		{
			rangeLow = 0;
		}

		TF1 gammaFunc("gammaFunc", "TMath::GammaDist(x,[0],0,[1])", rangeLow, rangeHigh);
		double shape = Ak * Ak / (sige * sige);
		double scale = sige * sige / Ak;
		gammaFunc.SetParameters(shape, scale);
		output = gammaFunc.GetRandom();

	}
	else
	{
		vare = KevinGausConstx1[detectorNumber - 1] * sige * sige + KevinGausConstx2[detectorNumber - 1] * pow(e_dep / e_ref, 2);
		vare = max(vare, 1e-14);
		sige = sqrt(vare);
		sige = KevinGausCal_Fac[detectorNumber - 1] * sige;

		if(e_dep < 20)
		{
			sige = 1.1139958 * sqrt(e_dep);
		}

		if(eDepModel == KHODYUK_EDEP_MODEL_DOUBLE || eDepModel == KHODYUK_PARAM_MODEL_DBLRES)
		{
			sige *= 2;
		}
		else if(eDepModel == KHODYUK_EDEP_MODEL_HALF)
		{
			sige *= .5;
		}

		output = rootRanGen->Gaus(Ak, sige);

		if(output < 0 && (eDepModel == KHODYUK_PARAM_MODEL_WGAUSS || eDepModel == KHODYUK_PARAM_MODEL_DBLRES || eDepModel == ALLEXP_PARAM_MODEL || eDepModel == ALLEXP_PARAM_MODEL_HIGH || eDepModel == ALLEXP_PARAM_MODEL_LOW))
		{
			output = 0;
		}

		if(eDepModel == KHODYUK_PARAM_MODEL_DBLGAUSS && output < Ak)
		{
			output = 2. * Ak - output;
		}

		if(eDepModel == KHODYUK_PARAM_MODEL_NORES)
		{
			output = Ak;
		}
	}
//    else
//    {
//
//        output=rootRanGen->PoissonD(Ak*1.2409866424);
//
//
//    }

	//cout << "Mean: " << Ak << "  Sigma: " << sige << endl;

	//Models for testing dependence on energy resolution

	return output;

}

void makeMultiplicityFile(TString firstGammaFile, int numEPGFiles) //includes Kevin's Energy Resolution
{
	TRandom3* rootRanGen = new TRandom3(0);
	TString nextResultFileName = RESULTS_DIR;
	nextResultFileName += firstGammaFile;
	TString nextResultMultiFileName = RESULTS_DIR;
	nextResultMultiFileName += addBeforeExtension(firstGammaFile, "_Multi");
	TFile* resultsMultiFile;
	TFile* resultsFile;
	TString fileNameWDir;
	TTree* resultsTree;
	TTree* resultsMultiTree;
	double BGOEDepG[12];
	TString gChainString;
	int multiplicity;
	double BGOEDepGMulti;

	for (int k = 0; k < numEPGFiles; k++)
	{

		resultsFile = new TFile(nextResultFileName, "READ");
		resultsTree = (TTree*) resultsFile->Get("g");

		resultsMultiFile = new TFile(nextResultMultiFileName, "RECREATE");
		resultsMultiFile->cd();
		resultsMultiTree = new TTree("resultsMultiTree", "Results Tree containing BGO multiplicty based data");

		for (int i = 0; i < 12; i++)
		{
			gChainString = "BGO" + int2str(i + 1) + "EDepTotal";
			resultsTree->SetBranchAddress(gChainString, &(BGOEDepG[i]));
		}

		resultsMultiTree->Branch("multiplicity", &multiplicity, "multiplicity/I");
		resultsMultiTree->Branch("BGOEDepGMulti", &BGOEDepGMulti, "BGOEDepGMulti/D");

		cout << "Processing " << nextResultMultiFileName;
		for (int j = 0; j < resultsTree->GetEntries(); j++)
		{
			resultsTree->GetEntry(j);

			BGOEDepGMulti = 0;
			multiplicity = 0;

			for (int i = 0; i < 12; i++)
			{
				if(BGOEDepG[i] > 0.1) BGOEDepG[i] = applyEDepModel(BGOEDepG[i], i + 1, rootRanGen);
				if(BGOEDepG[i] < BGO_MULTIPLICITY_CUT_ENERGY)
					BGOEDepG[i] = 0;
				else
					multiplicity++;
			}

			for (int i = 0; i < 12; i++)
			{
				if(BGOEDepG[i] > BGO_MULTIPLICITY_CUT_ENERGY)
				{
					BGOEDepGMulti = BGOEDepG[i];
					break;
				}
			}

			resultsMultiTree->Fill();

			if(j % 100000 == 0 && j != 0) cout << " . ";

		}
		cout << " completed!" << endl;

		nextResultFileName = getNextResultFileName(nextResultFileName);
		nextResultMultiFileName = getNextResultFileName(nextResultMultiFileName);
		resultsFile->Close();
		resultsMultiFile->Write();
		resultsMultiFile->Close();
		delete resultsMultiFile;
		delete resultsFile;

	}

}

void makeSingleBGOSpectrum(TString filePath, TString histName, bool useKevinResolution)
{
	TRandom3* rootRanGen = new TRandom3(0);
	TFile resultFile(filePath, "read");
	TTree* g = (TTree*) resultFile.Get("g");
	double BGO1EDepTotal;
	g->SetBranchAddress("BGO1EDepTotal", &BGO1EDepTotal);

	TH1D* theHist = new TH1D(histName, histName, 1000, 0., 1000.);

	for (int i = 0; i < g->GetEntries(); i++)
	{
		g->GetEntry(i);
		if(useKevinResolution) BGO1EDepTotal = applyEDepModel(BGO1EDepTotal, 1, rootRanGen);

		theHist->Fill(BGO1EDepTotal);

	}

	TString histPath = HISTS_DIR;
	histPath += histName;
	if(useKevinResolution) histPath += "_KevRes";
	histPath += ".txt";
	convertTH1ToTXT(theHist, histPath);
	resultFile.Close();
}

void makeBGOSim3DHists(int numEPGFiles, TString protonMidString, TString electronMidString, TString gammaMidString, TString protonHighString, TString electronHighString, TString gammaHighString, double epPer3BodyEvent, EDepModel eDepModel)
{
	const int numBGOs = 12;
	TRandom3* rootRanGen = new TRandom3(0);

	AnalysisCuts theCuts = { 50, 800., 2.02e-6, 24.98e-6, 10., 800. };

	double binCountEDepE = 200;
	double binLowEDepE = 0;
	double binHighEDepE = 1000;

	double binCountToFP = 41;
	double binLowToFP = 2.02e-6;
	double binHighToFP = 24.98e-6;

	double binCountEDepG = 200;
	double binLowEDepG = 0;
	double binHighEDepG = 1000;

	double branchingRatio = 0;

	TString histName, histTitle;

	const int numGammaLevels = 2;

	TH3D* BGOMult3DHists[numGammaLevels][numBGOs];
	TH3D* BGOSingles3DHists[numGammaLevels][numBGOs];

	TString gammaLevelString, eventDateString;
	eventDateString = "120208";
	TString protonString, electronString, gammaString;
	for (int gammaLevelNum = 0; gammaLevelNum < 2; gammaLevelNum++)
	{
		if(gammaLevelNum == 0)
		{
			gammaLevelString = "radmid";
			protonString = protonMidString;
			electronString = electronMidString;
			gammaString = gammaMidString;
		}
		else
		{
			gammaLevelString = "radhigh";
			protonString = protonHighString;
			electronString = electronHighString;
			gammaString = gammaHighString;
		}

		for (int j = 0; j < numBGOs; j++)
		{

			histName = "BGO3D" + gammaLevelString + "Mult" + int2str(j + 1) + "hist";
			histTitle = histName + " SBDEDepE,ToFP,BGOEDepG";
			BGOMult3DHists[gammaLevelNum][j] = new TH3D(histName, histTitle, binCountEDepE, binLowEDepE, binHighEDepE, binCountToFP, binLowToFP, binHighToFP, binCountEDepG, binLowEDepG, binHighEDepG);
			BGOMult3DHists[gammaLevelNum][j]->Sumw2();

			histName = "BGO3D" + gammaLevelString + "Singles" + int2str(j + 1) + "hist";
			histTitle = histName + " SBDEDepE,ToFP,BGOEDepG";
			BGOSingles3DHists[gammaLevelNum][j] = new TH3D(histName, histTitle, binCountEDepE, binLowEDepE, binHighEDepE, binCountToFP, binLowToFP, binHighToFP, binCountEDepG, binLowEDepG, binHighEDepG);
			BGOSingles3DHists[gammaLevelNum][j]->Sumw2();
		}

		histName = "BGO" + gammaLevelString + "GAngleHist";
		TH1D gAngleHist(histName, histName, 90, -90, 90);

		histName = "BGO" + gammaLevelString + "PAngleHist";
		TH1D pAngleHist(histName, histName, 90, -90, 90);

		histName = "BGO" + gammaLevelString + "EAngleHist";
		TH1D eAngleHist(histName, histName, 90, -90, 90);

		TString eventsPath = EVENTS_DIR;
		eventsPath += "MJB-Evn-" + gammaLevelString + "-n1.B4-0-20110711.root";
		TChain* eventsChain = makeEventChain(eventsPath, "t", numEPGFiles);

		double gAngle, pAngle, eAngle;
		double mxe0, mye0, mze0, mxp0, myp0, mzp0, mxg0, myg0, mzg0;

		eventsChain->SetBranchAddress("mxe0", &mxe0);
		eventsChain->SetBranchAddress("mye0", &mye0);
		eventsChain->SetBranchAddress("mze0", &mze0);

		eventsChain->SetBranchAddress("mxp0", &mxp0);
		eventsChain->SetBranchAddress("myp0", &myp0);
		eventsChain->SetBranchAddress("mzp0", &mzp0);

		eventsChain->SetBranchAddress("mxg0", &mxg0);
		eventsChain->SetBranchAddress("myg0", &myg0);
		eventsChain->SetBranchAddress("mzg0", &mzg0);

		TChain* resultsEChain;
		TChain* resultsPChain;
		TChain* resultsGammaChain;

		TString resultsEFileName = RESULTS_DIR;
		resultsEFileName += "Results" + eventDateString + "_" + gammaLevelString + "_E_Set0_" + electronString + ".root";

		TString resultsPFileName = RESULTS_DIR;
		resultsPFileName += "Results" + eventDateString + "_" + gammaLevelString + "_P_Set0_" + protonString + ".root";

		TString resultsGammaFileName = RESULTS_DIR;
		resultsGammaFileName += "Results" + eventDateString + "_" + gammaLevelString + "_G_Set0_" + gammaString + ".root";

		resultsEChain = makeResultChain(resultsEFileName, "g", numEPGFiles);
		resultsPChain = makeResultChain(resultsPFileName, "g", numEPGFiles);
		resultsGammaChain = makeResultChain(resultsGammaFileName, "g", numEPGFiles);

		double SBDEDepE, SBDEDepP, SBDToFP;
		double BGOEDep[numBGOs];
		double BGOKEFirst;
		resultsEChain->SetBranchAddress("SBDEDepTotal", &SBDEDepE);
		resultsPChain->SetBranchAddress("SBDEDepTotal", &SBDEDepP);
		resultsPChain->SetBranchAddress("SBDTimeFirst", &SBDToFP);

		TString gChainString;
		for (int i = 0; i < numBGOs; i++)
		{
			gChainString = "BGO" + int2str(i + 1) + "EDepTotal";
			resultsGammaChain->SetBranchAddress(gChainString, &(BGOEDep[i]));

		}

		if(eDepModel == KE_FIRST_MODEL || eDepModel == KE_FIRST_NOEP_MODEL)
		{
			resultsGammaChain->SetBranchAddress("BGOKEFirst", &BGOKEFirst);
		}

		int BGOMultiplicity;
		double BGOEDepGTotal;

//        int numEnergyViolators=0;

		//Check to see if model requests no EP
		bool checkForEP = true;

		if(eDepModel == KE_FIRST_NOEP_MODEL || eDepModel == KEVIN_GAUS_EDEP_NOEP)
		{
			checkForEP = false;
		}

		//Begin cutting/counting/histogram filling
		for (int j = 0; j < resultsGammaChain->GetEntries(); j++)
		{

			resultsEChain->GetEntry(j);
			resultsPChain->GetEntry(j);
			resultsGammaChain->GetEntry(j);
			eventsChain->GetEntry(j);

			BGOMultiplicity = 0;
			BGOEDepGTotal = 0;

			if(SBDEDepE > 0) SBDEDepE += rootRanGen->Gaus(0.0, SBD_ENERGY_SIGMA);

			//First is it an EP hit
			if(!checkForEP || (SBDToFP > theCuts.ToFPLow && SBDToFP < theCuts.ToFPHigh && SBDEDepP > PROTON_MIN_ENERGY && SBDEDepP < PROTON_MAX_ENERGY && SBDEDepE > theCuts.EDepELow && SBDEDepE < theCuts.EDepEHigh))
			{

				//BGO Loop
				for (int i = 0; i < numBGOs; i++)
				{
					if(BGOEDep[i] > 0 && eDepModel != NO_EDEP_MODEL)
					{
						BGOEDep[i] = applyEDepModel(BGOEDep[i], i + 1, rootRanGen, eDepModel);
					}

					if(eDepModel == KE_FIRST_MODEL)
					{
						BGOEDep[i] = 0;
					}

					if(BGOEDep[i] > theCuts.EDepGLow && BGOEDep[i] < theCuts.EDepGHigh)
					{

						BGOMultiplicity++;
						BGOEDepGTotal += BGOEDep[i];
					}
					else
						BGOEDep[i] = 0;

					BGOSingles3DHists[gammaLevelNum][i]->Fill(SBDEDepE, SBDToFP, BGOEDep[i]);

				}

				if(!checkForEP)
				{
					SBDEDepE = 300; //to make sure it survies a normal set of cuts when plotting
					SBDToFP = 4e-6;
				}

				if(eDepModel == KE_FIRST_MODEL || eDepModel == KE_FIRST_NOEP_MODEL)
				{
					if(BGOKEFirst > 10)
					{
						BGOEDepGTotal = BGOKEFirst;
						BGOMultiplicity = 1;
					}
					else
						BGOEDepGTotal = 0;
				}

				if(BGOEDepGTotal > theCuts.EDepGLow && BGOEDepGTotal < theCuts.EDepGHigh && BGOMultiplicity > 0)
				{
					BGOMult3DHists[gammaLevelNum][BGOMultiplicity - 1]->Fill(SBDEDepE, SBDToFP, BGOEDepGTotal);
					if(BGOMultiplicity == 1)
					{
						gAngle = atan2(sqrt(mxg0 * mxg0 + myg0 * myg0), mzg0) * 180. / PI;

						pAngle = atan2(sqrt(mxp0 * mxp0 + myp0 * myp0), mzp0) * 180. / PI;

						eAngle = atan2(sqrt(mxe0 * mxe0 + mye0 * mye0), mze0) * 180. / PI;

						gAngleHist.Fill(gAngle);
						eAngleHist.Fill(eAngle);
						pAngleHist.Fill(pAngle);
					}
				}

			}
		}

		delete resultsGammaChain;
		delete resultsEChain;
		delete resultsPChain;
	}

	TString histFilePath = HISTS_DIR;
	histFilePath += "mc3DHists_BGO_EPG_" + gammaHighString;
	if(eDepModel == MATT_EXP_EDEP) histFilePath += "_EXPEDEP";
	if(eDepModel == KEVIN_GAUS_EDEP_NOEP) histFilePath += "_NOEP";
	if(eDepModel == KE_FIRST_MODEL) histFilePath += "_KEFIRST";
	if(eDepModel == KE_FIRST_NOEP_MODEL) histFilePath += "_KEFIRSTNOEP";
	histFilePath += ".root";

	TFile* histFile = new TFile(histFilePath, "recreate");

	histFile->cd();

	TH3D* outMult3DHist[numBGOs];
	TH3D* outSingles3DHist[numBGOs];

	//Scale based on branching ratio and event count
	double scaleFactor;
	for (int gammaLevelNum = 0; gammaLevelNum < 2; gammaLevelNum++)
	{

		if(gammaLevelNum == 0)
		{
			branchingRatio = 0.00318548;
		}
		else
		{
			branchingRatio = 0.00252823;
		}
		scaleFactor = branchingRatio / (epPer3BodyEvent * numEPGFiles * 1000000);
		if(eDepModel == KEVIN_GAUS_EDEP_NOEP) scaleFactor *= epPer3BodyEvent;
		for (int i = 0; i < numBGOs; i++)
		{
			BGOMult3DHists[gammaLevelNum][i]->Scale(scaleFactor);
			BGOSingles3DHists[gammaLevelNum][i]->Scale(scaleFactor);
		}

	}

	for (int i = 0; i < numBGOs; i++)
	{

		histFile->cd();

		histName = "BGO3DMult" + int2str(i + 1) + "hist";
		histTitle = histName + " SBDEDepE,ToFP,BGOEDepG;Electron Energy (keV);Proton Time of Flight (s);BGO Energy Deposited (keV)";
		outMult3DHist[i] = new TH3D(histName, histTitle, binCountEDepE, binLowEDepE, binHighEDepE, binCountToFP, binLowToFP, binHighToFP, binCountEDepG, binLowEDepG, binHighEDepG);
		*(outMult3DHist[i]) = *(BGOMult3DHists[0][i]) + *(BGOMult3DHists[1][i]);
		outMult3DHist[i]->SetName(histName);
		outMult3DHist[i]->SetTitle(histTitle);

		histName = "BGO3DSingles" + int2str(i + 1) + "hist";
		histTitle = histName + " SBDEDepE,ToFP,BGOEDepG;Electron Energy (keV);Proton Time of Flight (s);BGO Energy Deposited (keV)";
		outSingles3DHist[i] = new TH3D(histName, histTitle, binCountEDepE, binLowEDepE, binHighEDepE, binCountToFP, binLowToFP, binHighToFP, binCountEDepG, binLowEDepG, binHighEDepG);
		*(outSingles3DHist[i]) = *(BGOSingles3DHists[0][i]) + *(BGOSingles3DHists[1][i]);
		outSingles3DHist[i]->SetName(histName);
		outSingles3DHist[i]->SetTitle(histTitle);
	}

	gROOT->cd();

	//TH1D* mc1DHist=(TH1D*) outMult3DHist[0]->Project3D("z");

	// histFilePath=HISTS_DIR;
	// histFilePath+="mcHistM1_BGO_EPG_"+analysisRunString+"_"+gammaHighString+"_BGOEDep.txt";
	// convertTH1ToTXT(mc1DHist,histFilePath);

	//  delete mc1DHist;

	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			delete BGOMult3DHists[j][i];
			delete BGOSingles3DHists[j][i];

		}

	}

	histFile->Write();
	histFile->Close();

}

void makeMCTextHistogramsFrom3D(TString analysisIDString, TString processingRunString, TString mcRIDString, AnalysisCuts theCuts, int multiplicity)
{

	gStyle->SetTitleW(0.95);
	if(multiplicity < 1 || multiplicity > 12)
	{
		cout << "Bad multiplicity entered" << endl;
		return;
	}

	//ExpHist
	TString histFilePath = HISTS_DIR;
	histFilePath += "expAnalysis3DHists_" + processingRunString + ".root";
	TFile expHistFile(histFilePath, "READ");

	TString histNameString = "BGOBGSubMult" + int2str(multiplicity) + "hist";

	TH3D* exp3DHist = (TH3D*) expHistFile.Get(histNameString);
	exp3DHist->SetTitle("Matt Exp. M=1 Data;Electron Energy Desposited (keV);Proton ToF (s);BGO Energy Deposited (keV)");
	exp3DHist->GetXaxis()->SetRangeUser(theCuts.EDepELow, theCuts.EDepEHigh);
	exp3DHist->GetYaxis()->SetRangeUser(theCuts.ToFPLow, theCuts.ToFPHigh);
	exp3DHist->GetZaxis()->SetRangeUser(theCuts.EDepGLow, theCuts.EDepGHigh);

	TString histPath = HISTS_DIR;

	TH1D* liveBGODetectorsHist = (TH1D*) expHistFile.Get("liveBGODetectorsHist");

	double numExpEPEventsTotal = liveBGODetectorsHist->Integral(2, 13);  //M=1-12

	double averageNumDetectors = 0;
	for (int i = 1; i <= 13; i++)
	{
		averageNumDetectors += liveBGODetectorsHist->GetBinContent(i) * (i - 1);

	}
	averageNumDetectors /= numExpEPEventsTotal;

	double numExpEPEventsTotalScaledByDetectors = numExpEPEventsTotal * averageNumDetectors / 12.;

//    double numExpEPGEventsTotal=exp3DHist->Integral();

	//MCHist
	histFilePath = HISTS_DIR;
	histFilePath += "mc3DHists_BGO_EPG_" + mcRIDString + ".root";
	TFile mcHistFile(histFilePath, "READ");

	histNameString = "BGO3DMult" + int2str(multiplicity) + "hist";

	TH3D* mc3DHist = (TH3D*) mcHistFile.Get(histNameString);

	mc3DHist->SetTitle("Geant4 MC M=1 Data;Electron Energy Desposited (keV);Proton ToF (s);BGO Energy Deposited (keV)");
	mc3DHist->Scale(numExpEPEventsTotalScaledByDetectors);
	mc3DHist->GetXaxis()->SetRangeUser(theCuts.EDepELow, theCuts.EDepEHigh);
	mc3DHist->GetYaxis()->SetRangeUser(theCuts.ToFPLow, theCuts.ToFPHigh);
	mc3DHist->GetZaxis()->SetRangeUser(theCuts.EDepGLow, theCuts.EDepGHigh);

	TH1D* mc1DHist = (TH1D*) mc3DHist->Project3D("z");
	TString txtFilePath = HISTS_DIR;
	txtFilePath += "MC_Matt_M" + int2str(multiplicity) + "_" + mcRIDString + "_" + analysisIDString + "EPG_EDepG.txt";
	convertTH1ToTXT(mc1DHist, txtFilePath);
	delete mc1DHist;

	mc1DHist = (TH1D*) mc3DHist->Project3D("x");
	txtFilePath = HISTS_DIR;
	txtFilePath += "MC_Matt_M" + int2str(multiplicity) + "_" + mcRIDString + "_" + analysisIDString + "EPG_EDepE.txt";
	convertTH1ToTXT(mc1DHist, txtFilePath);
	delete mc1DHist;

	mc1DHist = (TH1D*) mc3DHist->Project3D("y");
	txtFilePath = HISTS_DIR;
	txtFilePath += "MC_Matt_M" + int2str(multiplicity) + "_" + mcRIDString + "_" + analysisIDString + "EPG_ToFP.txt";
	convertTH1ToTXT(mc1DHist, txtFilePath);
	delete mc1DHist;
}

void makeBGOSingleEnergyGraphs()
{
	double peakSensitivity = .05;
	const int numEnergies = 45;
	const double energies[numEnergies] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 500, 600, 700, 800 };
	double peakCounts[2][numEnergies];
	double peakCountErrors[2][numEnergies];
	double countsAbove5keV[2][numEnergies];
	double countsAbove5keVErrors[2][numEnergies];
	double countsAbove200keV[2][numEnergies];
	double countsAbove200keVErrors[2][numEnergies];

	for (int runType = 0; runType < 2; runType++)
	{
		for (int i = 0; i < numEnergies; i++)
		{
			peakCounts[runType][i] = 0;
			countsAbove5keV[runType][i] = 0;
			countsAbove200keV[runType][i] = 0;
			TString fileName = MRKGEANTDATA_DIR;
			if(runType == 0)
				fileName += "Results_Geant4_SinglePhoton" + int2str((int) energies[i]) + "keV.root";
			else
				fileName += "Results_Geant4_SinglePhoton" + int2str((int) energies[i]) + "keVWWrap.root";

			TChain* resultsGammaChain = makeResultChain(fileName, "g", 1);

			double BGOEDepG[12];
			double BGOEDepGTotal;

			TString gChainString;
			for (int j = 0; j < 1; j++)
			{
				gChainString = "BGO" + int2str(j + 1) + "EDepTotal";
				resultsGammaChain->SetBranchAddress(gChainString, &(BGOEDepG[j]));

			}

			int numEntries = resultsGammaChain->GetEntries();
			for (int j = 0; j < numEntries; j++)
			{
				resultsGammaChain->GetEntry(j);
				BGOEDepGTotal = 0;
				for (int k = 0; k < 1; k++)
				{
					BGOEDepGTotal += BGOEDepG[k];
				}

				if(abs(energies[i] - BGOEDepGTotal) < peakSensitivity) peakCounts[runType][i] += 1.;
				if(BGOEDepGTotal > 5) countsAbove5keV[runType][i] += 1.;
				if(BGOEDepGTotal > 201) countsAbove200keV[runType][i] += 1.;
			}
			delete resultsGammaChain;
			peakCountErrors[runType][i] = sqrt(peakCounts[runType][i]) / numEntries;
			peakCounts[runType][i] = peakCounts[runType][i] / numEntries;
			countsAbove5keVErrors[runType][i] = sqrt(countsAbove5keV[runType][i]) / numEntries;
			countsAbove5keV[runType][i] = countsAbove5keV[runType][i] / numEntries;

			countsAbove200keVErrors[runType][i] = sqrt(countsAbove200keV[runType][i]) / numEntries;
			countsAbove200keV[runType][i] = countsAbove200keV[runType][i] / numEntries;

		}
	}

	TString imageNameString;
	TCanvas* c2 = getThesisPlotCanvas();

	TGraphErrors peakCountsGraph(numEnergies, energies, peakCounts[0], 0, peakCountErrors[0]);
	TGraphErrors peakCountsWWrapGraph(numEnergies, energies, peakCounts[1], 0, peakCountErrors[1]);
	peakCountsGraph.SetTitle(";Energy Incident (keV);Ratio of Counts in Peak");
	peakCountsGraph.GetYaxis()->SetTitleOffset(1.8);
	peakCountsGraph.Draw("AL");
	peakCountsWWrapGraph.SetMarkerColor(kRed);
	peakCountsWWrapGraph.SetLineColor(kRed);
	peakCountsWWrapGraph.Draw("L SAME");
	imageNameString = GRAPHS_DIR;
	imageNameString += "BGOSingle_PeakCountGraph140619.df";
	c2->SaveAs(imageNameString);

	TGraphErrors countsAbove5keVGraph(numEnergies, energies, countsAbove5keV[0], 0, countsAbove5keVErrors[0]);
	countsAbove5keVGraph.SetTitle(";Energy Incident (keV);Counts Detected over 5 keV per Counts Incident");
	countsAbove5keVGraph.GetYaxis()->SetTitleOffset(1.8);
	countsAbove5keVGraph.Draw("AL");
	imageNameString = GRAPHS_DIR;
	imageNameString += "BGOSingle_CountsAbove5keVGraph140619.png";
	c2->SaveAs(imageNameString);

	TGraphErrors countsAbove200keVGraph(numEnergies, energies, countsAbove200keV[0], 0, countsAbove200keVErrors[0]);
	countsAbove200keVGraph.SetTitle(";Energy Incident (keV);Counts Detected over 5 keV per Counts Incident");
	countsAbove200keVGraph.GetYaxis()->SetTitleOffset(1.8);
	countsAbove200keVGraph.Draw("AL");
	imageNameString = GRAPHS_DIR;
	imageNameString += "BGOSingle_CountsAbove200keVGraph140619.png";
	c2->SaveAs(imageNameString);

	c2->Close();

}

//examine positional light output

void plotQuadraticVsExpResponse()
{
	int detectorNumber = 1;
	const int numEnergies = 19;
	double energies[numEnergies] = { 10, 20, 30, 40, 50, 60, 70, 80, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
	double exponentialLO[numEnergies];
	double quadraticLO[numEnergies];

	for (int i = 0; i < numEnergies; i++)
	{
		quadraticLO[i] = applyEDepModel(energies[i], detectorNumber, nullptr, KEVIN_GAUS_EDEP) / energies[i];
		exponentialLO[i] = applyEDepModel(energies[i], detectorNumber, nullptr, MATT_EXP_EDEP) / energies[i];
	}

	TCanvas theCanvas("responseCanvas", "responseCanvas", 10, 10, 1600, 1400);
	theCanvas.SetLeftMargin(.14);
	theCanvas.SetRightMargin(.05);
	theCanvas.SetLogx();

	TGraph quadraticLOGraph(numEnergies, energies, quadraticLO);
	TGraph exponentialLOGraph(numEnergies, energies, exponentialLO);
	quadraticLOGraph.SetTitle("Energy Response Models;Energy (keV);LO relative to 511 keV");
	quadraticLOGraph.SetName("quadGraph");
	quadraticLOGraph.GetYaxis()->SetTitleOffset(1.8);
	quadraticLOGraph.GetYaxis()->SetRangeUser(.65, 1.15);
	quadraticLOGraph.GetXaxis()->SetRangeUser(8, 1100);

	exponentialLOGraph.SetMarkerColor(kRed);
	exponentialLOGraph.SetName("expGraph");
	exponentialLOGraph.SetLineColor(kRed);
	quadraticLOGraph.Draw("ALP");
	exponentialLOGraph.Draw("LP SAME");
	TLegend theLegend(0.65, 0.35, 0.9, 0.5);
	theLegend.SetFillColor(kWhite);
	theLegend.AddEntry("quadGraph", "Quadratic (Kevin)", "l");
	theLegend.AddEntry("expGraph", "Exponential (MOSZYNSKI)", "l");
	theLegend.Draw();

	TString imageName = GRAPHS_DIR;
	imageName += "EnergyDepositionModels.png";

	theCanvas.SaveAs(imageName);

	theCanvas.Close();

}

double getBGOHighToLowRatio(TH1* theHist, double& ratioError) // presumes 10 - 800 or 1000 keV range with 5 keV bins  Low is 10-200 keV high is 200-800 keV
{
	double highError, lowError;
	double highCounts = theHist->IntegralAndError(39, theHist->GetNbinsX(), highError);
	double lowCounts = theHist->IntegralAndError(1, 38, lowError);
	double ratio = highCounts / lowCounts;
	ratioError = ratio * sqrt(highError * highError / (highCounts * highCounts) + lowError * lowError / (lowCounts * lowCounts));
	cout << "Hist " << theHist->GetName() << "  Ratio: " << ratio << " +/- " << ratioError << endl;
	return ratio;
}

void getBGOHighToLowRatioList(int numHists, string* histFileList)
{
	for (int i = 0; i < numHists; i++)
	{
		gROOT->cd();
		TH1D* theHist = getTabSeperatedHist(histFileList[i], 158, 10, 800);
		double ratioError;
		getBGOHighToLowRatio(theHist, ratioError);
		delete theHist;

	}
}

TH1* combineAndScaleRadSetsToHist(RDK2Set* hmgSet, RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet, TString drawString, TCut epCut, TCut gCut, HistDim inpDimX, HistDim inpDimY)
{
	HistDim inpDimZ = { 0, 0, 0 };
	return combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, drawString, epCut, gCut, inpDimX, inpDimY, inpDimZ);
}

TH1* combineAndScaleRadSetsToHist(RDK2Set* hmgSet, RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet, TString drawString, TCut epCut, TCut gCut, HistDim inpDimX)
{
	HistDim inpDimY = { 0, 0, 0 };
	HistDim inpDimZ = { 0, 0, 0 };
	return combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, drawString, epCut, gCut, inpDimX, inpDimY, inpDimZ);
}

TH1* combineAndScaleRadSetsToHist(RDK2Set* hmgSet, RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet, TString drawString, TCut epCut, TCut gCut, HistDim inpDimX, HistDim inpDimY, HistDim inpDimZ)
{
	double epPer3BodyEvent;
	if(hmgSet == nullptr)
	{
		epPer3BodyEvent = 1;
	}
	else
	{
		epPer3BodyEvent = hmgSet->calcEPRate(epCut);
	}
	TH1* lowHist = nullptr;
	TH1* midHist = nullptr;
	TH1* highHist = nullptr;

	TCut theCut = gCut && epCut;

	if(lowSet != nullptr)
	{
		lowHist = lowSet->makeHist("lowHist", drawString, theCut, inpDimX, inpDimY, inpDimZ);
		lowHist->Scale(1. / (double) lowSet->GetNumFiles());
	}
	if(midSet != nullptr)
	{
		midHist = midSet->makeHist("midHist", drawString, theCut, inpDimX, inpDimY, inpDimZ);
		midHist->Scale(1. / (double) midSet->GetNumFiles());
	}
	if(highSet != nullptr)
	{
		highHist = highSet->makeHist("highHist", drawString, theCut, inpDimX, inpDimY, inpDimZ);
		highHist->Scale(1. / (double) highSet->GetNumFiles());
	}

	TH1* combinedHist = combineAndScaleRadHists("combinedHist", 1, lowHist, midHist, highHist, epPer3BodyEvent);  //NOTE the 1  I scaled before to ensure that different numbers of the radiative types could be used

	if(lowHist != nullptr) delete lowHist;
	if(midHist != nullptr) delete midHist;
	if(highHist != nullptr) delete highHist;

	return combinedHist;

}

TH1* combineAndScaleRadHists(TString newHistName, int numFiles, TH1* radLowHist, TH1* radMidHist, TH1* radHighHist, double epPer3BodyEvent)
{
	TH1* combinedHist = makeEmptyTH1Clone(radMidHist, "combinedHist");

	combinedHist->Sumw2();
	double branchingRatio, scaleFactor;

	if(radLowHist != 0)
	{
		branchingRatio = 0.00429224;
		scaleFactor = branchingRatio / (epPer3BodyEvent * numFiles * 1000000);
		combinedHist->Add(radLowHist, scaleFactor);
	}
	if(radMidHist != 0)
	{
		branchingRatio = 0.00318548;
		scaleFactor = branchingRatio / (epPer3BodyEvent * numFiles * 1000000);
		combinedHist->Add(radMidHist, scaleFactor);
	}
	if(radHighHist != 0)
	{
		branchingRatio = 0.00252823;
		scaleFactor = branchingRatio / (epPer3BodyEvent * numFiles * 1000000);
		combinedHist->Add(radHighHist, scaleFactor);
	}

	return combinedHist;

}

TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString drawString, TCut epCut, TCut gCut, HistDim inpDimX, HistDim inpDimY)
{
	HistDim inpDimZ = { 0, 0, 0 };
	return combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, epCut, gCut, inpDimX, inpDimY, inpDimZ);
}

TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString drawString, TCut epCut, TCut gCut, HistDim inpDimX)
{
	HistDim inpDimY = { 0, 0, 0 };
	HistDim inpDimZ = { 0, 0, 0 };
	return combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, epCut, gCut, inpDimX, inpDimY, inpDimZ);
}

TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString drawString, TCut epCut, TCut gCut, HistDim inpDimX, HistDim inpDimY, HistDim inpDimZ)
{
	TH1* outHist;
	TH1* histPrototype;

	if(inpDimZ.numBins > 0)
		histPrototype = makeTH3DFromDim("radHist", "radHist", inpDimX, inpDimY, inpDimZ);
	else if(inpDimY.numBins > 0)
		histPrototype = makeTH2DFromDim("radHist", "radHist", inpDimX, inpDimY);
	else
		histPrototype = makeTH1DFromDim("radHist", "radHist", inpDimX);
	outHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, epCut, gCut, histPrototype);
	delete histPrototype;
	return outHist;
}

TH1* combineAndScaleRadSetsToHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString drawString, TCut epCut, TCut gCut, TH1* histPrototype)
{
	double epPer3BodyEvent;
	if(threeBodySet == nullptr)
	{
		epPer3BodyEvent = 1;
	}
	else
	{
		epPer3BodyEvent = threeBodySet->calcEPRate(epCut);
	}
	TH1* radHist = nullptr;

	TCut theCut = gCut && epCut;
	radHist = fourBodySet->makeHist("radHist", drawString, theCut, histPrototype);
	if(fourBodySet->GetEventDateString() == "130429")  //Is it Fermi function based?
	{
		radHist->Scale(0.010277 / (epPer3BodyEvent * fourBodySet->GetNumFiles() * 1000000.)); //Now with fermi function
	}
	else
	{
		radHist->Scale((0.00252823 + 0.00318548 + 0.00429224) / (epPer3BodyEvent * fourBodySet->GetNumFiles() * 1000000.));  //Individual branching ratios from previous method add correctly
	}

	return radHist;
}

//From paper on BGO response.  All data is at the moment guessed from graph.
TGraph* getKhodyukBGOResponseGraph()
{
	const int numEntries = 16;
	const double energy[numEntries] = { 10, 13, 13.5, 15, 17, 20, 30, 40, 50, 60, 70, 80, 89, 91, 100, 120 };
	const double lightOutput[numEntries] = { 0.7, 0.74, 0.7, 0.71, 0.69, 0.72, 0.8, 0.845, 0.87, 0.88, 0.89, 0.9, 0.91, 0.875, 0.88, 0.9436 };

	TGraph* outGraph = new TGraph(numEntries, energy, lightOutput);
	outGraph->SetName("bgoResponseGraph");

	return outGraph;

}

//From paper on BGO response.  All data is at the moment guessed from graph.
TGraph* getKhodyukBGOResponseGraph2()
{
	const int numEntries = 15;
	const double energy[numEntries] = { 10, 13, 13.5, 15, 17, 20, 30, 40, 50, 60, 70, 80, 89, 91, 100 };
	const double lightOutput[numEntries] = { 0.7, 0.74, 0.7, 0.71, 0.69, 0.72, 0.8, 0.845, 0.87, 0.88, 0.89, 0.9, 0.91, 0.875, 0.88 };

	TGraph* outGraph = new TGraph(numEntries, energy, lightOutput);
	outGraph->SetName("bgoResponseGraph");

	return outGraph;

}

//From paper on BGO response.
TGraph* getKhodyukBGOResponseGraphOfficial()
{
	const int numEntries = 375;
	double energy[numEntries + 25];
	double lightOutput[numEntries + 25];
	TString filePath = string(HISTS_DIR) + "KhodyukOfficialEnergyResponse.txt";
	if(!FileExists(filePath)) return nullptr;
	ifstream inputFile(filePath);

	TString line;
	getline(inputFile, line);
	for (int i = 0; i < numEntries; i++)
	{
		inputFile >> energy[i] >> lightOutput[i];
		lightOutput[i] *= .01;
	}

//    for (int i=numEntries;i<numEntries+25;i++)
//    {
//        energy[i]=100+(i-numEntries)*900/25.;
//        lightOutput[i]=-.285656766*exp(energy[i]*-.008672992)+1.000916925;
//    }

	TGraph* outGraph = new TGraph(numEntries, energy, lightOutput);
	outGraph->SetName("bgoResponseGraph");

	return outGraph;

}

//Completely made up to suit my purposes
TGraph* getMattBGOResponseGraph()
{
	const int numEntries = 15;
	const double energy[numEntries] = { 10, 13, 13.5, 15, 17, 20, 30, 40, 50, 60, 70, 80, 89, 91, 100 };
	double lightOutput[numEntries] = { 0.7, 0.74, 0.7, 0.71, 0.69, 0.72, 0.8, 0.845, 0.87, 0.88, 0.89, 0.9, 0.91, 0.875, 0.88 };
	for (int i = 0; i < numEntries; i++)
	{
		lightOutput[i] *= .95;
	}

	TGraph* outGraph = new TGraph(numEntries, energy, lightOutput);
	outGraph->SetName("bgoResponseGraph");

	return outGraph;

}

void plotKhodyukBGOResponseGraph()
{

	TGraph* lowEnergyGraph = getMattBGOResponseGraph();
	const int numEntries = 790;
	double energy[numEntries];
	double lightOutput[numEntries];

	//const double A1=-.285656766;
	const double A1 = -0.3904193963;
	const double t1 = -.008672992;
	const double y0 = 1.000916925;

	for (int i = 0; i < numEntries; i++)
	{
		energy[i] = i + 10.;
		if(energy[i] > 100)
			lightOutput[i] = (A1 * exp(energy[i] * t1) + y0);
		else
			lightOutput[i] = lowEnergyGraph->Eval(energy[i]);
	}

	TGraph* fullGraph = new TGraph(numEntries, energy, lightOutput);

	TCanvas theCanvas("plotManyHistCanvas", "plotManyHistCanvas", 10, 10, 1024, 768);
	theCanvas.SetLeftMargin(.15);
	fullGraph->GetYaxis()->SetTitleOffset(1.6);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);
	//gPad->SetLeftMargin(.2);
	fullGraph->SetTitle("Matt's BGO Response Model Based on Khodyuk Data;Incident BGO Energy (keV);Light Output Ratio Relative to 662 keV");
	fullGraph->GetXaxis()->SetRangeUser(0., 800.);
	fullGraph->Draw("AL");

	TString imagePath = GRAPHS_DIR;
	imagePath += "BGOResponseGraphMatt.png";
	theCanvas.SaveAs(imagePath);

	delete fullGraph;
	delete lowEnergyGraph;
}

TH1D* makeBAPDCombinedHist(TString histName, RDK2Set* hmgSet, RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet, HistDim inpDim, TCut inpEPCut, double numEPHits, bool doKEFirst)
{
	TString drawString, gCutString, endString;

	if(doKEFirst)
		endString = "KEFirst";
	else
		endString = "EDepTotal";

	const int numBAPDS = 3;
	TH1D* individualBAPDHist;
	TH1D* bapdMCHist = makeTH1DFromDim(histName, histName, inpDim);
	bapdMCHist->Sumw2();

	for (int i = 0; i < numBAPDS; i++)
	{
		drawString = "gChn.BAPD" + int2str(i + 1) + endString;
		gCutString = "gChn.BAPD" + int2str(i + 1) + endString + " > .5";
		individualBAPDHist = (TH1D*) combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, drawString, inpEPCut, gCutString.Data(), inpDim);
		bapdMCHist->Add(individualBAPDHist);
		delete individualBAPDHist;
	}

	return bapdMCHist;
}

//Presumes histograms are really integers counts (represented as doubles)
TH1D* createHistWBGOResolution(TH1D* inpHist, int detectorNumber, TRandom3* rootRanGen, EDepModel eDepModel)
{

	TString outHistName = string(inpHist->GetName()) + "_WRes";
	TString outHistTitle = string(inpHist->GetTitle()) + " WRes";

	TH1D* outHist = new TH1D(outHistName, outHistTitle, inpHist->GetNbinsX(), getHistExtentLow(inpHist), getHistExtentHigh(inpHist));

	double binCenterEnergy, outputEnergy;
	for (int i = 1; i <= inpHist->GetNbinsX(); i++)
	{
		binCenterEnergy = inpHist->GetBinCenter(i);
		for (int j = 0; j < (int) inpHist->GetBinContent(i); j++)
		{
			outputEnergy = applyEDepModel(binCenterEnergy, detectorNumber, rootRanGen, eDepModel, nullptr);
			if(outputEnergy > getHistExtentLow(outHist) && outputEnergy < getHistExtentHigh(outHist)) outHist->Fill(outputEnergy);

		}
	}

	outHist->Sumw2();
	return outHist;
}

void makeAllMCHistograms(TString dateString, RDK2Set* hmgSet, RDK2Set* lowSet, RDK2Set* midSet, RDK2Set* highSet)
{
	TString HistsDir = HISTS_DIR;
	TH1* tempHist;
	TString baseHistName = "MC_Hist" + dateString;

	tempHist = combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, R_BAPDTOTAL, STD_EPCut, "gDChn.BAPDmultiplicity>0", BAPD_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepG_BAPD_EPG_Sum_" + lowSet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, "eDChn.SBDEDepBlur", STD_EPCut, "gDChn.BAPDmultiplicity>0", EDEPE_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepE_BAPD_EPG_Sum_" + lowSet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, "pChn.SBDTimeFirst", STD_EPCut, "gDChn.BAPDmultiplicity>0", TOFP_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_ToFP_BAPD_EPG_Sum_" + lowSet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, R_BGOTOTALEDGE, STD_EPCut, "gDChn.BGOmultiplicityEdge>0", BGO_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepG_BGO_EPG_Sum_" + highSet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, "eDChn.SBDEDepBlur", STD_EPCut, "gDChn.BGOmultiplicityEdge>0", EDEPE_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepE_BGO_EPG_Sum_" + highSet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(hmgSet, lowSet, midSet, highSet, "pChn.SBDTimeFirst", STD_EPCut, "gDChn.BGOmultiplicityEdge>0", TOFP_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_ToFP_BGO_EPG_Sum_" + highSet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = hmgSet->makeHist("tempHist", "eDChn.SBDEDepBlur", STD_EPHmgCut, EDEPE_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepE_EP_Sum_" + hmgSet->GetEResultsID() + ".txt");
	delete tempHist;

	tempHist = hmgSet->makeHist("tempHist", "pChn.SBDTimeFirst", STD_EPHmgCut, TOFP_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_ToFP_EP_Sum_" + hmgSet->GetEResultsID() + ".txt");
	delete tempHist;
}

void makeAllMCHistograms(TString dateString, RDK2Set* threeBodySet, RDK2Set* fourBodySet)
{
	TString HistsDir = HISTS_DIR;
	TH1* tempHist;
	TString baseHistName = "MC_Hist" + dateString;

	tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, R_BGOTOTALEDGE, STD_EPCut, "gDChn.BGOmultiplicityEdge>0", BGO_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepG_BGO_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, "eDChn.SBDEDepBlur", STD_EPCut, "gDChn.BGOmultiplicityEdge>0", EDEPE_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepE_BGO_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, "pChn.SBDTimeFirst", STD_EPCut, "gDChn.BGOmultiplicityEdge>0", TOFP_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_ToFP_BGO_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = fourBodySet->makeHist("tempHist", "eDChn.SBDEDepBlur", STD_EPCut, EDEPE_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepE_EP_Sum_" + fourBodySet->GetEResultsID() + ".txt");
	delete tempHist;

	tempHist = fourBodySet->makeHist("tempHist", "pChn.SBDTimeFirst", STD_EPCut, TOFP_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_ToFP_EP_Sum_" + fourBodySet->GetEResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, R_BAPDTOTAL, STD_EPCut, "gDChn.BAPDmultiplicity>0", BAPD_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepG_BAPD_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, "eDChn.SBDEDepBlur", STD_EPCut, "gDChn.BAPDmultiplicity>0", EDEPE_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_EDepE_BAPD_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt");
	delete tempHist;

	tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, "pChn.SBDTimeFirst", STD_EPCut, "gDChn.BAPDmultiplicity>0", TOFP_DIM);
	convertTH1ToTXT(tempHist, HistsDir + baseHistName + "_ToFP_BAPD_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt");
	delete tempHist;

}

void makeAllStandardComparisonGraphs(TString mcDateString, RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString expDateString, TString expPIDString, TString expGroupIDString)
{
	TString HistsDir = HISTS_DIR;
	TH1* mcEDepGHist[2] = { getTabSeperatedHist(HistsDir + "MC_Hist" + mcDateString + "_EDepG_BGO_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt", BGO_DIM), getTabSeperatedHist(HistsDir + "MC_Hist130116_EDepG_BGO_EPG_Sum_RID113.txt", BGO_DIM) };
	TH1* mcEDepEHist[2] = { getTabSeperatedHist(HistsDir + "MC_Hist" + mcDateString + "_EDepE_BGO_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt", EDEPE_DIM), getTabSeperatedHist(HistsDir + "MC_Hist130116_EDepE_BGO_EPG_Sum_RID113.txt", EDEPE_DIM) };
	TH1* mcToFPHist[2] = { getTabSeperatedHist(HistsDir + "MC_Hist" + mcDateString + "_ToFP_BGO_EPG_Sum_" + fourBodySet->GetGResultsID() + ".txt", TOFP_DIM), getTabSeperatedHist(HistsDir + "MC_Hist130116_ToFP_BGO_EPG_Sum_RID113.txt", TOFP_DIM) };

	mcEDepGHist[0]->SetTitle("MC_Matt;BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
	mcEDepEHist[0]->SetTitle("MC_Matt;Electron Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
	mcToFPHist[0]->SetTitle("MC_Matt;Proton Time of Flight (s);#frac{EPG}{EP} (count rate ratio)");

	mcEDepGHist[1]->SetTitle("MCOld_Matt;BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
	mcEDepEHist[1]->SetTitle("MCOld_Matt;Electron Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
	mcToFPHist[1]->SetTitle("MCOld_Matt;Proton Time of Flight (s);#frac{EPG}{EP} (count rate ratio)");

	TH1* expEDepGHist[3];
	TH1* expEDepEHist[1];
	TH1* expToFPHist[1];

	expEDepGHist[0] = getTabSeperatedHist(HistsDir + "EXPMatt_Hist" + expDateString + "_" + expPIDString + "_" + expGroupIDString + "_EDepG_BGO_EPG_Sum.txt", BGO_DIM);
	expEDepEHist[0] = getTabSeperatedHist(HistsDir + "EXPMatt_Hist" + expDateString + "_" + expPIDString + "_" + expGroupIDString + "_EDepE_BGO_EPG_Sum.txt", EDEPE_DIM);
	expToFPHist[0] = getTabSeperatedHist(HistsDir + "EXPMatt_Hist" + expDateString + "_" + expPIDString + "_" + expGroupIDString + "_ToFP_BGO_EPG_Sum.txt", TOFP_DIM);

	for (int i = 0; i < 1; i++)
	{
		expEDepGHist[i]->Rebin(5);
		expEDepEHist[i]->Rebin(5);
		expToFPHist[i]->Rebin(5);
	}

	mcEDepGHist[0]->Rebin(5);
	mcEDepEHist[0]->Rebin(5);
	mcToFPHist[0]->Rebin(5);

	mcEDepGHist[1]->Rebin(5);
	mcEDepEHist[1]->Rebin(5);
	mcToFPHist[1]->Rebin(5);

	expEDepGHist[0]->SetTitle("EXP_Matt(Jeff);BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
	expEDepEHist[0]->SetTitle("EXP_Matt(Jeff);Electron Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
	expToFPHist[0]->SetTitle("EXP_Matt(Jeff);Proton Time of Flight (s);#frac{EPG}{EP} (count rate ratio)");

	//Jeff
	HistDim JeffBGODim = { 160, 0, 800. };
	expEDepGHist[1] = getTabSeperatedHist(HistsDir + "Jeff_BGO_EPG_SpectrumInCounts121108.txt", JeffBGODim);
	double jeffScaleFactor = 1. / 5.12044e7;
	expEDepGHist[1]->Scale(jeffScaleFactor);
	expEDepGHist[1]->SetTitle("EXP_Jeff;BGO Energy (keV);#frac{EPG}{EP} (count rate ratio)");

	//Ben - Disabled since I dont' have getBenCSVHist anymore
//    expEDepGHist[2]=getBenCSVHist(HistsDir+"Ben_bgo_energy_12_10_23.csv",13, BGO_DIM);
	expEDepGHist[2]->Rebin(5);
	expEDepGHist[2]->SetTitle("EXP_Ben;BGO Energy (keV);#frac{EPG}{EP} (count rate ratio)");

	TH1* mcOldEDepGHist[1] = { mcEDepGHist[1] };

	plotResidualsToImage(3, expEDepGHist, 1, mcEDepGHist, "BGO Energy Normalized Residuals #frac{Exp-MC}{MC}", "EDepG_Sum_EXP_MC_" + expPIDString + "_Residuals" + mcDateString + ".png");
	plotResidualsToImage(1, expEDepEHist, 1, mcEDepEHist, "Electron Energy Normalized Residuals #frac{Exp-MC}{MC}", "EDepE_Sum_EXP_MC_" + expPIDString + "_Residuals" + mcDateString + ".png");
	plotResidualsToImage(1, expToFPHist, 1, mcToFPHist, "Proton ToF Normalized Residuals #frac{Exp-MC}{MC}", "ToFP_Sum_EXP_MC_" + expPIDString + "_Residuals" + mcDateString + ".png");
	plotResidualsToImage(1, mcEDepGHist, 1, mcOldEDepGHist, "BGO Energy Normalized Residuals #frac{New-Old}{Old}", "EDepG_Sum__MCComparison_Residuals" + mcDateString + ".png");

	stringstream theSStream;
	TString titleString;
	for (int i = 0; i < 3; i++)
	{
		theSStream << expEDepGHist[i]->GetTitle() << ": " << scientific << expEDepGHist[i]->Integral(1, 160);
		titleString = theSStream.str();
		expEDepGHist[i]->SetTitle(titleString);
		theSStream.str("");
	}

	for (int i = 0; i < 2; i++)
	{
		theSStream << mcEDepGHist[i]->GetTitle() << ": " << scientific << mcEDepGHist[i]->Integral(1, 160);
		titleString = theSStream.str();
		mcEDepGHist[i]->SetTitle(titleString);
		theSStream.str("");
	}

	plotExpVersusMCToImage(3, expEDepGHist, 2, mcEDepGHist, "BGO Energy", "EDepG_Sum_EXP_MC_" + expPIDString + "_" + mcDateString + ".png");
	plotExpVersusMCToImage(1, expEDepEHist, 1, mcEDepEHist, "Electron Energy", "EDepE_Sum_EXP_MC_" + expPIDString + "_" + mcDateString + ".png");
	plotExpVersusMCToImage(1, expToFPHist, 1, mcToFPHist, "Proton ToF", "ToFP_Sum_EXP_MC_" + expPIDString + "_" + mcDateString + ".png");
	//old MC vs new mc
	plotExpVersusMCToImage(0, nullptr, 2, mcEDepGHist, "BGO Energy", "EDepG_Sum_MCComparison_" + mcDateString + ".png");

	delete expEDepGHist[0];
	delete expEDepGHist[1];
	delete expEDepGHist[2];
	delete mcEDepGHist[0];
	delete mcEDepGHist[1];
	delete expEDepEHist[0];
	delete mcEDepEHist[0];
	delete expToFPHist[0];
	delete mcToFPHist[0];
}

TGraph* makeKevinResolutionGraph(int detectorNumber)
{
	double x[78];
	double y[78];

	double e_ref = 662;
	double FWHM = 2 * sqrt(2 * log(2));

	for (int i = 0; i < 78; i++)
	{
		double e_dep = 10 + i * 10;
		x[i] = e_dep;

		double sige = 0.1 / FWHM * sqrt(e_ref / e_dep) * e_dep;

		double vare = KevinGausConstx1[detectorNumber - 1] * sige * sige + KevinGausConstx2[detectorNumber - 1] * pow(e_dep / e_ref, 2);
		vare = max(vare, 1e-14);
		sige = sqrt(vare);
		sige = KevinGausCal_Fac[detectorNumber - 1] * sige;
		y[i] = sige;
		//double widthPercentage=sige*FWHM/e_dep*100;
		//y[i]=widthPercentage;
	}
	TString graphName = "ResGraphDet" + int2str(detectorNumber);
	TString graphTitle = "Resolution Graph for Detector " + int2str(detectorNumber);

	TGraph* returnGraph = new TGraph(78, x, y);

	returnGraph->SetName(graphName);
	returnGraph->SetTitle(graphTitle);

	return returnGraph;
}

double evalMattKhodyukParam(double energyIn)
{
	const int polyDegree = 6;
	double coefficients[polyDegree + 1] = { 0, 0, 0, 0, 0, 0, 0 };

	for (int i = 0; i < MattKhodyukParNumRanges; i++)
	{
		if(energyIn >= MattKhodyukParLow[i] && energyIn < MattKhodyukParHigh[i])
		{
			coefficients[0] = MattKhodyukParP0[i];
			coefficients[1] = MattKhodyukParP1[i];
			coefficients[2] = MattKhodyukParP2[i];
			coefficients[3] = MattKhodyukParP3[i];
			coefficients[4] = MattKhodyukParP4[i];
			coefficients[5] = MattKhodyukParP5[i];
			coefficients[6] = MattKhodyukParP6[i];
		}
	}

	double result = coefficients[0];
	for (int i = 0; i < polyDegree; i++)
	{
		result += coefficients[i + 1] * pow(energyIn, i + 1);
	}

	return result;

}

///130709
//double evalMattAllExpParamModel(double energyIn, EDepModel eDepModel)
//{
//    const double energyRegions[5]={.1,13.4,16.5,90,800};
//    double a[4]={0.988335,0.647738,1.02401,1.02005};
//    const double b[4]={4.41739,-1.18558,8.17012,0};
//
//    double result=0;
//    if(eDepModel==ALLEXP_PARAM_MODEL_HIGH)
//    {
//        a[0]*=1.03;
//        a[1]*=1.03;
//        a[2]*=1.03;
//        a[3]=1.01471;
//    }
//    if(eDepModel==ALLEXP_PARAM_MODEL_LOW)
//    {
//        a[0]*=.97;
//        a[1]*=.97;
//        a[2]*=.97;
//        a[3]=1.02578;
//    }
//    for (int i = 0;i< 3;i++)
//    {
//    	if(energyIn >= energyRegions[i] && energyIn <= energyRegions[i+1])
//        {
//            result = a[i]/(1+b[i]/energyIn);
//        }
//    }
//
//    if(energyIn >= energyRegions[3])
//    {
//        result = a[3]/(1+(662*(a[3]-1))/energyIn);
//    }
//
//
//    return result;
//
//}

///130710
double evalMattAllExpParamModel(double energyIn, EDepModel eDepModel)
{
	const double energyRegions[5] = { 0.1, 10., 13.4, 16.5, 90 };
	double a[5] = { 0.879065, 1.25406, 0.647738, 1.02546, 1.02021 };
	double b[5] = { -0.00785819, 9.1513, -1.18558, 8.36477, 0 };
	double c[5] = { -0.00153634, 0, 0, 0, 0 };

	double result = 0;

	if(eDepModel == ALLEXP_PARAM_MODEL_HIGH)
	{
		a[0] *= 1.03;
		a[1] *= 1.03;
		a[2] *= 1.03;
		a[3] *= 1.03;
		a[4] = 1.01471;
	}
	else if(eDepModel == ALLEXP_PARAM_MODEL_LOW)
	{

		a[1] *= .97;
		a[0] = a[1];
		b[0] = b[1];
		a[2] *= .97;
		a[3] *= .97;
		a[4] = 1.02578;
	}
	else
	{
		a[0] = a[1] / (1 + b[1] / energyRegions[1]);
		b[0] = 0;
	}

	for (int i = 0; i < 4; i++)
	{
		if(energyIn >= energyRegions[i] && energyIn <= energyRegions[i + 1])
		{
			result = a[i] / (1 + b[i] / energyIn);
		}
		if(eDepModel == ALLEXP_PARAM_MODEL_HIGH && i == 0)
		{
			result = a[0] + b[0] * energyIn + c[0] * energyIn * energyIn;
		}
	}

	if(energyIn >= energyRegions[4])
	{
		result = a[4] / (1 + (662 * (a[4] - 1)) / energyIn);
	}

	return result;

}

///131209 ends consistent
double evalMattAllExpParamModel2(double energyIn, EDepModel eDepModel)
{
	const double energyRegions[5] = { 0.1, 10., 13.4, 16.5, 90 };
	double a[5] = { 0.879065, 1.25406, 0.647738, 1.02546, 1.02021 };
	double b[5] = { -0.00785819, 9.1513, -1.18558, 8.36477, 0 };

	double result = 0;

	if(eDepModel == ALLEXP_PARAM_MODEL_HIGH)
	{
		a[1] *= 1.03;
		a[0] = a[1];
		b[0] = b[1];
		a[2] *= 1.03;
		a[3] *= 1.03;
		a[4] = 1.01471;
	}
	else if(eDepModel == ALLEXP_PARAM_MODEL_LOW)
	{

		a[1] *= .97;
		a[0] = a[1];
		b[0] = b[1];
		a[2] *= .97;
		a[3] *= .97;
		a[4] = 1.02578;
	}
	else
	{
		a[0] = a[1];
		b[0] = b[1];
	}

	for (int i = 0; i < 4; i++)
	{
		if(energyIn >= energyRegions[i] && energyIn <= energyRegions[i + 1])
		{
			result = a[i] / (1 + b[i] / energyIn);
		}

	}

	if(energyIn >= energyRegions[4])
	{
		result = a[4] / (1 + (662 * (a[4] - 1)) / energyIn);
	}

	return result;

}

///131209 ends consistent
double evalMattAllExpParamModel3(double energyIn, EDepModel eDepModel)
{
	const double energyRegions[5] = { 0.1, 10., 13.4, 16.5, 90 };
	double a[5] = { 0.7951731689, 0.7951731689, 0.647738, 1.02546, 1.0202 };
	double b[5] = { 0.8993068063, 0.8993068063, -1.18558, 8.36477, 0 };

	double result = 0;

	if(eDepModel == ALLEXP_PARAM_MODEL_HIGH)
	{
		a[0] = a[1] = 0.76751833;
		b[0] = b[1] = 0;
		a[2] *= 1.03;
		a[3] *= 1.03;
		a[4] = 1.01471;
	}
	else if(eDepModel == ALLEXP_PARAM_MODEL_LOW)
	{

		a[0] = a[1] = 1.1129584079;
		b[0] = b[1] = 7.2329092132;
		a[2] *= .97;
		a[3] *= .97;
		a[4] = 1.02578;
	}

	for (int i = 0; i < 4; i++)
	{
		if(energyIn >= energyRegions[i] && energyIn <= energyRegions[i + 1])
		{
			result = a[i] / (1 + b[i] / energyIn);
		}

	}

	if(energyIn >= energyRegions[4])
	{
		result = a[4] / (1 + (662 * (a[4] - 1)) / energyIn);
	}

	return result;

}

///Above 16.5 keV consistent
double evalMattAllExpParamModel4(double energyIn, EDepModel eDepModel)
{
	const double energyRegions[5] = { 0.1, 10., 13.4, 16.5, 90 };
	double a[5] = { 0.7951731689, 0.7951731689, 0.647738, 1.02546, 1.0202 };
	double b[5] = { 0.8993068063, 0.8993068063, -1.18558, 8.36477, 0 };

	double result = 0;

	if(eDepModel == ALLEXP_PARAM_MODEL_HIGH)
	{
		a[0] = a[1] = 0.76751833;
		b[0] = b[1] = 0;
		a[2] *= 1.03;

	}
	else if(eDepModel == ALLEXP_PARAM_MODEL_LOW)
	{

		a[0] = a[1] = 1.1129584079;
		b[0] = b[1] = 7.2329092132;
		a[2] *= .97;

	}

	for (int i = 0; i < 4; i++)
	{
		if(energyIn >= energyRegions[i] && energyIn <= energyRegions[i + 1])
		{
			result = a[i] / (1 + b[i] / energyIn);
		}

	}

	if(energyIn >= energyRegions[4])
	{
		result = a[4] / (1 + (662 * (a[4] - 1)) / energyIn);
	}

	return result;

}

///Above 16.5 keV consistent
double evalMattAllExpParamModel5(double energyIn, EDepModel eDepModel)
{
	const double energyRegions[5] = { 0.1, 10., 13.4, 16.5, 90 };
	double a[5] = { 0.7951731689, 0.7951731689, 0.647738, 1.02546, 1.0202 };
	double b[5] = { 0.8993068063, 0.8993068063, -1.18558, 8.36477, 0 };

	double result = 0;

	if(eDepModel == ALLEXP_PARAM_MODEL_HIGH)
	{
		a[0] = a[1] = 0.76751833;
		b[0] = b[1] = 0;

	}
	else if(eDepModel == ALLEXP_PARAM_MODEL_LOW)
	{

		a[0] = a[1] = 1.1129584079;
		b[0] = b[1] = 7.2329092132;

	}

	for (int i = 0; i < 4; i++)
	{
		if(energyIn >= energyRegions[i] && energyIn <= energyRegions[i + 1])
		{
			result = a[i] / (1 + b[i] / energyIn);
		}

	}

	if(energyIn >= energyRegions[4])
	{
		result = a[4] / (1 + (662 * (a[4] - 1)) / energyIn);
	}

	return result;

}

void makeIndividualDetectorHistograms(TString mcAnalysisDateString, RDK2Set* threeBodySet, RDK2Set* fourBodySet)
{

	TH1* mcIndividualHist;

	TString HistsDir = HISTS_DIR;

	TString drawString, cutString;

	TH1* mcGammaAvgHist;
	TH1* mcEDepEAvgHist;
	TH1* mcToFPAvgHist;

	mcGammaAvgHist = makeTH1DFromDim("mcGammaAvgHist", "mcAverage;BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)", BGO_DIM);
	mcGammaAvgHist->Sumw2();

	mcEDepEAvgHist = makeTH1DFromDim("mcEDepEAvgHist", "mcAverage;Electron Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)", EDEPE_DIM);
	mcEDepEAvgHist->Sumw2();

	mcToFPAvgHist = makeTH1DFromDim("mcToFPAvgHist", "mcAverage;Proton Time of Flight (s);#frac{EPG}{EP} (count rate ratio)", TOFP_DIM);
	mcToFPAvgHist->Sumw2();

	TString titleString;

	for (int i = 0; i < 12; i++)
	{
		cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurEdge > 10";

		//Gamma Hist
		drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurEdge";
		mcIndividualHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, STD_EPCut, TCut(cutString), BGO_DIM);
		convertTH1ToTXT(mcIndividualHist, HistsDir + "MC_Hist" + mcAnalysisDateString + "_EDepG_BGO_Det" + int2str(i + 1) + "_EPG_" + fourBodySet->GetGResultsID() + ".txt");
		mcGammaAvgHist->Add(mcIndividualHist);
		delete mcIndividualHist;

		//EDepE Hist
		drawString = "eDChn.SBDEDepBlur";
		mcIndividualHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, STD_EPCut, TCut(cutString), EDEPE_DIM);
		convertTH1ToTXT(mcIndividualHist, HistsDir + "MC_Hist" + mcAnalysisDateString + "_EDepE_BGO_Det" + int2str(i + 1) + "_EPG_" + fourBodySet->GetGResultsID() + ".txt");
		mcEDepEAvgHist->Add(mcIndividualHist);
		delete mcIndividualHist;

		//ToFP Hist
		drawString = "pChn.SBDTimeFirst";
		mcIndividualHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, STD_EPCut, TCut(cutString), TOFP_DIM);
		convertTH1ToTXT(mcIndividualHist, HistsDir + "MC_Hist" + mcAnalysisDateString + "_ToFP_BGO_Det" + int2str(i + 1) + "_EPG_" + fourBodySet->GetGResultsID() + ".txt");
		mcToFPAvgHist->Add(mcIndividualHist);
		delete mcIndividualHist;

	}
	mcGammaAvgHist->Scale(1. / 12.);
	mcEDepEAvgHist->Scale(1. / 12.);
	mcToFPAvgHist->Scale(1. / 12.);
	convertTH1ToTXT(mcGammaAvgHist, HistsDir + "MC_Hist" + mcAnalysisDateString + "_EDepG_BGO_AverageIndiv_EPG_" + fourBodySet->GetGResultsID() + ".txt");
	convertTH1ToTXT(mcEDepEAvgHist, HistsDir + "MC_Hist" + mcAnalysisDateString + "_EDepE_BGO_AverageIndiv_EPG_" + fourBodySet->GetGResultsID() + ".txt");
	convertTH1ToTXT(mcToFPAvgHist, HistsDir + "MC_Hist" + mcAnalysisDateString + "_ToFP_BGO_AverageIndiv_EPG_" + fourBodySet->GetGResultsID() + ".txt");
	delete mcGammaAvgHist;
	delete mcEDepEAvgHist;
	delete mcToFPAvgHist;

	return;
}

void createHistogramsForMCCalibraitonSourceData(int numFiles, TString runString)
{
	TString ResultsDir = RESULTS_DIR;
	TString HistsDir = HISTS_DIR;
	TChain resultsChain("g");

	TString fullPath;
	for (int i = 0; i < numFiles; i++)
	{
		TString fullPath = ResultsDir + "Results_" + runString + "_Set" + int2str(i) + ".root";
		resultsChain.Add(fullPath);
	}

	int numBGOs = 12;

	gROOT->cd();
	TH1D* individualHist;
	TH1D sumHist("SumHist", "SumHist", 10000, 0, 1000);
	sumHist.Sumw2();

	for (int i = 0; i < numBGOs; i++)
	{
		TString histName = "BGO_Det" + int2str(i + 1) + "_" + runString + ".txt";

		individualHist = new TH1D(histName, histName, 10000, 0, 1000);

		TString drawString = "BGO" + int2str(i + 1) + "EDepTotal >> " + histName;
		TString cutString = "BGO" + int2str(i + 1) + "EDepTotal > 0";
		resultsChain.Draw(drawString, cutString);
		convertTH1ToTXT(individualHist, HistsDir + histName);
		individualHist->Sumw2();
		sumHist.Add(individualHist);
		delete individualHist;

	}

	convertTH1ToTXT(&sumHist, HistsDir + "SumBGO_" + runString + ".txt");

}

void createDerivedResultsFile(int startingEvent, TTree* inpResultsTree, TString filePathOut, int numberEvents)
{
	TRandom3* rootRanGen = new TRandom3(0);

	const int numBGOs = 12;
	const int numBAPDs = 3;

	bool hasBGOs=inpResultsTree->GetBranchStatus("BGO1EDepTotal");
	bool hasBAPDs=inpResultsTree->GetBranchStatus("BAPD1EDepTotal");



	double SBDEDepTotal, BGOEDepTotal[numBGOs], BAPDEDepTotal[numBAPDs];

	inpResultsTree->SetBranchAddress("SBDEDepTotal", &SBDEDepTotal);
	TString gChainString;


	double SBDEDepBlur, BGOEDepBlur[numBGOs], BGOEDepBlurEdge[numBGOs], BGOEDepBlurLOAllExpParam[numBGOs], BGOEDepBlurLOAllExpParamLow[numBGOs], BGOEDepBlurLOAllExpParamHigh[numBGOs], BGOEDepResOnly[numBGOs];
//    double BAPDKEFirstNoSiO2[numBAPDS];
	int BGOmultiplicity, BGOmultiplicityEdge, BGOmultiplicityLOAllExpParamHigh, BGOmultiplicityLOAllExpParam, BGOmultiplicityLOAllExpParamLow, BGOmultiplicityResOnly;
	int BAPDmultiplicity;

//    int currentEventOnFluxTree[numBAPD]={0,0,0};

	TString branchName, branchNameWVar;

	TFile dervFile(filePathOut, "recreate");
	dervFile.cd();
	TTree* dervTree = new TTree(DEFAULT_RESULTSDERV_TREENAME, "MRK Results Derived");

	dervTree->Branch("SBDEDepBlur", &SBDEDepBlur, "SBDEDepBlur/D");

	TGraph* bgoResponseGraph=nullptr;

	if(hasBGOs)
	{
		bgoResponseGraph = getKhodyukBGOResponseGraph2();

		for (int i = 0; i < numBGOs; i++)
		{
			gChainString = "BGO" + int2str(i + 1) + "EDepTotal";
			inpResultsTree->SetBranchAddress(gChainString, &(BGOEDepTotal[i]));
		}
		dervTree->Branch("BGOmultiplicity", &BGOmultiplicity, "BGOmultiplicity/I");
		dervTree->Branch("BGOmultiplicityEdge", &BGOmultiplicityEdge, "BGOmultiplicityEdge/I");
		dervTree->Branch("BGOmultiplicityLOAllExpParam", &BGOmultiplicityLOAllExpParam, "BGOmultiplicityLOAllExpParam/I");
		dervTree->Branch("BGOmultiplicityLOAllExpParamLow", &BGOmultiplicityLOAllExpParam, "BGOmultiplicityLOAllExpParamLow/I");
		dervTree->Branch("BGOmultiplicityLOAllExpParamHigh", &BGOmultiplicityLOAllExpParamHigh, "BGOmultiplicityLOAllExpParamHigh/I");
		dervTree->Branch("BGOmultiplicityResOnly", &BGOmultiplicityResOnly, "BGOmultiplicityResOnly/I");
		for (int i = 0; i < numBGOs; i++)
		{
			branchName = "BGO" + int2str(i + 1) + "EDepBlur";
			branchNameWVar = branchName + "/D";
			dervTree->Branch(branchName, &BGOEDepBlur[i], branchNameWVar);

			branchName = "BGO" + int2str(i + 1) + "EDepBlurEdge";
			branchNameWVar = branchName + "/D";
			dervTree->Branch(branchName, &BGOEDepBlurEdge[i], branchNameWVar);

			branchName = "BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam";
			branchNameWVar = branchName + "/D";
			dervTree->Branch(branchName, &BGOEDepBlurLOAllExpParam[i], branchNameWVar);

			branchName = "BGO" + int2str(i + 1) + "EDepBlurLOAllExpParamLow";
			branchNameWVar = branchName + "/D";
			dervTree->Branch(branchName, &BGOEDepBlurLOAllExpParamLow[i], branchNameWVar);

			branchName = "BGO" + int2str(i + 1) + "EDepBlurLOAllExpParamHigh";
			branchNameWVar = branchName + "/D";
			dervTree->Branch(branchName, &BGOEDepBlurLOAllExpParamHigh[i], branchNameWVar);

			branchName = "BGO" + int2str(i + 1) + "EDepResOnly";
			branchNameWVar = branchName + "/D";
			dervTree->Branch(branchName, &BGOEDepResOnly[i], branchNameWVar);

		}
	}

	if(hasBAPDs)
	{
		for (int i = 0; i < numBAPDs; i++)
		{
			gChainString = "BAPD" + int2str(i + 1) + "EDepTotal";
			inpResultsTree->SetBranchAddress(gChainString, &(BAPDEDepTotal[i]));

		}
		dervTree->Branch("BAPDmultiplicity", &BAPDmultiplicity, "BAPDmultiplicity/I");
	}

	for (int i = startingEvent; i < startingEvent + numberEvents; i++)
	{
		inpResultsTree->GetEntry(i);


		if(SBDEDepTotal > 0)
		{
			SBDEDepBlur = SBDEDepTotal + rootRanGen->Gaus(0.0, SBD_ENERGY_SIGMA);
		}
		else
		{
			SBDEDepBlur = 0;
		}

		if(hasBGOs)
		{
			BGOmultiplicity = 0;
			BGOmultiplicityEdge = 0;
			BGOmultiplicityLOAllExpParam = 0;
			BGOmultiplicityLOAllExpParamLow = 0;
			BGOmultiplicityLOAllExpParamHigh = 0;
			BGOmultiplicityResOnly = 0;
			BAPDmultiplicity = 0;

			for (int k = 0; k < numBGOs; k++)
			{
				if(BGOEDepTotal[k] > 0)
				{
					BGOEDepBlur[k] = applyEDepModel(BGOEDepTotal[k], k + 1, rootRanGen, KEVIN_GAUS_EDEP);
					if(BGOEDepBlur[k] < 10)
					{
						// BGOEDepBlur[k]=0;
					}
					else
						BGOmultiplicity++;

					BGOEDepBlurEdge[k] = applyEDepModel(BGOEDepTotal[k], k + 1, rootRanGen, KHODYUK_PARAM_MODEL_WGAUSS, nullptr);
					if(BGOEDepBlurEdge[k] < 10)
					{
						// BGOEDepBlurEdge[k]=0;
					}
					else
						BGOmultiplicityEdge++;

					BGOEDepBlurLOAllExpParam[k] = applyEDepModel(BGOEDepTotal[k], k + 1, rootRanGen, ALLEXP_PARAM_MODEL, nullptr);
					if(BGOEDepBlurLOAllExpParam[k] < 10)
					{
						//BGOEDepBlurLOAllExpParam[k]=0;
					}
					else
						BGOmultiplicityLOAllExpParam++;

					BGOEDepBlurLOAllExpParamLow[k] = applyEDepModel(BGOEDepTotal[k], k + 1, rootRanGen, ALLEXP_PARAM_MODEL_LOW, nullptr);
					if(BGOEDepBlurLOAllExpParamLow[k] < 10)
					{
						// BGOEDepBlurLOAllExpParamLow[k]=0;
					}
					else
						BGOmultiplicityLOAllExpParamLow++;

					BGOEDepBlurLOAllExpParamHigh[k] = applyEDepModel(BGOEDepTotal[k], k + 1, rootRanGen, ALLEXP_PARAM_MODEL_HIGH, nullptr);
					if(BGOEDepBlurLOAllExpParamHigh[k] < 10)
					{
						// BGOEDepBlurLOAllExpParamHigh[k]=0;
					}
					else

						BGOmultiplicityLOAllExpParamHigh++;

					BGOEDepResOnly[k] = applyEDepModel(BGOEDepTotal[k], k + 1, rootRanGen, EDEP_WRESONLY, nullptr);
					if(BGOEDepResOnly[k] < 10)
					{
						//  BGOEDepResOnly[k]=0;
					}
					else
						BGOmultiplicityResOnly++;

				}
				else
				{
					BGOEDepBlur[k] = 0;
					BGOEDepBlurEdge[k] = 0;
					BGOEDepBlurLOAllExpParamHigh[k] = 0;
					BGOEDepBlurLOAllExpParam[k] = 0;
					BGOEDepBlurLOAllExpParamLow[k] = 0;
					BGOEDepResOnly[k] = 0;
				}
			}
		}

		if(hasBAPDs)
		{
			for (int k = 0; k < numBAPDs; k++)
			{
				if(BAPDEDepTotal[k] > 0.3)
				{
					BAPDmultiplicity++;
				}
			}
		}

		dervTree->Fill();

	}

	dervTree->Write("", TObject::kOverwrite);
	dervFile.Close();
	cout << "Made " << filePathOut << endl;

	delete rootRanGen;
	delete bgoResponseGraph;
	gROOT->cd();
	return;
}

void createDerivedResultsFilesForCal(int numFiles, TString runString)
{
	for (int i = 0; i < numFiles; i++)
	{
		TString resultsPath = string(RESULTS_DIR) + "Results_" + runString + "_Set" + int2str(i) + ".root";
		TFile resultsFile(resultsPath, "READ");
		TTree* resultsTree = (TTree*) resultsFile.Get(DEFAULT_RESULTS_TREENAME);

		TString derivedPath = string(RESULTS_DIR) + "Results_" + runString + "_Derv_Set" + int2str(i) + ".root";
		createDerivedResultsFile(0, resultsTree, derivedPath);
		resultsFile.Close();
	}

}

TH1D* createCalMCHistWithRes(int detectorNum, TString mcRunString, TString isotopeString, bool prebakeKhodyuk, double isotopeRatio)
{
	TH1D* mcOutHist;
	TH1D* mcHistNoRes;

	if(isotopeString == "Cs137Ba133")
	{
		TH1D* Cs137Hist = getTabSeperatedHist(TString(HISTS_DIR) + "BGO_Det" + int2str(detectorNum) + "_Cs137_" + mcRunString + ".txt", BGOCAL_DIM);
		TH1D* Ba133Hist = getTabSeperatedHist(TString(HISTS_DIR) + "BGO_Det" + int2str(detectorNum) + "_Ba133_" + mcRunString + ".txt", BGOCAL_DIM);
		Cs137Hist->Scale(isotopeRatio);

		mcHistNoRes = makeTH1DFromDim("tempHist", "tempHist", BGOCAL_DIM);

		*mcHistNoRes = *Cs137Hist + *Ba133Hist;
		delete Cs137Hist;
		delete Ba133Hist;
	}
	else
	{
		mcHistNoRes = getTabSeperatedHist(TString(HISTS_DIR) + "BGO_Det" + int2str(detectorNum) + "_" + isotopeString + "_" + mcRunString + ".txt", BGOCAL_DIM);
	}

	EDepModel theEDepModel;
	if(prebakeKhodyuk)
	{
		theEDepModel = KHODYUK_PARAM_MODEL_WGAUSS;
	}
	else
	{
		theEDepModel = EDEP_WRESONLY;
	}

	TRandom3 ranGen;
	mcOutHist = createHistWBGOResolution(mcHistNoRes, detectorNum, &ranGen, theEDepModel);
	delete mcHistNoRes;
	return mcOutHist;
}

int getMCGaussPeaksToFit(TString seriesName, bool prebakeKhodyuk, TF1** gaussPeakOut)
{
	TString guessMCAmplitudeList;
	TString guessMCMeanList;
	TString guessMCSigmaList;
	TString guessMCRangeList; //normally equivalent to sigma but could be different
	int numPeaksOut;

	//Determine series based information
	if(seriesName == "C504")
	{
		numPeaksOut = 4;
		guessMCAmplitudeList = "400 1300 400 400";
		guessMCMeanList = "33. 79.8 358. 662.";
		guessMCSigmaList = "10 10 30 40";
		guessMCRangeList = "8 15 22 40";

	}
	else if(seriesName == "C461" || seriesName == "C439")
	{
		numPeaksOut = 3;
		guessMCAmplitudeList = "400 1300 400";
		guessMCMeanList = "33. 79.8 358.";
		guessMCSigmaList = "10 10 30";
		guessMCRangeList = "8 15 22";

	}
	else if(seriesName == "C441")
	{
		numPeaksOut = 2;
		guessMCAmplitudeList = "200 400";
		guessMCMeanList = "190 662";
		guessMCSigmaList = "20 40";
		guessMCRangeList = "30 40";

	}
	else
	{
		cout << "Calibration series: " << seriesName << " not recognized." << endl;
		return 0;
	}

//    gaussPeakOut = new TF1*[numPeaksOut];

	for (int i = 0; i < numPeaksOut; i++)
	{
		double guessMCAmplitude, guessMCMean, guessMCSigma, guessMCRange;
		doubleFromList(guessMCAmplitudeList, i, guessMCAmplitude);
		doubleFromList(guessMCMeanList, i, guessMCMean);
		doubleFromList(guessMCSigmaList, i, guessMCSigma);
		doubleFromList(guessMCRangeList, i, guessMCRange);

		if(prebakeKhodyuk)
		{
			TRandom3 ranGen;
			double khodyukGuess = evalMattKhodyukParam(guessMCMean);
			guessMCMean *= khodyukGuess;
			guessMCSigma *= khodyukGuess;
			guessMCRange *= khodyukGuess;
		}

		TString peakName = "MCGauss" + seriesName + "Peak" + int2str(i + 1);
		gaussPeakOut[i] = new TF1(peakName, "gaus", guessMCMean - guessMCRange, guessMCMean + guessMCRange);
		gaussPeakOut[i]->SetParameters(guessMCAmplitude, guessMCMean, guessMCSigma);
	}
	return numPeaksOut;

}

//Returns deconvolution ratio from prebaking khodyuk into the MC (what to correct the mean by)
double deconvolutePrebakeKhodyuk(int numMCFiles, TString mcRunString, int BGONum, TF1* gaussPeakInOut)
{
	cout << "Deconvoluting mean from: " << mcRunString << ", Det " << BGONum << endl;
	double xMinEDepEdgeBlur, xMaxEDepEdgeBlur;

	gaussPeakInOut->GetRange(xMinEDepEdgeBlur, xMaxEDepEdgeBlur);

	TChain resultsChain("g");
	TChain resultsDervChain("rd");

	TString fullPath;
	for (int i = 0; i < numMCFiles; i++)
	{
		fullPath = string(RESULTS_DIR) + "Results_" + mcRunString + "_Set" + int2str(i) + ".root";
		resultsChain.Add(fullPath);
	}

	for (int i = 0; i < numMCFiles; i++)
	{
		fullPath = string(RESULTS_DIR) + "Results_" + mcRunString + "_Derv_Set" + int2str(i) + ".root";
		resultsDervChain.Add(fullPath);
	}
	resultsChain.AddFriend(&resultsDervChain, "rd");

	TString drawString = "BGO" + int2str(BGONum) + "EDepTotal>>tempHist";
	TString cutString = "rd.BGO" + int2str(BGONum) + "EDepBlurEdge > " + d2str(xMinEDepEdgeBlur) + "&& rd.BGO" + int2str(BGONum) + "EDepBlurEdge < " + d2str(xMaxEDepEdgeBlur);

	TH1D mcHist("tempHist", "tempHist", 400, 0, 800);

	resultsChain.Draw(drawString, cutString);

	double mcMeanRaw = mcHist.GetMean();
	double mcMeanWKhodyuk = gaussPeakInOut->GetParameter(1);

	cout << "Deconvolution complete -  MeanWithKhodyuk:  " << mcMeanWKhodyuk << "  MeanDeconvoluted: " << mcMeanRaw << endl;

	return mcMeanRaw / mcMeanWKhodyuk;

}

///Needs "mcHist" to be preloaded
Double_t fitFuncForMCHist(Double_t* v, Double_t* par)
{
	Double_t fitval;
	gROOT->cd();
	TH1D* mcHist;
	gDirectory->GetObject("mcHist", mcHist);
	int bin = mcHist->GetXaxis()->FindBin(v[0]);

	fitval = par[0] * mcHist->GetBinContent(bin);
//    fitval=par[0]*mcHist->Interpolate(v[0]);
	return fitval;
}

TF1* fitMCtoExp(TH1* expHist, TH1* mcHist, int binMin, int binMax)
{
	TString oldMCHistName = mcHist->GetName();
	mcHist->SetName("mcHist");

	double expIntegral = expHist->Integral(binMin, binMax);
	double mcIntegral = mcHist->Integral(binMin, binMax);
	if(expIntegral <= 0 || mcIntegral <= 0)
	{
		return 0;
	}

	TF1* histFitFunc = new TF1("histFitFunc", fitFuncForMCHist, expHist->GetBinLowEdge(binMin), expHist->GetBinLowEdge(binMax) + expHist->GetBinWidth(binMax), 1);
	double low, high;
	histFitFunc->GetRange(low, high);
	cout << low << " " << high << endl;
	histFitFunc->SetParameter(0, expIntegral / mcIntegral);

	cout << "Parameter default set to: " << histFitFunc->GetParameter(0) << endl;
	expHist->Fit("histFitFunc", "NV");
	mcHist->SetName(oldMCHistName);

	return histFitFunc;
}

double getRandom()
{
	return gRandom->Rndm();
}

TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString quantityToDraw, TCut epCut, TCut inpGCut, HistDim inpDimX)
{
	HistDim inpDimY = { 0, 0, 0 };
	HistDim inpDimZ = { 0, 0, 0 };
	return makeMCIndividualDetectorHist(threeBodySet, fourBodySet, quantityToDraw, epCut, inpGCut, inpDimX, inpDimY, inpDimZ);
}

TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString quantityToDraw, TCut epCut, TCut inpGCut, HistDim inpDimX, HistDim inpDimY)
{
	HistDim inpDimZ = { 0, 0, 0 };
	return makeMCIndividualDetectorHist(threeBodySet, fourBodySet, quantityToDraw, epCut, inpGCut, inpDimX, inpDimY, inpDimZ);
}

TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString quantityToDraw, TCut epCut, TCut inpGCut, HistDim inpDimX, HistDim inpDimY, HistDim inpDimZ)
{
	TH1* mcHist;
	TH1* histPrototype;
	if(inpDimZ.numBins == 0)
	{
		if(inpDimY.numBins == 0)
		{
			histPrototype = makeTH1DFromDim("prototype", quantityToDraw, inpDimX);
		}
		else
		{
			histPrototype = makeTH2DFromDim("prototype", quantityToDraw, inpDimX, inpDimY);
		}
	}
	else
	{
		histPrototype = makeTH3DFromDim("prototype", quantityToDraw, inpDimX, inpDimY, inpDimZ);
	}
	mcHist = makeMCIndividualDetectorHist(threeBodySet, fourBodySet, quantityToDraw, epCut, inpGCut, histPrototype);
	delete histPrototype;
	return mcHist;
}

TH1* makeMCIndividualDetectorHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TString quantityToDraw, TCut epCut, TCut inpGCut, TH1* histPrototype)
{
	if(quantityToDraw == "BGOEDepBlurEdge" || quantityToDraw == "SBDTimeFirstP" || quantityToDraw == "SBDEDepBlurE" || quantityToDraw == "BGOEDepBlurLOAllExpParam" || quantityToDraw == "BGOEDepBlurLOAllExpParamLow" || quantityToDraw == "BGOEDepBlurLOAllExpParamHigh" || quantityToDraw == "BGOEDepBlurLOAllExpParam_eg0" || quantityToDraw == "eg0" || quantityToDraw == "BGOEDepTotal" || quantityToDraw == "BAPD")
	{
	}
	else
	{
		cout << "Quantity to draw not recognized." << endl;
		return nullptr;
	}

	TString histName = quantityToDraw + "Hist";
	TH1* mcHist = (TH1*) histPrototype->Clone(histName);
	mcHist->Reset();
	mcHist->Sumw2();

	TString drawString;
	if(quantityToDraw == "BAPD")
	{
		for (int i = 0; i < 3; i++)
		{
			//EPG hists
			TString cutString = "";

			drawString = "gChn.BAPD" + int2str(i + 1) + "EDepTotal";
			cutString = "gChn.BAPD" + int2str(i + 1) + "EDepTotal" + "> .1";

			TCut theGCut = TCut(cutString) && inpGCut;

			TH1* tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, epCut, theGCut, histPrototype);
			mcHist->Add(tempHist);
			delete tempHist;

		}
		mcHist->Scale(1 / 3.);
	}
	else
	{
		for (int i = 0; i < 12; i++)
		{
			//EPG hists
			TString cutString = "";

			if(quantityToDraw == "BGOEDepBlurEdge")
			{
				drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurEdge";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurEdge" + "> 10";
			}
			else if(quantityToDraw == "BGOEDepBlurLOAllExpParam")
			{
				drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam" + "> 10";
			}
			else if(quantityToDraw == "BGOEDepBlurLOAllExpParamLow")
			{
				drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParamLow";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParamLow" + "> 10";
			}
			else if(quantityToDraw == "BGOEDepBlurLOAllExpParamHigh")
			{
				drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParamHigh";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParamHigh" + "> 10";
			}
			else if(quantityToDraw == "SBDTimeFirstP")
			{
				drawString = "pChn.SBDTimeFirst";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurEdge" + "> 10";
			}
			else if(quantityToDraw == "SBDEDepBlurE")
			{
				drawString = "eDChn.SBDEDepBlur";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurEdge" + "> 10";
			}
			else if(quantityToDraw == "BGOEDepBlurLOAllExpParam_eg0")
			{
				drawString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam:eg0";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam" + "> 10";
			}
			else if(quantityToDraw == "eg0")
			{
				drawString = "eg0";
				cutString = "gDChn.BGO" + int2str(i + 1) + "EDepBlurLOAllExpParam" + "> 10";
			}
			else if(quantityToDraw == "BGOEDepTotal")
			{
				drawString = "gChn.BGO" + int2str(i + 1) + "EDepTotal";
				//            cutString ="gChn.BGO"+int2str(i+1)+"EDepTotal" +"> 10";
			}
			TCut theGCut = TCut(cutString) && inpGCut;

			TH1* tempHist = combineAndScaleRadSetsToHist(threeBodySet, fourBodySet, drawString, epCut, theGCut, histPrototype);
			mcHist->Add(tempHist);
			delete tempHist;

		}
		mcHist->Scale(1 / 12.);
	}

	return mcHist;

}

TH1D* createEG0HitPercentHist(RDK2Set* threeBodySet, RDK2Set* fourBodySet, TCut epCut, const HistDim inpDimX)
{
	TH1D* eg0HitHist = (TH1D*) makeMCIndividualDetectorHist(threeBodySet, fourBodySet, "eg0", epCut, "", inpDimX);

	TH1D* eg0Hist = (TH1D*) fourBodySet->makeHist("eg0RawHist", "eg0", "", inpDimX);
	eg0Hist->Scale(0.010212 / (fourBodySet->GetNumFiles() * 1000000.));

	eg0HitHist->Divide(eg0Hist);

	return eg0HitHist;

}

TH1D* expSpecToTheorySpec(TH1D* expHist, RDK2Set* threeBodySet, RDK2Set* fourBodySet, TCut epCut, TString histPath2D, const HistDim inpDimX, const HistDim inpDimY)
{
	TH2D* eDetEG0Hist;
	if(!FileExists(histPath2D))
	{
		eDetEG0Hist = (TH2D*) makeMCIndividualDetectorHist(threeBodySet, fourBodySet, "BGOEDepBlurLOAllExpParam_eg0", epCut, "", inpDimX, inpDimY);

		TFile rootFile(histPath2D, "RECREATE");

		eDetEG0Hist->SetName("eDetEG0Hist");

		eDetEG0Hist->Write();

		rootFile.Close();
	}
	else
	{
		TFile rootFile(histPath2D, "READ");
		eDetEG0Hist = (TH2D*) rootFile.Get("eDetEG0Hist");
	}

	TH1D* percentHist = createEG0HitPercentHist(threeBodySet, fourBodySet, epCut, inpDimX);

	TCanvas* c3 = new TCanvas("c3", "c3", 10, 10, 1280, 1024);
	eDetEG0Hist->Draw("COLZ");

	TString imagePath = string(GRAPHS_DIR) + "EG0_EDet2DGraph.png";
	c3->SaveAs(imagePath);
	delete c3;

	TH1D* eg0ExpHist = makeTH1DFromDim("eg0ExpHist", "eg0ExpHist", inpDimX);
	eg0ExpHist->Sumw2();
	eg0ExpHist->Print();

	for (int i = 0; i < inpDimY.numBins; i++)
	{
		TH1D* eg0ForEDepHist = eDetEG0Hist->ProjectionX("projHist", i + 1, i + 1, "e");
		double integralOfProj = eg0ForEDepHist->Integral();
		if(integralOfProj > 0)
		{
			eg0ForEDepHist->Scale(1. / integralOfProj);

			double expBinContent = expHist->GetBinContent(i + 1);

			eg0ExpHist->Add(eg0ForEDepHist, expBinContent);
		}

		c3 = new TCanvas("c3", "c3", 10, 10, 1280, 1024);
		eg0ForEDepHist->Draw("L");
		imagePath = string(GRAPHS_DIR) + "EG0ForEDet" + d2str(expHist->GetBinCenter(i + 1)) + "_.png";
		c3->SaveAs(imagePath);
		delete c3;

		delete eg0ForEDepHist;
	}

	eg0ExpHist->Divide(percentHist);

	return eg0ExpHist;
}

TH1D* expSpecToTheorySpecWError(TH1D* expHist, RDK2Set* threeBodySet, RDK2Set* fourBodySet, TCut epCut, TString histPath2D)
{
	HistDim extractDim = { 80, 0, 800 };
	TH2D* eDetEG0Hist;
	if(!FileExists(histPath2D))
	{
		eDetEG0Hist = (TH2D*) makeMCIndividualDetectorHist(threeBodySet, fourBodySet, "BGOEDepBlurLOAllExpParam_eg0", epCut, "", extractDim, extractDim);

		TFile rootFile(histPath2D, "RECREATE");

		eDetEG0Hist->SetName("eDetEG0Hist");

		eDetEG0Hist->Write();

		rootFile.Close();
	}
	else
	{
		TFile rootFile(histPath2D, "READ");
		eDetEG0Hist = (TH2D*) rootFile.Get("eDetEG0Hist");
		eDetEG0Hist->SetDirectory(0);
		rootFile.Close();
	}

	TH1D* percentHist = createEG0HitPercentHist(threeBodySet, fourBodySet, epCut, extractDim);

	TCanvas* c3 = new TCanvas("c3", "c3", 10, 10, 1280, 1024);
	eDetEG0Hist->Draw("COLZ");

	TString imagePath = string(GRAPHS_DIR) + "EG0_EDet2DGraph.png";
	c3->SaveAs(imagePath);
	delete c3;

	const int numCurves = 100;
	TH1* eg0ExpHists[numCurves];

	const double scaleToMakeNearCounts = 1;

	for (int k = 0; k < numCurves; k++)
	{
		TString eg0ExpHistName = "eg0ExpHist" + int2str(k);
		eg0ExpHists[k] = makeTH1DFromDim(eg0ExpHistName, eg0ExpHistName, extractDim);

		for (int i = 0; i < extractDim.numBins; i++) //Skip first bin which is 10 keV
		{
			TH1D* eg0ForEDepHist = eDetEG0Hist->ProjectionX("projHist", i + 1, i + 1, "e");
			double integralOfProj = eg0ForEDepHist->Integral();
			if(integralOfProj > 0)
			{
				eg0ForEDepHist->Scale(1. / integralOfProj);

				double expBinContent = expHist->GetBinContent(i + 1) * scaleToMakeNearCounts;
				double expBinError = expHist->GetBinError(i + 1) * scaleToMakeNearCounts;

//                int randomExpBinContent=round(gRandom->Gaus(expBinContent,expBinError));
//
//                for (int j = 0;j< randomExpBinContent;j++)
//                {
//                    eg0ExpHists[k]->Fill(eg0ForEDepHist->GetRandom());
//
//                }

				double randomExpBinContent = gRandom->Gaus(expBinContent, expBinError);
				eg0ExpHists[k]->Add(eg0ForEDepHist, randomExpBinContent);

			}

//            c3 = new TCanvas("c3","c3",10,10,1280,1024);
//            eg0ForEDepHist->Draw("L");
//            imagePath=string(GRAPHS_DIR)+"EG0ForEDet"+d2str(expHist->GetBinCenter(i+1))+"_.png";
//            c3->SaveAs(imagePath);
//            delete c3;

			delete eg0ForEDepHist;
		}
	}

	TString averageHistName = "eg0ExpHistsAverage";
	TH1D* eg0ExpHistsAverage = (TH1D*) calcAvgHistWithErrorFromHists(numCurves, eg0ExpHists, averageHistName);

	eg0ExpHistsAverage->Divide(percentHist);

	eg0ExpHistsAverage->Scale(1. / scaleToMakeNearCounts);

	for (int i = 0; i < numCurves; i++)
	{
		delete eg0ExpHists[i];
	}

	return eg0ExpHistsAverage;
}

void createPlotOfRIDSequence(int numFiles, int beginningRID, TString dateString)
{
	if(dateString == "")
	{
		dateString = getCurrentDateString();
	}
	TString rids[5];
	for (int j = 0; j < 5; j++)
	{
		rids[j] = "RID" + int2str(beginningRID + j);
	}

	RDK2Set threeBodySet("threeBodySet", "STD NonRad W Fermi", numFiles, "130429", "STD_3B", rids[0], rids[1], "");
	RDK2Set fourBodySet("fourBodySet", "STD Rad W Fermi", numFiles, "130429", "STD_4B_100eV_Endpoint", rids[2], rids[3], rids[4]);

//    threeBodySet.makeDerivedResultsFiles(PROTON);
//    threeBodySet.makeDerivedResultsFiles(ELECTRON);
	double epgRate = fourBodySet.calcEPGRateFor12IndividualBGOs(STD_EPCut);
	double epgBAPDRate = fourBodySet.calcEPGRateFor3IndividualBAPDs(STD_EPCut);
//    fourBodySet.makeDerivedResultsFiles(PROTON);
//    fourBodySet.makeDerivedResultsFiles(ELECTRON);
//    fourBodySet.makeDerivedResultsFiles(GAMMA);

	TString histPath, titleString, imageString;
	TH1* mcHists[1];

	//PTOF
	mcHists[0] = threeBodySet.makeHist("tofHist", "pChn.SBDTimeFirst", STD_EPCut, TOFP_DIM);
	histPath = HISTS_DIR;
	histPath += "MC_Hist" + dateString + "_ToFP_EP_" + rids[4] + ".txt";
	convertTH1ToTXT(mcHists[0], histPath);

	double epRate = mcHists[0]->Integral();

	titleString = rids[4] + ";E-P Timing (us); Counts (AU)";
	mcHists[0]->SetTitle(titleString);
	imageString = "MC_ToFP_EP_" + rids[4] + "_" + dateString + ".png";
	plotExpVersusMCToImage(0, nullptr, 1, mcHists, "Proton Time of Flight " + rids[4] + " EP Rate: " + d2str(epRate), imageString);
	delete mcHists[0];

	//Electron Energy
	mcHists[0] = threeBodySet.makeHist("eEnergy", "eDChn.SBDEDepBlur", STD_EPCut, EDEPE_DIM);
	histPath = HISTS_DIR;
	histPath += "MC_Hist" + dateString + "_EDepE_EP_" + rids[4] + ".txt";
	convertTH1ToTXT(mcHists[0], histPath);

	titleString = rids[4] + ";SBD Electron Energy Detected (keV); Counts (AU)";
	mcHists[0]->SetTitle(titleString);
	imageString = "MC_EDepE_EP_" + rids[4] + "_" + dateString + ".png";
	plotExpVersusMCToImage(0, nullptr, 1, mcHists, "Electron Energy " + rids[4] + " EP Rate: " + d2str(epRate), imageString);
	delete mcHists[0];

	//BGO
	mcHists[0] = makeMCIndividualDetectorHist(&threeBodySet, &fourBodySet, "BGOEDepBlurLOAllExpParam", STD_EPCut, "", BGO_DIM);
	histPath = HISTS_DIR;
	histPath += "MC_Hist" + dateString + "_EDepGWLO_AverageIndiv_EPG_" + rids[4] + ".txt";
	convertTH1ToTXT(mcHists[0], histPath);

	titleString = rids[4] + ";BGO Energy Detected (keV); #frac{EPG}{EP}";
	mcHists[0]->SetTitle(titleString);
	imageString = "MC_EPG_" + rids[4] + "_" + dateString + ".png";
	plotExpVersusMCToImage(0, nullptr, 1, mcHists, "#frac{EPG}{EP} " + rids[4] + " EP Rate: " + d2str(epgRate), imageString);
	delete mcHists[0];

	//BAPD
	mcHists[0] = makeMCIndividualDetectorHist(&threeBodySet, &fourBodySet, "BAPD", STD_EPCut, "", BAPD_DIM);
	histPath = HISTS_DIR;
	histPath += "MC_Hist" + dateString + "_BAPD_AverageIndiv_EPG_" + rids[4] + ".txt";
	convertTH1ToTXT(mcHists[0], histPath);

	titleString = rids[4] + ";BAPD Energy Detected (keV); #frac{EPG}{EP}";
	mcHists[0]->SetTitle(titleString);
	imageString = "MC_EPGBAPD_" + rids[4] + "_" + dateString + ".png";
	plotExpVersusMCToImage(0, nullptr, 1, mcHists, "#frac{EPG}{EP} " + rids[4] + " EP Rate: " + d2str(epgBAPDRate), imageString);

	delete mcHists[0];
}

double BAPDCollectionEfficiencyModel(double inpEnergy, double depth)
{
	int DetNum = 1;
	double collectionEfficiency = 0.;

	if(depth > 60.001 || depth < 0)
	{

		typeAnythingToContinue("CE Error depth not in range.");
	}

	if(depth > 3.5355)
	{
		collectionEfficiency = 1.;
		if(DetNum == 2)
		{
			//Apply g factor as in paper
			if(inpEnergy < 2500)
				collectionEfficiency *= 1.12;
			else if(inpEnergy < 5900) collectionEfficiency *= -.000035294118 * inpEnergy + 1.208235; //Slope from 2.5 kev to 5.9 keV from 1.12 to 1.

		}

	}
	else
	{
		if(DetNum == 0)
		{
			if(depth > .6894)
				collectionEfficiency = 0.0495 * depth + 0.8259;
			else
				collectionEfficiency = 1.2162 * depth + .0215;
		}
		else if(DetNum == 1)
		{

			if(depth > .7248)
				collectionEfficiency = 0.0758 * depth + 0.73306;
			else if(depth > .1591)
				collectionEfficiency = 1.2162 * depth - 0.0935;
			else
				collectionEfficiency = .5657 * depth + 0.01;
		}
		else if(DetNum == 2)
		{
			if(depth > 3.5355)
				collectionEfficiency = 1.;
			else if(depth > 0.4066)
				collectionEfficiency = 0.0339 * depth + 0.8792;
			else
				collectionEfficiency = 1.2626 * depth + 0.37947;

			if(DetNum == 2)
			{
				//Apply g factor as in paper
				if(inpEnergy < 2500)
					collectionEfficiency *= 1.12;
				else if(inpEnergy < 5900) collectionEfficiency *= -.000035294118 * inpEnergy + 1.208235; //Slope from 2.5 kev to 5.9 keV from 1.12 to 1.

			}

		}
	}

	return collectionEfficiency;
}

void makeCombinedBGOAndBAPDPlot(TString expBGOPath, double bgoCorrection, double epBGOCount, TString expBAPDPath, double bapdCorrection, double epBAPDCount, int mcAIDNum, bool logx, bool logy, bool absScale, bool varWidth, bool usePDF)
{
	const TString HistsDir = HISTS_DIR;
	const TString GraphsDir = GRAPHS_DIR;

	double dataToResidRatio = 3;
	TString absScaleString = "";
	if(absScale)
	{
		absScaleString = "_scaled_";
	}
	TString varWidthString = "";
	if(varWidth) varWidthString = "_var_";
	TString logxString = "";
	if(logx) logxString = "_logx_";
	TString logyString = "";
	if(logy) logyString = "_logy_";

	TString imageString = "ExpPlot_AID" + int2str(mcAIDNum) + absScaleString + varWidthString + logxString + logyString;
	if(usePDF)
	{
		imageString += ".pdf";
	}
	else
	{
		imageString += ".root";
	}
	gStyle->SetErrorX(0);

	vector<double> var_bins = { 10, 12, 13, 14, 15, 17, 19, 21, 24, 27, 30, 34, 38, 43, 49, 56, 63, 71, 80, 92, 200, 780 };
	TH1D* expHistInit = getTabSeperatedHist(expBGOPath, BGO_DIM2);
	TH1D* mcHistInit = getTabSeperatedHist(HistsDir + "Hist_MC_AID" + int2str(mcAIDNum) + "_EPG_EDepGAvg.txt", BGO_DIM2);
	TH1D* theoryHist = getTabSeperatedHist(HistsDir + "EG0_BGODIM_Hist.txt", BGO_DIM2);

	cout << "BGO Counts: " << expHistInit->Integral() << endl;
	TH1* expHist = makeVarHistScaledByWidth(expHistInit, var_bins);
	TH1* mcHist = makeVarHistScaledByWidth(mcHistInit, var_bins);
	mcHist->SetName("theMCHist");

	double error;
	double integral = mcHistInit->IntegralAndError(1, mcHistInit->GetNbinsX(), error);
	cout << "MC: " << integral << " +/- " << error << endl;

	integral = expHistInit->IntegralAndError(1, expHistInit->GetNbinsX(), error);
	cout << "Exp Uncorrected: " << integral << " +/- " << error << endl;

	//Exp Scaled   Disabled corrections
//	expHist->Scale(1.+bgoCorrection);
//	expHistInit->Scale(1.+bgoCorrection);
//
//	integral=expHistInit->IntegralAndError(1,expHistInit->GetNbinsX(),error);
//	cout << "Exp Corrected: " << integral << " +/- " << error << endl;

	if(!absScale)
	{
		scaleHistogramsTogether(expHist, mcHist);
		scaleHistogramsTogether(expHistInit, mcHistInit);
	}
	else
	{
		mcHistInit->Scale(epBGOCount);
		mcHist->Scale(epBGOCount);
	}

	vector<double> var_bin_log = { 10, 12, 13, 14, 15, 17, 19, 21, 24, 27, 30, 34, 38, 43, 49, 56, 63, 71, 80, 92, 150, 200, 250, 300, 400, 500, 600, 700, 780 };
	TH1* mcHistCompare = makeVarHistScaledByWidth(mcHistInit, var_bin_log);
	TH1* theoryHistCompare = makeVarHistScaledByWidth(theoryHist, var_bin_log);
	scaleHistogramsTogether(mcHistCompare, theoryHistCompare);
	mcHistCompare->GetXaxis()->SetRangeUser(10, 1000);
	theoryHistCompare->GetXaxis()->SetRangeUser(10, 1000);

	vector<double> var_bins_bapd = { 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3, 1.5, 1.7, 2, 2.4, 3, 3.6, 4.3, 5.1, 6.1, 7.8, 14 };
	TH1D* expHistInitBAPD = getTabSeperatedHist(expBAPDPath, BAPD_DIM3);
	cout << "BAPD Counts: " << expHistInitBAPD->Integral() << endl;
	TH1D* mcHistInitBAPD = getTabSeperatedHist(HistsDir + "Hist_MC_AID" + int2str(mcAIDNum) + "_EPBG_EDepBGAvg.txt", BAPD_DIM3);
	TH1D* theoryHistBAPD = getTabSeperatedHist(HistsDir + "EG0_BAPDDIM_Hist.txt", BAPD_DIM3);

	TH1* expHistBAPD = makeVarHistScaledByWidth(expHistInitBAPD, var_bins_bapd);
	TH1* mcHistBAPD = makeVarHistScaledByWidth(mcHistInitBAPD, var_bins_bapd);
	mcHistBAPD->SetName("theMCHist");

	integral = mcHistInitBAPD->IntegralAndError(1, 140, error);
	cout << "MC Init: " << integral << " +/- " << error << endl;

//	//Exp Scaled
//	expHistBAPD->Scale(1.+bapdCorrection);
//	expHistInitBAPD->Scale(1.+bapdCorrection);

	integral = expHistInitBAPD->IntegralAndError(1, 140, error);
	cout << "Exp Uncorrected: " << integral << " +/- " << error << endl;

	if(!absScale)
	{
		scaleHistogramsTogether(expHistBAPD, mcHistBAPD, 1, 140);
		scaleHistogramsTogether(expHistInitBAPD, mcHistInitBAPD, 1, 140);
	}
	else
	{
		mcHistInitBAPD->Scale(epBAPDCount);
		mcHistBAPD->Scale(epBAPDCount);
	}

	scaleHistogramsTogether(mcHistInitBAPD, theoryHistBAPD, 1, 140);

	vector<double> var_bin_log_bad = { 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3, 1.5, 1.7, 2, 2.4, 3, 3.6, 4.3, 5.1, 6.1, 7.8, 14 };
	TH1* mcHistCompareBAPD = makeVarHistScaledByWidth(mcHistInitBAPD, var_bin_log_bad);
	TH1* theoryHistCompareBAPD = makeVarHistScaledByWidth(theoryHistBAPD, var_bin_log_bad);
	scaleHistogramsTogether(mcHistCompareBAPD, theoryHistCompareBAPD);

	TCanvas theCanvas("c2", "c2", 10, 10, 2400, 1200);
	gPad->SetFillColor(kWhite);

	TLatex Tla;
	Tla.SetTextFont(42);
	Tla.SetTextSize(.04);
	Tla.DrawLatexNDC(.03, .94, "(a)");

	double chi2 = expHist->Chi2Test(mcHist, "WWCHI2");

	TPaveText chi2Text(.20, .91, 0.45, .97);
	chi2Text.SetTextColor(1);
	chi2Text.SetTextFont(42);
	chi2Text.SetTextSize(.034);
	chi2Text.SetBorderSize(0);
	chi2Text.SetFillStyle(0);
	int ndf;

	ndf = expHist->GetNbinsX() - 1;
	TString latexString = "#chi^{2}/ndf = " + d2str(chi2, 1, 1) + "/" + int2str(ndf) + ", p = " + d2str(TMath::Prob(chi2, ndf), 2, 1);
	chi2Text.AddText(latexString);
	chi2Text.Draw("same");

	TLatex Tl;
	Tl.SetTextFont(42);
	Tl.SetTextSize(.04);
	Tl.DrawLatexNDC(.55, .94, "(b)");

	chi2 = expHistBAPD->Chi2Test(mcHistBAPD, "WWCHI2");

	TPaveText chi2TextBAPD(.70, .91, 0.95, .97);
	chi2TextBAPD.SetTextColor(1);
	chi2TextBAPD.SetTextFont(42);
	chi2TextBAPD.SetTextSize(.034);
	chi2TextBAPD.SetBorderSize(0);
	chi2TextBAPD.SetFillStyle(0);

	ndf = expHistBAPD->GetNbinsX() - 1;
	latexString = "#chi^{2}/ndf = " + d2str(chi2, 1, 1) + "/" + int2str(ndf) + ", p = " + d2str(TMath::Prob(chi2, ndf), 2, 1);
	chi2TextBAPD.AddText(latexString);
	chi2TextBAPD.Draw("same");

	gPad->SetLeftMargin(.0);
	gPad->SetRightMargin(.0);
	gPad->SetBottomMargin(0);

	if(logx && logy)
	{
		theCanvas.Divide(2, 0, 0, 0);
	}
	else
	{
		theCanvas.Divide(0, 2, 0, 0);
	}

	theCanvas.cd(1);
	gPad->SetLeftMargin(.0);
	gPad->SetRightMargin(.0);
	gPad->SetBottomMargin(0);
	gPad->SetTopMargin(0.05);
	TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
	TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
	pad1.SetFillColor(kWhite);
	pad2.SetFillColor(kWhite);
	pad1.SetLeftMargin(.15);
	pad2.SetLeftMargin(.15);
	pad1.SetRightMargin(.05);
	pad2.SetRightMargin(.05);
	pad1.SetBottomMargin(0);
	pad1.SetTopMargin(0);
	pad2.SetTopMargin(0);
	pad2.SetBottomMargin(.12);
	pad1.Draw();
	pad2.Draw();

	mcHistCompare->SetLineColor(kRed);
	theoryHistCompare->SetLineColor(kBlue);
	theoryHistCompare->SetLineStyle(2);
	theoryHistCompare->SetLineWidth(2);
	mcHistCompare->SetLineWidth(2);
	mcHistCompare->GetXaxis()->SetRangeUser(10, 600);
	theoryHistCompare->GetXaxis()->SetRangeUser(10, 600);
//
	expHist->SetMarkerStyle(8);
	expHist->SetMarkerSize(1.5);
	expHist->SetLineWidth(1);

	pad2.cd();
	if(logx) pad2.SetLogx(1);
	if(logy) pad2.SetLogy(1);
	pad2.SetGrid();

	if(logx) pad1.SetLogx(1);
	TH1F* hr = pad2.DrawFrame(9, 1, 1010, 900);

	hr->GetXaxis()->SetTitleSize(0.05);
	hr->GetYaxis()->SetTitleSize(0.05);

	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetXaxis()->SetMoreLogLabels(1);
//	hr->GetYaxis()->SetMoreLogLabels(1);
	hr->GetYaxis()->SetNoExponent(1);
	hr->GetXaxis()->SetNoExponent(1);
	hr->GetYaxis()->SetLabelSize(0.05);
	hr->GetXaxis()->SetTitleOffset(1.2);
	hr->GetYaxis()->SetTitleOffset(1.4);
	hr->GetYaxis()->CenterTitle();
	hr->GetYaxis()->SetNdivisions(3, 0, 0, false);
	hr->GetXaxis()->SetNdivisions(3, 0, 0, false);
	hr->SetTitle(";BGO photon pulse height (511=511 keV);BGO ep#gamma counts per keV");
	pad2.cd();
//	mcHist->Draw("HIST 9 same");

	theoryHistCompare->Draw("HIST C 9 same");
	mcHistCompare->Draw("HIST C 9 same");
	expHist->Draw("E 9 same");

	convertTH1ToTXT(theoryHistCompare, TString(HISTS_DIR) + "Hist_Paper_BGO_Theory.txt", true);
	convertTH1ToTXT(mcHistCompare, TString(HISTS_DIR) + "Hist_Paper_BGO_MC.txt", true);
	convertTH1ToTXT(expHist, TString(HISTS_DIR) + "Hist_Paper_BGO_Exp.txt", true);

	cout << "Exp bin 1: " << expHist->GetBinContent(3) << endl;
	cout << "theory bin 1: " << theoryHistCompare->GetBinCenter(1) << " " << mcHistCompare->GetBinContent(1) << endl;

	TLegend theLegend(0.19, .2, 0.66, 0.4);
	theLegend.SetTextSize(0.05);
	theLegend.SetBorderSize(1);
	theLegend.SetFillColor(kWhite);
	theLegend.AddEntry(theoryHistCompare, "Theory (scaled)", "L");
	theLegend.AddEntry(mcHistCompare, "Simulation (scaled)", "L");
	theLegend.AddEntry(expHist, "Experiment", "EP");

	theLegend.Draw("same");

	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
	residualHist->Sumw2();
	residualHist->Add(mcHist, -1);
	residualHist->Divide(mcHist);

	residualHist->SetMarkerColor(kBlack);
	residualHist->SetLineColor(kBlack);
	residualHist->GetXaxis()->SetLabelSize(0.2);
	residualHist->GetYaxis()->SetLabelSize(0.2);
	residualHist->GetXaxis()->SetTitleSize(0.2);
	residualHist->GetYaxis()->SetTitleSize(0.2);
	//residualHist->GetXaxis()->SetNdivisions(0,0,0);

	pad1.cd();
	pad1.SetLogy(0);
	pad1.SetGrid();
	TH1F* hr1 = pad1.DrawFrame(9, -0.25, 1010, 0.45);

	hr1->GetXaxis()->SetTitleSize(0.05);
	hr1->GetYaxis()->SetTitleSize(0.05 * dataToResidRatio);
	hr1->SetTitleOffset(.45, "y");
	TLine lineToDraw(9, 0, 1010, 0);
	lineToDraw.SetLineColor(kBlack);
	lineToDraw.SetLineWidth(2);
	lineToDraw.Draw("same");

	hr1->GetXaxis()->SetLabelSize(0.05);
	hr1->GetYaxis()->SetLabelSize(0.05 * dataToResidRatio);
	hr1->GetYaxis()->CenterTitle();
	hr1->SetTitle(";; #frac{Exp.- Sim.}{ Sim.}");
	hr1->SetTitle("");
	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
	hr1->GetXaxis()->SetNdivisions(3, 0, 0, false);
	residualHist->Draw("E 9 same");

	TLegend theLegend2(0.42, .06, 0.94, 0.32);
	theLegend2.SetTextSize(0.152);
	theLegend2.SetBorderSize(1);
	theLegend2.SetFillColor(kWhite);
	theLegend2.AddEntry(residualHist, "Normalized residuals", "ep");
	theLegend2.Draw("same");

//BAPDS
	theCanvas.cd(2);
	gPad->SetLeftMargin(.0);
	gPad->SetRightMargin(.0);
	gPad->SetBottomMargin(0);
	gPad->SetTopMargin(0.05);
	TPad pad1BAPD("pad1BAPD", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
	TPad pad2BAPD("pad2BAPD", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
	pad1BAPD.SetFillColor(kWhite);
	pad2BAPD.SetFillColor(kWhite);
	pad1BAPD.SetLeftMargin(.15);
	pad2BAPD.SetLeftMargin(.15);
	pad1BAPD.SetRightMargin(.05);
	pad2BAPD.SetRightMargin(.05);
	pad1BAPD.SetBottomMargin(0);
	pad1BAPD.SetTopMargin(0);
	pad2BAPD.SetTopMargin(0);
	pad2BAPD.SetBottomMargin(.12);
	pad1BAPD.Draw();
	pad2BAPD.Draw();

	mcHistCompareBAPD->SetLineColor(kGreen + 2);
	theoryHistCompareBAPD->SetLineColor(kBlue);
	theoryHistCompareBAPD->SetLineStyle(2);
	theoryHistCompareBAPD->SetLineWidth(2);
	mcHistCompareBAPD->SetLineWidth(2);
//
	expHistBAPD->SetMarkerStyle(8);
	expHistBAPD->SetMarkerSize(1.5);
	expHistBAPD->SetLineWidth(1);

	pad2BAPD.cd();
	if(logx) pad2BAPD.SetLogx(1);
	if(logy) pad2BAPD.SetLogy(1);
	pad2BAPD.SetGrid();

	if(logx) pad1BAPD.SetLogx(1);
	TH1F* hrBAPD = pad2BAPD.DrawFrame(0.3, 5, 14, 700);

	hrBAPD->GetXaxis()->SetTitleSize(0.05);
	hrBAPD->GetYaxis()->SetTitleSize(0.05);

	hrBAPD->GetXaxis()->SetLabelSize(0.05);
//	hrBAPD->GetXaxis()->SetMoreLogLabels(1);
//	hrBAPD->GetYaxis()->SetMoreLogLabels(1);
	hrBAPD->GetYaxis()->SetNoExponent(1);
	hrBAPD->GetXaxis()->SetNoExponent(1);
	hrBAPD->GetYaxis()->SetLabelSize(0.05);
	hrBAPD->GetXaxis()->SetTitleOffset(1.2);
	hrBAPD->GetYaxis()->SetTitleOffset(1.4);
	hrBAPD->GetXaxis()->SetNdivisions(3, 0, 0, false);
	hrBAPD->GetYaxis()->SetNdivisions(3, 0, 0, false);
	hrBAPD->GetYaxis()->CenterTitle();
	hrBAPD->SetTitle(";APD photon pulse height (5.9=5.9 keV);APD ep#gamma counts per keV");
	pad2BAPD.cd();

	theoryHistCompareBAPD->Draw("HIST C 9 same");
	mcHistCompareBAPD->Draw("HIST C 9 same");
	expHistBAPD->Draw("E 9 same");

	convertTH1ToTXT(theoryHistCompareBAPD, TString(HISTS_DIR) + "Hist_Paper_BAPD_Theory.txt", true);
	convertTH1ToTXT(mcHistCompareBAPD, TString(HISTS_DIR) + "Hist_Paper_BAPD_MC.txt", true);
	convertTH1ToTXT(expHistBAPD, TString(HISTS_DIR) + "Hist_Paper_BAPD_Exp.txt", true);

	TLegend theLegendBAPD(0.19, .2, 0.66, 0.4);
	theLegendBAPD.SetTextSize(0.05);
	theLegendBAPD.SetBorderSize(1);
	theLegendBAPD.SetFillColor(kWhite);
	theLegendBAPD.AddEntry(theoryHistCompareBAPD, "Theory (scaled)", "L");
	theLegendBAPD.AddEntry(mcHistCompareBAPD, "Simulation (scaled)", "L");
	theLegendBAPD.AddEntry(expHistBAPD, "Experiment", "EP");

	theLegendBAPD.Draw("same");

	TH1D* residualHistBAPD = (TH1D*) expHistBAPD->Clone("residualHistBAPD");
	residualHistBAPD->Sumw2();
	residualHistBAPD->Add(mcHistBAPD, -1);
	residualHistBAPD->Divide(mcHistBAPD);

	residualHistBAPD->SetMarkerColor(kBlack);
	residualHistBAPD->SetLineColor(kBlack);
	residualHistBAPD->GetXaxis()->SetLabelSize(0.2);
	residualHistBAPD->GetYaxis()->SetLabelSize(0.2);
	residualHistBAPD->GetXaxis()->SetTitleSize(0.2);
	residualHistBAPD->GetYaxis()->SetTitleSize(0.2);
	//residualHist->GetXaxis()->SetNdivisions(0,0,0);

	pad1BAPD.cd();
	pad1BAPD.SetLogy(0);
	pad1BAPD.SetGrid();
	TH1F* hr1BAPD = pad1BAPD.DrawFrame(0.3, -1.1, 15, 1.1);

	hr1BAPD->GetXaxis()->SetTitleSize(0.05);
	hr1BAPD->GetYaxis()->SetTitleSize(0.05 * dataToResidRatio);
	hr1BAPD->SetTitleOffset(.45, "y");
	TLine lineToDrawBAPD(0.3, 0, 15, 0);
	lineToDrawBAPD.SetLineColor(kBlack);
	lineToDrawBAPD.SetLineWidth(2);
	lineToDrawBAPD.Draw("same");

	hr1BAPD->GetXaxis()->SetLabelSize(0.05);
	hr1BAPD->GetYaxis()->SetLabelSize(0.05 * dataToResidRatio);
	hr1BAPD->GetYaxis()->CenterTitle();
	hr1BAPD->SetTitle(";; #frac{Exp.- Sim.}{ Sim.}");
	hr1BAPD->SetTitle("");
	hr1BAPD->GetYaxis()->SetNdivisions(5, 0, 0);
	hr1BAPD->GetXaxis()->SetNdivisions(3, 0, 0, false);
	residualHistBAPD->Draw("E 9 same");

	TLegend theLegend2BAPD(0.42, .06, 0.94, 0.32);
	theLegend2BAPD.SetTextSize(0.152);
	theLegend2BAPD.SetBorderSize(1);
	theLegend2BAPD.SetFillColor(kWhite);
	theLegend2BAPD.AddEntry(residualHistBAPD, "Normalized residuals", "ep");
	theLegend2BAPD.Draw("same");

	theCanvas.SaveAs(GraphsDir + imageString);

}

void makeEGTimingPlot(TString histPath, TString bhistPath)
{
	gStyle->SetPadTickY(0);
	TCanvas theCanvas("c2", "c2", 10, 10, 1200, 800);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(.10);
	gPad->SetRightMargin(.12);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.05);

	double sf = 40e-3;

	TH1D* avgHist = getTabSeperatedHist(histPath, 1800, -26, 46);
	TH1D* bavgHist = getTabSeperatedHist(bhistPath, 1800, -26, 46);

	cout << "Num:" << avgHist->Integral() << endl;

	avgHist->SetTitle(";Difference in arrival of electron & photon (#lower[0.05]{#scale[1.2]{#mu}}s);BGO Counts");

	avgHist->SetLineColor(kRed);
	avgHist->SetMarkerColor(kRed);
	bavgHist->SetLineColor(kGreen + 2);
	bavgHist->SetMarkerColor(kGreen + 2);
	avgHist->SetMarkerStyle(8);
	bavgHist->SetMarkerStyle(8);
	bavgHist->SetLineStyle(2);
	avgHist->SetLineWidth(2);
	bavgHist->SetLineWidth(2);
//    avgHist->Rebin(5);
//    bavgHist->Rebin(5);

	avgHist->GetXaxis()->SetRangeUser(-325 * sf, 575 * sf);
	avgHist->GetYaxis()->SetRangeUser(0, 1000);

	avgHist->GetXaxis()->SetTitleSize(0.05);
	avgHist->GetYaxis()->SetTitleSize(0.05);

	avgHist->GetXaxis()->SetLabelSize(0.05);
//	avgHist->GetXaxis()->SetMoreLogLabels(1);
//	avgHist->GetYaxis()->SetMoreLogLabels(1);
	avgHist->GetYaxis()->SetNoExponent(1);
	avgHist->GetXaxis()->SetNoExponent(1);
	avgHist->GetYaxis()->SetLabelSize(0.05);
	avgHist->GetXaxis()->SetTitleOffset(1.2);
	avgHist->GetYaxis()->SetTitleOffset(1.03);
	avgHist->GetYaxis()->CenterTitle();
	avgHist->GetYaxis()->SetNdivisions(509, true);

//	avgHist->Draw("P E1");
//	bavgHist->Draw("P E1 same");
//	avgHist->Draw("P E1 same");

	avgHist->Draw("HIST");
	theCanvas.Update();

	double textPosX = 0.15;
	double textWidth = 0.15;
	TPaveText textBox1(textPosX, 0.73, textPosX + textWidth, 0.91, "brNDC");
	textBox1.SetTextSize(0.045);
	textBox1.AddText("BGO");
	textBox1.AddText("Pre-peak");
	textBox1.AddText("Window");
	textBox1.AddText("(Bkgd.)");
	textBox1.SetTextFont(42);
	textBox1.SetFillColor(0);
	textBox1.SetFillStyle(0);
	textBox1.SetLineColor(0);
	textBox1.SetShadowColor(0);
	textBox1.SetLineWidth(0);
	textBox1.SetTextAlign(23);
	textBox1.Draw("same");

	TPaveText textBox2(0.5, 0.73, 0.82, 0.91, "brNDC");
	textBox2.SetTextSize(0.045);
	textBox2.AddText("BGO");
	textBox2.AddText("Post-peak");
	textBox2.AddText("Window");
	textBox2.AddText("(Bkgd.)");
	textBox2.SetTextFont(42);
	textBox2.SetFillColor(0);
	textBox2.SetFillStyle(0);
	textBox2.SetLineColor(0);
	textBox2.SetShadowColor(0);
	textBox2.SetLineWidth(0);
	textBox2.SetTextAlign(23);
	textBox2.Draw("same");

	textPosX = 0.456;
	textWidth = 0.11;
	TPaveText textBox3(textPosX, 0.73, textPosX + textWidth, 0.91, "brNDC");
	textBox3.SetTextSize(0.045);
	textBox3.AddText("BGO");
	textBox3.AddText("On-peak");
	textBox3.AddText("Window");
	textBox3.AddText("(Signal)");
	textBox3.SetTextFont(42);
	textBox3.SetFillColor(0);
	textBox3.SetFillStyle(0);
	textBox3.SetLineColor(0);
	textBox3.SetShadowColor(0);
	textBox3.SetLineWidth(0);
	textBox3.SetTextAlign(23);
	textBox3.Draw("same");

	TArrow ar2(3.3, 650, 5.5, 700, 0.03, "<|");
	ar2.SetAngle(27);
	ar2.SetLineWidth(4);
	ar2.Draw();

	TLine* line[2];
//	double xLine[4]={-325*sf,-75*sf,75*sf,575*sf};
	double xLine[2] = { -75 * sf, 75 * sf };
	for (int i = 0; i < 2; i++)
	{
		line[i] = new TLine(xLine[i], 0, xLine[i], 1000);
		line[i]->SetLineStyle(2);
		line[i]->Draw("same");
	}
	avgHist->Draw("HIST same");

	Float_t rightmax = 1.5 * bavgHist->GetMaximum();
	Float_t scale = gPad->GetUymax() / rightmax;
	bavgHist->Scale(scale);
	bavgHist->Draw("HIST same");
	TGaxis axis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax, 510, "+L");
	axis.SetTitle("APD Counts");
	axis.SetNoExponent(1);
	axis.SetNoExponent(1);
	axis.SetLabelSize(0.05);
	axis.SetTitleSize(0.05);
	axis.SetTitleOffset(1.2);
	axis.SetTitleOffset(1.3);
	axis.SetLabelFont(42);
	axis.SetTitleFont(42);
	axis.CenterTitle();
	axis.Draw();
	axis.SetTitleOffset(1.1);

	TLegend theLegend(0.67, .35, 0.82, 0.55);
	theLegend.SetTextSize(0.05);
	theLegend.SetBorderSize(1);
	theLegend.SetFillColor(kWhite);
	theLegend.AddEntry(avgHist, "BGO", "L");
	theLegend.AddEntry(bavgHist, "APD", "L");
	theLegend.Draw("same");

	theCanvas.SaveAs(TString(GRAPHS_DIR) + "ExpPlot_EP_EGTime.eps");
}

void makeZ0Hist()
{
	HistDim z0Dim = { 280, -.4, .72 };
	TString HistsDir = HISTS_DIR;

//	int numFiles=1;
//	RDK2CutSet STD_EPCutSetH("STD_EPCutSetH", "Standard Cuts 151015", 100., 757., 13., 31., 2.e-6, 25.e-6, 10, 780., 0.4, 30.);  // Final Cuts?
//	RDK2MCAnalysis theAnalysis( "AID86",  "Final???? (Geant4 version update and 25 nm)", numFiles,  "RID528", "RID529",  numFiles, "RID530",  "RID531", "RID532",  STD_EPCutSetH);
//	TH1* z0EPHist=theAnalysis.MakeStandardHist("epHist",DET_EP,"z0",z0Dim);
//	TH1* z0EPGHist=theAnalysis.MakeStandardHist("epgHist",DET_EPG,"z0",z0Dim);
//	TH1* z0EPBGHist=theAnalysis.MakeStandardHist("epbgHist",DET_EPBG,"z0",z0Dim);
//
//	convertTH1ToTXT(z0EPHist,HistsDir+"z0EPHist160201.txt");
//	convertTH1ToTXT(z0EPGHist,HistsDir+"z0EPGHist160201.txt");
//	convertTH1ToTXT(z0EPBGHist,HistsDir+"z0EPBGHist160201.txt");

	TH1* z0EPHist = getTabSeperatedHist(HistsDir + "z0EPHist151215.txt", z0Dim);
	TH1* z0EPGHist = getTabSeperatedHist(HistsDir + "z0EPGHist151215.txt", z0Dim);
	TH1* z0EPBGHist = getTabSeperatedHist(HistsDir + "z0EPBGHist151215.txt", z0Dim);

	TCanvas* theCanvas = new TCanvas("thesisPlot", "thesisPlot", 10, 10, 2400, 1000);
	theCanvas->SetLeftMargin(.085);
	theCanvas->SetRightMargin(.01);
	theCanvas->SetTopMargin(.05);
	theCanvas->SetBottomMargin(.18);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);

	for (int i = 1; i <= z0Dim.numBins; i++)
	{
		z0EPHist->SetBinError(i, 0);
		z0EPGHist->SetBinError(i, 0);
		z0EPBGHist->SetBinError(i, 0);

	}

	z0EPHist->Rebin(2);
	z0EPGHist->Rebin(2);
	z0EPBGHist->Rebin(2);

	scaleHistTo(z0EPHist, 800);
	scaleHistTo(z0EPGHist, 300);
	scaleHistTo(z0EPBGHist, 150);

	z0EPHist->SetTitle(";;");
//    z0EPHist->SetTitle(";;");

	z0EPHist->SetLineColor(kBlue);
	z0EPGHist->SetLineColor(kRed);
	z0EPBGHist->SetLineColor(kGreen + 2);

	z0EPHist->SetLineWidth(2);
	z0EPGHist->SetLineWidth(3);
	z0EPBGHist->SetLineWidth(3);

	z0EPGHist->SetLineStyle(2);
	z0EPBGHist->SetLineStyle(3);

	z0EPHist->GetXaxis()->SetRangeUser(-.12, .325);
	z0EPGHist->GetXaxis()->SetRangeUser(-.12, .325);
	z0EPBGHist->GetXaxis()->SetRangeUser(-.12, .325);
//
//	z0EPHist->GetYaxis()->SetRangeUser(0,.020*1000);

	TH1F* hr = theCanvas->DrawFrame(-.39, 0, .45, .020 * 1000);
	hr->GetYaxis()->SetTitleOffset(.6);
	hr->GetYaxis()->SetLabelSize(.08);
	hr->GetXaxis()->SetLabelSize(.08);
	hr->GetYaxis()->SetNdivisions(508);

	z0EPGHist->Draw("L same");
	z0EPBGHist->Draw("L same");
	z0EPHist->Draw("L same");

	TLegend legend(0.11, 0.4, 0.42, 0.9, nullptr, "brNDC");
	legend.SetTextSize(0.08);
//    legend.SetTextSize(1.4);
	legend.SetFillColor(kWhite);
	legend.AddEntry(z0EPHist, "#it{ep}", "L");
	legend.AddEntry(z0EPGHist, "#it{ep#gamma} (BGO)", "L");
	legend.AddEntry(z0EPBGHist, "#it{ep#gamma} (APD)", "L");
	legend.Draw();

	theCanvas->SaveAs(TString(GRAPHS_DIR) + "z0Plot160201.eps");

	delete theCanvas;
}

void makeOtherParameterPlot(TString expID, double correction, int mcAIDNum, TString coTypeStr, TString plotType, TString imageType, bool halfPlotNoResid)
{
	TString expPlotType;
	HistDim theDim;
	int rebinNum = 1;
	double xMin = 0, xMax = 1000;
	TString typeTitleString;
	if(plotType == "EDepE")
	{
		expPlotType = "eEn";
		typeTitleString = "Electron Energy Deposited (keV)";
		theDim = EDEPE_DIM2;
		rebinNum = 200;
		xMin = 0;
		xMax = 800;
		if(coTypeStr == "EPBG") rebinNum *= 4;
	}
	else if(plotType == "EDepP")
	{
		expPlotType = "pEn";
		typeTitleString = "Proton Energy Deposited (keV)";
		theDim = EDEPE_DIM2;
		rebinNum = 2;
		xMin = 10;
		xMax = 35;
		if(coTypeStr == "EPG") rebinNum *= 2;
		if(coTypeStr == "EPBG") rebinNum *= 4;
	}
	else if(plotType == "ToFP")
	{
		expPlotType = "epTOF";
		typeTitleString = "Proton Time of Flight (s)";
		theDim = TOFP_DIM5;
		rebinNum = 5;
		xMin = 0;
		xMax = 30e-6;
		if(coTypeStr == "EPBG") rebinNum = 25;
	}

	double dataToResidRatio = 3;
	TString imageString = "AID" + int2str(mcAIDNum) + "_" + plotType + "_" + coTypeStr + "_scaled";
	if(halfPlotNoResid)
	{
		imageString += "_half";
	}
	imageString += "." + imageType;
	gStyle->SetErrorX(0);

	TCanvas* theCanvas;
	TPad* pad1 = nullptr; //resididual hist
	TPad* pad2 = nullptr;
	if(halfPlotNoResid)
	{
		theCanvas = new TCanvas("largeCanvas", "largeCanvas", 10, 10, 600, 700);
		gPad->SetFillColor(kWhite);
		gPad->SetLeftMargin(.15);
		gPad->SetRightMargin(.05);
		gPad->SetBottomMargin(.12);
		gPad->SetTopMargin(0.05);
	}
	else
	{
		theCanvas = new TCanvas("largeCanvas", "largeCanvas", 10, 10, 1200, 1200);
		gPad->SetFillColor(kWhite);
		pad1 = new TPad("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
		pad2 = new TPad("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
		pad1->SetFillColor(kWhite);
		pad1->SetLeftMargin(.15);
		pad1->SetRightMargin(.05);
		pad1->SetBottomMargin(0);
		pad1->SetTopMargin(0);
		pad1->Draw();
		pad2->SetLeftMargin(.15);
		pad2->SetFillColor(kWhite);
		pad2->SetRightMargin(.05);
		pad2->SetTopMargin(0);
		pad2->SetBottomMargin(.12);
		pad2->Draw();

	}

	TString expDetStr = "BGO";
	if(coTypeStr == "EPBG") expDetStr = "bAPD";

	TString expCoTypeStr = "ep";
	if(coTypeStr == "EPG" || coTypeStr == "EPBG") expCoTypeStr = "epg";

	TH1D* expHist = getTabSeperatedHist(TString(HISTS_DIR) + expDetStr + "_" + expID + "_" + expCoTypeStr + "_" + expPlotType + ".txt", theDim);
	TH1D* mcHist = getTabSeperatedHist(TString(HISTS_DIR) + "Hist_MC_AID" + int2str(mcAIDNum) + "_" + coTypeStr + "_" + plotType + ".txt", theDim);

	if(coTypeStr == "EPG")
		mcHist->SetLineColor(kRed);
	else if(coTypeStr == "EPBG")
		mcHist->SetLineColor(kGreen + 2);
	else
		mcHist->SetLineColor(kBlue);

	mcHist->SetLineStyle(2);
	mcHist->SetLineWidth(4);

	expHist->SetMarkerStyle(20);
	expHist->SetMarkerSize(1);
	expHist->SetLineWidth(4);

	mcHist->Print();
	double error;
	double integral = mcHist->IntegralAndError(1, mcHist->GetNbinsX(), error);
	cout << "MC: " << integral << " +/- " << error << endl;
	expHist->Print();

	mcHist->Scale(1e6);
	expHist->Scale((1 + correction) * 1e6);
	mcHist->Rebin(rebinNum);
	expHist->Rebin(rebinNum);
	scaleHistogramsTogether(expHist, mcHist);

	TH1F* hr;
	if(halfPlotNoResid)
	{
		hr = theCanvas->DrawFrame(xMin, 0, xMax, expHist->GetMaximum() * 1.15);
		gPad->SetGrid();
	}
	else
	{
		hr = pad2->DrawFrame(xMin, 0, xMax, expHist->GetMaximum() * 1.15);
		pad2->cd();
		pad2->SetGrid();
	}

	hr->GetXaxis()->SetTitleSize(0.05);
	hr->GetYaxis()->SetTitleSize(0.05);
	hr->GetYaxis()->CenterTitle();
	hr->GetXaxis()->SetLabelSize(0.05);
	hr->GetYaxis()->SetLabelSize(0.05);
	hr->GetYaxis()->SetTitleOffset(1.6);
	TString titleString = ";" + typeTitleString;
	if(coTypeStr == "EP")
		titleString += ";ep counts";
	else if(coTypeStr == "EPG")
		titleString += ";BGO ep#gamma counts";
	else if(coTypeStr == "EPBG") titleString += ";APD ep#gamma counts";

	hr->SetTitle(titleString);
	mcHist->Draw("HIST 9 same");
	expHist->Draw("E1 9 same");

	TLegend theLegend(0.60, .65, 0.90, 0.9);
	theLegend.SetBorderSize(1);
	theLegend.SetFillColor(kWhite);
	theLegend.AddEntry(expHist, "Experiment", "ep");
	theLegend.AddEntry(mcHist, "Simulation", "L");
	theLegend.Draw("same");

	TH1D* residualHist = nullptr;
	if(!halfPlotNoResid)
	{
		residualHist = (TH1D*) expHist->Clone("residualHist");
		residualHist->Sumw2();
		residualHist->Add(mcHist, -1);
		residualHist->Divide(mcHist);

		residualHist->SetMarkerColor(kBlack);
		residualHist->SetLineColor(kBlack);
		residualHist->GetXaxis()->SetLabelSize(0.2);
		residualHist->GetYaxis()->SetLabelSize(0.2);
		residualHist->GetXaxis()->SetTitleSize(0.2);
		residualHist->GetYaxis()->SetTitleSize(0.2);

		pad1->cd();
		pad1->SetLogy(0);
		pad1->SetGrid();
		TH1F* hr1 = pad1->DrawFrame(xMin, -1.4, xMax, 1.4);

		hr1->GetXaxis()->SetTitleSize(0.05);
		hr1->GetYaxis()->SetTitleSize(0.05 * dataToResidRatio);
		hr1->SetTitleOffset(.32, "y");
		hr1->GetYaxis()->CenterTitle();
		hr1->GetXaxis()->SetLabelSize(0.05);
		hr1->GetYaxis()->SetLabelSize(0.05 * dataToResidRatio);

		hr1->SetTitle(";;#frac{Exp.-Sim.}{Sim.}");
		hr1->GetYaxis()->SetNdivisions(5, 0, 0);

		residualHist->Draw("E1 same");

		TLine lineToDraw(0, 0, 20, 0);
		lineToDraw.SetLineColor(kBlack);
		lineToDraw.Draw("same");

		TLegend theLegend2(0.50, .1, 0.90, 0.3);
		theLegend2.SetBorderSize(1);
		theLegend2.SetFillColor(kWhite);
		theLegend2.AddEntry(residualHist, "Normalized Residuals", "ep");
		theLegend2.Draw("same");
	}
	theCanvas->SaveAs(TString(GRAPHS_DIR) + imageString);

	delete expHist;
	delete mcHist;
	if(residualHist != nullptr)
	{
		delete residualHist;
	}
	delete theCanvas;
}

void makeOtherParameterPlots(TString bgoExpID, double bgoCorrection, TString apdExpID, double apdCorrection, int mcID, bool halfPlotNoResids)
{
	makeOtherParameterPlot(bgoExpID, 1, mcID, "EP", "EDepE", "eps", halfPlotNoResids);
	makeOtherParameterPlot(bgoExpID, bgoCorrection, mcID, "EPG", "EDepE", "eps", halfPlotNoResids);
	makeOtherParameterPlot(apdExpID, apdCorrection, mcID, "EPBG", "EDepE", "eps", halfPlotNoResids);

	makeOtherParameterPlot(bgoExpID, 1, mcID, "EP", "EDepP", "eps", halfPlotNoResids);
	makeOtherParameterPlot(bgoExpID, bgoCorrection, mcID, "EPG", "EDepP", "eps", halfPlotNoResids);
	makeOtherParameterPlot(apdExpID, apdCorrection, mcID, "EPBG", "EDepP", "eps", halfPlotNoResids);

	makeOtherParameterPlot(bgoExpID, 1, mcID, "EP", "ToFP", "eps", halfPlotNoResids);
	makeOtherParameterPlot(bgoExpID, bgoCorrection, mcID, "EPG", "ToFP", "eps", halfPlotNoResids);
	makeOtherParameterPlot(apdExpID, apdCorrection, mcID, "EPBG", "ToFP", "eps", halfPlotNoResids);
}

void outputParameterizedLOModelToText()
{
	ofstream outFile;

	outFile.open((TString(HISTS_DIR)+"ParamLOModel3_160530.txt"));

	outFile << "#energy(keV)\tlower yield curve\tbest estimate yield curve\tupper yield curve" << endl;
	double energy, lower,best,upper;
	for (int i = 0;i< BGOCAL_DIM.numBins;i++)
	{
		energy=i/10. + 0.05;
		lower=evalMattAllExpParamModel3(energy,ALLEXP_PARAM_MODEL_LOW);
		best=evalMattAllExpParamModel3(energy,ALLEXP_PARAM_MODEL);
		upper=evalMattAllExpParamModel3(energy,ALLEXP_PARAM_MODEL_HIGH);

		outFile << energy << "\t" << lower << "\t" << best << "\t" << upper;
		if(i != BGOCAL_DIM.numBins -1)
		{
			outFile <<endl;
		}

	}
	outFile.close();
}

void outputBGODataToFile()
{
	ofstream theFile(TString(HISTS_DIR)+"bgoExperimentsData.txt");
	theFile << "#experiment/energy deposit/incident/LO/LO uncertainty" << endl;

	const int numDewarPoints=14;
	double dewarEnergies[numDewarPoints]={122.1,81,59.5,31,99.8,88,23.3,122.1,99.8,88,23.3,122.1,185.7,107};
	double dewarLOs[numDewarPoints]={0.90033,0.91603,0.88772,0.81225,0.86602,0.90747,0.73221,0.8827,0.86363,0.90637,0.7364,0.89272,0.93333,0.8203};
	double dewarLOErrors[numDewarPoints]={0.010849,0.010122,0.010078,0.010542,0.012021,0.010652,0.03,0.010572,0.010722,0.010078,0.025,0.010022,0.015,0.02};
	for(int i=0; i <numDewarPoints; i++)
	{
		theFile << "0\t" << dewarEnergies[i] << "\t" << dewarLOs[i] << "\t" << dewarLOErrors[i] << endl;
	}

	const int numCalRDK2Points=6;
	double calRDK2Energies[numCalRDK2Points]={32.3,42.7,78.7,79.7,121.9,356};
	double calRDK2LOs[numCalRDK2Points]={0.81,0.86,0.94,0.93,0.92,0.986};
	double calRDK2LOErrors[numCalRDK2Points]={0.02,0.02,0.01,0.01,0.01,0.002};
	for(int i=0; i <numCalRDK2Points; i++)
	{
		theFile << "1\t" << calRDK2Energies[i] << "\t" << calRDK2LOs[i] << "\t" << calRDK2LOErrors[i] << endl;
	}


	//TGraph* khodyukGraph=getKhodyukBGOResponseGraphOfficial();
	//Assumed error was similar to moszynski and was 2%
	const int numKhodyukPoints=375;
	double khodyukEnergies[numKhodyukPoints]={11,11,12,12,13,13.025,13.05,13.075,13.1,13.125,13.15,13.175,13.2,13.225,13.25,13.275,13.3,13.325,13.35,13.375,
		13.4,13.425,13.45,13.475,13.5,13.525,13.55,13.575,13.6,13.625,13.65,13.675,13.7,13.725,13.75,13.775,13.8,13.825,13.85,13.875,13.9,13.925,13.95,13.975,14,
		14.025,14.05,14.075,14.1,14.125,14.15,14.175,14.2,14.225,14.25,14.275,14.3,14.325,14.35,14.375,14.4,14.425,14.45,14.475,14.5,14.525,14.55,14.575,14.6,14.625,
		14.65,14.675,14.7,14.725,14.75,14.775,14.8,14.825,14.85,14.875,14.9,14.925,14.95,14.975,15,15.025,15.05,15.075,15.1,15.125,15.15,15.175,15.2,15.225,15.25,
		15.275,15.3,15.325,15.35,15.375,15.4,15.425,15.45,15.475,15.5,15.525,15.55,15.575,15.6,15.625,15.65,15.675,15.7,15.725,15.75,15.775,15.8,15.825,15.85,15.875,
		15.9,15.925,15.95,15.975,16,16.025,16.05,16.075,16.1,16.125,16.15,16.175,16.2,16.225,16.25,16.275,16.3,16.325,16.35,16.375,16.4,16.425,16.45,16.475,16.5,16.525,
		16.55,16.575,16.6,16.625,16.65,16.675,16.7,16.725,16.75,16.775,16.8,16.825,16.85,16.875,16.9,16.925,16.95,16.975,17,17.025,17.05,17.075,17.1,17.125,17.15,17.175,
		17.2,17.225,17.25,17.275,17.3,17.325,17.35,17.375,17.4,17.425,17.45,17.475,17.5,17.525,17.55,17.575,17.6,17.625,17.65,17.675,17.7,17.725,17.75,17.775,17.8,17.825,
		17.85,17.875,17.9,17.925,17.95,17.975,18,18.025,18.05,18.075,18.1,18.125,18.15,18.175,18.2,18.225,18.25,18.275,18.3,18.325,18.35,18.375,18.4,18.425,18.45,18.475,
		18.5,18.525,18.55,18.575,18.6,18.625,18.65,18.675,18.7,18.725,18.75,18.775,18.8,18.825,18.85,18.875,18.9,18.925,18.95,18.975,19,19.025,19.05,19.075,19.1,19.125,
		19.15,19.175,19.2,19.225,19.25,19.275,19.3,19.325,19.35,19.375,19.4,19.425,19.45,19.475,19.5,19.525,19.55,19.575,19.6,19.625,19.65,19.675,19.7,19.725,19.75,19.775,
		19.8,19.825,19.85,19.875,19.9,19.925,19.95,19.975,20,20,25,25,30,30,30,30,35,35,40,40,45,45,45,45,50,50,55,55,60,60,65,65,70,70,75,75,80,80,85,85,89.5,89.6,89.7,
		89.8,89.9,90,90.1,90.2,90.3,90.4,90.5,90.6,90.7,90.8,90.9,91,91.1,91.2,91.3,91.4,91.5,91.6,91.7,91.8,91.9,92,92.1,92.2,92.3,92.4,92.5,92.6,92.7,92.8,92.9,93,93.1,
		93.2,93.3,93.4,93.5,93.6,93.7,93.8,93.9,94,94.1,94.2,94.3,94.4,94.5,94.6,94.7,94.8,94.9,95,95,100,100};
	double khodyukLOs[numKhodyukPoints]={0.700,0.707,0.725,0.724,0.740,0.742,0.737,0.743,0.743,0.741,0.737,0.740,0.738,0.736,0.743,0.746,0.743,0.739,0.742,0.742,
		0.736,0.722,0.703,0.704,0.705,0.695,0.703,0.701,0.698,0.699,0.702,0.702,0.704,0.695,0.701,0.696,0.705,0.696,0.706,0.701,0.696,0.693,0.697,0.697,0.694,0.700,0.699,0.698,
		0.699,0.702,0.697,0.692,0.693,0.696,0.697,0.689,0.697,0.683,0.684,0.693,0.687,0.698,0.695,0.702,0.699,0.697,0.697,0.699,0.697,0.699,0.694,0.691,0.699,0.691,0.692,0.694,
		0.694,0.702,0.697,0.699,0.701,0.695,0.698,0.699,0.702,0.698,0.697,0.700,0.702,0.698,0.697,0.698,0.701,0.708,0.703,0.697,0.699,0.698,0.701,0.705,0.704,0.706,0.702,0.700,
		0.704,0.704,0.706,0.712,0.701,0.705,0.704,0.704,0.702,0.698,0.697,0.702,0.703,0.698,0.699,0.693,0.701,0.700,0.701,0.700,0.694,0.697,0.698,0.698,0.699,0.695,0.701,0.704,
		0.698,0.701,0.697,0.701,0.697,0.700,0.694,0.697,0.691,0.687,0.684,0.690,0.692,0.685,0.689,0.687,0.693,0.686,0.694,0.697,0.690,0.691,0.690,0.688,0.695,0.690,0.690,0.690,
		0.690,0.694,0.693,0.693,0.690,0.688,0.692,0.687,0.689,0.692,0.687,0.694,0.695,0.692,0.693,0.693,0.698,0.694,0.695,0.692,0.699,0.694,0.693,0.695,0.695,0.696,0.692,0.694,
		0.690,0.694,0.697,0.698,0.700,0.700,0.700,0.696,0.701,0.696,0.696,0.698,0.704,0.705,0.696,0.694,0.698,0.701,0.702,0.703,0.697,0.701,0.702,0.701,0.698,0.708,0.702,0.702,
		0.701,0.702,0.699,0.706,0.705,0.700,0.714,0.699,0.710,0.705,0.704,0.702,0.704,0.703,0.703,0.703,0.705,0.706,0.703,0.701,0.712,0.701,0.710,0.707,0.709,0.704,0.709,0.713,
		0.706,0.711,0.710,0.709,0.719,0.709,0.705,0.711,0.712,0.711,0.708,0.709,0.709,0.710,0.710,0.716,0.716,0.712,0.709,0.717,0.713,0.715,0.711,0.716,0.716,0.717,0.714,0.715,
		0.711,0.720,0.714,0.718,0.715,0.711,0.720,0.716,0.717,0.722,0.722,0.724,0.720,0.718,0.768,0.770,0.803,0.800,0.800,0.800,0.827,0.826,0.844,0.845,0.861,0.864,0.858,0.858,
		0.875,0.874,0.880,0.880,0.885,0.885,0.890,0.889,0.896,0.893,0.897,0.898,0.901,0.895,0.906,0.909,0.911,0.911,0.912,0.910,0.911,0.910,0.909,0.910,0.908,0.908,0.899,0.885,
		0.885,0.884,0.885,0.881,0.883,0.880,0.882,0.881,0.879,0.879,0.878,0.881,0.879,0.879,0.880,0.880,0.882,0.882,0.877,0.878,0.878,0.880,0.876,0.879,0.880,0.881,0.879,0.878,
		0.878,0.879,0.878,0.880,0.880,0.879,0.878,0.878,0.877,0.878,0.880,0.877,0.878,0.878,0.878,0.881,0.880,0.882,0.883};
	double khodyukLOErrors[numKhodyukPoints]={0.014,0.014,0.015,0.014,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,
		0.015,0.015,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
		0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.015,0.015,0.016,0.016,0.016,0.016,0.017,0.017,0.017,0.017,0.017,0.017,0.017,
		0.017,0.018,0.017,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
		0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
		0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018};
	for(int i=0; i <numKhodyukPoints; i++)
	{
		theFile << "2\t" << khodyukEnergies[i] << "\t" << khodyukLOs[i] << "\t" << khodyukLOErrors[i] << endl;
	}

	const int numMoszynskiPoints=30;
	double moszynskiEnergies[numMoszynskiPoints]={14.4,16.6,22,31,32,45,52,59.6,74,81,84,88,122,244.1,302.9,320.1,335.8,356,411.8,435,493.7,511,569.7,586,662,758.8,834.8,1063.7,1173.2,1275};
	double moszynskiLOs[numMoszynskiPoints]={0.711,0.719,0.775,0.829,0.8412,0.853,0.9,0.893,0.91,0.926,0.937,0.933,0.915,0.964,0.979,0.964,0.97,0.986,0.97,0.995,1,0.996,1,0.999,1,1,1,1.01,1.013,1.01};
	double moszynskiLOErrors[numMoszynskiPoints]={0.01422,0.01438,0.0155,0.01658,0.016824,0.02559,0.018,0.01786,0.0182,0.01852,0.01874,0.01866,0.0183,0.02892,0.01958,0.01928,0.0291,0.01972,0.0194,0.02985,0.03,0.01992,0.02,0.02997,0.02,0.03,0.02,0.0202,0.02026,0.0202};
	for(int i=0; i <numMoszynskiPoints; i++)
	{
		theFile << "3\t" << moszynskiEnergies[i] << "\t" << moszynskiLOs[i] << "\t" << moszynskiLOErrors[i] << endl;
	}

	const int numAverkievPoints=20;
	double averkievEnergies[numAverkievPoints]={2.1,6.4,7.5,8.0,8.8,10.0,11.2,16.6,17.2,23.5,33.1,58.3,71.6,78.3,120.1,162.3,277.4,384.8,506.5,657.9};
	double averkievLOs[numAverkievPoints]={0.854,0.777,0.739,0.706,0.678,0.657,0.641,0.579,0.590,0.685,0.731,0.837,0.844,0.827,0.884,0.903,0.959,0.977,0.986,0.996};
	double averkievLOErrors[numAverkievPoints]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.035,0.035,0.035,0.035,0.035,0.02,0.015,0.01,0.005,0.005};
	for(int i=0; i <numAverkievPoints; i++)
	{
		theFile << "4\t" << averkievEnergies[i] << "\t" << averkievLOs[i] << "\t" << averkievLOErrors[i] << endl;
	}

	theFile.close();

}
