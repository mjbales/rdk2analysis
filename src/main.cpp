#include <RDK2AnalysisPlotter.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>

#include "TRandom3.h"
#include "TTreePlayer.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TVectorT.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"

#include "RDK2Constants.h"
#include "RDK2CutSet.h"
#include "RDK2Analysis.h"
#include "RDK2Set.h"
#include "RDK2MCAnalyzer.h"
#include "RDK2ExpAnalyzer.h"
#include "RDK2Events.h"
#include "RDK2Fierz.h"
#include "RDK2IO.h"
#include "LTimeAnalysis.h"

#include "MRGraphics.h"
#include "MRIO.h"

using namespace std;

int main(int argc, char * const argv[])
{
	delete gRandom;
	gRandom = new TRandom3(0);
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetTitleFont(62);
	gStyle->SetTitleFont(62, "xyz");
	gStyle->SetTitleBorderSize(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadColor(kWhite);
	gStyle->SetPadLeftMargin(.09);
	gStyle->SetPadRightMargin(.05);
	gStyle->SetLabelFont(62);
	gStyle->SetGridColor(12);

	const Int_t font = 42;
	gStyle->SetLegendFont(font);
	gStyle->SetLabelFont(font, "xyz");
	gStyle->SetTitleFont(font, "xyz");
	gStyle->SetTitleFont(font, "");
	gStyle->SetLabelFont(font, "");
	gStyle->SetStatFont(font);
//	gStyle->SetEndErrorSize(100);
//	gStyle->SetErrorX(10);
//	gROOT->ForceStyle();
	const TString HistsDir = HISTS_DIR;
	const TString GraphsDir = GRAPHS_DIR;
	const TString ResultsDir = RESULTS_DIR;
	const TString MRKRunDir = MRKGEANTDATA_DIR;
	const TString EventsDir = EVENTS_DIR;

	RDK2CutSet STD_EPCutSetA("STD_EPCutSetA", "Standard Cuts", 50., 800., 10., 31., 2.e-6, 10.e-6, 10., 800., 0.1, 20.);
	RDK2CutSet STD_EPCutSetB("STD_EPCutSetB", "Standard Cuts, 100 keV eE limit", 100., 800., 10., 31., 2.e-6, 10.e-6, 10., 800., 0.1, 20.);
	RDK2CutSet STD_EPCutSetC("STD_EPCutSetC", "Standard Cuts, 100 keV eE limit, 13 keV pE Limit ", 100., 800., 13., 31., 2.e-6, 10.e-6, 10, 800., 0.1, 20.);  // Higher electron cut
	RDK2CutSet STD_EPCutSetD("STD_EPCutSetD", "Standard Cuts, 100 keV eE limit, 13 keV pE Limit, 10-25us ", 100., 800., 13., 31., 10.e-6, 25.e-6, 10, 800., 0.1, 20.);  // Higher electron cut
	RDK2CutSet STD_EPCutSetE("STD_EPCutSetE", "Standard Cuts, 100 keV eE limit, 13 keV pE Limit, 2-25us ", 100., 800., 13., 31., 2.e-6, 25.e-6, 10, 800., 0.1, 20.);  // Higher electron cut
	RDK2CutSet STD_EPCutSetF("STD_EPCutSetF", "Standard Cuts,20 keV Gamma Cut, 100 keV eE limit, 13 keV pE Limit, 2-25us ", 100., 800., 13., 31., 2.e-6, 25.e-6, 20, 800., 0.1, 20.);  // Higher electron cut
	RDK2CutSet STD_EPCutSetG("STD_EPCutSetG", "Standard Cuts 150420", 100., 757., 13., 31., 2.e-6, 25.e-6, 10, 780., 0.1, 20.);  // Higher electron cut
	RDK2CutSet STD_EPCutSetH("STD_EPCutSetH", "Standard Cuts 151015", 100., 757., 13., 31., 2.e-6, 25.e-6, 10, 780., 0.4, 30.);  // Final Cuts?
	RDK2CutSet STD_EPCutSetI("STD_EPCutSetI", "Standard Cuts 151015", 100., 757., 0, 31., 2.e-6, 25.e-6, 10, 780., 0.4, 30.);  // Final Cuts? (no lower proton cut)
	RDK2CutSet STD_EPCutSetJ("STD_EPCutSetJ", "Standard Cuts 151015", 0.01, 757., 0.01, 31., 2.e-6, 25.e-6, 10, 780., 0.4, 30.);  // Final Cuts? (no lower proton cut or electron cut)
	RDK2CutSet STD_EPCutSetK("STD_EPCutSetK", "Standard Cuts 151015", 0.01, 757., 0.01, 31., 0, 50.e-6, 10, 780., 0.4, 30.);  // Final Cuts? (no lower proton cut or electron cut or time)
	RDK2CutSet STD_EPCutSetL("STD_EPCutSetH", "Standard Cuts 160214", 100., 757., 13., 31., 2.e-6, 25.e-6, 10, 780., 0.4, 14.);  // Final Cuts?
	RDK2CutSet STD_ECutSetA("STD_ECutSetA", "E Cut Only", 50., 800., -REALLY_BIG_DBL, REALLY_BIG_DBL, -REALLY_BIG_DBL,
	REALLY_BIG_DBL, 10., 800., 0.1, 20.);
	RDK2CutSet STD_NoEPCutSetA("STD_NoEPCutSetA", "No EP Cut", -REALLY_BIG_DBL, REALLY_BIG_DBL, -REALLY_BIG_DBL,
	REALLY_BIG_DBL, -REALLY_BIG_DBL, REALLY_BIG_DBL, 10., 800., 0.1, 20.);
	RDK2CutSet CutSetNoTOF("CutSetNoTOF", "No ToF, 100 keV eE limit, 13 keV pE Limit, no TOF ", 100., 800., 13., 31., 0, 30.e-6, 10, 800., 0.1, 20.);  // Higher electron cut

/*AID 83 for BAPDS*/
//
//	double dataToResidRatio = 3;
////	TString imageString = "AID83_BAPD_scaled_norm_var_loglog.pdf";
//	TString imageString = "AID83_BAPD_scaled_norm_var_loglog.root";
//	gStyle->SetErrorX(0);
//	TCanvas theCanvas("largeCanvas", "largeCanvas", 10, 10, 1200, 1200);
//	gPad->SetFillColor(kWhite);
//	TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
//	TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
//	pad1.SetFillColor(kWhite);
//	pad2.SetFillColor(kWhite);
//	pad1.SetLeftMargin(.15);
//	pad2.SetLeftMargin(.15);
//	pad1.SetRightMargin(.05);
//	pad2.SetRightMargin(.05);
//	pad1.SetBottomMargin(0);
//	pad1.SetTopMargin(0);
//	pad2.SetTopMargin(0);
//	pad2.SetBottomMargin(.12);
//	pad1.Draw();
//	pad2.Draw();
//
//	TH1D* expHistInit = getTabSeperatedHist(HistsDir + "bAPD_d236_d237_epg_gEn.txt", BAPD_DIM3);
//	TH1D* mcHistInit = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPBG_EDepBGAvg.txt", BAPD_DIM3);
//
//	vector<double> var_bins={0,.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,2,2.4,3,3.6,4.3,5.1,6.1,7.8,14};
//
//	TH1* expHist=makeVarHistScaledByWidth(expHistInit,var_bins);
//	TH1* mcHist=makeVarHistScaledByWidth(mcHistInit,var_bins);
//	expHistInit->Scale(1./2.51376e+07);
//	expHist->Scale(1./2.51376e+07);
//
//	double integral, error;
//	integral=expHistInit->IntegralAndError(1,expHistInit->GetNbinsX(),error);
//	cout << integral << " +/- " << error << endl;
//	integral=mcHistInit->IntegralAndError(1,mcHistInit->GetNbinsX(),error);
//	cout << integral << " +/- " << error << endl;
//
////	mcHistInit->Rebin(10);
//	mcHistInit->Scale(1e7);
//	mcHistInit->GetXaxis()->SetRangeUser(0.4,14);
//
//
////	expHist->Rebin(4);
////	mcHist->Rebin(4);
////	fitMCtoExp(expHist, mcHist,1,expHist->GetNbinsX());
//	scaleHistogramsTogether(expHist, mcHist);
//
//	mcHistInit->SetLineColor(kGreen+3);
//	//mcHistInit->SetLineStyle(2);
//	mcHistInit->SetLineWidth(3);
//
//	expHist->SetMarkerStyle(4);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(3);
//
//	mcHist->Print();
//	expHist->Print();
//
//
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//
//
////	mcHist->Rebin(10);
////	expHist->Rebin(10);
//
//	pad2.cd();
//	pad2.SetLogx(1);
//	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	pad1.SetLogx(1);
//
//
//	TH1F* hr = pad2.DrawFrame(.3, 0.095, 15, expHist->GetMaximum() * 1.5);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
////	hr->GetXaxis()->SetMoreLogLabels(1);
////	hr->GetYaxis()->SetMoreLogLabels(1);
//	hr->GetYaxis()->SetNoExponent(1);
//	hr->GetXaxis()->SetNoExponent(1);
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";APD photon peak height (5.9=5.9 keV);#frac{APD ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHistInit->Draw("HIST 9 same");
//
//	expHist->Draw("E1 9 same");
//	cout << mcHistInit->GetBinLowEdge(5) << endl;
//	cout << mcHistInit->GetBinContent(5) << endl;
//	cout << mcHistInit->GetBinError(5) << endl;
//
//
//
//	TLegend theLegend(0.23, .32, 0.665, 0.52);
//	theLegend.SetBorderSize(1);
//	theLegend.SetTextSize(0.06);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(mcHistInit, "Theory/simulation", "L");
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist);
//
//	residualHist->SetMarkerColor(kGreen+3);
//	residualHist->SetLineColor(kGreen+3);
//	residualHist->GetXaxis()->SetLabelSize(0.2);
//	residualHist->GetYaxis()->SetLabelSize(0.2);
//	residualHist->GetXaxis()->SetTitleSize(0.2);
//	residualHist->GetYaxis()->SetTitleSize(0.2);
//	//residualHist->GetXaxis()->SetNdivisions(0,0,0);
//
//	pad1.cd();
//	pad1.SetGrid();
//	TH1F* hr1 = pad1.DrawFrame(.3, -1.4, 15, 1.4);
//
//
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05 * dataToResidRatio);
//	hr1->SetTitleOffset(.41, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05 * dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
//	hr1->SetTitle(";;#frac{Exp. - Sim}{Sim.}");
//
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	residualHist->Draw("E1 same");
//
//	TLegend theLegend2(0.45, .06, 0.92, 0.32);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetTextSize(0.16);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Normalized residuals", "ep");
//	theLegend2.Draw("same");
//
//	theCanvas.SaveAs(GraphsDir + imageString);


/* BGOS */
//	double dataToResidRatio = 3;
//	bool absScale=true;
//	TString imageString = "AID83_BGO_scaled_norm_var_loglog_abs.root";
////	TString imageString = "AID83_BGO_scaled_norm_var_loglog_abs.pdf";
//	gStyle->SetErrorX(0);
//	TCanvas theCanvas("AID83_BGO_scaled_norm_var_loglog_abs", "AID83_BGO_scaled_norm_var_loglog_abs", 10, 10, 1200, 1200);
//	gPad->SetFillColor(kWhite);
//	TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
//	TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
//	pad1.SetFillColor(kWhite);
//	pad2.SetFillColor(kWhite);
//	pad1.SetLeftMargin(.15);
//	pad2.SetLeftMargin(.15);
//	pad1.SetRightMargin(.05);
//	pad2.SetRightMargin(.05);
//	pad1.SetBottomMargin(0);
//	pad1.SetTopMargin(0);
//	pad2.SetTopMargin(0);
//	pad2.SetBottomMargin(.12);
//	pad1.Draw();
//	pad2.Draw();
//
////	vector<double> var_bins={0,10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,400};
//	vector<double> var_bins={0,10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,200,800};
//	TH1D* expHistInit = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_gEn.txt", BGO_DIM2);
//	double correction=-0.014;
//	TH1D* mcHistInit = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_EDepGAvg.txt", BGO_DIM2);
//	TH1D* theoryHist = getTabSeperatedHist(HistsDir + "EG0_BGODIM_Hist.txt", BGO_DIM2);
//
//
//	TH1* expHist=makeVarHistScaledByWidth(expHistInit,var_bins);
//	TH1* mcHist=makeVarHistScaledByWidth(mcHistInit,var_bins);
//
//	//Exp Scaled
//	expHist->Scale((1.-correction)/(2.31112e+08));
//	expHist->Scale(1e6);
//	expHistInit->Scale((1.-correction)/(2.31112e+08));
//
//	double error;
//	double integral=mcHistInit->IntegralAndError(1,mcHistInit->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	integral=expHistInit->IntegralAndError(1,expHistInit->GetNbinsX(),error);
//	cout << "Exp: " << integral << " +/- " << error << endl;
//	expHistInit->Scale(1e6);
//
//	if(!absScale)
//	{
//		scaleHistogramsTogether(expHist,mcHist);
//		scaleHistogramsTogether(expHistInit,mcHistInit);
//	}
//	else
//	{
//		mcHistInit->Scale(1e6);
//		mcHist->Scale(1e6);
//	}
//	mcHistInit->Rebin(10);
//	mcHistInit->GetXaxis()->SetRangeUser(10,1000);
//
//	scaleHistogramsTogether(mcHistInit,theoryHist);
//	theoryHist->Rebin(10);
//	theoryHist->GetXaxis()->SetRangeUser(10,1000);
//
//
//	mcHistInit->SetLineColor(kRed);
//	theoryHist->SetLineColor(kBlue);
//	theoryHist->SetLineStyle(2);
//	theoryHist->SetLineWidth(3);
//	mcHistInit->SetLineWidth(3);
////
//	expHist->SetMarkerStyle(4);
//	expHist->SetMarkerSize(3);
//	expHist->SetLineWidth(3);
//
//	pad2.cd();
//	pad2.SetLogx(1);
//	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	pad1.SetLogx(1);
//	TH1F* hr = pad2.DrawFrame(9, 1e-3, 1010, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
////	hr->GetXaxis()->SetMoreLogLabels(1);
////	hr->GetYaxis()->SetMoreLogLabels(1);
//	hr->GetYaxis()->SetNoExponent(1);
//	hr->GetXaxis()->SetNoExponent(1);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetXaxis()->SetTitleOffset(1.2);
//	hr->GetYaxis()->SetTitleOffset(1.4);
//	hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";BGO photon peak height (511=511 keV);#frac{BGO ep#gamma counts}{1 million ep counts} (per keV)");
//	pad2.cd();
////	mcHist->Draw("HIST 9 same");
//	theoryHist->Draw("HIST 9 same");
//	mcHistInit->Draw("HIST 9 same");
//	expHist->Draw("E1 9 same");
//
//	TLegend theLegend(0.23, .32, 0.665, 0.52);
//	theLegend.SetTextSize(0.06);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(theoryHist, "Theory", "L");
//	theLegend.AddEntry(mcHistInit, "Theory+sim.", "L");
//	theLegend.AddEntry(expHist, "Experiment", "E1P");
//
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist);
//
//	residualHist->SetMarkerColor(kRed);
//	residualHist->SetLineColor(kRed);
//	residualHist->GetXaxis()->SetLabelSize(0.2);
//	residualHist->GetYaxis()->SetLabelSize(0.2);
//	residualHist->GetXaxis()->SetTitleSize(0.2);
//	residualHist->GetYaxis()->SetTitleSize(0.2);
//	//residualHist->GetXaxis()->SetNdivisions(0,0,0);
//
//	pad1.cd();
//	pad1.SetLogy(0);
//	pad1.SetGrid();
//	TH1F* hr1 = pad1.DrawFrame(9, -0.25, 1010, 0.45);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.45, "y");
//	TLine lineToDraw(9,0,1010,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.SetLineWidth(2);
//	lineToDraw.Draw("same");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
//	hr1->SetTitle(";; #frac{Exp.- Sim.}{ Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	residualHist->Draw("E1 9 same");
//
//
//
//	TLegend theLegend2(0.42, .06, 0.94, 0.32);
//	theLegend2.SetTextSize(0.16);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Normalized residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/*Make EG0 Hist */
//
//	RDK2MCAnalysis mcAnal("AID83",  "Final?? (fixed BAPD cuts)", 48,  "RID518", "RID519",  48, "RID520",  "RID521", "RID522",  STD_EPCutSetH);
//	TH1* theoryHist=mcAnal.MakeHist("eg0",FOURBODY,"eg0","",BGO_DIM2);
//	convertTH1ToTXT(theoryHist,HistsDir+"EG0_BGODIM_Hist.txt");
//	delete theoryHist;
//	theoryHist=mcAnal.MakeHist("eg0",FOURBODY,"eg0","",BAPD_DIM3);
//	convertTH1ToTXT(theoryHist,HistsDir+"EG0_BAPDDIM_Hist.txt");
//	delete theoryHist;

	/* BAPDs */
//	double dataToResidRatio = 3;
////	TString imageString = "AID83_BAPD_scaled_norm_var_loglog.root";
//	TString imageString = "AID83_BAPD_scaled_norm_var_loglog.pdf";
//	gStyle->SetErrorX(0);
//	TCanvas theCanvas("AID83_BAPD_scaled_norm_var_loglog", "AID83_BAPD_scaled_norm_var_loglog", 10, 10, 1200, 1200);
//	gPad->SetFillColor(kWhite);
//	TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
//	TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
//	pad1.SetFillColor(kWhite);
//	pad2.SetFillColor(kWhite);
//	pad1.SetLeftMargin(.15);
//	pad2.SetLeftMargin(.15);
//	pad1.SetRightMargin(.05);
//	pad2.SetRightMargin(.05);
//	pad1.SetBottomMargin(0);
//	pad1.SetTopMargin(0);
//	pad2.SetTopMargin(0);
//	pad2.SetBottomMargin(.12);
//	pad1.Draw();
//	pad2.Draw();
//
//	vector<double> var_bins={0,.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,2,2.4,3,3.6,4.3,5.1,6.1,7.8,14};
//	TH1D* expHistInit = getTabSeperatedHist(HistsDir + "bAPD_d236_d237_epg_gEn.txt", BAPD_DIM3);
//	double correction=-0.003;
//	double epCount=2.51376e+07;
//	TH1D* mcHistInit = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPBG_EDepBGAvg.txt", BAPD_DIM3);
//	TH1D* theoryHist = getTabSeperatedHist(HistsDir + "EG0_BAPDDIM_Hist.txt", BAPD_DIM3);
//
//
//	TH1* expHist=makeVarHistScaledByWidth(expHistInit,var_bins);
//	TH1* mcHist=makeVarHistScaledByWidth(mcHistInit,var_bins);
//
//	//Exp Scaled
//	expHist->Scale((1.+correction)/(epCount));
//	expHist->Scale(1e6);
//	expHistInit->Scale((1.+correction)/(epCount));
//
//	double error, integral;
//	integral=mcHist->IntegralAndError(1,140,error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	integral=expHist->IntegralAndError(1,140,error);
//	cout << "Exp: " << integral << " +/- " << error << endl;
//	integral=mcHistInit->IntegralAndError(1,140,error);
//	cout << "MC Init: " << integral << " +/- " << error << endl;
//	integral=expHistInit->IntegralAndError(1,140,error);
//	cout << "Exp Init: " << integral << " +/- " << error << endl;
//
//	expHistInit->Scale(1e6);
//
//
//	scaleHistogramsTogether(expHist,mcHist,1,140);
//	scaleHistogramsTogether(expHistInit,mcHistInit,1,140);
//	mcHistInit->Scale(10);
//	mcHistInit->GetXaxis()->SetRangeUser(0.4,14);
//
//	scaleHistogramsTogether(mcHistInit,theoryHist,1,140);
//	theoryHist->GetXaxis()->SetRangeUser(0.4,14);
//
//
//
//	mcHistInit->SetLineColor(kGreen+2);
//	theoryHist->SetLineColor(kBlue);
//	theoryHist->SetLineStyle(2);
//	theoryHist->SetLineWidth(3);
//	mcHistInit->SetLineWidth(3);
////
//	expHist->SetMarkerStyle(4);
//	expHist->SetMarkerSize(3);
//	expHist->SetLineWidth(3);
//
//	pad2.cd();
//	pad2.SetLogx(1);
//	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	pad1.SetLogx(1);
//	TH1F* hr = pad2.DrawFrame(0.3, 0.095, 15, theoryHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
////	hr->GetXaxis()->SetMoreLogLabels(1);
////	hr->GetYaxis()->SetMoreLogLabels(1);
//	hr->GetYaxis()->SetNoExponent(1);
//	hr->GetXaxis()->SetNoExponent(1);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetXaxis()->SetTitleOffset(1.2);
//	hr->GetYaxis()->SetTitleOffset(1.4);
//	hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";APD photon peak height (5.9=5.9 keV);#frac{APD ep#gamma counts}{1 million ep counts} (per keV)");
//	pad2.cd();
//
//	theoryHist->Draw("HIST 9 same");
//	mcHistInit->Draw("HIST 9 same");
//	expHist->Draw("E1 9 same");
//
//	TLegend theLegend(0.23, .32, 0.665, 0.52);
//	theLegend.SetTextSize(0.06);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(theoryHist, "Theory", "L");
//	theLegend.AddEntry(mcHistInit, "Theory+sim.", "L");
//	theLegend.AddEntry(expHist, "Experiment", "E1P");
//
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist);
//
//	residualHist->SetMarkerColor(kGreen+2);
//	residualHist->SetLineColor(kGreen+2);
//	residualHist->GetXaxis()->SetLabelSize(0.2);
//	residualHist->GetYaxis()->SetLabelSize(0.2);
//	residualHist->GetXaxis()->SetTitleSize(0.2);
//	residualHist->GetYaxis()->SetTitleSize(0.2);
//	//residualHist->GetXaxis()->SetNdivisions(0,0,0);
//
//	pad1.cd();
//	pad1.SetLogy(0);
//	pad1.SetGrid();
//	TH1F* hr1 = pad1.DrawFrame(0.3, -1.1, 15, 1.1);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.45, "y");
//	TLine lineToDraw(0.3,0,15,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.SetLineWidth(2);
//	lineToDraw.Draw("same");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
//	hr1->SetTitle(";; #frac{Exp.- Sim.}{ Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	residualHist->Draw("E1 9 same");
//
//
//
//	TLegend theLegend2(0.42, .06, 0.94, 0.32);
//	theLegend2.SetTextSize(0.16);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Normalized residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

/* Make Kevin graph without incorrect gamma cut*/
//	RDK2MCAnalysis* theAnal = new RDK2MCAnalysis( "AID83",  "Final?? (fixed BAPD cuts)", 48,  "RID518", "RID519",  48, "RID520",  "RID521", "RID522",  STD_EPCutSetH);
//	for (int i = 0; i < 12; ++i)
//	{
//		TH1* theHist=theAnal->MakeHist("BGO"+int2str(i+1),FOURBODY,"gChn.BGO"+int2str(i+1)+"EDepTotal",STD_EPCutSetC.CreateMCEPCut(1,0),BGOCAL_DIM);
//		convertTH1ToTXT(theHist,HistsDir+"HistAID83_BGO"+int2str(i+1)+"EDepTotal.txt");
//		delete theHist;
//	}
//	TCanvas theCanvas("theCanvas", "theCanvas", 1200, 800);
//	for (int i = 0; i < 12; ++i)
//	{
//		TH1D* theHist=getTabSeperatedHist(HistsDir+"HistAID83_BGO"+int2str(i+1)+"EDepTotal.txt",BGOCAL_DIM);
//		theHist->Rebin(10);
//		theHist->GetYaxis()->SetRangeUser(0,5000);
//		theHist->SetBinContent(1,0);
//		theHist->Draw("L");
//		theCanvas.SaveAs(GraphsDir+"HistAID83_BGO"+int2str(i+1)+"EDepTotal.png");
//
//		delete theHist;
//	}

/* Try Herbert's idea on x shifts */
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID61_EP_ToFP.txt", TOFP_DIM2);
//	TH1D* mcHistXPlus = getTabSeperatedHist(HistsDir + "Hist_MC_AID64_EP_ToFP.txt", TOFP_DIM2);
//	TH1D* mcHistXNeg = getTabSeperatedHist(HistsDir + "Hist_MC_AID65_EP_ToFP.txt", TOFP_DIM2);
//	mcHist->Sumw2();
//	mcHistXPlus->Sumw2();
//	mcHistXNeg->Sumw2();
//
//
//
//
//	mcHist->Scale(10./48.);
//
//	int rebin=9;
//	mcHist->Rebin(rebin);
//	mcHistXPlus->Rebin(rebin);
//	mcHistXNeg->Rebin(rebin);
//
//	TH1D* mcHistXPlusResid=(TH1D*) mcHistXPlus->Clone();
//	mcHistXPlusResid->Add(mcHist,-1);
//	mcHistXPlusResid->Divide(mcHist);
//
//	TH1D* mcHistXNegResid=(TH1D*) mcHistXNeg->Clone();
//	mcHistXNegResid->Add(mcHist,-1);
//	mcHistXNegResid->Divide(mcHist);
//
//	mcHistXPlusResid->SetTitle("Beam +1mm X Dir;Proton Time of Flight (s);Normalized Residuals");
//	mcHistXNegResid->SetTitle("Beam -1mm X Dir;Proton Time of Flight (s);Normalized Residuals");
//
//
//	TH1* expHists[2]={mcHistXPlusResid,mcHistXNegResid};
//
//	expHists[0]->GetYaxis()->SetRangeUser(-.06,.06);
//	expHists[0]->GetXaxis()->SetRangeUser(0,10.e-6);
//
//	TCanvas* theCanvas=plotExpVersusMC(2, expHists,0,NULL,"X Shifts");
//
//	theCanvas->SaveAs(GraphsDir+"HerbertToFXShiftIdea.png");

//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID523", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 20 nm"); //“Final” Run (1) with 20 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID524", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 20 nm"); //“Final” Run (1) with 20 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID525", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 20 nm"); //“Final” Run (1) with 20 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID526", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 20 nm"); //“Final” Run (1) with 20 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID527", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 20 nm"); //“Final” Run (1) with 20 nm dead layer

/*Mess with fitting */
//	vector<double> var_bins={0,10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,400,800};
//	TH1D* expHistInit = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_gEn.txt", BGO_DIM2);
//	double correction=-0.014;
//	TH1D* mcHistInit = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_EDepGAvg.txt", BGO_DIM2);
//	TH1D* theoryHist = getTabSeperatedHist(HistsDir + "EG0_BGODIM_Hist.txt", BGO_DIM2);
//
//
//	TH1* expHist=makeVarHistScaledByWidth(expHistInit,var_bins);
//	TH1* mcHist=makeVarHistScaledByWidth(mcHistInit,var_bins);
//
////	TH1* expHist = expHistInit;
////	TH1* mcHist = mcHistInit;
//	mcHist->SetName("mcHist");
//
//	//Exp Scaled
//	expHist->Scale((1.-correction)/(2.31112e+08));
////	expHist->Scale(1e6);
////	expHistInit->Scale((1.-correction)/(2.31112e+08));
//
//	double error;
//	double integral=mcHistInit->IntegralAndError(1,mcHistInit->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	integral=expHistInit->IntegralAndError(1,expHistInit->GetNbinsX(),error);
//	cout << "Exp: " << integral << " +/- " << error << endl;
////	expHistInit->Scale(1e6);
//
////
////	mcHistInit->Scale(1e6);
////	mcHist->Scale(1e6);
//
//	TF1* result=fitMCtoExp(expHist, mcHist,1,expHist->GetNbinsX());
//
//	cout << "NDF:" << result->GetNDF() << endl;
//	cout << "Chi Square:" << result->GetChisquare() << endl;
//	cout << "Prob:" << result->GetProb() << endl;


//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_d239_d240_epg_gEn.txt", -0.014, 2.31112e+08, HistsDir + "bAPD_d236_d237_epg_gEn.txt", -0.003, 2.51376e+07, 83,true,true,false,true,false);
//	makeEGTimingPlot(HistsDir+"ExpJeff_epg_egTOF_151122.txt",HistsDir+"ExpJeff_epbg_ebgTOF_151122.txt");

//	cout << gRandom->Rndm()<<endl;
//	cEvents theEvents;
//	double a=theEvents.calcBR(5,781,0,10000000,false, 885.7);
//	double b=theEvents.calcBR(5,781,0,10000000,true, 885.7);
//	double c=theEvents.calcBR(5,781,0,10000000,false,880.3);
//	double d=theEvents.calcBR(5,781,0,10000000,true, 880.3);
//	cout << "5 keV to 782 keV, Fermi Off, 885.7 s lifetime: " << a << endl;
//	cout << "5 keV to 782 keV, Fermi On, 885.7 s lifetime: " << b << endl;
//	cout << "5 keV to 782 keV, Fermi Off, 880.3 s lifetime: " << c << endl;
//	cout << "5 keV to 782 keV, Fermi On, 880.3 s lifetime: " << d << endl;
//
//	cout << a << endl;
//	cout << b << endl;
//	cout << c << endl;
//	cout << d << endl;
//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_d239_d240_epg_gEn.txt", -0.014, 2.31112e+08, HistsDir + "bAPD_d236_d237_epg_gEn.txt", -0.003, 2.51376e+07, 83,true,true,false,true,false);


/* Remake nonprop for most recent cuts and analysis */
//
//		 RDK2MCAnalysis mcAnalysis( "AID83",  "Final?? (fixed BAPD cuts)", 48,  "RID518", "RID519",  48, "RID520",  "RID521", "RID522",  STD_EPCutSetH);
//
//
//		 TH1* normHist=mcAnalysis.MakeRadAvgHist("Normal",DET_EPG,"",0,BGO_DIM2);
//		 TH1* highHist=mcAnalysis.MakeRadAvgHist("High",DET_EPG,"",2,BGO_DIM2);
//		 TH1* lowHist=mcAnalysis.MakeRadAvgHist("Low",DET_EPG,"",3,BGO_DIM2);
//
//		 convertTH1ToTXT(normHist,HistsDir+"ResultsParamModel3Norm.txt");
//		 convertTH1ToTXT(highHist,HistsDir+"ResultsParamModel3High.txt");
//		 convertTH1ToTXT(lowHist,HistsDir+"ResultsParamModel3Low.txt");
//
//	//	 TH1* normHist=getTabSeperatedHist(HistsDir+"ResultsParamModel2NormFixed.txt",BGO_DIM2);
//	//	 TH1* highHist=getTabSeperatedHist(HistsDir+"ResultsParamModel2HighFixed.txt",BGO_DIM2);
//	//	 TH1* lowHist=getTabSeperatedHist(HistsDir+"ResultsParamModel2LowFixed.txt",BGO_DIM2);
//
//		 normHist->Rebin(50);
//		 highHist->Rebin(50);
//		 lowHist->Rebin(50);
//
//		 TH1* theHists2[3]={normHist,highHist,lowHist};
//
//		 double integrals[3]={normHist->Integral(),highHist->Integral(),lowHist->Integral()	 };
//
//		 theHists2[0]->SetTitle("Norm: "+d2str(integrals[0]));
//		 theHists2[1]->SetTitle("High: +"+d2str((integrals[1]-integrals[0])*100/integrals[0])+"%");
//		 theHists2[2]->SetTitle("Low: "+d2str((integrals[2]-integrals[0])*100/integrals[0])+"%");
//
//		 plotExpVersusMCToImage(3,theHists2,0,NULL,"Third model of non-proportionality: 10 keV Limit","ResultsNonPropModel3_10keVLimit.pdf");


	/*Final???? analysis 151216 */
//	createMCAnalysisAndPlot( "AID85",  "Final??? (fixed SiDioxideLayer)", 48,  "RID523", "RID524",  48, "RID525",  "RID526", "RID527",  STD_EPCutSetH);
//	RDK2AnalysisPlotter comPlot("Final???","","AID85");
//	comPlot.makeAll();

	/*Make plots*/
//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	int mcID=83;
//	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID);

	/* Dead Layer figuring */
//	TH1D* mcHist = getTabSeperatedHist(TString(HISTS_DIR) + "BGO_d239_d240_ep_pEn.txt", EDEPE_DIM2);
//
////	RDK2MCAnalysis mcAnalysis( "AID85",  "Final??? (fixed SiDioxideLayer)", 1,  "RID528", "RID524",  10, "RID525",  "RID526", "RID527",  STD_EPCutSetH);
////	TH1* mcHist=mcAnalysis.MakeHist("mcHist",THREEBODY,"pDChn.SBDEDepBlur","pDChn.SBDEDepBlur > 0 ",EDEPE_DIM2);
////	TH1* mcHist=mcAnalysis.MakeStandardHist("mcHist",DET_EP,"pDChn.SBDEDepBlur",EDEPE_DIM2);
//
//
////	TFile f(ResultsDir+"Results_test2_P_STD_4B_100eV_Endpoint_130429_Set0.root");
////
////	TString resultsFilePath=ResultsDir+"Results_RID528_P_STD_3B_130429_Set0.root";
////	TString dervFilePath=ResultsDir+"Results_RID528_P_STD_3B_130429_Derv_Set0.root";
//////	TString pathToResultsFile=MRKRunDir+"Results_DEADLTest22c_P_STD_3B_130429_Derv_Set0.root";
////	TFile resultsFile(resultsFilePath);
////	TTree* resultsTree=(TTree*) resultsFile.Get("g");
////	createDerivedResultsFile(0, resultsTree,dervFilePath,10000);
////
////	TFile dervFile(dervFilePath);
////	TTree* dervTree=(TTree*) dervFile.Get("rd");
////	dervTree->Draw("SBDEDepBlur>>mcHist","SBDEDepBlur > 0");
////	TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//
//	cout << mcHist->GetMean() << endl;
//	TF1* func = new TF1("f1","gaus",17,23);
//	func->SetParameters(500,mcHist->GetMean(),mcHist->GetRMS());
//	mcHist->Fit("f1","V");


//	makeZ0Hist();
//
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID528", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 25 nm"); //Final run with 25 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID529", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 25 nm"); //Final run with 25 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID530", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 25 nm"); //Final run with 25 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID531", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 25 nm"); //Final run with 25 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID532", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 25 nm"); //Final run with 25 nm dead layer

/* AID86 */
	// createMCAnalysisAndPlot( "AID86",  "Final???? (Geant4 version update and 25 nm)", 48,  "RID528", "RID524",  48, "RID530",  "RID526", "RID532",  STD_EPCutSetH);
//	createMCAnalysisAndPlot( "AID86",  "Final???? (Geant4 version update and 25 nm)", 48,  "RID528", "RID529",  48, "RID530",  "RID531", "RID532",  STD_EPCutSetH);
/*Make plots*/
//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	double bgoEPCounts=2.31112e+08;
//	double apdEPCounts=2.51376e+07;
//	int mcID=86;
//	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID);
//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_"+bgoExpID+"_epg_gEn.txt", bgoCorrection, bgoEPCounts, HistsDir + "bAPD_"+apdExpID+"_epg_gEn.txt",apdCorrection, apdEPCounts, mcID,true,true,false,true,false);
//
//	makeEGTimingPlot(HistsDir+"ExpJeff_epg_egTOF_151122.txt",HistsDir+"ExpJeff_epbg_ebgTOF_151122.txt");

	/*Check Columb correction for BAPDs*/
//	cEvents theEvents;
//	double a=theEvents.calcBR(14.1,781,0,10000000,false,880.3);
//	double b=theEvents.calcBR(14.1,781,0,10000000,true, 880.3);
//	double c=theEvents.calcBR(0.4,14,0,10000000,false,880.3);
//	double d=theEvents.calcBR(0.4,14,0,10000000,true, 880.3);
//	cout << "14.1 keV to 782 keV, Fermi Off, 880.3s lifetime: " << a << endl;
//	cout << "14.1 keV to 782 keV, Fermi On, 880.3s lifetime: " << b << endl;
//	cout << "0.4 keV to 14 keV, Fermi Off, 880.3s lifetime: " << c << endl;
//	cout << "0.4 keV to 14 keV, Fermi On, 880.3 s lifetime: " << d << endl;

//	cEvents theEvents;
//	double a=theEvents.calcBR(0.1,781,0,10000000,true,880.3);
//	cout << "0.1 keV to 782 keV, Fermi Off, 880.3s lifetime: " << a << endl;


	/*Check new vs old analysis 160126*/
//	RDK2MCAnalysis* anals[2];
//	int numFiles=1;
//
//	anals[0]=new RDK2MCAnalysis( "AID83",  "Final?? (fixed BAPD cuts)", numFiles,  "RID518", "RID519",  numFiles, "RID520",  "RID521", "RID522",  STD_EPCutSetK);
//	anals[1]=new RDK2MCAnalysis("AID86",  "Final???? (Geant4 version update and 25 nm)", numFiles,  "RID528", "RID529",  numFiles, "RID530",  "RID531", "RID532",  STD_EPCutSetK);
//	double epCount[2],epCountError[2];
//	double ep4Count[2],ep4CountError[2];
//	double epgCount[2],epgCountError[2];
//	double epbgCount[2],epbgCountError[2];
//	double eCount[2],eCountError[2];
//	double pCount[2],pCountError[2];
//	double e4Count[2],e4CountError[2];
//	double p4Count[2],p4CountError[2];
//	double gCount[2],gCountError[2];
//	double bgCount[2],bgCountError[2];
//	for(int i=0;i<2;i++)
//	{
//		eCount[i]=anals[i]->CalcEPer3Decay(eCountError[i]);
//		pCount[i]=anals[i]->CalcPPer3Decay(pCountError[i]);
//		e4Count[i]=anals[i]->CalcEPer4Decay(e4CountError[i]);
//		p4Count[i]=anals[i]->CalcPPer4Decay(p4CountError[i]);
//		gCount[i]=anals[i]->CalcGPer4Decay(DET_G, gCountError[i]);
//		bgCount[i]=anals[i]->CalcGPer4Decay(DET_BG, bgCountError[i]);
//
//		epCount[i]=anals[i]->CalcEPPer3Decay(epCountError[i]);
//		ep4Count[i]=anals[i]->CalcEPPer4Decay(ep4CountError[i]);
//		epgCount[i]=anals[i]->CalcEPGPer4Decay(DET_EPG, epgCountError[i]);
//		epbgCount[i]=anals[i]->CalcEPGPer4Decay(DET_EPBG, epbgCountError[i]);
//	}
//
//	TString titles[2]={"Old","New"};
//	for(int i=0;i<2;i++)
//	{
//		cout << titles[i] + " E count Per 3 Decay:   " +d2str(eCount[i]) +" +/- "+d2str(eCountError[i]) << endl;
//		cout << titles[i] + " P count Per 3 Decay:   " +d2str(pCount[i]) +" +/- "+d2str(pCountError[i])<< endl;
//		cout << titles[i] + " EP count Per 3 Decay:   " +d2str(epCount[i]) +" +/- "+d2str(epCountError[i])<< endl;
//		cout << titles[i] + " E count Per 4 Decay:   " +d2str(e4Count[i]) +" +/- "+d2str(e4CountError[i]) << endl;
//		cout << titles[i] + " P count Per 4 Decay:   " +d2str(p4Count[i]) +" +/- "+d2str(p4CountError[i])<< endl;
//		cout << titles[i] + " BGO G count Per 4 Decay:   " +d2str(gCount[i]) +" +/- "+d2str(gCountError[i])<< endl;
//		cout << titles[i] + " BAPD G count Per 4 Decay:   " +d2str(bgCount[i]) +" +/- "+d2str(bgCountError[i])<< endl;
//		cout << titles[i] + " EP count Per 4 Decay:   " +d2str(ep4Count[i]) +" +/- "+d2str(ep4CountError[i])<< endl;
//		cout << titles[i] + " BGO EPG count Per 4 Decay:   " +d2str(epgCount[i]) +" +/- "+d2str(epgCountError[i])<< endl;
//		cout << titles[i] + " BAPD EPG count Per 4 Decay:   " +d2str(epbgCount[i]) +" +/- "+d2str(epbgCountError[i])<< endl;
//
//
//	}

	/*Make plots OLD*/
//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	double bgoEPCounts=2.31112e+08;
//	double apdEPCounts=2.51376e+07;
//	int mcID=83;
////	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID);
//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_"+bgoExpID+"_epg_gEn.txt", bgoCorrection, bgoEPCounts, HistsDir + "bAPD_"+apdExpID+"_epg_gEn.txt",apdCorrection, apdEPCounts, mcID,true,true,false,true,false);


//	makeZ0Hist();

	/* Dead Layer figuring */

//	const int numFiles=6;
//	int deadLayers[numFiles]={0,6,12,16,20,24};
//	double means[numFiles];
//	double fitMeans[numFiles];
//	double counts[numFiles];
//	TString geantVer="Geant4.10.02";
//	for(int i =0;i<numFiles;i++)
//	{
//		TString dervFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_"+geantVer+"_Derv.root";
//
//		TString resultsFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_"+geantVer+".root";
//		TFile resultsFile(resultsFilePath);
//		TTree* resultsTree=(TTree*) resultsFile.Get("g");
//		createDerivedResultsFile(0, resultsTree,dervFilePath,100000);
//
//		TFile dervFile(dervFilePath);
//		TTree* dervTree=(TTree*) dervFile.Get("rd");
//		dervTree->Draw("SBDEDepBlur>>mcHist","SBDEDepBlur > 0");
//		TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//		means[i]=mcHist->GetMean();
//		counts[i]=mcHist->Integral(1,mcHist->GetNbinsX());
//		TF1* func = new TF1("f1","gaus",17,23);
//		func->SetParameters(500,mcHist->GetMean(),mcHist->GetRMS());
//		mcHist->Fit("f1","V");
//		fitMeans[i]=func->GetParameter(1);
//	}
//
//	for(int i =0;i<numFiles;i++)
//	{
//		cout << deadLayers[i] << "\t" << counts[i] << "\t" << means[i] << "\t" << fitMeans[i] << endl;
//	}

//	const int numFiles=5;
//	int deadLayers[numFiles]={6,12,16,20,24};
//	TString geantVers[2]={"Geant4.9.6.p02","Geant4.9.6.p03"}
//	for(int i =0;i<numFiles;i++)
//	{
//
//		TFile dervFile(dervFilePath);
//		TTree* dervTree=(TTree*) dervFile.Get("rd");
//		dervTree->Draw("SBDEDepBlur>>mcHist","SBDEDepBlur > 0");
//		TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//		means[i]=mcHist->GetMean();
//		counts[i]=mcHist->Integral(1,mcHist->GetNbinsX());
//		TF1* func = new TF1("f1","gaus",17,23);
//		func->SetParameters(500,mcHist->GetMean(),mcHist->GetRMS());
//		mcHist->Fit("f1","V");
//		fitMeans[i]=func->GetParameter(1);
//	}
//
//	for(int i =0;i<numFiles;i++)
//	{
//		cout << deadLayers[i] << "\t" << counts[i] << "\t" << means[i] << "\t" << fitMeans[i] << endl;
//	}

//	cEvents theEvents;
//
//	theEvents.makeEventFile(MRKRunDir+"Events_Test.root",MRKRunDir+"Fluxmap_RDK2_Matt1.txt",100000,THREEBODY);
//
//	TString resultsFilePath[2]={MRKRunDir+"Results_160202DL20_Geant4.9.6.p02.root",
//								MRKRunDir+"Results_160202DL20_Geant4.9.6.p03.root"};
//	TFile* resultsFiles[2];
//	resultsFiles[0]=new TFile(resultsFilePath[0]);
//	resultsFiles[1]=new TFile(resultsFilePath[1]);
//
//	TTree* resultsTrees[2];
//	resultsTrees[0]=(TTree*) resultsFiles[0]->Get("g");
//	resultsTrees[1]=(TTree*) resultsFiles[1]->Get("g");
//
//	TString dervFilePath[2]={MRKRunDir+"Results_160202DL20_Geant4.9.6.p02_Derv.root",
//							MRKRunDir+"Results_160202DL20_Geant4.9.6.p03_Derv.root"};
//	TFile* dervFiles[2];
//	dervFiles[0]=new TFile(dervFilePath[0]);
//	dervFiles[1]=new TFile(dervFilePath[1]);
//
//	TTree* dervTrees[2];
//	dervTrees[0]=(TTree*) dervFiles[0]->Get("rd");
//	dervTrees[1]=(TTree*) dervFiles[1]->Get("rd");
//
//	TFile* eventFile=new TFile(MRKRunDir+"Events_Test.root");
//
//	TTree* eventTree=(TTree*) eventFile->Get("t");

//	eventTree->AddFriend(resultsTrees[0],"rp02");
//	eventTree->AddFriend(resultsTrees[1],"rp03");
//	eventTree->AddFriend(dervTrees[0],"dp02");
//	eventTree->AddFriend(dervTrees[1],"dp03");

//	TCanvas c2("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,800);
//
//	eventTree->Draw("sqrt(x0*x0+y0*y0)*1000");
//
//
//	c2.SaveAs(GraphsDir+"test.pdf");

//	cEvents theEvents;
//
//	theEvents.loadFluxFileMap(MRKRunDir+"Fluxmap_RDK2_Matt1.txt");
//
//	TH1D hist("test",";Radius;Counts",100,0,20);
//	double x,y,r;
//	for(int i=0;i<100000;i++)
//	{
//		theEvents.genPos();
//		x=theEvents.getx0();
//		y=theEvents.gety0();
//		r=sqrt(x*x+y*y)*1000;
//		hist.Fill(r);
//	}
//
//	TCanvas c2("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,800);
//	hist.Draw();
//	c2.SaveAs(GraphsDir+"test2.pdf");


	/*Make some macros maybe for last time ahahhahhahahhahhahahhahahhaha*/
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID533", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID534", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID535", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID536", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID537", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID538", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 22.4 nm"); //Geant4.9.6.p03 fixed beam ith 22.4 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID539", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 22.4 nm"); //Geant4.9.6.p03 fixed beam ith 22.4 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID540", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 22.4 nm"); //Geant4.9.6.p03 fixed beam ith 22.4 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID541", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 22.4 nm"); //Geant4.9.6.p03 fixed beam ith 22.4 nm dead layer
//	makeGeantMacroFilesForClusterEventGenerator(0,1, "RID542", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 22.4 nm"); //Geant4.9.6.p03 fixed beam ith 22.4 nm dead layer

	/* check no can/gold */
//	const int numFiles=1;
//	int deadLayers[numFiles]={0};
//	double means[numFiles];
//	double fitMeans[numFiles];
//	double counts[numFiles];
//	TString geantVer="Geant4.9.6.p03b";
//	for(int i =0;i<numFiles;i++)
//	{
//		TString dervFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_CanOff_"+geantVer+"_Derv.root";
//
//		TString resultsFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_CanOff_"+geantVer+".root";
//		TFile resultsFile(resultsFilePath);
//		TTree* resultsTree=(TTree*) resultsFile.Get("g");
//		createDerivedResultsFile(0, resultsTree,dervFilePath,100000);
//
//		TFile dervFile(dervFilePath);
//		TTree* dervTree=(TTree*) dervFile.Get("rd");
//		dervTree->Draw("SBDEDepBlur>>mcHist","SBDEDepBlur > 0");
//		TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//		means[i]=mcHist->GetMean();
//		counts[i]=mcHist->Integral(1,mcHist->GetNbinsX());
//		TF1* func = new TF1("f1","gaus",17,23);
//		func->SetParameters(500,mcHist->GetMean(),mcHist->GetRMS());
//		mcHist->Fit("f1","V");
//		fitMeans[i]=func->GetParameter(1);
//	}
//
//	for(int i =0;i<numFiles;i++)
//	{
//		cout << deadLayers[i] << "\t" << counts[i] << "\t" << means[i] << "\t" << fitMeans[i] << endl;
//	}
//
//	const int numFiles=1;
//		int deadLayers[numFiles]={20};
//		double means[numFiles];
//		double fitMeans[numFiles];
//		double counts[numFiles];
//		TString geantVer="Geant4.9.6.p02Expanded";
//		for(int i =0;i<numFiles;i++)
//		{
//			TString dervFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_"+geantVer+"_Derv.root";
//
//			TString resultsFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_"+geantVer+".root";
//			TFile resultsFile(resultsFilePath);
//			TTree* resultsTree=(TTree*) resultsFile.Get("g");
//			createDerivedResultsFile(0, resultsTree,dervFilePath,100000);
//
//			TFile dervFile(dervFilePath);
//			TTree* dervTree=(TTree*) dervFile.Get("rd");
//			dervTree->Draw("SBDEDepBlur>>mcHist","SBDEDepBlur > 0");
//			TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//			means[i]=mcHist->GetMean();
//			counts[i]=mcHist->Integral(1,mcHist->GetNbinsX());
//			TF1* func = new TF1("f1","gaus",17,23);
//			func->SetParameters(500,mcHist->GetMean(),mcHist->GetRMS());
//			mcHist->Fit("f1","V");
//			fitMeans[i]=func->GetParameter(1);
//		}
//
//		for(int i =0;i<numFiles;i++)
//		{
//			cout << deadLayers[i] << "\t" << counts[i] << "\t" << means[i] << "\t" << fitMeans[i] << endl;
//		}


	/* check full production run after beam change */
//		const int numFiles=4;
//		double means[numFiles];
//		double fitMeans[numFiles];
//		double counts[numFiles];
//		TString rids[numFiles]={"RID533","RID535","RID538","RID540"};
//		TString titleStrings[numFiles]={"MC:p02:ThreeBody","MC:p02:FourBody","MC:p03:ThreeBody","MC:p03:FourBody"};
//		TH1D* expHist = getTabSeperatedHist(TString(HISTS_DIR) + "BGO_d239_d240_ep_pEn.txt", EDEPE_DIM2);
//		expHist->SetTitle("Exp;Proton energy (keV);Counts");
//		expHist->GetXaxis()->SetRangeUser(0,32);
//
//		TH1* expHists[1]={expHist};
//
//
//		for(int i =0; i<numFiles; i++)
//		{
//			TString resultsFilePath=ResultsDir+"Results_"+rids[i]+"_P_STD_";
//			TString dervFilePath=ResultsDir+"Results_"+rids[i]+"_P_STD_";
//
//			if(i % 2 ==false)
//			{
//				resultsFilePath+="3B_130429_Set0.root";
//				dervFilePath+="3B_130429_Derv_Set0.root";
//			}
//			else
//			{
//				resultsFilePath+="4B_100eV_Endpoint_130429_Set0.root";
//				dervFilePath+="4B_100eV_Endpoint_130429_Derv_Set0.root";
//			}
//
////			TFile resultsFile(resultsFilePath);
////			TTree* resultsTree=(TTree*) resultsFile.Get("g");
////			createDerivedResultsFile(0, resultsTree,dervFilePath,1000000);
//
//			TFile dervFile(dervFilePath);
//			TTree* dervTree=(TTree*) dervFile.Get("rd");
//			dervTree->Draw("SBDEDepBlur >> mcHist(8000,0,800)","SBDEDepBlur > 13 && SBDEDepBlur < 31");
//			TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//			means[i]=mcHist->GetMean();
//			counts[i]=mcHist->Integral(1,mcHist->GetNbinsX());
//			TF1* func = new TF1("f1","gaus",17,23);
//			func->SetParameters(500,mcHist->GetMean(),mcHist->GetRMS());
//			mcHist->Fit("f1","V");
//			fitMeans[i]=func->GetParameter(1);
//			TH1* mcHists[1]={mcHist};
//			scaleHistogramsTogether(expHist,mcHist);
//			mcHist->SetTitle(titleStrings[i]);
//			expHist->SetName("Exp");
//
//			plotExpVersusMCToImage(1,expHists,1,mcHists,rids[i],"test_"+rids[i]+".png");
//
//			mcHist->Reset();
//		}
//
//		for(int i =0;i<numFiles;i++)
//		{
//			cout << rids[i] << "\t" << counts[i] << "\t" << means[i] << "\t" << fitMeans[i] << endl;
//		}

	/* check no can/gold */
//	const int numFiles=1;
//	int deadLayers[numFiles]={0};
//	double means[numFiles];
//	double fitMeans[numFiles];
//	double counts[numFiles];
//	TString geantVer="Geant4.9.6.p02b";
//	for(int i =0;i<numFiles;i++)
//	{
//		TString dervFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_GoldOff_"+geantVer+"_Derv.root";
//
//		TString resultsFilePath=MRKRunDir+"Results_160202DL"+int2str(deadLayers[i])+"_GoldOff_"+geantVer+".root";
//		TFile resultsFile(resultsFilePath);
//		TTree* resultsTree=(TTree*) resultsFile.Get("g");
//		createDerivedResultsFile(0, resultsTree,dervFilePath,100000);
//
//		TFile dervFile(dervFilePath);
//		TTree* dervTree=(TTree*) dervFile.Get("rd");
//		dervTree->Draw("SBDEDepBlur>>mcHist","SBDEDepBlur > 0");
//		TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//		means[i]=mcHist->GetMean();
//		counts[i]=mcHist->Integral(1,mcHist->GetNbinsX());
//		TF1* func = new TF1("f1","gaus",17,23);
//		func->SetParameters(500,mcHist->GetMean(),mcHist->GetRMS());
//		mcHist->Fit("f1","V");
//		fitMeans[i]=func->GetParameter(1);
//	}
//
//	for(int i =0;i<numFiles;i++)
//	{
//		cout << deadLayers[i] << "\t" << counts[i] << "\t" << means[i] << "\t" << fitMeans[i] << endl;
//	}


/* plot 10.02*/

//	const int numFiles=1;
//	TH1D* expHist = getTabSeperatedHist(TString(HISTS_DIR) + "BGO_d239_d240_ep_pEn.txt", EDEPE_DIM2);
//	expHist->SetTitle("Exp;Proton energy (keV);Counts");
//	expHist->GetXaxis()->SetRangeUser(0,32);
//
//	TH1* expHists[1]={expHist};
//
//
//	TString resultsFilePath=MRKRunDir+"Results_160202DL24_Geant4.10.02.root";
//	TString dervFilePath=MRKRunDir+"Results_160202DL24_Geant4.10.02_Derv.root";
//
//
//	TFile dervFile(dervFilePath);
//	TTree* dervTree=(TTree*) dervFile.Get("rd");
//	dervTree->Draw("SBDEDepBlur >> mcHist(8000,0,800)","SBDEDepBlur > 13 && SBDEDepBlur < 31");
//	TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//	TH1* mcHists[1]={mcHist};
//	scaleHistogramsTogether(expHist,mcHist);
//	mcHist->SetTitle("10.02");
//	expHist->SetName("Exp");
//
//	plotExpVersusMCToImage(1,expHists,1,mcHists,"","test_10.02.png");
//
//	mcHist->Reset();

	/*version analysis runs*/
//	createMCAnalysisAndPlot( "AID87",  "4.9.6.p02 with beam update", 1,  "RID533", "RID534",  1, "RID535",  "RID536", "RID537",  STD_EPCutSetH);
//	createMCAnalysisAndPlot( "AID88",  "4.9.6.p04 with beam update", 1,  "RID538", "RID539",  1, "RID540",  "RID541", "RID542",  STD_EPCutSetH);

	/*Check new vs old analysis 160208*/
//	RDK2MCAnalysis* anals[2];
//	int numFiles=1;
//
//	anals[0]=new RDK2MCAnalysis( "AID87",  "4.9.6.p02 with beam update", 1,  "RID533", "RID534",  1, "RID535",  "RID536", "RID537",  STD_EPCutSetH);
//	anals[1]=new RDK2MCAnalysis("AID88",  "4.9.6.p04 with beam update", 1,  "RID538", "RID539",  1, "RID540",  "RID541", "RID542",  STD_EPCutSetH);
//	double epCount[2],epCountError[2];
//	double ep4Count[2],ep4CountError[2];
//	double epgCount[2],epgCountError[2];
//	double epbgCount[2],epbgCountError[2];
//	double eCount[2],eCountError[2];
//	double pCount[2],pCountError[2];
//	double e4Count[2],e4CountError[2];
//	double p4Count[2],p4CountError[2];
//	double gCount[2],gCountError[2];
//	double bgCount[2],bgCountError[2];
//	for(int i=0;i<2;i++)
//	{
//		eCount[i]=anals[i]->CalcEPer3Decay(eCountError[i]);
//		pCount[i]=anals[i]->CalcPPer3Decay(pCountError[i]);
//		e4Count[i]=anals[i]->CalcEPer4Decay(e4CountError[i]);
//		p4Count[i]=anals[i]->CalcPPer4Decay(p4CountError[i]);
//		gCount[i]=anals[i]->CalcGPer4Decay(DET_G, gCountError[i]);
//		bgCount[i]=anals[i]->CalcGPer4Decay(DET_BG, bgCountError[i]);
//
//		epCount[i]=anals[i]->CalcEPPer3Decay(epCountError[i]);
//		ep4Count[i]=anals[i]->CalcEPPer4Decay(ep4CountError[i]);
//		epgCount[i]=anals[i]->CalcEPGPer4Decay(DET_EPG, epgCountError[i]);
//		epbgCount[i]=anals[i]->CalcEPGPer4Decay(DET_EPBG, epbgCountError[i]);
//	}
//
//	TString titles[2]={"Old","New"};
//	for(int i=0;i<2;i++)
//	{
//		cout << titles[i] + " E count Per 3 Decay:   " +d2str(eCount[i]) +" +/- "+d2str(eCountError[i]) << endl;
//		cout << titles[i] + " P count Per 3 Decay:   " +d2str(pCount[i]) +" +/- "+d2str(pCountError[i])<< endl;
//		cout << titles[i] + " EP count Per 3 Decay:   " +d2str(epCount[i]) +" +/- "+d2str(epCountError[i])<< endl;
//		cout << titles[i] + " E count Per 4 Decay:   " +d2str(e4Count[i]) +" +/- "+d2str(e4CountError[i]) << endl;
//		cout << titles[i] + " P count Per 4 Decay:   " +d2str(p4Count[i]) +" +/- "+d2str(p4CountError[i])<< endl;
//		cout << titles[i] + " BGO G count Per 4 Decay:   " +d2str(gCount[i]) +" +/- "+d2str(gCountError[i])<< endl;
//		cout << titles[i] + " BAPD G count Per 4 Decay:   " +d2str(bgCount[i]) +" +/- "+d2str(bgCountError[i])<< endl;
//		cout << titles[i] + " EP count Per 4 Decay:   " +d2str(ep4Count[i]) +" +/- "+d2str(ep4CountError[i])<< endl;
//		cout << titles[i] + " BGO EPG count Per 4 Decay:   " +d2str(epgCount[i]) +" +/- "+d2str(epgCountError[i])<< endl;
//		cout << titles[i] + " BAPD EPG count Per 4 Decay:   " +d2str(epbgCount[i]) +" +/- "+d2str(epbgCountError[i])<< endl;
//
//
//	}

//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	double bgoEPCounts=2.31112e+08;
//	double apdEPCounts=2.51376e+07;
//	int mcID=88;
//	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID);
//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_"+bgoExpID+"_epg_gEn.txt", bgoCorrection, bgoEPCounts, HistsDir + "bAPD_"+apdExpID+"_epg_gEn.txt",apdCorrection, apdEPCounts, mcID,true,true,false,true,false);


	/*Compare SRIM Calculations */
//	TH1D* expHist = getTabSeperatedHist(TString(HISTS_DIR) + "BGO_d239_d240_ep_pEn.txt", EDEPE_DIM2);
//	expHist->SetTitle("Exp;Proton energy (keV);Counts");
//	expHist->GetXaxis()->SetRangeUser(0,32);
//	TH1* expHists[1]={expHist};
//
//	int numSims=3;
//	TH1* theHists[numSims];
//
//	TString titles[numSims]={"SRIM","4.9.6.p02","4.9.6.p03"};
//
//	for(int i=1; i<numSims;i++)
//	{
//		TFile directFile(MRKGEANTDATA_DIR+TString("SRIMComparison"+titles[i]+".root"));
//		TTree* resultsTree=(TTree*) directFile.Get("g");
//		gROOT->cd();
//		theHists[i]=new TH1D("geantHist"+titles[i],"Geant"+titles[i]+";Energy(keV);Normalized counts",250,0,25);
//		resultsTree->Draw("SBDEDepTotal >> geantHist"+titles[i],"SBDEDepTotal > 0");
//	//	resultsTree->Draw("SBDKEFirst >> geantHist","SBDKEFirst > 0");
//		theHists[i]->Sumw2();
//
//		directFile.Close();
//		cout << theHists[i]->Integral() << endl;
////		scaleHistTo(theHists[i],1e6);
//	}
//	TH1D* srimHist=getTabSeperatedHist(HistsDir+"Hist_SRIM_ProtonAt0Deg_NoDead_Jeff140414.txt",250,0,25);
//	srimHist->SetTitle("SRIM");
//	theHists[0]=srimHist;
//
//
//
//	cout << srimHist->Integral() << endl;
//
//
//
//	plotExpVersusMCToImage(0,NULL,numSims,theHists,"25 keV protons at 0 degrees through 20 nm AU and no deadlayer","Geant4vsSRIM_0degProton20nmGoldNoDeadLayer.png");

	/*Beam bug post mortem*/


//		TChain* currentChain=makeEventChain(MRKRunDir+"Events_STD_3B_130429_Set0.root","t",1);
//		TChain* buggedChain=makeEventChain(EventsDir+"Events_STD_3B_130429_Set0.root","t",1);
//		TFile currentDervFile(MRKRunDir+"Results_160202DL20_Geant4.9.6.p02b_Derv.root");
//		TTree* currentDervTree=(TTree*) currentDervFile.Get("rd");
//		TFile buggedDervFile(MRKRunDir+"Results_160202DL20_Geant4.9.6.p02c_Derv.root");
//		TTree* buggedDervTree=(TTree*) buggedDervFile.Get("rd");
//
//		currentChain->AddFriend(currentDervTree,"rd");
//		buggedChain->AddFriend(buggedDervTree,"rd");
//		gROOT->cd();
//
//		TH1* buggedHist =new TH1D("buggedHist","Bugged;Radius (mm);Counts (AU)",130,0,13);
//		buggedChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> buggedHist","Entry$ < 100000");
//		buggedHist->Sumw2();
//
//		TH1* buggedHitsHist =new TH1D("buggedHitsHist","Bugged Hits;Radius (mm);Counts (AU)",130,0,13);
//		buggedChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> buggedHitsHist","rd.SBDEDepBlur > 10 && Entry$ < 100000");
//		buggedHitsHist->Sumw2();
//
//		TH1* currentHist =new TH1D("currentHist","Fixed;Radius (mm);Counts (AU)",130,0,13);
//		currentChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> currentHist","Entry$ < 100000");
//		currentHist->Sumw2();
//
//		TH1* currentHitsHist =new TH1D("currentHitsHist","Fixed Hits;Radius (mm);Counts (AU)",130,0,13);
//		currentChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> currentHitsHist","rd.SBDEDepBlur>10 && Entry$ < 100000");
//		currentHitsHist->Sumw2();
//
//
//
//
//
//
//		TH1* theHists[4]={currentHist,buggedHist,currentHitsHist,buggedHitsHist};
//
//		for(int i=0;i<4;i++)
//		{
//			scaleHistTo(theHists[i],1);
//		}
//
//		plotExpVersusMCToImage(0,NULL,4,theHists,"","BuggedBeamTest.pdf");
//
//		TCanvas c2("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,800);
//		buggedChain->Draw("ep0:SBDEDepBlur","SBDEDepBlur>0","lego0");
//		c2.SaveAs(GraphsDir+"test.pdf");

	/*Expanded beam*/
//	TChain* currentChain=makeEventChain(MRKRunDir+"Events_STD_3B_130429_Set0.root","t",1);
//	TChain* expandedChain=makeEventChain(MRKRunDir+"Events_STD_3B_Expanded_130429_Set0.root","t",1);
//	TFile currentDervFile(MRKRunDir+"Results_160202DL20_Geant4.9.6.p02b.root");
//	TTree* currentDervTree=(TTree*) currentDervFile.Get("g");
//	TFile expandedDervFile(MRKRunDir+"Results_160202DL20_Geant4.9.6.p02Expanded.root");
//	TTree* expandedDervTree=(TTree*) expandedDervFile.Get("g");
//
//	currentChain->AddFriend(currentDervTree,"g");
//	expandedChain->AddFriend(expandedDervTree,"g");
//	gROOT->cd();
//
//	TH1* expandedHist =new TH1D("expandedHist","Expanded;Radius (mm);Counts (AU)",100,0,50);
//	expandedChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> expandedHist","Entry$ < 100000");
//	expandedHist->Sumw2();
//
//	TH1* expandedHitsHist =new TH1D("expandedHitsHist","Expanded Hits;Radius (mm);Counts (AU)",100,0,50);
//	expandedChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> expandedHitsHist","g.SBDEDepTotal > 0 && Entry$ < 100000");
//	expandedHitsHist->Sumw2();
//
//	TH1* currentHist =new TH1D("currentHist","Fixed;Radius (mm);Counts (AU)",100,0,50);
//	currentChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> currentHist","Entry$ < 100000");
//	currentHist->Sumw2();
//
//	TH1* currentHitsHist =new TH1D("currentHitsHist","Fixed Hits;Radius (mm);Counts (AU)",100,0,50);
//	currentChain->Draw("sqrt(x0*x0+y0*y0)*1000 >> currentHitsHist","g.SBDEDepTotal>0 && Entry$ < 100000");
//	currentHitsHist->Sumw2();






//	TH1* theHists[4]={currentHist,expandedHist,currentHitsHist,expandedHitsHist};
//
//	for(int i=0;i<4;i++)
//	{
//		if(i%2==0)
//			scaleHistTo(theHists[i],1);
//		else
//			scaleHistTo(theHists[i],3);
//	}
//	scaleHistogramsTogether(theHists[1],theHists[3],1,10);
//	plotExpVersusMCToImage(0,NULL,4,theHists,"","expandedBeamTest.pdf");


	/*Final run part 5? 6?*/
//	createMCAnalysisAndPlot( "AID87",  "4.9.6.p02 with beam update", 48,  "RID533", "RID534",  48, "RID535",  "RID536", "RID537",  STD_EPCutSetL);
//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	double bgoEPCounts=2.31112e+08;
//	double apdEPCounts=2.51376e+07;
//	int mcID=87;
//	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID);
//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_"+bgoExpID+"_epg_gEn.txt", bgoCorrection, bgoEPCounts, HistsDir + "bAPD_"+apdExpID+"_epg_gEn.txt",apdCorrection, apdEPCounts, mcID,true,true,false,true,false);

//	createMCAnalysisAndPlot( "AID89",  "Final?? (fixed BAPD cuts)", 48,  "RID518", "RID519",  48, "RID520",  "RID521", "RID522",  STD_EPCutSetL);
//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	double bgoEPCounts=2.31112e+08;
//	double apdEPCounts=2.51376e+07;
//	int mcID=89;
//	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID);
//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_"+bgoExpID+"_epg_gEn.txt", bgoCorrection, bgoEPCounts, HistsDir + "bAPD_"+apdExpID+"_epg_gEn.txt",apdCorrection, apdEPCounts, mcID,true,true,false,true,false);

 /*Check all files*/
//	int numFiles=48;
//
//
//	for(int i =0;i<numFiles;i++)
//	{
//		TString dervFilePath=ResultsDir+"Results_RID537_G_STD_4B_100eV_Endpoint_130429_Derv_Set"+int2str(i)+".root";
//
//		TFile dervFile(dervFilePath);
//		TTree* dervTree=(TTree*) dervFile.Get("rd");
//		gROOT->cd();
//		dervTree->Draw("BGO1EDepBlurLOAllExpParam+BGO2EDepBlurLOAllExpParam+BGO3EDepBlurLOAllExpParam+BGO4EDepBlurLOAllExpParam+BGO5EDepBlurLOAllExpParam+BGO6EDepBlurLOAllExpParam+BGO7EDepBlurLOAllExpParam+BGO8EDepBlurLOAllExpParam+BGO9EDepBlurLOAllExpParam+BGO10EDepBlurLOAllExpParam+BGO11EDepBlurLOAllExpParam+BGO12EDepBlurLOAllExpParam>>mcHist","SBDEDepBlur > 0");
//		TH1F *mcHist = (TH1F*)gDirectory->Get("mcHist");
//
//		cout << mcHist->Integral() << endl;
//	}


//	RDK2MCAnalysis* anals[1];
//	int numFiles=1;
//
//	anals[0]=new RDK2MCAnalysis( "AID87",  "4.9.6.p02 with beam update", 48,  "RID533", "RID534",  48, "RID535",  "RID536", "RID537",  STD_EPCutSetH);
//	double epCount[2],epCountError[2];
//	double ep4Count[2],ep4CountError[2];
//	double epgCount[2],epgCountError[2];
//	double epbgCount[2],epbgCountError[2];
//	double eCount[2],eCountError[2];
//	double pCount[2],pCountError[2];
//	double e4Count[2],e4CountError[2];
//	double p4Count[2],p4CountError[2];
//	double gCount[2],gCountError[2];
//	double bgCount[2],bgCountError[2];
//	for(int i=0;i<1;i++)
//	{
//		eCount[i]=anals[i]->CalcEPer3Decay(eCountError[i]);
//		pCount[i]=anals[i]->CalcPPer3Decay(pCountError[i]);
//		e4Count[i]=anals[i]->CalcEPer4Decay(e4CountError[i]);
//		p4Count[i]=anals[i]->CalcPPer4Decay(p4CountError[i]);
//		gCount[i]=anals[i]->CalcGPer4Decay(DET_G, gCountError[i]);
//		bgCount[i]=anals[i]->CalcGPer4Decay(DET_BG, bgCountError[i]);
//
//		epCount[i]=anals[i]->CalcEPPer3Decay(epCountError[i]);
//		ep4Count[i]=anals[i]->CalcEPPer4Decay(ep4CountError[i]);
//		epgCount[i]=anals[i]->CalcEPGPer4Decay(DET_EPG, epgCountError[i]);
//		epbgCount[i]=anals[i]->CalcEPGPer4Decay(DET_EPBG, epbgCountError[i]);
//	}
//
//	TString titles[1]={"New"};
//	for(int i=0;i<1;i++)
//	{
//		cout << titles[i] + " E count Per 3 Decay:   " +d2str(eCount[i]) +" +/- "+d2str(eCountError[i]) << endl;
//		cout << titles[i] + " P count Per 3 Decay:   " +d2str(pCount[i]) +" +/- "+d2str(pCountError[i])<< endl;
//		cout << titles[i] + " EP count Per 3 Decay:   " +d2str(epCount[i]) +" +/- "+d2str(epCountError[i])<< endl;
//		cout << titles[i] + " E count Per 4 Decay:   " +d2str(e4Count[i]) +" +/- "+d2str(e4CountError[i]) << endl;
//		cout << titles[i] + " P count Per 4 Decay:   " +d2str(p4Count[i]) +" +/- "+d2str(p4CountError[i])<< endl;
//		cout << titles[i] + " BGO G count Per 4 Decay:   " +d2str(gCount[i]) +" +/- "+d2str(gCountError[i])<< endl;
//		cout << titles[i] + " BAPD G count Per 4 Decay:   " +d2str(bgCount[i]) +" +/- "+d2str(bgCountError[i])<< endl;
//		cout << titles[i] + " EP count Per 4 Decay:   " +d2str(ep4Count[i]) +" +/- "+d2str(ep4CountError[i])<< endl;
//		cout << titles[i] + " BGO EPG count Per 4 Decay:   " +d2str(epgCount[i]) +" +/- "+d2str(epgCountError[i])<< endl;
//		cout << titles[i] + " BAPD EPG count Per 4 Decay:   " +d2str(epbgCount[i]) +" +/- "+d2str(epbgCountError[i])<< endl;
//
//	}

//	makeEGTimingPlot(HistsDir+"ExpJeff_epg_egTOF_151122.txt",HistsDir+"ExpJeff_epbg_ebgTOF_151122.txt");
//	makeZ0Hist();

//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	double bgoEPCounts=2.31112e+08;
//	double apdEPCounts=2.51376e+07;
//	int mcID=87;
//	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID, true);

	/* Make beam shift macros */
//	makeGeantMacroFilesForClusterEventGenerator(0,10, "RID543", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: -1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,10, "RID544", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: -1mm x beam shift
//	makeGeantMacroFilesForhClusterEventGenerator(0,48, "RID545", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: -1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID546", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: -1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID547", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: -1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,10, "RID548", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: +1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,10, "RID549", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: +1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID550", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: +1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID551", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: +1mm x beam shift
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID552", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: +1mm x beam shift

	/* Make homogenous event set  and macros*/
//	createHomogenousEventSettings(500,"160323");
//	makeGeantMacroFilesForClusterEventGenerator(0,100, "RID553", PROTON, "Hmg_3B_160323", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: Homegenous
//	makeGeantMacroFilesForClusterEventGenerator(0,100, "RID554", ELECTRON, "Hmg_3B_160323", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSiDeadLayerLengthCmd 13.7 nm"); //Geant4.9.6.p02 fixed beam with 13.7 nm dead layer: Homegenous

	/* Testing homogenous events */
//	createHomogenousEventsUsingSettings(1, 1000000, "160323");
//
	/*Analyze and plot beam shifts again after beam radii problem fix 160331*/
//	createMCAnalysisAndPlot( "AID90",  "Beam shift 1mm x-", 10,  "RID543", "RID544",  48, "RID545",  "RID546", "RID547",  STD_EPCutSetL);
//	createMCAnalysisAndPlot( "AID91",  "Beam shift 1mm x+", 10,  "RID548", "RID549",  48, "RID550",  "RID551", "RID552",  STD_EPCutSetL);
//
//	RDK2AnalysisPlotter thePlotter("BAPD Field from COMSOL","","AID87 AID90 AID91"); thePlotter.makeAll();

	/*Create Paper Hists */
//	double bgoCorrection=-.014;
//	double apdCorrection=-.003;
//	TString bgoExpID="d239_d240";
//	TString apdExpID="d236_d237";
//	double bgoEPCounts=2.31112e+08;
//	double apdEPCounts=2.51376e+07;
//	int mcID=89;
////	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID);
//	makeCombinedBGOAndBAPDPlot(HistsDir+"BGO_"+bgoExpID+"_epg_gEn.txt", bgoCorrection, bgoEPCounts, HistsDir + "bAPD_"+apdExpID+"_epg_gEn.txt",apdCorrection, apdEPCounts, mcID,true,true,false,true,false);



//	makeGeantMacroFilesForClusterEventGeneratorMakeEventsOnly(0,100, "Hmg_3B_160323");

	/* Analyze hmg runs */
//	createLilBTreesFromHmgEventFiles(100, "Hmg_3B_160323", 21, -1, 0.1);
//	makeDWCutHists(21, "Hmg02", "First try Hmg Geant4.9.6.p02", 99,  "RID553", "RID554",  0, "",  "", "",  STD_EPCutSetL);
//	plotLilBExpFitToMC(TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg02", EDEPE_DIM2, 21, -1, 0.1, 100., 700.);
//	makeEE0Plots(21, "Hmg02", "First try Hmg Geant4.9.6.p02", 99,  "RID553", "RID554",  0, "",  "", "",  STD_EPCutSetL);

	/*Lifetime backscattering*/
	//	analyzeProtonBackscattering("LifetimeProtonBackScattering");
//	makePosXYPlotFromReflect("LifetimeProtonBackScattering");

	/*Provide Kevin Graph of Parameterized LO Model*/
//	outputParameterizedLOModelToText();
//	outputBGODataToFile();

	/*Lifetime backscattering Version differences*/
//	analyzeProtonBackscattering("LifetimeProtonBackScattering_9.6.p02");
//	analyzeProtonBackscattering("LifetimeProtonBackScattering_10.02");
//	makePosXYPlotFromReflect("LifetimeProtonBackScattering_9.6.p02","Protons detected entering SBD Active Layer (Geant4.9.6.p02 with deadlayer/gold RDK2 Geo)");
//	makePosXYPlotFromReflect("LifetimeProtonBackScattering_10.02","Protons detected entering SBD Active Layer (Geant4.10.02 with deadlayer/gold RDK2 Geo)");
//	analyzeProtonBackscattering("LifetimeProtonBackScattering_9.6.p03");
//	analyzeProtonBackscattering("LifetimeProtonBackScattering_9.6.p04");
//	makePosXYPlotFromReflect("LifetimeProtonBackScattering_9.6.p03","Protons detected entering SBD Active Layer (Geant4.9.6.p03 with deadlayer/gold RDK2 Geo)");
//	makePosXYPlotFromReflect("LifetimeProtonBackScattering_9.6.p04","Protons detected entering SBD Active Layer (Geant4.9.6.p04 with deadlayer/gold RDK2 Geo)");

	/*Electron backscattering */
//	analyzeAllEBackscatteringSBDOnly();
//	makeBSPlotsSBDOnly();
//	makeBSPlotsSBDOnly2();

	/*Electron backscattering through secondary events*/
//	makeAllIncidentEventsFromResults(100, "RID554");
//	vector<TString> rids={"Liv","GS","WVI","SingleScatter","Penelope"};
//	for(unsigned int i=0;i<rids.size();i++)
//	{kwatch
//		makeSecondaryMacros(100, "RID554", "RID"+TString::Itoa(i+555,10),rids[i]);
//	}
//	makeDWCutHists(21, "Hmg03", "Geant4.10.02 Livermore Secondary Events", 99,  "RID553", "RID555",  0, "",  "", "",  STD_EPCutSetL);
//	makeDWCutHists(21, "Hmg04", "Geant4.10.02 WVI Secondary Events", 99,  "RID553", "RID557",  0, "",  "", "",  STD_EPCutSetL);
//	plotLilBExpFitToMC("Geant4.9.6.p02 Livermore (small range) Added constant error",TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg01", EDEPE_DIM2, 21, -0.1, 0.01, 100., 700.);
//	plotLilBExpFitToMC("Geant4.9.6.p02 Livermore (large range) Added constant error",TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg02", EDEPE_DIM2, 21, -1, 0.1, 100., 700.);
//	plotLilBExpFitToMC("Geant4.10.02 Livermore Secondary Added constant error",TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg03", EDEPE_DIM2, 21, -1, 0.1, 100., 700.);
//	plotLilBExpFitToMC("Geant4.10.02 WVI Secondary Added constant error",TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg04", EDEPE_DIM2, 21, -1, 0.1, 100., 700.);
//	createLilBTreesFromHmgEventFiles(100, "Hmg_3B_160323", 41, -0.2, 0.01,1);
//    makeDWCutHists(41, "Hmg05", "Geant4.9.6.p02 Livermore (medium range)", 99,  "RID553", "RID554",  0, "",  "", "",  STD_EPCutSetL);
//	makeDWCutHists(41, "Hmg06", "Geant4.10.02 WVI Secondary Events (medium range)", 99,  "RID553", "RID557",  0, "",  "", "",  STD_EPCutSetL);
//	makeDWCutHists(41, "Hmg07", "Geant4.10.02 Livermore Secondary Events (medium range)", 99,  "RID553", "RID555",  0, "",  "", "",  STD_EPCutSetL);
//	plotLilBExpFitToMC("Geant4.9.6.p02 Livermore (medium range)",TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg05", EDEPE_DIM2, 41, -.2, 0.01, 100., 700.);
//	plotLilBExpFitToMC("Geant4.10.02 WVI Secondary Added constant error",TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg06", EDEPE_DIM2, 41, -.2, 0.01, 100., 700.);
//	plotLilBExpFitToMC("Geant4.10.02 Livermore Secondary Added constant error",TStringkwatch(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", "Hmg07", EDEPE_DIM2, 41, -.2, 0.01, 100., 700.);

	/*Analyze little b through secondary events*/
	RDK2MCAnalyzer theMCAnalyzer("Hmg10", "Geant4.10.02 Livermore Secondary Events", 1,  "RID553", "RID555",  0, "",  "", "",  STD_EPCutSetL,true);
//	RDK2MCAnalyzer theMCAnalyzer("Hmg11", "Geant4.10.02 Livermore Secondary Events", 1,  "RID533", "RID534",  0, "",  "", "",  STD_EPCutSetL,false);
	makeBootStrapLittleB( TString(HISTS_DIR)+"BGO_d239_d240_ep_eEn.txt", EDEPE_DIM2, theMCAnalyzer,   3, -0.1, 0.1, 100., 700.);

	/*Checking hists*/
	TH1* hists[1]={getTabSeperatedHist(HistsDir+"Hist_MC_Hmg10_dwcut1_EP_EDepEProb.txt",EDEPE_DIM2)};
	hists[0]->Rebin(100);
	plotExpVersusMCToImage(1,hists,0,nullptr,"test",GraphsDir+"test.png",false);

//	RDK2MCAnalyzer theMCAnalyzer("Hmg10", "Geant4.10.02 Livermore Secondary Events", 1,  "RID553", "RID555",  0, "",  "", "",  STD_EPCutSetL,true);
//	theMCAnalyzer.MakeAllStandardHists();

	return 0;  //end int main()
}
