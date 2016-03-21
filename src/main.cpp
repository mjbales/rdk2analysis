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

//#include "littlebanalysis.h"

//#include "cParticle.h"
//#include "MCanalysis2.h"
//#include "SBD_Cal.h"
//#include "mattexpanalysis.h"
//#include "RDK2Vars.h"
//#include "RDK2Set.h"

//#include "cintProblemFunctions.h"
//#include "constants.h"
//#include "FermiRCooper.h"
//#include "RDK2MCAnalysis.h"
//#include "RDK2ExpAnalysis.h"
//#include "RDK2Deconv.h"
//#include "cEvents.h"
//#include "cMRKText.h"
//#include "mattmisc.h"
//#include "mattio.h"
//#include "mattanalysis.h"
//#include "mattgraphics.h"
//#include "ROOTField.h"

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

///Coplot weird amplitude issues
//    createBGSubCalibrationSourceHistogramsForBa133ToFiles("C461",CALPEAKTIMINGCUT_ONPEAK);
//    createBGSubCalibrationSourceHistogramsForBa133ToFiles("C461",CALPEAKTIMINGCUT_ONPEAKEXT);
//    analyzeCs137Ba133Peaks("C504", "DoubleSource_Center_130310", true,CALPEAKTIMINGCUT_ONPEAK);
//    const int numHist =3;
//    TH1* expHists[numHist];
//    TH1* mcHists[1];
//    expHists[0]=getTabSeperatedHist(TString(HISTS_DIR)+"CalBGSubHist_C504_Det11_MJB.txt",BGOEXPCAL_DIM);
//    expHists[1]=getTabSeperatedHist(TString(HISTS_DIR)+"CalBGSubHist_C504_Det12_MJB.txt",BGOEXPCAL_DIM);
//    expHists[2]=getTabSeperatedHist(TString(HISTS_DIR)+"CalBGSubHist_C504OnPeakExt_Det12_MJB.txt",BGOEXPCAL_DIM);
//    expHists[3]=getTabSeperatedHist(TString(HISTS_DIR)+"CalBGSubHist_C461_Det11_MJB.txt",BGOEXPCAL_DIM);
//    expHists[4]=getTabSeperatedHist(TString(HISTS_DIR)+"CalBGSubHist_C461OnPeakExt_Det11_MJB.txt",BGOEXPCAL_DIM);
//    mcHists[0]=getTabSeperatedHist(TString(HISTS_DIR)+"BGO_Det11_Cs137Ba133_DoubleSource_Center_130310.txt",BGOCAL_DIM);
//
//    mcHists[0]->GetXaxis()->Set(10000,0.,28700.);
//    mcHists[0]->Rebin(20);
//    mcHists[0]->Scale(.25);
//    mcHists[0]->GetYaxis()->SetRangeUser(0,3400);
//    mcHists[0]->GetXaxis()->SetRangeUser(100,22000);
//
//    expHists[0]->SetTitle("C504_Det11_OnPeak;Energy Channel (corrected to 19000=662);Counts (AU)");
//    expHists[1]->SetTitle("C504_Det12_OnPeak;Energy Channel (corrected to 19000=662);Counts (AU)");
//    expHists[2]->SetTitle("C504_Det12_OnPeakExt;Energy Channel (corrected to 19000=662);Counts (AU)");
//    expHists[3]->SetTitle("C461_Det11_OnPeak;Energy Channel (corrected to 19000=662);Counts (AU)");
//    expHists[4]->SetTitle("C461_Det11_OnPeakExt;Energy Channel (corrected to 19000=662);Counts (AU)");
//    mcHists[0]->SetTitle("MC_WKhod_Cs137Ba133_Det11;Energy Channel (corrected to 19000=662);Counts (AU)");
//
//    for (int i = 0;i< numHist;i++)
//    {
//    	expHists[i]->Rebin(4);
//    	expHists[i]->GetYaxis()->SetRangeUser(0,3400);
//    	expHists[i]->GetXaxis()->SetRangeUser(100,22000);
//    }
//
//
//    plotExpVersusMCToImage(numHist,expHists,1,mcHists,"Relative amplitude irregularities and timing in Cs137 + Ba-133","CalAmpIrregularities_130404.png");

///Creating timing with OnPeakExt
//    const int numDetectors=12;
//    TString seriesList="C439";
//    TString triggerList="11 11 5 5 5 1 5 1";
//    int numSeries= numItemsInStringList(seriesList);
//    TString seriesString,triggerString;
//    TH1* expHists[numDetectors-1];
//    TString filePath;
//    HistDim chnDim={200,0,1500};
//    for (int j = 0;j< numSeries;j++)
//    {
//
//        stringFromList(seriesList,j,seriesString);
//        stringFromList(triggerList,j,triggerString);
//        createRawChnCalibrationSourceHistogramsToFiles(seriesString,chnDim);
//
//        double yMax=1;
//        double yMin=1;
//        gROOT->cd();
//        const int colorList[12]= {kBlack, kBlue, kRed, kGreen,kYellow , kCyan,kMagenta, kGray, kPink, kOrange, kViolet,kAzure};
//        const int markerList[14]={21,20,22,23,25,24,32,33,3,5,28,26,30,29};
//        gStyle->SetTitleColor(kWhite);
//        TCanvas* c2 = new TCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,800);
//        c2->SetFillColor(kGray);
//        c2->cd();
//        c2->SetLeftMargin(.12);
//        c2->SetRightMargin(.06);
//        c2->SetLogy();
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
//        TLegend *theLegend;
//        double legendSize=.025*(numDetectors-1);
//        if(legendSize < .15)
//            legendSize=.15;
//        theLegend = new TLegend(0.68,.9-legendSize,0.93,0.9);
//
//        for (int i = 0;i< numDetectors-1;i++)
//        {
//            int currentDetector=i+1;
//            if(currentDetector>6)
//            {
//                currentDetector=i+2;
//            }
//            filePath=HistsDir;
//            filePath+="CalRawPkChnHist_"+seriesString+"_Det"+int2str(currentDetector)+"_MJB.txt";
//            expHists[i]=getTabSeperatedHist(filePath,chnDim);
//            TString histNameString=seriesString+"Det"+int2str(currentDetector)+"Hist";
//            expHists[i]->SetName(histNameString);
//            TString titleString="Det"+int2str(currentDetector)+";PkChn (channels);counts(AU)";
//            expHists[i]->SetTitle(titleString);
//            if(expHists[i]->GetMaximum() > yMax)
//                yMax=expHists[i]->GetMaximum();
//            expHists[i]->SetLineColor(colorList[(i)%12]);
//            expHists[i]->SetMarkerColor(colorList[(i)%12]);
//            expHists[i]->SetMarkerStyle(markerList[(i)%14]);
//            expHists[i]->SetLineWidth(2);
//            expHists[i]->SetMarkerSize(1);
//            expHists[i]->GetYaxis()->SetTitleOffset(1.8);
//            theLegend->AddEntry(expHists[i],expHists[i]->GetTitle(),"PE");
//        }
//        expHists[0]->GetYaxis()->SetRangeUser(yMin,yMax);
//        TString imageName="CalRawPkChnHist_OnPeakExt_"+seriesString+".png";
//        TString imageTitle="Timing for " + seriesString;
//        if(seriesString[0]=='S')
//            imageTitle+= " Trigger on Channel "+triggerString;
//
//
//
//        expHists[0]->SetTitle("");
//        expHists[0]->Draw("P E1 X0");
//        for (int i = 1;i< numDetectors-1;i++)
//        {
//            expHists[i]->Draw("P E1 X0 same");
//        }
//
//        c2->Update();
//
//        TString funcString=d2str(yMax);
//        TF1 onPeakRegion("onPeakRegion",funcString,0.,800.);
//        TF1 offPeakRegion("offPeakRegion",funcString,800.,1300.);
//
//        onPeakRegion.SetFillStyle(3013);
//        onPeakRegion.SetFillColor(kGreen);
//        onPeakRegion.SetLineColor(kBlack);
//        onPeakRegion.Draw("f same");
//
//        offPeakRegion.SetFillStyle(3013);
//        offPeakRegion.SetFillColor(kRed);
//        offPeakRegion.SetLineColor(kBlack);
//        offPeakRegion.Draw("f same");
//
//        theLegend->Draw();
//        theLegend->SetBorderSize(1);
//        theLegend->SetFillColor(kWhite);
//
//        TPaveText* textBox = new TPaveText(.05,.925,.95,.995,"NDC");
//        textBox->SetTextSize(.03);
//        textBox->AddText(imageTitle);
//        textBox->SetFillColor(kWhite);
//        textBox->SetLineColor(kWhite);
//
//        textBox->SetShadowColor(kWhite);
//        textBox->SetFillStyle(4000);//transparent
//        textBox->SetLineWidth(0);
//        textBox->Draw();
//
//        TString fullPath=GRAPHS_DIR;
//        fullPath+=imageName;
//        c2->SaveAs(fullPath);
//
//        delete c2;
//        for (int i = 0;i< numDetectors-1;i++)
//        {
//        	delete expHists[i];
//        }
//
//
//    }

///Testing new calibration function

//    analyzeLOofCalPeaks("C504", "Center_130310", true,CALPEAKTIMINGCUT_ONPEAK);
//   createHistogramsForMCCalibraitonSourceData(10,"Ba133_Center_130310");
//    createBGSubCalibrationSourceHistogramsForBa133ToFiles("C461",CALPEAKTIMINGCUT_ONPEAK);
//    analyzeLOofCalPeaks("C461", "Center_130310", false,CALPEAKTIMINGCUT_ONPEAKEXT);
	//  createBGSubCalibrationSourceHistogramsForBa133ToFiles("C461",CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C461", "Offset2_130310", false,CALPEAKTIMINGCUT_OFFPEAK);

//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Center2_130408");
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Center_130310");
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Offset2_130310");
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_DoubleSource_Center_130310");
//    TH1* expHists[4];
//    expHists[0]=createCalMCHistWithRes(1, "Center_130310", "Ba133",false,1);
//    expHists[1]=createCalMCHistWithRes(1, "Offset2_130310", "Ba133",false,1);
//    expHists[2]=createCalMCHistWithRes(1, "DoubleSource_Center_130310", "Ba133",false,1);
//    expHists[3]=createCalMCHistWithRes(1, "Center2_130408", "Ba133",false,1);
//
//    for (int i = 0;i< 4;i++)
//    {
//        expHists[i]->Rebin(20);
//
//    }
//    expHists[0]->GetYaxis()->SetRangeUser(0,13000);
//    expHists[0]->GetXaxis()->SetRangeUser(0,500);
//
//    plotExpVersusMCToImage(4,expHists,0,NULL,"","test.png");
//    analyzeLOofCalPeaks("C461", "Center2_130408", false,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C461", "Center2_130408", false,CALPEAKTIMINGCUT_ONPEAKEXT);

//    createBGSubCalibrationSourceHistogramsForBa133ToFiles("C439",CALPEAKTIMINGCUT_OFFPEAK);

//    createBGSubCalibrationSourceHistogramsForBa133ToFiles("C439",CALPEAKTIMINGCUT_ONPEAKEXT);

//    analyzeLOofCalPeaks("C439", "Center2_130408", false,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C439", "Center2_130408", false,CALPEAKTIMINGCUT_ONPEAKEXT);

	//To Be Run
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Center2_130408");
//    createHistogramsForMCCalibraitonSourceData(10,"Co57_Center2_130408");
//    createHistogramsForMCCalibraitonSourceData(10,"Cs137_Center2_130408");
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_DoubleSource_Center2_130408");
//    createHistogramsForMCCalibraitonSourceData(10,"Cs137_DoubleSource_Center2_130408");
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Offset3_130408");
//    analyzeLOofCalPeaks("C461", "Center2_130408", false,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C461", "Center2_130408", false,CALPEAKTIMINGCUT_ONPEAKEXT);
//    analyzeLOofCalPeaks("C439", "Center2_130408", false,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C439", "Center2_130408", false,CALPEAKTIMINGCUT_ONPEAKEXT);
//    analyzeLOofCalPeaks("C439", "Center2_130408", false,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C439", "Center2_130408", false,CALPEAKTIMINGCUT_ONPEAKEXT);
//    analyzeLOofCalPeaks("C504", "DoubleSource_Center2_130408", false,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C504", "DoubleSource_Center2_130408", false,CALPEAKTIMINGCUT_ONPEAKEXT);

///Create co-LO plot
//    const int numHists=6;
//    const int colorList[12]= {kBlack, kBlack, kBlack, kBlue,kBlue ,kBlue, kRed, kRed, kRed, kGreen,kGreen,kGreen};
//    const int markerList[14]={21,20,22,21,20,22,32,33,3,5,28,26,30,29};
//    TString legendString[numHists];
//    for (int peakNum=0;peakNum<1;peakNum++)
//    {
//        TH1* loHists[numHists];
//
//
//        loHists[0]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Cs137Ba133_C504OnPeakExt_DoubleSource_Center_130310_Peak1.txt",13,.5,13.5);
//        legendString[0]="C504_OnPeakExt";
//
//        loHists[1]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C461OnPeakExt_Center2_130408_Peak1.txt",13,.5,13.5);
//        legendString[1]="C461_OnPeakExt";
//
//        loHists[2]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C439OnPeakExt_Center2_130408_Peak1.txt",13,.5,13.5);
//        legendString[2]="C439_OnPeakExt";
//
//        loHists[3]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Cs137Ba133_C504OffPeak_DoubleSource_Center_130310_Peak1.txt",13,.5,13.5);
//        legendString[3]="C504_OffPeak";
//
//        loHists[4]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C461OffPeak_Center2_130408_Peak1.txt",13,.5,13.5);
//        legendString[4]="C461_OffPeak";
//
//        loHists[5]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C439OffPeak_Center2_130408_Peak1.txt",13,.5,13.5);
//        legendString[5]="C439_OffPeak";
//
//        TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1600,1000);
//        theCanvas.SetFillColor(kGray);
//        theCanvas.cd();
//        theCanvas.SetLeftMargin(.12);
//        theCanvas.SetRightMargin(.2);
//        theCanvas.SetBottomMargin(.06);
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
//        double legendSize=.05*numHists;
//        TLegend* theLegend = new TLegend(0.8,.9-legendSize,0.99,0.9);
//        theLegend->SetBorderSize(1);
//        theLegend->SetFillColor(kWhite);
//
//        for (int i = 0;i< numHists;i++)
//        {
//            loHists[i]->GetXaxis()->SetLabelOffset(99);
//            loHists[i]->SetLineColor(colorList[i%12]);
//            loHists[i]->SetMarkerColor(colorList[i%12]);
//            loHists[i]->SetMarkerStyle(markerList[i%14]);
//            //loHists[i]->SetMarkerStyle(21);
//            loHists[i]->GetYaxis()->SetTitle("Light Output");
//            loHists[i]->SetMarkerSize(1.2);
//            loHists[i]->GetXaxis()->Set(13,.45+i*.03,13.45+i*.03);
//            theLegend->AddEntry(loHists[i],legendString[i],"PE");
//        }
//
//
//
//        double xAxisTitleYPosition;
//        if(peakNum==0)
//        {
//            loHists[0]->GetYaxis()->SetRangeUser(0.7,.9);
//            xAxisTitleYPosition=.69;
//
//        }
//        else if(peakNum==1)
//        {
//            loHists[0]->GetYaxis()->SetRangeUser(0.9,1.0);
//            xAxisTitleYPosition=.895;
//        }
//        else
//        {
//            loHists[0]->GetYaxis()->SetRangeUser(0.96,1.01);
//            xAxisTitleYPosition=.958;
//        }
//        loHists[0]->GetXaxis()->SetRangeUser(0,16);
//        TString plotTitleString="Comparison of LO Analysis for Peak "+int2str(peakNum+1);
//        loHists[0]->SetTitle(plotTitleString);
//        loHists[0]->GetXaxis()->SetTimeOffset(2.2);
//
//        gStyle->SetErrorX(0);
//        loHists[0]->Draw("P E1");
//
//        for (int i = 1;i< numHists;i++)
//        {
//        	loHists[i]->Draw("P E1 same");
//        }
//
//        theLegend->Draw();
//
//        //Create a text based X axis
//        TText xAxisTText;
//        //xAxisTText.SetTextAngle(90);
//        xAxisTText.SetTextSize(0.03);
//        xAxisTText.SetTextAlign(21);
//        for (int i=0;i<13;i++)
//        {
//            Float_t x = loHists[0]->GetXaxis()->GetBinCenter(i+1);
//            TString axisTitle=int2str(i+1);
//            if(i==12)
//            {
//              axisTitle="Avg.";
//            }
//            xAxisTText.DrawText(x,xAxisTitleYPosition,axisTitle);
//        }
//
//        TString graphPathString=GraphsDir+"LOAnalysisComparisonPeak"+int2str(peakNum+1)+"130411.png";
//        theCanvas.Draw();
//
//        theCanvas.SaveAs(graphPathString);
//
//        for (int i = 0;i< numHists;i++)
//        {
//        	delete loHists[i];
//        }
//
//
//    }

///Creating timing with OnPeakExt
//    const int numDetectors=12;
//    TString seriesList="C504";
//    TString triggerList="11 11 5 5 5 1 5 1";
//    int numSeries= numItemsInStringList(seriesList);
//    TString seriesString,triggerString;
//    TH1* expHists[numDetectors-1];
//    TString filePath;
//    HistDim chnDim={200,0,1500};
//    for (int j = 0;j< numSeries;j++)
//    {
//
//        stringFromList(seriesList,j,seriesString);
//        stringFromList(triggerList,j,triggerString);
//        createRawChnCalibrationSourceHistogramsToFiles(seriesString,chnDim);
//
//        double yMax=1;
//        double yMin=1;
//        gROOT->cd();
//        const int colorList[12]= {kBlack, kBlue, kRed, kGreen,kYellow , kCyan,kMagenta, kGray, kPink, kOrange, kViolet,kAzure};
//        const int markerList[14]={21,20,22,23,25,24,32,33,3,5,28,26,30,29};
//        gStyle->SetTitleColor(kWhite);
//        TCanvas* c2 = new TCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,800);
//        c2->SetFillColor(kGray);
//        c2->cd();
//        c2->SetLeftMargin(.12);
//        c2->SetRightMargin(.06);
//        c2->SetLogy();
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
//        TLegend *theLegend;
//        double legendSize=.025*(numDetectors-1);
//        if(legendSize < .15)
//            legendSize=.15;
//        theLegend = new TLegend(0.68,.9-legendSize,0.93,0.9);
//
//        for (int i = 0;i< numDetectors-1;i++)
//        {
//            int currentDetector=i+1;
//            if(currentDetector>6)
//            {
//                currentDetector=i+2;
//            }
//            filePath=HistsDir;
//            filePath+="CalRawPkChnHist_"+seriesString+"_Det"+int2str(currentDetector)+"_MJB_NoCuts.txt";
//            expHists[i]=getTabSeperatedHist(filePath,chnDim);
//            TString histNameString=seriesString+"Det"+int2str(currentDetector)+"Hist";
//            expHists[i]->SetName(histNameString);
//            TString titleString="Det"+int2str(currentDetector)+";PkChn (channels);counts(AU)";
//            expHists[i]->SetTitle(titleString);
//            if(expHists[i]->GetMaximum() > yMax)
//                yMax=expHists[i]->GetMaximum();
//            expHists[i]->SetLineColor(colorList[(i)%12]);
//            expHists[i]->SetMarkerColor(colorList[(i)%12]);
//            expHists[i]->SetMarkerStyle(markerList[(i)%14]);
//            expHists[i]->SetLineWidth(2);
//            expHists[i]->SetMarkerSize(1);
//            expHists[i]->GetYaxis()->SetTitleOffset(1.8);
//            theLegend->AddEntry(expHists[i],expHists[i]->GetTitle(),"PE");
//        }
//        expHists[0]->GetYaxis()->SetRangeUser(yMin,yMax);
//        TString imageName="CalRawPkChnHist_OnPeakExt_"+seriesString+"_NoCuts.png";
//        TString imageTitle="Timing for " + seriesString;
//        if(seriesString[0]=='S')
//            imageTitle+= " Trigger on Channel "+triggerString;
//
//
//
//        expHists[0]->SetTitle("");
//        expHists[0]->Draw("P E1 X0");
//        for (int i = 1;i< numDetectors-1;i++)
//        {
//            expHists[i]->Draw("P E1 X0 same");
//        }
//
//        c2->Update();
//
//        TString funcString=d2str(yMax);
//        TF1 onPeakRegion("onPeakRegion",funcString,0.,800.);
//        TF1 offPeakRegion("offPeakRegion",funcString,800.,1300.);
//
//        onPeakRegion.SetFillStyle(3013);
//        onPeakRegion.SetFillColor(kGreen);
//        onPeakRegion.SetLineColor(kBlack);
//        onPeakRegion.Draw("f same");
//
//        offPeakRegion.SetFillStyle(3013);
//        offPeakRegion.SetFillColor(kRed);
//        offPeakRegion.SetLineColor(kBlack);
//        offPeakRegion.Draw("f same");
//
//        theLegend->Draw();
//        theLegend->SetBorderSize(1);
//        theLegend->SetFillColor(kWhite);
//
//        TPaveText* textBox = new TPaveText(.05,.925,.95,.995,"NDC");
//        textBox->SetTextSize(.03);
//        textBox->AddText(imageTitle);
//        textBox->SetFillColor(kWhite);
//        textBox->SetLineColor(kWhite);
//
//        textBox->SetShadowColor(kWhite);
//        textBox->SetFillStyle(4000);//transparent
//        textBox->SetLineWidth(0);
//        textBox->Draw();
//
//        TString fullPath=GRAPHS_DIR;
//        fullPath+=imageName;
//        c2->SaveAs(fullPath);
//
//        delete c2;
//        for (int i = 0;i< numDetectors-1;i++)
//        {
//        	delete expHists[i];
//        }
//
//
//    }

///Analyze with Khodyuk prebake (working correctly with deconvolution)

//    createDerivedResultsFilesForCal(10,"Ba133_Center2_130408");
//    createDerivedResultsFilesForCal(10,"Ba133_DoubleSource_Center2_130408");
//    createDerivedResultsFilesForCal(10,"Cs137_DoubleSource_Center2_130408");

//    analyzeLOofCalPeaks("C461", "Center2_130408", true,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C461", "Center2_130408", true,CALPEAKTIMINGCUT_ONPEAKEXT);
//    analyzeLOofCalPeaks("C439", "Center2_130408", true,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C439", "Center2_130408", true,CALPEAKTIMINGCUT_ONPEAKEXT);

	//  analyzeLOofCalPeaks("C504", "DoubleSource_Center2_130408", true,CALPEAKTIMINGCUT_OFFPEAK);  //Not working yet
//    analyzeLOofCalPeaks("C504", "DoubleSource_Center2_130408", true,CALPEAKTIMINGCUT_ONPEAKEXT);

///Create co-LO plot with Khodyuk prebake
//    const int numHists=4;
//    const int colorList[12]= {kBlue, kBlue, kRed,kRed ,kBlack,kBlue, kRed, kRed, kRed, kGreen,kGreen,kGreen};
//    const int markerList[14]={21,20,21,20,20,22,32,33,3,5,28,26,30,29};
//    TString legendString[numHists];
//    for (int peakNum=0;peakNum<1;peakNum++)
//    {
//        TH1* loHists[numHists];
//
//
//        loHists[0]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C461OffPeak_Center2_130408_Peak1.txt",13,.5,13.5);
//        legendString[0]="C461_OffPeak";
//
//        loHists[1]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C439OffPeak_Center2_130408_Peak1.txt",13,.5,13.5);
//        legendString[1]="C439_OffPeak";
//
//        loHists[2]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C461OffPeak_Center2_130408_KhodyukPrebake_Peak1.txt",13,.5,13.5);
//        legendString[2]="C461_OffPeak_KhodyukPrebake";
//
//        loHists[3]=getTabSeperatedHist(TString(HISTS_DIR)+"LOForIndivDets_Ba133_C439OffPeak_Center2_130408_KhodyukPrebake_Peak1.txt",13,.5,13.5);
//        legendString[3]="C439_OffPeak_KhodyukPrebake";
//
//        TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1600,1000);
//        theCanvas.SetFillColor(kGray);
//        theCanvas.cd();
//        theCanvas.SetLeftMargin(.12);
//        theCanvas.SetRightMargin(.2);
//        theCanvas.SetBottomMargin(.06);
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
//        double legendSize=.05*numHists;
//        TLegend* theLegend = new TLegend(0.8,.9-legendSize,0.99,0.9);
//        theLegend->SetBorderSize(1);
//        theLegend->SetFillColor(kWhite);
//
//        for (int i = 0;i< numHists;i++)
//        {
//            loHists[i]->GetXaxis()->SetLabelOffset(99);
//            loHists[i]->SetLineColor(colorList[i%12]);
//            loHists[i]->SetMarkerColor(colorList[i%12]);
//            loHists[i]->SetMarkerStyle(markerList[i%14]);
//            //loHists[i]->SetMarkerStyle(21);
//            loHists[i]->GetYaxis()->SetTitle("Light Output");
//            loHists[i]->SetMarkerSize(1.2);
//            loHists[i]->GetXaxis()->Set(13,.45+i*.03,13.45+i*.03);
//            theLegend->AddEntry(loHists[i],legendString[i],"PE");
//        }
//
//
//
//        double xAxisTitleYPosition;
//        if(peakNum==0)
//        {
//            loHists[0]->GetYaxis()->SetRangeUser(0.7,.9);
//            xAxisTitleYPosition=.69;
//
//        }
//        else if(peakNum==1)
//        {
//            loHists[0]->GetYaxis()->SetRangeUser(0.9,1.0);
//            xAxisTitleYPosition=.895;
//        }
//        else
//        {
//            loHists[0]->GetYaxis()->SetRangeUser(0.96,1.01);
//            xAxisTitleYPosition=.958;
//        }
//        loHists[0]->GetXaxis()->SetRangeUser(0,16);
//        TString plotTitleString="Comparison of LO Analysis for Peak "+int2str(peakNum+1);
//        loHists[0]->SetTitle(plotTitleString);
//        loHists[0]->GetXaxis()->SetTimeOffset(2.2);
//
//        gStyle->SetErrorX(0);
//        loHists[0]->Draw("P E1");
//
//        for (int i = 1;i< numHists;i++)
//        {
//        	loHists[i]->Draw("P E1 same");
//        }
//
//        theLegend->Draw();
//
//        //Create a text based X axis
//        TText xAxisTText;
//        //xAxisTText.SetTextAngle(90);
//        xAxisTText.SetTextSize(0.03);
//        xAxisTText.SetTextAlign(21);
//        for (int i=0;i<13;i++)
//        {
//            Float_t x = loHists[0]->GetXaxis()->GetBinCenter(i+1);
//            TString axisTitle=int2str(i+1);
//            if(i==12)
//            {
//              axisTitle="Avg.";
//            }
//            xAxisTText.DrawText(x,xAxisTitleYPosition,axisTitle);
//        }
//
//        TString graphPathString=GraphsDir+"LOAnalysisComparisonPeak"+int2str(peakNum+1)+"130417.png";
//        theCanvas.Draw();
//
//        theCanvas.SaveAs(graphPathString);
//
//        for (int i = 0;i< numHists;i++)
//        {
//        	delete loHists[i];
//        }
//
//
//    }

///Plot LO from different detector setups
//    const int numGraphs=4;
//
//    const int numDewarPoints=14;
//    double dewarEnergies[numDewarPoints]={122.1,81,59.5,31,99.8,88,23.3,122.1,99.8,88,23.3,122.1,185.7,107};
//    double dewarLOs[numDewarPoints]={0.90033,0.91603,0.88772,0.81225,0.86602,0.90747,0.73221,0.8827,0.86363,0.90637,0.7364,0.89272,0.93333,0.8203};
//    double dewarLOErrors[numDewarPoints]={0.010849,0.010122,0.010078,0.010542,0.012021,0.010652,0.03,0.010572,0.010722,0.010078,0.025,0.010022,0.015,0.02};
//    TGraphErrors dewarGraph(numDewarPoints,dewarEnergies,dewarLOs,NULL,dewarLOErrors);
//    dewarGraph.SetName("DewarGraph");
//    dewarGraph.SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    TString dewarLegendString="RDK2 BGO in Dewar [MCA] (Various Sources)";
//    dewarGraph.SetLineColor(kRed);
//    dewarGraph.SetMarkerColor(kRed);
//    dewarGraph.SetMarkerStyle(kFullCircle);
//    dewarGraph.SetMarkerSize(1.5);
//
//
//    const int numCalRDK2Points=3;
//    double calRDK2Energies[numCalRDK2Points]={32.3,79.7,356};
//    double calRDK2LOs[numCalRDK2Points]={0.81,0.93,0.986};
//    double calRDK2LOErrors[numCalRDK2Points]={0.02,0.01,0.002};
//    TGraphErrors calRDK2Graph(numCalRDK2Points,calRDK2Energies,calRDK2LOs,NULL,calRDK2LOErrors);
//    calRDK2Graph.SetName("CalRDK2Graph");
//    TString calRDK2LegendString="RDK2 BGO Array [DAQ] (Ba-133 Cal Data)";
//    calRDK2Graph.SetLineColor(kViolet);
//    calRDK2Graph.SetMarkerColor(kViolet);
//    calRDK2Graph.SetMarkerStyle(kFullSquare);
//    calRDK2Graph.SetMarkerSize(1.5);
//
//    TGraph* khodyukGraph=getKhodyukBGOResponseGraphOfficial();
//    TString khodyukLegendString="I. Khodyuk and P. Dorenbos (monochromator)";
//    khodyukGraph->SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    khodyukGraph->SetLineColor(kBlack);
//    khodyukGraph->SetMarkerColor(kBlack);
//    khodyukGraph->SetMarkerStyle(kFullCircle);
//
//    const int numMoszynskiPoints=30;
//    double moszynskiEnergies[numMoszynskiPoints]={14.4,16.6,22,31,32,45,52,59.6,74,81,84,88,122,244.1,302.9,320.1,335.8,356,411.8,435,493.7,511,569.7,586,662,758.8,834.8,1063.7,1173.2,1275};
//    double moszynskiLOs[numMoszynskiPoints]={0.711,0.719,0.775,0.829,0.8412,0.853,0.9,0.893,0.91,0.926,0.937,0.933,0.915,0.964,0.979,0.964,0.97,0.986,0.97,0.995,1,0.996,1,0.999,1,1,1,1.01,1.013,1.01};
//    TGraph moszynskiGraph(numMoszynskiPoints,moszynskiEnergies,moszynskiLOs);
//    moszynskiGraph.SetName("MoszynskiGraph");
//    TString moszynskiLegendString="Moszynski et. al.";
//    moszynskiGraph.SetLineColor(kCyan);
//    moszynskiGraph.SetMarkerColor(kCyan);
//    moszynskiGraph.SetMarkerStyle(kFullSquare);
//    moszynskiGraph.SetMarkerSize(1.5);
//
//    int detectorNumber=1;
//    const int numQuadraticPoints=19;
//    double quadraticEnergies[numQuadraticPoints]={10,20,30,40,50,60,70,80,100,150,200,300,400,500,600,700,800,900,1000};
//    double quadraticLOs[numQuadraticPoints];
//
//    for (int i=0;i<numQuadraticPoints;i++)
//    {
//        quadraticLOs[i]=applyEDepModel(quadraticEnergies[i],detectorNumber, NULL, KEVIN_GAUS_EDEP)/quadraticEnergies[i];
//    }
//    TGraph quadraticGraph(numQuadraticPoints,quadraticEnergies,quadraticLOs);
//    quadraticGraph.SetName("QuadraticGraph");
//    TString quadraticLegendString="Old Quadratic Model";
//    quadraticGraph.SetLineColor(kGreen);
//    quadraticGraph.SetMarkerColor(kGreen);
//    quadraticGraph.SetMarkerStyle(kFullSquare);
//    quadraticGraph.SetMarkerSize(1.5);
//
//
//
//
//    double legendSize=.05*numGraphs;
//    TLegend* theLegend = new TLegend(0.8,.9-legendSize,0.99,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
//    theLegend->AddEntry(&quadraticGraph,quadraticLegendString,"L");
//    theLegend->AddEntry(khodyukGraph,khodyukLegendString,"P");
//    //theLegend->AddEntry(&moszynskiGraph,moszynskiLegendString,"P");
//    theLegend->AddEntry(&dewarGraph,dewarLegendString,"P");
//    theLegend->AddEntry(&calRDK2Graph,calRDK2LegendString,"P");
//
//    TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,800,600);
//    theCanvas.SetFillColor(kGray);
//    theCanvas.cd();
//    theCanvas.SetLeftMargin(.12);
//    theCanvas.SetRightMargin(.1);
//    theCanvas.SetBottomMargin(.1);
//    theCanvas.SetGridx();
//    theCanvas.SetGridy();
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//
//    khodyukGraph->GetXaxis()->SetLimits(0,400.);
//    khodyukGraph->GetYaxis()->SetRangeUser(.65,1.);
//
//
//    khodyukGraph->Draw("A P");
//    quadraticGraph.Draw("L same");
//   // moszynskiGraph.Draw("P same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//    theLegend->Draw("same");
//
//
//    TString imagePath=string(GRAPHS_DIR)+ "LOExperimentsGraph_130418.root";
//
//    theCanvas.SaveAs(imagePath);
//
//    delete khodyukGraph;
//    delete theLegend;

///Create Cumulative Graphs for Collab 130318
//cEvents theEvents;
//double energyIn=14.397;
//double lowerLimit=energyIn*evalMattKhodyukParam(energyIn);
//double lowEnergyBR=theEvents.calcBR(20,EKEMAX);

//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM);
//
//    double theoryBRForExpRange=0.002960299; //14.397-endpoint which is the experimental data range
//
////    double theoryBRForMCEventsRange=0.00252823+0.00318548+0.00429224; //.2 - endpoint (unused at the moment
//
//    TH1D* lowToHighCumulativeIntegralHist=makeTH1DFromDim("lowToHighCumulativeIntegralHist","Integral Low to High;BGO Energy Detected (keV);Exp-MC/MC",BGO_DIM);
//    TH1D* highToLowCumulativeIntegralHist=makeTH1DFromDim("highToLowCumulativeIntegralHist","Integral High to Low;BGO Energy Detected (keV);Exp-MC/MC",BGO_DIM);
//    TH1D* lowToHighCumulativeSpecFitHist=makeTH1DFromDim("lowToHighCumulativeSpecFitHist","SpecFit Low to High;BGO Energy Detected (keV);Extracted Branching Ratio for 14.39 keV to the Endpoint",BGO_DIM);
//    TH1D* highToLowCumulativeSpecFitHist=makeTH1DFromDim("highToLowCumulativeSpecFitHist","SpecFit High to Low;BGO Energy Detected (keV);Extracted Branching Ratio",BGO_DIM);
//
//    double expIntegral,expError;
//    double mcIntegral, mcError;
//    double spectralFit;
//    for (int i = 0;i< 800;i++)
//    {
//        //Low to High
//        expIntegral=expHist->IntegralAndError(0,i+1,expError);
//        mcIntegral=mcHist->IntegralAndError(0,i+1,mcError);
//        spectralFit=fitMCtoExp(expHist,mcHist,0,i+1);
//        if(mcIntegral != 0)
//        {
//            lowToHighCumulativeIntegralHist->SetBinContent(i+1,(expIntegral/mcIntegral)-1);
//            lowToHighCumulativeIntegralHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION));
//        }
//
//        if(spectralFit > 0)
//        {
//            lowToHighCumulativeSpecFitHist->SetBinContent(i+1,spectralFit -1.);
//            lowToHighCumulativeSpecFitHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION)); //using integral error for now
//        }
//
//        //High to low
//        expIntegral=expHist->IntegralAndError(i+1,800,expError);
//        mcIntegral=mcHist->IntegralAndError(i+1,800,mcError);
//        spectralFit=fitMCtoExp(expHist,mcHist,i+1,800);
//        if(mcIntegral != 0)
//        {
//            highToLowCumulativeIntegralHist->SetBinContent(i+1,(expIntegral/mcIntegral)-1.);
//            highToLowCumulativeIntegralHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION));
//        }
//
//        if(spectralFit > 0)
//        {
//            highToLowCumulativeSpecFitHist->SetBinContent(i+1,spectralFit-1.);
//            highToLowCumulativeSpecFitHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION)); //using integral error for now
//        }
//
//
//    }
//
//    convertTH1ToTXT(lowToHighCumulativeIntegralHist,HistsDir+"CumulativeExpMCComparison_LowToHighIntegral_130419.txt");
//    convertTH1ToTXT(lowToHighCumulativeSpecFitHist,HistsDir+"CumulativeExpMCComparison_LowToHighSpecFit_130419.txt");
//    convertTH1ToTXT(highToLowCumulativeIntegralHist,HistsDir+"CumulativeExpMCComparison_HighToLowIntegral_130419.txt");
//    convertTH1ToTXT(highToLowCumulativeSpecFitHist,HistsDir+"CumulativeExpMCComparison_HighToLowSpecFit_130419.txt");

//    TH1D* lowToHighCumulativeIntegralHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_LowToHighIntegral_130419.txt",BGO_DIM);
//    TH1D* lowToHighCumulativeSpecFitHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_LowToHighSpecFit_130419.txt",BGO_DIM);
//    TH1D* highToLowCumulativeIntegralHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_HighToLowIntegral_130419.txt",BGO_DIM);
//    TH1D* highToLowCumulativeSpecFitHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_HighToLowSpecFit_130419.txt",BGO_DIM);
//
//    double x[2]={-1,801};
//    double y[2]={0,0};
//    TGraph zeroAxis(2,x,y);
//    zeroAxis.SetLineWidth(3);
//
//
//    TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,800,600);
//    theCanvas.SetFillColor(kGray);
//    theCanvas.cd();
//    theCanvas.SetLeftMargin(.12);
//    //theCanvas.SetRightMargin(.2);
//    theCanvas.SetBottomMargin(.1);
//    theCanvas.SetGridx();
//    theCanvas.SetGridy();
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//
//
//    lowToHighCumulativeIntegralHist->GetYaxis()->SetRangeUser(-.05, .2);
//    lowToHighCumulativeIntegralHist->GetYaxis()->SetTitleOffset(1.3);
//    highToLowCumulativeIntegralHist->GetYaxis()->SetRangeUser(-.05, 2);
//    highToLowCumulativeIntegralHist->GetYaxis()->SetTitleOffset(1.3);
//
//    lowToHighCumulativeIntegralHist->SetLineColor(kBlue);
//    highToLowCumulativeIntegralHist->SetLineColor(kBlue);
//    lowToHighCumulativeSpecFitHist->SetLineColor(kRed);
//    highToLowCumulativeSpecFitHist->SetLineColor(kRed);
//
//    lowToHighCumulativeIntegralHist->SetTitle("Cumulative Comparison of Exp & MC: 10 keV to X keV;Upper Limit on BGO Energy Detected (keV);#frac{Exp-MC}{MC}");
//    highToLowCumulativeIntegralHist->SetTitle("Cumulative Comparison of Exp & MC: X keV to 800 keV;Lower Limit on BGO Energy Detected (keV);#frac{Exp-MC}{MC}");
//
//
////    lowToHighCumulativeIntegralHist->Draw("P E1");
//    highToLowCumulativeIntegralHist->Draw("P E1");
////    lowToHighCumulativeSpecFitHist->Draw("P E1 same");
//    highToLowCumulativeSpecFitHist->Draw("P E1 same");
//    zeroAxis.Draw("same");
//
//    double legendSize=.05*3;
//    TLegend* theLegend = new TLegend(0.6,.5-legendSize,0.79,0.5);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
//
//
////    theLegend->AddEntry(lowToHighCumulativeIntegralHist,"Integral","LE");
//    theLegend->AddEntry(highToLowCumulativeIntegralHist,"Integral","LE");
////    theLegend->AddEntry(lowToHighCumulativeSpecFitHist,"Spectral Fit","LE");
//    theLegend->AddEntry(highToLowCumulativeSpecFitHist,"Spectral Fit","LE");
//    theLegend->Draw("L same");
//
//    TString imagePath=GraphsDir+"test2.root";
//
//    theCanvas.SaveAs(imagePath);
//
//
////    delete mcHist;
////    delete expHist;
//    delete lowToHighCumulativeIntegralHist;
//    delete highToLowCumulativeIntegralHist;

	///Plot Jeff Average det vs Matt Average det
//    const int numExp=3;
//    const TString jeffPath=HistsDir+"EXPJeff_Hist130308_EDepG_BGO_IndivCounts_EPG.txt";
//    const HistDim jeffDim={160,0,800};
//    TH1* expHists[numExp]={getTabSeperatedHist(jeffPath,jeffDim),getTabSeperatedHist(HistsDir+"EXPBen_Hist130319_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM),getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM)};
//    expHists[1]->Rebin(5);
//    expHists[2]->Rebin(5);
//    expHists[0]->Scale(8.683e-5/(double) expHists[0]->Integral()); //Jeff's hist is in average counts per BGO
//    expHists[0]->SetTitle("Exp (Jeff);BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
//    expHists[2]->SetTitle("Exp (Matt);BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
//    expHists[1]->SetTitle("Exp (Ben);BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)");
//
//    expHists[0]->SetName("jeffHist");
//    expHists[1]->SetName("benHist");
//    expHists[2]->SetName("mattHist");
//
//
//    TH1* mcHists[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM)};
//    mcHists[0]->Rebin(5);
//    mcHists[0]->SetTitle("MC (Matt);BGO Energy Detected (keV);#frac{EPG}{EP}   (count rate ratio x 10^{-6})");
//    mcHists[0]->SetName("mcHist");
//
//
////    for (int i = 0;i< 40;i++)
////    {
////        expHists[0]->SetBinContent(i+1,0);
////        expHists[1]->SetBinContent(i+1,0);
////        expHists[2]->SetBinContent(i+1,0);
////        mcHists[0]->SetBinContent(i+1,0);
////
////        expHists[0]->SetBinError(i+1,0);
////        expHists[1]->SetBinError(i+1,0);
////        expHists[2]->SetBinError(i+1,0);
////        mcHists[0]->SetBinError(i+1,0);
////
////    }
//
//    expHists[0]->GetXaxis()->SetRangeUser(0,800);
//    expHists[1]->GetXaxis()->SetRangeUser(0,800);
//    mcHists[0]->GetXaxis()->SetRangeUser(0,800);
//
//
//    TH1* benHist[1]={expHists[1]};
//
//
//    plotResidualsToImage(2, expHists, 1, mcHists,"Average of Individual Detectors - Residuals (Exp - MC)","EDepG_AverageIndiv_EXP_MC_Comparison_Residuals_130419.root",-2e-6,2e-6,false);
//
//
//    plotResidualsToImage(2, expHists, 1, mcHists,"Average of Individual Detectors - Normalized Residuals (#frac{Exp - MC}{MC})","EDepG_AverageIndiv_EXP_MC_Comparison_NormResiduals_130419.root",-1,1,true);
//    TString titleString;
//
//
//    double integral, integralError;
//
//    integral=mcHists[0]->IntegralAndError(1,160,integralError);
//    titleString="Monte Carlo: " + d2str(integral,4) +" #pm " + d2str(integralError,3);
//    mcHists[0]->SetTitle(titleString);
//
//    integral=expHists[0]->IntegralAndError(1,160,integralError);
//    titleString="Jeff (Exp): " + d2str(integral,3) +" #pm " + d2str(integralError,3);
//    expHists[0]->SetTitle(titleString);
//
//    integral=expHists[1]->IntegralAndError(1,160,integralError);
//    titleString="Ben (Exp): " + d2str(integral,3) +" #pm " + d2str(integralError,3);
//    expHists[1]->SetTitle(titleString);
//    plotExpVersusMCToImage(2, expHists, 1, mcHists,"Average of Individual Detectors - Comparison","EDepG_AverageIndiv_EXP_MC_Comparison_130419.root");

///Make single detector histograms
///Make all MC Hists 130213 - but now with Kevin Gauss
////    RDK2Set threeBodySet(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
////    RDK2Set fourBodySet(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
////
////    TH1* mcIndividualHists[12];
////
////
////    TString drawString;
////    TString PIDString="PID5";
////    TString groupIDString="G6";
////
////
////
////    TH1* mcAverageHist[1];
////
////    mcAverageHist[0]=makeTH1DFromDim("mcAverageHist","mcAverage;BGO Energy Detected (keV);#frac{EPG}{EP} (count rate ratio)",BGO_DIM);
////    mcAverageHist[0]->Sumw2();
////
////
////    TString titleString;
////
////    for (int i = 0;i< 12;i++)
////    {
////
////        //MC
////        TString drawString="gDChn.BGO"+int2str(i+1)+"EDepBlur";
////        TString cutString="gDChn.BGO"+int2str(i+1)+"EDepBlur > 10";
////        mcIndividualHists[i]=combineAndScaleRadSetsToHist(&threeBodySet,&fourBodySet, drawString,STD_EPCut,TCut(cutString),BGO_DIM);
//////        convertTH1ToTXT(mcIndividualHists[i],HistsDir+"MC_Hist130213_EDepG_BGO_Det"+int2str(i+1)+"_EPG_"+fourBodySet.GetGResultsID()+".txt");
//////        mcIndividualHists[i]=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_Det"+int2str(i+1)+"_EPG_"+fourBodySet.GetGResultsID()+".txt",BGO_DIM);
////        mcAverageHist[0]->Add(mcIndividualHists[i]);
////        delete mcIndividualHists[i];
////
////    }
////    mcAverageHist[0]->Scale(1./12.);
////    mcAverageHist[0]->Rebin(5);
////    convertTH1ToTXT(mcAverageHist[0],HistsDir+"MC_Hist130213_EDepGQuadratic_BGO_AverageIndiv_EPG_"+fourBodySet.GetGResultsID()+".txt");
//    TH1* mcAverageHist[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepGQuadratic_BGO_AverageIndiv_EPG_RID130.txt",160,0,800)};
//    mcAverageHist[0]->SetTitle("Quadratic;BGO Energy Detected (keV);#frac{EPG}{EP}   (count rate ratio x 10^{-6})");
//
//    TH1* khodyukHists[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM)};
//    khodyukHists[0]->Rebin(5);
//    khodyukHists[0]->SetTitle("Khodyuk;BGO Energy Detected (keV);#frac{EPG}{EP}   (count rate ratio x 10^{-6})");
//    khodyukHists[0]->SetName("khodyukHist");
//
//    TString titleString;
//    double integral,integralError;
//    integral=khodyukHists[0]->IntegralAndError(1,800,integralError);
//    titleString="Khodyuk: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    khodyukHists[0]->SetTitle(titleString);
//
//    integral=mcAverageHist[0]->IntegralAndError(1,800,integralError);
//    titleString="Quadratic: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    mcAverageHist[0]->SetTitle(titleString);
//
//
//    plotExpVersusMCToImage(1,khodyukHists,1,mcAverageHist,"Comparison of Quadratic and Khodyuk models","QuadraticKhodyukComparison_Comparison.root");
//    khodyukHists[0]->SetTitle("Total Normalized Residual");
//    plotResidualsToImage(1,khodyukHists,1,mcAverageHist,"Comparison of Quadratic and Khodyuk models #frac{khodyuk-quadratic}{quadratic}","QuadraticKhodyukComparison_NormalizedResiduals.root");

///Make mult1 hist with new geometry
////    RDK2Set threeBodySet(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
////    RDK2Set fourBodySet(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
////
////
////
////    TH1D* mult1Hist=makeTH1DFromDim("mult1Hist","",BGO_DIM);
////    mult1Hist->Sumw2();
////    for (int i = 0;i< 12;i++)
////    {
////        TString drawString="gDChn.BGO"+int2str(i+1)+"EDepBlurEdge";
////        TString cutString=drawString +"> 10 && gDChn.BGOmultiplicityEdge==1";
////        //TString cutString="";
////        TCut theGCut=TCut(cutString);
////    	TH1* tempHist=combineAndScaleRadSetsToHist(&threeBodySet,&fourBodySet,drawString,STD_EPCut,theGCut,BGO_DIM);
////    	mult1Hist->Add(tempHist);
////    	delete tempHist;
////    }
////
////    convertTH1ToTXT(mult1Hist,HistsDir+"MC_Hist130213_EDepG_BGO_M1_EPG_RID130.txt");
//    TH1D* mult1Hist={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_M1_EPG_RID130.txt",BGO_DIM)};
//
//    TH1D* averageHist={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM)};
//    averageHist->Scale(12);
//
//    TH1* histsA[1]={averageHist};
//    TH1* histsB[2]={mult1Hist,averageHist};
//
//    TString titleString;
//    double integral,integralError;
//    integral=histsA[0]->IntegralAndError(1,800,integralError);
//    titleString="Individual Detector Average x12: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsA[0]->SetTitle(titleString);
//
//    integral=histsB[0]->IntegralAndError(1,800,integralError);
//    titleString="Multiplicity 1 of Detector Array: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsB[0]->SetTitle(titleString);
//
//
//
//    plotExpVersusMCToImage(0,NULL,2,histsB,"Comparison of individual detector and M1 spectra","Mult1AndIndividualDetectors_130421.root");
//    histsA[0]->SetTitle("Total Norm Residual");
//    plotResidualsToImage(1,histsA,1,histsB,"Normalized Residuals of individual detector and M1 spectra #frac{indiv-mult1}{mult1}","Mult1AndIndividualDetectors_NormResiduals_130421.root");
//
//    delete histsA[0];
//    delete histsB[0];
///Make mult1 hist with template analysis
////    TH1D* mult1Hist=makeTH1DFromDim("mult1Hist","",BGO_DIM);
////    mult1Hist->Sumw2();
////    for (int i = 0;i< 12;i++)
////    {
////        TString drawString="BGOEDep"+int2str(i+1);
////        TString cutString=drawString+" > 10";
////        TCut theGCut=TCut(cutString);
////    	TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, theGCut, 1,BGO_DIM);
////    	mult1Hist->Add(tempHist);
////    	delete tempHist;
////    }
//
////    convertTH1ToTXT(mult1Hist,HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepG_BGO_M1_EPG.txt");
//    TH1D* mult1Hist={getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepG_BGO_M1_EPG.txt",BGO_DIM)};
//    TH1D* averageHist={getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM)};
//    averageHist->Scale(12);
//
//    TH1* histsA[1]={averageHist};
//    TH1* histsB[2]={mult1Hist,averageHist};
//
//    TString titleString;
//    double integral,integralError;
//    integral=histsA[0]->IntegralAndError(1,800,integralError);
//    titleString="Individual Detector Average x12: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsA[0]->SetTitle(titleString);
//
//    integral=histsB[0]->IntegralAndError(1,800,integralError);
//    titleString="Multiplicity 1 of Detector Array: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsB[0]->SetTitle(titleString);
//
//
//
//    plotExpVersusMCToImage(0,NULL,2,histsB,"Exp Comparison of individual detector and M1 spectra","Mult1AndIndividualDetectors_Exp_130421.root");
//    histsA[0]->SetTitle("Total Norm Residual");
//    plotResidualsToImage(1,histsA,1,histsB,"Exp Normalized Residuals of individual detector and M1 spectra #frac{indiv-mult1}{mult1}","Mult1AndIndividualDetectors_Exp_NormResiduals_130421.root");
//
//    delete histsA[0];
//    delete histsB[0];
///Plot new geometry vs old geometry
//    RDK2Set threeBodySet(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
//    RDK2Set fourBodySet(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//
//
//
//    TH1D* newGeomM1Hist=makeTH1DFromDim("newGeomM1Hist","",BGO_DIM);
//    newGeomM1Hist->Sumw2();
//    for (int i = 0;i< 12;i++)
//    {
//        TString drawString="gDChn.BGO"+int2str(i+1)+"EDepBlur";
//        TString cutString=drawString +"> 10 && gDChn.BGOmultiplicity==1";
//        //TString cutString="";
//        TCut theGCut=TCut(cutString);
//    	TH1* tempHist=combineAndScaleRadSetsToHist(&threeBodySet,&fourBodySet,drawString,STD_EPCut,theGCut,BGO_DIM);
//    	newGeomM1Hist->Add(tempHist);
//    	delete tempHist;
//    }
//
//    convertTH1ToTXT(newGeomM1Hist,HistsDir+"MC_Hist130213_EDepGQuadratic_BGO_M1_EPG_RID130.txt");
//
//
//    //TH1D* newGeomM1Hist=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepGQuadratic_BGO_M1_EPG_RID130.txt",BGO_DIM);
//    newGeomM1Hist->Rebin(5);
//    TString histFilePath=HistsDir+"hist_BGO_M1_RID55_2.04us_Combined.root";
//    TFile theFile(histFilePath);
//    TH1D* oldGeomM1Hist=(TH1D*) theFile.Get("hist_BGO_M1_Combined");
//    convertTH1ToTXT(oldGeomM1Hist,HistsDir+"hist_BGO_M1_Combined_RID55_2.04us.txt");
//    theFile.Close();
//    oldGeomM1Hist=getTabSeperatedHist(HistsDir+"hist_BGO_M1_Combined_RID55_2.04us.txt",160,0,800);
//    oldGeomM1Hist->Scale(1./55660500.);
//
//    TH1* histsA[1]={newGeomM1Hist};
//    TH1* histsB[2]={oldGeomM1Hist,newGeomM1Hist};
//
//    TString titleString;
//    double integral,integralError;
//    integral=histsA[0]->IntegralAndError(1,160,integralError);
//    titleString="New Geometry M1 Hist: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsA[0]->SetTitle(titleString);
//
//    integral=histsB[0]->IntegralAndError(1,160,integralError);
//    titleString="Old Geometry: "+d2str(integral,4)+" #pm "+d2str(integralError,4)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsB[0]->SetTitle(titleString);
//
//
//
//    plotExpVersusMCToImage(0,NULL,2,histsB,"Comparison of new and old geometry M=1","NewOldGeom_130421.root");
//    histsA[0]->SetTitle("Total Norm Residual");
//    plotResidualsToImage(1,histsA,1,histsB,"Normalized Residuals of new and old geometry M1 spectra #frac{new-old}{old}","NewOldGeom_NormResiduals_130421.root");
//
//    delete histsA[0];
//    delete histsB[0];

///Make Electron energy, proton energy, and time of flight hists for MC
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//
//    TH1D* tofEPGHist=makeTH1DFromDim("tofEPGHist","",TOFP_DIM);
//    tofEPGHist->Sumw2();
//
//    TH1D* eDepEEPGHist=makeTH1DFromDim("eDepEEPGHist","",EDEPE_DIM);
//    eDepEEPGHist->Sumw2();
//
//    TString drawString;
//    for (int i = 0;i< 12;i++)
//    {
//        //EPG hists
//        TString cutString="gDChn.BGO"+int2str(i+1)+"EDepBlurEdge" +"> 10";
//        TCut theGCut=TCut(cutString);
//        TH1* tempHist;
//
//        drawString="pChn.SBDTimeFirst";
//        tempHist=combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,drawString,STD_EPCut,theGCut,TOFP_DIM);
//        tofEPGHist->Add(tempHist);
//        delete tempHist;
//
//        drawString="eDChn.SBDEDepBlur";
//        tempHist=combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,drawString,STD_EPCut,theGCut,EDEPE_DIM);
//        eDepEEPGHist->Add(tempHist);
//        delete tempHist;
//
//    }
//    delete fourBodySet;
//
//    tofEPGHist->Scale(1/12.);
//    eDepEEPGHist->Scale(1/12.);
//    convertTH1ToTXT(tofEPGHist,HistsDir+"MC_Hist130213_ToFP_AverageIndiv_EPG_RID130.txt");
//    convertTH1ToTXT(eDepEEPGHist,HistsDir+"MC_Hist130213_EDepE_AverageIndiv_EPG_RID130.txt");
//
//
//    //EPHists
//    drawString="pChn.SBDTimeFirst";
//    TH1* tofEPHist=threeBodySet->makeHist("tempHist",drawString,STD_EPCut,TOFP_DIM);
//    convertTH1ToTXT(tofEPHist,HistsDir+"MC_Hist130213_ToFP_EP_RID130.txt");
//
//    drawString="eDChn.SBDEDepBlur";
//    TH1* eDepEEPHist=threeBodySet->makeHist("tempHist",drawString,STD_EPCut,EDEPE_DIM);
//    convertTH1ToTXT(eDepEEPHist,HistsDir+"MC_Hist130213_EDepE_EP_RID130.txt");
//
//    delete eDepEEPHist;
//    delete tofEPHist;
//    delete tofEPGHist;
//    delete eDepEEPGHist;
//    delete threeBodySet;

///Make Electron energy, proton energy, and time of flight hists for Exp
//    TH1D* tofEPGHist=makeTH1DFromDim("tofEPGHist","",TOFP_DIM);
//    tofEPGHist->Sumw2();
//
//    TH1D* eDepEEPGHist=makeTH1DFromDim("eDepEEPGHist","",EDEPE_DIM);
//    eDepEEPGHist->Sumw2();
//
//    TString drawString;
//
//    for (int i = 0;i< 12;i++)
//    {
//        TString cutString="BGOEDep"+int2str(i+1) + " >10";
//        TCut theGCut=TCut(cutString);
//
//        //TofP
//        drawString="SBDToFP";
//    	TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, theGCut, -(i+1),TOFP_DIM);
//    	tofEPGHist->Add(tempHist);
//    	delete tempHist;
//
//    	//EDepE
//        drawString="SBDEDepE";
//    	tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, theGCut, -(i+1),EDEPE_DIM);
//    	eDepEEPGHist->Add(tempHist);
//    	delete tempHist;
//
//    }
//
//    tofEPGHist->Scale(1/12.);
//    eDepEEPGHist->Scale(1/12.);
//    convertTH1ToTXT(tofEPGHist,HistsDir+"EXPMatt_Hist130213_PID5_G6_ToFP_BGO_AverageIndiv_EPG.txt");
//    convertTH1ToTXT(eDepEEPGHist,HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepE_BGO_AverageIndiv_EPG.txt");
//
//
//    //EPHists
//    drawString="SBDToFP";
//    tofEPHist=makeExpEPHist("PID5", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString,"",TOFP_DIM);
//    convertTH1ToTXT(tofEPHist,HistsDir+"EXPMatt_Hist130213_PID5_G6_ToFP_EP.txt");
//
//    drawString="SBDEDepE";
//    eDepEEPHist=makeExpEPHist("PID5", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString,"",EDEPE_DIM);
//    convertTH1ToTXT(eDepEEPHist,HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepE_EP.txt");
//
//    delete eDepEEPHist;
//    delete tofEPHist;
//    delete tofEPGHist;
//    delete eDepEEPGHist;

///Make plots of exp vs MC for ToF and Electron energy 130422
//    TH1* expHists[1];
//    TH1* mcHists[1];
//
//    //EPG
//
//    //ToFP
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_ToFP_BGO_AverageIndiv_EPG.txt",TOFP_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130213_ToFP_AverageIndiv_EPG_RID130.txt",TOFP_DIM);
//
//    expHists[0]->Rebin(5);
//    mcHists[0]->Rebin(5);
//
//    mcHists[0]->GetYaxis()->SetRangeUser(0,5.2e-6);
//
//    expHists[0]->SetTitle("Exp;Proton Time of Flight (s);EPG/EP (count rate ratio)");
//    mcHists[0]->SetTitle("MC;Proton Time of Flight (s);EPG/EP (count rate ratio)");
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Exp & Simulation Comparison: Proton Time of Flight","ToFP_EPG_AverageIndiv_EXP_MC_Comparison_130419.png");
//    expHists[0]->SetTitle("Normalized Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Exp & Simulation Normalized Residuals #frac{Exp - MC}{MC}: Proton Time of Flight","ToFP_EPG__AverageIndiv_EXP_MC_NormResiduals_130419.png");
//
//    cout << "ToFP:" << endl;
//    cout << "Exp: " << expHists[0]->Integral(1,TOFP_DIM.numBins) << endl;
//    cout << "MC: " << mcHists[0]->Integral(1,TOFP_DIM.numBins) << endl << endl;
//    delete expHists[0];
//    delete mcHists[0];
//
//    //EDepE
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepE_BGO_AverageIndiv_EPG.txt",EDEPE_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepE_AverageIndiv_EPG_RID130.txt",EDEPE_DIM);
//
//    expHists[0]->Rebin(10);
//    mcHists[0]->Rebin(10);
//
//    mcHists[0]->GetYaxis()->SetRangeUser(0,2.6e-6);
//
//    expHists[0]->SetTitle("Exp;Electron Energy (keV);EPG/EP (count rate ratio)");
//    mcHists[0]->SetTitle("MC;Electron Energy (keV);EPG/EP (count rate ratio)");
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Exp & Simulation Comparison: Electron Energy","EDepE_EPG_AverageIndiv_EXP_MC_Comparison_130419.png");
//    expHists[0]->SetTitle("Normalized Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Exp & Simulation Normalized Residuals #frac{Exp - MC}{MC}: Electron Energy","EDepE_EPG__AverageIndiv_EXP_MC_NormResiduals_130419.png");
//
//    cout << "EDepE:" << endl;
//    cout << "Exp: " << expHists[0]->Integral(1,EDEPE_DIM.numBins) << endl;
//    cout << "MC: " << mcHists[0]->Integral(1,EDEPE_DIM.numBins) << endl << endl;
//    delete expHists[0];
//    delete mcHists[0];
//
//    //EP
//
//    //ToFP
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_ToFP_EP.txt",TOFP_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130213_ToFP_EP_RID130.txt",TOFP_DIM);
//
//    expHists[0]->Rebin(5);
//    mcHists[0]->Rebin(5);
//
//    expHists[0]->Scale(1./expHists[0]->Integral(1,expHists[0]->GetNbinsX()));
//    mcHists[0]->Scale(1./mcHists[0]->Integral(1,mcHists[0]->GetNbinsX()));
//
//    expHists[0]->SetTitle("Exp;Proton Time of Flight (s);Normalized Counts (AU)");
//    mcHists[0]->SetTitle("MC;Proton Time of Flight (s);Normalized Counts (AU)");
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Exp & Simulation Comparison: Proton Time of Flight","ToFP_EP_EXP_MC_Comparison_130419.png");
//    expHists[0]->SetTitle("Normalized Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Exp & Simulation Normalized Residuals #frac{Exp - MC}{MC}: Proton Time of Flight","ToFP_EP_EXP_MC_NormResiduals_130419.png",-.2,.2);
//
//    cout << "ToFP:" << endl;
//    cout << "Exp: " << expHists[0]->Integral(1,TOFP_DIM.numBins) << endl;
//    cout << "MC: " << mcHists[0]->Integral(1,TOFP_DIM.numBins) << endl << endl;
//    delete expHists[0];
//    delete mcHists[0];
//
//    //EDepE
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepE_EP.txt",EDEPE_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepE_EP_RID130.txt",EDEPE_DIM);
//
//    expHists[0]->Rebin(5);
//    mcHists[0]->Rebin(5);
//
//    expHists[0]->Scale(1./expHists[0]->Integral(1,expHists[0]->GetNbinsX()));
//    mcHists[0]->Scale(1./mcHists[0]->Integral(1,mcHists[0]->GetNbinsX()));
//
//    expHists[0]->SetTitle("Exp;Electron Energy (keV);Normalized Counts (AU)");
//    mcHists[0]->SetTitle("MC;Electron Energy (keV);Normalized Counts (AU)");
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Exp & Simulation Comparison: Electron Energy","EDepE_EP_EXP_MC_Comparison_130419.png");
//    expHists[0]->SetTitle("Normalized Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Exp & Simulation Normalized Residuals #frac{Exp - MC}{MC}: Electron Energy","EDepE_EP_EXP_MC_NormResiduals_130419.png",-.2,.2);
//
//    cout << "EDepE:" << endl;
//    cout << "Exp: " << expHists[0]->Integral(1,EDEPE_DIM.numBins) << endl;
//    cout << "MC: " << mcHists[0]->Integral(1,EDEPE_DIM.numBins) << endl << endl;
//    delete expHists[0];
//    delete mcHists[0];

///Make BAPD Hists MC 130423
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//
//    TH1D* mcHist=makeTH1DFromDim("mcHist","",BAPD_DIM);
//    mcHist->Sumw2();
//
//    for (int i = 0;i< 3;i++)
//    {
//        //EPG hists
//        TString cutString="gChn.BAPD"+int2str(i+1)+"EDepTotal" +"> .3";
//        TCut theGCut=TCut(cutString);
//        TH1* tempHist;
//
//        TString drawString="gChn.BAPD"+int2str(i+1)+"EDepTotal";
//        tempHist=combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,drawString,STD_EPCut,theGCut,BAPD_DIM);
//        convertTH1ToTXT(tempHist,HistsDir+"MC_Hist130213_EDepG_BAPD_"+int2str(i+1)+"_EPG_RID130.txt");
//        mcHist->Add(tempHist);
//        delete tempHist;
//
//    }
//    delete fourBodySet;
//
//    mcHist->Scale(1/3.);
//
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130213_EDepG_BAPD_AverageIndiv_EPG_RID130.txt");
//
///Make BAPD Hists Exp
//    TH1D* expHist=makeTH1DFromDim("expHist","",BAPD_DIM);
//    expHist->Sumw2();
//    for (int i = 0;i< 3;i++)
//    {
//        TString drawString="BAPDEDep"+int2str(i+1);
//        TString cutString="";
//        TCut theGCut=TCut(cutString);
//    	TH1* tempHist=makeExpEPGHist(true,"PID5", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, theGCut, -(i+1),BAPD_DIM);
//    	convertTH1ToTXT(tempHist,HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepG_BAPD_Det"+int2str(i+1)+"_EPG.txt");
//    	expHist->Add(tempHist);
//    	delete tempHist;
//    }
//    expHist->Scale(1./3.);
//
//    convertTH1ToTXT(expHist,HistsDir+"EXPMatt_Hist130213_PID5_G6_EDepG_BAPD_AverageIndiv_EPG.txt");

///Create Cumulative Graphs for Collab 130318
//////cEvents theEvents;
//////double energyIn=14.397;
//////double lowerLimit=energyIn*evalMattKhodyukParam(energyIn);
//////double lowEnergyBR=theEvents.calcBR(20,EKEMAX);
////
////
////    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
////    TH1D* expHist=getTabSeperatedHist(HistsDir+"EXPJeff_Hist130422_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM);
////
////    expHist->Scale(8.683e-5/(double) expHist->Integral()); //Jeff's hist is the sum of the 12 individual detectors and he did not calculate EPG/EP
////
//////    double theoryBRForExpRange=0.002960299; //14.397-endpoint which is the experimental data range
////
//////    double theoryBRForMCEventsRange=0.00252823+0.00318548+0.00429224; //.2 - endpoint (unused at the moment
////
////    TH1D* lowToHighCumulativeIntegralHist=makeTH1DFromDim("lowToHighCumulativeIntegralHist","Integral Low to High;BGO Energy Detected (keV);Exp-MC/MC",BGO_DIM);
////    TH1D* highToLowCumulativeIntegralHist=makeTH1DFromDim("highToLowCumulativeIntegralHist","Integral High to Low;BGO Energy Detected (keV);Exp-MC/MC",BGO_DIM);
////    TH1D* lowToHighCumulativeSpecFitHist=makeTH1DFromDim("lowToHighCumulativeSpecFitHist","SpecFit Low to High;BGO Energy Detected (keV);Extracted Branching Ratio for 14.39 keV to the Endpoint",BGO_DIM);
////    TH1D* highToLowCumulativeSpecFitHist=makeTH1DFromDim("highToLowCumulativeSpecFitHist","SpecFit High to Low;BGO Energy Detected (keV);Extracted Branching Ratio",BGO_DIM);
////
////    double expIntegral,expError;
////    double mcIntegral, mcError;
////    double spectralFit;
////    for (int i = 0;i< 800;i++)
////    {
////        //Low to High
////        expIntegral=expHist->IntegralAndError(0,i+1,expError);
////        mcIntegral=mcHist->IntegralAndError(0,i+1,mcError);
////        spectralFit=fitMCtoExp(expHist,mcHist,0,i+1);
////        if(mcIntegral != 0)
////        {
////            lowToHighCumulativeIntegralHist->SetBinContent(i+1,(expIntegral/mcIntegral)-1);
////            lowToHighCumulativeIntegralHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION));
////        }
////
////        if(spectralFit > 0)
////        {
////            lowToHighCumulativeSpecFitHist->SetBinContent(i+1,spectralFit -1.);
////            lowToHighCumulativeSpecFitHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION)); //using integral error for now
////        }
////
////        //High to low
////        expIntegral=expHist->IntegralAndError(i+1,800,expError);
////        mcIntegral=mcHist->IntegralAndError(i+1,800,mcError);
////        spectralFit=fitMCtoExp(expHist,mcHist,i+1,800);
////        if(mcIntegral != 0)
////        {
////            highToLowCumulativeIntegralHist->SetBinContent(i+1,(expIntegral/mcIntegral)-1.);
////            highToLowCumulativeIntegralHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION));
////        }
////
////        if(spectralFit > 0)
////        {
////            highToLowCumulativeSpecFitHist->SetBinContent(i+1,spectralFit-1.);
////            highToLowCumulativeSpecFitHist->SetBinError(i+1,getBasicErrorProgation(expIntegral, mcIntegral, expError, mcError,DIVISION)); //using integral error for now
////        }
////
////
////    }
//
////    convertTH1ToTXT(lowToHighCumulativeIntegralHist,HistsDir+"CumulativeExpMCComparison_LowToHighIntegral_130423.txt");
////    convertTH1ToTXT(lowToHighCumulativeSpecFitHist,HistsDir+"CumulativeExpMCComparison_LowToHighSpecFit_130423.txt");
////    convertTH1ToTXT(highToLowCumulativeIntegralHist,HistsDir+"CumulativeExpMCComparison_HighToLowIntegral_130423.txt");
////    convertTH1ToTXT(highToLowCumulativeSpecFitHist,HistsDir+"CumulativeExpMCComparison_HighToLowSpecFit_130423.txt");
//
//    TH1D* lowToHighCumulativeIntegralHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_LowToHighIntegral_130423.txt",BGO_DIM);
//    TH1D* lowToHighCumulativeSpecFitHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_LowToHighSpecFit_130423.txt",BGO_DIM);
//    TH1D* highToLowCumulativeIntegralHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_HighToLowIntegral_130423.txt",BGO_DIM);
//    TH1D* highToLowCumulativeSpecFitHist=getTabSeperatedHist(HistsDir+"CumulativeExpMCComparison_HighToLowSpecFit_130423.txt",BGO_DIM);
//
//    double x[2]={-1,801};
//    double y[2]={0,0};
//    TGraph zeroAxis(2,x,y);
//    zeroAxis.SetLineWidth(3);
//
//
//    TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,800,600);
//    theCanvas.SetFillColor(kGray);
//    theCanvas.cd();
//    theCanvas.SetLeftMargin(.12);
//    //theCanvas.SetRightMargin(.2);
//    theCanvas.SetBottomMargin(.1);
//    theCanvas.SetGridx();
//    theCanvas.SetGridy();
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//
//
//    lowToHighCumulativeIntegralHist->GetYaxis()->SetRangeUser(-.05, .2);
//    lowToHighCumulativeIntegralHist->GetYaxis()->SetTitleOffset(1.3);
//    highToLowCumulativeIntegralHist->GetYaxis()->SetRangeUser(-.05, 2);
//    highToLowCumulativeIntegralHist->GetYaxis()->SetTitleOffset(1.3);
//
//    lowToHighCumulativeIntegralHist->SetLineColor(kBlue);
//    highToLowCumulativeIntegralHist->SetLineColor(kBlue);
//    lowToHighCumulativeSpecFitHist->SetLineColor(kRed);
//    highToLowCumulativeSpecFitHist->SetLineColor(kRed);
//
//    lowToHighCumulativeIntegralHist->SetTitle("Cumulative Comparison of Exp & MC: 10 keV to X keV;Upper Limit on BGO Energy Detected (keV);#frac{Exp-MC}{MC}");
//    highToLowCumulativeIntegralHist->SetTitle("Cumulative Comparison of Exp & MC: X keV to 800 keV;Lower Limit on BGO Energy Detected (keV);#frac{Exp-MC}{MC}");
//
//
////    lowToHighCumulativeIntegralHist->Draw("P E1");
//    highToLowCumulativeIntegralHist->Draw("P E1");
////    lowToHighCumulativeSpecFitHist->Draw("P E1 same");
//    highToLowCumulativeSpecFitHist->Draw("P E1 same");
//    zeroAxis.Draw("same");
//
//    double legendSize=.05*3;
//    TLegend* theLegend = new TLegend(0.6,.5-legendSize,0.79,0.5);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
//
//
// //   theLegend->AddEntry(lowToHighCumulativeIntegralHist,"Integral","LE");
//   theLegend->AddEntry(highToLowCumulativeIntegralHist,"Integral","LE");
////    theLegend->AddEntry(lowToHighCumulativeSpecFitHist,"Spectral Fit","LE");
//    theLegend->AddEntry(highToLowCumulativeSpecFitHist,"Spectral Fit","LE");
//    theLegend->Draw("L same");
//
//    TString imagePath=GraphsDir+"test2.root";
//
//    theCanvas.SaveAs(imagePath);
//
//
////    delete mcHist;
////    delete expHist;
//    delete lowToHighCumulativeIntegralHist;
//    delete highToLowCumulativeIntegralHist;
//

///Make plots of exp vs MC for BAPD Energy (Ben) 130423
//    TH1* expHists[1];
//    TH1* mcHists[1];
//
//    //EPG
//
//    //ToFP
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPBen_Hist130321_PID5_G6_EDepG_BAPD_AverageIndiv_EPG.txt",BAPD_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BAPD_AverageIndiv_EPG_RID130.txt",BAPD_DIM);
//
//   // expHists[0]->Rebin(5);
//   // mcHists[0]->Rebin(5);
//
//    //mcHists[0]->GetYaxis()->SetRangeUser(0,5.2e-6);
//    TString titleString;
//    double integral, integralError;
//
//    integral=expHists[0]->IntegralAndError(1,BAPD_DIM.numBins,integralError);
//    titleString= "Exp: "+d2str(integral,3)+" #pm "+d2str(integralError,1)+";BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
//    expHists[0]->SetTitle(titleString);
//
//    integral=mcHists[0]->IntegralAndError(1,BAPD_DIM.numBins,integralError);
//    titleString= "MC: "+d2str(integral,3)+" #pm "+d2str(integralError,1)+";BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
//    mcHists[0]->SetTitle(titleString);
//
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Exp & Simulation Comparison: BAPD Energy Detected","EDep_BAPD_EPG_AverageIndiv_EXP_MC_Comparison_130419.png");
//    expHists[0]->SetTitle("Normalized Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Exp & Simulation Normalized Residuals #frac{Exp - MC}{MC}: BAPD Energy Detected","EDepBAPD_EPG__AverageIndiv_EXP_MC_NormResiduals_130419.png");
//
//    delete expHists[0];
//    delete mcHists[0];

///Electron Energy + Gamma Energy for energy conservation
//
//    HistDim theHistDim={1000,0.,1000};
////    makeTH1DFromDim("expHist","",theHistDim);
////    TString drawString="SPECIAL_EDepG+SBDEDepE";
////    TString cutString="";
////    TCut theGCut=TCut(cutString);
////    TH1* expHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, theGCut, 0,theHistDim);
////    convertTH1ToTXT(expHist,HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck.txt");
//    TH1* expHist=getTabSeperatedHist(HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck.txt",theHistDim);
//
////    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
////    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
////
////    TH1* mcHist=combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,R_BGOTOTALEDGE+"+eDChn.SBDEDepBlur",STD_EPCut,"gDChn.BGOmultiplicityEdge > 0",theHistDim);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_EDepGPlusEDepE_EnergyConservationCheck.txt");
//    TH1* mcHist=getTabSeperatedHist(HistsDir+"MC_EDepGPlusEDepE_EnergyConservationCheck.txt",theHistDim);
//    TH1* mcHists[1]={mcHist};
//
//
//    TH1* expHists[1]={expHist};
//
//    double integral, integralError;
//
//    ofstream energyConservationOut;
//
//    TString outPath=HistsDir+"EnergyConservationIntegratedData.txt";
//
//    energyConservationOut.open(outPath);
//
//
//    energyConservationOut << "0-800 keV" << endl;
//    integral = expHist->IntegralAndError(11,800, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(11,800, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//    energyConservationOut << "0-1000 keV" << endl;
//    integral = expHist->IntegralAndError(11,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(11,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//
//    energyConservationOut << "800-1000 keV" << endl;
//    integral = expHist->IntegralAndError(801,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(801,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//    energyConservationOut << "850-1000 keV" << endl;
//    integral = expHist->IntegralAndError(851,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(851,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//
//    energyConservationOut << "900-1000 keV" << endl;
//    integral = expHist->IntegralAndError(901,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(901,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//    energyConservationOut.close();
//    expHist->Rebin(20);
//    mcHist->Rebin(20);
//    mcHist->GetYaxis()->SetRangeUser(-5e-6,50e-6);
//
//    mcHist->SetTitle("MC;Total Energy Detected (keV); EPG/EP (count rate ratio)");
//    expHist->SetTitle("Exp (Matt);Total Energy Detected (keV); EPG/EP (count rate ratio)");
//
//    plotExpVersusMCToImage(1, expHists, 1, mcHists, "Test for energy conservation violation in Exp (BGOEDepAll + SBDEDepE)", "EDepGPlusEDepE_EnergyConservationCheck.png");
//    expHists[0]->SetTitle("Total Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Test for energy conservation violation in Exp (BGOEDepAll + SBDEDepE) - Residuals","EDepGPlusEDepE_EnergyConservationCheck_Residuals.png",-2e-6,10e-6,false);
//    expHists[0]->SetTitle("Normalized Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Test for energy conservation violation in Exp (BGOEDepAll + SBDEDepE) - #frac{Exp-MC}{Exp}","EDepGPlusEDepE_EnergyConservationCheck_NormResiduals.png");

///Create data cycles 3/4 for Jeff

//    TH1* jeffHist[2]={getTabSeperatedHist(HistsDir+"Jeff_DataCycles3_4_130426.txt",BGO_DIM),getTabSeperatedHist(HistsDir+"Jeff_DataCycles5_8_130426.txt",BGO_DIM)};
//    jeffHist[0]->Scale(1./1.73e7);
//    jeffHist[1]->Scale(1./3.45e7);
//
//
//
//
////    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
////    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
////    TH1* mcHist[1]={combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,R_BGOTOTALEDGE,STD_EPCut,"gDChn.BGOmultiplicityEdge>0",BGO_DIM)};
////    convertTH1ToTXT(mcHist[0],HistsDir+"MC_Hist130213_EDepG_BGO_Sum_EPG_RID130.txt");
//    TH1* mcHist[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_Sum_EPG_RID130.txt",BGO_DIM)};
//
//    jeffHist[0]->Rebin(10);
//    jeffHist[1]->Rebin(10);
//    mcHist[0]->Rebin(10);
//
//    jeffHist[0]->SetTitle("Cycles 3-4;BGO Energy Detected (keV); #frac{EPG}{EP} (count rate ratio)");
//    jeffHist[1]->SetTitle("Cycles 5-8;BGO Energy Detected (keV); #frac{EPG}{EP} (count rate ratio)");
//    mcHist[0]->SetTitle("MC;BGO Energy Detected (keV); #frac{EPG}{EP} (count rate ratio)");
//
//
//    plotResidualsToImage(2,jeffHist,1, mcHist,"Normalized Residuals", "test.png");

///Make tof by series
//
//
//    TString drawString;
//
//    TH1* expEPGHists[NUM_EXP_SERIES_G6];
//    TH1* expEPHists[NUM_EXP_SERIES_G6];
//
//    for (int j = 0;j< NUM_EXP_SERIES_G6;j++)
//    {
//        int currentSeries=EXP_SERIES_ID_LIST_G6[j];
//        int currentSeriesArray[1]={currentSeries};
//
//        TString titleString="Series"+int2str(currentSeries)+"_EP";
//
//        drawString="SBDToFP";
//        expEPHists[j]=makeExpEPHist("PID5", 1, currentSeriesArray, drawString,"",TOFP_DIM);
//        expEPHists[j]->SetName(titleString);
//        expEPHists[j]->SetTitle(titleString);
//        expEPHists[j]->Scale(1./expEPHists[j]->Integral());
//        expEPHists[j]->Rebin(5);
//
////        expEPGHists[j]=makeTH1DFromDim(titleString+"_EPG",titleString+"_EPG",TOFP_DIM);
////        expEPGHists[j]->Sumw2();
////
////        for (int i = 0;i< 12;i++)
////        {
////            TString cutString="BGOEDep"+int2str(i+1) + " >10";
////            TCut theGCut=TCut(cutString);
////
////            //TofP
////            drawString="SBDToFP";
////            TH1* tempHist=makeExpEPGHist(false,"PID5", 1, currentSeriesArray, drawString, theGCut, -(i+1),TOFP_DIM);
////            expEPGHists[j]->Add(tempHist);
////            delete tempHist;
////        }
////
////        expEPGHists[j]->Scale(1/12.);
//    }
//
//
//
//    TH1* mcEPGHists[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_ToFP_AverageIndiv_EPG_RID130.txt",TOFP_DIM)};
//    TH1* mcEPHists[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_ToFP_EP_RID130.txt",TOFP_DIM)};
//    mcEPHists[0]->Scale(1./mcEPHists[0]->Integral());
//    mcEPHists[0]->Rebin(5);
//
//    mcEPHists[0]->SetTitle("mcHist;Proton Time of Flight (s); Counts (AU)");
//
////    plotExpVersusMCToImage(NUM_EXP_SERIES_G6,expEPGHists,1,mcEPGHists,"Proton Time of Flight-EPG","EPG_PTOF_bySeriesWMC.png");
//    plotExpVersusMCToImage(NUM_EXP_SERIES_G6,expEPHists,1,mcEPHists,"Proton Time of Flight EP - Colors = Series, MC = Black line","EP_PTOF_bySeriesWMC.png",true);

///Make tof by series for EPG

//    TString drawString;
//
//    TH1* expEPGHists[2];
//
//    TString titleStrings[2]={"Cycle3_4","Cycle5_8"};
//    int seriesArrayCounts[2]={NUM_EXP_SERIES_G7_Cycle3_4,NUM_EXP_SERIES_G7_Cycle5_8};
//    const int* seriesArrays[2]={EXP_SERIES_ID_LIST_G7_Cycle3_4,EXP_SERIES_ID_LIST_G7_Cycle5_8};
//
//    for (int j = 0;j< 2;j++)
//    {
//
//        expEPGHists[j]=makeTH1DFromDim(titleStrings[j],titleStrings[j],TOFP_DIM);
//        expEPGHists[j]->Sumw2();
//
//        for (int i = 0;i< 12;i++)
//        {
//            TString cutString="BGOEDep"+int2str(i+1) + " >10";
//            TCut theGCut=TCut(cutString);
//
//            //TofP
//            drawString="SBDToFP";
//            TH1* tempHist=makeExpEPGHist(false,"PID5", seriesArrayCounts[j], seriesArrays[j], drawString, theGCut, -(i+1),TOFP_DIM);
//            expEPGHists[j]->Add(tempHist);
//            delete tempHist;
//        }
//
//        expEPGHists[j]->Scale(1/12.);
//        TString histPath=HistsDir+"ExpMatt_Hists130427_G7_ToFP_EPG_"+titleStrings[j]+".png";
//        convertTH1ToTXT(expEPGHists[j],histPath);
//
//
////                expEPGHists[j]=getTabSeperatedHist(histPath,TOFP_DIM);
//                expEPGHists[j]->Rebin(25);
//                expEPGHists[j]->SetTitle(titleStrings[j]);
//
//    }
//
//
//
//
//    TH1* mcEPGHists[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_ToFP_AverageIndiv_EPG_RID130.txt",TOFP_DIM)};
//
//
//    mcEPGHists[0]->Rebin(25);
//
//    mcEPGHists[0]->SetTitle("mcHist;Proton Time of Flight (s); Counts (AU)");
//    mcEPGHists[0]->GetYaxis()->SetRangeUser(0,23e-6);
//
//    plotExpVersusMCToImage(2,expEPGHists,1,mcEPGHists,"Proton Time of Flight-EPG","EPG_PTOF_bySeriesWMC_NewSeries.png");
//    plotResidualsToImage(2,expEPGHists,1,mcEPGHists,"Proton Time of Flight-EPG: Normalized Residuals","EPG_PTOF_bySeriesWMC_NormResiduals_NewSeries.png");
	//plotExpVersusMCToImage(NUM_EXP_SERIES_G6,expEPHists,1,mcEPHists,"Proton Time of Flight EP - Colors = Series, MC = Black line","EP_PTOF_bySeriesWMC.png",true);

///Make events from cluster to match runs completed 130429
//    makeMacrosForEventsOnCluster(0,48,"MAKEEVENTS_3B",THREEBODY);
//    makeMacrosForEventsOnCluster(0,48,"MAKEEVENTS_4B",FOURBODY);

///Create Fermi Function modified runs
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID171",PROTON,THREEBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID172",ELECTRON,THREEBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID173",PROTON,FOURBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID174",ELECTRON,FOURBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID175",GAMMA,FOURBODY);

///Create no fields runs
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID176",PROTON,THREEBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID177",ELECTRON,THREEBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID178",PROTON,FOURBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID179",ELECTRON,FOURBODY);
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID180",GAMMA,FOURBODY);

///Create 0 mirror runs
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID181",PROTON,THREEBODY,"/MRK/det/setMirrorFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID182",ELECTRON,THREEBODY,"/MRK/det/setMirrorFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID183",PROTON,FOURBODY,"/MRK/det/setMirrorFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID184",ELECTRON,FOURBODY,"/MRK/det/setMirrorFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID185",GAMMA,FOURBODY,"/MRK/det/setMirrorFieldScale 0");

///Double check overlaps did nothing
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID186",GAMMA,FOURBODY,"");

///Create 0 HV SBD runs
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID187",PROTON,THREEBODY,"/MRK/det/setSBDHVFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID188",ELECTRON,THREEBODY,"/MRK/det/setSBDHVFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID189",PROTON,FOURBODY,"/MRK/det/setSBDHVFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID190",ELECTRON,FOURBODY,"/MRK/det/setSBDHVFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID191",GAMMA,FOURBODY,"/MRK/det/setSBDHVFieldScale 0");

///Create 0 Magnetic Field runs
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID192",PROTON,THREEBODY,"/MRK/det/setMagneticFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID193",ELECTRON,THREEBODY,"/MRK/det/setMagneticFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID194",PROTON,FOURBODY,"/MRK/det/setMagneticFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID195",ELECTRON,FOURBODY,"/MRK/det/setMagneticFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID196",GAMMA,FOURBODY,"/MRK/det/setMagneticFieldScale 0");

///Analyze high energy originating spectrum

//    TString titleString;
//
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//    TH1* histsB[1]={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM)};
//    histsB[0]->Rebin(20);
//    titleString="Standard;BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsB[0]->SetTitle(titleString);
//
//    TH1* expHists[1]={getTabSeperatedHist(HistsDir+"EXPJeff_Hist130422_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM)};
//    expHists[0]->Scale(8.683e-5/(double) expHists[0]->Integral()); //Jeff's hist is the sum of the 12 individual detectors and he did not calculate EPG/EP
//    expHists[0]->Rebin(20);
//    titleString="Exp;BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    expHists[0]->SetTitle(titleString);
//
//    for (int i = 0;i < 7;i++)
//    {
//
//
//
//
//        double fractionLessLowEnergyEvents=0.3+i*.1;
//
////        TH1D* mcHist=makeTH1DFromDim("mcHist","",BGO_DIM);
////        mcHist->Sumw2();
////        for (int i = 0;i< 12;i++)
////        {
////            //EPG hists
////            TString cutString="gChn.BGO"+int2str(i+1)+"EDepTotal" +"> 10";
////            cutString+=" && rndm < ((t.eg0/781.571)*(1.-"+d2str(fractionLessLowEnergyEvents)+") +" +d2str(fractionLessLowEnergyEvents)+")";
////            TCut theGCut=TCut(cutString);
////            TH1* tempHist;
////
////            TString drawString="gChn.BGO"+int2str(i+1)+"EDepTotal";
////            tempHist=combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,drawString,STD_EPCut,theGCut,BGO_DIM);
////            //convertTH1ToTXT(tempHist,HistsDir+"MC_Hist130213_EDepG_BGO_"+int2str(i+1)+"_EPG_RID130.txt");
////            mcHist->Add(tempHist);
////            delete tempHist;
////
////        }
////
////
////        mcHist->Scale(1./12.);
////
////        convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_LinearTheoryMod"+d2str(fractionLessLowEnergyEvents)+"_EPG_RID130.txt");
////
////        TH1* histsA[2]={histsB[0],mcHist};
//        TH1* histsA[2]={histsB[0],getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_LinearTheoryMod"+d2str(fractionLessLowEnergyEvents)+"_EPG_RID130.txt",BGO_DIM)};
//
//        histsA[1]->Rebin(20);
//        titleString="LinearMod"+d2str(fractionLessLowEnergyEvents)+";BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//        histsA[1]->SetTitle(titleString);
//
//        histsA[1]->Scale(histsB[0]->Integral()/histsA[1]->Integral());  //Scale modified MC to previous MC
//
//        titleString="LinearMod-Standard: Linear Mod "+d2str(fractionLessLowEnergyEvents);
//        TString imageFileName="LinearMod"+d2str(fractionLessLowEnergyEvents)+"_Residual.png";
//        plotResidualsToImage(1,&(histsA[1]),1,histsB,titleString,imageFileName,-3e-6,1e-6,false);
//
//        titleString="#frac{LinearMod-Standard}{Standard}: Linear Mod "+d2str(fractionLessLowEnergyEvents);
//        imageFileName="LinearMod"+d2str(fractionLessLowEnergyEvents)+"_NormResidual.png";
//        plotResidualsToImage(1,&(histsA[1]),1,histsB,titleString,imageFileName);
//
//        titleString="#frac{Exp-MC}{MC}: Linear Mod "+d2str(fractionLessLowEnergyEvents);
//        imageFileName="LinearMod"+d2str(fractionLessLowEnergyEvents)+"_ExpNormResidual.png";
//        plotResidualsToImage(1,expHists,2,histsA,titleString,imageFileName);
//
//        titleString="Exp-MC: Linear Mod "+d2str(fractionLessLowEnergyEvents);
//        imageFileName="LinearMod"+d2str(fractionLessLowEnergyEvents)+"_ExpResidual.png";
//        plotResidualsToImage(1,expHists,2,histsA,titleString,imageFileName,-0.5e-6,4e-6,false);
//
//        delete histsA[1];
//    }
//
//    delete expHists[0];
//    delete histsB[0];
//    delete threeBodySet;
//    delete fourBodySet;

///Double check that we are correctly creating a higher energy originating spectrum 130501
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//
//    double fractionLessLowEnergyEvents=1;
//    TString titleString;
//    TString drawString="t.eg0";
//    TString cutString="";
//    cutString="rndm < ((t.eg0/781.571)*(1.-"+d2str(fractionLessLowEnergyEvents)+") +" +d2str(fractionLessLowEnergyEvents)+")";
//
//    TH1* histsB[1]={combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,drawString,"",cutString,BGO_DIM)};
//    histsB[0]->Rebin(20);
//    titleString="Standard;BGO Energy Detected (keV);EPG/EP (count rate ratio)";
//    histsB[0]->SetTitle(titleString);
//
//    for (int i = 0;i < 7;i++)
//    {
//
//        fractionLessLowEnergyEvents=0.3+i*.1;
//
//        //EPG hists
//        cutString="rndm < ((t.eg0/781.571)*(1.-"+d2str(fractionLessLowEnergyEvents)+") +" +d2str(fractionLessLowEnergyEvents)+")";
//
//
//        TH1* histsA[1]={combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,drawString,"",cutString,BGO_DIM)};
//        histsA[0]->Rebin(20);
//        titleString="LinearMod_EG0_"+d2str(fractionLessLowEnergyEvents)+";Initial Radiative Energy (keV);Counts (Arbitrary)";
//        histsA[0]->SetTitle(titleString);
//        //histsA[0]->Scale(histsB[0]->Integral()/histsA[0]->Integral());  //Scale modified MC to previous MC
//
//        titleString="#frac{LinearMod-Standard}{Standard}: Linear Mod "+d2str(fractionLessLowEnergyEvents);
//        TString imageFileName="LinearMod_EG0_"+d2str(fractionLessLowEnergyEvents)+"_NormResidual.png";
//        plotResidualsToImage(1,histsA,1,histsB,titleString,imageFileName);
//
//        delete histsA[0];
//    }
//
//    delete histsB[0];
//    delete threeBodySet;
//    delete fourBodySet;

///Calculate new Fermi'd BR
//    cEvents theEvents;
////    cout << theEvents.calcBR(.1,EKEMAX) << endl; //0.0102123
//    cout << theEvents.calcBR(15,340) << endl;

///Make Fermi Functioned based results and compare
//    double newBranchingRatio=0.010212;
//    double oldBranchingRatio=0.010006;  //From 0.00429224 +0.00318548 +0.00252823 from old method of low, mid, and high energy event sets

//    double newToOldRatio=newBranchingRatio/oldBranchingRatio;
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD Rad W Fermi");
//    RDK2Set* oldFourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//
////    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut,BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130501_EDepG_AverageIndiv_EPG_RID175.txt");
//
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130501_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1* mcHists[2]={getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM),mcHist};
//
//    TH1* expHists[1]={getTabSeperatedHist(HistsDir+"EXPJeff_Hist130422_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM)};
//
//    expHists[0]->Scale(8.683e-5/(double) expHists[0]->Integral()); //Jeff's hist is the sum of the 12 individual detectors and he did not calculate EPG/EP
//
//    mcHists[0]->SetTitle("NoFermi;BGO Energy Detected (keV);EPG/EP (count rate ratio)");
//    mcHists[1]->SetTitle("WFermi;BGO Energy Detected (keV);EPG/EP (count rate ratio)");
//
//    TH1* eg0Hists[2]={oldFourBodySet->makeHist("oldEG0Hist","t.eg0","t.ee0 > 50",BGO_DIM),fourBodySet->makeHist("newEG0Hist","t.eg0","t.ee0 > 50",BGO_DIM)};
//    convertTH1ToTXT(eg0Hists[0],HistsDir+"MC_Hist130501_EG0_NoFermi_130429.txt");
//    convertTH1ToTXT(eg0Hists[1],HistsDir+"MC_Hist130501_EG0_Fermi_130429.txt");
//
////    TH1* eg0Hists[2]={getTabSeperatedHist(HistsDir+"MC_Hist130501_EG0_NoFermi_130429.txt",BGO_DIM),getTabSeperatedHist(HistsDir+"MC_Hist130501_EG0_Fermi_130429.txt",BGO_DIM)};
//    eg0Hists[0]->SetTitle("NoFermi;Photon Energy Emitted (keV);Counts (AU)");
//    eg0Hists[1]->SetTitle("WFermi;Photon Energy Emitted (keV);Counts (AU)");
//    eg0Hists[1]->Scale(newToOldRatio);
//
//    eg0Hists[0]->Rebin(10);
//    eg0Hists[1]->Rebin(10);
//    mcHists[0]->Rebin(10);
//    mcHists[1]->Rebin(10);
//    expHists[0]->Rebin(10);
//
//    TString titleString,imageFileName;
//
//    titleString="Fermi Function BGO Comparison#frac{Exp-MC}{MC}";
//    imageFileName="FermiFunctionEvents_ExpVsMCComparison_1305041_NormalizedResiduals.png";
//    plotResidualsToImage(1,expHists,2,mcHists,titleString,imageFileName);
//
//    titleString="Fermi Function BGO Comparison #frac{WithFermi-NoFermi}{NoFermi}";
//    imageFileName="FermiFunctionEvents_MCComparison_1305041_NormalizedResiduals.png";
//    plotResidualsToImage(1,mcHists,1,&(mcHists[1]),titleString,imageFileName);
//
//    titleString="Fermi Function Photon Events Comparison #frac{WithFermi-NoFermi}{NoFermi}";
//    imageFileName="FermiFunctionEvents_eg0Comparison_1305041_NormalizedResiduals.png";
//    plotResidualsToImage(1,&(eg0Hists[1]),1,eg0Hists,titleString,imageFileName,-.1,.1);
//
//    delete threeBodySet;
//    delete fourBodySet;
//    delete mcHists[0];
//    delete mcHists[1];
//    delete eg0Hists[0];
//    delete eg0Hists[1];
//    delete expHists[0];

///B Field off

//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID192","RID193","","NoMag NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID194","RID195","RID196","NoMag Rad W Fermi");
//
//    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut,BGO_DIM);
//
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130501_EDepG_AverageIndiv_EPG_RID196.txt");
//
//    delete threeBodySet;
//    delete fourBodySet;
//    delete mcHist;

///Remake Electron Energy + Gamma Energy for energy conservation 130501
//
//    makeExpEPCoincidenceFiles("PID6",NUM_EXP_SERIES_G6,EXP_SERIES_ID_LIST_G6); //Fixed electron energy max limits prebuilt into exp analysis
//    HistDim theHistDim={2000,0.,2000};
////    makeTH1DFromDim("expHist","",theHistDim);
////    TString drawString="SPECIAL_EDepG+SBDEDepE";
////    TString cutString="";
////    TCut theGCut=TCut(cutString);
////    TH1* expHist=makeExpEPGHist(false,"PID6", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, theGCut, 0,theHistDim);
////    convertTH1ToTXT(expHist,HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck.txt");
//    TH1* expHist=getTabSeperatedHist(HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck.txt",theHistDim);
//
////    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
////    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//
////    TH1* mcHist=combineAndScaleRadSetsToHist(threeBodySet,fourBodySet,R_BGOTOTALEDGE+"+eDChn.SBDEDepBlur",STD_EPCut,"gDChn.BGOmultiplicityEdge > 0",theHistDim);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_EDepGPlusEDepE_EnergyConservationCheck.txt");
//    TH1* mcHist=getTabSeperatedHist(HistsDir+"MC_EDepGPlusEDepE_EnergyConservationCheck.txt",theHistDim);
//    TH1* mcHists[1]={mcHist};
//
//
//    TH1* expHists[1]={expHist};
//
//    double integral, integralError;
//
//    ofstream energyConservationOut;
//
//    TString outPath=HistsDir+"EnergyConservationIntegratedData.txt";
//
//    energyConservationOut.open(outPath);
//
//
//    energyConservationOut << "0-800 keV" << endl;
//    integral = expHist->IntegralAndError(11,800, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(11,800, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//    energyConservationOut << "0-1000 keV" << endl;
//    integral = expHist->IntegralAndError(11,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(11,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//
//    energyConservationOut << "800-1000 keV" << endl;
//    integral = expHist->IntegralAndError(801,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(801,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//    energyConservationOut << "850-1000 keV" << endl;
//    integral = expHist->IntegralAndError(851,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(851,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//
//    energyConservationOut << "900-1000 keV" << endl;
//    integral = expHist->IntegralAndError(901,1000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(901,1000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//    energyConservationOut << "1000-2000 keV" << endl;
//    integral = expHist->IntegralAndError(1001,2000, integralError);
//    energyConservationOut << "Exp: " << integral << " +/- " << integralError << endl;
//    integral = mcHist->IntegralAndError(1001,2000, integralError);
//    energyConservationOut << "MC: " << integral << " +/- " << integralError << endl << endl;
//
//    energyConservationOut.close();
//    expHist->Rebin(20);
//    mcHist->Rebin(20);
//    mcHist->GetYaxis()->SetRangeUser(-5e-6,50e-6);
//
//    mcHist->SetTitle("MC;Total Energy Detected (keV); EPG/EP (count rate ratio)");
//    expHist->SetTitle("Exp (Matt);Total Energy Detected (keV); EPG/EP (count rate ratio)");
//
//    plotExpVersusMCToImage(1, expHists, 1, mcHists, "Test for energy conservation violation in Exp (BGOEDepAll + SBDEDepE)", "EDepGPlusEDepE_EnergyConservationCheck_Extended.png");
//    expHists[0]->SetTitle("Total Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Test for energy conservation violation in Exp (BGOEDepAll + SBDEDepE) - Residuals","EDepGPlusEDepE_EnergyConservationCheck_Residuals_Extended.png",-2e-6,10e-6,false);
//    expHists[0]->SetTitle("Normalized Residual");
//    plotResidualsToImage(1,expHists,1,mcHists,"Test for energy conservation violation in Exp (BGOEDepAll + SBDEDepE) - #frac{Exp-MC}{Exp}","EDepGPlusEDepE_EnergyConservationCheck_NormResiduals_Extended.png");
//    delete expHists[0];
//    delete mcHists[0];

///Plotting individual spectra of energy conservation check (incomplete...need to modify makeExpEPGHist
//    TH1* expHists[1];
//    TString drawString, cutString,titleString;
////
//    HistDim theHistDim={1000,0.,2000};
////
//    //BGOEDepAll
//    makeTH1DFromDim("expHist","",theHistDim);
//    drawString="SPECIAL_EDepG";
//    cutString="SPECIAL_EnergyConservationCheck";
//    expHists[0]=makeExpEPGHist(false,"PID6", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, cutString, 0,theHistDim);
//    convertTH1ToTXT(expHists[0],HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck_BGOEDepAll.txt");
////    expHists[0]=getTabSeperatedHist(HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck_BGOEDepAll.txt",theHistDim);
//    expHists[0]->SetTitle("ExpMatt;BGO Energy Detected (keV);#frac{EPG}{EP} (cout rate ratio)");
//    expHists[0]->Rebin(20);
//    plotExpVersusMCToImage(1,expHists,0,NULL,"BGOEDepAll of 800 < BGOEDepAll + SBDEDepE < 2000","BGOEDepAll_EDepGPlusEDepE_EnergyConservationCheck.png");
//    delete expHists[0];
//
//    //SBDEDepE
//    makeTH1DFromDim("expHist","",theHistDim);
//    drawString="SBDEDepE";
//    cutString="SPECIAL_EnergyConservationCheck";
//    expHists[0]=makeExpEPGHist(false,"PID6", NUM_EXP_SERIES_G6, EXP_SERIES_ID_LIST_G6, drawString, cutString, 0,theHistDim);
//    convertTH1ToTXT(expHists[0],HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck_SBDEDepE.txt");
////    expHists[0]=getTabSeperatedHist(HistsDir+"ExpMatt_EDepGPlusEDepE_EnergyConservationCheck_SBDEDepE.txt",theHistDim);
//    expHists[0]->SetTitle("ExpMatt;Electron Energy Detected (keV);#frac{EPG}{EP} (cout rate ratio)");
//    expHists[0]->Rebin(20);
//    plotExpVersusMCToImage(1,expHists,0,NULL,"SBDEDepE of 800 < SBDEDepE + SBDEDepE < 2000","SBDEDepE_EDepGPlusEDepE_EnergyConservationCheck.png");
//    delete expHists[0];

///Becchetti spectrometer work 130104
////
//    TString distanceString="18cm";
//    TLegend *theLegend;
//    theLegend = new TLegend(0.8,0.5,0.93,0.9);
//    TCanvas* c2 = new TCanvas("plotManyResolution","plotManyResolution",10,10,1024,768);
//    c2->SetFillColor(kGray);
//    c2->cd();
//    c2->SetLeftMargin(.12);
//    c2->SetRightMargin(.12);
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//    TString filePath=string(MRKGEANTDATA_DIR)+"spectrometer130104_"+distanceString+".root";
//
//    TFile specFile(filePath);
//
//    TTree* fluxTree=(TTree*) specFile.Get("SBDFluxTree");
//
//    HistDim energyDim={100,0.,1000.};
//    HistDim positionDim={300,0,30.0};
//
//    TString titleString="Single Magnet Detector D_{O} = "+distanceString+";Energy (keV);D_{I}=dM-3.53 (cm)";
//
//    TH2D* the2DHist=makeTH2DFromDim("SingleHist",titleString, energyDim,positionDim);
//
//
//
//    fluxTree->Draw("100 * PosZ + 1.27:KE >> SingleHist","KE > 5","COLZ");
//
//
//
//    TString imagePath=GraphsDir+"BecchettiSingleMagnetMC130104_"+distanceString+".png";
//    c2->SaveAs(imagePath);
//    delete c2;
//
//    ofstream outFile;
//
//    TString outFilePath=HistsDir+"SpectrometerKEDistance"+distanceString+"130104.txt";
//    outFile.open(outFilePath);
//    outFile << "zi\tPeak1\tPeak2" << endl;
//
//    const int numDistances=70;
//    TH1* distanceHists[numDistances];
//    for (int i = 10;i< numDistances+1;i++)
//    {
//
//        c2 = new TCanvas("plotManyResolution","plotManyResolution",10,10,1024,768);
//        c2->SetFillColor(kGray);
//        c2->cd();
//        c2->SetLeftMargin(.12);
//        c2->SetRightMargin(.05);
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//        TString nameString="distancehist"+int2str(i+1)+"cm";
//        TString titleString="D_{O} = "+distanceString+", D_{I} = "+int2str(i+1)+" cm;Energy (keV);counts (AU)";
//        distanceHists[i]=makeTH1DFromDim(nameString,titleString,energyDim);
//        the2DHist->ProjectionX(nameString,i+1,i+1);
//        distanceHists[i]->GetYaxis()->SetTitleOffset(1.8);
//
//        TF1 firstPeak("firstPeak","gaus",0.,120.);
//        firstPeak.SetParameter(0,1000);
//        firstPeak.SetParameter(1,50);
//        firstPeak.SetParameter(2,100);
//        TF1 secondPeak("secondPeak","gaus",120.,1000.);
//        secondPeak.SetParameter(0,1000);
//        secondPeak.SetParameter(1,150);
//        secondPeak.SetParameter(2,100);
//
//
//        distanceHists[i]->Fit(&firstPeak,"VMRN");
//        distanceHists[i]->Fit(&secondPeak,"VMR+");
//        double Peak1Location=firstPeak.GetParameter(1);
//        double Peak2Location=secondPeak.GetParameter(1);
//        outFile << (i+1)/10. << "\t" << Peak1Location << "\t" << Peak2Location << endl;
//        titleString="D_{O} = "+distanceString+", D_{I} = "+int2str(i+1)+" cm, Peak1="+d2str(Peak1Location,3)+" keV, Peak2="+d2str(Peak2Location,3)+" keV;Energy (keV);counts (AU)";
//        distanceHists[i]->SetTitle(titleString);
//        distanceHists[i]->Draw();
//        firstPeak.Draw("same");
//        secondPeak.Draw("same");
//        TString imagePath=GraphsDir+"SpectrometerDO_"+distanceString+"_DI_"+int2str(i+1)+"cmGraph.png";
//        //c2->SaveAs(imagePath);
//        delete c2;
//
//    }

//    outFile.close();

///Check if input of COMSOL field files works 130515

//    cVField3D theComsolField;
//    theComsolField.loadField("/home/mjbales/school/rdk/comsol/COMSOL_CompleteElectricField_Fine130515.txt",1);
//    convertCVField3DToROOT(&theComsolField,"/home/mjbales/school/rdk/comsol/COMSOL_CompleteElectricField_Fine130515.root","E");

//    cVRootField3D theRootField;
//
//    theRootField.loadField("/home/mjbales/school/rdk/comsol/COMSOL_CompleteElectricField_Fine130515.root","E",1);
//
//    TString response="";
//    while(response != "end")
//    {
//        cout << "Enter z position: ";
//        cin >> response;
//        double zPos=str2double(response);
//        cVector3D pos(0,0,zPos);
//        cVector3D result(0,0,0);
//        theRootField.linearInterp(pos,result);
//        cout << endl << "Field at 0 0 " << response << ": " << result.x << " " << result.y << " " << result.z << " V/m"  << endl;
//    }

//    createHistogramsForMCCalibraitonSourceData(10,"Co57_Offset3_130310");
//    createHistogramsForMCCalibraitonSourceData(10,"Cs137_Offset3_130310");

///Compare performance of ROOTField to CVField 130520
//    TString fieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_SBD_Matt4.txt";
//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_SBD_Matt5.txt";
////    cVField3D oldField;
////    oldField.setSymmetryInY(true);
////    oldField.loadField(fieldToTest,1);
//    cVField2D oldField;
//    oldField.setTransposedField(true);
//    oldField.loadField(fieldToTest,1);
////    ROOTField newField(fieldToTest,"solenoidField",1,3,3);
//    ROOTField newField2(newFieldToTest,"solenoidField",1,2,2);
////    newField.setSymmetry(1,true);
////    newField2.setSymmetry(1,true);
//    newField2.setRotation(0,-9.5*PI/180.,0);
//
//    TVector3 pos,newFieldOut,newFieldOut2;
//    cVector3D pos2,oldFieldOut;
//    for (int i = 0;i< 100;i++)
//    {
//
//        oldFieldOut.setVal(0,0,0);
//        newFieldOut.SetXYZ(0,0,0);
//        newFieldOut2.SetXYZ(0,0,0);
//        double x,y,z;
//        cout << "x: ";
//        cin >> x;
//        cout << endl << "y:";
//        cin >> y;
//        cout << endl << "z:";
//        cin >> z;
//        cout << endl;
//        pos.SetX(x);
//        pos.SetY(y);
//        pos.SetZ(z);
//
////        pos.SetX(gRandom->Rndm()*newField.getXLength()+newField.getXStart());
////        pos.SetY(gRandom->Rndm()*newField.getYLength()+newField.getYStart());
////        pos.SetZ(gRandom->Rndm()*newField.getZLength()+newField.getZStart());
//
//        pos2.setVal(pos.x(),pos.y(),pos.z());
//
//        oldField.linearInterpAs3D(pos2,oldFieldOut);
////        newField.linearInterp3D(pos,newFieldOut);
//        newField2.linearInterp3D(pos,newFieldOut2);
//
//
//        cout << "For location " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
//        cout << "Old Field: " << oldFieldOut.x<< " " << oldFieldOut.y<<" " << oldFieldOut.z << endl;
////        cout << "New Field: " << newFieldOut.x()<< " " << newFieldOut.y()<<" " << newFieldOut.z() << endl << endl;
//        cout << "New Field2: " << newFieldOut2.x()<< " " << newFieldOut2.y()<<" " << newFieldOut2.z() << endl << endl;
//
//    }
////    int numEvents=50000000;
////
////    time_t startTime, endTime;
////
////    startTime = time(NULL);
////
////    for (int i = 0;i< numEvents;i++)
////    {
////        oldFieldOut.setVal(0,0,0);
////        newFieldOut.SetXYZ(0,0,0);
////
////        pos.SetX(gRandom->Rndm()*newField2.getXLength()*1.1+newField2.getXStart()*.95);
////        pos.SetY(gRandom->Rndm()*newField2.getYLength()*1.1+newField2.getYStart()*.95);
////        pos.SetZ(gRandom->Rndm()*newField2.getZLength()*1.1+newField2.getZStart()*.95);
////
////        pos2.setVal(pos.x(),pos.y(),pos.z());
////
////        oldField.linearInterpAs3D(pos2,oldFieldOut);
////        //newField.linearInterp3D(pos,newFieldOut);
////
////    }
////    endTime=time(NULL);
////    cout << difftime(endTime,startTime) <<" seconds run time." << endl;
////
//////    startTime = time(NULL);
//////    for (int i = 0;i< numEvents;i++)
//////    {
//////        oldFieldOut.setVal(0,0,0);
//////        newFieldOut.SetXYZ(0,0,0);
//////
//////        pos.SetX(gRandom->Rndm()*newField.getXLength()*1.1+newField.getXStart()*.95);
//////        pos.SetY(gRandom->Rndm()*newField.getYLength()*1.1+newField.getYStart()*.95);
//////        pos.SetZ(gRandom->Rndm()*newField.getZLength()*1.1+newField.getZStart()*.95);
//////
//////        pos2.setVal(pos.x(),pos.y(),pos.z());
//////
//////        //oldField.linearInterp(pos2,oldFieldOut);
//////        newField.linearInterp3D(pos,newFieldOut);
//////
//////    }
//////    endTime=time(NULL);
//////    cout << difftime(endTime,startTime) <<" seconds run time." << endl;
////
////
////    startTime = time(NULL);
////    for (int i = 0;i< numEvents;i++)
////    {
////        oldFieldOut.setVal(0,0,0);
////        newFieldOut.SetXYZ(0,0,0);
////
////        pos.SetX(gRandom->Rndm()*newField2.getXLength()*1.1+newField2.getXStart()*.95);
////        pos.SetY(gRandom->Rndm()*newField2.getYLength()*1.1+newField2.getYStart()*.95);
////        pos.SetZ(gRandom->Rndm()*newField2.getZLength()*1.1+newField2.getZStart()*.95);
////
////        pos2.setVal(pos.x(),pos.y(),pos.z());
////
////        //oldField.linearInterp(pos2,oldFieldOut);
////        newField2.linearInterp3D(pos,newFieldOut);
////
////    }
////    endTime=time(NULL);
////    cout << difftime(endTime,startTime) <<" seconds run time." << endl;

///test field txt files to root files
//    TString fieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Bfield_RDK2_Matt2.txt";
////    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_SBD_Matt5.txt";
//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Bfield_RDK2_Matt2.root";
////    TString fieldToSaveAs="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_SBD_Matt5.root";
//    cVField3D oldField;
//    oldField.loadField(fieldToTest,1);
////    oldField.setSymmetryInY(true);
//    ROOTField newField(newFieldToTest,"solenoid",1.,3,3);
//    newField.setSymmetry(false,true,false);
////    newField.setRotation(0,-9.5*PI/180.,0);
////    newField.saveFieldToFile(fieldToSaveAs,"solenoid");
//
////    ROOTField newField2(fieldToSaveAs,"solenoid",1,3,3);
////    newField2.setSymmetry(false,true,false);
////    newField2.setRotation(0,-9.5*PI/180.,0);
//
//
//    TVector3 pos,newFieldOut,newFieldOut2;
//    TVector3 posInField;
//    cVector3D pos2,oldFieldOut;
//    double magOld,magNew;
//    int numInField=0;
//    for (int i = 0;i< 100000;i++)
//    {
//
//        oldFieldOut.setVal(0,0,0);
//        newFieldOut.SetXYZ(0,0,0);
////        newFieldOut2.SetXYZ(0,0,0);
//        double x,y,z;
//        cout << "x: ";
//        cin >> x;
//        cout << endl << "y:";
//        cin >> y;
//        cout << endl << "z:";
//        cin >> z;
//        cout << endl;
//        pos.SetX(x);
//        pos.SetY(y);
//        pos.SetZ(z);
//
////        pos.SetX(gRandom->Rndm()*.1-0.05);
////        pos.SetY(gRandom->Rndm()*.1-0.05);
////        pos.SetZ(gRandom->Rndm()-.5);
//
//        pos2.setVal(pos.x(),pos.y(),pos.z());
//
//        oldField.linearInterp(pos2,oldFieldOut);
//        newField.linearInterp3D(pos,newFieldOut);
////        newField2.linearInterp3D(pos,newFieldOut2);
//
////        posInField=pos;
////        posInField.RotateY(-9.5*PI/180.);
//
////        cout << "Old: " << -pos.Z()*SIN_BEND+pos.X()*COS_BEND << " " << pos.Y() << " " << pos.Z()*COS_BEND+pos.X()*SIN_BEND << " " << endl;
////        cout << "New: " << posInField.X() << " " << posInField.Y() << " " << posInField.Z() << " " << endl;
//
//
//
//
//        cout << "For location " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
//        cout << "Old Field: " << oldFieldOut.x<< " " << oldFieldOut.y<<" " << oldFieldOut.z << endl;
//        cout << "New Field: " << newFieldOut.x()<< " " << newFieldOut.y()<<" " << newFieldOut.z() << endl << endl;
////        cout << "New Field2: " << newFieldOut2.x()<< " " << newFieldOut2.y()<<" " << newFieldOut2.z() << endl << endl;
//        magOld=oldFieldOut.mag();
//        magNew=newFieldOut.Mag();
//        if(magOld > 0)
//        {
//
//            numInField++;
//            if(abs((magOld-magNew)/magOld) > 0.001)
//            {
//                cout << "Detected difference between fields" << endl;
//            }
//        }
//
//
//    }
//    cout << numInField << endl;

///Convert txt to root
	///SBD

//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_SBD_Matt5.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_SBD_Matt5.root";
//    ROOTField newField(newFieldToTest,"sbd",1,2,2);
//    newField.saveFieldToFile(fieldToSaveAs,"sbd");

	///Mirror

//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_Mirror_Matt4.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_Mirror_Matt4.root";
//    ROOTField newField(newFieldToTest,"mirror",1,2,2);
//    newField.saveFieldToFile(fieldToSaveAs,"mirror");

	///Solenoid
//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/Bfield_RDK2_Matt2.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/mrkgeantrun/Bfield_RDK2_Matt2.root";
//    ROOTField newField(newFieldToTest,"solenoid",1,3,3);
//    newField.saveFieldToFile(fieldToSaveAs,"solenoid");

	///BAPD

//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130528.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130528.root";
//    ROOTField newField(newFieldToTest,"bapd",1,3,3);
//    newField.saveFieldToFile(fieldToSaveAs,"bapd");

///Test BAPD field
////
//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130528.root";
//    ROOTField newField(newFieldToTest,"bapd",1,3,3);
//    TVector3 newFieldOut;
//    newField.linearInterp3D(TVector3(-.023,-.023,-0.039),newFieldOut);
//    newFieldOut.Print();

///Create 0 BAPD field with mrkgeant1.023
//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID197",PROTON,THREEBODY,"/MRK/det/setBAPDFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID198",ELECTRON,THREEBODY,"/MRK/det/setBAPDFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID199",PROTON,FOURBODY,"/MRK/det/setBAPDFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID200",ELECTRON,FOURBODY,"/MRK/det/setBAPDFieldScale 0");
//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID201",GAMMA,FOURBODY,"/MRK/det/setBAPDFieldScale 0");

///Create full BAPD field with mrkgeant1.023
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID202",PROTON,THREEBODY,"");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID203",ELECTRON,THREEBODY,"");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID204",PROTON,FOURBODY,"");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID205",ELECTRON,FOURBODY,"");
//makeGeantMacroFilesForClusterEventGenerator(0,48,"RID206",GAMMA,FOURBODY,"");

///Make EDepGTotal Hists for Kevin 130611
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
//
//    TString cutString, drawString;
//    for (int i = 0;i< 12;i++)
//    {
//        cutString="gChn.BGO"+int2str(i+1)+"EDepTotal > .1";
//
//        //Gamma Hist
//        drawString="gChn.BGO"+int2str(i+1)+"EDepTotal";
//        TH1* mcIndividualHist=combineAndScaleRadSetsToHist(threeBodySet,fourBodySet, drawString,STD_EPCut,TCut(cutString),BGOCAL_DIM);
//        convertTH1ToTXT(mcIndividualHist,HistsDir+"MC_Hist130611_EDepG_BGO_Det"+int2str(i+1)+"_EPG_"+fourBodySet->GetGResultsID()+".txt");
//        delete mcIndividualHist;
//
//
//    }

///Calc EP/EPG rates for timing regeimes
// RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad");
// RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD Rad");
// threeBodySet->calcEPRate(STD_EPCut && TCut("pChn.SBDTimeFirst < 10.e-6"));
// fourBodySet->calcEPGRateFor12IndividualBGOs(STD_EPCut);
// fourBodySet->calcEPGRateFor12IndividualBGOs(STD_EPCut && TCut("pChn.SBDTimeFirst < 10.e-6"));

///Create hists for different timing regimes with new cuts 130613
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD NonRad W Fermi");
////
////
////    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut,BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130617_EDepG_AverageIndiv_EPG_RID175.txt");
////
////    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_ECut&&TCut("pChn.SBDEDepTotal > 10. && pChn.SBDEDepTotal < 31. && pChn.SBDTimeFirst > 10.0e-6 && pChn.SBDTimeFirst < 25.e-6"),BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130617_EDepG_AverageIndiv_10us_25us_EPG_RID175.txt");
////    delete mcHist;
//
//    TH1* mcHists[3] = {getTabSeperatedHist(HistsDir+"MC_Hist130617_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM),getTabSeperatedHist(HistsDir+"MC_Hist130617_EDepG_AverageIndiv_10us_25us_EPG_RID175.txt",BGO_DIM),getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM)};
//
//    plotExpVersusMCToImage(0,NULL,3,mcHists,"test","test.png");

///Plot LO from different detector setups 130617
//
//    const int numGraphs=4;
//
//    const int numDewarPoints=14;
//    double dewarEnergies[numDewarPoints]={122.1,81,59.5,31,99.8,88,23.3,122.1,99.8,88,23.3,122.1,185.7,107};
//    double dewarLOs[numDewarPoints]={0.90033,0.91603,0.88772,0.81225,0.86602,0.90747,0.73221,0.8827,0.86363,0.90637,0.7364,0.89272,0.93333,0.8203};
//    double dewarLOErrors[numDewarPoints]={0.010849,0.010122,0.010078,0.010542,0.012021,0.010652,0.03,0.010572,0.010722,0.010078,0.025,0.010022,0.015,0.02};
//    TGraphErrors dewarGraph(numDewarPoints,dewarEnergies,dewarLOs,NULL,dewarLOErrors);
//    dewarGraph.SetName("DewarGraph");
//    dewarGraph.SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    TString dewarLegendString="RDK2 BGO in Dewar [MCA] (Various Sources)";
//    dewarGraph.SetLineColor(kRed);
//    dewarGraph.SetMarkerColor(kRed);
//    dewarGraph.SetMarkerStyle(kFullCircle);
//    dewarGraph.SetMarkerSize(1.5);
//
//
//    const int numCalRDK2Points=6;
//    double calRDK2Energies[numCalRDK2Points]={32.3,42.7,78.7,79.7,121.9,356};
//    double calRDK2LOs[numCalRDK2Points]={0.81,0.86,0.94,0.93,0.92,0.986};
//    double calRDK2LOErrors[numCalRDK2Points]={0.02,0.02,0.01,0.01,0.01,0.002};
//    TGraphErrors calRDK2Graph(numCalRDK2Points,calRDK2Energies,calRDK2LOs,NULL,calRDK2LOErrors);
//    calRDK2Graph.SetName("CalRDK2Graph");
//    TString calRDK2LegendString="RDK2 BGO Array [DAQ] (Matt's Cal & Jeff's 80keV)";
//    calRDK2Graph.SetLineColor(kViolet);
//    calRDK2Graph.SetMarkerColor(kViolet);
//    calRDK2Graph.SetMarkerStyle(kFullSquare);
//    calRDK2Graph.SetMarkerSize(1.5);
//
//    //TGraph* khodyukGraph=getKhodyukBGOResponseGraphOfficial();
//    //Assumed error was similar to moszynski and was 2%
//    const int numKhodyukPoints=375;
//    double khodyukEnergies[numKhodyukPoints]={11,11,12,12,13,13.025,13.05,13.075,13.1,13.125,13.15,13.175,13.2,13.225,13.25,13.275,13.3,13.325,13.35,13.375,
//        13.4,13.425,13.45,13.475,13.5,13.525,13.55,13.575,13.6,13.625,13.65,13.675,13.7,13.725,13.75,13.775,13.8,13.825,13.85,13.875,13.9,13.925,13.95,13.975,14,
//        14.025,14.05,14.075,14.1,14.125,14.15,14.175,14.2,14.225,14.25,14.275,14.3,14.325,14.35,14.375,14.4,14.425,14.45,14.475,14.5,14.525,14.55,14.575,14.6,14.625,
//        14.65,14.675,14.7,14.725,14.75,14.775,14.8,14.825,14.85,14.875,14.9,14.925,14.95,14.975,15,15.025,15.05,15.075,15.1,15.125,15.15,15.175,15.2,15.225,15.25,
//        15.275,15.3,15.325,15.35,15.375,15.4,15.425,15.45,15.475,15.5,15.525,15.55,15.575,15.6,15.625,15.65,15.675,15.7,15.725,15.75,15.775,15.8,15.825,15.85,15.875,
//        15.9,15.925,15.95,15.975,16,16.025,16.05,16.075,16.1,16.125,16.15,16.175,16.2,16.225,16.25,16.275,16.3,16.325,16.35,16.375,16.4,16.425,16.45,16.475,16.5,16.525,
//        16.55,16.575,16.6,16.625,16.65,16.675,16.7,16.725,16.75,16.775,16.8,16.825,16.85,16.875,16.9,16.925,16.95,16.975,17,17.025,17.05,17.075,17.1,17.125,17.15,17.175,
//        17.2,17.225,17.25,17.275,17.3,17.325,17.35,17.375,17.4,17.425,17.45,17.475,17.5,17.525,17.55,17.575,17.6,17.625,17.65,17.675,17.7,17.725,17.75,17.775,17.8,17.825,
//        17.85,17.875,17.9,17.925,17.95,17.975,18,18.025,18.05,18.075,18.1,18.125,18.15,18.175,18.2,18.225,18.25,18.275,18.3,18.325,18.35,18.375,18.4,18.425,18.45,18.475,
//        18.5,18.525,18.55,18.575,18.6,18.625,18.65,18.675,18.7,18.725,18.75,18.775,18.8,18.825,18.85,18.875,18.9,18.925,18.95,18.975,19,19.025,19.05,19.075,19.1,19.125,
//        19.15,19.175,19.2,19.225,19.25,19.275,19.3,19.325,19.35,19.375,19.4,19.425,19.45,19.475,19.5,19.525,19.55,19.575,19.6,19.625,19.65,19.675,19.7,19.725,19.75,19.775,
//        19.8,19.825,19.85,19.875,19.9,19.925,19.95,19.975,20,20,25,25,30,30,30,30,35,35,40,40,45,45,45,45,50,50,55,55,60,60,65,65,70,70,75,75,80,80,85,85,89.5,89.6,89.7,
//        89.8,89.9,90,90.1,90.2,90.3,90.4,90.5,90.6,90.7,90.8,90.9,91,91.1,91.2,91.3,91.4,91.5,91.6,91.7,91.8,91.9,92,92.1,92.2,92.3,92.4,92.5,92.6,92.7,92.8,92.9,93,93.1,
//        93.2,93.3,93.4,93.5,93.6,93.7,93.8,93.9,94,94.1,94.2,94.3,94.4,94.5,94.6,94.7,94.8,94.9,95,95,100,100};
//    double khodyukLOs[numKhodyukPoints]={0.700,0.707,0.725,0.724,0.740,0.742,0.737,0.743,0.743,0.741,0.737,0.740,0.738,0.736,0.743,0.746,0.743,0.739,0.742,0.742,
//        0.736,0.722,0.703,0.704,0.705,0.695,0.703,0.701,0.698,0.699,0.702,0.702,0.704,0.695,0.701,0.696,0.705,0.696,0.706,0.701,0.696,0.693,0.697,0.697,0.694,0.700,0.699,0.698,
//        0.699,0.702,0.697,0.692,0.693,0.696,0.697,0.689,0.697,0.683,0.684,0.693,0.687,0.698,0.695,0.702,0.699,0.697,0.697,0.699,0.697,0.699,0.694,0.691,0.699,0.691,0.692,0.694,
//        0.694,0.702,0.697,0.699,0.701,0.695,0.698,0.699,0.702,0.698,0.697,0.700,0.702,0.698,0.697,0.698,0.701,0.708,0.703,0.697,0.699,0.698,0.701,0.705,0.704,0.706,0.702,0.700,
//        0.704,0.704,0.706,0.712,0.701,0.705,0.704,0.704,0.702,0.698,0.697,0.702,0.703,0.698,0.699,0.693,0.701,0.700,0.701,0.700,0.694,0.697,0.698,0.698,0.699,0.695,0.701,0.704,
//        0.698,0.701,0.697,0.701,0.697,0.700,0.694,0.697,0.691,0.687,0.684,0.690,0.692,0.685,0.689,0.687,0.693,0.686,0.694,0.697,0.690,0.691,0.690,0.688,0.695,0.690,0.690,0.690,
//        0.690,0.694,0.693,0.693,0.690,0.688,0.692,0.687,0.689,0.692,0.687,0.694,0.695,0.692,0.693,0.693,0.698,0.694,0.695,0.692,0.699,0.694,0.693,0.695,0.695,0.696,0.692,0.694,
//        0.690,0.694,0.697,0.698,0.700,0.700,0.700,0.696,0.701,0.696,0.696,0.698,0.704,0.705,0.696,0.694,0.698,0.701,0.702,0.703,0.697,0.701,0.702,0.701,0.698,0.708,0.702,0.702,
//        0.701,0.702,0.699,0.706,0.705,0.700,0.714,0.699,0.710,0.705,0.704,0.702,0.704,0.703,0.703,0.703,0.705,0.706,0.703,0.701,0.712,0.701,0.710,0.707,0.709,0.704,0.709,0.713,
//        0.706,0.711,0.710,0.709,0.719,0.709,0.705,0.711,0.712,0.711,0.708,0.709,0.709,0.710,0.710,0.716,0.716,0.712,0.709,0.717,0.713,0.715,0.711,0.716,0.716,0.717,0.714,0.715,
//        0.711,0.720,0.714,0.718,0.715,0.711,0.720,0.716,0.717,0.722,0.722,0.724,0.720,0.718,0.768,0.770,0.803,0.800,0.800,0.800,0.827,0.826,0.844,0.845,0.861,0.864,0.858,0.858,
//        0.875,0.874,0.880,0.880,0.885,0.885,0.890,0.889,0.896,0.893,0.897,0.898,0.901,0.895,0.906,0.909,0.911,0.911,0.912,0.910,0.911,0.910,0.909,0.910,0.908,0.908,0.899,0.885,
//        0.885,0.884,0.885,0.881,0.883,0.880,0.882,0.881,0.879,0.879,0.878,0.881,0.879,0.879,0.880,0.880,0.882,0.882,0.877,0.878,0.878,0.880,0.876,0.879,0.880,0.881,0.879,0.878,
//        0.878,0.879,0.878,0.880,0.880,0.879,0.878,0.878,0.877,0.878,0.880,0.877,0.878,0.878,0.878,0.881,0.880,0.882,0.883};
//    double khodyukLOErrors[numKhodyukPoints]={0.014,0.014,0.015,0.014,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,
//        0.015,0.015,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.015,0.015,0.016,0.016,0.016,0.016,0.017,0.017,0.017,0.017,0.017,0.017,0.017,
//        0.017,0.018,0.017,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
//        0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
//        0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018};
//    TGraphErrors khodyukGraph(numKhodyukPoints,khodyukEnergies,khodyukLOs,NULL,khodyukLOErrors);
//    TString khodyukLegendString="I. Khodyuk and P. Dorenbos (monochromator)";
//    khodyukGraph.SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    khodyukGraph.SetLineColor(kBlack);
//    khodyukGraph.SetMarkerColor(kBlack);
//    khodyukGraph.SetMarkerStyle(kFullCircle);
//
//
//    const int numMoszynskiPoints=30;
//    double moszynskiEnergies[numMoszynskiPoints]={14.4,16.6,22,31,32,45,52,59.6,74,81,84,88,122,244.1,302.9,320.1,335.8,356,411.8,435,493.7,511,569.7,586,662,758.8,834.8,1063.7,1173.2,1275};
//    double moszynskiLOs[numMoszynskiPoints]={0.711,0.719,0.775,0.829,0.8412,0.853,0.9,0.893,0.91,0.926,0.937,0.933,0.915,0.964,0.979,0.964,0.97,0.986,0.97,0.995,1,0.996,1,0.999,1,1,1,1.01,1.013,1.01};
//    double moszynskiLOErrors[numMoszynskiPoints]={0.01422,0.01438,0.0155,0.01658,0.016824,0.02559,0.018,0.01786,0.0182,0.01852,0.01874,0.01866,0.0183,0.02892,0.01958,0.01928,0.0291,0.01972,0.0194,0.02985,0.03,0.01992,0.02,0.02997,0.02,0.03,0.02,0.0202,0.02026,0.0202};
//    TGraphErrors moszynskiGraph(numMoszynskiPoints,moszynskiEnergies,moszynskiLOs,NULL,moszynskiLOErrors);
//    moszynskiGraph.SetName("MoszynskiGraph");
//    TString moszynskiLegendString="Moszynski et. al.";
//    moszynskiGraph.SetLineColor(kCyan);
//    moszynskiGraph.SetMarkerColor(kCyan);
//    moszynskiGraph.SetMarkerStyle(kFullSquare);
//    moszynskiGraph.SetMarkerSize(1.5);
//
////    int detectorNumber=1;
//////    const int numQuadraticPoints=19;
//////    double quadraticEnergies[numQuadraticPoints]={10,20,30,40,50,60,70,80,100,150,200,300,400,500,600,700,800,900,1000};
//////    double quadraticLOs[numQuadraticPoints];
////
////    for (int i=0;i<numQuadraticPoints;i++)
////    {
////        quadraticLOs[i]=applyEDepModel(quadraticEnergies[i],detectorNumber, NULL, KEVIN_GAUS_EDEP)/quadraticEnergies[i];
////    }
////    TGraph quadraticGraph(numQuadraticPoints,quadraticEnergies,quadraticLOs);
////    quadraticGraph.SetName("QuadraticGraph");
////    TString quadraticLegendString="Old Quadratic Model";
////    quadraticGraph.SetLineColor(kGreen);
////    quadraticGraph.SetMarkerColor(kGreen);
////    quadraticGraph.SetMarkerStyle(kFullSquare);
////    quadraticGraph.SetMarkerSize(1.5);
//
//
//
//
//    double legendSize=.05*numGraphs;
//    TLegend* theLegend = new TLegend(0.6,.4-legendSize,0.99,0.4);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
////    theLegend->AddEntry(&quadraticGraph,quadraticLegendString,"L");
//    theLegend->AddEntry(&khodyukGraph,khodyukLegendString,"P");
//    theLegend->AddEntry(&moszynskiGraph,moszynskiLegendString,"P");
//    theLegend->AddEntry(&dewarGraph,dewarLegendString,"P");
//    theLegend->AddEntry(&calRDK2Graph,calRDK2LegendString,"P");
//
//    TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,1024);
//    theCanvas.SetFillColor(kGray);
//    theCanvas.cd();
//    theCanvas.SetLeftMargin(.12);
//    theCanvas.SetRightMargin(.1);
//    theCanvas.SetBottomMargin(.1);
//    theCanvas.SetGridx();
//    theCanvas.SetGridy();
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//
//    khodyukGraph.GetYaxis()->SetTitleOffset(1.3);
//    khodyukGraph.GetXaxis()->SetLimits(0,400.);
//    khodyukGraph.GetYaxis()->SetRangeUser(.65,1.);
//
//
//    khodyukGraph.Draw("A P E1");
////    quadraticGraph.Draw("L same");
//    moszynskiGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//    theLegend->Draw("same");
//
//    TString imagePath=string(GRAPHS_DIR)+ "LOExperimentsGraph_130617.png";
//
//    theCanvas.SaveAs(imagePath);
//
//    ///create master graph that has all the data in it with errors scaled to correspond with confidence and number of data points for fitting
//    const int numCumulativePoints=numKhodyukPoints+numMoszynskiPoints+numDewarPoints+numCalRDK2Points;
//    double cumulativeEnergies[numCumulativePoints];
//    double cumulativeLOs[numCumulativePoints];
//    double cumulativeLOErrors[numCumulativePoints];
//
//    double khodyukSigmaFactor=2;
//    double moszynskiSigmaFactor=2;
//    double dewarSigmaFactor=2;
//    double calRDK2SigmaFactor=1;
//
//    //Add Khodyuk graph
//    for (int i = 0;i< numKhodyukPoints;i++)
//    {
//        cumulativeEnergies[i]=khodyukEnergies[i];
//        cumulativeLOs[i]=khodyukLOs[i];
//        cumulativeLOErrors[i]=khodyukLOErrors[i]*sqrt(numKhodyukPoints)*khodyukSigmaFactor;
//    }
//
//    //Add moszynski graph
//    for (int i = 0;i< numMoszynskiPoints;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints]=moszynskiEnergies[i];
//        cumulativeLOs[i+numKhodyukPoints]=moszynskiLOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints]=moszynskiLOErrors[i]*sqrt(numMoszynskiPoints)*moszynskiSigmaFactor;
//    }
//
//    //Add dewar graph
//    for (int i = 0;i< numDewarPoints;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints+numMoszynskiPoints]=dewarEnergies[i];
//        cumulativeLOs[i+numKhodyukPoints+numMoszynskiPoints]=dewarLOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints+numMoszynskiPoints]=dewarLOErrors[i]*sqrt(numDewarPoints)*dewarSigmaFactor;
//    }
//
//    //Add RDK2 apparatus graph
//    for (int i = 0;i< numCalRDK2Points;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2Energies[i];
//        cumulativeLOs[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2LOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2LOErrors[i]*sqrt(numCalRDK2Points)*calRDK2SigmaFactor;
//    }
//
//
//    TGraphErrors cumulativeGraph(numCumulativePoints,cumulativeEnergies,cumulativeLOs,NULL,cumulativeLOErrors);
//
//    TCanvas* c2 = new TCanvas("fitFunc","fitFunc",10,10,1280,1024);
//    TH2F* hpx = new TH2F("hpx","Matt's Fit and Extrapolation to Energy Response Data;Energy (keV);Light Output Relative to 662",10,0,30,10,.3,1.05); // axis range
//    hpx->GetYaxis()->SetTitleOffset(1.3);
//    hpx->SetStats(kFALSE);   // no statistics
//    hpx->Draw();
//    cumulativeGraph.SetMarkerStyle(21);
//    cumulativeGraph.SetMarkerSize(.8);
//
//    khodyukGraph.Draw("P E1 same");
//    moszynskiGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//    theLegend->Draw("same");
//
//
//
//    int polynomialDegree=2;
//    const int numRegions=4;
//    TF1* regionFuncs[numRegions];
//    double regionEnergyDivisions[numRegions+1]={9,13.5,16,90,400};
//
//    int numRepeatFits=8;
//    TString fitString;
//
//    TString filePath=HistsDir+"EnergyResponseOutputParametersForCumulativeFit130617.txt";
//    ofstream outFile(filePath);
//
//    for (int i = 0;i< numRegions;i++)
//    {
//        TString regionName="region"+int2str(i)+"Func";
////        TString regionFunctionString="pol"+int2str(polynomialDegree);
//        TString regionFunctionString="-[0]/(1+[1]/x)";
//        regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i],regionEnergyDivisions[i+1]);
//        regionFuncs[i]->SetLineColor(kGreen);
//        regionFuncs[i]->SetParameter(0,100);
//        regionFuncs[i]->SetParameter(1,16);
////        if(i==0)
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[0],regionEnergyDivisions[i+2]);
////        }
////        else if(i==numRegions-1)
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i-1],regionEnergyDivisions[i+1]);
////        }
////        else
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i-1],regionEnergyDivisions[i+2]);
////        }
//    	for (int j = 0;j< numRepeatFits-1;j++)
//    	{
//            fitString="VMRN";
//            cumulativeGraph.Fit(regionFuncs[i],fitString);
//    	}
//
//        fitString="VMRN+";
//    	cumulativeGraph.Fit(regionFuncs[i],fitString);
//    	regionFuncs[i]->SetRange(regionEnergyDivisions[i],regionEnergyDivisions[i+1]);
//    	regionFuncs[i]->Draw("SAME");
//    	for (int j = 0;j< polynomialDegree+1;j++)
//    	{
//            outFile << regionFuncs[i]->GetParameter(j) << "\t";
//    	}
//    	outFile << endl;
//    }
//
//
//
//    regionFuncs[0]->SetRange(0,regionEnergyDivisions[1]);
////    regionFuncs[numRegions-1]->SetRange(regionEnergyDivisions[numRegions-1],120);
//
//    //Fit exponential high energy region
////    TString regionName="regionExpFunc";
////    TString regionFunctionString="-[0]*exp([1]*x+[2])+[3]";
////    TF1* expRegionFunc=new TF1(regionName,regionFunctionString,120*.95,800*1.05);
////    expRegionFunc->SetLineColor(kGreen);
////    expRegionFunc->SetParameter(0,-.285656766);
////    expRegionFunc->SetParameter(1,-.008672992);
////    expRegionFunc->SetParameter(2,0);
////    expRegionFunc->SetParameter(3,1.000916925);
////    for (int j = 0;j< numRepeatFits-1;j++)
////    {
////        fitString="VMRN";
////        cumulativeGraph.Fit(expRegionFunc,fitString);
////    }
////
////    fitString="VMRN+";
////    cumulativeGraph.Fit(expRegionFunc,fitString);
////    expRegionFunc->SetRange(120,800);
////    expRegionFunc->Draw("same");
////    for (int j = 0;j< 4;j++)
////    {
////        outFile << expRegionFunc->GetParameter(j) << "\t";
////    }
////    outFile << endl;
//
//
////    TF1 endFunc=TF1("endFunc","(-.285656766*exp(x*-.008672992)+1.000916925)",100,1000);
////    endFunc.Draw("SAME");
//
//    TString imagePath2=GraphsDir+"EnergyResponseOverallFit130617.png";
//    c2->SaveAs(imagePath2);
//    delete c2;
//
////    expRegionFunc->SetRange(120,800);
//
//    for (int i = 0;i< numRegions;i++)
//    {
//        delete regionFuncs[i];
//
//    }
//
//    outFile.close();
//
//
//
//
//    delete theLegend;

///Examine 0V data 130617

//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID181","RID182","","STD NonRad W Fermi 0V Mirror");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID183","RID184","RID185","STD NonRad W Fermi 0V Mirror");
//
//
//    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut,BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130613_EDepG_AverageIndiv_EPG_RID185.txt");
//
//    TH1* mcHists[2] = {getTabSeperatedHist(HistsDir+"MC_Hist130617_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM),getTabSeperatedHist(HistsDir+"MC_Hist130613_EDepG_AverageIndiv_EPG_RID185.txt",BGO_DIM)};
//
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"test","test.png");

///Convert COMSOL files to root

	///SBD
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_SBDElectricField_Fine130515.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/comsol/COMSOL_SBDElectricField_Fine130515.root";
//    ROOTField newField(newFieldToTest,"sbd",1,3,3);
//    newField.saveFieldToFile(fieldToSaveAs,"sbd");

	///Mirror
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_MirrorElectricField_Fine130515.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/comsol/COMSOL_MirrorElectricField_Fine130515.root";
//    ROOTField newField(newFieldToTest,"mirror",1,3,3);
//    newField.saveFieldToFile(fieldToSaveAs,"mirror");

///Create COMSOL electric field runs
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID207",PROTON,THREEBODY,"/MRK/det/setBAPDFieldScale 0\n/MRK/det/setMirrorFieldFilePath COMSOL_MirrorElectricField_Fine130515.root");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID208",ELECTRON,THREEBODY,"/MRK/det/setBAPDFieldScale 0\n/MRK/det/setMirrorFieldFilePath COMSOL_MirrorElectricField_Fine130515.root");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID209",PROTON,THREEBODY,"/MRK/det/setBAPDFieldScale 0\n/MRK/det/setSBDFieldFilePath COMSOL_SBDElectricField_Fine130515.root");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID210",ELECTRON,THREEBODY,"/MRK/det/setBAPDFieldScale 0\n/MRK/det/setSBDFieldFilePath COMSOL_SBDElectricField_Fine130515.root");

///Create macros to test Geant4.9.6.p02 130619
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID211",PROTON,THREEBODY,"/MRK/det/setBAPDFieldScale 0");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID212",ELECTRON,THREEBODY,"/MRK/det/setBAPDFieldScale 0");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID213",PROTON,FOURBODY,"/MRK/det/setBAPDFieldScale 0");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID214",ELECTRON,FOURBODY,"/MRK/det/setBAPDFieldScale 0");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID215",GAMMA,FOURBODY,"/MRK/det/setBAPDFieldScale 0");

///Create macros to test physics libraries in Geant4.9.6.p02 130621
//
//    TString standardLines="/MRK/det/setBAPDFieldScale 0";
//    standardLines+="\n/MRK/useManualPos true";
//    standardLines+="\n/MRK/useManualDir true";
//    standardLines+="\n/gun/position 0 0 22.4 cm";
//    standardLines+="\n/gun/direction 0.2588190451 0.9659258263 0";  //towards BGO 1
//
//    TString physicsList="emlivermore empenelope emstandard_opt0 emstandard_opt1 emstandard_opt2 emstandard_opt3 emstandard_opt4 standardSS standardSSM standardWVI standardGS";
//
//    int RIDNumber=216;
//    for (int i = 0;i< numItemsInStringList(physicsList);i++)
//    {
//        TString thePhysics;
//        stringFromList(physicsList,i,thePhysics);
//    	makeGeantMacroFilesForClusterEventGenerator(0,1,"RID"+int2str(RIDNumber),GAMMA,FOURBODY,standardLines+"\n/MRK/phys/addPhysics "+ thePhysics);
//    	RIDNumber++;
//    }

///Create central light response curve with correct region by region experiment weighting
//    const int numGraphs=4;
//
//    const int numDewarPoints=14;
//    double dewarEnergies[numDewarPoints]={122.1,81,59.5,31,99.8,88,23.3,122.1,99.8,88,23.3,122.1,185.7,107};
//    double dewarLOs[numDewarPoints]={0.90033,0.91603,0.88772,0.81225,0.86602,0.90747,0.73221,0.8827,0.86363,0.90637,0.7364,0.89272,0.93333,0.8203};
//    double dewarLOErrors[numDewarPoints]={0.010849,0.010122,0.010078,0.010542,0.012021,0.010652,0.03,0.010572,0.010722,0.010078,0.025,0.010022,0.015,0.02};
//    TGraphErrors dewarGraph(numDewarPoints,dewarEnergies,dewarLOs,NULL,dewarLOErrors);
//    dewarGraph.SetName("DewarGraph");
//    dewarGraph.SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    TString dewarLegendString="RDK2 BGO in Dewar [MCA] (Various Sources)";
//    dewarGraph.SetLineColor(kRed);
//    dewarGraph.SetMarkerColor(kRed);
//    dewarGraph.SetMarkerStyle(kFullCircle);
//    dewarGraph.SetMarkerSize(1.5);
//
//
//    const int numCalRDK2Points=6;
//    double calRDK2Energies[numCalRDK2Points]={32.3,42.7,78.7,79.7,121.9,356};
//    double calRDK2LOs[numCalRDK2Points]={0.81,0.86,0.94,0.93,0.92,0.986};
//    double calRDK2LOErrors[numCalRDK2Points]={0.02,0.02,0.01,0.01,0.01,0.002};
//    TGraphErrors calRDK2Graph(numCalRDK2Points,calRDK2Energies,calRDK2LOs,NULL,calRDK2LOErrors);
//    calRDK2Graph.SetName("CalRDK2Graph");
//    TString calRDK2LegendString="RDK2 BGO Array [DAQ] (Matt's Cal & Jeff's 80keV)";
//    calRDK2Graph.SetLineColor(kViolet);
//    calRDK2Graph.SetMarkerColor(kViolet);
//    calRDK2Graph.SetMarkerStyle(kFullSquare);
//    calRDK2Graph.SetMarkerSize(1.5);
//
//    //TGraph* khodyukGraph=getKhodyukBGOResponseGraphOfficial();
//    //Assumed error was similar to moszynski and was 2%
//    const int numKhodyukPoints=375;
//    double khodyukEnergies[numKhodyukPoints]={11,11,12,12,13,13.025,13.05,13.075,13.1,13.125,13.15,13.175,13.2,13.225,13.25,13.275,13.3,13.325,13.35,13.375,
//        13.4,13.425,13.45,13.475,13.5,13.525,13.55,13.575,13.6,13.625,13.65,13.675,13.7,13.725,13.75,13.775,13.8,13.825,13.85,13.875,13.9,13.925,13.95,13.975,14,
//        14.025,14.05,14.075,14.1,14.125,14.15,14.175,14.2,14.225,14.25,14.275,14.3,14.325,14.35,14.375,14.4,14.425,14.45,14.475,14.5,14.525,14.55,14.575,14.6,14.625,
//        14.65,14.675,14.7,14.725,14.75,14.775,14.8,14.825,14.85,14.875,14.9,14.925,14.95,14.975,15,15.025,15.05,15.075,15.1,15.125,15.15,15.175,15.2,15.225,15.25,
//        15.275,15.3,15.325,15.35,15.375,15.4,15.425,15.45,15.475,15.5,15.525,15.55,15.575,15.6,15.625,15.65,15.675,15.7,15.725,15.75,15.775,15.8,15.825,15.85,15.875,
//        15.9,15.925,15.95,15.975,16,16.025,16.05,16.075,16.1,16.125,16.15,16.175,16.2,16.225,16.25,16.275,16.3,16.325,16.35,16.375,16.4,16.425,16.45,16.475,16.5,16.525,
//        16.55,16.575,16.6,16.625,16.65,16.675,16.7,16.725,16.75,16.775,16.8,16.825,16.85,16.875,16.9,16.925,16.95,16.975,17,17.025,17.05,17.075,17.1,17.125,17.15,17.175,
//        17.2,17.225,17.25,17.275,17.3,17.325,17.35,17.375,17.4,17.425,17.45,17.475,17.5,17.525,17.55,17.575,17.6,17.625,17.65,17.675,17.7,17.725,17.75,17.775,17.8,17.825,
//        17.85,17.875,17.9,17.925,17.95,17.975,18,18.025,18.05,18.075,18.1,18.125,18.15,18.175,18.2,18.225,18.25,18.275,18.3,18.325,18.35,18.375,18.4,18.425,18.45,18.475,
//        18.5,18.525,18.55,18.575,18.6,18.625,18.65,18.675,18.7,18.725,18.75,18.775,18.8,18.825,18.85,18.875,18.9,18.925,18.95,18.975,19,19.025,19.05,19.075,19.1,19.125,
//        19.15,19.175,19.2,19.225,19.25,19.275,19.3,19.325,19.35,19.375,19.4,19.425,19.45,19.475,19.5,19.525,19.55,19.575,19.6,19.625,19.65,19.675,19.7,19.725,19.75,19.775,
//        19.8,19.825,19.85,19.875,19.9,19.925,19.95,19.975,20,20,25,25,30,30,30,30,35,35,40,40,45,45,45,45,50,50,55,55,60,60,65,65,70,70,75,75,80,80,85,85,89.5,89.6,89.7,
//        89.8,89.9,90,90.1,90.2,90.3,90.4,90.5,90.6,90.7,90.8,90.9,91,91.1,91.2,91.3,91.4,91.5,91.6,91.7,91.8,91.9,92,92.1,92.2,92.3,92.4,92.5,92.6,92.7,92.8,92.9,93,93.1,
//        93.2,93.3,93.4,93.5,93.6,93.7,93.8,93.9,94,94.1,94.2,94.3,94.4,94.5,94.6,94.7,94.8,94.9,95,95,100,100};
//    double khodyukLOs[numKhodyukPoints]={0.700,0.707,0.725,0.724,0.740,0.742,0.737,0.743,0.743,0.741,0.737,0.740,0.738,0.736,0.743,0.746,0.743,0.739,0.742,0.742,
//        0.736,0.722,0.703,0.704,0.705,0.695,0.703,0.701,0.698,0.699,0.702,0.702,0.704,0.695,0.701,0.696,0.705,0.696,0.706,0.701,0.696,0.693,0.697,0.697,0.694,0.700,0.699,0.698,
//        0.699,0.702,0.697,0.692,0.693,0.696,0.697,0.689,0.697,0.683,0.684,0.693,0.687,0.698,0.695,0.702,0.699,0.697,0.697,0.699,0.697,0.699,0.694,0.691,0.699,0.691,0.692,0.694,
//        0.694,0.702,0.697,0.699,0.701,0.695,0.698,0.699,0.702,0.698,0.697,0.700,0.702,0.698,0.697,0.698,0.701,0.708,0.703,0.697,0.699,0.698,0.701,0.705,0.704,0.706,0.702,0.700,
//        0.704,0.704,0.706,0.712,0.701,0.705,0.704,0.704,0.702,0.698,0.697,0.702,0.703,0.698,0.699,0.693,0.701,0.700,0.701,0.700,0.694,0.697,0.698,0.698,0.699,0.695,0.701,0.704,
//        0.698,0.701,0.697,0.701,0.697,0.700,0.694,0.697,0.691,0.687,0.684,0.690,0.692,0.685,0.689,0.687,0.693,0.686,0.694,0.697,0.690,0.691,0.690,0.688,0.695,0.690,0.690,0.690,
//        0.690,0.694,0.693,0.693,0.690,0.688,0.692,0.687,0.689,0.692,0.687,0.694,0.695,0.692,0.693,0.693,0.698,0.694,0.695,0.692,0.699,0.694,0.693,0.695,0.695,0.696,0.692,0.694,
//        0.690,0.694,0.697,0.698,0.700,0.700,0.700,0.696,0.701,0.696,0.696,0.698,0.704,0.705,0.696,0.694,0.698,0.701,0.702,0.703,0.697,0.701,0.702,0.701,0.698,0.708,0.702,0.702,
//        0.701,0.702,0.699,0.706,0.705,0.700,0.714,0.699,0.710,0.705,0.704,0.702,0.704,0.703,0.703,0.703,0.705,0.706,0.703,0.701,0.712,0.701,0.710,0.707,0.709,0.704,0.709,0.713,
//        0.706,0.711,0.710,0.709,0.719,0.709,0.705,0.711,0.712,0.711,0.708,0.709,0.709,0.710,0.710,0.716,0.716,0.712,0.709,0.717,0.713,0.715,0.711,0.716,0.716,0.717,0.714,0.715,
//        0.711,0.720,0.714,0.718,0.715,0.711,0.720,0.716,0.717,0.722,0.722,0.724,0.720,0.718,0.768,0.770,0.803,0.800,0.800,0.800,0.827,0.826,0.844,0.845,0.861,0.864,0.858,0.858,
//        0.875,0.874,0.880,0.880,0.885,0.885,0.890,0.889,0.896,0.893,0.897,0.898,0.901,0.895,0.906,0.909,0.911,0.911,0.912,0.910,0.911,0.910,0.909,0.910,0.908,0.908,0.899,0.885,
//        0.885,0.884,0.885,0.881,0.883,0.880,0.882,0.881,0.879,0.879,0.878,0.881,0.879,0.879,0.880,0.880,0.882,0.882,0.877,0.878,0.878,0.880,0.876,0.879,0.880,0.881,0.879,0.878,
//        0.878,0.879,0.878,0.880,0.880,0.879,0.878,0.878,0.877,0.878,0.880,0.877,0.878,0.878,0.878,0.881,0.880,0.882,0.883};
//    double khodyukLOErrors[numKhodyukPoints]={0.014,0.014,0.015,0.014,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,
//        0.015,0.015,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.015,0.015,0.016,0.016,0.016,0.016,0.017,0.017,0.017,0.017,0.017,0.017,0.017,
//        0.017,0.018,0.017,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
//        0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
//        0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018};
//    TGraphErrors khodyukGraph(numKhodyukPoints,khodyukEnergies,khodyukLOs,NULL,khodyukLOErrors);
//    TString khodyukLegendString="I. Khodyuk and P. Dorenbos (monochromator)";
//    khodyukGraph.SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    khodyukGraph.SetLineColor(kGreen);
//    khodyukGraph.SetMarkerColor(kGreen);
//    khodyukGraph.SetMarkerStyle(kFullCircle);
//
//
//    const int numMoszynskiPoints=30;
//    double moszynskiEnergies[numMoszynskiPoints]={14.4,16.6,22,31,32,45,52,59.6,74,81,84,88,122,244.1,302.9,320.1,335.8,356,411.8,435,493.7,511,569.7,586,662,758.8,834.8,1063.7,1173.2,1275};
//    double moszynskiLOs[numMoszynskiPoints]={0.711,0.719,0.775,0.829,0.8412,0.853,0.9,0.893,0.91,0.926,0.937,0.933,0.915,0.964,0.979,0.964,0.97,0.986,0.97,0.995,1,0.996,1,0.999,1,1,1,1.01,1.013,1.01};
//    double moszynskiLOErrors[numMoszynskiPoints]={0.01422,0.01438,0.0155,0.01658,0.016824,0.02559,0.018,0.01786,0.0182,0.01852,0.01874,0.01866,0.0183,0.02892,0.01958,0.01928,0.0291,0.01972,0.0194,0.02985,0.03,0.01992,0.02,0.02997,0.02,0.03,0.02,0.0202,0.02026,0.0202};
//    TGraphErrors moszynskiGraph(numMoszynskiPoints,moszynskiEnergies,moszynskiLOs,NULL,moszynskiLOErrors);
//    moszynskiGraph.SetName("MoszynskiGraph");
//    TString moszynskiLegendString="Moszynski et. al.";
//    moszynskiGraph.SetLineColor(kCyan);
//    moszynskiGraph.SetMarkerColor(kCyan);
//    moszynskiGraph.SetMarkerStyle(kFullSquare);
//    moszynskiGraph.SetMarkerSize(1.5);
//
//
//
////    int detectorNumber=1;
//////    const int numQuadraticPoints=19;
//////    double quadraticEnergies[numQuadraticPoints]={10,20,30,40,50,60,70,80,100,150,200,300,400,500,600,700,800,900,1000};
//////    double quadraticLOs[numQuadraticPoints];
////
////    for (int i=0;i<numQuadraticPoints;i++)
////    {
////        quadraticLOs[i]=applyEDepModel(quadraticEnergies[i],detectorNumber, NULL, KEVIN_GAUS_EDEP)/quadraticEnergies[i];
////    }
////    TGraph quadraticGraph(numQuadraticPoints,quadraticEnergies,quadraticLOs);
////    quadraticGraph.SetName("QuadraticGraph");
////    TString quadraticLegendString="Old Quadratic Model";
////    quadraticGraph.SetLineColor(kGreen);
////    quadraticGraph.SetMarkerColor(kGreen);
////    quadraticGraph.SetMarkerStyle(kFullSquare);
////    quadraticGraph.SetMarkerSize(1.5);
//
//
//
//
//    double legendSize=.05*numGraphs;
//    TLegend* theLegend = new TLegend(0.6,.4-legendSize,0.99,0.4);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
////    theLegend->AddEntry(&quadraticGraph,quadraticLegendString,"L");
//    theLegend->AddEntry(&khodyukGraph,khodyukLegendString,"P");
//    theLegend->AddEntry(&moszynskiGraph,moszynskiLegendString,"P");
//    theLegend->AddEntry(&dewarGraph,dewarLegendString,"P");
//    theLegend->AddEntry(&calRDK2Graph,calRDK2LegendString,"P");
//
//    TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,1024);
//    theCanvas.SetFillColor(kGray);
//    theCanvas.cd();
//    theCanvas.SetLeftMargin(.12);
//    theCanvas.SetRightMargin(.1);
//    theCanvas.SetBottomMargin(.1);
//    theCanvas.SetGridx();
//    theCanvas.SetGridy();
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//
//    khodyukGraph.GetYaxis()->SetTitleOffset(1.3);
//    khodyukGraph.GetXaxis()->SetLimits(0,400.);
//    khodyukGraph.GetYaxis()->SetRangeUser(.65,1.);
//
//
//    khodyukGraph.Draw("A P E1");
////    quadraticGraph.Draw("L same");
//    moszynskiGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//    theLegend->Draw("same");
//
//    TString imagePath=string(GRAPHS_DIR)+ "LOExperimentsGraph_130626.png";
//
//    theCanvas.SaveAs(imagePath);
//
//    ///create master graph that has all the data in it with errors scaled to correspond with confidence and number of data points for fitting;
//
//
//    const int numCumulativePoints=numKhodyukPoints+numMoszynskiPoints+numDewarPoints+numCalRDK2Points;
//    double cumulativeEnergies[numCumulativePoints];
//    double cumulativeLOs[numCumulativePoints];
//    double cumulativeLOErrors[numCumulativePoints];
//
//    double khodyukSigmaFactor=2;
//    double moszynskiSigmaFactor=2;
//    double dewarSigmaFactor=2;
//    double calRDK2SigmaFactor=1;
//
//    //Add Khodyuk graph
//    for (int i = 0;i< numKhodyukPoints;i++)
//    {
//        cumulativeEnergies[i]=khodyukEnergies[i];
//        cumulativeLOs[i]=khodyukLOs[i];
//        cumulativeLOErrors[i]=khodyukLOErrors[i]*sqrt(numKhodyukPoints)*khodyukSigmaFactor;
//    }
//
//    //Add moszynski graph
//    for (int i = 0;i< numMoszynskiPoints;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints]=moszynskiEnergies[i];
//        cumulativeLOs[i+numKhodyukPoints]=moszynskiLOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints]=moszynskiLOErrors[i]*sqrt(numMoszynskiPoints)*moszynskiSigmaFactor;
//    }
//
//    //Add dewar graph
//    for (int i = 0;i< numDewarPoints;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints+numMoszynskiPoints]=dewarEnergies[i];
//        cumulativeLOs[i+numKhodyukPoints+numMoszynskiPoints]=dewarLOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints+numMoszynskiPoints]=dewarLOErrors[i]*sqrt(numDewarPoints)*dewarSigmaFactor;
//    }
//
//    //Add RDK2 apparatus graph
//    for (int i = 0;i< numCalRDK2Points;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2Energies[i];
//        cumulativeLOs[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2LOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2LOErrors[i]*sqrt(numCalRDK2Points)*calRDK2SigmaFactor;
//    }
//
//
//    TGraphErrors cumulativeGraph(numCumulativePoints,cumulativeEnergies,cumulativeLOs,NULL,cumulativeLOErrors);
//
//    TCanvas* c2 = new TCanvas("fitFunc","fitFunc",10,10,1280,1024);
//    TH2F* hpx = new TH2F("hpx","Fit and Extrapolation to Energy Response Data;Energy Incident (keV);Light Output Relative to 662 keV",10,0,662,10,.3,1.05); // axis range
//    hpx->GetYaxis()->SetTitleOffset(1.3);
//    hpx->SetStats(kFALSE);   // no statistics
//    hpx->Draw();
//
//
//    khodyukGraph.Draw("P E1 same");
//    moszynskiGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//
//    int polynomialDegree=1; //or num parameters-1
//    const int numRegions=4;
//    TF1* regionFuncs[numRegions];
//    double regionEnergyDivisions[numRegions+1]={9,13.4,16.5,90,400};
//
//    int numRepeatFits=8;
//    TString fitString;
//
//    TString filePath=HistsDir+"EnergyResponseOutputParametersForCumulativeFit130617.txt";
//    ofstream outFile(filePath);
//    outFile << "EnergyLowLimit(keV)" << "\t";
//    outFile << "EnergyHighLimit(keV)" << "\t";
//    for (int i = 0;i< polynomialDegree+1;i++)
//    {
//        outFile << "Parameter" << i+1 << "\t" << "Parameter" << i+1 << "Error\t";
//    }
//    outFile << endl;
//
//
//
//    TGraphErrors* theGraphs[numGraphs]={&khodyukGraph,&moszynskiGraph,&dewarGraph,&calRDK2Graph};
//    double graphSigmaFactors[numGraphs]={2.,2.,2.,1.};
//    for (int i = 0;i< numRegions;i++)
//    {
//        ///Create graph for region
//
//        vector<double> regionEnergies;
//        vector<double> regionLOs;
//        vector<double> regionErrors;
//        for (int j = 0;j< numGraphs;j++)
//        {
//            int expPointsInRegion=0;
//            double* currentEnergies=theGraphs[j]->GetX();
//            double* currentLOs=theGraphs[j]->GetY();
//            double* currentLOErrors=theGraphs[j]->GetEY();
//            for (int k = 0;k< theGraphs[j]->GetN();k++)
//            {
//
//            	if(currentEnergies[k] >= regionEnergyDivisions[i] && currentEnergies[k] <= regionEnergyDivisions[i+1])
//            	{
//            	    regionEnergies.push_back(currentEnergies[k]);
//            	    regionLOs.push_back(currentLOs[k]);
//                    expPointsInRegion++;
//
//            	}
//            }
//
//            for (int k = 0;k< theGraphs[j]->GetN();k++)
//            {
//                if(currentEnergies[k] >= regionEnergyDivisions[i] && currentEnergies[k] <= regionEnergyDivisions[i+1])
//            	{
//                    regionErrors.push_back(currentLOErrors[k]*sqrt(expPointsInRegion)*graphSigmaFactors[j]);
//                }
//            }
//
//
//        }
//        TGraphErrors regionGraph(regionEnergies.size(),regionEnergies,regionLOs,NULL,regionErrors);
//        TString regionName="region"+int2str(i)+"Func";
////        TString regionFunctionString="pol"+int2str(polynomialDegree);
//        TString regionFunctionString;
//        if(i==numRegions-1)
//        {
//            regionFunctionString="[0]/(1.+(662.*([0]-1.))/x)";
//        }
//        else
//        {
//            regionFunctionString="[0]/(1.+[1]/x)";
//        }
//        regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i],regionEnergyDivisions[i+1]);
//        regionFuncs[i]->SetLineColor(kBlack);
//        regionFuncs[i]->SetParameter(0,100);
//        if(i!=numRegions-1)
//        {
//            regionFuncs[i]->SetParameter(1,16);
//        }
////        if(i==0)
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[0],regionEnergyDivisions[i+2]);
////        }
////        else if(i==numRegions-1)
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i-1],regionEnergyDivisions[i+1]);
////        }
////        else
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i-1],regionEnergyDivisions[i+2]);
////        }
//    	for (int j = 0;j< numRepeatFits-1;j++)
//    	{
//            fitString="VMRN";
//            regionGraph.Fit(regionFuncs[i],fitString);
//    	}
//
//        fitString="VMRN+";
//    	regionGraph.Fit(regionFuncs[i],fitString);
//    	regionFuncs[i]->SetRange(regionEnergyDivisions[i],regionEnergyDivisions[i+1]);
//    	regionFuncs[i]->Draw("SAME");
//    	outFile << regionEnergyDivisions[i] << "\t";
//    	outFile << regionEnergyDivisions[i+1] << "\t";
//    	for (int j = 0;j< polynomialDegree+1;j++)
//    	{
//            outFile << regionFuncs[i]->GetParameter(j) << "\t" << regionFuncs[i]->GetParError(j) << "\t" ;
//    	}
//    	outFile << endl;
//    }
//
//
//
//    regionFuncs[0]->SetRange(0,regionEnergyDivisions[1]);
//   regionFuncs[numRegions-1]->SetRange(regionEnergyDivisions[numRegions-1],800);
//
//   TF1* regionFuncsLow[numRegions];
//   TF1* regionFuncsHigh[numRegions];
//   double lowFactor=.97;
//   double highFactor=1.03;
//   for (int i = 0;i< numRegions-1;i++)
//   {
//       TString newName=string(regionFuncs[i]->GetName())+"_Low";
//       regionFuncsLow[i]=(TF1*) regionFuncs[i]->Clone(newName);
//       regionFuncsLow[i]->SetParameter(0,regionFuncs[i]->GetParameter(0)*lowFactor);
//       regionFuncsLow[i]->Draw("SAME");
//       regionFuncsLow[i]->SetLineColor(kBlue);
//
//       newName=string(regionFuncs[i]->GetName())+"_High";
//       regionFuncsHigh[i]=(TF1*) regionFuncs[i]->Clone(newName);
//       regionFuncsHigh[i]->SetParameter(0,regionFuncs[i]->GetParameter(0)*highFactor);
//       regionFuncsHigh[i]->Draw("SAME");
//       regionFuncsHigh[i]->SetLineColor(kBlue);
//   }
//
//   TString newName=string(regionFuncs[numRegions-1]->GetName())+"_Low";
//   regionFuncsLow[numRegions-1]=(TF1*) regionFuncs[numRegions-1]->Clone(newName);
//   double LOOfCentered=regionFuncs[numRegions-1]->Eval(regionEnergyDivisions[numRegions-1]);
//   double newParam=(regionEnergyDivisions[numRegions-1]-662.)/(regionEnergyDivisions[numRegions-1]/(LOOfCentered*lowFactor)-662.);
//   cout << "Alow: " << newParam << endl;
//   regionFuncsLow[numRegions-1]->SetParameter(0,newParam);
//   regionFuncsLow[numRegions-1]->Draw("SAME");
//   regionFuncsLow[numRegions-1]->SetLineColor(kBlue);
//
//   newName=string(regionFuncs[numRegions-1]->GetName())+"_High";
//   regionFuncsHigh[numRegions-1]=(TF1*) regionFuncs[numRegions-1]->Clone(newName);
//   newParam=(regionEnergyDivisions[numRegions-1]-662.)/(regionEnergyDivisions[numRegions-1]/(LOOfCentered*highFactor)-662.);
//   cout << "Ahigh: " << newParam << endl;
//   regionFuncsHigh[numRegions-1]->SetParameter(0,newParam);
//   regionFuncsHigh[numRegions-1]->Draw("SAME");
//   regionFuncsHigh[numRegions-1]->SetLineColor(kBlue);
//
//    ///connect points at region intersections visually
//
//    for (int i = 0;i< numRegions-1;i++)
//    {
//        double x1,x2,y1,y2;
//
//        //center
//        x1=x2=regionEnergyDivisions[i+1];
//        y1=regionFuncs[i]->Eval(x1);
//        y2=regionFuncs[i+1]->Eval(x1);
//        TLine* lineToDraw=new TLine(x1,y1,x2,y2);
//        lineToDraw->SetLineColor(kBlack);
//        lineToDraw->Draw("SAME");
//
//        //Low
//        x1=x2=regionEnergyDivisions[i+1];
//        y1=regionFuncsLow[i]->Eval(x1);
//        y2=regionFuncsLow[i+1]->Eval(x1);
//        lineToDraw=new TLine(x1,y1,x2,y2);
//        lineToDraw->SetLineColor(kBlue);
//        lineToDraw->Draw("SAME");
//
//        //High
//        x1=x2=regionEnergyDivisions[i+1];
//        y1=regionFuncsHigh[i]->Eval(x1);
//        y2=regionFuncsHigh[i+1]->Eval(x1);
//        lineToDraw=new TLine(x1,y1,x2,y2);
//        lineToDraw->SetLineColor(kBlue);
//        lineToDraw->Draw("SAME");
//
//    }
//
//    //Add old Matt Param model
//    const int numOldPoints=1000;
//    double xOld[numOldPoints];
//    double yOld[numOldPoints];
//    for (int i = 0;i< 1000;i++)
//    {
//        xOld[i]=i;
//        yOld[i]=evalMattKhodyukParam(xOld[i]);
//    }
//    TGraph oldGraph(numOldPoints,xOld,yOld);
//    oldGraph.SetLineColor(kPink);
//    oldGraph.Draw("LSame");
//
//    theLegend->AddEntry("region0Func","Fit","L");
//    theLegend->AddEntry("region0Func_Low","#pm 3% of fit","L");
//    theLegend->AddEntry(&oldGraph,"oldKhodyuk Model","L");
//
//
//    theLegend->Draw("same");
//
//
//
//    TString imagePath2=GraphsDir+"EnergyResponseOverallFit130626.png";
//    c2->SaveAs(imagePath2);
//    delete c2;
//
//
//    for (int i = 0;i< numRegions;i++)
//    {
//        delete regionFuncs[i];
//
//    }
//
//    outFile.close();
//
//    delete theLegend;

///make graphs of new all exp LO param model 130627
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD NonRad W Fermi");
////    fourBodySet->makeDerivedResultsFiles(GAMMA);
//
////    TH1D* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut&&"Entry$ % 4 == 0",BGO_DIM);
////    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut&&"Entry$ % 4 == 1",BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130627_EDepGWLO_AverageIndiv_EPG_RID175.txt");
//
//    TH1D* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",STD_EPCut&&"Entry$ % 4 == 2",BGO_DIM);
//    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist130627_EDepGWLOLow_AverageIndiv_EPG_RID175.txt");
//
//    TH1D* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",STD_EPCut&&"Entry$ % 4== 3",BGO_DIM);
//    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist130627_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt");
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->Se tTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp Param Model Comparison","ParamModelComparison130627.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison130627_Residuals.png");

///Check that field class changes worked in mrkgeant
////    RDK2Set* oldThreeBodySet=new RDK2Set(1,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
////    RDK2Set* oldTourBodySet=new RDK2Set(1,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD NonRad W Fermi");
//
//    RDK2Set* threeBodySet=new RDK2Set(1,"130429","STD_3B","RID209","RID210","","STD NonRad W Fermi newFieldClass");
//    RDK2Set* fourBodySet=new RDK2Set(1,"130429","STD_4B_100eV_Endpoint","RID213","RID214","RID215","STD NonRad W Fermi newFieldClass");
//
//    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut,BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID215.txt");
//
//    TH1* mcHists[2] = {getTabSeperatedHist(HistsDir+"MC_Hist130617_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM),mcHist};
//
//    mcHists[0]->SetTitle("Old field method;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New field method;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"compareNewFieldMethod130627","compareNewFieldMethod130627.png");

///Compare geant4 physics libraries

//    TString physicsList="emlivermore empenelope emstandard_opt0 emstandard_opt1 emstandard_opt2 emstandard_opt3 emstandard_opt4 standardSS standardSSM standardWVI standardGS";
//
//    int RIDNumber=216;
//    const int numHists=numItemsInStringList(physicsList);
//    TH1* mcHists[numHists];
//    for (int i = 0;i< numItemsInStringList(physicsList);i++)
//    {
//        TString thePhysics;
//        stringFromList(physicsList,i,thePhysics);
//        TString RIDString="RID"+int2str(RIDNumber+i);
//
//        RDK2Set fourBodySet(1,"130429","STD_4B_100eV_Endpoint","RID213","RID214",RIDString,"special");
//        mcHists[i]=makeMCIndividualDetectorHist(NULL,&fourBodySet,"BGOEDepBlurEdge","",BGO_DIM);
//        convertTH1ToTXT(mcHists[i],HistsDir+"MC_Hist130627_DirectEDepG_AverageIndiv_EPG_"+RIDString+".txt");
//        TString titleString=thePhysics+";Energy Detected (keV);Counts Detected (AU)";
//        mcHists[i]->SetTitle(titleString);
//        mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,numHists,mcHists,"Geant4.9.6.p02 Physics Libraries","PhysicsLibraryTest_BGO_Geant4.9.6.p02.png");
//    plotResidualsToImage(numHists-1,&(mcHists[1]),1,mcHists,"Geant4.9.6.p02 Physics Libraries #frac{X - emlivermore}{emlivermore}","PhysicsLibraryTest_BGO_Geant4.9.6.p02_Residuals.png");

///make graphs of new all exp LO param model 130709
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD NonRad W Fermi");
//////    fourBodySet->makeDerivedResultsFiles(GAMMA);
//////
//    TH1* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut&&"Entry$ % 4 == 0","",BGO_DIM);
//    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist130627b_EDepG_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut&&"Entry$ % 4 == 1",BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130709_EDepGWLO_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130709_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1D* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",STD_EPCut&&"Entry$ % 4 == 2",BGO_DIM);
//    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist130709_EDepGWLOLow_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130709_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1D* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",STD_EPCut&&"Entry$ % 4== 3",BGO_DIM);
//    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist130709_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130709_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130709.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp. Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison130709_Residuals.png");
//
//      TString drawString;
//
//    TH1* expHist[1];
//
//
//    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
//    expHist[0]->Sumw2();
//
//    for (int i = 0;i< 12;i++)
//    {
//        TString cutString="SBDToFP < 10";
//        TCut EPCut=TCut(cutString);
//
//        drawString="SPECIAL_EDepG";
//        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, EPCut, -(i+1),BGO_DIM);
//        expHist[0]->Add(tempHist);
//        delete tempHist;
//    }
//    expHist[0]->Scale(1/12.);
//
//    TString histPath=HistsDir+"ExpMatt_Hists130709_G7_BGOEDepAvg_EPG.png";
//    convertTH1ToTXT(expHist[0],histPath);
////  expEPGHists[j]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130709.png");
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison #frac{EXP-MC}{MC}","ParamModelComparison130709_ResidualsVsExp.png");

///Create central light response curve with correct region by region experiment weighting with Averkiev
//    const int numGraphs=5;
//
//    const int numDewarPoints=14;
//    double dewarEnergies[numDewarPoints]={122.1,81,59.5,31,99.8,88,23.3,122.1,99.8,88,23.3,122.1,185.7,107};
//    double dewarLOs[numDewarPoints]={0.90033,0.91603,0.88772,0.81225,0.86602,0.90747,0.73221,0.8827,0.86363,0.90637,0.7364,0.89272,0.93333,0.8203};
//    double dewarLOErrors[numDewarPoints]={0.010849,0.010122,0.010078,0.010542,0.012021,0.010652,0.03,0.010572,0.010722,0.010078,0.025,0.010022,0.015,0.02};
//    TGraphErrors dewarGraph(numDewarPoints,dewarEnergies,dewarLOs,NULL,dewarLOErrors);
//    dewarGraph.SetName("DewarGraph");
//    dewarGraph.SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    TString dewarLegendString="RDK2 BGO in Dewar [MCA] (Various Sources)";
//    dewarGraph.SetLineColor(kRed);
//    dewarGraph.SetMarkerColor(kRed);
//    dewarGraph.SetMarkerStyle(kFullCircle);
//    dewarGraph.SetMarkerSize(1.5);
//
//
//    const int numCalRDK2Points=6;
//    double calRDK2Energies[numCalRDK2Points]={32.3,42.7,78.7,79.7,121.9,356};
//    double calRDK2LOs[numCalRDK2Points]={0.81,0.86,0.94,0.93,0.92,0.986};
//    double calRDK2LOErrors[numCalRDK2Points]={0.02,0.02,0.01,0.01,0.01,0.002};
//    TGraphErrors calRDK2Graph(numCalRDK2Points,calRDK2Energies,calRDK2LOs,NULL,calRDK2LOErrors);
//    calRDK2Graph.SetName("CalRDK2Graph");
//    TString calRDK2LegendString="RDK2 BGO Array [DAQ] (Matt's Cal & Jeff's 80keV)";
//    calRDK2Graph.SetLineColor(kViolet);
//    calRDK2Graph.SetMarkerColor(kViolet);
//    calRDK2Graph.SetMarkerStyle(kFullSquare);
//    calRDK2Graph.SetMarkerSize(1.5);
//
//    //TGraph* khodyukGraph=getKhodyukBGOResponseGraphOfficial();
//    //Assumed error was similar to moszynski and was 2%
//    const int numKhodyukPoints=375;
//    double khodyukEnergies[numKhodyukPoints]={11,11,12,12,13,13.025,13.05,13.075,13.1,13.125,13.15,13.175,13.2,13.225,13.25,13.275,13.3,13.325,13.35,13.375,
//        13.4,13.425,13.45,13.475,13.5,13.525,13.55,13.575,13.6,13.625,13.65,13.675,13.7,13.725,13.75,13.775,13.8,13.825,13.85,13.875,13.9,13.925,13.95,13.975,14,
//        14.025,14.05,14.075,14.1,14.125,14.15,14.175,14.2,14.225,14.25,14.275,14.3,14.325,14.35,14.375,14.4,14.425,14.45,14.475,14.5,14.525,14.55,14.575,14.6,14.625,
//        14.65,14.675,14.7,14.725,14.75,14.775,14.8,14.825,14.85,14.875,14.9,14.925,14.95,14.975,15,15.025,15.05,15.075,15.1,15.125,15.15,15.175,15.2,15.225,15.25,
//        15.275,15.3,15.325,15.35,15.375,15.4,15.425,15.45,15.475,15.5,15.525,15.55,15.575,15.6,15.625,15.65,15.675,15.7,15.725,15.75,15.775,15.8,15.825,15.85,15.875,
//        15.9,15.925,15.95,15.975,16,16.025,16.05,16.075,16.1,16.125,16.15,16.175,16.2,16.225,16.25,16.275,16.3,16.325,16.35,16.375,16.4,16.425,16.45,16.475,16.5,16.525,
//        16.55,16.575,16.6,16.625,16.65,16.675,16.7,16.725,16.75,16.775,16.8,16.825,16.85,16.875,16.9,16.925,16.95,16.975,17,17.025,17.05,17.075,17.1,17.125,17.15,17.175,
//        17.2,17.225,17.25,17.275,17.3,17.325,17.35,17.375,17.4,17.425,17.45,17.475,17.5,17.525,17.55,17.575,17.6,17.625,17.65,17.675,17.7,17.725,17.75,17.775,17.8,17.825,
//        17.85,17.875,17.9,17.925,17.95,17.975,18,18.025,18.05,18.075,18.1,18.125,18.15,18.175,18.2,18.225,18.25,18.275,18.3,18.325,18.35,18.375,18.4,18.425,18.45,18.475,
//        18.5,18.525,18.55,18.575,18.6,18.625,18.65,18.675,18.7,18.725,18.75,18.775,18.8,18.825,18.85,18.875,18.9,18.925,18.95,18.975,19,19.025,19.05,19.075,19.1,19.125,
//        19.15,19.175,19.2,19.225,19.25,19.275,19.3,19.325,19.35,19.375,19.4,19.425,19.45,19.475,19.5,19.525,19.55,19.575,19.6,19.625,19.65,19.675,19.7,19.725,19.75,19.775,
//        19.8,19.825,19.85,19.875,19.9,19.925,19.95,19.975,20,20,25,25,30,30,30,30,35,35,40,40,45,45,45,45,50,50,55,55,60,60,65,65,70,70,75,75,80,80,85,85,89.5,89.6,89.7,
//        89.8,89.9,90,90.1,90.2,90.3,90.4,90.5,90.6,90.7,90.8,90.9,91,91.1,91.2,91.3,91.4,91.5,91.6,91.7,91.8,91.9,92,92.1,92.2,92.3,92.4,92.5,92.6,92.7,92.8,92.9,93,93.1,
//        93.2,93.3,93.4,93.5,93.6,93.7,93.8,93.9,94,94.1,94.2,94.3,94.4,94.5,94.6,94.7,94.8,94.9,95,95,100,100};
//    double khodyukLOs[numKhodyukPoints]={0.700,0.707,0.725,0.724,0.740,0.742,0.737,0.743,0.743,0.741,0.737,0.740,0.738,0.736,0.743,0.746,0.743,0.739,0.742,0.742,
//        0.736,0.722,0.703,0.704,0.705,0.695,0.703,0.701,0.698,0.699,0.702,0.702,0.704,0.695,0.701,0.696,0.705,0.696,0.706,0.701,0.696,0.693,0.697,0.697,0.694,0.700,0.699,0.698,
//        0.699,0.702,0.697,0.692,0.693,0.696,0.697,0.689,0.697,0.683,0.684,0.693,0.687,0.698,0.695,0.702,0.699,0.697,0.697,0.699,0.697,0.699,0.694,0.691,0.699,0.691,0.692,0.694,
//        0.694,0.702,0.697,0.699,0.701,0.695,0.698,0.699,0.702,0.698,0.697,0.700,0.702,0.698,0.697,0.698,0.701,0.708,0.703,0.697,0.699,0.698,0.701,0.705,0.704,0.706,0.702,0.700,
//        0.704,0.704,0.706,0.712,0.701,0.705,0.704,0.704,0.702,0.698,0.697,0.702,0.703,0.698,0.699,0.693,0.701,0.700,0.701,0.700,0.694,0.697,0.698,0.698,0.699,0.695,0.701,0.704,
//        0.698,0.701,0.697,0.701,0.697,0.700,0.694,0.697,0.691,0.687,0.684,0.690,0.692,0.685,0.689,0.687,0.693,0.686,0.694,0.697,0.690,0.691,0.690,0.688,0.695,0.690,0.690,0.690,
//        0.690,0.694,0.693,0.693,0.690,0.688,0.692,0.687,0.689,0.692,0.687,0.694,0.695,0.692,0.693,0.693,0.698,0.694,0.695,0.692,0.699,0.694,0.693,0.695,0.695,0.696,0.692,0.694,
//        0.690,0.694,0.697,0.698,0.700,0.700,0.700,0.696,0.701,0.696,0.696,0.698,0.704,0.705,0.696,0.694,0.698,0.701,0.702,0.703,0.697,0.701,0.702,0.701,0.698,0.708,0.702,0.702,
//        0.701,0.702,0.699,0.706,0.705,0.700,0.714,0.699,0.710,0.705,0.704,0.702,0.704,0.703,0.703,0.703,0.705,0.706,0.703,0.701,0.712,0.701,0.710,0.707,0.709,0.704,0.709,0.713,
//        0.706,0.711,0.710,0.709,0.719,0.709,0.705,0.711,0.712,0.711,0.708,0.709,0.709,0.710,0.710,0.716,0.716,0.712,0.709,0.717,0.713,0.715,0.711,0.716,0.716,0.717,0.714,0.715,
//        0.711,0.720,0.714,0.718,0.715,0.711,0.720,0.716,0.717,0.722,0.722,0.724,0.720,0.718,0.768,0.770,0.803,0.800,0.800,0.800,0.827,0.826,0.844,0.845,0.861,0.864,0.858,0.858,
//        0.875,0.874,0.880,0.880,0.885,0.885,0.890,0.889,0.896,0.893,0.897,0.898,0.901,0.895,0.906,0.909,0.911,0.911,0.912,0.910,0.911,0.910,0.909,0.910,0.908,0.908,0.899,0.885,
//        0.885,0.884,0.885,0.881,0.883,0.880,0.882,0.881,0.879,0.879,0.878,0.881,0.879,0.879,0.880,0.880,0.882,0.882,0.877,0.878,0.878,0.880,0.876,0.879,0.880,0.881,0.879,0.878,
//        0.878,0.879,0.878,0.880,0.880,0.879,0.878,0.878,0.877,0.878,0.880,0.877,0.878,0.878,0.878,0.881,0.880,0.882,0.883};
//    double khodyukLOErrors[numKhodyukPoints]={0.014,0.014,0.015,0.014,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,
//        0.015,0.015,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,
//        0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.015,0.015,0.016,0.016,0.016,0.016,0.017,0.017,0.017,0.017,0.017,0.017,0.017,
//        0.017,0.018,0.017,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
//        0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,
//        0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018};
//    TGraphErrors khodyukGraph(numKhodyukPoints,khodyukEnergies,khodyukLOs,NULL,khodyukLOErrors);
//    TString khodyukLegendString="I. Khodyuk and P. Dorenbos (monochromator)";
//    khodyukGraph.SetTitle("Light Output from Experimental Sources; Energy Incident/Absorbed (keV); Light Output Relative to 662 keV");
//    khodyukGraph.SetLineColor(kGreen);
//    khodyukGraph.SetMarkerColor(kGreen);
//    khodyukGraph.SetMarkerStyle(kFullCircle);
//
//
//    const int numMoszynskiPoints=30;
//    double moszynskiEnergies[numMoszynskiPoints]={14.4,16.6,22,31,32,45,52,59.6,74,81,84,88,122,244.1,302.9,320.1,335.8,356,411.8,435,493.7,511,569.7,586,662,758.8,834.8,1063.7,1173.2,1275};
//    double moszynskiLOs[numMoszynskiPoints]={0.711,0.719,0.775,0.829,0.8412,0.853,0.9,0.893,0.91,0.926,0.937,0.933,0.915,0.964,0.979,0.964,0.97,0.986,0.97,0.995,1,0.996,1,0.999,1,1,1,1.01,1.013,1.01};
//    double moszynskiLOErrors[numMoszynskiPoints]={0.01422,0.01438,0.0155,0.01658,0.016824,0.02559,0.018,0.01786,0.0182,0.01852,0.01874,0.01866,0.0183,0.02892,0.01958,0.01928,0.0291,0.01972,0.0194,0.02985,0.03,0.01992,0.02,0.02997,0.02,0.03,0.02,0.0202,0.02026,0.0202};
//    TGraphErrors moszynskiGraph(numMoszynskiPoints,moszynskiEnergies,moszynskiLOs,NULL,moszynskiLOErrors);
//    moszynskiGraph.SetName("MoszynskiGraph");
//    TString moszynskiLegendString="Moszynski et. al.";
//    moszynskiGraph.SetLineColor(kCyan);
//    moszynskiGraph.SetMarkerColor(kCyan);
//    moszynskiGraph.SetMarkerStyle(kFullSquare);
//    moszynskiGraph.SetMarkerSize(1.5);
//
//    const int numAverkievPoints=20;
//    double averkievEnergies[numAverkievPoints]={2.1,6.4,7.5,8.0,8.8,10.0,11.2,16.6,17.2,23.5,33.1,58.3,71.6,78.3,120.1,162.3,277.4,384.8,506.5,657.9};
//    double averkievLOs[numAverkievPoints]={0.854,0.777,0.739,0.706,0.678,0.657,0.641,0.579,0.590,0.685,0.731,0.837,0.844,0.827,0.884,0.903,0.959,0.977,0.986,0.996};
//    double averkievLOErrors[numAverkievPoints]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.035,0.035,0.035,0.035,0.035,0.02,0.015,0.01,0.005,0.005};
//    TGraphErrors averkievGraph(numAverkievPoints,averkievEnergies,averkievLOs,NULL,averkievLOErrors);
//    averkievGraph.SetName("AverkievGraph");
//    TString averkievLegendString="Averkiev et. al.";
//    averkievGraph.SetLineColor(kMagenta);
//    averkievGraph.SetMarkerColor(kMagenta);
//    averkievGraph.SetMarkerStyle(kOpenCircle);
//    averkievGraph.SetMarkerSize(1.5);
//
//
//
////    int detectorNumber=1;
//////    const int numQuadraticPoints=19;
//////    double quadraticEnergies[numQuadraticPoints]={10,20,30,40,50,60,70,80,100,150,200,300,400,500,600,700,800,900,1000};
//////    double quadraticLOs[numQuadraticPoints];
////
////    for (int i=0;i<numQuadraticPoints;i++)
////    {
////        quadraticLOs[i]=applyEDepModel(quadraticEnergies[i],detectorNumber, NULL, KEVIN_GAUS_EDEP)/quadraticEnergies[i];
////    }
////    TGraph quadraticGraph(numQuadraticPoints,quadraticEnergies,quadraticLOs);
////    quadraticGraph.SetName("QuadraticGraph");
////    TString quadraticLegendString="Old Quadratic Model";
////    quadraticGraph.SetLineColor(kGreen);
////    quadraticGraph.SetMarkerColor(kGreen);
////    quadraticGraph.SetMarkerStyle(kFullSquare);
////    quadraticGraph.SetMarkerSize(1.5);
//
//
//
//
//    double legendSize=.05*numGraphs;
//    TLegend* theLegend = new TLegend(0.6,.4-legendSize,0.99,0.4);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
////    theLegend->AddEntry(&quadraticGraph,quadraticLegendString,"L");
//    theLegend->AddEntry(&khodyukGraph,khodyukLegendString,"P");
//    theLegend->AddEntry(&moszynskiGraph,moszynskiLegendString,"P");
//    theLegend->AddEntry(&averkievGraph,averkievLegendString,"P");
//    theLegend->AddEntry(&dewarGraph,dewarLegendString,"P");
//    theLegend->AddEntry(&calRDK2Graph,calRDK2LegendString,"P");
//
//
//    TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1280,1024);
//    theCanvas.SetFillColor(kGray);
//    theCanvas.cd();
//    theCanvas.SetLeftMargin(.12);
//    theCanvas.SetRightMargin(.1);
//    theCanvas.SetBottomMargin(.1);
//    theCanvas.SetGridx();
//    theCanvas.SetGridy();
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//
//    khodyukGraph.GetYaxis()->SetTitleOffset(1.3);
//    khodyukGraph.GetXaxis()->SetLimits(0,150.);
//    khodyukGraph.GetYaxis()->SetRangeUser(.55,1.);
//
//
//    khodyukGraph.Draw("A P E1");
////    quadraticGraph.Draw("L same");
//    moszynskiGraph.Draw("P E1 same");
//    averkievGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//    theLegend->Draw("same");
//
//    TString imagePath=string(GRAPHS_DIR)+ "LOExperimentsGraph_130710.png";
//
//    theCanvas.SaveAs(imagePath);
//
//    ///create master graph that has all the data in it with errors scaled to correspond with confidence and number of data points for fitting;
//
//
//    const int numCumulativePoints=numKhodyukPoints+numMoszynskiPoints+numDewarPoints+numCalRDK2Points+numAverkievPoints;
//    double cumulativeEnergies[numCumulativePoints];
//    double cumulativeLOs[numCumulativePoints];
//    double cumulativeLOErrors[numCumulativePoints];
//
//    double khodyukSigmaFactor=2;
//    double moszynskiSigmaFactor=2;
//    double dewarSigmaFactor=2;
//    double calRDK2SigmaFactor=1;
//    double averkievSigmaFactor=4;
//
//    //Add Khodyuk graph
//    for (int i = 0;i< numKhodyukPoints;i++)
//    {
//        cumulativeEnergies[i]=khodyukEnergies[i];
//        cumulativeLOs[i]=khodyukLOs[i];
//        cumulativeLOErrors[i]=khodyukLOErrors[i]*sqrt(numKhodyukPoints)*khodyukSigmaFactor;
//    }
//
//    //Add moszynski graph
//    for (int i = 0;i< numMoszynskiPoints;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints]=moszynskiEnergies[i];
//        cumulativeLOs[i+numKhodyukPoints]=moszynskiLOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints]=moszynskiLOErrors[i]*sqrt(numMoszynskiPoints)*moszynskiSigmaFactor;
//    }
//
//    //Add dewar graph
//    for (int i = 0;i< numDewarPoints;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints+numMoszynskiPoints]=dewarEnergies[i];
//        cumulativeLOs[i+numKhodyukPoints+numMoszynskiPoints]=dewarLOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints+numMoszynskiPoints]=dewarLOErrors[i]*sqrt(numDewarPoints)*dewarSigmaFactor;
//    }
//
//    //Add RDK2 apparatus graph
//    for (int i = 0;i< numCalRDK2Points;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2Energies[i];
//        cumulativeLOs[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2LOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints]=calRDK2LOErrors[i]*sqrt(numCalRDK2Points)*calRDK2SigmaFactor;
//    }
//
//    //Add Averkiev graph
//    for (int i = 0;i< numAverkievPoints;i++)
//    {
//        cumulativeEnergies[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints+numCalRDK2Points]=averkievEnergies[i];
//        cumulativeLOs[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints+numCalRDK2Points]=averkievLOs[i];
//        cumulativeLOErrors[i+numKhodyukPoints+numMoszynskiPoints+numDewarPoints+numCalRDK2Points]=averkievLOErrors[i]*sqrt(numAverkievPoints)*averkievSigmaFactor;
//    }
//
//
//    TGraphErrors cumulativeGraph(numCumulativePoints,cumulativeEnergies,cumulativeLOs,NULL,cumulativeLOErrors);
//
//    TCanvas* c2 = new TCanvas("fitFunc","fitFunc",10,10,1280,1024);
//    TH2F* hpx = new TH2F("hpx","Fit and Extrapolation to Energy Response Data;Energy Incident (keV);Light Output Relative to 662 keV",10,0,662,10,.3,1.05); // axis range
//    hpx->GetYaxis()->SetTitleOffset(1.3);
//    hpx->SetStats(kFALSE);   // no statistics
//    hpx->Draw();
//
//
//    khodyukGraph.Draw("P E1 same");
//    moszynskiGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//    averkievGraph.Draw("P E1 same");
//
//    int polynomialDegree=2; //or num parameters-1
//    const int numRegions=5;
//    TF1* regionFuncs[numRegions];
//    double regionEnergyDivisions[numRegions+1]={1.9,10,13.4,16.5,90,400};
//
//    int numRepeatFits=8;
//    TString fitString;
//
//    TString filePath=HistsDir+"EnergyResponseOutputParametersForCumulativeFit130710.txt";
//    ofstream outFile(filePath);
//    outFile << "EnergyLowLimit(keV)" << "\t";
//    outFile << "EnergyHighLimit(keV)" << "\t";
//    for (int i = 0;i< polynomialDegree+1;i++)
//    {
//        outFile << "Parameter" << i+1 << "\t" << "Parameter" << i+1 << "Error\t";
//    }
//    outFile << endl;
//
//
//
//    TGraphErrors* theGraphs[numGraphs]={&khodyukGraph,&moszynskiGraph,&dewarGraph,&calRDK2Graph,&averkievGraph};
//    double graphSigmaFactors[numGraphs]={2.,2.,2.,1.,2};
//    for (int i = 0;i< numRegions;i++)
//    {
//        ///Create graph for region
//
//        vector<double> regionEnergies;
//        vector<double> regionLOs;
//        vector<double> regionErrors;
//        for (int j = 0;j< numGraphs;j++)
//        {
//            int expPointsInRegion=0;
//            double* currentEnergies=theGraphs[j]->GetX();
//            double* currentLOs=theGraphs[j]->GetY();
//            double* currentLOErrors=theGraphs[j]->GetEY();
//            for (int k = 0;k< theGraphs[j]->GetN();k++)
//            {
//
//            	if(currentEnergies[k] >= regionEnergyDivisions[i] && currentEnergies[k] <= regionEnergyDivisions[i+1])
//            	{
//            	    regionEnergies.push_back(currentEnergies[k]);
//            	    regionLOs.push_back(currentLOs[k]);
//                    expPointsInRegion++;
//
//            	}
//            }
//
//            for (int k = 0;k< theGraphs[j]->GetN();k++)
//            {
//                if(currentEnergies[k] >= regionEnergyDivisions[i] && currentEnergies[k] <= regionEnergyDivisions[i+1])
//            	{
//                    regionErrors.push_back(currentLOErrors[k]*sqrt(expPointsInRegion)*graphSigmaFactors[j]);
//                }
//            }
//
//
//        }
//        TGraphErrors regionGraph(regionEnergies.size(),regionEnergies,regionLOs,NULL,regionErrors);
//        TString regionName="region"+int2str(i)+"Func";
////        TString regionFunctionString="pol"+int2str(polynomialDegree);
//        TString regionFunctionString;
//        if(i==numRegions-1)
//        {
//            regionFunctionString="[0]/(1.+(662.*([0]-1.))/x)";
//        }
//        else if(i==0)
//        {
//            regionFunctionString="pol"+int2str(polynomialDegree);
//        }
//        else
//        {
//            regionFunctionString="[0]/(1.+[1]/x)";
//        }
//        regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i],regionEnergyDivisions[i+1]);
//        regionFuncs[i]->SetLineColor(kBlack);
//        regionFuncs[i]->SetParameter(0,100);
//        if(i!=numRegions-1)
//        {
//            regionFuncs[i]->SetParameter(1,16);
//        }
////        if(i==0)
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[0],regionEnergyDivisions[i+2]);
////        }
////        else if(i==numRegions-1)
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i-1],regionEnergyDivisions[i+1]);
////        }
////        else
////        {
////            regionFuncs[i]=new TF1(regionName,regionFunctionString,regionEnergyDivisions[i-1],regionEnergyDivisions[i+2]);
////        }
//    	for (int j = 0;j< numRepeatFits-1;j++)
//    	{
//            fitString="VMRN";
//            regionGraph.Fit(regionFuncs[i],fitString);
//    	}
//
//        fitString="VMRN+";
//    	regionGraph.Fit(regionFuncs[i],fitString);
//    	regionFuncs[i]->SetRange(regionEnergyDivisions[i],regionEnergyDivisions[i+1]);
//    	regionFuncs[i]->Draw("SAME");
//    	outFile << regionEnergyDivisions[i] << "\t";
//    	outFile << regionEnergyDivisions[i+1] << "\t";
//    	for (int j = 0;j< polynomialDegree+1;j++)
//    	{
//            outFile << regionFuncs[i]->GetParameter(j) << "\t" << regionFuncs[i]->GetParError(j) << "\t" ;
//    	}
//
//    	outFile << endl;
//    }
//
//
//
//   regionFuncs[0]->SetRange(0,regionEnergyDivisions[1]);
//   regionFuncs[numRegions-1]->SetRange(regionEnergyDivisions[numRegions-1],800);
//
//   TF1* regionFuncsLow[numRegions];
//   TF1* regionFuncsHigh[numRegions];
//   double lowFactor=.97;
//   double highFactor=1.03;
//   for (int i = 0;i< numRegions-1;i++)
//   {
//       TString newName=string(regionFuncs[i]->GetName())+"_Low";
//
//       if(i==0)
//       {
//            regionFuncsLow[0]=(TF1*) regionFuncs[1]->Clone(newName);
//            regionFuncsLow[0]->SetParameter(0,regionFuncs[1]->GetParameter(0)*lowFactor);
//            regionFuncsLow[0]->SetParameter(1,regionFuncs[1]->GetParameter(1));
//            regionFuncsLow[0]->SetRange(0,regionEnergyDivisions[1]);
//       }
//       else
//       {
//           regionFuncsLow[i]=(TF1*) regionFuncs[i]->Clone(newName);
//           regionFuncsLow[i]->SetParameter(0,regionFuncs[i]->GetParameter(0)*lowFactor);
//       }
//       regionFuncsLow[i]->Draw("SAME");
//       regionFuncsLow[i]->SetLineColor(kBlue);
//
//       newName=string(regionFuncs[i]->GetName())+"_High";
//       regionFuncsHigh[i]=(TF1*) regionFuncs[i]->Clone(newName);
//       regionFuncsHigh[i]->SetParameter(0,regionFuncs[i]->GetParameter(0)*highFactor);
//       regionFuncsHigh[i]->Draw("SAME");
//       regionFuncsHigh[i]->SetLineColor(kBlue);
//   }
//
//
//   regionFuncs[0]->SetParameter(0,regionFuncs[1]->Eval(regionEnergyDivisions[1]));
//   for (int i = 0;i< polynomialDegree;i++)
//   {
//        regionFuncs[0]->SetParameter(i+1,0);
//   }
//
//   regionFuncs[0]->Draw("SAME");
//
//   TString newName=string(regionFuncs[numRegions-1]->GetName())+"_Low";
//   regionFuncsLow[numRegions-1]=(TF1*) regionFuncs[numRegions-1]->Clone(newName);
//   double LOOfCentered=regionFuncs[numRegions-1]->Eval(regionEnergyDivisions[numRegions-1]);
//   double newParam=(regionEnergyDivisions[numRegions-1]-662.)/(regionEnergyDivisions[numRegions-1]/(LOOfCentered*lowFactor)-662.);
//   cout << "Alow: " << newParam << endl;
//   regionFuncsLow[numRegions-1]->SetParameter(0,newParam);
//   regionFuncsLow[numRegions-1]->Draw("SAME");
//   regionFuncsLow[numRegions-1]->SetLineColor(kBlue);
//
//   newName=string(regionFuncs[numRegions-1]->GetName())+"_High";
//   regionFuncsHigh[numRegions-1]=(TF1*) regionFuncs[numRegions-1]->Clone(newName);
//   newParam=(regionEnergyDivisions[numRegions-1]-662.)/(regionEnergyDivisions[numRegions-1]/(LOOfCentered*highFactor)-662.);
//   cout << "Ahigh: " << newParam << endl;
//   regionFuncsHigh[numRegions-1]->SetParameter(0,newParam);
//   regionFuncsHigh[numRegions-1]->Draw("SAME");
//   regionFuncsHigh[numRegions-1]->SetLineColor(kBlue);
//
//    ///connect points at region intersections visually
//
//    for (int i = 0;i< numRegions-1;i++)
//    {
//        double x1,x2,y1,y2;
//
//        //center
//        x1=x2=regionEnergyDivisions[i+1];
//        y1=regionFuncs[i]->Eval(x1);
//        y2=regionFuncs[i+1]->Eval(x1);
//        TLine* lineToDraw=new TLine(x1,y1,x2,y2);
//        lineToDraw->SetLineColor(kBlack);
//        lineToDraw->Draw("SAME");
//
//        //Low
//        x1=x2=regionEnergyDivisions[i+1];
//        y1=regionFuncsLow[i]->Eval(x1);
//        y2=regionFuncsLow[i+1]->Eval(x1);
//        lineToDraw=new TLine(x1,y1,x2,y2);
//        lineToDraw->SetLineColor(kBlue);
//        lineToDraw->Draw("SAME");
//
//        //High
//        x1=x2=regionEnergyDivisions[i+1];
//        y1=regionFuncsHigh[i]->Eval(x1);
//        y2=regionFuncsHigh[i+1]->Eval(x1);
//        lineToDraw=new TLine(x1,y1,x2,y2);
//        lineToDraw->SetLineColor(kBlue);
//        lineToDraw->Draw("SAME");
//
//    }
//
//    //Add old Matt Param model
//    const int numOldPoints=1000;
//    double xOld[numOldPoints];
//    double yOld[numOldPoints];
//    for (int i = 0;i< 1000;i++)
//    {
//        xOld[i]=i;
//        yOld[i]=evalMattKhodyukParam(xOld[i]);
//    }
//    TGraph oldGraph(numOldPoints,xOld,yOld);
//    oldGraph.SetLineColor(kPink);
//    oldGraph.Draw("LSame");
//
//    theLegend->AddEntry("region0Func","Fit","L");
//    theLegend->AddEntry("region0Func_Low","#pm 3% of fit","L");
//    theLegend->AddEntry(&oldGraph,"oldKhodyuk Model","L");
//
//
//    theLegend->Draw("same");
//
//
//
//    TString imagePath2=GraphsDir+"EnergyResponseOverallFit130710.png";
//    c2->SaveAs(imagePath2);
//    delete c2;
//
//
//    for (int i = 0;i< numRegions;i++)
//    {
//        delete regionFuncs[i];
//
//    }
//
//    outFile.close();
//
//    delete theLegend;

///make graphs of new all exp LO param model 130710
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD NonRad W Fermi");
//    fourBodySet->makeDerivedResultsFiles(GAMMA);
//
////    TH1D* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut&&"Entry$ % 4 == 0",BGO_DIM);
////    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1D* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut&&"Entry$ % 4 == 1",BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1D* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",STD_EPCut&&"Entry$ % 4 == 2",BGO_DIM);
//    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist130710_EDepGWLOLow_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1D* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",STD_EPCut&&"Entry$ % 4== 3",BGO_DIM);
//    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist130710_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130710.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp. Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison130710_Residuals.png");
//
//      TString drawString;
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130709_G7_BGOEDepAvg_EPG.txt";
//
//
////    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
////    expHist[0]->Sumw2();
////
////    for (int i = 0;i< 12;i++)
////    {
////        TString cutString="SBDToFP < 10";
////        TCut EPCut=TCut(cutString);
////
////        drawString="SPECIAL_EDepG";
////        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, EPCut, -(i+1),BGO_DIM);
////        expHist[0]->Add(tempHist);
////        delete tempHist;
////    }
////    expHist[0]->Scale(1/12.);
////    convertTH1ToTXT(expHist[0],histPath);
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130710.png");
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison #frac{EXP-MC}{MC}","ParamModelComparison130710_ResidualsVsExp.png");

///Check new LO fit function

//int numPoints=1000;
//double x[numPoints];
//double y[numPoints];
//double yLow[numPoints];
//double yHigh[numPoints];
//
//for (int i = 0;i< 1000;i++)
//{
//    x[i]=i;
//    y[i]= evalMattAllExpParamModel(x[i], ALLEXP_PARAM_MODEL);
//    yLow[i]= evalMattAllExpParamModel(x[i], ALLEXP_PARAM_MODEL_LOW);
//    yHigh[i]= evalMattAllExpParamModel(x[i], ALLEXP_PARAM_MODEL_HIGH);
//}
//TGraph loGraph(numPoints,x,y);
//TGraph loGraphLow(numPoints,x,yLow);
//TGraph loGraphHigh(numPoints,x,yHigh);
//TCanvas* c2 = new TCanvas("fitFunc","fitFunc",10,10,1280,1024);
//loGraph.GetXaxis()->SetRangeUser(0,200);
//loGraph.GetYaxis()->SetRangeUser(0.3,1);
//loGraph.Draw("AL");
//loGraphLow.Draw("L");
//loGraphHigh.Draw("L");
//
//TString imagePath=GraphsDir+"test.png";
//c2->SaveAs(imagePath);

///Testing electric fields
//    ROOTField eField;
//    eField.loadFieldFromFile(MRKRunDir+"COMSOL_SBDElectricField_Fine130515.root", "sbd", 1, 3,  3);
//
//    TVector3 outVect(0,0,0);
//    eField.linearInterp3D(TVector3(-4.6/100.,0,-27.4/100.),outVect);
//    outVect.Print();

///Convert COMSOL files to root

	///SBD
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_SBDElectricField_Fine130515.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/comsol/COMSOL_SBDElectricField_Fine130515.root";
//    ROOTField newField(newFieldToTest,"sbd",1,3,3);
//    newField.saveFieldToFile(fieldToSaveAs,"sbd");

	///Mirror
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_MirrorElectricField_Fine130515.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/comsol/COMSOL_MirrorElectricField_Fine130515.root";
//    ROOTField newField(newFieldToTest,"mirror",1,3,3);
//    newField.saveFieldToFile(fieldToSaveAs,"mirror");

	///BAPD
//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130528.txt";
//    TString fieldToSaveAs="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130528.root";
//    ROOTField newField(newFieldToTest,"bapd",1,3,3);
//    newField.saveFieldToFile(fieldToSaveAs,"bapd");

///Create Tim's EG0_Exp extraction
//    HistDim extractDim={80,0,800};
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD NonRad W Fermi");
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"ExpMatt_Hists130709_G7_BGOEDepAvg_EPG.txt",BGO_DIM);
//    expHist->Rebin(10);
//
//    TH1D* eg0ExpHist=expSpecToTheorySpec( expHist,  threeBodySet, fourBodySet,STD_EPCut, extractDim, extractDim);
//    convertTH1ToTXT(eg0ExpHist,HistsDir+"expEG0Hist130716.txt");
////    TH1D* eg0ExpHist=getTabSeperatedHist(HistsDir+"expEG0Hist130716.txt",extractDim);
//
//    eg0ExpHist->SetTitle("Exp;Energy of Decay (keV);AU");
//
//
//
////    TH1D* eg0Hist=(TH1D*) fourBodySet->makeHist("eg0RawHist","eg0","",extractDim);
////    eg0Hist->Scale(0.010212/(fourBodySet->GetNumFiles()*1000000.));
////    convertTH1ToTXT(eg0Hist,HistsDir+"test2.txt");
//////    TH1D* eg0Hist=getTabSeperatedHist(HistsDir+"test2.txt",extractDim);
//
//    HistDim oldDim={800,0,800};
//    TH1D* eg0Hist=getTabSeperatedHist(HistsDir+"MC_Hist130501_EG0_Fermi_130429.txt",oldDim);
//    eg0Hist->Rebin(10);
//    eg0Hist->Scale(0.010212/(48*1000000.));
//
//
//    eg0Hist->SetTitle("Theory;Energy of Decay (keV);AU");
//    eg0Hist->GetYaxis()->SetRangeUser(0,0.001);
//
//
//    TH1* theHists[2]={eg0Hist,eg0ExpHist};
//    plotExpVersusMCToImage(0,NULL,2,theHists,"eg0","eg0Exp130716.png");

///Check two different versions of the same MC 130712
//
//    TH1* mcHists[2]={getTabSeperatedHist(HistsDir+"MC_Hist130617_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM),getTabSeperatedHist(HistsDir+"MC_Hist130627b_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM)};
//    mcHists[0]->SetTitle("A");
//    mcHists[1]->SetTitle("B");
//    plotResidualsToImage(1,&(mcHists[1]),1,mcHists,"test","test.png");

///make graphs of new all exp LO param model but for 2-25 us this time 130712
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD NonRad W Fermi");
////    fourBodySet->makeDerivedResultsFiles(GAMMA);
//
//    TCut EPCut=STD_ECut && "pChn.SBDEDepTotal > 10. && pChn.SBDEDepTotal < 31. && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 25.e-6";
//
//    TH1* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",EPCut&&"Entry$ % 4 == 0",BGO_DIM);
//    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",EPCut&&"Entry$ % 4 == 1",BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130712_EDepGWLO_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130712_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",EPCut&&"Entry$ % 4 == 2",BGO_DIM);
//    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist130712_EDepGWLOLow_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130712_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",EPCut&&"Entry$ % 4== 3",BGO_DIM);
//    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist130712_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt");
////    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130712_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130712.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp. Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison130712_Residuals.png");
//
//      TString drawString;
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130712_G7_BGOEDepAvg_EPG.txt";
//
//
//    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
//    expHist[0]->Sumw2();
//
//    for (int i = 0;i< 12;i++)
//    {
////        TString cutString="SBDToFP < 10";
////        TCut EPCut=TCut(cutString);
//
//        drawString="SPECIAL_EDepG";
//        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "", -(i+1),BGO_DIM);
//        expHist[0]->Add(tempHist);
//        delete tempHist;
//    }
//    expHist[0]->Scale(1/12.);
//    convertTH1ToTXT(expHist[0],histPath);
//
////    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130712.png");
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison #frac{EXP-MC}{MC}","ParamModelComparison130712_ResidualsVsExp.png");

///Create cycle based 10-25us histograms to understand Jeff's post
//    TH1* expHistSingle[6];
//    TH1* expHistSum[6];
//    for (int i = 0;i< 6;i++)
//    {
//
//
//
//
//        TString histName = "Cycle"+int2str(i+3)+"Single";
//        expHistSingle[i]=makeTH1DFromDim(histName,histName,BGO_DIM);
//        expHist[i]->Sumw2();
//
//        for (int i = 0;i< 12;i++)
//        {
//    //        TString cutString="SBDToFP < 10";
//    //        TCut EPCut=TCut(cutString);
//
//            drawString="SPECIAL_EDepG";
//            TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "", -(i+1),BGO_DIM);
//            expHist[i]->Add(tempHist);
//            delete tempHist;
//        }
//        expHist[i]->Scale(1/12.);
//        TString histPath=HistsDir+"ExpMatt_Hists130712_G7_Cycle"+int2str(i+3)+"_BGOEDepAvg_EPG.txt";
//        convertTH1ToTXT(expHist[i],histPath);
//    }

///make graphs of new all exp LO param model 130712 - 2-10us but with no Fermi Func
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD NonRad W Fermi");
//    fourBodySet->makeDerivedResultsFiles(GAMMA);
//
//    TH1* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut&&"Entry$ % 4 == 0",BGO_DIM);
//    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//
//    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut&&"Entry$ % 4 == 1",BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130712b_EDepGWLO_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130712b_EDepGWLO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    TH1* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",STD_EPCut&&"Entry$ % 4 == 2",BGO_DIM);
//    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist130712b_EDepGWLOLow_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130712b_EDepGWLOLow_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    TH1* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",STD_EPCut&&"Entry$ % 4== 3",BGO_DIM);
//    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist130712b_EDepGWLOHigh_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130712b_EDepGWLOHigh_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130712b.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp. Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison130712b_Residuals.png");
//
//      TString drawString;
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130709_G7_BGOEDepAvg_EPG.txt"; //Don't need to update this
//
//
////    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
////    expHist[0]->Sumw2();
////
////    for (int i = 0;i< 12;i++)
////    {
////        TString cutString="SBDToFP < 10";
////        TCut EPCut=TCut(cutString);
////
////        drawString="SPECIAL_EDepG";
////        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, EPCut, -(i+1),BGO_DIM);
////        expHist[0]->Add(tempHist);
////        delete tempHist;
////    }
////    expHist[0]->Scale(1/12.);
////    convertTH1ToTXT(expHist[0],histPath);
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130712b.png");
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison #frac{EXP-MC}{MC}","ParamModelComparison130712b_ResidualsVsExp.png");

///make graphs of new all exp LO param model 130712 - 2-25us and old proton energy cut both with no Fermi Func
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD NonRad W Fermi");
////    fourBodySet->makeDerivedResultsFiles(GAMMA);
//
//    TH1* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_OLDEPCut&&"Entry$ % 4 == 0",BGO_DIM);
//    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//
//    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_OLDEPCut&&"Entry$ % 4 == 1",BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist130712c_EDepGWLO_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130712c_EDepGWLO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    TH1* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",STD_OLDEPCut&&"Entry$ % 4 == 2",BGO_DIM);
//    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist130712c_EDepGWLOLow_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130712c_EDepGWLOLow_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    TH1* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",STD_OLDEPCut&&"Entry$ % 4== 3",BGO_DIM);
//    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist130712c_EDepGWLOHigh_AverageIndiv_EPG_RID130.txt");
////    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130712c_EDepGWLOHigh_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130712c.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp. Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison130712c_Residuals.png");
//
//      TString drawString;
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130712c_G7_BGOEDepAvg_EPG.txt"; //Don't need to update this
//
//
//    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
//    expHist[0]->Sumw2();
//
//    for (int i = 0;i< 12;i++)
//    {
////        TString cutString="SBDToFP < 10";
////        TCut EPCut=TCut(cutString);
//
//        drawString="SPECIAL_EDepG";
//        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "", -(i+1),BGO_DIM);
//        expHist[0]->Add(tempHist);
//        delete tempHist;
//    }
//    expHist[0]->Scale(1/12.);
//    convertTH1ToTXT(expHist[0],histPath);
//
////    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130712c.png");
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison #frac{EXP-MC}{MC}","ParamModelComparison130712c_ResidualsVsExp.png");

///Start comparing graphs to figure out why things are different 130712
//    const int numMCHists=6;
//    TH1* mcHists[numMCHists];
//    //2-25us
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//    mcHists[0]->SetTitle("RID130old");
//    mcHists[1]=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//    mcHists[1]->SetTitle("RID130new");
//    mcHists[2]=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//    mcHists[2]->SetTitle("RID175(wFermi)");
//
//    //2-10us
//    mcHists[3]=getTabSeperatedHist(HistsDir+"MC_Hist130617_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//    mcHists[3]->SetTitle("RID175Khodyuk");
//    mcHists[4]=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//    mcHists[4]->SetTitle("RID175Param");
//    //10-25us
//    mcHists[5]=getTabSeperatedHist(HistsDir+"MC_Hist130617_EDepG_AverageIndiv_10us_25us_EPG_RID175.txt",BGO_DIM);
//    mcHists[5]->SetTitle("RID175Khodyuk_10_25us");
//
//
//     const int numExpHists=4;
//     TH1* expHists[numExpHists];
//     expHists[0]=getTabSeperatedHist(HistsDir+"EXPJeff_Hist130422_EDepG_BGO_AverageIndiv_EPG.txt",BGO_DIM);
//     expHists[0]->Scale(8.683e-5/(double) expHists[0]->Integral()); //Jeff's hist is the sum of the 12 individual detectors and he did not calculate EPG/EP
//     expHists[0]->SetTitle("JeffOld");
//
//     expHists[1]=getTabSeperatedHist(HistsDir+"ExpJeff_Hists130712_2_10_BGOEDepAvg_EPG.txt",BGO_DIM);
//     expHists[1]->SetTitle("JeffNew2_10");
//
//     expHists[2]=getTabSeperatedHist(HistsDir+"ExpJeff_Hists130712_10_25_BGOEDepAvg_EPG.txt",BGO_DIM);
//     expHists[2]->SetTitle("JeffNew10_25");
//
//     expHists[3]=getTabSeperatedHist(HistsDir+"ExpJeff_Hists130712_2_25_BGOEDepAvg_EPG.txt",BGO_DIM);
//     expHists[3]->SetTitle("JeffNew2_25");
//
//
//     int rebin=50;
//     for (int i = 0;i< numMCHists;i++)
//     {
//     	mcHists[i]->Rebin(rebin);
//     	mcHists[i]->GetXaxis()->SetRangeUser(0,800);
//     }
//
//     for (int i = 0;i< numExpHists;i++)
//     {
//     	expHists[i]->Rebin(rebin);
//     	expHists[i]->GetXaxis()->SetRangeUser(0,800);
//     }

//    plotResidualsAndComparisonToImage(1,expHists,3,mcHists,"Old Exp/Old Cuts #frac{Exp-MC}{MC}","OldExp_OldCuts130715.png");
//     plotResidualsAndComparisonToImage(1,expHists,1,&(mcHists[2]),"Old Exp/Old Cuts/Fermi  #frac{Exp-MC}{MC}","OldExp_OldCuts130715_ZOOM.png");
//     plotResidualsAndComparisonToImage(1,&(expHists[1]),2,&(mcHists[3]),"New Exp/New Cuts/Fermi #frac{Exp-MC}{MC}","NewExp_NewCuts130715_Rebin.png");
//     plotResidualsAndComparisonToImage(1,&(expHists[1]),1,&(mcHists[3]),"New Exp/New Cuts/Fermi  #frac{Exp-MC}{MC}","NewExp_NewCuts130715_ZOOM.png");
//
//    plotResidualsAndComparisonToImage(1,&(expHists[3]),numMCHists,mcHists,"New Jeff Hist #frac{Exp-MC}{MC}","test2.png");
//    plotResidualsAndComparisonToImage(1,&(expHists[3]),1,expHists,"New Jeff Hist #frac{Exp-MC}{MC}","test3.png");

//    plotResidualsAndComparisonToImage(2,&(expHists[0]),1,&(mcHists[3]),"#frac{Exp-MC}{MC}"","test.png");

///Provide Kevin Graph of Parameterized LO Model

//    ofstream outFile;
//
//    outFile.open((HistsDir+"ParamLOModel130830.txt"));
//
//    outFile << "#energy(keV)\tlower yield curve\tbest estimate yield curve\tupper yield curve" << endl;
//    double energy, lower,best,upper;
//    for (int i = 0;i< BGOCAL_DIM.numBins;i++)
//    {
//        energy=i/10. + 0.05;
//        lower=evalMattAllExpParamModel(energy,ALLEXP_PARAM_MODEL_LOW);
//        best=evalMattAllExpParamModel(energy,ALLEXP_PARAM_MODEL);
//        upper=evalMattAllExpParamModel(energy,ALLEXP_PARAM_MODEL_HIGH);
//
//        outFile << energy << "\t" << lower << "\t" << best << "\t" << upper;
//        if(i != BGOCAL_DIM.numBins -1)
//            outFile <<endl;
//
//    }
//    outFile.close();

///Compare geant4 physics libraries 130716

//    TString physicsList="emlivermore empenelope emstandard_opt0 emstandard_opt1 emstandard_opt2 emstandard_opt3 emstandard_opt4 standardSS standardWVI geant4_9_5p1_Livermore";
//
//    int RIDNumber=216;
//    const int numHists=numItemsInStringList(physicsList);
//    TH1* mcHists[numHists];
//    for (int i = 0;i< numHists;i++)
//    {
//        TString RIDString;
//        if(i<8)
//        {
//            RIDString="RID"+int2str(RIDNumber+i);
//        }
//        else
//        {
//            RIDString="RID"+int2str(RIDNumber+i+1);
//        }
//        TString histPath=HistsDir+"MC_Hist130627_DirectEDepG_AverageIndiv_EPG_"+RIDString+".txt";
//        if(!FileExists(histPath))
//        {
//
//            RDK2Set fourBodySet(1,"130429","STD_4B_100eV_Endpoint","RID213","RID214",RIDString,"special");  //RID213 and 214 are dummies since we dont' care about coincidences
//            mcHists[i]=makeMCIndividualDetectorHist(NULL,&fourBodySet,"BGOEDepBlurEdge","",BGO_DIM);
//            convertTH1ToTXT(mcHists[i],histPath);
//        }
//        else
//        {
//            mcHists[i]=getTabSeperatedHist(histPath,BGO_DIM);
//        }
//        TString thePhysics;
//        stringFromList(physicsList,i,thePhysics);
//        TString titleString=thePhysics+";Energy Detected (keV);Counts Detected (AU)";
//        mcHists[i]->SetTitle(titleString);
//        mcHists[i]->Rebin(20);
//        mcHists[i]->GetYaxis()->SetRangeUser(0,400);
//    }
//
//
//    plotExpVersusMCToImage(numHists-1,mcHists,1,&(mcHists[numHists-1]),"Geant4.9.6.p02 Physics Libraries compared to Geant4.9.5.p01","PhysicsLibraryTest_BGO_Geant4.9.6.p02.png");
//    plotResidualsToImage(numHists-1,mcHists,1,&(mcHists[numHists-1]),"Geant4.9.6.p02 Physics Libraries #frac{X - geant4_9_5p1_Livermore}{geant4_9_5p1_Livermore}","PhysicsLibraryTest_BGO_Geant4.9.6.p02.png",-.2,.2);

///Check Geant4.9.6.p02 overall 130717
//  RDK2Set* threeBodySet=new RDK2Set(1,"130429","STD_3B","RID211","RID212","","STD NonRad W Fermi");
//  RDK2Set* fourBodySet=new RDK2Set(1,"130429","STD_4B_100eV_Endpoint","RID213","RID214","RID215","STD NonRad W Fermi");
//  threeBodySet->makeDerivedResultsFiles(PROTON);
//  threeBodySet->makeDerivedResultsFiles(ELECTRON);
//  fourBodySet->makeDerivedResultsFiles(PROTON);
//  fourBodySet->makeDerivedResultsFiles(ELECTRON);
//  fourBodySet->makeDerivedResultsFiles(GAMMA);
//
//
//    TH1* mcHistNew=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,BGO_DIM);
//
//
//
//    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130709_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//
//
//    TH1* mcHists[2] = {mcHistOld,mcHistNew};
//
//    mcHists[0]->SetTitle("Old Geant;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Geant;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//
//    for (int i = 0;i< 2;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"Geant Version Comparison","GeantVersionComparison130717.png");
//    plotResidualsToImage(1,&(mcHists[1]),1,mcHists,"Geant Version Comparison #frac{New-OLD}{OLD}","GeantVersionComparison130717_Residuals.png");

///Create macros to test physics libraries in Geant4.9.6.p02 130717 for electrons

//    TString standardLines="/MRK/det/setBAPDFieldScale 0";
//    TString physicsList="empenelope emstandard_opt0 emstandard_opt1 emstandard_opt2 emstandard_opt3 emstandard_opt4 standardSS standardSSM standardWVI standardGS emDNA";
//
//    int RIDNumber=227;
//    for (int i = 0;i< numItemsInStringList(physicsList);i++)
//    {
//        TString thePhysics;
//        stringFromList(physicsList,i,thePhysics);
//    	makeGeantMacroFilesForClusterEventGenerator(0,1,"RID"+int2str(RIDNumber),ELECTRON,THREEBODY,standardLines+"\n/MRK/phys/addPhysics "+ thePhysics);
//    	RIDNumber++;
//    }

///Create macros for Geant4.9.6.p02 48 million full run of Livermore 130717
//    standardLines="/MRK/det/setBAPDFieldScale 0";
//
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID239",PROTON,THREEBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID240",ELECTRON,THREEBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID241",PROTON,FOURBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID242",ELECTRON,FOURBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID243",GAMMA,FOURBODY,standardLines);

///Create Tim's EG0_Exp extraction
//    HistDim extractDim={80,0,800};
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD Rad W Fermi");
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"ExpMatt_Hists130709_G7_BGOEDepAvg_EPG.txt",BGO_DIM);
//    expHist->Rebin(10);
//
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130709_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//    mcHist->Rebin(10);
//
////    TH1D* eg0ExpHist=expSpecToTheorySpecWError( expHist, threeBodySet,fourBodySet,STD_EPCut,HistsDir+"EG0_EDet2DGraph_130717.root");
////    convertTH1ToTXT(eg0ExpHist,HistsDir+"expEG0Hist130717.txt");
//    TH1D* eg0ExpHist=getTabSeperatedHist(HistsDir+"expEG0Hist130716.txt",extractDim);
//
//    eg0ExpHist->SetTitle("Exp;Energy of Decay (keV);Probability per decay per bin");
//
////    TH1* percentHist[1]={createEG0HitPercentHist(threeBodySet,fourBodySet,STD_EPCut, extractDim)};
////    convertTH1ToTXT(percentHist[0],HistsDir+"eg0HitPercentage130722.txt");
////    plotExpVersusMCToImage(0,NULL,1,percentHist,"EG0HistPercentage","eg0HitPercentage130722.png");
//
//
//
//
//
//
//
////    TH1D* eg0Hist=(TH1D*) fourBodySet->makeHist("eg0RawHist","eg0","",extractDim);
////    eg0Hist->Scale(0.010212/(fourBodySet->GetNumFiles()*1000000.));
////    convertTH1ToTXT(eg0Hist,HistsDir+"test2.txt");
//////    TH1D* eg0Hist=getTabSeperatedHist(HistsDir+"test2.txt",extractDim);
//
//    HistDim oldDim={800,0,800};
//    TH1D* eg0Hist=getTabSeperatedHist(HistsDir+"MC_Hist130501_EG0_Fermi_130429.txt",oldDim);
//    eg0Hist->Rebin(10);
//    eg0Hist->Scale(0.010212/(48*1000000.));
//
//    eg0ExpHist->SetBinContent(1,0);
//    eg0Hist->SetBinContent(1,0);
//    eg0Hist->SetBinError(1,0);
//    eg0Hist->SetBinError(1,0);
//
//
//    eg0Hist->SetTitle("Theory;Energy of Decay (keV);Probability per decay per bin");
//    eg0Hist->GetYaxis()->SetRangeUser(0,0.001);
//
//
//    eg0Hist->GetYaxis()->SetTitleOffset(1.8);
//    eg0ExpHist->GetYaxis()->SetTitleOffset(1.8);
////    gStyle->SetTitleYOffset(1.2);
//
//    TH1* mcHists[1]={eg0Hist};
//    TH1* expHists[1]={eg0ExpHist};
//    plotResidualsToImage(1,expHists,1,mcHists,"Residuals (Exp - Theory)","eg0Exp130717_UnNormalizedResidual.png",-0.00001,.000025,false);
//// plotExpVersusMCToImage(0,NULL,1,mcHists,"Initial Gamma Energy of Neutron Decay","eg0130722.png");
//
//    eg0Hist->GetXaxis()->SetRangeUser(0,300);
//    eg0ExpHist->GetXaxis()->SetRangeUser(0,300);
//
//
//
////    plotResidualsAndComparisonToImage(1,expHists,1,mcHists,"Initial Gamma Energy of Neutron Decay","eg0Exp130717_Zoom.png");
//
//    expHists[0]=expHist;
//    mcHists[0]=mcHist;
//
//    expHist->SetTitle("Exp;Energy Detected (keV);Counts");
//    mcHist->SetTitle("MC+Theory;Energy Detected (keV);#frac{EPG}{EP}");
//
//    expHist->Scale(4.485e8);
//    expHist->GetYaxis()->SetRangeUser(0,9500);
//    plotExpVersusMCToImage(1,expHists,0,mcHists,"#Experimental Spectrum after Background Subtraction","exp130722.png");
//
//
////    expHist->GetXaxis()->SetRangeUser(0,800);
////    mcHist->GetXaxis()->SetRangeUser(0,800);
////    plotResidualsAndComparisonToImage(1,expHists,1,mcHists,"#frac{EPG}{EP} Comparison","expVsMC130722.png",-.6,1.5);

///Compare geant4 physics libraries: electrons 130722

//    TString physicsList="emlivermore empenelope em_opt4 standardSS OLD_Livermore";
//    TString RIDList="212 227 231 232 172";
//
//    const int numHists=numItemsInStringList(physicsList);
//    TH1* mcHists[numHists];
//    for (int i = 0;i< numHists;i++)
//    {
//        TString RIDString;
//        stringFromList(RIDList,i,RIDString);
//        RIDString="RID"+RIDString;
//        TString histPath=HistsDir+"MC_Hist130722_EDepE_AverageIndiv_"+RIDString+".txt";
//        if(!FileExists(histPath))
//        {
//
//            RDK2Set threeBodySet(1,"130429","STD_3B","RID171",RIDString,"","special");  //RID171 is a dummy since we dont' care about coincidences
//            TString histName=RIDString+"Hist";
//            mcHists[i]=threeBodySet.makeHist(histName,"eChn.SBDEDepTotal","",EDEPE_DIM);
//            convertTH1ToTXT(mcHists[i],histPath);
//        }
//        else
//        {
//            mcHists[i]=getTabSeperatedHist(histPath,BGO_DIM);
//        }
//        TString thePhysics;
//        stringFromList(physicsList,i,thePhysics);
//        TString titleString=thePhysics+";Electron Energy Absorbed (keV);Counts Detected (AU)";
//        mcHists[i]->SetTitle(titleString);
//        mcHists[i]->SetBinContent(1,0);
//        mcHists[i]->Rebin(50);
//        mcHists[i]->GetYaxis()->SetRangeUser(0,800*50);
//
//    }
//
//
//    plotExpVersusMCToImage(numHists-1,mcHists,1,&(mcHists[numHists-1]),"Geant4.9.6.p02 Physics Libraries compared to Geant4.9.5.p01","PhysicsLibraryTest_SBDE_Geant4.9.6.p02.png");
//    plotResidualsToImage(numHists-1,mcHists,1,&(mcHists[numHists-1]),"Geant4.9.6.p02 Physics Libraries #frac{X - geant4_9_5p1_Livermore}{geant4_9_5p1_Livermore}","PhysicsLibraryTest_SBDE_Geant4.9.6.p02_Residual.png",-.15,.15);

	///Compare geant4 physics libraries: Gamma 130722

//    TString physicsList="emlivermore empenelope em_opt4 OLD_Livermore";
//    TString RIDList="216 217 222 226";
//
//    const int numHists=numItemsInStringList(physicsList);
//    TH1* mcHists[numHists];
//    for (int i = 0;i< numHists;i++)
//    {
//        TString RIDString;
//        stringFromList(RIDList,i,RIDString);
//        RIDString="RID"+RIDString;
//        TString histPath=HistsDir+"MC_Hist130722_DirectEDepG_AverageIndiv_"+RIDString+".txt";
//        if(!FileExists(histPath))
//        {
//            TString histName=RIDString+"Hist";
//            RDK2Set fourBodySet(1,"130429","STD_4B_100eV_Endpoint","RID213","RID214",RIDString,"special");  //RID213 and 214 are dummies since we dont' care about coincidences
//            mcHists[i]=fourBodySet.makeHist(histName,"gChn.BGO1EDepTotal","",BGO_DIM);
//            convertTH1ToTXT(mcHists[i],histPath);
//        }
//        else
//        {
//            mcHists[i]=getTabSeperatedHist(histPath,BGO_DIM);
//        }
//        TString thePhysics;
//        stringFromList(physicsList,i,thePhysics);
//        TString titleString=thePhysics+";Gamma Energy Absorbed (keV);Counts Detected (AU)";
//        mcHists[i]->SetTitle(titleString);
//        mcHists[i]->SetBinContent(1,0);
//        mcHists[i]->Rebin(20);
////        mcHists[i]->GetYaxis()->SetRangeUser(0,0.00025);
//
//    }
//
//
//    plotExpVersusMCToImage(numHists-1,mcHists,1,&(mcHists[numHists-1]),"Geant4.9.6.p02 Physics Libraries compared to Geant4.9.5.p01","PhysicsLibraryTest_BGOE_Geant4.9.6.p02.png");
//    plotResidualsToImage(numHists-1,mcHists,1,&(mcHists[numHists-1]),"Geant4.9.6.p02 Physics Libraries #frac{X - geant4_9_5p1_Livermore}{geant4_9_5p1_Livermore}","PhysicsLibraryTest_BGOE_Geant4.9.6.p02_Residual.png",-.15,.15);

///Testing and checking new field files 130723
//    TString filePath;
//
//    filePath=MRKRunDir+"NewSBDField.root";
//    TFile newFile(filePath,"READ");
//    TTree* newTree=(TTree*) newFile.Get("g");
//    gROOT->cd();
//    TH1D* newHist=makeTH1DFromDim("newHist","New Field;Proton Energy Absorbed (keV);Counts (AU)",EDEPP_DIM);
//    newTree->Draw("SBDEDepTotal >> newHist","SBDEDepTotal >0");
//    newFile.Close();
//    newTree=NULL;
//
//    filePath=MRKRunDir+"OldSBDField.root";
//    TFile oldFile(filePath,"READ");
//    TTree* oldTree=(TTree*) oldFile.Get("g");
//    gROOT->cd();
//    TH1D* oldHist=makeTH1DFromDim("oldHist","Old Field;Proton Energy Absorbed (keV);Counts (AU)",EDEPP_DIM);
//    oldTree->Draw("SBDEDepTotal >> oldHist","SBDEDepTotal >0");
//    oldFile.Close();
//    oldTree=NULL;
//
//    TH1* mcHists[2]={newHist,oldHist};
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"Old/New SBD Field","test.png");
//
//    delete newHist;
//    delete oldHist;

	///SBD
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_SBDElectricField_Fine130515.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_SBDElectricField_Fine130515.root";
//    ROOTField sbdField(newFieldToTest,"sbd",1,3,3);
//    sbdField.saveFieldToFile(fieldToSaveAs,"sbd");

	///Mirror
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_MirrorElectricField_Fine130515.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_MirrorElectricField_Fine130515.root";
//    ROOTField mirrorField(newFieldToTest,"mirror",1,3,3);
//    mirrorField.saveFieldToFile(fieldToSaveAs,"mirror");

	///BAPD
//    TString newFieldToTest="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130528.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_BAPDElectricField_Fine130528.root";
//    ROOTField bapdField(newFieldToTest,"bapd",1,3,3);
//    bapdField.saveFieldToFile(fieldToSaveAs,"bapd");

//    double distanceDownBend;
//    double offsetX=0;
//    double angle=9.5;
//
//
//    ROOTField newField;
//    newField.loadFieldFromFile(MRKRunDir+"COMSOL_SBDElectricField_Fine130515.root", "sbd", 1, 3,  3);
//
//    ROOTField oldField;
//    oldField.loadFieldFromFile(MRKRunDir+"Efield_RDK2_SBD_Matt5.root", "sbd", 1, 2,  2);
//    oldField.setRotation(0,-9.5*PI/180.,0);
//
//    for (int i = 0;i< 100;i++)
//    {
//        cout << "Distance (m) down bend: ";
//        cin >> distanceDownBend;
//
//        double xPos=offsetX*cos(angle*PI/180.)+distanceDownBend*sin(angle*PI/180.);
//        double yPos=0;
//        double zPos=-offsetX*sin(angle*PI/180.)+distanceDownBend*cos(angle*PI/180.);
//
//
//        TVector3 outVect(0,0,0);
//        newField.linearInterp3D(TVector3(xPos,yPos,zPos),outVect);
//        outVect.Print();
//
//        outVect.SetXYZ(0,0,0);
//        oldField.linearInterp3D(TVector3(xPos,yPos,zPos),outVect);
//        outVect.Print();
//    }

//    TString filePath;
//
//    filePath=MRKRunDir+"NewMirrorField.root";
//    TFile newFile(filePath,"READ");
//    TTree* newTree=(TTree*) newFile.Get("g");
//    gROOT->cd();
//    TH1D* newHist=makeTH1DFromDim("newHist","New Field;Proton Energy Absorbed (keV);Counts (AU)",EDEPP_DIM);
//    newTree->Draw("SBDEDepTotal >> newHist","SBDEDepTotal >0");
//    newFile.Close();
//    newTree=NULL;
//
//    filePath=MRKRunDir+"OldMirrorField.root";
//    TFile oldFile(filePath,"READ");
//    TTree* oldTree=(TTree*) oldFile.Get("g");
//    gROOT->cd();
//    TH1D* oldHist=makeTH1DFromDim("oldHist","Old Field;Proton Energy Absorbed (keV);Counts (AU)",EDEPP_DIM);
//    oldTree->Draw("SBDEDepTotal >> oldHist","SBDEDepTotal >0");
//    oldFile.Close();
//    oldTree=NULL;
//
//    TH1* mcHists[2]={newHist,oldHist};
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"Old/New Mirror Field","test2.png");
//
//    delete newHist;
//    delete oldHist;

///Create macros for Geant4.9.6.p02 48 million full run of Livermore with new fields 130723
//    TString standardLines="/MRK/det/setMirrorFieldFilePath COMSOL_MirrorElectricField_Fine130515.root";
//    standardLines+="\n/MRK/det/setSBDFieldFilePath COMSOL_SBDElectricField_Fine130515.root";
//
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID244",PROTON,THREEBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID245",ELECTRON,THREEBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID246",PROTON,FOURBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID247",ELECTRON,FOURBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID248",GAMMA,FOURBODY,standardLines);

///Make exp set that includes new EP cuts and 1000 keV gamma cut - PID7
//
////    makeExpEPCoincidenceFiles("PID7",NUM_EXP_SERIES_G7,EXP_SERIES_ID_LIST_G7);
//
//    TH1* expHists[2];
//    TString histPath=HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt";
//
//
////    expHists[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
////    expHists[0]->Sumw2();
////
////    for (int i = 0;i< 12;i++)
////    {
////        TString drawString="BGOEDep"+int2str(i+1);
////        TH1* tempHist=makeExpEPGHist(false,"PID7", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "", -(i+1),BGO_DIM);
////        expHists[0]->Add(tempHist);
////        delete tempHist;
////    }
////    expHists[0]->Scale(1/12.);
////    convertTH1ToTXT(expHists[0],histPath);
//    expHists[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHists[1]=getTabSeperatedHist(HistsDir+"ExpJeff_Hists130712_2_10_BGOEDepAvg_EPG.txt",BGO_DIM);
//
////    histPath=HistsDir+"ExpMatt_Hists130709_G7_BGOEDepAvg_EPG.txt";
////    expHists[1]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHists[0]->Rebin(10);
//    expHists[1]->Rebin(10);
//    expHists[0]->SetTitle("Matt;Energy Detected (keV);#frac{EPG}{EP}");
//    expHists[1]->SetTitle("Jeff;Energy Detected (keV);#frac{EPG}{EP}");
//
//    plotExpVersusMCToImage(2,expHists,0,NULL,"Exp. Comparison","ExpWUpdatedCuts.png");
//    plotResidualsToImage(1,expHists,1,&(expHists[1]),"Exp. Comparison #frac{Matt-Jeff}{Jeff}","ExpWUpdatedCuts_Residuals.png");

	///Convert Bryce's full electric field to root format
//    TString newFieldToTest="/media/mjbexternal/comsol/Ef_Total.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_SBD_Mirror_BAPDElectricField_SuperFine130724.root";
//    ROOTField totalField(newFieldToTest,"total",1,3,3);
//    totalField.saveFieldToFile(fieldToSaveAs,"total");

///Create Tim's EG0_Exp extraction with corrected Exp Hist
//    HistDim extractDim={80,0,800};
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD Rad W Fermi");
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt",BGO_DIM);
//    expHist->Rebin(10);
//
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130709_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//    mcHist->Rebin(10);
//
////    TString histPath2D=HistsDir+"EG0_EDet2DGraph_130724.root";
////    TFile rootFile(histPath2D,"READ");
////    TH2D* eDetEG0Hist=(TH2D*) rootFile.Get("eDetEG0Hist");
////    eDetEG0Hist->SetDirectory(0);
////    rootFile.Close();
////
////    TCanvas* c3 = new TCanvas("c3","c3",10,10,1280,1024);
////    eDetEG0Hist->SetTitle("Energy Detected versus Energy of Decay;Energy of Decay (keV);Energy Detected (keV)");
////    eDetEG0Hist->Draw("SCAT");
////
////    TString imagePath=string(GRAPHS_DIR)+"EG0_EDet2DGraph_Scat.png";
////    c3->SaveAs(imagePath);
////    delete c3;
//
////    TH1D* eg0ExpHist=expSpecToTheorySpecWError( expHist, threeBodySet,fourBodySet,STD_EPCut,HistsDir+"EG0_EDet2DGraph_130724.root");
////    convertTH1ToTXT(eg0ExpHist,HistsDir+"expEG0Hist130724.txt");
//    TH1D* eg0ExpHist=getTabSeperatedHist(HistsDir+"expEG0Hist130716.txt",extractDim);
//
//    eg0ExpHist->SetTitle("Exp+MC;Energy of Decay (keV);Probability per decay per bin");
//
////    TH1* percentHist[1]={createEG0HitPercentHist(threeBodySet,fourBodySet,STD_EPCut, extractDim)};
////    convertTH1ToTXT(percentHist[0],HistsDir+"eg0HitPercentage130722.txt");
////    plotExpVersusMCToImage(0,NULL,1,percentHist,"EG0HistPercentage","eg0HitPercentage130722.png");
//
//
//
//
//
//
//
////    TH1D* eg0Hist=(TH1D*) fourBodySet->makeHist("eg0RawHist","eg0","",extractDim);
////    eg0Hist->Scale(0.010212/(fourBodySet->GetNumFiles()*1000000.));
////    convertTH1ToTXT(eg0Hist,HistsDir+"test2.txt");
//////    TH1D* eg0Hist=getTabSeperatedHist(HistsDir+"test2.txt",extractDim);
//
//    HistDim oldDim={800,0,800};
//    TH1D* eg0Hist=getTabSeperatedHist(HistsDir+"MC_Hist130501_EG0_Fermi_130429.txt",oldDim);
//    eg0Hist->Rebin(10);
//    eg0Hist->Scale(0.010212/(48*1000000.));
//
//    eg0ExpHist->SetBinContent(1,0);
//    eg0Hist->SetBinContent(1,0);
//    eg0Hist->SetBinError(1,0);
//    eg0Hist->SetBinError(1,0);
//
//
//    eg0Hist->SetTitle("Theory;Energy of Decay (keV);Probability per decay per bin");
//    eg0Hist->GetYaxis()->SetRangeUser(0,0.001);
//
//
//    eg0Hist->GetYaxis()->SetTitleOffset(1.8);
//    eg0ExpHist->GetYaxis()->SetTitleOffset(1.8);
////    gStyle->SetTitleYOffset(1.2);
//
//    TH1* mcHists[1]={eg0Hist};
//    TH1* expHists[1]={eg0ExpHist};
////    plotResidualsToImage(1,expHists,1,mcHists,"Residuals (Exp - Theory)","eg0Exp130724_UnNormalizedResidual.png",-0.00001,.000025,false);
//// plotExpVersusMCToImage(0,NULL,1,mcHists,"Initial Gamma Energy of Neutron Decay","eg0130724.png");
//
//    eg0Hist->GetXaxis()->SetRangeUser(0,800);
//    eg0ExpHist->GetXaxis()->SetRangeUser(0,800);
//
//
//
//    plotResidualsAndComparisonToImage(1,expHists,1,mcHists,"Initial Gamma Energy of Neutron Decay","eg0Exp130724.png",-1.,5.);
//
////    expHists[0]=expHist;
////    mcHists[0]=mcHist;
////
////    expHist->SetTitle("Exp;Energy Detected (keV);Counts");
////    mcHist->SetTitle("MC+Theory;Energy Detected (keV);#frac{EPG}{EP}");
//
////    expHist->Scale(4.485e8);
////    expHist->GetYaxis()->SetRangeUser(0,9500);
////    plotExpVersusMCToImage(1,expHists,0,mcHists,"#Experimental Spectrum after Background Subtraction","exp130724.png");
//
//
////    expHist->GetXaxis()->SetRangeUser(0,300);
////    mcHist->GetXaxis()->SetRangeUser(0,300);
////    mcHist->GetYaxis()->SetRangeUser(0,23e-6);
////    plotResidualsAndComparisonToImage(1,expHists,1,mcHists,"#frac{EPG}{EP} Comparison","expVsMC130724_Zoom.png",-.6,1.5);

///Replot post on parameterized LO
//
//    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist130627_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    	mcHists[i]->GetXaxis()->SetRangeUser(0,300);
//    }
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt";
//
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//    expHist[0]->GetXaxis()->SetRangeUser(0,300);
//
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison130724Zoom.png");
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison EXP-MC","ParamModelComparison130724_UnnormalizedResidualsVsExpZoom.png",-2e-7,2e-6,false);

///Make exp set that includes new EP cuts and 1000 keV gamma cut - PID7 - but now using larger pre peak and new post peak window

////    makeExpEPCoincidenceFiles("PID7",NUM_EXP_SERIES_G7,EXP_SERIES_ID_LIST_G7);
//
//    TH1* expHists[2];
//    TString histPath=HistsDir+"ExpMatt_Hists130725_G7_BGOEDepAvg_EPG.txt";
//
//
//    expHists[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
//    expHists[0]->Sumw2();
//
//    for (int i = 0;i< 12;i++)
//    {
//        TString drawString="BGOEDep"+int2str(i+1);
//        TH1* tempHist=makeExpEPGHist(false,"PID7", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "", -(i+1),BGO_DIM);
//        expHists[0]->Add(tempHist);
//        delete tempHist;
//    }
//    expHists[0]->Scale(1/12.);
//    convertTH1ToTXT(expHists[0],histPath);
////    expHists[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHists[1]=getTabSeperatedHist(HistsDir+"ExpJeff_Hists130712_2_10_BGOEDepAvg_EPG.txt",BGO_DIM);
//
////    histPath=HistsDir+"ExpMatt_Hists130709_G7_BGOEDepAvg_EPG.txt";
////    expHists[1]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHists[0]->Rebin(10);
//    expHists[1]->Rebin(10);
//    expHists[0]->SetTitle("Matt;Energy Detected (keV);#frac{EPG}{EP}");
//    expHists[1]->SetTitle("Jeff;Energy Detected (keV);#frac{EPG}{EP}");
//
//    expHists[0]->GetXaxis()->SetRangeUser(0,300);
//    expHists[1]->GetXaxis()->SetRangeUser(0,300);
//
//    plotExpVersusMCToImage(2,expHists,0,NULL,"Exp. Comparison","ExpWUpdatedCuts130725.png");
//    plotResidualsToImage(1,expHists,1,&(expHists[1]),"Exp. Comparison #frac{Matt-Jeff}{Jeff}","ExpWUpdatedCuts_Residuals130725.png");

	///Make exp set that includes new EP cuts and 1000 keV gamma cut - PID7 - but now using larger pre peak and new post peak window and SBDEPkChn

//    makeExpEPCoincidenceFiles("PID8",NUM_EXP_SERIES_G7,EXP_SERIES_ID_LIST_G7);  //Remake with SBDEPkCh
//
	///BGO Energies
//    HistDim PkChnDim={100,0,700};
//
//
//    TH2D* exp2DHist=makeTH2DFromDim("expHist2D","expHist2D;E Peak Channel (1ch=40ns);BGO Energy Detected (keV)",BGO_DIM,PkChnDim);
//    exp2DHist->Sumw2();
//
//    for (int i = 0;i< 12;i++)
//    {
//        TString drawString="SBDEPkCh:BGOEDep"+int2str(i+1);
//        TH2D* tempHist=(TH2D*) makeExpEPGHist(false,"PID7", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "", -(i+1),BGO_DIM,PkChnDim);
//        exp2DHist->Add(tempHist);
//        delete tempHist;
//    }
//    exp2DHist->Scale(1/12.);
//
//    TString hist2DPath=HistsDir+"SBDEPkChn_BGOEDep_2DHist.root";
//    TFile outFile(hist2DPath,"RECREATE");
//    exp2DHist->Write("",TObject::kOverwrite);
//    outFile.Close();

//    TString hist2DPath=HistsDir+"SBDEPkChn_BGOEDep_2DHist.root";
//    TFile inpFile(hist2DPath,"READ");
//    TH2D* exp2DHist=(TH2D*) inpFile.Get("expHist2D");
//    exp2DHist->SetDirectory(0);
//    inpFile.Close();
//
//    TCanvas* c3 = new TCanvas("c3","c3",10,10,1280,1024);
//    exp2DHist->Draw("COLZ");
//    TString imagePath=string(GRAPHS_DIR)+"SBDEPkChn_BGOEDep2DGraph.png";
//    c3->SaveAs(imagePath);
//    exp2DHist->Draw("SCAT");
//    imagePath=string(GRAPHS_DIR)+"SBDEPkChn_BGOEDep_Scat2DGraph.png";
//    c3->SaveAs(imagePath);
//    delete c3;
//
//    TH1D* highEPkChn=exp2DHist->ProjectionX("highPkChnHist",88,100,"e");
//    TH1D* lowEPkChn=exp2DHist->ProjectionX("lowPkChnHist",1,87,"e");
//
//    highEPkChn->SetTitle("#geq 616 PkChn;BGO Energy Detected (keV);#frac{EPG}{EP}");
//    lowEPkChn->SetTitle("< 616 PkChn;BGO Energy Detected (keV);#frac{EPG}{EP}");
//
//    TH1* expHists[2]={highEPkChn,lowEPkChn};
//
//    plotExpVersusMCToImage(2,expHists,0,NULL,"BGO Energies for high/low Electron timing regions","Exp_BGOEDep_PkChnEComparison.png",true);

	///Electron Energies
//    HistDim PkChnDim={700,0,700};
//    TString hist2DPath=HistsDir+"SBDEPkChn_SBDEEDep_2DHist_WCuts.root";
////
////    TString drawString="SBDEPkCh:SBDEDepE";
////    TH2D* exp2DHist= (TH2D*) makeExpEPHist("PID8", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "",EDEPE_DIM,PkChnDim);
////    exp2DHist->SetTitle("expHist2D;Electron Energy Detected (keV);E Peak Channel (1ch=40ns)");
////
////
////    TFile outFile(hist2DPath,"RECREATE");
////    exp2DHist->Write("",TObject::kOverwrite);
////    outFile.Close();
//
//    TFile inpFile(hist2DPath,"READ");
//    TH2D* exp2DHist=(TH2D*) inpFile.Get("PID8hist");
//    exp2DHist->SetDirectory(0);
//    inpFile.Close();
//
////    TCanvas* c3 = new TCanvas("c3","c3",10,10,1280,1024);
////    exp2DHist->Draw("COLZ");
////    TString imagePath=string(GRAPHS_DIR)+"SBDEPkChn_SBDEEDep2DGraph.png";
////    c3->SaveAs(imagePath);
////    exp2DHist->Draw("SCAT");
////    imagePath=string(GRAPHS_DIR)+"SBDEPkChn_SBDEEDep_Scat2DGraph.png";
////    c3->SaveAs(imagePath);
////    delete c3;
//
//    TH1D* highEPkChn=exp2DHist->ProjectionX("highPkChnHist",600,700,"e");
//    TH1D* lowEPkChn=exp2DHist->ProjectionX("lowPkChnHist",1,600,"e");
//
//    double highIntegral=highEPkChn->Integral(1,800);
//    double lowIntegral=lowEPkChn->Integral(1,800);
//    TString titleString;
//
//    titleString="#geq 600 PkChns ("+d2str(highIntegral)+");Electron Energy Detected (keV);counts";
//    highEPkChn->SetTitle(titleString);
//    titleString="< 600 PkChns ("+d2str(lowIntegral)+");Electron Energy Detected (keV);counts";
//    lowEPkChn->SetTitle(titleString);
//
//    highEPkChn->GetYaxis()->SetRangeUser(1,2e5);
//
//    TH1* expHists[2]={highEPkChn,lowEPkChn};
//
//    plotExpVersusMCToImage(2,expHists,0,NULL,"Electron Energies for high/low Electron timing regions (with cuts )","Exp_SBDEEDep_PkChnEComparison_WCuts.png",true);
//
////    TH1D* ePkChnHist=exp2DHist->ProjectionY("ePkChnHist",1,700,"e");
////    cout << "1D Integral:"<< ePkChnHist->Integral(1,700) << endl;
////
////    ePkChnHist->SetTitle("Electron PkChn;Timing (channel);counts");
////    ePkChnHist->GetYaxis()->SetRangeUser(.1,1e7);
////
////    TH1* expHists[1]={ePkChnHist};
////
////    plotExpVersusMCToImage(1,expHists,0,NULL,"Electron Peak Channel Timing","Exp_PkChnE130728.png",true);
//
////    TString hist2DPath=HistsDir+"SBDEPkChn_SBDEEDep_2DHist.root";
////    TFile inpFile(hist2DPath,"READ");
////    TH2D* exp2DHist=(TH2D*) inpFile.Get("PID8hist");
////    exp2DHist->SetDirectory(0);
////    inpFile.Close();
////
////    TH1D* ePkChnHist=exp2DHist->ProjectionY("ePkChnHist",1,700,"e");
////    cout << "1D Integral:"<< ePkChnHist->Integral(1,700) << endl;
////
////    TH1* expHists[1]={ePkChnHist,ePkChnHistsWCuts};
////
////    plotExpVersusMCToImage(2,expHists,0,NULL,"Electron Peak Channel Timing","Exp_PkChnE130728.png",true);

///Convert Bryce's full electric field to root format - divide and conquer test
//    TString newFieldToTest="/media/mjbexternal/comsol/Ef4_solo.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_SBD_Mirror_BAPDElectricField_SuperFine130724_4.root";
//    ROOTField totalField(newFieldToTest,"total4",1,3,3);
//    totalField.saveFieldToFile(fieldToSaveAs,"total4");

//    TString newFieldToTest="/media/mjbexternal/comsol/Ef_Total.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_SBD_Mirror_BAPDElectricField_SuperFine130724.root";
//    ROOTField totalField(newFieldToTest,"total",1,3,3);
//    totalField.saveFieldToFile(fieldToSaveAs,"total");

///Check full new geant version and new electric field run 130807
//    TString histPath;
//    TH1D* mcOldHist=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID239","RID240","","STD NonRad W Fermi New Geant");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID241","RID242","RID243","STD NonRad W Fermi New Geant");
//    histPath=HistsDir+"MC_Hist130807_EDepG_AverageIndiv_EPG_RID243.txt";
//
////    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,BGO_DIM);
////    convertTH1ToTXT(mcHist,histPath);
//    TH1D* mcHist=getTabSeperatedHist(histPath,BGO_DIM);
//
//    RDK2Set* threeBodyNewSet=new RDK2Set(48,"130429","STD_3B","RID244","RID245","","STD NonRad W Fermi New Geant, New Fields");
//    RDK2Set* fourBodyNewSet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID246","RID247","RID248","STD NonRad W Fermi New Geant, New Fields");
//    histPath=HistsDir+"MC_Hist130807_EDepG_AverageIndiv_EPG_RID248.txt";
//
////    TH1* mcNewHist=makeMCIndividualDetectorHist(threeBodyNewSet,fourBodyNewSet,"BGOEDepBlurLOAllExpParam",STD_EPCut,BGO_DIM);
////    convertTH1ToTXT(mcHist,histPath);
//    TH1D* mcNewHist=getTabSeperatedHist(histPath,BGO_DIM);
//
//    TH1* mcHists[3]={mcOldHist,mcHist,mcNewHist};
//
//    mcOldHist->SetTitle("Old Geant/Old Fields;BGO Energy Detected (keV);#frac{EPG}{EP}");
//    mcHist->SetTitle("New Geant/Old Fields;BGO Energy Detected (keV);#frac{EPG}{EP}");
//    mcNewHist->SetTitle("New Geant/New Fields;BGO Energy Detected (keV);#frac{EPG}{EP}");
//
//    for (int i = 0;i< 3;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,3,mcHists,"MC Comparison","GeantVersionNewFieldsComparison.png");
//    plotResidualsToImage(2, &(mcHists[1]),1,mcHists,"#frac{New-Old}{Old}","GeantVersionNewFieldsComparison_Residuals.png");

///Multiplicity 2 to investigate 130820
//    TH1* expHist;
//    TH1* mcHist;
//    TH1* expHists[1];
//    TH1* mcHists[1];
//    TString histFileName;
//
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID244","RID245","","STD NonRad W Fermi New Geant, New Fields");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID246","RID247","RID248","STD NonRad W Fermi New Geant, New Fields");
//
//    const int numTypes=3;
//    TString expDrawStrings[numTypes]={"SPECIAL_EDepG","SBDEDepE","SBDToFP"};
//    TString mcDrawStrings[numTypes]={R_BGOSUM,"eDChn.SBDESBDEDepBlur","pChn.SBDTimeFirst"};
//    TString nameStrings[numTypes]={"BGOSum","EDepE","ToFP"};
//    TString titleStrings[numTypes]={"BGO Energy Detected (keV)","Electron Energy Detected (keV)","Proton Time of Flight (us)"};
//    TString dateString="130820";
//    HistDim histDims[3]={BGO_DIM,EDEPE_DIM,TOFP_DIM};
//
//    for (int i = 0;i< 1;i++)
//    {
//        histFileName="ExpMatt_Hists"+dateString+"_G7_"+nameStrings[i]+"_M2_EPG.txt";
//
//
////        expHist=makeExpEPGHist(false,"PID7", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, expDrawStrings[i], "SBDEDepE > 50. && SBDEDepE < 800. && SBDToFP >  2e-6 && SBDToFP < 10e-6 && SBDEDepP > 10. && SBDEDepP < 31.", 2,histDims[i]);
//        expHist=makeExpEPGHist(false,"PID7", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, expDrawStrings[i], "SBDEDepE > 50. && SBDEDepE < 800. && SBDEDepP > 10. && SBDEDepP < 31.", 2,histDims[i]);
////        expHist=makeExpEPGHist(false,"PID7", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, expDrawStrings[i], "", 2,histDims[i]);
//        convertTH1ToTXT(expHist,HistsDir+histFileName);
//
////        expHist=getTabSeperatedHist(HistsDir+histFileName,histDims[i]);
//
//        histFileName="MC_Hist"+dateString+"_"+nameStrings[i]+"_Mult2_EPG_RID248.txt";
//
////        mcHist=combineAndScaleRadSetsToHist( threeBodySet, fourBodySet,mcDrawStrings[i],STD_EPCut,"gDChn.BGOmultiplicityLOAllExpParam==2",histDims[i]);
////        convertTH1ToTXT(mcHist,HistsDir+histFileName);
//
//        mcHist=getTabSeperatedHist(HistsDir+histFileName,histDims[i]);
//
//        TString titleString="Exp;"+titleStrings[i]+";#frac{EPG}{EP}";
//        expHist->SetTitle(titleString);
//        titleString="MC;"+titleStrings[i]+";#frac{EPG}{EP}";
//        mcHist->SetTitle(titleString);
//
//        if(i==0 || i==1)
//        {
//            expHist->Rebin(10);
//            mcHist->Rebin(10);
//        }
//
//        expHists[0]=expHist;
//        mcHists[0]=mcHist;
//        TString imageFileName="Mult2Hist_"+nameStrings[i]+"_"+dateString+".png";
//        plotExpVersusMCToImage(1,expHists,1,mcHists,"Multiplicity 2",imageFileName);
//        imageFileName=addBeforeExtension(imageFileName,"_Residuals");
//        plotResidualsToImage(1,expHists,1,mcHists,"Multiplicity 2",imageFileName,-2,2);
//
//        delete expHist;
//        delete mcHist;
//    }
//
//    delete threeBodySet;
//    delete fourBodySet;

///Multiplicity 2 a different way to investigate 130822
//    TH1* expHist;
//    TH1* mcHist;
//    TH1* expHists[1];
//    TH1* mcHists[1];
//    TString histFileName;
//
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID244","RID245","","STD NonRad W Fermi New Geant, New Fields");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID246","RID247","RID248","STD NonRad W Fermi New Geant, New Fields");
//
//
//    TString nameString="BGOIndiv";
//    TString titleString="BGO Energy Detected (keV)";
//    TString dateString="130821";
//    HistDim histDim=BGO_DIM;
//
//
//    histFileName="ExpMatt_Hists"+dateString+"_G7_"+nameString+"_M2Indiv_EPG.txt";
//
////    expHist=makeTH1DFromDim("expAvg","",histDim);
////    for (int i = 0;i< 12;i++)
////    {
////    	TString expDrawString="(BGOEDep"+int2str(i+1)+" >10 && BGODeltaT"+int2str(i+1)+" > -75.01 && BGODeltaT"+int2str(i+1)+" < 74.99)*BGOEDep"+int2str(i+1)+"";
////    	expDrawString+=";(BGOEDep"+int2str(i+1)+" >10 && BGODeltaT"+int2str(i+1)+" > -450.01 && BGODeltaT"+int2str(i+1)+" < -300.99)*BGOEDep"+int2str(i+1)+"";
////
////    	TString expCutString="((BGOEDep1 >10 && BGODeltaT1 > -75.01 && BGODeltaT1 < 74.99) + (BGOEDep2 >10 && BGODeltaT2 > -75.01 && BGODeltaT2 < 74.99) + (BGOEDep3 >10 && BGODeltaT3 > -75.01 && BGODeltaT3 < 74.99) + (BGOEDep4 >10 && BGODeltaT4 > -75.01 && BGODeltaT4 < 74.99) + (BGOEDep5 >10 && BGODeltaT5 > -75.01 && BGODeltaT5 < 74.99) + (BGOEDep6 >10 && BGODeltaT6 > -75.01 && BGODeltaT6 < 74.99) + (BGOEDep7 >10 && BGODeltaT7 > -75.01 && BGODeltaT7 < 74.99) + (BGOEDep8 >10 && BGODeltaT8 > -75.01 && BGODeltaT8 < 74.99) + (BGOEDep9 >10 && BGODeltaT9 > -75.01 && BGODeltaT9 < 74.99)+ (BGOEDep10 >10 && BGODeltaT10 > -75.01 && BGODeltaT10 < 74.99) + (BGOEDep11 >10 && BGODeltaT11 > -75.01 && BGODeltaT11 < 74.99) + (BGOEDep12 >10 && BGODeltaT12 > -75.01 && BGODeltaT12 < 74.99)) == 2";
////    	expCutString +=" && SBDEDepE > 50. && SBDEDepE < 800. && SBDEDepP > 10. && SBDEDepP < 31.";
////    	expCutString +=";((BGOEDep1 >10 && BGODeltaT1 > -450.01 && BGODeltaT1 < -300.99) + (BGOEDep2 >10 && BGODeltaT2 > -450.01 && BGODeltaT2 < -300.99) + (BGOEDep3 >10 && BGODeltaT3 > -450.01 && BGODeltaT3 < -300.99) + (BGOEDep4 >10 && BGODeltaT4 > -450.01 && BGODeltaT4 < -300.99) + (BGOEDep5 >10 && BGODeltaT5 > -450.01 && BGODeltaT5 < -300.99) + (BGOEDep6 >10 && BGODeltaT6 > -450.01 && BGODeltaT6 < -300.99) + (BGOEDep7 >10 && BGODeltaT7 > -450.01 && BGODeltaT7 < -300.99) + (BGOEDep8 >10 && BGODeltaT8 > -450.01 && BGODeltaT8 < -300.99) + (BGOEDep9 >10 && BGODeltaT9 > -450.01 && BGODeltaT9 < -300.99)+ (BGOEDep10 >10 && BGODeltaT10 > -450.01 && BGODeltaT10 < -300.99) + (BGOEDep11 >10 && BGODeltaT11 > -450.01 && BGODeltaT11 < -300.99) + (BGOEDep12 >10 && BGODeltaT12 > -450.01 && BGODeltaT12 < -300.99)) == 2";
////    	expCutString +=" && SBDEDepE > 50. && SBDEDepE < 800. && SBDEDepP > 10. && SBDEDepP < 31.";
////
////        TH1* tempHist= makeExpEPGHist(false,"PID7", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, expDrawString, expCutString, -i,histDim);
////        expHist->Add(tempHist);
////        delete tempHist;
////    }
////    expHist->Scale(1./12.);
////    expHist->SetBinContent(1,0);
////    expHist->SetBinError(1,0);
////    convertTH1ToTXT(expHist,HistsDir+histFileName);
//
//
//    expHist=getTabSeperatedHist(HistsDir+histFileName,histDim);
//
//    histFileName="MC_Hist"+dateString+"_"+nameString+"_Mult2Indiv_EPG_RID248.txt";
//
////    TString mcDrawString="BGOEDepBlurLOAllExpParam";
////    mcHist=makeMCIndividualDetectorHist(threeBodySet, fourBodySet,mcDrawString,STD_EPCut,"gDChn.BGOmultiplicityLOAllExpParam == 2",histDim);
////    mcHist->Scale(12.);
////    convertTH1ToTXT(mcHist,HistsDir+histFileName);
//
//    mcHist=getTabSeperatedHist(HistsDir+histFileName,histDim);
//    mcHist->Scale(12.);
//
//
//    TString fullTitleString="Exp;"+titleString+";#frac{EPG}{EP}";
//    expHist->SetTitle(fullTitleString);
//    fullTitleString="MC;"+titleString+";#frac{EPG}{EP}";
//    mcHist->SetTitle(fullTitleString);
//
//
//    expHist->Rebin(10);
//    mcHist->Rebin(10);
//
//    mcHist->GetYaxis()->SetRangeUser(0,.8e-6);
//
//    expHists[0]=expHist;
//    mcHists[0]=mcHist;
//    TString imageFileName="Mult2IndivHist_"+nameString+"_"+dateString+".png";
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Multiplicity 2, Indiv",imageFileName);
//    imageFileName=addBeforeExtension(imageFileName,"_Residuals");
//    plotResidualsToImage(1,expHists,1,mcHists,"Multiplicity 2, Indiv",imageFileName,-2,2);
//
//    delete expHist;
//    delete mcHist;
//
//    delete threeBodySet;
//    delete fourBodySet;

///Convert Bryce BAPD field to ROOT format
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_BAPDElectricField_Fine130826.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_BAPDElectricField_Fine130826.root";
//    ROOTField totalField(newFieldToTest,"bapd",1,3,3);
//    totalField.saveFieldToFile(fieldToSaveAs,"bapd");

///Create macros for Geant4.9.6.p02 48 million full run of Livermore with new fields 130827 + new BAPD Field
//    TString standardLines="/MRK/det/setMirrorFieldFilePath COMSOL_MirrorElectricField_Fine130515.root";
//    standardLines+="\n/MRK/det/setSBDFieldFilePath COMSOL_SBDElectricField_Fine130515.root";
//    standardLines+="\n/MRK/det/setBAPDFieldScale 1350";
//
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID249",PROTON,THREEBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID250",ELECTRON,THREEBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID251",PROTON,FOURBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID252",ELECTRON,FOURBODY,standardLines);
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID253",GAMMA,FOURBODY,standardLines);

///Make event settings files to create energy deposit histograms per energy initial for deconvolution of exp spectrum 130827
//    TString  filePath=MRKRunDir+"4BRandomNumbers.txt";
//    ifstream randomNumbersFile(filePath);
//    cEvents theEvents;
//
//    for (int i = 0;i< 156;i++)
//    {
//        double e1=i*5;
//        double e2=(i+1)*5;
//        if(i==0)
//        {
//            e1=.1;
//        }
//
//        TString filePath=MRKRunDir+"EventSettings_STD_4B_"+d2str(e1)+"keV_"+d2str(e2)+"keV_130827.root";
//        TFile ranGenFile(filePath,"recreate");
//        ranGenFile.cd();
//
//
//        TVectorD fourBody(1);
//        TVectorD homogeneousEvents(1);
//        TVectorD littleb(1);
//        TVectorD normalizationConstant(1);
//        TVectorD minimumGammaEnergy(1);
//        TVectorD maximumGammaEnergy(1);
//
//        fourBody[0]=1;
//        homogeneousEvents[0]=0;
//        littleb[0]=0;
//        normalizationConstant[0]=theEvents.setGoodNormConstant(FOURBODY,0,e1,e2);
//
//        minimumGammaEnergy[0]=e1;
//        maximumGammaEnergy[0]=e2;
//
//        fourBody.Write("fourBody",TObject::kOverwrite);
//        homogeneousEvents.Write("homogeneousEvents",TObject::kOverwrite);
//        littleb.Write("littleb",TObject::kOverwrite);
//        normalizationConstant.Write("normalizationConstant",TObject::kOverwrite);
//        minimumGammaEnergy.Write("minimumGammaEnergy",TObject::kOverwrite);
//        maximumGammaEnergy.Write("maximumGammaEnergy",TObject::kOverwrite);
//
//        int randomSeed;
//        randomNumbersFile >> randomSeed;
//        TRandom3* ranGen=new TRandom3(randomSeed);
//        TString ranGenName="ranGen0";
//        ranGen->SetName(ranGenName);
//        ranGen->Write("",TObject::kOverwrite);
//        ranGenFile.Write("",TObject::kOverwrite);
//        ranGenFile.Close();
//
//    }
//    randomNumbersFile.close();
///Check EP rate for BAPD field run 130903
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID249","RID245","","STD NonRad W Fermi New Geant, New Fields with protons w BAPD field");
//    RDK2Set* threeBodySetOld=new RDK2Set(48,"130429","STD_3B","RID244","RID245","","STD NonRad W Fermi New Geant, New Fields with no BAPD field");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID246","RID247","RID248","STD NonRad W Fermi New Geant, New Fields");
//
//
//
//    double epRateNew = threeBodySet->calcEPRate(STD_EPCut);
//    double epRateOld = threeBodySetOld->calcEPRate(STD_EPCut);
//
//
//    double totalChange=(epRateNew-epRateOld)/epRateOld;
//
//    cout << "Total change (NEW-OLD)/OLD: " << totalChange << endl;
//    typeAnythingToContinue("");

///Check BAPD field
//    TString fieldPath="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130826b.root";
//    ROOTField bryceField(fieldPath,"bapd",1,3,3);
//
//    fieldPath="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130528.root";
//    ROOTField mattField(fieldPath,"bapd",1,3,3);
//
//
//    TVector3 vecOut(0,0,0);
//    bryceField.linearInterp3D(TVector3(0,0,.0321),vecOut);
//    cout << "Bryce field" << endl;
//    vecOut.Print();
//
//    vecOut.SetXYZ(0.,0.,0.);
//    cout << "Matt field" << endl;
//    mattField.linearInterp3D(TVector3(0,0,.032),vecOut);
//    vecOut.Print();
//
//    vecOut.SetXYZ(0,0,0);
//    bryceField.linearInterp3D(TVector3(.02399,-.024,-.04),vecOut);
//    cout << "Bryce field head" << endl;
//    vecOut.Print();
//
//    fieldPath="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_MirrorElectricField_Fine130515.root";
//    ROOTField mirrorField(fieldPath,"mirror",1,3,3);
//
//     vecOut.SetXYZ(0.,0.,0.);
//    cout << "Mirror field" << endl;
//    mirrorField.linearInterp3D(TVector3(.02399,-.02399,.220001),vecOut);
//    vecOut.Print();

///Convert Bryce BAPD field to ROOT format but do better
//    TString newFieldToTest="/home/mjbales/school/rdk/comsol/COMSOL_BAPDElectricField_Fine130826.txt";
//    TString fieldToSaveAs=MRKRunDir+"COMSOL_BAPDElectricField_Fine130826b.root";
//    ROOTField totalField(newFieldToTest,"bapd",1,3,3);
//    totalField.saveFieldToFile(fieldToSaveAs,"bapd");

///Calculate BAPD graph

	//Latest
//    RDK2Set* threeBodyNewSet=new RDK2Set(48,"130429","STD_3B","RID244","RID245","","STD NonRad W Fermi New Geant, New Fields");
//    RDK2Set* fourBodyNewSet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID246","RID247","RID248","STD NonRad W Fermi New Geant, New Fields");
//    TString histPath=HistsDir+"MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248.txt";

	//Fermi added, no geometry fix
//    RDK2Set* threeBodyNewSet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","");
//    RDK2Set* fourBodyNewSet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","");
//    TString histPath=HistsDir+"MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID175.txt";

//    //Geometry fix, w Fermi
//    RDK2Set* threeBodyNewSet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","");
//    RDK2Set* fourBodyNewSet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID186","");
//    TString histPath=HistsDir+"MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID186.txt";
//
//    mcHists[0]=makeMCIndividualDetectorHist(threeBodyNewSet,fourBodyNewSet,"BAPD",STD_EPCut,"",BAPD_DIM);
//    convertTH1ToTXT(mcHists[0],histPath);
////    mcHists[0]=getTabSeperatedHist(histPath,BAPD_DIM);
//    mcHists[1]=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BAPD_AverageIndiv_EPG_RID130.txt",BAPD_DIM);

//    const int numExpHists=2;
//    TH1* expHists[numExpHists];
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPBen_Hist121023_EDepG_BAPD_Sum_EPG.txt",BAPD_DIM);
//    expHists[1]=getTabSeperatedHist(HistsDir+"EXPBen_Hist130321_PID5_G6_EDepG_BAPD_AverageIndiv_EPG.txt",BAPD_DIM);
//
//    const int numMCHists=5;
//    const int numMCHists=1;
//    TH1* mcHists[numMCHists];
////    TString RIDList="RID107 RID130 RID175 RID186 RID248";
//    TString RIDList="RID248";
//    TString nameList="MC"
////    TString mcFileList="MC_Hist130116_EDepG_BAPD_EPG_Sum_RID107.txt MC_Hist130213_EDepG_BAPD_AverageIndiv_EPG_RID130.txt MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID175.txt MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID186.txt MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248.txt";
//    TString mcFileList="MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248.txt";
//    TString titleString;
//    double integral, integralError;
//
////    integral=expHists[0]->IntegralAndError(3,BAPD_DIM.numBins,integralError);
////    titleString= "Ben Old: "+d2str(integral,3)+" #pm "+d2str(integralError,1)+";BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
////    expHists[0]->SetTitle(titleString);
//
//    const int numExpHists=1;
//    integral=expHists[0]->IntegralAndError(3,BAPD_DIM.numBins,integralError);
//    titleString= "Ben New: "+d2str(integral,3)+" #pm "+d2str(integralError,1)+";BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
//    expHists[0]->SetTitle(titleString);
//
//    for (int i = 0;i< numMCHists;i++)
//    {
//        TString filePath;
//        stringFromList(mcFileList,i,filePath);
//        filePath=HistsDir+filePath;
//        mcHists[i]=getTabSeperatedHist(filePath,BAPD_DIM);
//
//        TString RIDString;
//        stringFromList(RIDList,i,RIDString);
//        if(RIDString=="RID107")
//        {
//            mcHists[i]->Scale(1/3.);
//        }
//        integral=mcHists[i]->IntegralAndError(3,BAPD_DIM.numBins,integralError);
//        titleString= "MC "+RIDString+": "+d2str(integral,3)+" #pm "+d2str(integralError,1)+";BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
//        mcHists[i]->SetTitle(titleString);
//
//    }
//
//    mcHists[0]->GetYaxis()->SetRangeUser(0,3e-6);
//
//    plotExpVersusMCToImage(numExpHists,expHists,numMCHists,mcHists,"Exp & Simulation Comparison: BAPD Energy Detected","EDep_BAPD_EPG_AverageIndiv_EXP_MC_Comparison_13010.root");
//
//
//    for (int i = 0;i< numMCHists;i++)
//    {
//    	delete mcHists[i];
//    }
//
//    for (int i = 0;i< numExpHists;i++)
//    {
//    	delete expHists[i];
//    }

	///Check EP rate for BAPD field run 130909
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID249","RID250","","STD NonRad W Fermi New Geant, New Fields with protons w BAPD field");
//    RDK2Set* threeBodySetOld=new RDK2Set(48,"130429","STD_3B","RID244","RID245","","STD NonRad W Fermi New Geant, New Fields with no BAPD field");
//    RDK2Set* fourBodySet=new RDK2Set(25,"130429","STD_4B_100eV_Endpoint","RID251","RID252","RID248","STD Rad W Fermi New Geant, New Fields");
//    RDK2Set* fourBodySetOld=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID246","RID247","RID248","STD Rad W Fermi New Geant, New Fields with no BAPD field");
//
//
////    double epRateNew = threeBodySet->calcEPRate(STD_EPCut);
////    double epRateOld = threeBodySetOld->calcEPRate(STD_EPCut);
////
////
////    double totalChange=(epRateNew-epRateOld)/epRateOld;it
////
////    cout << "Total change (NEW-OLD)/OLD: " << totalChange << endl;  //Answer: -6.2% 130909
//    HistDim channelDim={130,0,650};
//    TH1* mcHists[2];
////    mcHists[0]=threeBodySet->makeHist("newHist","pChn.SBDTimeFirst/40e-9",STD_ECut&&"pChn.SBDEDepTotal > 10. && pChn.SBDEDepTotal < 31.",channelDim);
////    mcHists[1]=threeBodySetOld->makeHist("oldHist","pChn.SBDTimeFirst/40e-9",STD_ECut&&"pChn.SBDEDepTotal > 10. && pChn.SBDEDepTotal < 31.",channelDim);
////
////
////
////    convertTH1ToTXT(mcHists[0],HistsDir+"MC_Hist130910_TOFP_EP_RID249.txt");
////    convertTH1ToTXT(mcHists[1],HistsDir+"MC_Hist130910_TOFP_EP_RID244.txt");
//
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130910_TOFP_EP_RID249.txt",channelDim);
//    mcHists[1]=getTabSeperatedHist(HistsDir+"MC_Hist130910_TOFP_EP_RID244.txt",channelDim);
//
//    double epRateNew = mcHists[0]->Integral();
//    double epRateOld = mcHists[1]->Integral();
//    double totalChange=(epRateNew-epRateOld)/epRateOld;
//    cout << "Total change (NEW-OLD)/OLD: " << totalChange << endl;  //Answer: -6.1% 130915
//    mcHists[0]->SetTitle("With BAPD Field;Proton Time of Flight (1 channel = 40ns);Counts (AU)");
//    mcHists[1]->SetTitle("Without BAPD Field;Proton Time of Flight (1 channel = 40ns);Counts (AU)");
//    mcHists[0]->GetXaxis()->SetRangeUser(0,650);
//    mcHists[0]->GetYaxis()->SetRangeUser(100,.7e6);
//
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"Effect of BAPD electric field on proton time of flight","BAPD_EField_Effect_On_PTOF_130910.png",true);
////    delete mcHists[0];
////    delete mcHists[1];
////
////    mcHists[0]=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,"",BGO_DIM);
////    convertTH1ToTXT(mcHists[0],HistsDir+"MC_Hist130915_EDepG_AverageIndiv_EPG_RID248b.txt");
////
////    mcHists[1]=getTabSeperatedHist(HistsDir+"MC_Hist130807_EDepG_AverageIndiv_EPG_RID248.txt",BGO_DIM);
////
////    mcHists[0]->SetTitle("With BAPD Field;BGO Energy Detected (keV);Counts (AU)");
////    mcHists[1]->SetTitle("Without BAPD Field;BGO Energy Detected (keV);Counts (AU)");
////
////    plotExpVersusMCToImage(0,NULL,2,mcHists,"Effect of BAPD electric field on EPG/EP","BAPD_EField_Effect_On_BGOEDepG_130915.png");
////    plotExpVersusMCToImage(0,NULL,2,mcHists,"Effect of BAPD electric field on EPG/EP","BAPD_EField_Effect_On_BGOEDepG_Residuals_130915.png");

///Plot BAPD field from Bryce

//    TString fieldPath="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_BAPDElectricField_Fine130826.root";
//    ROOTField bryceField(fieldPath,"bapd",1350,3,3);
//
//    fieldPath="/home/mjbales/school/rdk/mrkgeantrun/COMSOL_MirrorElectricField_Fine130515.root";
//    ROOTField mirrorField(fieldPath,"mirror",1400,3,3);
//
//    const int numPoints=1000;
//    double x[numPoints];
//    double y[numPoints];
//    double y2[numPoints];
//    TVector3 pos(0,0,0);
//    double length=.6;
//    double startingZ=-.1;
//    for (int i = 0;i< numPoints;i++)
//    {
//        x[i]=-.1+.6*i/numPoints;
//        pos.SetZ(x[i]);
//        TVector3 vecOut(0,0,0);
//        bryceField.linearInterp3D(pos,vecOut);
//        mirrorField.linearInterp3D(pos,vecOut);
//        y[i]=vecOut.Mag();
//        if(i!=0)
//        {
//            y2[i]=-vecOut.Z()*.6/numPoints+y2[i-1];
//        }
//
//
//    }
//
//
//
//
//    TGraph eFieldGraph(numPoints,x,y);
//    eFieldGraph.SetTitle("Total E-Field Magnitude of BAPD field and mirror field along beam axis; MC Z Coordinate (m);Total E-Field Magnitude (V/M)");
//
//    TCanvas* c2 = new TCanvas("graphCanvas","graphCanvas",10,10,1024,768);
//    eFieldGraph.Draw("AL");
//    TString imagePath=GraphsDir+"BAPD_Mirror_EFieldMag_130915.png";
//    c2->SaveAs(imagePath);
//
//    TGraph potentialGraph(numPoints,x,y2);
//    potentialGraph.SetTitle("Potential of BAPD field and mirror field along beam axis; MC Z Coordinate (m);Potential (V)");
//
//    potentialGraph.Draw("AL");
//    imagePath=GraphsDir+"BAPD_Mirror_Potential_130915.png";
//    c2->SaveAs(imagePath);

///Make BAPD Graph for Presentations

//    TH1* mcHists[1]={getTabSeperatedHist(HistsDir+"MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248.txt",BAPD_DIM)};
//
//
//    //    TH1* expHists[1]={getTabSeperatedHist(HistsDir+"EXPBen_Hist130321_PID5_G6_EDepG_BAPD_AverageIndiv_EPG.txt",BAPD_DIM)};
//    TH1* expHists[1]={getTabSeperatedHist(HistsDir+"EXPBen_Hist131009_BAPD_AverageIndiv_EPG.txt",BAPD_DIM)};
//
//
//    for (int i = 0;i< 2;i++)
//    {
//    	mcHists[0]->SetBinContent(i+1,0);
//    	mcHists[0]->SetBinError(i+1,0);
//    	expHists[0]->SetBinContent(i+1,0);
//    	expHists[0]->SetBinError(i+1,0);
//
//    }
//
//
//    double stdError;
//    double stdIntegral=expHists[0]->IntegralAndError(1,80,stdError,"");
//
//    double error;
//    double integral=mcHists[0]->IntegralAndError(1,80,error,"");
//
//    TString titleString="THEORY+MC: "+d2str(integral,3)+" #pm "+d2str(error,2)+";BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
////    TString titleString="THEORY+MC;BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
//    mcHists[0]->SetTitle(titleString);
//
//    titleString="EXP Ben: "+d2str(stdIntegral,3)+" #pm "+d2str(stdError,2)+";BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
////    titleString="EXP;BAPD Energy Detected (keV);EPG/EP (count rate ratio)";
//    expHists[0]->SetTitle(titleString);
//
//
//
//
//
//
//    mcHists[0]->Rebin(2);
//    expHists[0]->Rebin(2);
//
//
//    mcHists[0]->GetYaxis()->SetRangeUser(0,5e-6);
//    expHists[0]->GetYaxis()->SetRangeUser(0,5e-6);
//
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Bare APD Energy Detected","EDep_BAPD_EPG_AverageIndiv_EXP_MC_Comparison_131101.png");
//    mcHists[0]->SetTitle("MC");
//    expHists[0]->SetTitle("Exp");
//
//    plotResidualsToImage(1,expHists,1,mcHists,"Bare APD Energy Detected","EDep_BAPD_EPG_AverageIndiv_EXP_MC_Residuals_131101.png");
//
//
//
//    for (int i = 0;i< 1;i++)
//    {
//    	delete mcHists[i];
//    }
//
//    for (int i = 0;i< 1;i++)
//    {
//    	delete expHists[i];
//    }

///Make BGO Graph for Presentations
//    double integral,error;
//    TString titleString;
//    const int numMC=1;
//    const int numExp=2;
//    TH1* mcHists[numMC];
//    TH1* expHists[numExp];
//    TString xAxisTitle=";BGO Energy Detected (keV)";
//    TString yAxisTitle=";#frac{EPG}{EP} (count rate ratio)";
//
//
//
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//    integral=mcHists[0]->IntegralAndError(1,800,error,"");
//    titleString="THEORY+MC ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    mcHists[0]->SetTitle(titleString);
//
//    //    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
////
////
////    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    expHists[0]=getTabSeperatedHist(HistsDir+"ExpJeff_Hists130712_2_10_BGOEDepAvg_EPG.txt",BGO_DIM);
//    integral=expHists[0]->IntegralAndError(1,800,error,"");
//    titleString="EXP Jeff ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    expHists[0]->SetTitle(titleString);
//
//
//    expHists[1]=getTabSeperatedHist(HistsDir+"Ben_bgo_energy_20131009.txt",BGO_DIM);
//    integral=expHists[1]->IntegralAndError(1,800,error,"");
//    titleString="EXP Ben ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    expHists[1]->SetTitle(titleString);
////    expHists[1]->Scale(1./12.);
//
//
//
//    TH1* theHists[3]={expHists[0],expHists[1],mcHists[0]};
//    const int numRanges=6;
//    int lowEnd[numRanges] = {11,11,21,31,201,401};
//    int highEnd[numRanges] = {800,400,800,800,800,800};
//    TString nameString[3]={"Jeff","Ben","MC"};
//
//    for (int i = 0;i< 3;i++)
//    {
//        for (int j = 0;j< numRanges;j++)
//        {
//            integral = theHists[i]->IntegralAndError(lowEnd[j],highEnd[j],error,"");
//            cout << nameString[i] << " " << lowEnd[j]-1 << " " << highEnd[j] << " " << integral << " " << error << endl;
//
//        }
//
//    }
//
//    expHists[0]->Rebin(25);
//    expHists[1]->Rebin(25);
//    mcHists[0]->Rebin(25);
//
//    mcHists[0]->GetYaxis()->SetRangeUser(0,14e-6);
//    mcHists[0]->GetXaxis()->SetRangeUser(0,800);
//
//    plotExpVersusMCToImage(2,expHists,1,mcHists,"Average of Individual BGO Detectors","EDep_BGO_EPG_AverageIndiv_EXP_MC_Comparison_131014b.root");
//    expHists[0]->SetTitle("EXP Jeff");
//    expHists[1]->SetTitle("EXP Ben");
//    plotResidualsToImage(2,expHists,1,mcHists,"Average of Individual BGO Detectors","EDep_BGO_EPG_AverageIndiv_EXP_MC_NormalizedResidual_131014b.root");
//    plotResidualsToImage(2,expHists,1,mcHists,"Average of Individual BGO Detectors","EDep_BGO_EPG_AverageIndiv_EXP_MC_RawResidual_131014b.root",-2e-6,2e-6,false);
//
//
//    for (int i = 0;i< 1;i++)
//    {
//    	delete mcHists[i];
//    }
//
//    for (int i = 0;i< 2;i++)
//    {
//    	delete expHists[i];
//    }

///Check what's going wrong with newer MC files
//    TString histPath;
//    TH1D* mcOldHist=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID239","RID240","","STD NonRad W Fermi New Geant");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID241","RID242","RID243","STD NonRad W Fermi New Geant");
//    histPath=HistsDir+"MC_Hist130807_EDepG_AverageIndiv_EPG_RID243.txt";
//
//    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,BGO_DIM);
//    convertTH1ToTXT(mcHist,histPath);
//    TH1D* mcHist=getTabSeperatedHist(histPath,BGO_DIM);
//
//    RDK2Set* threeBodyNewSet=new RDK2Set(48,"130429","STD_3B","RID244","RID245","","STD NonRad W Fermi New Geant, New Fields");
//    RDK2Set* fourBodyNewSet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID246","RID247","RID248","STD NonRad W Fermi New Geant, New Fields");
//    histPath=HistsDir+"MC_Hist130807_EDepG_AverageIndiv_EPG_RID248.txt";
//
////    TH1* mcNewHist=makeMCIndividualDetectorHist(threeBodyNewSet,fourBodyNewSet,"BGOEDepBlurLOAllExpParam",STD_EPCut,BGO_DIM);
////    convertTH1ToTXT(mcHist,histPath);
//    TH1D* mcNewHist=getTabSeperatedHist(histPath,BGO_DIM);
//
//    TH1* mcHists[3]={mcOldHist,mcHist,mcNewHist};
//
//    mcOldHist->SetTitle("Old Geant/Old Fields;BGO Energy Detected (keV);#frac{EPG}{EP}");
//    mcHist->SetTitle("New Geant/Old Fields;BGO Energy Detected (keV);#frac{EPG}{EP}");
//    mcNewHist->SetTitle("New Geant/New Fields;BGO Energy Detected (keV);#frac{EPG}{EP}");
//
//    for (int i = 0;i< 3;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,3,mcHists,"MC Comparison","GeantVersionNewFieldsComparison131008.png");
//    plotResidualsToImage(2, &(mcHists[1]),1,mcHists,"#frac{New-Old}{Old}","GeantVersionNewFieldsComparison_Residuals131008.png");

///Check beam shifts

//    TH1* stdHist=getTabSeperatedHist(HistsDir+"MC_Hist130213_EDepG_BGO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//    double stdError;
//    double stdIntegral=stdHist->IntegralAndError(11,800,stdError,"");
//    delete stdHist;
//
//    const int numHists=9;
//    TString ridNumbers[numHists]={"RID130","RID135","RID140","RID145","RID150","RID155","RID160","RID165","RID170"};
//
//    for (int i = 0;i< numHists;i++)
//    {
//        TH1* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130218_EDepG_BGO_AverageIndiv_EPG_"+ridNumbers[i]+".txt",BGO_DIM);
//        if(mcHist != NULL)
//        {
//
//            double error;
//            double integral=mcHist->IntegralAndError(11,800,error,"");
//
//            delete mcHist;
//
//            double ratio=(stdIntegral-integral)/stdIntegral;
//            double ratioError=sqrt(stdError*stdError+error*error);
//            ratioError=ratio*sqrt(ratioError*ratioError/(ratio*ratio)+stdError*stdError/(stdIntegral*stdIntegral));
//
//            cout << ridNumbers[i] << ": " << ratio << " " << ratioError << endl;
//        }
//    }

///Check Nonlinearity
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist130710_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//     double stdError;
//    double stdIntegral=mcHist->IntegralAndError(11,800,stdError,"");
//    delete mcHist;
//
//    double error;
//    double integral=mcHistLow->IntegralAndError(11,800,error,"");
//
//    delete mcHistLow;
//
//    double ratio=(stdIntegral-integral)/stdIntegral;
//    double ratioError=sqrt(stdError*stdError+error*error);
//    ratioError=ratio*sqrt(ratioError*ratioError/(ratio*ratio)+stdError*stdError/(stdIntegral*stdIntegral));
//
//    cout << "Low: " << ratio << " " << ratioError << endl;
//
//    integral=mcHistHigh->IntegralAndError(11,800,error,"");
//
//    delete mcHistHigh;
//
//    ratio=(stdIntegral-integral)/stdIntegral;
//    ratioError=sqrt(stdError*stdError+error*error);
//    ratioError=ratio*sqrt(ratioError*ratioError/(ratio*ratio)+stdError*stdError/(stdIntegral*stdIntegral));
//
//    cout << "High: " << ratio << " " << ratioError << endl;

///Create BGO Spectrum with Cycles5-8 and then with Electron Energy Cut >100 keV 131030

///Cycles 3-8
//    TString titleString;
//    TString xAxisTitle=";Energy Detected (keV)";
//    TString yAxisTitle=";#Frac{EPG}{EP} Ratio";
//    double integral, error;
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD NonRad W Fermi");
//    const int numMCHists=1;
//    TH1* mcHists[numMCHists];
//    TH1* expHist[1];
//    TString histPath,drawString;
//
//
////    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,"",BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist131030_EDepGWLO_AverageIndiv_EPG_RID130.txt");
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist131030_EDepGWLO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    integral=mcHists[0]->IntegralAndError(1,800,error,"");
//    titleString="THEORY+MC ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    mcHists[0]->SetTitle(titleString);
//
//    for (int i = 0;i< numMCHists;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    histPath=HistsDir+"ExpMatt_Hists131030_G7_Cycle3_8_BGOEDepAvg_EPG.txt";
//
//
//    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
//    expHist[0]->Sumw2();
//    for (int i = 0;i< 12;i++)
//    {
//
//        drawString="BGOEDep"+int2str(i+1);
//        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, "", -(i+1),BGO_DIM);
//        expHist[0]->Add(tempHist);
//        delete tempHist;
//    }
//    expHist[0]->Scale(1/12.);
//    convertTH1ToTXT(expHist[0],histPath);
//
////    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    integral=expHist[0]->IntegralAndError(1,800,error,"");
//    titleString="EXP_Matt ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    expHist[0]->SetTitle(titleString);
//
//    expHist[0]->Rebin(10);
//
//
//    plotExpVersusMCToImage(1,expHist,1,mcHists,"Cycle 3-8 Comparison","Cycles38Comparison131030.png");
//    mcHists[0]->SetTitle("THEORY+MC");
//    expHist[0]->SetTitle("Exp_Matt");
//    plotResidualsToImage(1,expHist,1,mcHists,"Cycle 3-8 Comparison #frac{EXP-MC}{MC}","Cycles38Comparison131030_Residuals.png");

///Cycles 3-4
//    TString titleString;
//    TString xAxisTitle=";Energy Detected (keV)";
//    TString yAxisTitle=";#Frac{EPG}{EP} Ratio";
//    double integral, error;
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD NonRad W Fermi");
//    const int numMCHists=1;
//    TH1* mcHists[numMCHists];
//    TH1* expHist[1];
//    TString histPath,drawString;
//
//
////    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,"",BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist131030_EDepGWLO_AverageIndiv_EPG_RID130.txt");
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist131030_EDepGWLO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    integral=mcHists[0]->IntegralAndError(1,800,error,"");
//    titleString="THEORY+MC ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    mcHists[0]->SetTitle(titleString);
//
//    for (int i = 0;i< numMCHists;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    histPath=HistsDir+"ExpMatt_Hists131030_G7_Cycle3_4_BGOEDepAvg_EPG.txt";
//
//
//    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
//    expHist[0]->Sumw2();
//
//    for (int i = 0;i< 12;i++)
//    {
//
//        drawString="BGOEDep"+int2str(i+1);
//        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7_Cycle3_4, EXP_SERIES_ID_LIST_G7_Cycle3_4, drawString, "", -(i+1),BGO_DIM);
//        expHist[0]->Add(tempHist);
//        delete tempHist;
//    }
//    expHist[0]->Scale(1/12.);
//    convertTH1ToTXT(expHist[0],histPath);
//
////    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    integral=expHist[0]->IntegralAndError(1,800,error,"");
//    titleString="EXP_Matt ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    expHist[0]->SetTitle(titleString);
//
//    expHist[0]->Rebin(10);
//
//
//    plotExpVersusMCToImage(1,expHist,1,mcHists,"Cycle 3-4 Comparison","Cycles34Comparison131030.png");
//    mcHists[0]->SetTitle("THEORY+MC");
//    expHist[0]->SetTitle("Exp_Matt");
//    plotResidualsToImage(1,expHist,1,mcHists,"Cycle 3-4 Comparison #frac{EXP-MC}{MC}","Cycles34Comparison131030_Residuals.png");

///Cycles 5-8
//    TString titleString;
//    TString xAxisTitle=";Energy Detected (keV)";
//    TString yAxisTitle=";#Frac{EPG}{EP} Ratio";
//    double integral, error;
//    RDK2Set* threeBodySet=new RDK2Set(48,"130116","STD_3B","RID126","RID127","","STD NonRad W Fermi");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","RID128","RID129","RID130","STD NonRad W Fermi");
//    const int numMCHists=1;
//    TH1* mcHists[numMCHists];
//    TH1* expHist[1];
//    TString histPath,drawString;
//
//
////    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,"",BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist131030_EDepGWLO_AverageIndiv_EPG_RID130.txt");
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist131030_EDepGWLO_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    integral=mcHists[0]->IntegralAndError(1,800,error,"");
//    titleString="THEORY+MC ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    mcHists[0]->SetTitle(titleString);
//
//    for (int i = 0;i< numMCHists;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    histPath=HistsDir+"ExpMatt_Hists131030_G7_Cycle5_8_BGOEDepAvg_EPG.txt";
//
//
////    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
////    expHist[0]->Sumw2();
////
////    for (int i = 0;i< 12;i++)
////    {
////
////        drawString="BGOEDep"+int2str(i+1);
////        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7_Cycle5_8, EXP_SERIES_ID_LIST_G7_Cycle5_8, drawString, "", -(i+1),BGO_DIM);
////        expHist[0]->Add(tempHist);
////        delete tempHist;
////    }
////    expHist[0]->Scale(1/12.);
////    convertTH1ToTXT(expHist[0],histPath);
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    integral=expHist[0]->IntegralAndError(1,800,error,"");
//    titleString="EXP_Matt ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    expHist[0]->SetTitle(titleString);
//
//    expHist[0]->Rebin(10);
//
//
//    plotExpVersusMCToImage(1,expHist,1,mcHists,"Cycle 5-8 Comparison","Cycles58Comparison131030.png");
//    mcHists[0]->SetTitle("THEORY+MC");
//    expHist[0]->SetTitle("Exp_Matt");
//    plotResidualsToImage(1,expHist,1,mcHists,"Cycle 5-8 Comparison #frac{EXP-MC}{MC}","Cycles58Comparison131030_Residuals.png");

///Cycles 5-8 100 keV Cut

//
////   mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCutB,"",BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist131030_EDepGWLO_SBEDepE100_AverageIndiv_EPG_RID130.txt");
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist131030_EDepGWLO_SBEDepE100_AverageIndiv_EPG_RID130.txt",BGO_DIM);
//
//    integral=mcHists[0]->IntegralAndError(1,800,error,"");
//    titleString="THEORY+MC ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    mcHists[0]->SetTitle(titleString);
//
//
//    for (int i = 0;i< numMCHists;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//
//
//    histPath=HistsDir+"ExpMatt_Hists131030_G7_Cycle5_8_SBEDepE100_BGOEDepAvg_EPG.txt";
//
////    expHist[0]=makeTH1DFromDim("expHist2","expHist2",BGO_DIM);
////    expHist[0]->Sumw2();
////    for (int i = 0;i< 12;i++)
////    {
////
////        drawString="BGOEDep"+int2str(i+1);
////        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7_Cycle5_8, EXP_SERIES_ID_LIST_G7_Cycle5_8, drawString, "SBDEDepE > 100", -(i+1),BGO_DIM);
////        expHist[0]->Add(tempHist);
////        delete tempHist;
////    }
////    expHist[0]->Scale(1/12.);
////    convertTH1ToTXT(expHist[0],histPath);
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    integral=expHist[0]->IntegralAndError(1,800,error,"");
//    titleString="EXP_Matt ("+d2str(integral,3)+" #pm "+d2str(error,2)+")"+xAxisTitle+yAxisTitle;
//    expHist[0]->SetTitle(titleString);
//
//    expHist[0]->Rebin(10);
//
//
//    plotExpVersusMCToImage(1,expHist,1,mcHists,"Cycle 5-8 Comparison","Cycles58SBEDepE100Comparison131030.png");
//    mcHists[0]->SetTitle("THEORY+MC");
//    expHist[0]->SetTitle("Exp_Matt");
//    plotResidualsToImage(1,expHist,1,mcHists,"Cycle 5-8 Comparison #frac{EXP-MC}{MC}","Cycles58_SBEDepE100Comparison131030_Residuals.png");

///Check BR up to approximate
//  cEvents theEvents;
//    cout << theEvents.calcBR(14,EKEMAX) << endl; //0.0102123
//    cout << theEvents.calcBR(14,78) << endl;

///Checking SBD field
//    TString fieldToSaveAs="/home/mjbales/school/rdk/mrkgeantrun/Efield_RDK2_SBD_Matt5.root";
//    ROOTField sbdField(fieldToSaveAs,"sbd",1,2,2);
//
//    sbdField.linearInterp3D(

///Calc EPG Rate 131105
//    TCut PCutSlow=TCut("pChn.SBDEDepTotal > 10. && pChn.SBDEDepTotal < 31. && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 25.e-6");
//    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","");
//    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID186","");
////    fourBodySet->calcEPGRateFor3IndividualBAPDs(STD_EPCut,1);
////    fourBodySet->calcEPGRateFor12IndividualBGOs(STD_EPCut);
////    threeBodySet->calcEPRate(STD_EPCut);
//    threeBodySet->calcEPRate(STD_ECut && PCutSlow);
//    fourBodySet->calcEPGRateFor3IndividualBAPDs(STD_EPCut);
//    fourBodySet->calcEPGRateFor12IndividualBGOs(STD_ECut && PCutSlow);

///Test latest Geant4 version with changes to SBD cuts and simplified field classes
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID253",PROTON,THREEBODY,"");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID254",ELECTRON,THREEBODY,"");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID255",PROTON,FOURBODY,"");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID256",ELECTRON,FOURBODY,"");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID257",GAMMA,FOURBODY,"");

///Shift SBD around.
	//Shift SBD Detector by 1 cm towards bend
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID258",PROTON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 1 cm;/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID259",ELECTRON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 1 cm;/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID260",PROTON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 1 cm;/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID261",ELECTRON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 1 cm;/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID262",GAMMA,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 1 cm;/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");

	//Shift SBD Detector 1 cm away from bend
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID263",PROTON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm;/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID264",ELECTRON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm;/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID265",PROTON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm;/MRK/field/setFieldOffset -0.1650476059 0 -0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID266",ELECTRON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm;/MRK/field/setFieldOffset -0.1650476059 0 -0.9862856015 cm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID267",GAMMA,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm;/MRK/field/setFieldOffset -0.1650476059 0 -0.9862856015 cm");

	//Shift SBD Detector by 1 mm side to side
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID268",PROTON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset 1 0 0 mm;/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID269",ELECTRON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset 1 0 0 mm;/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID270",PROTON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 1 0 0 mm;/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID271",ELECTRON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 1 0 0 mm;/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID272",GAMMA,FOURBODY,"/MRK/det/setSBDDetectorBendOffset 1 0 0 mm;/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");

	//Shift SBD Detector -1 mm side to side
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID273",PROTON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset -1 0 0 mm;/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID274",ELECTRON,THREEBODY,"/MRK/det/setSBDDetectorBendOffset -1 0 0 mm;/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID275",PROTON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset -1 0 0 mm;/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID276",ELECTRON,FOURBODY,"/MRK/det/setSBDDetectorBendOffset -1 0 0 mm;/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");
	//makeGeantMacroFilesForClusterEventGenerator(0,1,"RID277",GAMMA,FOURBODY,"/MRK/det/setSBDDetectorBendOffset -1 0 0 mm;/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");

///Test Matt new fields in this version - 1.026
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID278",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID279",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID280",PROTON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID281",ELECTRON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID282",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");

/////Analyze test of radial originating point to detection radius down the bend - 131117
//
//    TString listOfResultFiles="radialTest20NoField.root radialTest20WField.root radialTest25NoField.root radialTest25WField.root radialTest27.8NoField.root radialTest27.8WField.root";
//    TString listOfEventFiles="Events_radialTest20NoField.root Events_radialTest20WField.root Events_radialTest25NoField.root Events_radialTest25WField.root Events_radialTest27.8NoField.root Events_radialTest27.8WField.root";
//    TString listOfDistances="20 20 25 25 27.8 27.8";
//    TString listOfEfield="No With No With No With";
//
//
//    for (int i = 0;i< 6;i++)
//    {
//        TString resultFileString;
//        stringFromList(listOfResultFiles,i,resultFileString);
//        TString resultsFilePath=MRKRunDir+resultFileString;
//        TFile resultsFile(resultsFilePath);
//        TTree* resultsTree=(TTree*) resultsFile.Get("g");
//
//        TString eventFileString;
//        stringFromList(listOfEventFiles,i,eventFileString);
//        TString eventFilePath=MRKRunDir+eventFileString;
//        TFile eventFile(eventFilePath);
//        TTree* eventTree=(TTree*) eventFile.Get("t");
//        resultsTree->AddFriend(eventTree,"t");
//
//        TCanvas c2("canvas","plotManyHistCanvas",1928,269,799,591);
//        c2.SetFillColor(kGray);
//        c2.cd();
//        c2.SetLeftMargin(.09);
//        c2.SetRightMargin(.05);
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
//        TH2F axisHist("axisHist","axisHist",2,0,0.02,2,0,0.026);
//        axisHist.GetYaxis()->SetRangeUser(0,0.026);
//        axisHist.GetXaxis()->SetRangeUser(0,0.02);
//        axisHist.Draw();
//
//
//
//        TString drawString="sqrt(pow(SBDPosZFirst*"+d2str(SIN_BEND)+"-SBDPosXFirst*"+d2str(COS_BEND)+",2)+SBDPosYFirst*SBDPosYFirst):sqrt(t.x0*t.x0+t.y0*t.y0)";
//        resultsTree->Draw(drawString,"SBDKEFirst > 0.1 && t.z0 > 0.03 && t.z0 < 0.40","same");
//
//        TH2F* htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes
//        axisHist.GetXaxis()->SetTitle("Emission Radius (m)");
//        axisHist.GetYaxis()->SetTitle("Detection Radius (m)");
//        TString distString,eFieldString;
//        stringFromList(listOfDistances,i,distString);
//        stringFromList(listOfEfield,i,eFieldString);
//        TString titleString=eFieldString+" -25 kV SBD Field: 3 cm < z0 < 40 cm @ SBD at "+distString+" cm from bend";
//        axisHist.SetTitle(titleString);
//        axisHist.GetYaxis()->SetTitleOffset(1.14);
//
//
//        TF1 silicon600mm("silicon600mm","0.01382", 0, 0.06);
//        silicon600mm.SetLineColor(kGreen);
//        silicon600mm.Draw("same");
//
//        TF1 case600mm("case600mm","0.01805", 0, 0.06);
//        case600mm.SetLineColor(kRed);
//        case600mm.Draw("same");
//
//        TF1 silicon300mm("silicon300mm","0.01382/sqrt(2)", 0, 0.06);
//        silicon300mm.SetLineColor(kCyan);
//        silicon300mm.Draw("same");
//
//        TF1 case300mm("case300mm","0.01805/sqrt(2)", 0, 0.06);
//        case300mm.SetLineColor(kOrange);
//        case300mm.Draw("same");
//
//        TF1 electrode("electrode","0.1", 0.01235, 0.01245);
//        electrode.SetLineColor(kBlue);
//        electrode.SetFillColor(kBlue);
//        electrode.SetFillStyle(1001);
//        electrode.Draw("f same");
//
//        TLegend theLegend(0.7,0.1,0.93,0.5);
//
//
//        theLegend.AddEntry(htemp,"Proton Events","P");
//        theLegend.AddEntry(&silicon600mm,"600 mm^{2} silicon","L");
//        theLegend.AddEntry(&case600mm,"600 mm^{2} case","L");
//        theLegend.AddEntry(&silicon300mm,"300 mm^{2} silicon","L");
//        theLegend.AddEntry(&case300mm,"300 mm^{2} case","L");
//        theLegend.AddEntry(&electrode,"Electrode inner radius","L");
//
//        theLegend.Draw();
//        theLegend.SetBorderSize(1);
//        theLegend.SetFillColor(kWhite);
//
//        TString imagePath=GraphsDir+resultFileString+".png";
//        c2.SaveAs(imagePath);
//        eventFile.Close();
//        resultsFile.Close();
//
//
//    }

///Recreate origin graph for RDK2 posting
///Analyze test of radial originating point to detection radius down the bend - 131117

//    TString listOfResultFiles="radialTest27.8NoField.root radialTest27.8WField.root";
//    TString listOfEventFiles="Events_radialTest27.8NoField.root Events_radialTest27.8WField.root";
//    TString listOfDistances="27.8 27.8";
//    TString listOfEfield="No With";
//
//
//    for (int i = 0;i< 2;i++)
//    {
//        TString resultFileString;
//        stringFromList(listOfResultFiles,i,resultFileString);
//        TString resultsFilePath=MRKRunDir+resultFileString;
//        TFile resultsFile(resultsFilePath);
//        TTree* resultsTree=(TTree*) resultsFile.Get("g");
//
//        TString eventFileString;
//        stringFromList(listOfEventFiles,i,eventFileString);
//        TString eventFilePath=MRKRunDir+eventFileString;
//        TFile eventFile(eventFilePath);
//        TTree* eventTree=(TTree*) eventFile.Get("t");
//        resultsTree->AddFriend(eventTree,"t");
//
//        TCanvas c2("canvas","plotManyHistCanvas",1928,269,799,591);
//        c2.SetFillColor(kGray);
//        c2.cd();
//        c2.SetLeftMargin(.09);
//        c2.SetRightMargin(.05);
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
//        TH2F axisHist("axisHist","axisHist",2,0,0.02,2,0,0.026);
//        axisHist.GetYaxis()->SetRangeUser(0,0.026);
//        axisHist.GetXaxis()->SetRangeUser(0,0.02);
//        axisHist.Draw();
//
//
//
//        TString drawString="sqrt(pow(SBDPosZFirst*"+d2str(SIN_BEND)+"-SBDPosXFirst*"+d2str(COS_BEND)+",2)+SBDPosYFirst*SBDPosYFirst):sqrt(t.x0*t.x0+t.y0*t.y0)";
//        resultsTree->Draw(drawString,"SBDKEFirst > 0.1 && t.z0 > 0.03 && t.z0 < 0.32","same");
//
//        TH2F* htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes
//        axisHist.GetXaxis()->SetTitle("Emission Radius (m)");
//        axisHist.GetYaxis()->SetTitle("Detection Radius (m)");
//        TString distString,eFieldString;
//        stringFromList(listOfDistances,i,distString);
//        stringFromList(listOfEfield,i,eFieldString);
//        TString titleString=eFieldString+" -25 kV SBD Field: 3 cm < z0 < 32 cm @ SBD at "+distString+" cm from bend";
//        axisHist.SetTitle(titleString);
//        axisHist.GetYaxis()->SetTitleOffset(1.14);
//
//
//        TF1 silicon600mm("silicon600mm","0.01382", 0, 0.06);
//        silicon600mm.SetLineColor(kGreen);
//        silicon600mm.Draw("same");
//
//        TF1 case600mm("case600mm","0.01805", 0, 0.06);
//        case600mm.SetLineColor(kRed);
//        case600mm.Draw("same");
//
//        TF1 bapdHolder("bapdHolder","0.1", 0.01745, 0.01755);
//        bapdHolder.SetLineColor(kBlue);
//        bapdHolder.SetFillColor(kBlue);
//        bapdHolder.SetFillStyle(1001);
//        bapdHolder.Draw("f same");
//
//        TF1 beam("beam","0.1", 0.01195, 0.01205);
//        beam.SetLineColor(kViolet);
//        beam.SetFillColor(kViolet);
//        beam.SetFillStyle(1001);
//        beam.Draw("f same");
//
//        TLegend theLegend(0.7,0.1,0.93,0.5);
//
//
//        theLegend.AddEntry(htemp,"Proton Events","P");
//        theLegend.AddEntry(&silicon600mm,"600 mm^{2} silicon OR","L");
//        theLegend.AddEntry(&case600mm,"600 mm^{2} case OR","L");
//        theLegend.AddEntry(&beam,"Rough Beam OR","L");
//        theLegend.AddEntry(&bapdHolder,"BAPD Holder IR","L");
//
//
//        theLegend.Draw();
//        theLegend.SetBorderSize(1);
//        theLegend.SetFillColor(kWhite);
//
//        TString imagePath=GraphsDir+resultFileString+"_RDK2.png";
//        c2.SaveAs(imagePath);
//        eventFile.Close();
//        resultsFile.Close();
//    }

///Make origin graph for RDK2 W/Beam
///Analyze test of radial originating point to detection radius down the bend - 131117

//    TString listOfResultFiles="radialTestRDK2Beam27.8NoField.root radialTestRDK2Beam27.8WField.root";
//    TString listOfEventFiles="Events_STD_3B_130429_Set0.root Events_STD_3B_130429_Set0.root";
//    TString listOfDistances="27.8 27.8";
//    TString listOfEfield="No With";
//
//
//    for (int i = 0;i< 1;i++)
//    {
//        TString resultFileString;
//        stringFromList(listOfResultFiles,i,resultFileString);
//        TString resultsFilePath=MRKRunDir+resultFileString;
//        TFile resultsFile(resultsFilePath);
//        TTree* resultsTree=(TTree*) resultsFile.Get("g");
//
//        TString eventFileString;
//        stringFromList(listOfEventFiles,i,eventFileString);
//        TString eventFilePath=MRKRunDir+eventFileString;
//        TFile eventFile(eventFilePath);
//        TTree* eventTree=(TTree*) eventFile.Get("t");
//        resultsTree->AddFriend(eventTree,"t");
//
//        TCanvas c2("canvas","plotManyHistCanvas",1928,269,799,591);
//        c2.SetFillColor(kGray);
//        c2.cd();
//        c2.SetLeftMargin(.09);
//        c2.SetRightMargin(.05);
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
////        TH2F axisHist("axisHist","axisHist",2,0,0.02,2,0,0.026);
////        axisHist.GetYaxis()->SetRangeUser(0,0.026);
////        axisHist.GetXaxis()->SetRangeUser(0,0.02);
////        axisHist.Draw();
//
//
//
//
//        TString drawString="sqrt(pow(SBDPosZFirst*"+d2str(SIN_BEND)+"-SBDPosXFirst*"+d2str(COS_BEND)+",2)+SBDPosYFirst*SBDPosYFirst):sqrt(t.x0*t.x0+t.y0*t.y0)";
////        TString drawString="0.0095:sqrt(t.x0*t.x0+t.y0*t.y0)";
////        TString drawString="SBDKEFirst:t.ep0";
//        resultsTree->Draw(drawString,"SBDKEFirst>0.1","");
//
//        TH2F* htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes
////        axisHist.GetXaxis()->SetTitle("Emission Radius (m)");
////        axisHist.GetYaxis()->SetTitle("Detection Radius (m)");
//        TString distString,eFieldString;
//        stringFromList(listOfDistances,i,distString);
//        stringFromList(listOfEfield,i,eFieldString);
//        TString titleString=eFieldString+" -25 kV SBD Field: From Beam @ SBD at "+distString+" cm from bend";
////        axisHist.SetTitle(titleString);
////        axisHist.GetYaxis()->SetTitleOffset(1.14);
//
//
//        TF1 silicon600mm("silicon600mm","0.01382", 0, 0.06);
//        silicon600mm.SetLineColor(kGreen);
//        silicon600mm.Draw("same");
//
//        TF1 case600mm("case600mm","0.01805", 0, 0.06);
//        case600mm.SetLineColor(kRed);
//        case600mm.Draw("same");
//
//        TF1 bapdHolder("bapdHolder","0.1", 0.01745, 0.01755);
//        bapdHolder.SetLineColor(kBlue);
//        bapdHolder.SetFillColor(kBlue);
//        bapdHolder.SetFillStyle(1001);
//        bapdHolder.Draw("f same");
//
//        TF1 beam("beam","0.1", 0.01195, 0.01205);
//        beam.SetLineColor(kViolet);
//        beam.SetFillColor(kViolet);
//        beam.SetFillStyle(1001);
//        beam.Draw("f same");
//
//        TLegend theLegend(0.7,0.1,0.93,0.5);
//
//
//        theLegend.AddEntry(htemp,"Proton Events","P");
//        theLegend.AddEntry(&silicon600mm,"600 mm^{2} silicon OR","L");
//        theLegend.AddEntry(&case600mm,"600 mm^{2} case OR","L");
//        theLegend.AddEntry(&beam,"Rough Beam OR","L");
//        theLegend.AddEntry(&bapdHolder,"BAPD Holder IR","L");
//
//
//        theLegend.Draw();
//        theLegend.SetBorderSize(1);
//        theLegend.SetFillColor(kWhite);
//
//        TString imagePath=GraphsDir+resultFileString+"_RDK2beam.png";
//        c2.SaveAs(imagePath);
//        eventFile.Close();
//        resultsFile.Close();
//    }

///make graphs of new all exp LO param model 131209
///First remaking to double check then temporarily changing LO curves so they all agree at lower energies
////    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
////    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD Rad W Fermi");
////    fourBodySet->makeDerivedResultsFiles(GAMMA);
//
////    TH1* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut&&"Entry$ % 4 == 0","",BGO_DIM);
////    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist131209_EDepG_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist131209_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
////    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut&&"Entry$ % 4 == 1","",BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist131209_EDepGWLO_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist131209_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
////    TH1* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",STD_EPCut&&"Entry$ % 4 == 2","",BGO_DIM);
////    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist131209_EDepGWLOLow_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist131209_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
////    TH1* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",STD_EPCut&&"Entry$ % 4== 3","",BGO_DIM);
////    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist131209_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist131209_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//
//    TString titleString;
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison131209.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp. Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison131209_Residuals.png");
//
//      TString drawString;
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt";
//
//
////    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
////    expHist[0]->Sumw2();
////
////    for (int i = 0;i< 12;i++)
////    {
////        TString cutString="SBDToFP < 10";
////        TCut EPCut=TCut(cutString);
////
////        drawString="SPECIAL_EDepG";
////        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, EPCut, -(i+1),BGO_DIM);
////        expHist[0]->Add(tempHist);
////        delete tempHist;
////    }
////    expHist[0]->Scale(1/12.);
////    convertTH1ToTXT(expHist[0],histPath);
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//
//    titleString=expHist[0]->GetTitle();
//    expHist[0]->SetTitle(titleString);
//
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison #frac{EXP-MC}{MC}","ParamModelComparison131209_ResidualsVsExp.png");
//     for (int i = 0;i< 4;i++)
//    {
//        titleString=mcHists[i]->GetTitle();
//        titleString+=": "+d2str(mcHists[i]->Integral(1,1000));
//        mcHists[i]->SetTitle(titleString);
//    }
//    titleString+=": "+d2str(expHist[0]->Integral(1,1000));
//    expHist[0]->SetTitle(titleString);
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison131209.png");

///Testing modified EDep formula
//    const int n=100;
//    double x[3][n];
//    double y[3][n];
//    EDepModel theModel[3]={ALLEXP_PARAM_MODEL,ALLEXP_PARAM_MODEL_HIGH,ALLEXP_PARAM_MODEL_LOW};
//    const int colorList[12]= {kBlack, 632, 416, 600, 616, 432,800, 820, 840, 860, 870,880};
//    TString modelString[3]={"Norm","High","Low"};
//
//    TGraph theGraphs[3];
//    TCanvas* c2 = new TCanvas("canvas","canvas",10,10,1024,768);
//    c2->SetFillColor(kGray);
//    c2->cd();
//    c2->SetLeftMargin(.12);
//    c2->SetRightMargin(.06);
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//    TLegend *theLegend;
//    theLegend = new TLegend(0.8,0.5,0.93,0.9);
//
//    for (int i = 0;i< 3;i++)
//    {
//        for (int j = 0;j< n;j++)
//        {
//            x[i][j]=(double) j+1;
//            y[i][j]=evalMattAllExpParamModel2(x[i][j], theModel[i]);
//        }
//
//
//    	theGraphs[i]=TGraph(n,x[i],y[i]);
//
//    	theGraphs[i].SetLineColor(colorList[i]);
//    	if(i==0)
//    	{
//    	    theGraphs[i].Draw("AL");
//    	    theGraphs[i].GetYaxis()->SetRangeUser(0.3,1.0);
//    	    theGraphs[i].SetTitle("Nonlinearity Model (Alternate);BGO Energy EDep (keV);Energy Relative to 662 keV");
//    	}
//    	else
//    	{
//    	    theGraphs[i].Draw("Lsame");
//    	}
//    	theLegend->AddEntry(&theGraphs[i],modelString[i],"L");
//    }
//
//    theLegend->Draw();
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
//    TString imagePath=GraphsDir+"ParamModelComparison131218alternate.png";
//    c2->SaveAs(imagePath);
//    delete c2;

///ok now with new files
////    RDK2Set* threeBodySet=new RDK2Set(48,"130429","STD_3B","RID171","RID172","","STD NonRad W Fermi");
////    RDK2Set* fourBodySet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","RID173","RID174","RID175","STD Rad W Fermi");
////    fourBodySet->makeDerivedResultsFiles(GAMMA);
//
////    TH1* mcHistOld=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurEdge",STD_EPCut&&"Entry$ % 4 == 0","",BGO_DIM);
////    convertTH1ToTXT(mcHistOld,HistsDir+"MC_Hist131209b_EDepG_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist131209b_EDepG_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
////
////    TH1* mcHist=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut&&"Entry$ % 4 == 1","",BGO_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist131209b_EDepGWLO_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist131209b_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
////    TH1* mcHistLow=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamLow",STD_EPCut&&"Entry$ % 4 == 2","",BGO_DIM);
////    convertTH1ToTXT(mcHistLow,HistsDir+"MC_Hist131209b_EDepGWLOLow_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistLow=getTabSeperatedHist(HistsDir+"MC_Hist131209b_EDepGWLOLow_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
////    TH1* mcHistHigh=makeMCIndividualDetectorHist(threeBodySet,fourBodySet,"BGOEDepBlurLOAllExpParamHigh",STD_EPCut&&"Entry$ % 4== 3","",BGO_DIM);
////    convertTH1ToTXT(mcHistHigh,HistsDir+"MC_Hist131209b_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt");
//    TH1D* mcHistHigh=getTabSeperatedHist(HistsDir+"MC_Hist131209b_EDepGWLOHigh_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//
//
//
//    TH1* mcHists[4] = {mcHistOld,mcHist,mcHistLow,mcHistHigh};
//
//    mcHists[0]->SetTitle("Old Khodyuk Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Param Model;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[2]->SetTitle("New Param Model - 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[3]->SetTitle("New Param Model + 3%;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< 4;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison131209b.png");
//    plotResidualsToImage(3,&(mcHists[1]),1,mcHists,"All Exp. Param Model Comparison #frac{X-OLD}{OLD}","ParamModelComparison131209b_Residuals.png");
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt";
//
//
////    expHist[0]=makeTH1DFromDim("expHist","expHist",BGO_DIM);
////    expHist[0]->Sumw2();
////
////    for (int i = 0;i< 12;i++)
////    {
////        TString cutString="SBDToFP < 10";
////        TCut EPCut=TCut(cutString);
////
////        drawString="SPECIAL_EDepG";
////        TH1* tempHist=makeExpEPGHist(false,"PID5", NUM_EXP_SERIES_G7, EXP_SERIES_ID_LIST_G7, drawString, EPCut, -(i+1),BGO_DIM);
////        expHist[0]->Add(tempHist);
////        delete tempHist;
////    }
////    expHist[0]->Scale(1/12.);
////    convertTH1ToTXT(expHist[0],histPath);
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//
//    TString titleString;
//    plotResidualsToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison #frac{EXP-MC}{MC}","ParamModelComparison131209_ResidualsVsExpb.png");
//     for (int i = 0;i< 4;i++)
//    {
//        titleString=mcHists[i]->GetTitle();
//        titleString+=": "+d2str(mcHists[i]->Integral(1,1000));
//        mcHists[i]->SetTitle(titleString);
//    }
//    titleString+=": "+d2str(expHist[0]->Integral(1,1000));
//    expHist[0]->SetTitle(titleString);
//    plotExpVersusMCToImage(1,expHist,4,mcHists,"All Exp. Param Model Comparison","ParamModelComparison131209b.png");

///Running files done in 11/11/13
//int beginningRID=254;
//int numberOfSets=5;
//
//for (int i = 0;i< numberOfSets;i++)
//{
//    createPlotOfRIDSequence(1,beginningRID+i*numberOfSets);
//}

///Figuring out why I had issues with figuring out beam radii to detection radii
///131228

//    TString listOfResultFiles="radialTestRDK2Beam27.8NoField.root radialTestRDK2Beam27.8WField.root";
//    TString listOfEventFiles="Events_STD_3B_130429_Set0.root Events_STD_3B_130429_Set0.root";
//    TString listOfDistances="27.8 27.8";
//    TString listOfEfield="No With";
//
//
//    for (int i = 0;i< 1;i++)
//    {
//        TString resultFileString;
//        stringFromList(listOfResultFiles,i,resultFileString);
//        TString resultsFilePath=MRKRunDir+resultFileString;
//        TFile resultsFile(resultsFilePath);
//        TTree* resultsTree=(TTree*) resultsFile.Get("g");
//
//        TString eventFileString;
//        stringFromList(listOfEventFiles,i,eventFileString);
//        TString eventFilePath=MRKRunDir+eventFileString;
//        TFile eventFile(eventFilePath);
//        TTree* eventTree=(TTree*) eventFile.Get("t");
//        resultsTree->AddFriend(eventTree,"t");
//
//        TCanvas c2("canvas","plotManyHistCanvas",1928,269,799,591);
//        c2.SetFillColor(kGray);
//        c2.cd();
//        c2.SetLeftMargin(1.8);
//        c2.SetRightMargin(.05);
//        gPad->SetTickx(1);
//        gPad->SetTicky(1);
//        gPad->SetFillColor(kWhite);
//
//        TH2F axisHist("axisHist","axisHist",2,0,0.022,2,0,0.022);
//        axisHist.Draw();
//
//
//
//
//        TString drawString="sqrt(pow(SBDPosZFirst*"+d2str(SIN_BEND)+"-SBDPosXFirst*"+d2str(COS_BEND)+",2)+SBDPosYFirst*SBDPosYFirst):sqrt((t.x0*t.x0)+(t.y0*t.y0))";
////        TString drawString="0.0095:sqrt(t.x0*t.x0+t.y0*t.y0)";
////        TString drawString="SBDKEFirst:t.ep0";
//        resultsTree->Draw(drawString,"SBDKEFirst>0.1","same");
//
////        TH2F* htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes
////        htemp->GetXaxis()->SetTitle("Emission Radius (m)");
////        htemp->GetYaxis()->SetTitle("Detection Radius (m)");
////        htemp->GetYaxis()->SetRangeUser(0,0.022);
////        htemp->GetXaxis()->SetRangeUser(0,0.022);
//        TString distString,eFieldString;
//        stringFromList(listOfDistances,i,distString);
//        stringFromList(listOfEfield,i,eFieldString);
//        TString titleString=eFieldString+" -25 kV SBD Field: From Beam @ SBD at "+distString+" cm from bend";
//        axisHist.SetTitle(titleString);
//        axisHist.GetYaxis()->SetTitleOffset(1.5);
//        axisHist.GetXaxis()->SetTitle("Emission Radius (m)");
//        axisHist.GetYaxis()->SetTitle("Detection Radius (m)");
//        axisHist.GetYaxis()->SetRangeUser(0,0.022);
//        axisHist.GetXaxis()->SetRangeUser(0,0.022);
//
//
//        TF1 silicon600mm("silicon600mm","0.01382", 0, 0.06);
//        silicon600mm.SetLineColor(kGreen);
//        silicon600mm.Draw("same");
//
//        TF1 case600mm("case600mm","0.01805", 0, 0.06);
//        case600mm.SetLineColor(kRed);
//        case600mm.Draw("same");
//
//        TF1 bapdHolder("bapdHolder","0.1", 0.01745, 0.01755);
//        bapdHolder.SetLineColor(kBlue);
//        bapdHolder.SetFillColor(kBlue);
//        bapdHolder.SetFillStyle(1001);
//        bapdHolder.Draw("f same");
//
//        TF1 beam("beam","0.1", 0.01195, 0.01205);
//        beam.SetLineColor(kViolet);
//        beam.SetFillColor(kViolet);
//        beam.SetFillStyle(1001);
//        beam.Draw("f same");
//
//        TLegend theLegend(0.7,0.1,0.93,0.5);
//
//
////        theLegend.AddEntry(htemp,"Proton Events","P");
//        theLegend.AddEntry(&axisHist,"Proton Events","P");
//        theLegend.AddEntry(&silicon600mm,"600 mm^{2} silicon OR","L");
//        theLegend.AddEntry(&case600mm,"600 mm^{2} case OR","L");
//        theLegend.AddEntry(&beam,"Rough Beam OR","L");
//        theLegend.AddEntry(&bapdHolder,"BAPD Holder IR","L");
//
//
//        theLegend.Draw();
//        theLegend.SetBorderSize(1);
//        theLegend.SetFillColor(kWhite);
//
//        TString imagePath=GraphsDir+resultFileString+"_RDK2beam131228.png";
//        c2.SaveAs(imagePath);
//        eventFile.Close();
//        resultsFile.Close();
//    }

	///Test latest Geant4 version with correction of beam radii 1.027 131230
//    makeGeantMacroFilesForClusterEventGenerator(1,47,"RID253",PROTON,THREEBODY,"/MRK/field/addNamedFields OLD_RDK2");
//    makeGeantMacroFilesForClusterEventGenerator(1,47,"RID254",ELECTRON,THREEBODY,"/MRK/field/addNamedFields OLD_RDK2");
//    makeGeantMacroFilesForClusterEventGenerator(1,47,"RID255",PROTON,FOURBODY,"/MRK/field/addNamedFields OLD_RDK2");
//    makeGeantMacroFilesForClusterEventGenerator(1,47,"RID256",ELECTRON,FOURBODY,"/MRK/field/addNamedFields OLD_RDK2");
//    makeGeantMacroFilesForClusterEventGenerator(1,47,"RID257",GAMMA,FOURBODY,"/MRK/field/addNamedFields OLD_RDK2");

///Test field settings (one last time?) 140102
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID278",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID279",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID280",PROTON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID281",ELECTRON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID282",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");

///Check latest Geant4 version with correction of beam radii
//   createPlotOfRIDSequence(48,253);
//   createPlotOfRIDSequence(48,171);

//    const int numMCHist=2;
//    TH1D* mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist131209_EDepGWLO_AverageIndiv_EPG_RID175.txt",BGO_DIM);
//    TH1D* mcHistNew=getTabSeperatedHist(HistsDir+"MC_Hist131230_EDepGWLO_AverageIndiv_EPG_RID257.txt",BGO_DIM);
////
//
//
//    TH1* mcHists[numMCHist] = {mcHistOld,mcHistNew};
//
//    mcHists[0]->SetTitle("Small Beam;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("Correct Beam;Energy Detected (keV);#Frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< numMCHist;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,numMCHist,mcHists,"Beam size correction","BeamSizeCorrection131230b.png");
//    plotResidualsToImage(numMCHist-1,&(mcHists[1]),1,mcHists,"Beam size correction #frac{X-OLD}{OLD}","BeamSizeCorrection131230b_Residuals.png");
//
//    TH1* expHist[1];
//    TString histPath=HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt";
//
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//    expHist[0]->Rebin(10);
//    expHist[0]->SetTitle("expHist");
//
//    TString titleString;
//    plotResidualsToImage(1,expHist,numMCHist,mcHists,"Beam size correction","BeamSizeCorrection131230b_ResidualsVsExpb.png");
//     for (int i = 0;i< numMCHist;i++)
//    {
//        titleString=mcHists[i]->GetTitle();
//        titleString+=": "+d2str(mcHists[i]->Integral(1,1000));
//        mcHists[i]->SetTitle(titleString);
//    }
//    titleString+=": "+d2str(expHist[0]->Integral(1,1000));
//    expHist[0]->SetTitle(titleString);
//    plotExpVersusMCToImage(1,expHist,numMCHist,mcHists,"Beam size correction","BeamSizeCorrection131230b.png");
//
//    delete mcHistOld;
//    delete mcHistNew;
//
//
//
//
//
//
//
/////Create plot for new fields
//  //  createPlotOfRIDSequence(48,278);
//
//
//    mcHistOld=getTabSeperatedHist(HistsDir+"MC_Hist131230_EDepGWLO_AverageIndiv_EPG_RID257.txt",BGO_DIM);
//    mcHistNew=getTabSeperatedHist(HistsDir+"MC_Hist140113_EDepGWLO_AverageIndiv_EPG_RID282.txt",BGO_DIM);
//
//
//
//    mcHists[0] = mcHistOld;
//    mcHists[1] =mcHistNew;
//
//    mcHists[0]->SetTitle("Old Field, Correct Beam;Energy Detected (keV);#frac{EPG}{EP} Ratio");
//    mcHists[1]->SetTitle("New Field, Correct Beam;Energy Detected (keV);#frac{EPG}{EP} Ratio");
//
//    for (int i = 0;i< numMCHist;i++)
//    {
//    	mcHists[i]->Rebin(10);
//    }
//
//    plotExpVersusMCToImage(0,NULL,numMCHist,mcHists,"Beam size correction","NewFieldCheck140113.png");
//    plotResidualsToImage(numMCHist-1,&(mcHists[1]),1,mcHists,"Beam size correction #frac{X-OLD}{OLD}","NewFieldCheck140113_Residuals.png");
//
//    plotResidualsToImage(1,expHist,numMCHist,mcHists,"Beam size correction","NewFieldCheck140113_ResidualsVsExp.png");
//     for (int i = 0;i< numMCHist;i++)
//    {
//        titleString=mcHists[i]->GetTitle();
//        titleString+=": "+d2str(mcHists[i]->Integral(1,1000));
//        mcHists[i]->SetTitle(titleString);
//    }
//    titleString+=": "+d2str(expHist[0]->Integral(1,1000));
//    expHist[0]->SetTitle(titleString);
//    plotExpVersusMCToImage(1,expHist,numMCHist,mcHists,"Beam size correction","NewFieldCheck140113vsExp.png");

///Analyze no PLO in LO extraction for Ba133 (Center3 = no PLO)
//    createDerivedResultsFilesForCal(10,"Ba133_Center3_140113");
//
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Center3_140113");
//    analyzeLOofCalPeaks("C461", "Center3_140113", true,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C439", "Center3_140113", true,CALPEAKTIMINGCUT_OFFPEAK);

///Make some blank thesis plots
//createBlankThesisPlot("tp_rdk1edetg_0.pdf", "RDK1 Gamma Energy Spectrum;Energy Detected (keV);Counts", "Exp. Sim.");
//createBlankThesisPlot("tp_ploresponse_0.pdf", "BGO Scintillation Based on Incident Position;Position (m);Relative Energy Detected", "BGO1 BGO6");
//createBlankThesisPlot("tp_bgonlexp_0.pdf", "BGO Non-proportionality Measurements;Energy Incident (keV);Energy Detected Relative to 662 keV", "Khodyuk Moszynski Gentile Averkiev Fit");
//createBlankThesisPlot("tp_beam_0.pdf", "Density Plot of Beam Flux;x (mm);y (mm)", "DensityKey");
//createBlankThesisPlot("tp_calbgo_0.pdf", "Calibration Source Example: Ba-133+Cs-137;Scintillation Counts (Channel);Counts", "Ba-133+Cs-137");
//createBlankThesisPlot("tp_calbapd_0.pdf", "Calibration Examples from Lightsource;BAPD Channel;Counts", "100eV 300eV 500eV 1000eV");
//createBlankThesisPlot("tp_bapdnlexp_0.pdf", "BAPD Non-proportionality Model for BAPDs;Depth (um);Relative Response", "BAPD1 BAPD2 BAPD3");
//createBlankThesisPlot("tp_bfieldcenter_0.pdf", "Magnetic Field Along Centerline;Position Along Centerline (m);Magnetic Field Strength (T)", "Measured Simulation");
//createBlankThesisPlot("tp_efieldcenter_0.pdf", "Electric Field Along Centerline;Position Along Centerline (m);Electric Field Strength (T)", "Measured Simulation");
//createBlankThesisPlot("tp_ep0_0.pdf", "Initial Proton Event Energy; Energy (eV);Events (AU)", "Threebody Fourbody");
//createBlankThesisPlot("tp_ee0_0.pdf", "Initial Electron Event Energy; Energy (eV);Events (AU)", "Threebody Fourbody");
//createBlankThesisPlot("tp_eg0_0.pdf", "Initial Gamma Event Energy; Energy (eV);Events (AU)", "Fourbody");
//createBlankThesisPlot("tp_tofp_0.pdf", "Electron - Proton Timing;Time Difference (us);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_edepp_0.pdf", "Proton Energy Desposited;Energy Desposited (keV);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_edepe_0.pdf", "Electron Energy Desposited;Energy Desposited (keV);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_bgoedepgfull_0.pdf", "Full BGO Gamma Energy Desposited;Energy Desposited (keV);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_bgoedepglow_0.pdf", "Low BGO Gamma Energy Desposited;Energy Desposited (keV);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_bgoedepgresid_0.pdf", "Residuals BGO Gamma Energy Desposited #frac{Exp-MC}{MC};Energy Desposited (keV);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_bapdedepgfull_0.pdf", "Full BAPD Gamma Energy Desposited;Energy Desposited (keV);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_bapdedepglow_0.pdf", "Low BAPD Gamma Energy Desposited;Energy Desposited (keV);#frac{EPG}{EP} Counts", "Measured Simulation");
//createBlankThesisPlot("tp_bapdedepgresid_0.pdf", "Residuals BGO Gamma Energy Desposited #frac{Exp-MC}{MC};#frac{EPG}{EP} Counts", "Measured Simulation");

///Make BGO length contraction gamma run 140120
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID283",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");

///Calc MC hit ratios for flux calcs 140120
//    RDK2Set threeBodySet(48,"130429","STD_3B","RID278","RID279","","STD NonRad W Fermi");
//    RDK2Set fourBodySet(48,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID282","STD Rad W Fermi");
//
//
//
//
////    threeBodySet.calcEPRate(STD_EPCut);
//    fourBodySet.calcEPGRateFor12IndividualBGOs(STD_EPCut); //0.0149523
//    fourBodySet.calcEPGRateFor3IndividualBAPDs(STD_EPCut);
//    fourBodySet.calcEPGRateFor3IndividualBAPDs(STD_EPCut,1);

///Check length contraction results 140122
//    RDK2Set threeBodySet(48,"130429","STD_3B","RID278","RID279","","STD NonRad W Fermi");
//    RDK2Set fourBodySet(48,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID283","STD Rad W Fermi");

//    fourBodySet.makeDerivedResultsFiles(GAMMA);

//    fourBodySet.calcEPGRateFor12IndividualBGOs(STD_EPCut); //0.0148503

///Interpolation checking
//    ROOTField bField;
//    bField.loadFieldFromFile(MRKRunDir+"Bfield_RDK2_Matt2.root", "solenoid", 1, 3,  3);
//    bField.setSymmetry(false,true,false);
//
//    TVector3 outVect(0,0,0);
//    for (int i = 0;i< 10;i++)
//    {
//        TVector3 pos(gRandom->Rndm()*bField.getXLength()+bField.getXStart(),gRandom->Rndm()*bField.getYLength()+bField.getYStart(),gRandom->Rndm()*bField.getZLength()+bField.getZStart());
//        pos.Print();
//        cout << endl;
//        bField.linearInterp3D(pos,outVect);
//        cout << "Linear result: " << endl;
//        outVect.Print();
//        outVect.SetXYZ(0,0,0);
//
//        bField.cubicInterp3D(pos,outVect);
//        cout << "Cubic result: " << endl;
//        outVect.Print();
//        outVect.SetXYZ(0,0,0);
//
//        bField.cubicInterp3D2(pos,outVect);
//        cout << "Cubic Other result: " << endl;
//        outVect.Print();
//        outVect.SetXYZ(0,0,0);
//        cout << endl << endl;
//
//    }

//    ROOTField testField;
//    testField.loadFieldFromFile(MRKRunDir+"Field_Test.txt", "", 1, 3,  3);
//
//    TVector3 outVect(0,0,0);
//    for (int i = 0;i< 10;i++)
//    {
//        TVector3 pos(gRandom->Rndm()*testField.getXLength()+testField.getXStart(),gRandom->Rndm()*testField.getYLength()+testField.getYStart(),gRandom->Rndm()*testField.getZLength()+testField.getZStart());
//        pos.Print();
//        cout << endl;
//        testField.linearInterp3D(pos,outVect);
//        cout << "Linear result: " << endl;
//        outVect.Print();
//        outVect.SetXYZ(0,0,0);
//
//        testField.cubicInterp3D(pos,outVect);
//        cout << "Cubic result: " << endl;
//        outVect.Print();
//        outVect.SetXYZ(0,0,0);
//
//        testField.cubicInterp3D2(pos,outVect);
//        cout << "Cubic Other result: " << endl;
//        outVect.Print();
//        outVect.SetXYZ(0,0,0);
//        cout << endl << endl;
//    }
//
//    TVector3 pos(1.5,1.5,1.5);
//    pos.Print();
//    cout << endl;
//    testField.linearInterp3D(pos,outVect);
//    cout << "Linear result: " << endl;
//    outVect.Print();
//    outVect.SetXYZ(0,0,0);
//
//    testField.cubicInterp3D(pos,outVect);
//    cout << "Cubic result: " << endl;
//    outVect.Print();
//    outVect.SetXYZ(0,0,0);
//
//    testField.cubicInterp3D2(pos,outVect);
//    cout << "Cubic Other result: " << endl;
//    outVect.Print();
//    outVect.SetXYZ(0,0,0);
//    cout << endl << endl;

///Test timing differences for linear, old cubic, and new cubic interpolation
//    ROOTField bField;
//    bField.loadFieldFromFile(MRKRunDir+"Bfield_RDK2_Matt2.root", "solenoid", 1, 3,  3);
//    bField.setSymmetry(false,true,false);
//
//    int numEvents=5000000;
//
//    time_t startTime, endTime;
//    TVector3 newFieldOut,pos;
//
//    //Linear
//    startTime = time(NULL);
//    for (int i = 0;i< numEvents;i++)
//    {
//        newFieldOut.SetXYZ(0,0,0);
//
//        pos.SetX(gRandom->Rndm()*bField.getXLength()*1.1+bField.getXStart()*.95);
//        pos.SetY(gRandom->Rndm()*bField.getYLength()*1.1+bField.getYStart()*.95);
//        pos.SetZ(gRandom->Rndm()*bField.getZLength()*1.1+bField.getZStart()*.95);
//
//
//        bField.linearInterp3D(pos,newFieldOut);
//    }
//    endTime=time(NULL);
//    cout <<"Linear Interp: " << difftime(endTime,startTime) <<" seconds run time." << endl; //~9 seconds
//
//    //Cubic - old
//    startTime = time(NULL);
//    for (int i = 0;i< numEvents;i++)
//    {
//        newFieldOut.SetXYZ(0,0,0);
//
//        pos.SetX(gRandom->Rndm()*bField.getXLength()*1.1+bField.getXStart()*.95);
//        pos.SetY(gRandom->Rndm()*bField.getYLength()*1.1+bField.getYStart()*.95);
//        pos.SetZ(gRandom->Rndm()*bField.getZLength()*1.1+bField.getZStart()*.95);
//
//
//        bField.cubicInterp3D(pos,newFieldOut);
//    }
//    endTime=time(NULL);
//    cout <<"Old Cubic Interp: " << difftime(endTime,startTime) <<" seconds run time." << endl; //~164 seconds
//
//    //Cubic - new
//    startTime = time(NULL);
//    for (int i = 0;i< numEvents;i++)
//    {
//        newFieldOut.SetXYZ(0,0,0);
//
//        pos.SetX(gRandom->Rndm()*bField.getXLength()*1.1+bField.getXStart()*.95);
//        pos.SetY(gRandom->Rndm()*bField.getYLength()*1.1+bField.getYStart()*.95);
//        pos.SetZ(gRandom->Rndm()*bField.getZLength()*1.1+bField.getZStart()*.95);
//
//
//        bField.cubicInterp3D2(pos,newFieldOut);
//    }
//    endTime=time(NULL);
//    cout <<"New Cubic Interp: " << difftime(endTime,startTime) <<" seconds run time." << endl;  //~74 seconds

///Main.cpp code for package for Peter 140205
//
//
//
//    //Load test field - stored on disk as txt file in cm
//    ROOTField testField;
//    testField.loadFieldFromFile(MRKRunDir+"Field_Test.txt", "", 1, 3,  3);
//
//
//    //Load RDK2 magnetic field - stored on disk as a TH3D histogram in a root file in meters
//    ROOTField bField;
//    bField.loadFieldFromFile(MRKRunDir+"Bfield_RDK2_Matt2.root", "solenoid", 1, 3,  3);
//    bField.setSymmetry(false,true,false);
//
//        int numEvents=5000000;
//
//    time_t startTime, endTime;
//    TVector3 newFieldOut,pos;
//
//    //Linear
//    cout << "Beginning Timing test" << endl;
//    startTime = time(NULL);
//    for (int i = 0;i< numEvents;i++)
//    {
//        newFieldOut.SetXYZ(0,0,0);
//
//        //Get random positions inside the field
//
//        pos.SetX(gRandom->Rndm()*bField.getXLength()*1.1+bField.getXStart()*.95);
//        pos.SetY(gRandom->Rndm()*bField.getYLength()*1.1+bField.getYStart()*.95);
//        pos.SetZ(gRandom->Rndm()*bField.getZLength()*1.1+bField.getZStart()*.95);
//
//
//        bField.linearInterp3D(pos,newFieldOut);
//    }
//    endTime=time(NULL);
//    cout <<"Linear Interp: " << difftime(endTime,startTime) <<" seconds run time." << endl;
//
//    //Cubic - old
//    startTime = time(NULL);
//    for (int i = 0;i< numEvents;i++)
//    {
//        newFieldOut.SetXYZ(0,0,0);
//
//        //Get random positions inside the field
//        pos.SetX(gRandom->Rndm()*bField.getXLength()*1.1+bField.getXStart()*.95);
//        pos.SetY(gRandom->Rndm()*bField.getYLength()*1.1+bField.getYStart()*.95);
//        pos.SetZ(gRandom->Rndm()*bField.getZLength()*1.1+bField.getZStart()*.95);
//
//
//        bField.cubicInterp3D(pos,newFieldOut);
//    }
//    endTime=time(NULL);
//    cout <<"Cubic Interp: " << difftime(endTime,startTime) <<" seconds run time." << endl;

///Recalculate length contraction in radial direction
//    RDK2Set threeBodySet(48,"130429","STD_3B","RID278","RID279","","STD NonRad W Fermi");
//    RDK2Set fourBodySet(48,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID283","STD Rad W Fermi");
//
//    fourBodySet.makeDerivedResultsFiles(GAMMA);
//
//    fourBodySet.calcEPGRateFor12IndividualBGOs(STD_EPCut);

///Analyze no PLO in LO extraction for Ba133 (Center3 = no PLO, Center2=PLO)
//    createDerivedResultsFilesForCal(10,"Ba133_Center3_140113");
//    createDerivedResultsFilesForCal(10,"Ba133_Center2_140124");
//    createDerivedResultsFilesForCal(10,"Cs137_Center3_140124");
//    createDerivedResultsFilesForCal(10,"Cs137_Center2_140124");
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Center3_140113");
//    createHistogramsForMCCalibraitonSourceData(10,"Ba133_Center2_140124");
//    createHistogramsForMCCalibraitonSourceData(10,"Cs137_Center3_140124");
//    createHistogramsForMCCalibraitonSourceData(10,"Cs137_Center2_140124");
//
//    analyzeLOofCalPeaks("C461", "Center3_140113", true,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C439", "Center3_140113", true,CALPEAKTIMINGCUT_OFFPEAK);
//
//    analyzeLOofCalPeaks("C461", "Center2_140124", true,CALPEAKTIMINGCUT_OFFPEAK);
//    analyzeLOofCalPeaks("C439", "Center2_140124", true,CALPEAKTIMINGCUT_OFFPEAK);

///Fit Cs137 peak
//    const int numDets=12;
//    TH1D* mcHists[12];
//
//    TString mcRunString="Center3_140124";
//    TString isotopeString="Cs137";
//    bool prebakeKhodyuk=true;
//    double isotopeRatio=0;
//    double energyCorrectionFromPLO=662./653.;
////    double energyCorrectionFromPLO=1;
//
//    TF1* theMCEnergyPeaks[5]={NULL,NULL,NULL,NULL,NULL};
//
//    const int numPeaks=getMCGaussPeaksToFit( "C441",  prebakeKhodyuk, theMCEnergyPeaks);
//
//    double mcMean[numDets][numPeaks];
//    double mcMeanError[numDets][numPeaks];
//    double mcSigma[numDets][numPeaks];
//    double mcAmplitude[numDets][numPeaks];
//    double averageMCAmplitude[numPeaks];
//    double averageMCMean[numPeaks];
//    double averageMCMeanError[numPeaks];
//    double averageMCSigma[numPeaks];
//
//
//
//    TRandom3 ranGen;
//
//    for (int j = 0;j< numPeaks;j++)
//    {
//        averageMCAmplitude[j]=0;
//        averageMCMean[j]=0;
//        averageMCMeanError[j]=0;
//        averageMCSigma[j]=0;
//
//          }
//
//    for (int i = 0;i< numDets;i++)
//    {
//        if(i!=6)
//        {
//
//
//            mcHists[i]=createCalMCHistWithRes(i+1, mcRunString, isotopeString,prebakeKhodyuk,isotopeRatio);
//
//            if(isotopeString=="Cs137Ba133")
//            {
//                convertTH1ToTXT(mcHists[i],TString(HISTS_DIR)+"BGO_Det"+int2str(i+1)+"_"+isotopeString+"_"+mcRunString+".txt");
//            }
//
//            mcHists[i]->Rebin(20);
//
//            TString titleString="Det"+int2str(i+1)+";Energy Deposited[";
//
//            if(prebakeKhodyuk)
//            {
//                titleString+="with";
//            }
//            else
//            {
//                titleString+="no";
//            }
//            titleString+= "LO response](keV);Counts";
//            mcHists[i]->SetTitle(titleString);
//            if(isotopeString=="Ba133")
//            {
//                mcHists[i]->GetXaxis()->SetRangeUser(0,500);
//            }
//            else
//            {
//                mcHists[i]->GetXaxis()->SetRangeUser(0,750);
//            }
//
//
//            for (int j = 0;j< numPeaks;j++)
//            {
//
//                fitHistWithGaus(mcHists[i],theMCEnergyPeaks[j]);
//
//                double deconvolutionRatio=1;
//                if(prebakeKhodyuk)
//                {
//                    if(isotopeString=="Cs137Ba133")
//                    {
//                        double csDeconvolutionRatio=deconvolutePrebakeKhodyuk(10, "Cs137_"+mcRunString, i+1, theMCEnergyPeaks[j]);
//                        csDeconvolutionRatio*=isotopeRatio;
//                        double baDeconvolutionRatio=deconvolutePrebakeKhodyuk(10, "Ba133_"+mcRunString, i+1, theMCEnergyPeaks[j]);
//                        deconvolutionRatio=(csDeconvolutionRatio+baDeconvolutionRatio)/(1.+isotopeRatio);
//
//                    }
//                    else
//                    {
//                        deconvolutionRatio=deconvolutePrebakeKhodyuk(10, isotopeString+"_"+mcRunString, i+1, theMCEnergyPeaks[j]);
//                    }
//                }
//
//                mcAmplitude[i][j] = theMCEnergyPeaks[j]->GetParameter(0);
//
//                mcMean[i][j] = theMCEnergyPeaks[j]->GetParameter(1)*energyCorrectionFromPLO*deconvolutionRatio;
//                mcMeanError[i][j] = theMCEnergyPeaks[j]->GetParError(1)*energyCorrectionFromPLO*deconvolutionRatio;
//                mcSigma[i][j] = theMCEnergyPeaks[j]->GetParameter(2)*energyCorrectionFromPLO*deconvolutionRatio;
//
//
//
//
//                //add code for prebake to include finding original mean
//
//                averageMCAmplitude[j]+=mcAmplitude[i][j]/11.;
//                averageMCMean[j]+=mcMean[i][j]/11.;
//                averageMCSigma[j]+=mcSigma[i][j]/11.;
//            }
//
//            TCanvas theCanvas("graphsCanvas","graphsCanvas",10,10,1280,800);
//            theCanvas.SetFillColor(kGray);
//            theCanvas.cd();
//            theCanvas.SetLeftMargin(.12);
//            theCanvas.SetRightMargin(.06);
//            gPad->SetTickx(1);
//            gPad->SetTicky(1);
//            gPad->SetFillColor(kWhite);
//
//            mcHists[i]->Draw("HIST");
//            for (int j = 0;j< numPeaks;j++)
//            {
//            	theMCEnergyPeaks[j]->Draw("L SAME");
//            }
//
//            TString imagePath=string(GRAPHS_DIR)+"CalMC_"+isotopeString+"_"+mcRunString;
//            if(prebakeKhodyuk)
//            {
//                imagePath+="_KhodyukPrebake";
//            }
//            imagePath+="_Det"+int2str(i+1)+"_CalData.png";
//
//            theCanvas.SaveAs(imagePath);
//
//            delete mcHists[i];
//        }
//    }
//
//    cout << "Average Mean:" << averageMCMean[1] << endl;

///Macros for E,P check with cubic interpolation (will change duration by hand.
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID284",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setUseCubicInterpolation 1\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setUseCubicInterpolation 1\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setUseCubicInterpolation 1");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID285",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setUseCubicInterpolation 1\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setUseCubicInterpolation 1\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setUseCubicInterpolation 1");

///Plot latest and greatest BGO data for Tom 140213
//    const int numMCHist=1;
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist140113_EDepGWLO_AverageIndiv_EPG_RID282.txt",BGO_DIM);
//
//    TH1* mcHists[numMCHist] = {mcHist};
//
//    mcHists[0]->SetTitle("Theory+MC;Energy Detected (keV);#frac{EPG}{EP} Ratio");
//
//
//    for (int i = 0;i< numMCHist;i++)
//    {
//    	mcHists[i]->Rebin(5);
//    }
//
//
//
//    TH1* expHist[1];
//
//    TString histPath=HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt";
//
//
//    expHist[0]=getTabSeperatedHist(histPath,BGO_DIM);
//
//
//    expHist[0]->Rebin(5);
//    expHist[0]->SetTitle("Exp.;Energy Detected (keV);#frac{EPG}{EP}");
//    expHist[0]->GetYaxis()->SetRangeUser(0,14e-6);
//    mcHists[0]->GetYaxis()->SetRangeUser(0,14e-6);
//
//
//    TString titleString;
//    plotExpVersusMCToImage(1,expHist,numMCHist,mcHists,"Average of BGO Crystal Detectors","NSACMeetingBGOPlot140213.png");

///The Case of the Mysterious Fermi Function - by M.J. Bales 140213
////    int numIterationsPerBin=1000;
////    TH1D* noFermiHist=calcBRHist(numIterationsPerBin,690,10.,700.,false);
////    noFermiHist->SetName("noFermiHist");
////    noFermiHist->SetTitle("No Fermi");
////    TH1D* fermiHist=calcBRHist(numIterationsPerBin,690,10,700,true);
////    fermiHist->SetName("fermiHist");
////    fermiHist->SetTitle("With Fermi");
//
//    TH1D* fermiHist=getTabSeperatedHist(HistsDir+"fermiHist140213.txt",690,10,700);
//    TH1D* noFermiHist=getTabSeperatedHist(HistsDir+"noFermiHist140213.txt",690,10,700);
//
////    convertTH1ToTXT(noFermiHist,HistsDir+"noFermiHist140213.txt");
////    convertTH1ToTXT(fermiHist,HistsDir+"fermiHist140213.txt");
//
//    TH1D* ratioHist=(TH1D*) fermiHist->Clone("ratioHist");
//    ratioHist->Divide(noFermiHist);
//    ratioHist->SetTitle("Ratio of Fermi to No Fermi BR");
////    convertTH1ToTXT(ratioHist,HistsDir+"fermiRatioHist140213.txt");
//
//
//    TH1* mcHists[2]={noFermiHist,fermiHist};
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"BR for 1 keV bins, With and Without Fermi Function","fermiComparison140213.png");
//    mcHists[0]=ratioHist;
//    mcHists[0]->GetYaxis()->SetRangeUser(1.028,1.058);
//    plotExpVersusMCToImage(0,NULL,1,mcHists,"Ratio of Fermi to No Fermi BR","fermiBRRatio140213.png");

//    cEvents theEvents;
//    theEvents.makeEventFiles(EventsDir+"
//    double  nFermi=theEvents.calcBR(390,410,1000000,false);
//    double  wFermi=theEvents.calcBR(390,410,1000000,true);
//    double ratio=wFermi/nFermi;
//    cout << "BR no Fermi: " << nFermi << endl;
//    cout << "BR with Fermi: " << wFermi  << endl;
//    cout << "BR ratio: " <<ratio << endl;

//    int numIterationsPerBin=1000;
//    TH1D* noFermiHist=calcBRHist(numIterationsPerBin,690,10.,700.,false);
//    noFermiHist->SetName("noFermiHist");
//    noFermiHist->SetTitle("No Fermi");
//    TH1D* fermiHist=calcBRHist(numIterationsPerBin,690,10,700,true);
//    fermiHist->SetName("fermiHist");
//    fermiHist->SetTitle("With Fermi");
//
////    TH1D* fermiHist=getTabSeperatedHist(HistsDir+"fermiHist140213.txt",690,10,700);
////    TH1D* noFermiHist=getTabSeperatedHist(HistsDir+"noFermiHist140213.txt",690,10,700);
//
//    convertTH1ToTXT(noFermiHist,HistsDir+"noFermiHist140214.txt");
//    convertTH1ToTXT(fermiHist,HistsDir+"fermiHist140214.txt");
//
//    TH1D* ratioHist=(TH1D*) fermiHist->Clone("ratioHist");
//    ratioHist->Divide(noFermiHist);
//    ratioHist->SetTitle("Ratio of Fermi to No Fermi BR");
//    convertTH1ToTXT(ratioHist,HistsDir+"fermiRatioHist140214.txt");
//
//
//    TH1* mcHists[2]={noFermiHist,fermiHist};
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"BR for 1 keV bins, With and Without Fermi Function","fermiComparison140214.png");
//    mcHists[0]=ratioHist;
//    mcHists[0]->GetYaxis()->SetRangeUser(1.028,1.058);
//    plotExpVersusMCToImage(0,NULL,1,mcHists,"Ratio of Fermi to No Fermi BR","fermiBRRatio140214.png");

	///make norm constant fit to speed things up
//    cEvents theEvents;
//    double temp=theEvents.getFitNormConstant(FOURBODY,100);
//    TH1D* normConstHist=new TH1D("normConstHist","Hist of Norm Constant for Event Gen",10,300,400);
//    for (int i = 0;i< normConstHist->GetNbinsX();i++)
//    {
//    	normConstHist->SetBinContent(i+1,theEvents.calcNormConstant(FOURBODY,0,normConstHist->GetBinLowEdge(i+1)));
//    }
//
//   convertTH1ToTXT(normConstHist,HistsDir+"normConstHist140214c.txt");
////
////    TH1D* normConstHist=getTabSeperatedHist(HistsDir+"normConstHist140214c.txt",10,300,400);
////
//    for (int i = 0;i< normConstHist->GetNbinsX();i++)
//    {
//            normConstHist->SetBinError(i+1,normConstHist->GetBinContent(i+1)*0.01);
//    }
////
//    TF1 fitFunc("fitFunc","pol2",normConstHist->GetBinLowEdge(2),getHistExtentHigh(normConstHist));
//    fitFunc.SetParameter(0,-5.51e-14);
//    fitFunc.SetParameter(1,1.11e-12);
//    fitFunc.SetParameter(2,5.11e-15);
////    fitFunc.SetParameter(3,-3.39e-10);
//    TH1* mcHists[1]={normConstHist};
//    normConstHist->Fit("fitFunc","VMRN");
//    TCanvas* theCanvas=plotExpVersusMC(0, NULL,1,mcHists,"Norm Const");
//    fitFunc.Draw("same");
//    theCanvas->SetLogy();
//    TString imagePath=GraphsDir+"normConst140214d.png";
//    theCanvas->SaveAs(imagePath);

///Try Rob's code:
//    FermiRCooper fermiClass;
//    TH1D* robHist= fermiClass.Fermi(1000, 1234567, .050+fermiClass.Me, fermiClass.eEmax);
//    convertTH1ToTXT(robHist,HistsDir+"fermiRatioHistRob140216.txt");
//
////    TH1D* robHist=getTabSeperatedHist(HistsDir+"fermiRatioHistRob140216.txt",690,.010,.700);
//
//    int numIterationsPerBin=10000;
//    TH1D* noFermiHist=calcBRHist(numIterationsPerBin,690,10.,700.,50,false);
//    noFermiHist->SetName("noFermiHist");
//    noFermiHist->SetTitle("No Fermi");
//    TH1D* fermiHist=calcBRHist(numIterationsPerBin,690,10,700,50,true);
//    fermiHist->SetName("fermiHist");
//    fermiHist->SetTitle("With Fermi");
//
//    convertTH1ToTXT(noFermiHist,HistsDir+"noFermiHist140216.txt");
//    convertTH1ToTXT(fermiHist,HistsDir+"fermiHist140216.txt");
//
////    TH1D* fermiHist=getTabSeperatedHist(HistsDir+"fermiHist140213.txt",690,10,700);
////    TH1D* noFermiHist=getTabSeperatedHist(HistsDir+"noFermiHist140213.txt",690,10,700);
//
//    TH1D* ratioHist=(TH1D*) fermiHist->Clone("ratioHist");
//    ratioHist->Divide(noFermiHist);
//    ratioHist->SetTitle("Ratio of Fermi to No Fermi BR");
//    convertTH1ToTXT(ratioHist,HistsDir+"fermiRatioHist140216.txt");
//
//    TH1D* mattHist=getTabSeperatedHist(HistsDir+"fermiRatioHist140216.txt",690,.010,.700);
//
//    robHist->SetTitle("Rob");
//    mattHist->SetTitle("Matt;MeV;BR ratio");
//
//    TH1* mcHists[2]={mattHist,robHist};
//    mcHists[0]->GetYaxis()->SetRangeUser(1.028,1.058);
//    plotExpVersusMCToImage(0,NULL,2,mcHists,"Ratio of Fermi to No Fermi BR: Rob EE>50 keV","fermiBRRatioRob140214.png");

///Ok all figured out now I just need to see about the event method.
//    cEvents theEvents;
//    double gLow=10;
//    int numEventFiles=48;
//    int numEventsPerFile=1000000;
//    int numEvents=numEventFiles*numEventsPerFile;
//
////    theEvents.makeEventFiles(EventsDir+"Events_STD_4B_100eV_Endpoint_140217_Set0.root",MRKRunDir+"Fluxmap_RDK2_Matt1.txt",numEventFiles,numEventsPerFile,FOURBODY,0,false,gLow,EKEMAX,true);
////    theEvents.makeEventFiles(EventsDir+"Events_STD_4B_100eV_Endpoint_140217NoFermi_Set0.root",MRKRunDir+"Fluxmap_RDK2_Matt1.txt",numEventFiles,numEventsPerFile,FOURBODY,0,false,gLow,EKEMAX,false);
//    double fermiBR = theEvents.calcBR(gLow,EKEMAX,1000000,true);
//    double noFermiBR = theEvents.calcBR(gLow,EKEMAX,1000000,false);
//    cout << "Fermi on BR:" << fermiBR << endl; //0.0107331
//    cout << "Fermi off BR:" << noFermiBR<< endl; //0.0103071
//
//    //Old attempt
//    RDK2Set* oldFermiSet=new RDK2Set(48,"130429","STD_4B_100eV_Endpoint","","","","STD Rad W Fermi");
//    RDK2Set* oldNoFermiSet=new RDK2Set(48,"130116","STD_4B_100eV_Endpoint","","","","STD Rad");
//
//    RDK2Set* fermiSet=new RDK2Set(numEventFiles,"140218","STD_4B_100eV_Endpoint","","","","STD Rad W Fermi");
//    RDK2Set* noFermiSet=new RDK2Set(1,"14021NoFermi","STD_4B_100eV_Endpoint","","","","STD Rad");
//
//    HistDim brDim={690,10,700};
//
//    TH1* eg0
//
//    TH1* eg0Hists[4]={getTabSeperatedHist(HistsDir+"noFermiHist140216.txt",brDim),getTabSeperatedHist(HistsDir+"fermiHist140216.txt",brDim),oldNoFermiSet->makeHist("oldEG0Hist","t.eg0","t.ee0 > 50",brDim),oldFermiSet->makeHist("newEG0Hist","t.eg0","t.ee0 > 50",brDim)};
//    //TH1* eg0Hists[4]={getTabSeperatedHist(HistsDir+"noFermiHist140216.txt",brDim),getTabSeperatedHist(HistsDir+"fermiHist140216.txt",brDim),getTabSeperatedHist(HistsDir+"noFermiEventMethodHist140217b.txt",brDim),getTabSeperatedHist(HistsDir+"fermiEventMethodHist140217b.txt",brDim)};
//
//    convertTH1ToTXT(eg0Hists[2],HistsDir+"noFermiEventMethodHist140217b.txt");
//    convertTH1ToTXT(eg0Hists[3],HistsDir+"fermiEventMethodHist140217b.txt");
//
//    eg0Hists[2]->SetTitle("Events_NoFermi;Gamma Energy (keV);BR");
//    eg0Hists[3]->SetTitle("Events_Fermi;Gamma Energy (keV);BR");
//    eg0Hists[0]->SetTitle("Direct_NoFermi;Gamma Energy (keV);BR");
//    eg0Hists[1]->SetTitle("Direct_Fermi;Gamma Energy (keV);BR");
//
//    eg0Hists[2]->Scale(noFermiBR/numEvents);
//    eg0Hists[3]->Scale(fermiBR/numEvents);
//    //eg0Hists[0]->GetXaxis()->SetRangeUser(0,0.004);
//
//    for (int i = 0;i< 4;i++)
//    {
//        eg0Hists[i]->Rebin(30);
//    }
//
//    plotExpVersusMCToImage(0,NULL,4,eg0Hists,"Comparison of BR methods","BRMethods140219.png",false);
//    plotExpVersusMCToImage(0,NULL,4,eg0Hists,"Comparison of BR methods","BRMethods140219_Log.png",true);
//    plotResidualsToImage(3,&(eg0Hists[1]),1,&(eg0Hists[0]),"#frac{X-Direct_NoFermi}{Direct_NoFermi}","BRMethods140219_Resid.png",-.1,.1);
//
////    plotResidualsToImage(1,&(eg0Hists[1]),1,eg0Hists,"#frac{Fermi-NoFermi}{NoFermi}","fermiRatioEventMethod.png",-.1,.1);
////
////    makeMacrosForEventsOnCluster(0,48,"MAKEEVENTS_4B",FOURBODY); //New files to compare to old ones

///Make other NSAC Talk Plots
	///PID8
////    createPlotOfRIDSequence(48,278);
//   // makeAllStandardExpHists(getCurrentDateString(),"PID8","G7_Cycle5_8","EPBGO",NUM_EXP_SERIES_G7_Cycle5_8,EXP_SERIES_ID_LIST_G7_Cycle5_8);
//    //makeAllStandardExpHists(getCurrentDateString(),"PID8","G7","EPBGO",NUM_EXP_SERIES_G7,EXP_SERIES_ID_LIST_G7);
// //   makeAllStandardExpHists(getCurrentDateString(),"PID8","G7_Cycle5_8","EP",NUM_EXP_SERIES_G7_Cycle5_8,EXP_SERIES_ID_LIST_G7_Cycle5_8);
//
//    TH1* expHists[4];
//    TH1* mcHists[1];
//
//    ///Plot EP TOFP
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140221_PID8_G7_Cycle5_8_ToFP_EP.txt",TOFP_DIM2);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist140220_ToFP_EP_RID282.txt",TOFP_DIM2);
//
//    expHists[0]->SetTitle("Exp;E-P Time Delta (s);Experimental Counts");
//    mcHists[0]->SetTitle("Theory+MC;E-P Time Delta (s);Experimental Counts");
//
//    scaleHistogramsTogether(expHists[0],mcHists[0]);
//
//    expHists[0]->Rebin(3);
//    mcHists[0]->Rebin(3);
//
//    expHists[0]->GetXaxis()->SetRangeUser(0,10e-6);
//    mcHists[0]->GetXaxis()->SetRangeUser(0,10e-6);
//
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Proton Time of Flight","ExpVsMC_140221_PID8_G7_Cycle5_8_RID278_EP_ToFP.png");
//    delete expHists[0];
//    delete mcHists[0];
//
//   ///Plot EP Electron Energy
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140221_PID8_G7_Cycle5_8_EDepE_EP.txt",EDEPE_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist140220_EDepE_EP_RID282.txt",EDEPE_DIM);
//
//    expHists[0]->SetTitle("Exp;Electron Energy Detected (keV);Experimental Counts");
//    mcHists[0]->SetTitle("Theory+MC;Electron Energy Detected (keV);Experimental Counts");
//
//    scaleHistogramsTogether(expHists[0],mcHists[0]);
//
////    expHists->GetXaxis()->SetRangeUser(0,10e-6);
//
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Electron Energy Detected","ExpVsMC_140221_PID8_G7_Cycle5_8_RID278_EP_EDepE.png");
//
//    expHists[0]->Rebin(5);
//    mcHists[0]->Rebin(5);
//
//    plotResidualsToImage(1,expHists,1,mcHists,"Electron Energy Detected","ExpVsMC_140221_PID8_G7_Cycle5_8_RID278_EP_EDepE_Residuals.png",-.2,.2);
//    delete expHists[0];
//    delete mcHists[0];
//
//    ///Plot EPG BGO Average
//    int bgoRebin=5;
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140221_PID8_G7_Cycle5_8_EDepG_BGOAverage_EPG.txt",BGO_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist140220_EDepGWLO_AverageIndiv_EPG_RID282.txt",BGO_DIM);
//
//    expHists[0]->SetTitle("Exp;BGO Pulse Height (511 = 511 keV);#frac{EPG}{EP}");
//    mcHists[0]->SetTitle("Theory+MC;BGO Pulse Height (511 = 511 keV);#frac{EPG}{EP}");
//
//    expHists[0]->Rebin(bgoRebin);
//    mcHists[0]->Rebin(bgoRebin);
//    expHists[0]->GetYaxis()->SetRangeUser(0,14e-6);
//    mcHists[0]->GetYaxis()->SetRangeUser(0,14e-6);
//    mcHists[0]->GetXaxis()->SetRangeUser(0,800);
//
////    expHists[1]=getTabSeperatedHist(HistsDir+"ExpMatt_Hists130724_G7_BGOEDepAvg_EPG.txt",BGO_DIM);
////    expHists[1]->SetTitle("Old;Energy Detected (keV);#frac{EPG}{EP}");
////    expHists[1]->Rebin(bgoRebin);
////
////    expHists[2]=getTabSeperatedHist(HistsDir+"ExpJeff_Hists130712_2_10_BGOEDepAvg_EPG.txt",BGO_DIM);
//////    expHists[2]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140221_PID7_G7_Cycle5_8_EDepG_BGOAverage_EPG.txt",BGO_DIM);
////    expHists[2]->SetTitle("Jeff;Energy Detected (keV);#frac{EPG}{EP}");
////    expHists[2]->Rebin(bgoRebin);
////
////    expHists[3]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140221_PID8_G7_EDepG_BGOAverage_EPG.txt",BGO_DIM);
////    expHists[3]->SetTitle("PID_ALL;Energy Detected (keV);#frac{EPG}{EP}");
////    expHists[3]->Rebin(bgoRebin);
//
//
//
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Average of 12 BGO Photon Detectors","ExpVsMC_140221_PID8_G7_Cycle5_8_RID278_EPG_EDepG_BGOAverage.png");
//    expHists[0]->Rebin(4);
//    mcHists[0]->Rebin(4);
//    plotResidualsToImage(1,expHists,1,mcHists,"Average of 12 BGO Photon Detectors","ExpVsMC_140221_PID8_G7_Cycle5_8_RID278_EPG_EDepG_BGOAverage_Residuals.png");
//    delete expHists[0];
////    delete expHists[1];
////    delete expHists[2];
//    delete mcHists[0];
//
//    ///Plot EPG BAPD Average
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPBen_Hist131009_BAPD_AverageIndiv_EPG.txt",BAPD_DIM);
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist140220_BAPD_AverageIndiv_EPG_RID282.txt",BAPD_DIM);
//
//    expHists[0]->SetTitle("Exp;BAPD Pulse Height (5.9 = 5.9 keV);#frac{EPG}{EP}");
//    mcHists[0]->SetTitle("Theory+MC;BAPD Pulse Height (5.9 = 5.9 keV);#frac{EPG}{EP}");
//
//    for (int i = 0;i< 2;i++)
//    {
//    	mcHists[0]->SetBinContent(i+1,0);
//    	mcHists[0]->SetBinError(i+1,0);
//    	expHists[0]->SetBinContent(i+1,0);
//    	expHists[0]->SetBinError(i+1,0);
//
//    }
//
//    expHists[0]->Rebin(2);
//    mcHists[0]->Rebin(2);
//
//
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Average of 3 BAPD Photon Detectors","ExpVsMC_140221_PID8_G7_Cycle5_8_RID278_EPG_EDepG_BAPDAverage.png");
//
//    delete expHists[0];
//    delete mcHists[0];

///Check BR for event vs direct method 140221
//    cEvents theEvents;
//    double gLow=0.1;
//    int numEventFiles=48;
//    int numEventsPerFile=1000000;
//    int numEvents=numEventFiles*numEventsPerFile;
//    TH1* eg0Hists[3];
//    HistDim brDim={690,10,700};
//
////    double fermiBR = theEvents.calcBR(gLow,EKEMAX,1000000,true);
////    cout << "Fermi on (no ecut) BR:" << fermiBR << endl; //0.0102108
////    RDK2Set* fermiSet=new RDK2Set(numEventFiles,"140218","STD_4B_100eV_Endpoint","","","","STD Rad W Fermi");
////    eg0Hists[1]=fermiSet->makeHist("newEG0Hist","t.eg0","t.ee0 > 50",brDim);
////    eg0Hists[1]->Scale(fermiBR/numEvents);
////    convertTH1ToTXT(eg0Hists[1],HistsDir+"fermiEventMethodHist140221.txt");
//
//    eg0Hists[1]=getTabSeperatedHist(HistsDir+"fermiEventMethodHist140221.txt",brDim);
//
//    int numIterationsPerBin=100000;
//    eg0Hists[2]=calcBRHist(numIterationsPerBin,690,10.,700.,50.,false);
//    eg0Hists[2]->SetName("noFermiHist");
//    eg0Hists[2]->SetTitle("No Fermi");
//    eg0Hists[0]=calcBRHist(numIterationsPerBin,690,10,700,50,true);
//    eg0Hists[0]->SetName("fermiHist");
//    eg0Hists[0]->SetTitle("With Fermi");
//
//    convertTH1ToTXT(eg0Hists[2],HistsDir+"noFermiHist140222.txt");
//    convertTH1ToTXT(eg0Hists[0],HistsDir+"fermiHist140222.txt");
//
////    eg0Hists[0]=getTabSeperatedHist(HistsDir+"fermiHist140222.txt",brDim);
////    eg0Hists[2]=getTabSeperatedHist(HistsDir+"noFermiHist140222.txt",brDim);
//
//    for (int i = 1;i<=brDim.numBins ;i++)
//    {
//    	eg0Hists[0]->SetBinError(i,eg0Hists[0]->GetBinContent(i)*0.001);
//    	eg0Hists[2]->SetBinError(i,eg0Hists[2]->GetBinContent(i)*0.001);
//    }
//
//
//
//    eg0Hists[1]->SetTitle("Events_Fermi;Gamma Energy (keV);BR");
//    eg0Hists[2]->SetTitle("Direct_NoFermi;Gamma Energy (keV);BR");
//    eg0Hists[0]->SetTitle("Direct_Fermi;Gamma Energy (keV);BR");
//
//
//    for (int i = 0;i< 3;i++)
//    {
//        eg0Hists[i]->Rebin(5);
//    }
//
//    plotExpVersusMCToImage(0,NULL,3,eg0Hists,"Comparison of BR methods","BRMethods140221.png",false);
//    plotExpVersusMCToImage(0,NULL,3,eg0Hists,"Comparison of BR methods","BRMethods140221_Log.png",true);
//    plotResidualsToImage(2,eg0Hists,1,&(eg0Hists[2]),"#frac{X-Direct_NoFermi}{Direct_NoFermi}","BRMethods140221_Resid.png",-.1,.1);

///Make better graph for Tom 130204
//    RDK2Set hmgSet(44,"100421","hmg","RID114","RID115","","HMG NonRad");
//    RDK2Set lowSet(44,"120208","radlow","RID105","RID106","RID107","RadLow Set");
//    RDK2Set midSet(44,"120208","radmid","RID108","RID109","RID110","RadMid Set");
//    RDK2Set highSet(44,"120208","radhigh","RID111","RID112","RID113","RadHigh Set");
//
////    lowSet.makeDerivedResultsFiles(GAMMA);
////    midSet.makeDerivedResultsFiles(GAMMA);
////    highSet.makeDerivedResultsFiles(GAMMA);
//
//
//    //Make other cuts
//    TH1* mcHists[2];
//    TH1* mcHists2[1];
//    TH1* expHists[1];
////    mcHists[0]=combineAndScaleRadSetsToHist(&hmgSet,&lowSet, &midSet, &highSet,"((gDChn.BGO1EDepResOnly > 10.)*gDChn.BGO1EDepResOnly + (gDChn.BGO2EDepResOnly > 10.)*gDChn.BGO2EDepResOnly + (gDChn.BGO3EDepResOnly > 10.)*gDChn.BGO3EDepResOnly + (gDChn.BGO4EDepResOnly > 10.)*gDChn.BGO4EDepResOnly + (gDChn.BGO5EDepResOnly > 10.)*gDChn.BGO5EDepResOnly + (gDChn.BGO6EDepResOnly > 10.)*gDChn.BGO6EDepResOnly + (gDChn.BGO7EDepResOnly > 10.)*gDChn.BGO7EDepResOnly + (gDChn.BGO8EDepResOnly > 10.)*gDChn.BGO8EDepResOnly + (gDChn.BGO9EDepResOnly > 10.)*gDChn.BGO9EDepResOnly + (gDChn.BGO10EDepResOnly > 10.)*gDChn.BGO10EDepResOnly + (gDChn.BGO11EDepResOnly > 10.)*gDChn.BGO11EDepResOnly + (gDChn.BGO12EDepResOnly > 10.)*gDChn.BGO12EDepResOnly )",STD_EPCut,"((gDChn.BGO1EDepResOnly > 10.)*1 + (gDChn.BGO2EDepResOnly > 10.)*1 + (gDChn.BGO3EDepResOnly > 10.)*1 + (gDChn.BGO4EDepResOnly > 10.)*1 + (gDChn.BGO5EDepResOnly > 10.)*1 + (gDChn.BGO6EDepResOnly > 10.)*1 + (gDChn.BGO7EDepResOnly > 10.)*1 + (gDChn.BGO8EDepResOnly > 10.)*1 + (gDChn.BGO9EDepResOnly > 10.)*1 + (gDChn.BGO10EDepResOnly > 10.)*1 + (gDChn.BGO11EDepResOnly > 10.)*1 + (gDChn.BGO12EDepResOnly > 10.)*1 ) > 0",BGO_DIM);
////    convertTH1ToTXT(mcHists[0],HistsDir+"MC_Hist130116_EDepGRawWRes_BGO_EPG_Sum_RID113.txt");
//
//
//    //Jeff
//    HistDim JeffBGODim={160,0,800.};
//    expHists[0]=getTabSeperatedHist(HistsDir+"Jeff_BGO_EPG_SpectrumInCounts121108.txt",JeffBGODim);
//    double jeffScaleFactor=1./5.12044e7;
//    expHists[0]->Scale(jeffScaleFactor);
//    expHists[0]->SetTitle("Exp.;BGO Pulse Height (511 = 511 keV);#frac{EPG}{EP} (count rate ratio)");
//
//// Need to fix the filenames
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist130116_EDepGRawWRes_BGO_EPG_Sum_RID113.txt",BGO_DIM);
//     mcHists[1]=getTabSeperatedHist(HistsDir+"MC_Hist130116_EDepG_BGO_EPG_Sum_RID113.txt",BGO_DIM);
//
//
//    mcHists[0]->SetTitle("Proportional MC;BGO Pulse Height (511 = 511 keV);#frac{EPG}{EP} (count rate ratio)");
//    mcHists[1]->SetTitle("Nonproportional MC;BGO Pulse Height (511 = 511 keV);#frac{EPG}{EP} (count rate ratio)");
//    mcHists2[0]=mcHists[1];
//
//    mcHists[0]->Rebin(5);
//    mcHists[1]->Rebin(5);
//    expHists[0]->Rebin(1);
//
//    mcHists[0]->GetXaxis()->SetRangeUser(0,300);
//    mcHists[1]->GetXaxis()->SetRangeUser(0,300);
//    expHists[0]->GetXaxis()->SetRangeUser(0,300);
//    plotExpVersusMCToImage(1,expHists,2,mcHists,"With and Without Nonproportionality","EDepG_NonpropComparison_MC_140221.png");
//    mcHists[0]->Rebin(2);
//    mcHists[1]->Rebin(2);
//    expHists[0]->Rebin(2);
//    mcHists[0]->GetXaxis()->SetRangeUser(0,300);
//    mcHists[1]->GetXaxis()->SetRangeUser(0,300);
//    expHists[0]->GetXaxis()->SetRangeUser(0,300);
//    plotResidualsToImage(1, mcHists2,1,mcHists,"With and Without Nonproportionality #frac{Prop.-Nonprop}{Nonprop}","EDepG_NonpropComparisonResiduals_MC_140221.png",-.4,.4);
//    mcHists[0]->GetXaxis()->SetRangeUser(0,300);
//    mcHists[1]->GetXaxis()->SetRangeUser(0,300);
//    expHists[0]->GetXaxis()->SetRangeUser(0,300);
//    plotResidualsToImage(1,expHists,2,mcHists,"With and Without Nonproportionality #frac{Exp-MC}{MC}","EDepG_NonpropComparisonResidualsWExp_MC_140221.png",-.4,1.);

///Make noFermiEvents on Cluster
	//makeMacrosForEventsOnCluster(0,48,"MAKEEVENTS_4B_NoFermi",FOURBODY,"140218_NoFermi"); //New files to compare to old ones
///Make
//      cEvents theEvents;
//    double gLow=10;
//    int numEventFiles=48;
//    int numEventsPerFile=1000000;
//    int numEvents=numEventFiles*numEventsPerFile;
//    TString filePath;
//
//    theEvents.setGoodNormConstant(FOURBODY,0,gLow);
//    for (int i = 0;i< numEventFiles;i++)
//    {
//        filePath=EventsDir+"Events_STD_4B_1keV_Endpoint_140222_Set"+int2str(i)+".root";
//    	theEvents.makeEventFile(filePath,MRKRunDir+"Fluxmap_RDK2_Matt1.txt",numEventsPerFile,FOURBODY,0,false,gLow,EKEMAX,true);
//    	filePath=EventsDir+"Events_STD_4B_1keV_Endpoint_140222NoFermi_Set"+int2str(i)+".root";
//    	theEvents.makeEventFile(filePath,MRKRunDir+"Fluxmap_RDK2_Matt1.txt",numEventsPerFile,FOURBODY,0,false,gLow,EKEMAX,false);
//    }

///Check BR for event vs direct method 140221
//    cEvents theEvents;
//    int numEventFiles=48;
//    int numEventsPerFile=1000000;
//    int numEvents=numEventFiles*numEventsPerFile;
//    TH1* fermiHists[3];
//    TH1* noFermiHists[3];
//    TH1* directHists[2];
//    TH1* eventHists[4];
//    TH1* highEventHists[2];
//    TH1* lowEventHists[2];
//    TH1* allHists[6];
//    double branchingRatio;
//    RDK2Set* eventSet;
//
//    HistDim brDim={690,10,700};
//
//
//
////    ///Make low fermi hist
//////    branchingRatio = theEvents.calcBR(0.1,EKEMAX,0,1000000,true);
////    branchingRatio = 0.0102767;
////    eventSet=new RDK2Set(numEventFiles,"140218","STD_4B_100eV_Endpoint","","","","STD Rad W Fermi");
////    lowEventHists[0]=eventSet->makeHist("a","t.eg0","t.ee0 > 50",brDim);
////    lowEventHists[0]->Scale(branchingRatio/numEvents);
////    convertTH1ToTXT(lowEventHists[0],HistsDir+"fermiEventMethodHist140222.txt");
////    delete eventSet;
////
////    ///Make low no fermi hist
//////    branchingRatio = theEvents.calcBR(0.1,EKEMAX,0,1000000,false);
////    branchingRatio = 0.0099645;
////    eventSet=new RDK2Set(numEventFiles,"140218_NoFermi","STD_4B_100eV_Endpoint","","","","STD Rad W Fermi");
////    lowEventHists[1]=eventSet->makeHist("b","t.eg0","t.ee0 > 50",brDim);
////    lowEventHists[1]->Scale(branchingRatio/numEvents);
////    convertTH1ToTXT(lowEventHists[1],HistsDir+"noFermiEventMethodHist140222.txt");
////    delete eventSet;
////
////    ///Make high fermi hist
//////    branchingRatio = theEvents.calcBR(10,EKEMAX,0,1000000,true);
////    branchingRatio = 0.0035493;
////    eventSet=new RDK2Set(numEventFiles,"140222","STD_4B_1keV_Endpoint","","","","STD Rad W Fermi");
////    highEventHists[0]=eventSet->makeHist("c","t.eg0","t.ee0 > 50",brDim);
////    highEventHists[0]->Scale(branchingRatio/numEvents);
////    convertTH1ToTXT(highEventHists[0],HistsDir+"fermiEventHighMethodHist140222.txt");
////    delete eventSet;
////
////    ///Make high no fermi hist
//////    branchingRatio = theEvents.calcBR(10,EKEMAX,0,1000000,false);
////    branchingRatio = 0.003445;
////    eventSet=new RDK2Set(numEventFiles,"140222NoFermi","STD_4B_1keV_Endpoint","","","","STD Rad W Fermi");
////    highEventHists[1]=eventSet->makeHist("d","t.eg0","t.ee0 > 50",brDim);
////    highEventHists[1]->Scale(branchingRatio/numEvents);
////    convertTH1ToTXT(highEventHists[1],HistsDir+"noFermiHighEventMethodHist140222.txt");
////    delete eventSet;
//
//    ///Make BR Hists
////    int numIterationsPerBin=100000;
////    eg0Hists[2]=calcBRHist(numIterationsPerBin,690,10.,700.,50.,false);
////    eg0Hists[2]->SetName("noFermiHist");
////    eg0Hists[2]->SetTitle("No Fermi");
////    eg0Hists[0]=calcBRHist(numIterationsPerBin,690,10,700,50,true);
////    eg0Hists[0]->SetName("fermiHist");
////    eg0Hists[0]->SetTitle("With Fermi");
////
////    convertTH1ToTXT(eg0Hists[2],HistsDir+"noFermiHist140222.txt");
////    convertTH1ToTXT(eg0Hists[0],HistsDir+"fermiHist140222.txt");
//
//    ///Load Hists from Text
//    directHists[0]=getTabSeperatedHist(HistsDir+"fermiHist140222.txt",brDim);
//    directHists[1]=getTabSeperatedHist(HistsDir+"noFermiHist140222.txt",brDim);
//    lowEventHists[0]=getTabSeperatedHist(HistsDir+"fermiEventMethodHist140222.txt",brDim);
//    lowEventHists[1]=getTabSeperatedHist(HistsDir+"noFermiEventMethodHist140222.txt",brDim);
//    highEventHists[0]=getTabSeperatedHist(HistsDir+"fermiEventHighMethodHist140222.txt",brDim);
//    highEventHists[1]=getTabSeperatedHist(HistsDir+"noFermiHighEventMethodHist140222.txt",brDim);
//
//    directHists[0]->SetTitle("Direct_Fermi;Gamma Energy (keV);BR");
//    directHists[1]->SetTitle("Direct_NoFermi;Gamma Energy (keV);BR");
//    lowEventHists[0]->SetTitle("EventLow_Fermi;Gamma Energy (keV);BR");
//    lowEventHists[1]->SetTitle("EventLow_NoFermi;Gamma Energy (keV);BR");
//    highEventHists[0]->SetTitle("EventHigh_Fermi;Gamma Energy (keV);BR");
//    highEventHists[1]->SetTitle("EventHigh_NoFermi;Gamma Energy (keV);BR");
//
//    eventHists[0]=lowEventHists[0];
//    eventHists[1]=lowEventHists[1];
//    eventHists[2]=highEventHists[0];
//    eventHists[3]=highEventHists[1];
//
//    fermiHists[0]=directHists[0];
//    fermiHists[1]=lowEventHists[0];
//    fermiHists[2]=highEventHists[0];
//
//    noFermiHists[0]=directHists[1];
//    noFermiHists[1]=lowEventHists[1];
//    noFermiHists[2]=highEventHists[1];
//
//    allHists[0]=fermiHists[0];
//    allHists[1]=fermiHists[1];
//    allHists[2]=fermiHists[2];
//    allHists[3]=noFermiHists[0];
//    allHists[4]=noFermiHists[1];
//    allHists[5]=noFermiHists[2];
//
//    for (int i = 1;i<=brDim.numBins ;i++)
//    {
//    	directHists[0]->SetBinError(i,directHists[0]->GetBinContent(i)*0.001);
//    	directHists[1]->SetBinError(i,directHists[1]->GetBinContent(i)*0.001);
//    }
//
//
//    for (int i = 0;i< 6;i++)
//    {
//        allHists[i]->Rebin(30);
//    }
//
//    plotExpVersusMCToImage(0,NULL,6,allHists,"Comparison of BR methods","BRMethods140221.png",false);
//    plotExpVersusMCToImage(0,NULL,6,allHists,"Comparison of BR methods","BRMethods140221_Log.png",true);
//    plotResidualsToImage(3,fermiHists,3,noFermiHists,"#frac{Fermi-NoFermi}{NoFermi}","BRMethods140221_Resid.png",-.1,.1);

	///Examine BRs
//    cEvents theEvents;

//    cout << "Fermi BR: " << theEvents.calcBR(0.1,EKEMAX,0,1000000,true) << endl;  //0.0102108
//    cout << "No Fermi BR: " << theEvents.calcBR(0.1,EKEMAX,0,10000000,false) << endl; //0.00989184
//    cout << "No Fermi BR: " << theEvents.calcBR(20,EKEMAX,0,1000000,false) << endl; //0.00989184
//    cout << "No Fermi BR: " << theEvents.calcBR(0.1,EKEMAX,0,10000000,false) << endl; //0.00989184
//
//    const int numBR=1000;
//    double br[numBR];
//
//
//    double gammaLow=5;
//    double gammaHigh=EKEMAX;
//
//    for (int j = 1;j< 2;j++)
//    {
//
//
//        TH1D theHist("disHist","BR Distribution;BR;Measurements",100,0.0045*.9,0.0045*1.1);
//
//        for (int i = 0;i< numBR;i++)
//        {
//            br[i]=theEvents.calcBR(gammaLow,gammaHigh,0,100000,j,885.7);
//            theHist.Fill(br[i]);
//            cout << "BR: " << br[i] << endl;
//        }
//
//
//
//        TCanvas* c2 = new TCanvas("brPlot","brPlot",10,10,1024,768);
//
//        TString titleString="BR "+d2str(gammaLow)+"-"+d2str(gammaHigh,3)+" keV: " + d2str(theHist.GetMean(),7,1)+" #pm "+d2str(theHist.GetMeanError(),7,1)+" ("+d2str(theHist.GetMeanError()*100./theHist.GetMean(),4)+"%)";
//
//        theHist.SetTitle(titleString);
//
//        theHist.Draw();
//
//        TString imageName;
//        if(j)
//        {
//            imageName=GraphsDir+"BRTestFermi1403035eV_885.7.png";
//        }
//        else
//        {
//            imageName=GraphsDir+"BRTestNoFermi1403035keV_885.7.png";
//        }
//        c2->SaveAs(imageName);
//        delete c2;
//    }

///Exp analysis of electron triggers only 140305
//    doFunctionToAllFilesInDirectory("/media/mjbexternal/analysis/conversion/",&createRootFileFromExpDataITX);

//    makeExpEPCoincidenceFiles("PID8",NUM_EXP_SERIES_G8,EXP_SERIES_ID_LIST_G8);
//
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G8","EP",NUM_EXP_SERIES_G8,EXP_SERIES_ID_LIST_G8);
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G7_Cycle5_8","EP",NUM_EXP_SERIES_G7_Cycle5_8,EXP_SERIES_ID_LIST_G7_Cycle5_8);
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G8","E",NUM_EXP_SERIES_G8,EXP_SERIES_ID_LIST_G8);
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G7_Cycle5_8","E",NUM_EXP_SERIES_G7_Cycle5_8,EXP_SERIES_ID_LIST_G7_Cycle5_8);
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G8","EBGO",NUM_EXP_SERIES_G8,EXP_SERIES_ID_LIST_G8);
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G7_Cycle5_8","EBGO",NUM_EXP_SERIES_G7_Cycle5_8,EXP_SERIES_ID_LIST_G7_Cycle5_8);

//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G9","EBGO",NUM_EXP_SERIES_G9,EXP_SERIES_ID_LIST_G9);
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G9","E",NUM_EXP_SERIES_G9,EXP_SERIES_ID_LIST_G9);
//    makeAllStandardExpHists(getCurrentDateString(),"PID8","G9","EP",NUM_EXP_SERIES_G9,EXP_SERIES_ID_LIST_G9);

///Make E triggerMC

//    RDK2Set threeBodySet(48,"130429","STD_3B","RID278","RID279","","STD NonRad W Fermi");
//    RDK2Set fourBodySet(48,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID282","STD Rad W Fermi");
////    TH1* mcHist=threeBodySet.makeHist("eEnergy","eDChn.SBDEDepBlur",STD_EPCut,EDEPE_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist"+getCurrentDateString()+"_EDepE_EP_282.txt");
////    delete mcHist;
////    mcHist=threeBodySet.makeHist("eEnergy","eDChn.SBDEDepBlur",STD_ECut,EDEPE_DIM);
////    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist"+getCurrentDateString()+"_EDepE_E_282.txt");
////    delete mcHist;
//
//    TH1* mcHist=makeMCIndividualDetectorHist(&threeBodySet,&fourBodySet,"BGOEDepBlurLOAllExpParam",STD_ECut,"",BGO_DIM);
//    convertTH1ToTXT(mcHist,HistsDir+"MC_Hist"+getCurrentDateString()+"_EDepGWLO_AverageIndiv_EG_282.txt");
//    delete mcHist;

///Plot e trigger results

///EP cuts
//    TH1* expEEHists[2];
//    TH1* mcEEHists[1];
//    expEEHists[1]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140305_PID8_G7_Cycle5_8_EDepE_EP.txt",EDEPE_DIM);
//    expEEHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140306_PID8_G9_EDepE_EP.txt",EDEPE_DIM);
//
//    expEEHists[1]->SetTitle("EP Triggers (scal);Electron Energy Peak Channel (782 Channel=782 keV);E Trigger Counts");
//    expEEHists[0]->SetTitle("E Triggers;Electron Energy Peak Channel (782 Channel=782 keV);E Trigger Counts");
//
//    mcEEHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist140305_EDepE_EP_282.txt",EDEPE_DIM);
//
//    mcEEHists[0]->SetTitle("EP MC (scal);Electron Energy Peak Channel (782 Channel=782 keV);E Trigger Counts");
//
//    scaleHistogramsTogether(expEEHists[0],expEEHists[1],401,801);  //Scale hist2 to have the same number within the range as hist1
//    scaleHistogramsTogether(expEEHists[0],mcEEHists[0],401,801);  //Scale hist2 to have the same number within the range as hist1
//
//    mcEEHists[0]->GetYaxis()->SetRangeUser(0,4100.);
//
//    for (int i = 0;i< 2;i++)
//    {
//    	expEEHists[i]->Rebin(5);
//    }
//    mcEEHists[0]->Rebin(5);
//
//    plotExpVersusMCToImage(2,expEEHists,1,mcEEHists,"Electron Energy with EP Cut","EDepE_ETriggers_EPCut_140306.png");
//    delete expEEHists[0];
//    delete expEEHists[1];
//    delete mcEEHists[0];

	///E cut only
////    TH1* expEEHists[2];
////    TH1* mcEEHists[1];
//    expEEHists[1]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140305_PID8_G7_Cycle5_8_EDepE_E.txt",EDEPE_DIM);
//    expEEHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140306_PID8_G9_EDepE_E.txt",EDEPE_DIM);
//
//
//
//    expEEHists[1]->SetTitle("EP Triggers (scal);Electron Energy Peak Channel (782 Channel=782 keV);E Trigger Counts");
//    expEEHists[0]->SetTitle("E Triggers;Electron Energy Peak Channel (782 Channel=782 keV);E Trigger Counts");
//
//    double error;
//    double integral=expEEHists[0]->IntegralAndError(51,801,error);
//
//    mcEEHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist140305_EDepE_E_282.txt",EDEPE_DIM);
//
//
//
//    mcEEHists[0]->SetTitle("E MC (scal);Electron Energy Peak Channel (782 Channel=782 keV);E Trigger Counts");
//
//    scaleHistogramsTogether(expEEHists[0],expEEHists[1],401,801);  //Scale hist2 to have the same number within the range as hist1
//    scaleHistogramsTogether(expEEHists[0],mcEEHists[0],401,801);  //Scale hist2 to have the same number within the range as hist1
//
//    mcEEHists[0]->GetYaxis()->SetRangeUser(0,4100.);
//
////    for (int i = 0;i< 2;i++)
////    {
////    	expEEHists[i]->Rebin(5);
////    }
////    mcEEHists[0]->Rebin(5);
//
//    plotExpVersusMCToImage(2,expEEHists,1,mcEEHists,"Electron Energy with E Cut","EDepE_ETriggers_ECut_140306.png");
//    delete expEEHists[0];
//    delete expEEHists[1];
//    delete mcEEHists[0];

	///EG cut only
//    TH1* expHists[2];
//    TH1* mcHists[1];
//    expHists[1]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140305_PID8_G7_Cycle5_8_EDepG_BGOAverage_EG.txt",BGO_DIM);
//    expHists[0]=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140306_PID8_G9_EDepG_BGOAverage_EG.txt",BGO_DIM);
//
//    expHists[1]->SetTitle("EP Triggers;Electron Energy Peak Channel (782 Channel=782 keV);#frac(EPG}{EP}");
//    expHists[0]->SetTitle("E Triggers;Electron Energy Peak Channel (782 Channel=782 keV);#frac(EPG}{EP}");
//
//    mcHists[0]=getTabSeperatedHist(HistsDir+"MC_Hist140305_EDepGWLO_AverageIndiv_EG_282.txt",BGO_DIM);
//
//    mcHists[0]->SetTitle("E MC;BGO Energy Peak Channel (511 Channel=511 keV);#frac{EPG}{EP}");
//
//
//    for (int i = 0;i< 2;i++)
//    {
//    	expHists[i]->Rebin(5);
//    }
//    mcHists[0]->Rebin(5);
//
//    plotResidualsToImage(2,expHists,1,mcHists,"BGO Energy with EG Cut","BGO_ETriggers_EGCut_Resid_140306.png");
//
//    plotExpVersusMCToImage(2,expHists,1,mcHists,"BGO Energy with EG Cut","BGO_ETriggers_EGCut_140306.png");

///Convert BAPD Field 140309
//    TString txtPath="/home/mjbales/school/rdk/comsol/COMSOL_BAPDOnlyField_A_1350V_140307.txt";
//    TString rootPath="/home/mjbales/school/rdk/comsol/COMSOL_BAPDOnlyField_A_1350V_140307.root";
//    ROOTField newField(txtPath,"bapd",1,3,3);
//    newField.saveFieldToFile(rootPath,"bapd");

///Create runs with new BAPD field
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID286",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD");

///Create Mirror Field shifts
//    ///Shift Mirror +1 mm in z direction
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID287",PROTON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0 0 1 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID288",ELECTRON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0 0  1 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID289",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0 0  1 mm");
//
//    ///Shift Mirror -1 mm in z direction
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID290",PROTON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0 0 -1 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID291",ELECTRON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0 0  -1 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID292",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0 0  -1 mm");
//
//    ///Shift Mirror +1 mm in x direction
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID293",PROTON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1 0 0 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID294",ELECTRON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1 0 0 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID295",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1 0 0 mm");
//
//    ///Shift Mirror -1 mm in x direction
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID296",PROTON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1 0 0 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID297",ELECTRON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1 0 0 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID298",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1 0 0 mm");

///Create SBD Shifts (again)
//    ///Shift SBD Detector by 1 cm towards bend
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID299",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID300",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
//
//
//    ///Shift SBD Detector 1 cm away from bend
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID301",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID302",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm");
//
//
//   ///Shift SBD Detector by 1 mm side to side
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID303",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID304",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
//
//
//    ///Shift SBD Detector -1 mm side to side
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID305",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset -1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1,"RID306",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset -1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");
//
///Create 1.0 mm SBD run
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID307",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/use1_5mmSBD 0");

///Create Partial Length Contraction runs
//    ///No Al length contraction
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID308",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setAlLengthContractionFactor 1");
//
//    ///No BGO length contraction
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID309",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setBGOLengthContractionFactor 1");

///Check old length contraction no EP cuts
//    RDK2Set fourBodySet(48,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID282","STD Rad W Fermi");
//    fourBodySet.calcEPGRateFor12IndividualBGOs(STD_EPCut); //0.0149523
//    fourBodySet.calcEPGRateFor12IndividualBGOs(""); //0.0495034

//    RDK2Set fourBodySetB(48,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID283","STD Rad W Fermi");
//    fourBodySetB.calcEPGRateFor12IndividualBGOs(STD_EPCut); //0.0148503
//    fourBodySetB.calcEPGRateFor12IndividualBGOs(""); //0.0489479

///Make 2D graphs for Kevin for deconvolution
//    const int numFiles=48;
//    const int numEvents=numFiles*1000000;
//    RDK2Set fourBodySet(numFiles,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID283","STD Rad W Fermi");
//    TString histPath=HistsDir+"Hist2D_EDepGBGO_vs_EG0_140310.root";
//    TFile histFile(histPath,"RECREATE");
//    histFile.cd();
//    TH2D* hist2DAvg=makeTH2DFromDim("eDepAvgvseg0Hist","eDepAvgvseg0Hist",BGO_DIM,BGO_DIM);
//    hist2DAvg->Sumw2();
//    TH2D* hists2D[12];
//    for (int i = 0;i< 12;i++)
//    {
//        TString histName="eDepvseg0Hist"+int2str(i+1);
//    	hists2D[i]=fourBodySet.makeHist2D(histName,"(gChn.BGO"+int2str(i+1)+"EDepTotal)*(eDChn.SBDEDepBlur > 50. && eDChn.SBDEDepBlur < 800. && pChn.SBDEDepTotal > 10. && pChn.SBDEDepTotal < 31. && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 10.e-6):eg0","",BGO_DIM,BGO_DIM);
//    	hists2D[i]->Scale(1./numEvents);
//        hist2DAvg->Add(hists2D[i],1./12.);
//        hists2D[i]->Write("",TObject::kOverwrite);
//    }
//    hist2DAvg->Write("",TObject::kOverwrite);
//
////    TString histPath=HistsDir+"Hist2D_EDepGBGO_vs_EG0_140310.root";
////    TFile histFile(histPath,"READ");
////    TH2D* hist2DAvg;
////    TH2D* hists2D[12];
//
//    hist2DAvg=(TH2D*) histFile.Get("eDepAvgvseg0Hist");
//    convertTH2ToTXT(hist2DAvg,HistsDir+"Hist2D_Average_EDepGBGO_vs_EG0_140310.txt");
//
//    for (int i = 0;i< 12;i++)
//    {
//        TString histName="eDepvseg0Hist"+int2str(i+1);
//    	hists2D[i]=(TH2D*) histFile.Get(histName);
//    	convertTH2ToTXT(hists2D[i],HistsDir+"Hist2D_BGO"+int2str(i+1)+"_EDepGBGO_vs_EG0_140310.txt");
//    }
//    histFile.Close();

///Make proton energy cut histograms 140312
//    RDK2Set threeBodySet("stdset","stdset",1,"130429","STD_3B","RID278","RID279","");
//
//    const int numHists=6;
//    TString histName,drawString,cutString,histFilePath;
//
//
//    TH1* tempHist;
//    for (int i = 0;i< numHists;i++)
//    {
//
//        double lowEn=10+i*3;
//        double highEn=lowEn+3;
//        histName="Hist_p_En"+int2str(lowEn)+"_"+int2str(highEn);
//
////
////        ///pUnblur
////        cutString="eDChn.SBDEDepBlur > 50. && eDChn.SBDEDepBlur < 800. && pChn.SBDEDepTotal > "+int2str(lowEn)+". && pChn.SBDEDepTotal < "+int2str(highEn)+" && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 10.e-6";
////        //pE
////    	tempHist=threeBodySet.makeHist("temp","pChn.SBDEDepTotal",cutString,EDEPP_DIM);
////        histFilePath=HistsDir+"MC"+histName+"_EpUnblur_EPunblur_RID278_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////        //eE
////    	tempHist=threeBodySet.makeHist("temp","eDChn.SBDEDepBlur",cutString,EDEPE_DIM);
////        histFilePath=HistsDir+"MC"+histName+"_Ee_EPunblur_RID278_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////        //pToF
////    	tempHist=threeBodySet.makeHist("temp","pChn.SBDTimeFirst/(40.0e-9)",cutString,TOFPCCHAN_DIM);
////        histFilePath=HistsDir+"MC"+histName+"_pToF_EPunblur_RID278_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////        ///pBlur
////        cutString="eDChn.SBDEDepBlur > 50. && eDChn.SBDEDepBlur < 800. && pDChn.SBDEDepBlur > "+int2str(lowEn)+" && pDChn.SBDEDepBlur < "+int2str(highEn)+" && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 10.e-6";
////
////         //pE
////    	tempHist=threeBodySet.makeHist("temp","pDChn.SBDEDepBlur",cutString,EDEPP_DIM);
////        histFilePath=HistsDir+"MC"+histName+"_Ep_EP_RID278_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////        //eE
////    	tempHist=threeBodySet.makeHist("temp","eDChn.SBDEDepBlur",cutString,EDEPE_DIM);
////        histFilePath=HistsDir+"MC"+histName+"_Ee_EP_RID278_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////        //pToF
////    	tempHist=threeBodySet.makeHist("temp","pChn.SBDTimeFirst/(40.0e-9)",cutString,TOFPCCHAN_DIM);
////        histFilePath=HistsDir+"MC"+histName+"_pToF_EP_RID278_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////
////
////        ///ExpHists
////        cutString="SBDEDepE > 50 && SBDEDepE < 800 && SBDToFP >= 2e-6 && SBDToFP <= 10e-6 && SBDEDepP > "+int2str(lowEn)+" && SBDEDepP < "+int2str(highEn);
////        int seriesArray[1]={76};
////
////        //pE
////        tempHist=makeExpEPHist("PID8", 1, seriesArray,"SBDEDepP", cutString, EDEPP_DIM);
////        histFilePath=HistsDir+"Exp"+histName+"_Ep_EP_PID8_S76_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////        //eE
////        tempHist=makeExpEPHist("PID8", 1, seriesArray,"SBDEDepE", cutString, EDEPE_DIM);
////        histFilePath=HistsDir+"Exp"+histName+"_Ee_EP_PID8_S76_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
////
////        pToF
////        tempHist=makeExpEPHist("PID8", 1, seriesArray,"SBDToFP/(40e-9)", cutString, TOFPCCHAN_DIM);
////        histFilePath=HistsDir+"Exp"+histName+"_pToF_EP_PID8_S76_140312.txt";
////        convertTH1ToTXT(tempHist,histFilePath);
////        delete tempHist;
//
//        ///NoGoldLayer
//        RDK2Set threeBodyNoGoldSet("nogoldset","nogoldset",1,"130429","STD_3B","RID310","RID312","");
//        cutString="eDChn.SBDEDepBlur > 50. && eDChn.SBDEDepBlur < 800. && pChn.SBDKEFirst > "+int2str(lowEn)+" && pDChn.SBDKEFirst < "+int2str(highEn)+" && pChn.SBDKEFirst > 2.0e-6 && pChn.SBDTimeFirst < 10.e-6";
//        TCut theCut(cutString);
//         //pE
//    	tempHist=threeBodyNoGoldSet.makeHist("temp","pDChn.SBDEDepBlur",theCut,EDEPP_DIM);
//        histFilePath=HistsDir+"MC"+histName+"_Ep_EP_RID310_140312.txt";
//        convertTH1ToTXT(tempHist,histFilePath);
//        delete tempHist;
//
//        //eE
//    	tempHist=threeBodyNoGoldSet.makeHist("temp","eDChn.SBDEDepBlur",theCut,EDEPE_DIM);
//        histFilePath=HistsDir+"MC"+histName+"_Ee_EP_RID310_140312.txt";
//        convertTH1ToTXT(tempHist,histFilePath);
//        delete tempHist;
//
//        //pToF
//    	tempHist=threeBodyNoGoldSet.makeHist("temp","pChn.SBDTimeFirst/(40.0e-9)",theCut,TOFPCCHAN_DIM);
//        histFilePath=HistsDir+"MC"+histName+"_pToF_EP_RID310_140312.txt";
//        convertTH1ToTXT(tempHist,histFilePath);
//        delete tempHist;
//
//
//    }

///Make Graphs: Unblured proton energies
//    const int numHists=6;
//    TString histName,histFilePath;
//
//    TH1* mcPeHists[numHists];
//    TH1* mcEeHists[numHists];
//    TH1* mcTHists[numHists];
//    TH1* expPeHists[numHists];
//    TH1* expEeHists[numHists];
//    TH1* expTHists[numHists];
//    TH1* mcBlurPeHists[numHists];
//    TH1* mcBlurEeHists[numHists];
//    TH1* mcBlurTHists[numHists];
//
//    const int colorList[numHists]= {kBlack, kRed, kGreen, kBlue, kViolet, kGray};
//
//
//     ///GetHists
//     for (int i = 0;i< numHists;i++)
//     {
//        double lowEn=10+i*3;
//        double highEn=lowEn+3;
//        histName="Hist_p_En"+int2str(lowEn)+"_"+int2str(highEn);
//
//        mcPeHists[i]=getTabSeperatedHist(HistsDir+"MC"+histName+"_EpUnblur_EPunblur_RID278_140312.txt",EDEPP_DIM);
//        mcEeHists[i]=getTabSeperatedHist(HistsDir+"MC"+histName+"_Ee_EPunblur_RID278_140312.txt",EDEPE_DIM);
//        mcTHists[i]=getTabSeperatedHist(HistsDir+"MC"+histName+"_pToF_EPunblur_RID278_140312.txt",TOFPCCHAN_DIM);
//
//
//        mcBlurPeHists[i]=getTabSeperatedHist(HistsDir+"MC"+histName+"_Ep_EP_RID278_140312.txt",EDEPP_DIM);
//        mcBlurEeHists[i]=getTabSeperatedHist(HistsDir+"MC"+histName+"_Ee_EP_RID278_140312.txt",EDEPE_DIM);
//        mcBlurTHists[i]=getTabSeperatedHist(HistsDir+"MC"+histName+"_pToF_EP_RID278_140312.txt",TOFPCCHAN_DIM);
//
//        expPeHists[i]=getTabSeperatedHist(HistsDir+"Exp"+histName+"_Ep_EP_PID8_S76_140312.txt",EDEPP_DIM);
//        expEeHists[i]=getTabSeperatedHist(HistsDir+"Exp"+histName+"_Ee_EP_PID8_S76_140312.txt",EDEPE_DIM);
//        expTHists[i]=getTabSeperatedHist(HistsDir+"Exp"+histName+"_pToF_EP_PID8_S76_140312.txt",TOFPCCHAN_DIM);
//
//        mcPeHists[i]->SetLineColor(colorList[i]);
//        mcEeHists[i]->SetLineColor(colorList[i]);
//        mcTHists[i]->SetLineColor(colorList[i]);
//        expPeHists[i]->SetLineColor(colorList[i]);
//        expEeHists[i]->SetLineColor(colorList[i]);
//        expTHists[i]->SetLineColor(colorList[i]);
//        mcBlurPeHists[i]->SetLineColor(colorList[i]);
//        mcBlurEeHists[i]->SetLineColor(colorList[i]);
//        mcBlurTHists[i]->SetLineColor(colorList[i]);
//
//
//        expPeHists[i]->SetMarkerColor(colorList[i]);
//        expEeHists[i]->SetMarkerColor(colorList[i]);
//        expTHists[i]->SetMarkerColor(colorList[i]);
//
//        mcPeHists[i]->SetTitle(";Proton Energy (keV);Counts");
//        mcEeHists[i]->SetTitle(";Electron Energy (keV);Counts");
//        mcTHists[i]->SetTitle(";e-p TOF (1ch=40ns);Counts");
//        mcBlurPeHists[i]->SetTitle(";Proton Energy (keV);Counts");
//        mcBlurEeHists[i]->SetTitle(";Electron Energy (keV);Counts");
//        mcBlurTHists[i]->SetTitle(";e-p TOF (1ch=40ns);Counts");
//
//
//        expPeHists[i]->SetTitle(";Proton Energy (keV);Counts");
//        expEeHists[i]->SetTitle(";Electron Energy (keV);Counts");
//        expTHists[i]->SetTitle(";e-p TOF (1ch=40ns);Counts");
//        expEeHists[i]->Rebin(10);
//        expTHists[i]->Rebin(4);
//
//     }
//
	///Make duplicate Exp graph as Jeff did
//    TCanvas* theCanvas = new TCanvas("c","c",1024,2000);
//    theCanvas->Divide(1,3);
//
//    TLegend* theLegends[3];
//    TH1** theHists[3]={expPeHists,expEeHists,expTHists};
//
//
//    for (int j = 0;j< 3;j++) //graph in layout loop
//    {
//        double max=0;
//        theCanvas->cd(j+1);
//        theLegends[j] = new TLegend(0.8,0.5,0.93,0.9);
//        theLegends[j]->SetBorderSize(1);
//        theLegends[j]->SetFillColor(kWhite);
//        for (int i = 0;i< numHists;i++)  //PEnergy group loop
//        {
//            double lowEn=10+i*3;
//            double highEn=lowEn+3;
//            histName="Hist_p_En"+int2str(lowEn)+"_"+int2str(highEn);
//
//            //PEnergy
//
//            theLegends[j]->AddEntry(theHists[j][i],histName,"L");
//            if(i==0)
//            {
//                theHists[j][i]->Draw("hist");
//            }
//            else
//            {
//                theHists[j][i]->Draw("hist same");
//            }
//            if(j!=0 && i != 0)
//            {
//                scaleHistogramsTogether(theHists[j][0],theHists[j][i]);
//            }
//            if(theHists[j][i]->GetMaximum() > max)
//            {
//                max=theHists[j][i]->GetMaximum();
//            }
//        }
//        theLegends[j]->Draw();
//        theHists[j][0]->GetYaxis()->SetRangeUser(0,max);
//        theCanvas->Draw();
//
//    }
//
//    TString imagePath=GraphsDir+"Exp_pEnergyGroups_140312.pdf";
//    theCanvas->Print(imagePath);

	///Make MC graph Hist as Jeff did
//    TCanvas* theCanvas = new TCanvas("c","c",1024,2000);
//    theCanvas->Divide(1,3);
//
//    TLegend* theLegends[3];
//    TH1** theHists[3]={mcPeHists,mcEeHists,mcTHists};
//
//
//    for (int j = 0;j< 3;j++) //graph in layout loop
//    {
//        double max=0;
//        theCanvas->cd(j+1);
//        theLegends[j] = new TLegend(0.8,0.5,0.93,0.9);
//        theLegends[j]->SetBorderSize(1);
//        theLegends[j]->SetFillColor(kWhite);
//        for (int i = 0;i< numHists;i++)  //PEnergy group loop
//        {
//            double lowEn=10+i*3;
//            double highEn=lowEn+3;
//            histName="Hist_p_En"+int2str(lowEn)+"_"+int2str(highEn);
//
//            //PEnergy
//
//            theLegends[j]->AddEntry(theHists[j][i],histName,"L");
//            if(i==0)
//            {
//                theHists[j][i]->Draw("hist");
//            }
//            else
//            {
//                theHists[j][i]->Draw("hist same");
//            }
//            if(j!=0 && i != 0)
//            {
//                scaleHistogramsTogether(theHists[j][0],theHists[j][i]);
//            }
//            if(theHists[j][i]->GetMaximum() > max)
//            {
//                max=theHists[j][i]->GetMaximum();
//            }
//
//        }
//
//
//        theLegends[j]->Draw();
//        theHists[j][0]->GetYaxis()->SetRangeUser(0,max);
//        theCanvas->Draw();
//
//    }
//
//    TString imagePath=GraphsDir+"MC_pEnergyGroups_normalized_140312.pdf";
//    theCanvas->Print(imagePath);

	///Make MC with blur on SBD graph Hist as Jeff did
//    TCanvas* theCanvas = new TCanvas("c","c",1024,2000);
//    theCanvas->Divide(1,3);
//
//    TLegend* theLegends[3];
//    TH1** theHists[3]={mcBlurPeHists,mcBlurEeHists,mcBlurTHists};
//
//
//    for (int j = 0;j< 3;j++) //graph in layout loop
//    {
//        double max=0;
//        theCanvas->cd(j+1);
//        theLegends[j] = new TLegend(0.8,0.5,0.93,0.9);
//        theLegends[j]->SetBorderSize(1);
//        theLegends[j]->SetFillColor(kWhite);
//        for (int i = 0;i< numHists;i++)  //PEnergy group loop
//        {
//            double lowEn=10+i*3;
//            double highEn=lowEn+3;
//            histName="Hist_p_En"+int2str(lowEn)+"_"+int2str(highEn);
//
//            //PEnergy
//
//            theLegends[j]->AddEntry(theHists[j][i],histName,"L");
//            if(i==0)
//            {
//                theHists[j][i]->Draw("hist");
//            }
//            else
//            {
//                theHists[j][i]->Draw("hist same");
//            }
//
//            if(j!=0 && i != 0)
//            {
//                scaleHistogramsTogether(theHists[j][0],theHists[j][i]);
//            }
//            if(theHists[j][i]->GetMaximum() > max)
//            {
//                max=theHists[j][i]->GetMaximum();
//            }
//        }
//        theLegends[j]->Draw();
//        theHists[j][0]->GetYaxis()->SetRangeUser(0,max);
//        theCanvas->Draw();
//
//    }
//
//    TString imagePath=GraphsDir+"MC_BlurpEnergyGroups_scaled_140312.pdf";
//    theCanvas->Print(imagePath);

///Turn exp series back into asci for jeff

//    TString expPath=EXPANALYSIS_DIR;
//    expPath+="O_v8d10_76++.itx.root";
//    TFile* inpFile=new TFile(expPath,"read");
//    TTree* expTree=(TTree*) inpFile->Get("expTree");
//
//    double SBDEDepP;
//    expTree->SetBranchAddress("p_Area",&SBDEDepP);
//    ofstream outFile;
//    outFile.precision(15);
//    TString outPath=EXPANALYSIS_DIR;
//    outPath+="O_v8d10_76++_Simple.txt";
//    outFile.open(outPath);
//    outFile << "#EventNum\tp_En" << endl;
//    for (int i = 0;i< expTree->GetEntries();i++)
//    {
//        expTree->GetEntry(i);
//        outFile << i << "\t" << scientific << SBDEDepP << endl;
//
//    }
//
//    outFile.close();
//    inpFile->Close();

//    TString outPath=EXPANALYSIS_DIR;
//    outPath+="O_v8d10_76++.txt";
//    TTreePlayer* expTreePlayer=(TTreePlayer*) expTree->GetPlayer();
//    expTreePlayer->SetScanRedirect(kTRUE);
//    expTreePlayer->SetScanFileName(outPath);
//    expTree->Scan("*");
//    inpFile->Close();

///Test variable bin width histograms
////    TH1* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist140220_EDepGWLO_AverageIndiv_EPG_RID282.txt",BGO_DIM);
////
////    const int numBins=20;
////    int varBinLow[numBins];
////    double totalCounts=mcHist->Integral();
////    double goal=totalCounts/(double)numBins;
////    varBinLow[0]=1;
////
////    int currentBottomBin=0;
////    for (int i = 1;i<mcHist->GetNbinsX();i++)
////    {
////        if(mcHist->Integral(varBinLow[currentBottomBin],i) > goal)
////        {
////            currentBottomBin++;
////            if(currentBottomBin == numBins)
////            break;
////            varBinLow[currentBottomBin]=i;
////            cout << i << ",";
////        }
////
////    }
//
////    double BGOBinLowEdges10[11]={10,13,16,20,26,34,44,57,73,96,800};
//    double BGOBinLowEdges20[21]={10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,800};
//    TH1D* varBGOHistPrototype=new TH1D("varProtoType","varProtoType",20,BGOBinLowEdges20);
//
//    ///Make MC
////    RDK2Set threeBodySet(48,"130429","STD_3B","RID278","RID279","","STD NonRad W Fermi");
////    RDK2Set fourBodySet(48,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID282","STD Rad W Fermi");
////    TH1* mcHists[1]={makeMCIndividualDetectorHist(&threeBodySet,&fourBodySet,"BGOEDepBlurLOAllExpParam",STD_EPCut,"",varBGOHistPrototype)};
////    convertTH1ToTXT(mcHists[0],HistsDir+"mcVarHistTest.txt");unsigned int m_Counter
//
//
//    ///MakeExp
////    TH1* expHists[1]={makeExpAverageBGODetectorHist("PID8", NUM_EXP_SERIES_G7_Cycle5_8,EXP_SERIES_ID_LIST_G7_Cycle5_8,varBGOHistPrototype)};
////    convertTH1ToTXT(expHists[0],HistsDir+"expVarHistTest.txt");
//
//    TH1* mcHists[1]={getTabSeperatedHist(HistsDir+"mcVarHistTest.txt",varBGOHistPrototype)};
//    TH1* expHists[1]={getTabSeperatedHist(HistsDir+"expVarHistTest.txt",varBGOHistPrototype)};
////
//    plotResidualsToImage(1,expHists,1,mcHists,"Variable Bin Width Test", "VarBinWidthTest.png");

///Make no gold layer run
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID310",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/SBDGoldLayerOff 1");

///Testing user created TObject I/O  140316
///MC AID1 - Standard
//    TString mcFilePath=MISC_DIR;
//    mcFilePath+="AID1.root";
//
//    TFile mcFile(mcFilePath,"recreate");
//    mcFile.cd();
//    RDK2MCAnalysis* mcAnalysis=new RDK2MCAnalysis("AID1","Previous Standard",48,"RID278","RID279",48,"RID280","RID281","RID283",STD_EPCutSetA);
//    mcAnalysis->MakeAllStandardHists();
//    mcAnalysis->GetEP_pTHist()->Print();
//    mcFile.cd();
//    mcAnalysis->Write();
//    mcFile.Close();
//    delete mcAnalysis;

//    TFile testFile2(filePath,"read");
//    RDK2MCAnalysis* theMCAnalysis=(RDK2MCAnalysis*) testFile2.Get("AID1");
//    theMCAnalysis->PrintRates();
//    cout << theMCAnalysis->GetEP_pTHist()->GetTitle()<< ":"<<theMCAnalysis->GetEP_pTHist()->Integral() <<endl;

	///EID0 - S317
//    TString expFilePath=MISC_DIR;
//    TString EIDString="EID0";
//    TString groupIDString="S317";
//    expFilePath+=EIDString+".root";
//    TFile expFile(expFilePath,"recreate");
//    expFile.cd();
//    RDK2ExpAnalysis* expAnalysis=new RDK2ExpAnalysis(EIDString,"Previous Standard","PID8",groupIDString,STD_EPCutSetA,STD_BGO_WINDOWS_A,STD_BAPD_WINDOWS_A);
////    expAnalysis->CalcRates();
////    expAnalysis->PrintRates();
//    expAnalysis->MakeAllStandardHists();
//    expFile.cd();
//    expAnalysis->Write();
//    expFile.Close();
//    delete expAnalysis;

	///EID1 - G7_Cycle5_8
//    TString expFilePath=MISC_DIR;
//    TString EIDString="EID1";
//    TString groupIDString="G7_Cycle5_8";
//    expFilePath+=EIDString+".root";
//    TFile expFile(expFilePath,"recreate");
//    expFile.cd();
//    RDK2ExpAnalysis* expAnalysis=new RDK2ExpAnalysis(EIDString,"Previous Standard","PID8",groupIDString,STD_EPCutSetA,STD_BGO_WINDOWS_A,STD_BAPD_WINDOWS_A);
//    expAnalysis->MakeAllStandardHists();
//    expFile.cd();
//    expAnalysis->Write();
//    expFile.Close();
//    delete expAnalysis;

//    TFile testFile2(expFilePath,"read");
//    RDK2ExpAnalysis* theExpAnalysis=(RDK2ExpAnalysis*) testFile2.Get(EIDString);
//    double temp;
//    theExpAnalysis->CalcRates();
//    theExpAnalysis->PrintRates();
//

//    TH1* expHist=getTabSeperatedHist(HistsDir+"EXPMatt_Hist140221_PID8_G7_Cycle5_8_EDepG_BGOAverage_EPG.txt",BGO_DIM);
//    TH1* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist140220_EDepGWLO_AverageIndiv_EPG_RID282.txt",BGO_DIM);
//    cout << "Exp: "<< expHist->Integral() << endl;
//    cout << "MC: "<<mcHist->Integral() << endl;

///Make MC
//    TString mcFilePath,AIDString,descriptionString;
//    TFile* mcFile;
//    RDK2MCAnalysis* mcAnalysis;
/////MC AID1 - Standard
//    mcFilePath=MISC_DIR;
//    AIDString="AID1";
//    descriptionString="Standard";
//    mcFilePath+=AIDString+".root";
//
//    mcFile=new TFile(mcFilePath,"recreate");
//    mcFile->cd();
//    mcAnalysis=new RDK2MCAnalysis(AIDString,descriptionString,48,"RID278","RID279",48,"RID280","RID281","RID283",STD_EPCutSetA);
//    mcAnalysis->MakeAllStandardHists();
//    mcFile->cd();
//    mcAnalysis->Write();
//    mcFile->Close();
//    delete mcAnalysis;
//    delete mcFile;
//
/////MC AID2 - No Gold on SBD
//    mcFilePath=MISC_DIR;
//    AIDString="AID2";
//    descriptionString="No Gold on SBD";
//    mcFilePath+=AIDString+".root";
//
//    mcFile=new TFile(mcFilePath,"recreate");
//    mcFile->cd();
//    mcAnalysis=new RDK2MCAnalysis(AIDString,descriptionString,48,"RID310","RID279",48,"RID280","RID281","RID283",STD_EPCutSetA);
//    mcAnalysis->MakeAllStandardHists();
//    mcFile->cd();
//    mcAnalysis->Write();
//    mcFile->Close();
//    delete mcAnalysis;
//    delete mcFile;

	///Test plotting
//    TString EIDString="EID1";
//    TFile expFile(TString(MISC_DIR)+EIDString+".root","read");
//    RDK2ExpAnalysis* theExp=(RDK2ExpAnalysis*) expFile.Get(EIDString);
//    RDK2ExpAnalysis* theExps[2]={theExp,theExp};

//    double inp,error;
//    double epCounts=theExp->GetEPCounts();
//    inp=theExp->GetEPG_eEHist()->IntegralAndError(1,800,error);
//    theExp->SetEPGCounts(inp*epCounts,error*epCounts);
//
//
//    theExp->CalcEPGPerEP(DET_EPG,error);
//
//    theExp->PrintRates();
//
//    TFile exp2File(TString(MISC_DIR)+EIDString+"b.root","recreate");
//    exp2File.cd();
//    theExp->Write();
//    exp2File.Close();

//    TString AIDString="AID1";
//    TFile mcFile(TString(MISC_DIR)+AIDString+".root","read");
//    RDK2MCAnalysis* theMC=(RDK2MCAnalysis*) mcFile.Get(AIDString);
//    AIDString="AID2";
//    TFile mcFile2(TString(MISC_DIR)+AIDString+".root","read");
//    RDK2MCAnalysis* theMC2=(RDK2MCAnalysis*) mcFile2.Get(AIDString);
//    RDK2MCAnalysis* theMCs[2]={theMC,theMC2};

//    TH1* expHists[1]={theExp->GetEPG_gEAvgVarHist()};
//    TH1* mcHists[1]={theMC->GetEPG_gEAvgVarHist()};
//    plotResidualsToImage(1,expHists,1,mcHists,"test","test.png");

//    RDK2AnalysisPlotter theComparison("Plotting Test",0,theExps,1,theMCs);
//    TCanvas theCanvas("EPExpMCAnalysisComparisonPlot","EPExpMCAnalysisComparisonPlot",10,10,1200,800);
//    theComparison.makeAll();
//    expFile.Close();
//    delete theExp;
//    theExp=NULL;
//    theExps[1]=NULL;

	///Make a standard set of RID after event gen changes 3/19/14
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID311",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID312",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID313",PROTON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID314",ELECTRON,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID315",GAMMA,FOURBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT");

///Make analysis files
//    TString mcFilePath,AIDString,descriptionString;
//    TFile* mcFile;
//    RDK2MCAnalysis* mcAnalysis;
//    TString expFilePath, EIDString,groupIDString;
//    TFile* expFile;
//    RDK2ExpAnalysis* expAnalysis;
//
//    ///EID2 - G7_Cycle3_4
//    EIDString="EID2";
//    groupIDString="G7_Cycle3_4";
//    descriptionString="Previous Standard";
//
//    expFilePath=MISC_DIR;
//    expFilePath+=EIDString+".root";
//    expFile=new TFile(expFilePath,"recreate");
//    expFile->cd();
//    expAnalysis=new RDK2ExpAnalysis(EIDString,descriptionString,"PID8",groupIDString,STD_EPCutSetA,STD_BGO_WINDOWS_A,STD_BAPD_WINDOWS_A);
//    expAnalysis->MakeAllStandardHists();
//    expFile->cd();
//    expAnalysis->Write();
//    expFile->Close();
//    delete expAnalysis;
//    delete expFile;
//
//    ///MC AID3 - Electron Cut >100
//    AIDString="AID3";
//    descriptionString="Standard eE>100";
//
//    mcFilePath=MISC_DIR;
//    mcFilePath+=AIDString+".root";
//    mcFile=new TFile(mcFilePath,"recreate");
//    mcFile->cd();
//    mcAnalysis=new RDK2MCAnalysis(AIDString,descriptionString,48,"RID278","RID279",48,"RID280","RID281","RID283",STD_EPCutSetB);
//    mcAnalysis->MakeAllStandardHists();
//    mcFile->cd();
//    mcAnalysis->Write();
//    mcFile->Close();
//    delete mcAnalysis;
//    delete mcFile;
//
//
//
/////MC AID4 - Old Fields, no contraction
//    AIDString="AID4";
//    descriptionString="Old Fields, no contract.";
//
//    mcFilePath=MISC_DIR;
//    mcFilePath+=AIDString+".root";
//    mcFile=new TFile(mcFilePath,"recreate");
//    mcFile->cd();
//    mcAnalysis=new RDK2MCAnalysis(AIDString,descriptionString,48,"RID253","RID254",48,"RID255","RID256","RID257",STD_EPCutSetA);
//    mcAnalysis->MakeAllStandardHists();
//    mcFile->cd();
//    mcAnalysis->Write();
//    mcFile->Close();
//    delete mcAnalysis;
//    delete mcFile;
//
//      ///MC AID5 - No BGO length contraction
//    AIDString="AID5";
//    descriptionString="No BGO contraction";
//
//    mcFilePath=MISC_DIR;
//    mcFilePath+=AIDString+".root";
//    mcFile=new TFile(mcFilePath,"recreate");
//    mcFile->cd();
//    mcAnalysis=new RDK2MCAnalysis(AIDString,descriptionString,48,"RID278","RID279",48,"RID280","RID281","RID282",STD_EPCutSetA);
//    mcAnalysis->MakeAllStandardHists();
//    mcFile->cd();
//    mcAnalysis->Write();
//    mcFile->Close();
//    delete mcAnalysis;
//    delete mcFile;

///Make plots for existing ids 140312
//    TString expList="EID1 EID2";
//    TString mcList="AID1 AID3 AID4 AID5 AID14 AID15";
//
//    for (int i = 0;i< numItemsInStringList(expList);i++)
//    {
//    	RDK2AnalysisPlotter thePlotter(stringFromList(expList,i),stringFromList(expList,i),"");
//        thePlotter.makeAll();
//    }
//    for (int i = 0;i< numItemsInStringList(mcList);i++)
//    {
//    	RDK2AnalysisPlotter thePlotter(stringFromList(mcList,i),"",stringFromList(mcList,i));
//        thePlotter.makeAll();
//    }

///Big Analyze and plot job
///Done
//    createMCAnalysisAndPlot( "AID14",  "SBD shift +x 1mm (3 body only)", 1,  "RID303", "RID304",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID15",  "SBD shift -x 1mm (3 body only)", 1,  "RID305", "RID306",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createExpAnalysisAndPlot( "EID3",  "Standard", "PID8", "G7_Cycle5_8", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//    createMCAnalysisAndPlot( "AID2",  "No Gold on SBD", 1,  "RID310", "RID312",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID6",  "Standard (new event sequence)", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID8",  "Mirror Field +z 1mm (4 body only)", 1,  "RID311", "RID312",  1, "RID287",  "RID288", "RID289",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID9",  "Mirror Field -z 1mm (4 body only)", 1,  "RID311", "RID312",  1, "RID290",  "RID291", "RID292",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID10",  "Mirror Field +x 1mm (4 body only)", 1,  "RID311", "RID312",  1, "RID293",  "RID294", "RID295",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID11",  "Mirror Field -x 1mm (4 body only)", 1,  "RID311", "RID312",  1, "RID296",  "RID297", "RID298",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID16",  "No Al Length Contraction (G only)", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID308",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID17",  "No BGO Length Contraction (G only)", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID309",  STD_EPCutSetB);

///Todo
//
///Problems

///Fix Experiments EP counts
//    TString expList="EID3";
//    for (int i = 0;i< numItemsInStringList(expList);i++)
//    {
//        TString analString=stringFromList(expList,i);
//        TString filePath=ANALYSIS_DIR;
//        filePath+=analString+".root";
//        TFile f(filePath);
//        RDK2ExpAnalysis* theAnalysis=(RDK2ExpAnalysis*) f.Get(analString);
//        double tempError;
//        theAnalysis->CalcEPGPerEP(DET_EPG,tempError);
//        theAnalysis->PrintRates();
//        TFile f2(filePath+".grrr","recreate");
//        f2.cd();
//        theAnalysis->Write();
//        f2.Close();
//        f.Close();
//    }

//    RDK2AnalysisPlotter thePlotter5("EID1","EID1",""); thePlotter5.makeAll();
//RDK2ExpAnalysis theAnalysis( "S76",  "Standard", "PID8", "S76", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//double tempError;
//theAnalysis.CalcRates();
//theAnalysis.PrintRates();

///Plot pE section hists 140326

//    RDK2Set threeBodySet("standard","standard",48,"130429","STD_3B","RID311","RID304","");
//    createAndPlotPESections("pE",&threeBodySet, 6, 17, 25.75);
//    createAndPlotPESections("pEWRes",&threeBodySet, 6, 10, 31,1);
//    RDK2Set threeBodyNoGoldSet("nogoldset","nogoldset",1,"130429","STD_3B","RID310","RID312","");
//    createAndPlotPESections("pE_No_Gold",&threeBodyNoGoldSet, 6, 24.75, 25.75);
//    createAndPlotPESections("pEWRes_No_Gold",&threeBodyNoGoldSet, 6, 18, 31,1);
//    createAndPlotPESections("pEKE_No_Gold",&threeBodyNoGoldSet, 6, 24.75,25.75,2);

//     createAndPlotPESections("ExpCycles5_8",&threeBodySet, 6, 10, 31,0,true); //Make Exp

///Future to Be Analyzed
//    createMCAnalysisAndPlot( "AID7",  "BAPD Field (proton three body only)", 1,  "RID286", "RID312",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID12",  "SBD shift upstream 1cm (3 body only)", 1,  "RID299", "RID300",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID13",  "SBD shift downstream 1cm (3 body only)", 1,  "RID301", "RID302",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//RDK2AnalysisPlotter thePlotter3("BAPD Field from COMSOL","EID3","AID3 AID6"); thePlotter3.makeEPPlot(DET_EP);

///Checkunusual EPG/EP rate calc 140327
//RDK2ExpAnalysis theAnalysis( "S87",  "Standard", "PID8", "S87", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//theAnalysis.MakeAllStandardHists();

///Redo Experimental Analysis
//    createExpAnalysisAndPlot( "EID1",  "Standard", "PID8", "G7_Cycle5_8", STD_EPCutSetA, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//    createExpAnalysisAndPlot( "EID2",  "Standard", "PID8", "G7_Cycle3_4", STD_EPCutSetA, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//    createExpAnalysisAndPlot( "EID3",  "Standard", "PID8", "G7_Cycle5_8", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);

///Comparisons
//    RDK2AnalysisPlotter thePlotter1("SBD shift x direction","EID3","AID3 AID10 AID11"); thePlotter1.makeEPGPlot(DET_EPG);
//    RDK2AnalysisPlotter thePlotter2("BGO Length Contractions (Standard=With Both)","EID3","AID3 AID16 AID17"); thePlotter2.makeEPGPlot(DET_EPG);
//    RDK2AnalysisPlotter thePlotter4("No Gold for Protons","EID3","AID3 AID2"); thePlotter4.makeEPPlot(DET_EP);
//    RDK2AnalysisPlotter thePlotter5("Electron Energy Cut Change MC","EID3","AID1 AID3"); thePlotter5.makeAll();
//    RDK2AnalysisPlotter thePlotter6("Electron Energy Cut 50 keV","EID1","AID1"); thePlotter6.makeAll();
//    RDK2AnalysisPlotter thePlotter7("Electron Energy Cut 100 keV","EID3","AID3"); thePlotter7.makeAll();
//
//    createMCAnalysisAndPlot( "AID2",  "No Gold on SBD", 48,  "RID310", "RID312",  48, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID6",  "Standard (new event sequence)", 48,  "RID311", "RID312",  48, "RID313",  "RID314", "RID315",  STD_EPCutSetB);

///S233 only Etrigger analysis
//    createExpAnalysisAndPlot( "S233_E",  "S233_E", "PID8", "S233", STD_ECutSetA, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//    createExpAnalysisAndPlot( "S233_EP",  "S233_EP", "PID8", "S233", STD_EPCutSetA, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//    createMCAnalysisAndPlot( "AID1_E",  "Standard E TriggerOnly", 1,"RID278","RID279",1,"RID280","RID281","RID283",STD_ECutSetA);

//    RDK2AnalysisPlotter thePlotter1("S233 E Triggers E","S233_E","AID1_E"); thePlotter1.makeEPPlot(DET_EP);
//    RDK2AnalysisPlotter thePlotter2("S233 E Triggers EP","S233_EP","AID1"); thePlotter2.makeEPPlot(DET_EP);

///Process
//    createMCAnalysisAndPlot( "AID18",  "No BGO Length Contraction in anything", 48,  "RID278", "RID279",  48, "RID280",  "RID281", "RID282",  STD_EPCutSetB);

///Recreate BAPD runs with longer time 140330
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID286",PROTON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD");

///Create some plot comparisons for elogs
//RDK2AnalysisPlotter thePlotter1("Standard","EID3","AID3"); thePlotter1.makeAll();

///Recreate plots of Exp analysis
//TString expPlotList="EID1 EID2 EID3";
//for (int i = 0;i< numItemsInStringList(expPlotList);i++)
//{
//    TString analString=stringFromList(expPlotList,i);
//    RDK2AnalysisPlotter thePlotter1(analString,analString,""); thePlotter1.makeAll();
//}

///Check Length contraction results
//RDK2MCAnalysis* theAnalysis[4];
//theAnalysis[0]=new RDK2MCAnalysis("AID18noep",  "No Contract", 1,  "RID278", "RID279",  1, "RID280",  "RID281", "RID282",  STD_NoEPCutSetA);
//theAnalysis[1]=new RDK2MCAnalysis("AID16noep",  "BGO Contract", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID308",  STD_NoEPCutSetA);
//theAnalysis[2]=new RDK2MCAnalysis("AID17noep",  "Alum Contract", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID309",  STD_NoEPCutSetA);
//theAnalysis[3]=new RDK2MCAnalysis("AID6noep",  "BGO+Al Contract", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID315",  STD_NoEPCutSetA);
//
//TH1* theHists[4];
//for (int i = 3;i< 4;i++)
//{
//	theHists[i]=theAnalysis[i]->MakeRadAvgHist("hist"+int2str(i),DET_EPG,"eDChn.SBDEDepBlur",0,EDEPE_DIM);
//
//	cout<< theAnalysis[i]->GetTitle()<<": " << theHists[i]->Integral() << endl;
//}

//for (int i = 0;i< 4;i++)
//{
//	cout<< theAnalysis[i]->GetTitle()<<": " << theHists[i]->Integral() << endl;
//	delete theHists[i];
//	delete theAnalysis[i];
//}

//Results:
//No Contract: 4.22393e-05
//BGO Contract: 4.04982e-05
//Alum Contract: 4.41954e-05
//BGO+Al Contract: 4.21974e-05

///Now with ep's just to confirm
//RDK2MCAnalysis* theAnalysis[4];
//theAnalysis[0]=new RDK2MCAnalysis("AID18",  "No Contract", 1,  "RID278", "RID279",  1, "RID280",  "RID281", "RID282",  STD_EPCutSetB);
//theAnalysis[1]=new RDK2MCAnalysis("AID16",  "BGO Contract", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID308",  STD_EPCutSetB);
//theAnalysis[2]=new RDK2MCAnalysis("AID17",  "Alum Contract", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID309",  STD_EPCutSetB);
//theAnalysis[3]=new RDK2MCAnalysis("AID6",  "BGO+Al Contract", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//
//TH1* theHists[4];
//for (int i = 0;i< 4;i++)
//{
//	theHists[i]=theAnalysis[i]->MakeRadAvgHist("hist"+int2str(i),DET_EPG,"eDChn.SBDEDepBlur",0,EDEPE_DIM);
//
//	cout<< theAnalysis[i]->GetTitle()<<": " << theHists[i]->Integral() << endl;
//}
//
//for (int i = 0;i< 4;i++)
//{
//	cout<< theAnalysis[i]->GetTitle()<<": " << theHists[i]->Integral() << endl;
//	delete theHists[i];
//	delete theAnalysis[i];
//}

//Results
//No Contract: 8.1747e-05
//BGO Contract: 7.79928e-05
//Alum Contract: 8.48474e-05
//BGO+Al Contract: 8.14293e-05

///Comparisons
//    RDK2AnalysisPlotter thePlotter1("Mirror shift x direction","EID3","AID6 AID10 AID11"); thePlotter1.makeEPGPlot(DET_EPG);
//    RDK2AnalysisPlotter thePlotter2("Mirror shift z direction","EID3","AID6 AID8 AID9"); thePlotter2.makeEPGPlot(DET_EPG);
//    RDK2AnalysisPlotter thePlotter3("Standard","EID3","AID6"); thePlotter3.makeAll();

///Create fixed 1.0 mm SBD run
//    makeGeantMacroFilesForClusterEventGenerator(0,48,"RID307",ELECTRON,THREEBODY,"/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/use15mmSBD 0");

///Process SBD Shift and plot
//    createMCAnalysisAndPlot( "AID12",  "SBD shift upstream 1cm (3 body only)", 1,  "RID299", "RID300",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID13",  "SBD shift downstream 1cm (3 body only)", 1,  "RID301", "RID302",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    RDK2AnalysisPlotter thePlotter1("SBD shift x' direction","EID3","AID6 AID14 AID15"); thePlotter1.makeEPGPlot(DET_EPG);
//    RDK2AnalysisPlotter thePlotter2("SBD shift z' direction","EID3","AID6 AID12 AID13"); thePlotter2.makeEPGPlot(DET_EPG);

///Make four body SBD Shift runs
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID316", PROTON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID317", ELECTRON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID318", PROTON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID319", ELECTRON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID320", PROTON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID321", ELECTRON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID322", PROTON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset -1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID323", ELECTRON, "STD_B4_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset -1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm");

///Check BAPD Field Run
//createMCAnalysisAndPlot( "AID7",  "BAPD Field (proton three body only)", 1,  "RID286", "RID312",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//RDK2AnalysisPlotter thePlotter3("BAPD Field (Three Body Only)","EID3","AID6 AID7"); thePlotter3.makeAll(); //2/3 of protons are stopped.....great.....

///Make Comparison Graph with SRIM
//    TFile directFile(MRKGEANTDATA_DIR+TString("directProtonsNoFieldTest.root"));
//    TTree* resultsTree=(TTree*) directFile.Get("g");
//    gROOT->cd();
//    TH1D* geantHist=new TH1D("geantHist","Geant4;Energy(keV);Normalized counts",250,0,25);
//    resultsTree->Draw("SBDKEFirst >> geantHist");
//    geantHist->Sumw2();
//    directFile.Close();
//    TH1D* srimHist=getTabSeperatedHist(HistsDir+"Hist_SRIM_ProtonAt0Deg_Jeff140313.txt",250,0,25);
//    srimHist->SetTitle("SRIM");
//    TH1* theHists[2]={geantHist,srimHist};
//    scaleHistogramsTogether(geantHist,srimHist,2,251);
//    geantHist->GetXaxis()->SetRangeUser(10,25);
//    plotExpVersusMCToImage(0,NULL,2,theHists,"25 keV protons at 0 degrees through 20 nm AU","Geant4vsSRIM_0degProton20nmGold_2.png");

///Make high stat pE Section hist
//    RDK2Set threeBodySet("standard","standard",48,"130429","STD_3B","RID311","RID312","");
//    createAndPlotPESections("pE",&threeBodySet, 6, 17, 25.75);
//    createAndPlotPESections("pEWRes",&threeBodySet, 6, 10, 31,1);

///Check EDepBlur 140404
//RDK2Set threeBodySet("standard","standard",1,"130429","STD_3B","RID310","RID312","");
//TCanvas c2;
//c2.SetGrid(1,1);
//gPad->SetTickx(1);
//gPad->SetTicky(1);
//TH1* hist=threeBodySet.makeHist("hist","pDChn.SBDEDepBlur:pChn.SBDEDepTotal",TCut(),EDEPP_DIM,EDEPP_DIM);
//hist->SetTitle("Proton Energy;NoRes (keV);WithRes (keV)");
//hist->Draw();
//c2.SaveAs(GraphsDir+"NoGold_WResvsNoRes.png");

///Process exp with correct calibration on P and E (PID9)
//    doFunctionToAllFilesInDirectory("/media/mjbexternal/analysis/conversion/",&createRootFileFromExpDataITX);
//    makeExpEPCoincidenceFiles("PID9",NUM_EXP_SERIES_G7,EXP_SERIES_ID_LIST_G7);

//createExpAnalysisAndPlot( "EID4",  "New PE cal From Jeff", "PID9", "G7_Cycle5_8", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//RDK2AnalysisPlotter thePlotter3("New Exp","EID4 EID3","AID6"); thePlotter3.makeAll();

///Make series by series 10-13.5 keV plots

//    RDK2CutSet LowPeCut("LowPeCut","Low Pe Cuts",100.,800.,10.,13.5,2.e-6,10.e-6,10.,800.,0.1,20.);
//    TCanvas c2("plotManyHistCanvas","plotManyHistCanvas",1600,1200);
//    int rainbow[6]={kRed,kOrange,kYellow,kGreen,kBlue,kViolet};
//    c2.SetGrid(1,1);
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    TH1* theHists[NUM_EXP_SERIES_G7];
//    double max=0;
//    int currentCycle=0;
//    TLegend theLegend(0.55,.6,0.87,0.8);  //Top Right
//
//    for (int i = 0;i< NUM_EXP_SERIES_G7 ;i++)
//    {
//
//        TString seriesName="S"+int2str(EXP_SERIES_ID_LIST_G7[i]);
//        TString histPath=HistsDir+seriesName+"Specialb.txt";
//
//        RDK2ExpAnalysis expAnal(seriesName, seriesName, "PID9",seriesName, LowPeCut, STD_BGO_WINDOWS_A,STD_BAPD_WINDOWS_A);
//        theHists[i]=expAnal.MakeHist(seriesName+"Hist","SBDEDepE",LowPeCut.CreateExpEPCut(),EDEPE_DIM);
//        convertTH1ToTXT(theHists[i],histPath);
//
//        theHists[i]=getTabSeperatedHist(histPath,EDEPE_DIM);
//
//
//        if(getCycleSeriesIn(EXP_SERIES_ID_LIST_G7[i])!=currentCycle)
//        {
//            currentCycle=getCycleSeriesIn(EXP_SERIES_ID_LIST_G7[i]);
//            theLegend.AddEntry(theHists[i],"Cycle"+int2str(currentCycle),"L");
//
//        }
//        theHists[i]->SetLineColor(rainbow[getCycleSeriesIn(EXP_SERIES_ID_LIST_G7[i])-3]);
//        theHists[i]->SetTitle("Electron Energy for Proton Energies 10-13.5 keV;Electron Energy (keV);Relative Counts (AU)");
//        theHists[i]->Rebin(50);
//        scaleHistTo(theHists[i],1);
//        if(theHists[i]->GetMaximum() > max) max=theHists[i]->GetMaximum();
//    }
//    theHists[0]->GetYaxis()->SetRangeUser(0,max);
//    for (int i = 0;i< NUM_EXP_SERIES_G7;i++)
//    {
//        if(i==0)
//        {
//            theHists[i]->Draw("l hist");
//        }
//        else
//        {
//            theHists[i]->Draw("l hist same");
//        }
//    }
//    TH1* mcHist=getTabSeperatedHist(HistsDir+"MCHist_p_En_0_Ee_EPWRes_RID311_140402.txt",EDEPE_DIM);
//    mcHist->Rebin(50);
//    mcHist->SetBinContent(1,0);
//    mcHist->SetBinContent(2,0);
//    scaleHistTo(mcHist,1);
//    mcHist->SetLineColor(kBlack);
//    mcHist->Draw("l hist same");
//    theLegend.AddEntry(mcHist,"MC","L");
//    theLegend.SetBorderSize(1);
//    theLegend.SetFillColor(kWhite);
//    theLegend.Draw();
//
//    c2.SaveAs(GraphsDir+"SeriesBySeriesLowPE_EESpectrumb.pdf");
//
//    for (int i = 0;i< NUM_EXP_SERIES_G7;i++)
//    {
//        delete theHists[i];
//    }

///Create systematic runs 140411
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID324", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //Corrected Length contraction coefficients
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID325", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID326", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID327", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID328", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID329", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID330", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID331", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID332", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID333", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.  0. 0. mm"); //1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID334", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.  0. 0. mm"); //1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID335", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.  0. 0. mm"); //1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID336", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.  0. 0. mm"); //1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID337", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.  0. 0. mm"); //1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID338", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.  0. 0. mm"); //1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID339", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.  0. 0. mm"); //1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID340", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.  0. 0. mm"); //1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID341", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. 1. mm"); //1 mm shift z+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID342", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. 1. mm"); //1 mm shift z+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID343", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. 1. mm"); //1 mm shift z+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID344", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. 1. mm"); //1 mm shift z+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID345", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. -1. mm"); //1 mm shift z- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID346", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. -1. mm"); //1 mm shift z- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID347", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. -1. mm"); //1 mm shift z- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID348", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. -1. mm"); //1 mm shift z- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID349", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID350", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID351", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID352", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID353", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID354", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID355", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID356", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID316", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm"); //1 cm shift SBD down stream
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID317", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.1650476059 0 0.9862856015 cm"); //1 cm shift SBD down stream
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID318", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm"); //1 cm shift SBD up stream
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID319", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1 cm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.1650476059 -0 0.9862856015 cm"); //1 cm shift SBD up stream
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID320", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm"); //1 mm shift SBD +x
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID321", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.9862856015 0 -0.1650476059 mm"); //1 mm shift SBD +x
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID322", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset -1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm"); //1 mm shift SBD -x
//    makeGeantMacroFilesForClusterEventGenerator(0,1, "RID323", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset -1 0 0 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.9862856015 0 0.1650476059 mm"); //1 mm shift SBD -x
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID357", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/useBGOPositionalLightOut 0"); //Positional Light Output Off

///Check on 1mm SBD run
//    createMCAnalysisAndPlot( "AID19",  "1mm SBD", 1,  "RID311", "RID307",  1, "RID313",  "RID314", "RID315",  STD_EPCutSetB);
//    RDK2AnalysisPlotter thePlotter3("BAPD Field (Three Body Only)","EID3","AID6 AID19"); thePlotter3.makeAll();

///Dead layer
///Make Comparison Graph with SRIM
//    TFile directFile(MRKGEANTDATA_DIR+TString("directProtonsNoFieldTest.root"));
//    TTree* resultsTree=(TTree*) directFile.Get("g");
//    gROOT->cd();
//    TH1D* geantHist=new TH1D("geantHist","Geant4;Energy(keV);Normalized counts",250,0,25);
//    resultsTree->Draw("SBDKEFirst >> geantHist");
//    geantHist->Sumw2();
//    directFile.Close();
//    TH1D* srimHist=getTabSeperatedHist(HistsDir+"Hist_SRIM_ProtonAt0Deg_Jeff140313.txt",250,0,25);
//    srimHist->SetTitle("SRIM");
//    TH1* theHists[2]={geantHist,srimHist};
//    scaleHistogramsTogether(geantHist,srimHist,2,251);
//    geantHist->GetXaxis()->SetRangeUser(10,25);
//    plotExpVersusMCToImage(0,NULL,2,theHists,"25 keV protons at 0 degrees through 20 nm AU and 20 nm deadlayer","Geant4vsSRIM_0degProton20nmGold20nmDeadLayer.png");

///Make a whole bunch of analysis 140414

//    createMCAnalysisAndPlot( "AID19",  "1mm SBD", 1,  "RID311", "RID307",  1, "RID313",  "RID314", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID21",  "SBD shift upstream 1cm", 1,  "RID299", "RID300",  1, "RID316",  "RID317", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID22",  "SBD shift downstream 1cm", 1,  "RID301", "RID302",  1, "RID318",  "RID319", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID23",  "SBD shift +x 1mm", 1,  "RID303", "RID304",  1, "RID320",  "RID321", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID24",  "SBD shift -x 1mm", 1,  "RID305", "RID306",  1, "RID322",  "RID323", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID25",  "1 mm shift z+ B Field", 1,  "RID325", "RID326",  1, "RID327",  "RID328", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID26",  "1 mm shift z- B Field", 1,  "RID329", "RID330",  1, "RID331",  "RID332", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID27",  "1 mm shift x+ B Field", 1,  "RID333", "RID334",  1, "RID335",  "RID336", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID28",  "1 mm shift x- B Field", 1,  "RID337", "RID338",  1, "RID339",  "RID340", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID29",  "1 mm shift z+ neutronbeam", 1,  "RID341", "RID342",  1, "RID343",  "RID344", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID30",  "1 mm shift z- neutronbeam", 1,  "RID345", "RID346",  1, "RID347",  "RID348", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID31",  "1 mm shift x+ neutronbeam", 1,  "RID349", "RID350",  1, "RID351",  "RID352", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID32",  "1 mm shift x- neutronbeam", 1,  "RID353", "RID354",  1, "RID355",  "RID356", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID33",  "Positional Light Output Off", 48,  "RID311", "RID312",  48, "RID313",  "RID314", "RID357",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID20",  "Corrected Length Contraction", 48,  "RID311", "RID312",  48, "RID313",  "RID314", "RID324",  STD_EPCutSetB);

//    RDK2AnalysisPlotter* thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("SBD shift x' direction","EID3","AID20 AID23 AID24"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("SBD shift z' direction","EID3","AID20 AID21 AID22"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("B Field shift z' direction","EID3","AID20 AID25 AID26"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("B Field shift x' direction","EID3","AID20 AID27 AID28"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Neutron Beam shift z' direction","EID3","AID20 AID29 AID30"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Neutron Beam shift x' direction","EID3","AID20 AID31 AID32"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Positional Lightoutput effect","EID3","AID20 AID33"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Corrected Length Contraction","EID3","AID17 AID20"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("1 mm SBD","EID3","AID20 AID19"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("SBD shift x' direction","EID3","AID20 AID23 AID24"); thePlotter->makeEPPlot(DET_EP); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("SBD shift z' direction","EID3","AID20 AID21 AID22"); thePlotter->makeEPPlot(DET_EP); delete thePlotter;

///Dead layer runs
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID358", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //20 nm dead layer
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID359", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //20 nm dead layer
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID360", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //20 nm dead layer
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID361", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //20 nm dead layer

///Need to run gamma for neutron beam shift
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID362", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. 1. mm"); //1 mm shift z+ neutronbeam
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID363", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 0. 0. -1. mm"); //1 mm shift z- neutronbeam
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID364", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID365", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam

///Calc false EPG
////RDK2MCAnalysis* theAnal=new RDK2MCAnalysis( "AID20",  "Corrected Length Contraction", 48,  "RID311", "RID312",  48, "RID313",  "RID314", "RID324",  STD_EPCutSetB);
//double error;
//double answer=theAnal->CalcEPPer3Decay(error);
//cout << "EP Rate: " << answer << " +/- " << error << endl;
//answer=theAnal->CalcEPFalseGPerEP(DET_EPG,error);
//cout << "False EPG: " << answer << " +/- " << error << endl;
//delete theAnal;

/////Make Event for Deconvolution
//    cEvents theEvents;
//    int numEventsStart=200;
//    int numEventFiles=200;
//    int numEventsPerFile=250000;
//    int numEvents=numEventFiles*numEventsPerFile;
//    TString filePath;
//
//
//    for (int i = numEventsStart;i< numEventFiles;i++)
//    {
//        double gLow=(i+1);
//        double gHigh=(i+2);
//        theEvents.setGoodNormConstant(FOURBODY,0,gLow);
//        filePath=EventsDir+"Events_Decov_4B_140415_Set"+int2str(i)+".root";
//    	theEvents.makeEventFile(filePath,MRKRunDir+"Fluxmap_RDK2_Matt1.txt",numEventsPerFile,FOURBODY,0,false,gLow,gHigh,true,-.12,.32);
//    }

///Make events on cluster (move old copies of 130429 to old folder)
//    makeMacrosForEventsOnCluster(0,48,"MAKEEVENTS_3B",THREEBODY,"130429");
//    makeMacrosForEventsOnCluster(0,48,"MAKEEVENTS_4B",FOURBODY,"130429");

///Quick originating point hist
//    RDK2Set threeBodySet("std3","std3",48,"130429","STD_3B","RID278","RID279","");
//    RDK2Set fourBodySet("std4","std4",1,"130429","STD_4B_100eV_Endpoint","RID280","RID281","RID282");
//
//
//    HistDim z0Dim={200,-13,40};
////    TH1* theHist=fourBodySet.makeHist("z0epg","z0",STD_EPGCut,z0Dim);
//    TH1* theHist2=threeBodySet.makeHist("z0ep","z0",STD_EPCut,z0Dim);
//    TFile f(HistsDir+"eventsz0.root","recreate");
//    f.cd();
////    theHist->Write();
//    theHist2->Write();
//    f.Close();

///Create analysis of PLO off
//    createMCAnalysisAndPlot( "AID33",  "Positional Light Output Off", 1,  "RID311", "RID312",  1, "RID313",  "RID314", "RID357",  STD_EPCutSetB);
//    RDK2AnalysisPlotter* thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Effect of PLO","EID3","AID20 AID33"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;

///Create analysis of neutron beam shifts with gamma's run
//    createMCAnalysisAndPlot( "AID29",  "1 mm shift z+ neutronbeam", 1,  "RID341", "RID342",  1, "RID343",  "RID344", "RID362",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID30",  "1 mm shift z- neutronbeam", 1,  "RID345", "RID346",  1, "RID347",  "RID348", "RID363",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID31",  "1 mm shift x+ neutronbeam", 1,  "RID349", "RID350",  1, "RID351",  "RID352", "RID364",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID32",  "1 mm shift x- neutronbeam", 1,  "RID353", "RID354",  1, "RID355",  "RID356", "RID365",  STD_EPCutSetB);
//    RDK2AnalysisPlotter* thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Z shifts of Neutron beam","EID3","AID20 AID29 AID30"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Z shifts of Neutron beam","EID3","AID20 AID31 AID32"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;

///Convert BAPD Field 140422
//    TString txtPath="/media/mjbexternal/comsol/BAPDFieldNormal_140418.txt";
//    TString rootPath="/home/mjbales/school/rdk/comsol/COMSOL_BAPDOnlyField_A_1350V_140418.root";
//    ROOTField newField(txtPath,"bapd",1,3,3);
//    newField.saveFieldToFile(rootPath,"bapd");

///Make BAPD macro
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID369", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD"); //New BAPD Field attempt #1329384728937

///Test Deconv functions
//    const double BGODECONV_VARHIST_LOWEDGES[22]={0,10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,800};
//    TH1D* varHistPrototype=new TH1D("varHistPrototype","varHistPrototype",21,BGODECONV_VARHIST_LOWEDGES);
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"VarSpecialHist_Exp_EID3_EPG_EDepGAvg.txt",varHistPrototype);
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"VarSpecialHist_MC_AID20_EPG_EDepGAvg.txt",varHistPrototype);
////    TH1* expHists[1]={expHist};
////    TH1* mcHists[1]={mcHist};
////    plotExpVersusMCToImage(1,expHists,1,mcHists,"VarHistsComparison","VarHistsComparison140427.png");
//    RDK2Set theSet("DeconvSet","DeconvSet",48,"130429","STD_4B_100eV_Endpoint","RID313","RID314","RID324");
//
//    //Make theoryHist
////    TH1D* theoryHist=(TH1D*) theSet.makeHist("theoryHist","eg0","",varHistPrototype);
////    theoryHist->Scale(0.010277/(double) theSet.GetNumEvents());
////    convertTH1ToTXT(theoryHist,HistsDir+"EG0_140427.txt");
//    TH1D* theoryHist=getTabSeperatedHist(HistsDir+"EG0_140427.txt",varHistPrototype);
//
//    double epRate=0.13993;
////    make2DHistForDeconv("DC1", &theSet,epRate,  varHistPrototype);
//    performDeconvFit("DC1", 0,epRate,expHist);
//    makeComparisonPlotFromFit("DC1", 0,theoryHist);
//
//    delete expHist;
//    delete varHistPrototype;
//    delete theoryHist;
//    delete mcHist;

///Analyze some runs
//    createMCAnalysisAndPlot( "AID29",  "1 mm shift z+ neutronbeam", 1,  "RID341", "RID342",  1, "RID343",  "RID344", "RID362",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID30",  "1 mm shift z- neutronbeam", 1,  "RID345", "RID346",  1, "RID347",  "RID348", "RID363",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID31",  "1 mm shift x+ neutronbeam", 1,  "RID349", "RID350",  1, "RID351",  "RID352", "RID364",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID32",  "1 mm shift x- neutronbeam", 1,  "RID353", "RID354",  1, "RID355",  "RID356", "RID365",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID34",  "20 nm Si Dead Layer on SBD", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID35",  "New BAPD Field Change (3 body only)", 1,  "RID369", "RID359",  1, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
//    RDK2AnalysisPlotter* thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Z shifts of Neutron beam","EID3","AID20 AID29 AID30"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("X shifts of Neutron beam","EID3","AID20 AID31 AID32"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("20nm Si Dead Layer","EID3","AID20 AID34"); thePlotter->makeAll(); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("New BAPD Field","EID3","AID34 AID35"); thePlotter->makeEPPlot(DET_EP); delete thePlotter;

///Make More Event for Deconvolution
//    cEvents theEvents;
//    int numEventsStart=600;
//    int numEventFiles=100;
//    int numEventsPerFile=250000;
//    int numEvents=numEventFiles*numEventsPerFile;
//    TString filePath;
//
//
//    for (int i = numEventsStart;i< numEventsStart+numEventFiles;i++)
//    {
//        double gLow=(i+1);
//        double gHigh=(i+2);
//        theEvents.setGoodNormConstant(FOURBODY,0,gLow);
//        filePath=EventsDir+"Events_Deconv_4B_140415_Set"+int2str(i)+".root";
//    	theEvents.makeEventFile(filePath,MRKRunDir+"Fluxmap_RDK2_Matt1.txt",numEventsPerFile,FOURBODY,0,false,gLow,gHigh,true,-.12,.32);
//    }

///Create 20 nm si Dead layer run with energy resolution
//    createMCAnalysisAndPlot( "AID36",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
//    RDK2AnalysisPlotter* thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("20nm Si Dead Layer","EID3","AID20 AID34 AID36"); thePlotter->makeAll(); delete thePlotter;

///Make Fields along center line
//    TString sbdFieldPath=MRKRunDir+"COMSOL_SBDElectricField_Fine130515.root";
//    TString mirrorFieldPath=MRKRunDir+"COMSOL_MirrorElectricField_Fine130515.root";
//    TString bapdFieldPath=MRKRunDir+"COMSOL_BAPDOnlyField_A_1350V_140418.root";
//    TString solenoidFieldPath=MRKRunDir+"Bfield_RDK2_Matt2.root";
//
//    ROOTField sbdField(sbdFieldPath, "sbd", -25000,3,3);
//    ROOTField mirrorField(mirrorFieldPath, "mirror", 1400,3,3);
//    ROOTField bapdField(bapdFieldPath, "bapd", 1,3,3); bapdField.setSymmetry(false,true,false);
//    ROOTField solenoidField(solenoidFieldPath, "solenoid", 1,3,3); solenoidField.setSymmetry(false,true,false);
//
//    const int numberOfPoints=10000;
//    const double plotDistance=.705; //in meters
//    double distanceBetweenPoints=plotDistance/((double) numberOfPoints);
//    double x[numberOfPoints];
//    double eField[numberOfPoints];
//    double bField[numberOfPoints];
//    double potentials[numberOfPoints];
////    TH1D eFieldAlongCenter("eFieldAlongCenter","Electric Field Along Centerline;m;V/m",numberOfPoints,-.5*distanceBetweenPoints,plotDistance-.5*distanceBetweenPoints);
////    TH1D bFieldAlongCenter("bFieldAlongCenter","Magnetic Field Along Centerline;m;T",numberOfPoints,-.5*distanceBetweenPoints,plotDistance-.5*distanceBetweenPoints);
////    TH1D potentialAlongCenter("potentialAlongCenter","Potential Along Centerline;m;V",numberOfPoints,-distanceBetweenPoints,plotDistance-distanceBetweenPoints);
//
////    TVector3 position(0,0,solenoidField.getZFinal()-.001);
//    TVector3 position(0,0,.424);
//    double potential=0;
//    for (int i = 0;i< numberOfPoints;i++)
//    {
//        x[i]=.424 - i*distanceBetweenPoints;
//        TVector3 eFieldVec(0,0,0);
//        TVector3 bFieldVec(0,0,0);
//
//
//        sbdField.linearInterp3D(position,eFieldVec);
//        eFieldVec*=-1.;
//        mirrorField.linearInterp3D(position,eFieldVec);
//        bapdField.linearInterp3D(position,eFieldVec);
//        solenoidField.linearInterp3D(position,bFieldVec);
//
//        TVector3 dL;
//
//        if(bFieldVec.Mag() > 0)
//            dL =bFieldVec.Unit()*distanceBetweenPoints; //Magnetic field direction times the distance step
//        else
//            dL = TVector3(0,0,-distanceBetweenPoints);
//
//        potential+=-eFieldVec.Dot(dL);
//
//        eField[i]=eFieldVec.Mag();
//        bField[i]=bFieldVec.Mag();
//        potentials[i]=potential;
//
//
//        position+=dL;
//    }
//
////    TString eFieldHistPath=HistsDir+"eFieldAlongCenter140508.txt";
////    TString bFieldHistPath=HistsDir+"bFieldAlongCenter140508.txt";
////    TString potentialHistPath=HistsDir+"potentialAlongCenter140508.txt";
////
////    convertTH1ToTXT(&eFieldAlongCenter,eFieldHistPath);
////    convertTH1ToTXT(&bFieldAlongCenter,bFieldHistPath);
////    convertTH1ToTXT(&potentialAlongCenter,potentialHistPath);
////
////    TH1* theHists[3]={&eFieldAlongCenter,&bFieldAlongCenter,&potentialAlongCenter};
////    scaleHistogramsTogether(&eFieldAlongCenter,&bFieldAlongCenter);
////    scaleHistogramsTogether(&eFieldAlongCenter,&potentialAlongCenter);
////
////    plotExpVersusMCToImage(0,NULL,3,theHists,"Fields","test.png");
//
//    TGraph eFieldGraph(numberOfPoints,x,eField);
//    TGraph bFieldGraph(numberOfPoints,x,bField); bFieldGraph.SetLineColor(kRed); bFieldGraph.SetLineWidth(2);
//    TGraph potentialGraph(numberOfPoints,x,potentials); potentialGraph.SetLineColor(kBlue); potentialGraph.SetLineWidth(2);
//
////    TCanvas theCanvas("plotManyHistCanvas","plotManyHistCanvas",10,10,1200,700);
////    theCanvas.SetLeftMargin(.09);
////    theCanvas.SetRightMargin(.05);
////    theCanvas.SetTopMargin(.05);
////    gPad->SetTickx(1);
////    gPad->SetTicky(1);
////    gPad->SetFillColor(kWhite);
////    TCanvas* theCanvas=getThesisPlotCanvas();
//    TCanvas* theCanvas=getHalfPresentationPlot2();
//
////
////    blankHist.SetTitle("");
//
//    eFieldGraph.SetTitle(";Distance Along Centerline (m);Electric Field Strength (V/m)");
//    bFieldGraph.SetTitle(";Distance Along Centerline (m);Magnetic Field Strength (T)");
//    potentialGraph.SetTitle(";Distance Along Centerline (m);Electric Potential (V)");
//
//    eFieldGraph.GetXaxis()->SetRangeUser(-2.8,4.2);
////    eFieldGraph.Draw("AL"); theCanvas.SaveAs(GraphsDir+"eField140508.pdf");
////    eFieldGraph.Draw("AL"); theCanvas.SaveAs(GraphsDir+"eField140508.pdf");
//
//    bFieldGraph.GetYaxis()->SetRangeUser(3.2,4.7);
//    bFieldGraph.GetXaxis()->SetRangeUser(-.28,.42);
////    bFieldGraph.Draw("AL"); theCanvas->SaveAs(GraphsDir+"bField140508.pdf");
//    bFieldGraph.Draw("AL"); theCanvas->SaveAs(GraphsDir+"bField140508b.png");
//    theCanvas->SetLogy(); potentialGraph.Draw("AL");
//
//
//
//    TPaveText textBox1(0.75,0.6706767,0.85,0.7398496,"brNDC");
//    textBox1.SetTextSize(.06);
//    textBox1.AddText("Mirror Field");
//    textBox1.SetFillColor(kWhite);
//    textBox1.SetLineColor(kWhite);
//    textBox1.SetShadowColor(kWhite);
//    textBox1.SetLineWidth(0);
//    textBox1.Draw("same");
//
//    TPaveText textBox2(0.45,0.2285714,0.5,0.275188,"brNDC");
//    textBox2.SetTextSize(.06);
//    textBox2.AddText("Bare APD Field");
//    textBox2.SetFillColor(kWhite);
//    textBox2.SetLineColor(kWhite);
//    textBox2.SetShadowColor(kWhite);
//    textBox2.SetLineWidth(0);
//    textBox2.Draw("same");
//
//    TPaveText textBox3(0.1103679,0.3,0.25,0.4766917,"brNDC");
//    textBox3.SetTextSize(.06);
//    textBox3.AddText("SBD Field");
//    textBox3.AddText("(negative)");
//    textBox3.SetFillColor(kWhite);
//    textBox3.SetLineColor(kWhite);
//    textBox3.SetShadowColor(kWhite);
//    textBox3.SetLineWidth(0);
//    textBox3.Draw("same");
//
//    potentialGraph.GetXaxis()->SetRangeUser(-.28,.42);
//    potentialGraph.GetYaxis()->SetRangeUser(1.,30000.);
////    theCanvas->SaveAs(GraphsDir+"potential140508.pdf");
//    theCanvas->SaveAs(GraphsDir+"potential140508b.png");

///Plot waveforms for thesis

//    TCanvas* theCanvas=getThesisPlotCanvas();
//    theCanvas->SetLeftMargin(.12);
//    TString waveformList="EP1";
////    TString waveformList="BGO1 BGO2 BGO3 BGO4 BGO5";
////    TString waveformList="BAPD1 BAPD2 BAPD3 BAPD4 BAPD5";
//    for (int i = 0;i< numItemsInStringList(waveformList);i++)
//    {
//        TString waveformName;
//        stringFromList(waveformList,i,waveformName);
//
//        TString waveformPath=HistsDir+"WaveForm"+waveformName+".txt";
//        HistDim waveformDim={2048,0.5,2048.5};
//        TH1* waveformHist=getTabSeperatedHist(waveformPath,waveformDim);
//        waveformHist->SetTitle(";Channel Number (1 ch = 40 ns);Signal Channel");
//
//        waveformHist->GetYaxis()->SetTitleOffset(1.3);
//        waveformHist->SetLineColor(kBlue);
//        waveformHist->Draw("L");
//
//
//        TPaveText textBox2(0.2700669,0.7789474,0.3578261,0.8225564,"brNDC");
//        textBox2.SetTextSize(.03);
//        textBox2.AddText("Electron");
//        textBox2.SetFillColor(kWhite);
//        textBox2.SetLineColor(kWhite);
//        textBox2.SetShadowColor(kWhite);
//        textBox2.SetLineWidth(0);
//        textBox2.Draw("same");
//
//        TPaveText textBox3(0.4,0.3,0.48,0.37,"brNDC");
//        textBox3.SetTextSize(.03);
//        textBox3.AddText("Proton");
//        textBox3.SetFillColor(kWhite);
//        textBox3.SetLineColor(kWhite);
//        textBox3.SetShadowColor(kWhite);
//        textBox3.SetLineWidth(0);
//        textBox3.Draw("same");
//
//        theCanvas->SaveAs(GraphsDir+"WaveForm"+waveformName+".pdf");
//        delete waveformHist;
//    }

///Make money graph
//    HistDim z0Dim={200,-.18, .62};
//
////    RDK2MCAnalysis theAnalysis("AID36",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
////    TH1* z0EPHist=theAnalysis.MakeStandardHist("epHist",DET_EP,"z0",z0Dim);
////    TH1* z0EPGHist=theAnalysis.MakeStandardHist("epgHist",DET_EPG,"z0",z0Dim);
////    TH1* z0EPBGHist=theAnalysis.MakeStandardHist("epbgHist",DET_EPBG,"z0",z0Dim);
////
////    convertTH1ToTXT(z0EPHist,HistsDir+"z0EPHist140511.txt");
////    convertTH1ToTXT(z0EPGHist,HistsDir+"z0EPGHist140511.txt");
////    convertTH1ToTXT(z0EPBGHist,HistsDir+"z0EPBGHist140511.txt");
//
//    TCanvas* theCanvas=new TCanvas("thesisPlot","thesisPlot",10,10,2400,700);
//    theCanvas->SetLeftMargin(.06);
//    theCanvas->SetRightMargin(.01);
//    theCanvas->SetTopMargin(.05);
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//    TH1* z0EPHist=getTabSeperatedHist(HistsDir+"z0EPHist140511.txt",z0Dim);
//    TH1* z0EPGHist=getTabSeperatedHist(HistsDir+"z0EPGHist140511.txt",z0Dim);
//    TH1* z0EPBGHist=getTabSeperatedHist(HistsDir+"z0EPBGHist140511.txt",z0Dim);
//
//    for (int i = 1;i<=z0Dim.numBins ;i++)
//    {
//        z0EPHist->SetBinError(i,0);
//        z0EPGHist->SetBinError(i,0);
//        z0EPBGHist->SetBinError(i,0);
//
//    }
//
//    z0EPHist->Rebin(2);
//    z0EPGHist->Rebin(2);
//    z0EPBGHist->Rebin(2);
//
//    scaleHistTo(z0EPHist,1.);
//    scaleHistTo(z0EPGHist,.3);
//    scaleHistTo(z0EPBGHist,.20);
//
//
//    z0EPHist->SetTitle(";Decay Position in Beam (m);Relative Detection Probability (AU)");
//
//    z0EPHist->SetLineColor(kBlue);
//    z0EPGHist->SetLineColor(kRed);
//    z0EPBGHist->SetLineColor(kGreen);
//
//    z0EPHist->GetXaxis()->SetRangeUser(-.3,.35);
//    z0EPHist->GetYaxis()->SetRangeUser(0,.025);
//
//    z0EPHist->GetYaxis()->SetTitleOffset(.6);
//
//
//    z0EPHist->Draw("L");
//    z0EPGHist->Draw("L same");
//    z0EPBGHist->Draw("L same");
//
//    TLegend legend(0.1,0.7,0.25,0.9,NULL,"brNDC");
//    legend.SetFillColor(kWhite);
//    legend.AddEntry(z0EPHist,"EP Coincidence");
//    legend.AddEntry(z0EPGHist,"EPG Coincidence (BGO/APD)");
//    legend.AddEntry(z0EPBGHist,"EPG Coincidence (APD)");
//    legend.Draw();
//
//    theCanvas->SaveAs(GraphsDir+"z0Plot140518.pdf");
//
//    delete theCanvas;

///Make various exp graphs
///EG Delta T
//    HistDim egDeltaTDim={900,-450.5,449.5};
//    RDK2ExpAnalysis theAnalysis( "EID4",  "New PE cal From Jeff", "PID9", "G7_Cycle5_8", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//
//    TH1* bgoHist=makeTH1DFromDim("bgoHist",";E-G Time Difference (s);Counts",egDeltaTDim);
//    TH1* bapdHist=makeTH1DFromDim("bapdHist",";E-G Time Difference (s);Counts",egDeltaTDim);
//
//    //BGO
//    for (int i = 0;i< 12;i++)
//    {
//    	TH1* tempHist=theAnalysis.MakeHist("tempHist","BGODeltaT"+int2str(i+1),STD_EPCutSetB.CreateExpEPGCut(0,0,0,i+1,DET_EPG),egDeltaTDim);
//    	bgoHist->Add(tempHist);
//    	delete tempHist;
//    }
//
//    convertTH1ToTXT(bgoHist,HistsDir+"BGOEGDeltaT_EID4_140511.txt");
//    convertTH1ToTXT(bapdHist,HistsDir+"BAPDEGDeltaT_EID4_140511.txt");
//
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
//    theCanvas->SetLeftMargin(.12);
//    bgoHist->SetLineColor(kBlue);
//    bgoHist->Draw("p");
//    theCanvas->SaveAs(GraphsDir+"egDeltaTPlot140511.pdf");

//    delete theCanvas;

//    RDK2AnalysisPlotter* thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("SBD shift x' direction","","AID33 AID23 AID24"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("SBD shift z' direction","","AID33 AID21 AID22"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("B Field shift z' direction","","AID33 AID25 AID26"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("B Field shift x' direction","","AID33 AID27 AID28"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Neutron Beam shift z' direction","","AID33 AID29 AID30"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Neutron Beam shift x' direction","","AID33 AID31 AID32"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
////    thePlotter=new RDK2AnalysisPlotter("Neutron Beam shift x' direction","","AID20 AID33 AID34 AID35 AID36"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    delete thePlotter;

///Create eventfile gamma cap for deconvolution events
//    cEvents theEvents;
//    int numEventsStart=600;
//    int numEventFiles=100;
//    int numEventsPerFile=250000;
//    double gLow=400;
//    double gHigh=EKEMAX;
//    theEvents.setGoodNormConstant(FOURBODY,0,gLow);
//    TString filePath=MRKRunDir+"Events_Deconv_4B_140415_Set400CAP.root";
//    theEvents.makeEventFile(filePath,MRKRunDir+"Fluxmap_RDK2_Matt1.txt",numEventsPerFile,FOURBODY,0,false,gLow,gHigh,true,-.12,.32);

///Gotta keep that cluster going....lets get some systematic statistics!
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID370", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID371", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID372", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID373", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam

///More Runs!!!! 140516
///No length contraction,bapd length contraction, and gamma det+mirror offsets
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID374", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/det/setBGOLengthContractionFactor 1.0\n/MRK/det/setAlLengthContractionFactor 1.0"); //No Length Contraction
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID375", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //Normal with BAPD length contraction
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID376", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm"); //1 mm shift z+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID377", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm"); //1 mm shift z+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID378", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm"); //1 mm shift z+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID379", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm"); //1 mm shift z+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID380", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm"); //1 mm shift z+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID381", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm"); //1 mm shift z- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID382", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm"); //1 mm shift z- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID383", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm"); //1 mm shift z- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID384", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm"); //1 mm shift z- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID385", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm"); //1 mm shift z- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID386", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm"); //1 mm shift x+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID387", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm"); //1 mm shift x+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID388", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm"); //1 mm shift x+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID389", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm"); //1 mm shift x+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID390", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm"); //1 mm shift x+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID391", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm"); //1 mm shift x- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID392", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm"); //1 mm shift x- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID393", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm"); //1 mm shift x- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID394", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm"); //1 mm shift x- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID395", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm"); //1 mm shift x- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID396", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm"); //1 mm shift y+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID397", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm"); //1 mm shift y+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID398", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm"); //1 mm shift y+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID399", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm"); //1 mm shift y+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID400", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm"); //1 mm shift y+ gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID401", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm"); //1 mm shift y- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID402", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm"); //1 mm shift y- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID403", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm"); //1 mm shift y- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID404", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm"); //1 mm shift y- gamma det+electrostatic mirror
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID405", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm"); //1 mm shift y- gamma det+electrostatic mirror

///RDK I correction calculations 140516
//    double calEnergy=59.5;
//    double calNonLinearFactor=evalMattAllExpParamModel(60,ALLEXP_PARAM_MODEL);
//    cout << calEnergy << " keV has a nonlinear factor of " << calNonLinearFactor << endl;
//    double lowLimit=15;
//    double highLimit=340;
//    double lowLimitGuess=19;
//    double highLimitGuess=313;
//    cout << "Lower limit was 15 and now is "<< lowLimitGuess*evalMattAllExpParamModel(lowLimitGuess,ALLEXP_PARAM_MODEL)/calNonLinearFactor << endl;
//    cout << "Upper limit was 340 and now is "<< highLimitGuess*evalMattAllExpParamModel(highLimitGuess,ALLEXP_PARAM_MODEL)/calNonLinearFactor << endl;
//
//    cEvents theEvents;
//
//    double noFermiBROld = theEvents.calcBR(lowLimit,highLimit,0,10000000,false);
//    double noFermiBROldNewLife = theEvents.calcBR(lowLimit,highLimit,0,10000000,false,880.);
//    double noFermiBR = theEvents.calcBR(lowLimitGuess,highLimitGuess,0,10000000,false);
//    double fermiBR = theEvents.calcBR(lowLimitGuess,highLimitGuess,0,10000000,true);
//    double fermiBRNewLife = theEvents.calcBR(lowLimitGuess,highLimitGuess,0,10000000,true,880.);
//
//    cout << "No Fermi function rate with old range was: " << noFermiBROld << endl;
//    cout << "No Fermi function rate with old range but new lifetime is: " << noFermiBROldNewLife << endl;
//    cout << "No Fermi function rate with new range is: " << noFermiBR << endl;
//    cout << "With Fermi function rate with new range is: " << fermiBR << endl;
//    cout << "With Fermi function rate with new range and new lifetime is: " << fermiBRNewLife << endl;

///RDKII paper range

//    double lowLimitGuess=15;
//    double highLimitGuess=781;
//
//    cEvents theEvents;
//    double fermiBRNewLife = theEvents.calcBR(lowLimitGuess,highLimitGuess,0,10000000,true,880.);
//    cout << "With Fermi function rate with new range and new lifetime is: " << fermiBRNewLife << endl;

///New z0 hist for paper
//HistDim z0Dim={280,-.4, .72};
//
////    RDK2MCAnalysis theAnalysis("AID36",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
////    TH1* z0EPHist=theAnalysis.MakeStandardHist("epHist",DET_EP,"z0",z0Dim);
////    TH1* z0EPGHist=theAnalysis.MakeStandardHist("epgHist",DET_EPG,"z0",z0Dim);
////    TH1* z0EPBGHist=theAnalysis.MakeStandardHist("epbgHist",DET_EPBG,"z0",z0Dim);
////
////    convertTH1ToTXT(z0EPHist,HistsDir+"z0EPHist140511.txt");
////    convertTH1ToTXT(z0EPGHist,HistsDir+"z0EPGHist140511.txt");
////    convertTH1ToTXT(z0EPBGHist,HistsDir+"z0EPBGHist140511.txt");
//
//    TCanvas* theCanvas=new TCanvas("thesisPlot","thesisPlot",10,10,2400,1000);
//    theCanvas->SetLeftMargin(.085);
//    theCanvas->SetRightMargin(.01);
//    theCanvas->SetTopMargin(.05);
//    theCanvas->SetBottomMargin(.18);
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//    TH1* z0EPHist=getTabSeperatedHist(HistsDir+"z0EPHist140518.txt",z0Dim);
//    TH1* z0EPGHist=getTabSeperatedHist(HistsDir+"z0EPGHist140518.txt",z0Dim);
//    TH1* z0EPBGHist=getTabSeperatedHist(HistsDir+"z0EPBGHist140518.txt",z0Dim);
//
//    for (int i = 1;i<=z0Dim.numBins ;i++)
//    {
//        z0EPHist->SetBinError(i,0);
//        z0EPGHist->SetBinError(i,0);
//        z0EPBGHist->SetBinError(i,0);
//
//    }
//
//    z0EPHist->Rebin(2);
//    z0EPGHist->Rebin(2);
//    z0EPBGHist->Rebin(2);
//
//    scaleHistTo(z0EPHist,1.*50);
//    scaleHistTo(z0EPGHist,.3*50);
//    scaleHistTo(z0EPBGHist,.15*50);
//
//
//    z0EPHist->SetTitle(";Decay Position in Neutron Beam (m);Detection Prob.(AU)");
////    z0EPHist->SetTitle(";;");
//
//    z0EPHist->SetLineColor(kBlue);
//    z0EPGHist->SetLineColor(kRed);
//    z0EPBGHist->SetLineColor(kGreen);
//
//    z0EPHist->SetLineWidth(2);
//    z0EPGHist->SetLineWidth(2);
//    z0EPBGHist->SetLineWidth(2);
//
//    z0EPHist->GetXaxis()->SetRangeUser(-.39,.45);
//    z0EPHist->GetXaxis()->SetLabelSize(.08);
//    z0EPHist->GetXaxis()->SetTitleSize(.08);
//
//    z0EPHist->GetYaxis()->SetRangeUser(0,.020*50);
//    z0EPHist->GetYaxis()->SetTitleOffset(.5);
//    z0EPHist->GetYaxis()->SetLabelSize(.08);
//    z0EPHist->GetYaxis()->SetTitleSize(.08);
//    z0EPHist->GetYaxis()->SetNdivisions(508);
//
//
//
//    z0EPHist->Draw("L");
//    z0EPGHist->Draw("L same");
//    z0EPBGHist->Draw("L same");
//    z0EPHist->Draw("L same");
//
//    TLegend legend(0.11,0.4,0.42,0.9,NULL,"brNDC");
////    legend.SetTextSize(1.4);
//    legend.SetFillColor(kWhite);
//    legend.AddEntry(z0EPHist,"EP");
//    legend.AddEntry(z0EPGHist,"EPG (BGO/APD)");
//    legend.AddEntry(z0EPBGHist,"EPG (Direct APD)");
//    legend.Draw();
//
//    theCanvas->SaveAs(GraphsDir+"z0Plot140518.pdf");
//
//    delete theCanvas;

///BGO true energy range
//double calNonLinearFactor=evalMattAllExpParamModel(511,ALLEXP_PARAM_MODEL);
//double lowLimitGuess=14.06;
//cout << "Lower limit was 10 and now is "<< lowLimitGuess*evalMattAllExpParamModel(lowLimitGuess,ALLEXP_PARAM_MODEL)/calNonLinearFactor << endl;
//cEvents theEvents;
//    double fermiBRNewLife = theEvents.calcBR(lowLimitGuess,EKEMAX,0,10000000,true,880.);
//    cout << fermiBRNewLife << endl;

///Deconv 140525
///DC2
//    double binBoundaries[402];
//    for (int i = 0;i< 401;i++)
//    {
//        binBoundaries[i]=i;
//    }
//    binBoundaries[401]=800;
//
//    TH1D* varHistPrototype=new TH1D("varHistPrototype","varHistPrototype",401,binBoundaries);
////    TH1D* expHist=getTabSeperatedHist(HistsDir+"VarSpecialHist_Exp_EID3_EPG_EDepGAvg.txt",varHistPrototype);
////    TH1D* mcHist=getTabSeperatedHist(HistsDir+"VarSpecialHist_MC_AID20_EPG_EDepGAvg.txt",varHistPrototype);
////    TH1* expHists[1]={expHist};
////    TH1* mcHists[1]={mcHist};
////    plotExpVersusMCToImage(1,expHists,1,mcHists,"VarHistsComparison","VarHistsComparison140427.png");
//    RDK2Set theSet("DeconvSet","DeconvSet",400,"140415","Deconv_4B","RID366","RID367","RID368");
//
//    //Make theoryHist
//    TH1D* theoryHist=(TH1D*) theSet.makeHist("theoryHist","eg0","",varHistPrototype);
//    scaleHistTo(theoryHist,0.010277);
//    convertTH1ToTXT(theoryHist,HistsDir+"EG0_140525.txt");
////    TH1D* theoryHist=getTabSeperatedHist(HistsDir+"EG0_140525.txt",varHistPrototype);
//
//    make2DHistForDeconv("DC2", &theSet,1,  varHistPrototype);
////    performDeconvFit("DC1", 0,1,expHist);
////    makeComparisonPlotFromFit("DC1", 0,theoryHist);
//
////    delete expHist;
//    delete varHistPrototype;
//    delete theoryHist;
////    delete mcHist;
///DC3 - 1-92 keV with Deconv
//    const double binBoundaries[22]={0,10,12,13,14,15,17,19,21,24,27,30,34,38,43,49,56,63,71,80,92,800};
//    TH1D* varHistPrototype=new TH1D("varHistPrototype","varHistPrototype",21,binBoundaries);
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"VarSpecialHist_Exp_EID3_EPG_EDepGAvg.txt",varHistPrototype);
////    TH1D* mcHist=getTabSeperatedHist(HistsDir+"VarSpecialHist_MC_AID20_EPG_EDepGAvg.txt",varHistPrototype);
////    TH1* expHists[1]={expHist};
////    TH1* mcHists[1]={mcHist};
////    plotExpVersusMCToImage(1,expHists,1,mcHists,"VarHistsComparison","VarHistsComparison140427.png");
//    RDK2Set theSet("DeconvSet","DeconvSet",90,"140415","Deconv_4B","RID366","RID367","RID368");
//
//    TH1D* theoryHist=getTabSeperatedHist(HistsDir+"EG0_140427.txt",varHistPrototype);
//
//    make2DHistForDeconv("DC3", &theSet,1,  varHistPrototype);
////    performDeconvFit("DC1", 0,1,expHist);
////    makeComparisonPlotFromFit("DC1", 0,theoryHist);
//
////    delete expHist;
//    delete varHistPrototype;
//    delete theoryHist;
////    delete mcHist;

///SBD Corrected position based on length contraction
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID406", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID407", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID408", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)
//makeGeantMacroFilesForClusterEventGenerator(0,48, "RID409", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)
//makeGeantMacroFilesForClusterEventGenerator(0,10, "RID410", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -3 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift from bend
//makeGeantMacroFilesForClusterEventGenerator(0,10, "RID411", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -3 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift from bend
//makeGeantMacroFilesForClusterEventGenerator(0,10, "RID412", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -3 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift from bend
//makeGeantMacroFilesForClusterEventGenerator(0,10, "RID413", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -3 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift from bend

///Analyze some runs!
//    createMCAnalysisAndPlot( "AID37",  "Standard (no SBD Length Contr)", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID375",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID38",  "No Length Contraction", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID374",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID39",  "1 mm shift x+ neutronbeam (EP only)", 48,  "RID370", "RID371",  1, "RID360",  "RID361", "RID375",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID40",  "1 mm shift x- neutronbeam (EP only)", 48,  "RID372", "RID373",  1, "RID360",  "RID361", "RID375",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID41",  "1 mm shift z+ gamma det", 10,  "RID376", "RID377",  10, "RID378",  "RID379", "RID380",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID42",  "1 mm shift z- gamma det", 10,  "RID381", "RID382",  10, "RID383",  "RID384", "RID385",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID43",  "1 mm shift x+ gamma det", 10,  "RID386", "RID387",  10, "RID388",  "RID389", "RID390",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID44",  "1 mm shift x- gamma det", 10,  "RID391", "RID392",  10, "RID393",  "RID394", "RID395",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID45",  "1 mm shift y+ gamma det", 10,  "RID396", "RID397",  10, "RID398",  "RID399", "RID400",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID46",  "1 mm shift y- gamma det", 10,  "RID401", "RID402",  10, "RID403",  "RID404", "RID405",  STD_EPCutSetB);
//
//    RDK2AnalysisPlotter* thePlotter;
////    thePlotter=new RDK2AnalysisPlotter("Gamma Det Length Contraction","EID4","AID37 AID38"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Neutron Beam Shift X Higher Stats (Three body Only)","EID4","AID37 AID39 AID40"); thePlotter->makeEPPlot(DET_EP); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Gamma Det+MirrorShift Z","EID4","AID37 AID41 AID42"); thePlotter->makeAll(); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Gamma Det+MirrorShift X","EID4","AID37 AID43 AID44"); thePlotter->makeAll(); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Gamma Det+MirrorShift Y","EID4","AID37 AID45 AID46"); thePlotter->makeAll(); delete thePlotter;

///Make crap ton of runs
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID414", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID415", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID416", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset 1. 0. 0. mm"); //1 mm shift x+ neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID417", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID418", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID419", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/setBeamOffset -1. 0. 0. mm"); //1 mm shift x- neutronbeam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID420", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID421", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID422", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID423", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID424", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID425", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID426", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID427", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID428", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID429", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID430", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID431", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID432", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID433", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID434", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID435", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID436", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID437", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID438", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID439", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID440", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID441", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID442", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID443", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID444", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID445", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID446", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID447", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID448", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID449", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID450", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)+BAPD Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID451", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)+BAPD Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID452", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)+BAPD Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID453", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)+BAPD Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID454", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //Length contraction corrected SBD (2 mm from bend)+BAPD Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID455", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.6726951503 0 -2.03899024875 mm"); //Length contraction corrected SBD (2 mm from bend)+1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID456", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.6726951503 0 -2.03899024875 mm"); //Length contraction corrected SBD (2 mm from bend)+1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID457", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.6726951503 0 -2.03899024875 mm"); //Length contraction corrected SBD (2 mm from bend)+1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID458", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.6726951503 0 -2.03899024875 mm"); //Length contraction corrected SBD (2 mm from bend)+1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID459", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.6726951503 0 -2.03899024875 mm"); //Length contraction corrected SBD (2 mm from bend)+1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID460", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset -1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -1.2998760527 0 -1.70889503695 mm"); //Length contraction corrected SBD (2 mm from bend)-1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID461", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset -1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -1.2998760527 0 -1.70889503695 mm"); //Length contraction corrected SBD (2 mm from bend)-1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID462", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset -1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -1.2998760527 0 -1.70889503695 mm"); //Length contraction corrected SBD (2 mm from bend)-1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID463", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset -1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -1.2998760527 0 -1.70889503695 mm"); //Length contraction corrected SBD (2 mm from bend)-1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID464", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset -1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -1.2998760527 0 -1.70889503695 mm"); //Length contraction corrected SBD (2 mm from bend)-1 mm x SBD

///Make thesis plots
///BR Plot
////    TH1D* brHist=calcBRHist(10000,7000,0.1,700,0,true,true);
////    convertTH1ToTXT(brHist,HistsDir+"brHist140608.txt");
//    TH1D* brHist=getTabSeperatedHist(HistsDir+"brHist140608.txt",7000,0.1,700);
//    TCanvas* theCanvas = getThesisPlotCanvas();
//    brHist->Smooth(10);
//    brHist->SetLineColor(kBlue);
//
//    brHist->SetTitle(";Lower Photon Energy Limit (keV);Branching Ratio");
////    brHist->GetYaxis()->SetTitleOffset(1.5);
//    brHist->GetXaxis()->SetTitleOffset(1.3);
//    brHist->GetSumw2()->Set(0);
//    brHist->GetXaxis()->SetRangeUser(0.2,600);
//    theCanvas->SetLogy();
//    theCanvas->SetLogx();
//    brHist->Draw("C");
//    theCanvas->SaveAs(GraphsDir+"brHist140608.pdf");

///Make events hists
///egcostheta
//    //What I choose isn't very important so I grabbed RID's at random
//    int numEventFiles=48;
//    RDK2MCAnalysis theAnalysis("AID36",  "20 nm Si Dead Layer on SBD with P Resolution", numEventFiles,  "RID358", "RID359",  numEventFiles, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
//
//    HistDim cosDim={100,-1,1};
////    TH1* egcosthetaHist=theAnalysis.MakeHist("egcosthetaHist",FOURBODY,"(mxe0 * mxg0 + mye0 * myg0 + mze0 * mzg0)/(sqrt(mxe0 * mxe0 + mye0 * mye0 + mze0 * mze0) + sqrt(mxg0 * mxg0 + myg0 * myg0 + mzg0 * mzg0))","eg0 > 14.1",cosDim);
////    egcosthetaHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(egcosthetaHist,HistsDir+"egcosthetaHist140610.txt");
//    TH1D* egcosthetaHist=getTabSeperatedHist(HistsDir+"egcosthetaHist140610.txt",cosDim);
//    TCanvas* theCanvas = getThesisPlotCanvas();
////    egcosthetaHist->Smooth(10);
//    egcosthetaHist->SetLineColor(kBlue);
//
//    egcosthetaHist->SetTitle(";cos #theta_{e#gamma};Relative Probability (AU)");
//    egcosthetaHist->GetYaxis()->SetTitleOffset(1.2);
//    egcosthetaHist->GetXaxis()->SetTitleOffset(1.3);
//    egcosthetaHist->GetSumw2()->Set(0);
//    egcosthetaHist->GetXaxis()->SetRangeUser(-1,1);
//    egcosthetaHist->Draw("l");
//    theCanvas->SaveAs(GraphsDir+"egcosthetaHist140610.pdf");
///E_p0

///Check v1.035 vs 1.036
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID465", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //Check vs 1.035
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID466", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //Check vs 1.035
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID467", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //Check vs 1.035
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID468", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //Check vs 1.035
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID469", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //Check vs 1.035

///Process systematic runs
	//createMCAnalysisAndPlot( "AID47",  "V1.036 Old SBD Pos", 10,  "RID465", "RID466",  10, "RID467",  "RID468", "RID469",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID48",  "V1.036 Correct SBD Pos", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID49",  "V1.036 Correct SBD Pos +2mm shift from bend", 10,  "RID410", "RID411",  10, "RID412",  "RID413", "RID375",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID50",  "V1.036 Old SBD +1 mm x neutron beam", 10,  "RID370", "RID371",  10, "RID414",  "RID415", "RID416",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID51",  "V1.036 Old SBD -1 mm x neutron beam", 10,  "RID372", "RID373",  10, "RID417",  "RID418", "RID419",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID52",  "V1.036 Correct SBD Pos 1 mm shift z+ gamma det", 10,  "RID420", "RID421",  10, "RID422",  "RID423", "RID424",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID53",  "V1.036 Correct SBD Pos 1 mm shift z- gamma det", 10,  "RID425", "RID426",  10, "RID427",  "RID428", "RID429",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID54",  "V1.036 Correct SBD Pos 1 mm shift x+ gamma det", 10,  "RID430", "RID431",  10, "RID432",  "RID433", "RID434",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID55",  "V1.036 Correct SBD Pos 1 mm shift x- gamma det", 10,  "RID435", "RID436",  10, "RID437",  "RID438", "RID439",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID56",  "V1.036 Correct SBD Pos 1 mm shift y+ gamma det", 10,  "RID440", "RID441",  10, "RID442",  "RID443", "RID444",  STD_EPCutSetB);
	//createMCAnalysisAndPlot( "AID57",  "V1.036 Correct SBD Pos 1 mm shift y- gamma det", 10,  "RID445", "RID446",  10, "RID447",  "RID448", "RID449",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID58",  "V1.036 Correct SBD Pos+BAPD Field", 10,  "RID450", "RID451",  10, "RID452",  "RID453", "RID454",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID59",  "V1.036 Correct SBD Pos+1 mm x SBD", 10,  "RID455", "RID456",  10, "RID457",  "RID458", "RID459",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID60",  "V1.036 Correct SBD Pos-1 mm x SBD", 10,  "RID460", "RID461",  10, "RID462",  "RID463", "RID464",  STD_EPCutSetB);

//    RDK2AnalysisPlotter* thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Version 1.035 vs Version 1.036","EID4","AID37 AID47"); thePlotter->makeAll(); delete thePlotter; //whoops AID37 uses version 1.036 in the photon
//    thePlotter=new RDK2AnalysisPlotter("Version 1.035 vs Version 1.036","EID4","AID36 AID47"); thePlotter->makeAll(); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Corrected SBD Position +/- 2 mm","","AID47 AID48 AID49"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Corrected SBD Position +/- 2 mm","","AID47 AID48 AID49"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm z","","AID48 AID52 AID53"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm z","","AID48 AID52 AID53"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm x","","AID48 AID54 AID55"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//    thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm x","","AID48 AID54 AID55"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;

///Make a bunch of thesis plots 140619

///z0 Hist for thesis
//    HistDim z0Dim={280,-.4, .72};
//
////    RDK2MCAnalysis theAnalysis("AID36",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
////    TH1* z0EPHist=theAnalysis.MakeStandardHist("epHist",DET_EP,"z0",z0Dim);
////    TH1* z0EPGHist=theAnalysis.MakeStandardHist("epgHist",DET_EPG,"z0",z0Dim);
////    TH1* z0EPBGHist=theAnalysis.MakeStandardHist("epbgHist",DET_EPBG,"z0",z0Dim);
////
////    convertTH1ToTXT(z0EPHist,HistsDir+"z0EPHist140511.txt");
////    convertTH1ToTXT(z0EPGHist,HistsDir+"z0EPGHist140511.txt");
////    convertTH1ToTXT(z0EPBGHist,HistsDir+"z0EPBGHist140511.txt");
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
//    TH1* z0EPHist=getTabSeperatedHist(HistsDir+"z0EPHist140518.txt",z0Dim);
//    TH1* z0EPGHist=getTabSeperatedHist(HistsDir+"z0EPGHist140518.txt",z0Dim);
//    TH1* z0EPBGHist=getTabSeperatedHist(HistsDir+"z0EPBGHist140518.txt",z0Dim);
//
//    for (int i = 1;i<=z0Dim.numBins ;i++)
//    {
//        z0EPHist->SetBinError(i,0);
//        z0EPGHist->SetBinError(i,0);
//        z0EPBGHist->SetBinError(i,0);
//
//    }
//
//    z0EPHist->Rebin(2);
//    z0EPGHist->Rebin(2);
//    z0EPBGHist->Rebin(2);
//
//    scaleHistTo(z0EPHist,1.*50);
//    scaleHistTo(z0EPGHist,.3*50);
//    scaleHistTo(z0EPBGHist,.15*50);
//
//
//    z0EPHist->SetTitle(";Decay Position in Neutron Beam (m);Detection Probability (AU)");
////    z0EPHist->SetTitle(";;");
//
//    z0EPHist->SetLineColor(kBlue);
//    z0EPGHist->SetLineColor(kRed);
//    z0EPBGHist->SetLineColor(kGreen+2);
//
//
//    z0EPGHist->SetLineStyle(2);
//    z0EPBGHist->SetLineWidth(2);
//    z0EPBGHist->SetLineStyle(3);
//    z0EPBGHist->SetLineWidth(2);
//
//    z0EPHist->GetXaxis()->SetRangeUser(-.39,.45);
////    z0EPHist->GetXaxis()->SetLabelSize(.08);
////    z0EPHist->GetXaxis()->SetTitleSize(.08);
//
//    z0EPHist->GetYaxis()->SetRangeUser(0,.020*50);
////    z0EPHist->GetYaxis()->SetTitleOffset(.5);
////    z0EPHist->GetYaxis()->SetLabelSize(.08);
////    z0EPHist->GetYaxis()->SetTitleSize(.08);
//    z0EPHist->GetYaxis()->SetNdivisions(508);
//
//
//
//    z0EPHist->Draw("L");
//    z0EPGHist->Draw("L same");
//    z0EPBGHist->Draw("L same");
//    z0EPHist->Draw("L same");
//
//    TLegend legend(0.11,0.6,0.42,0.9,NULL,"brNDC");
////    legend.SetTextSize(1.4);
//    legend.SetFillColor(kWhite);
//    legend.AddEntry(z0EPHist,"EP");
//    legend.AddEntry(z0EPGHist,"EPG (BGO)");
//    legend.AddEntry(z0EPBGHist,"EPG (LAAPD)");
//    legend.Draw();
//
//    theCanvas->SaveAs(GraphsDir+"z0Plot140619.pdf");
//
//    delete theCanvas;

///Beam profile

//cSField2D fluxMap;
//fluxMap.loadField(MRKRunDir+"Fluxmap_RDK2_Matt1.txt",1);
//TH2D flux2DHist("flux2DHist",";X Direction (cm);Y Direction (cm)",150,-1.51,1.51,150,-1.51,1.51);
//
//for (int i = 0;i< 150;i++)
//{
//    for (int j = 0;j< 150;j++)
//    {
//        double flux=fluxMap.getValue(i,j);
//        flux2DHist.SetBinContent(i+1,j+1,flux);
//    }
//}
//
////    double min=flux2DHist.GetMinimum();
//    flux2DHist.SetMaximum(1);
//    flux2DHist.GetYaxis()->SetTitleOffset(1.2);
//
//    TCanvas* theCanvas=new TCanvas("thesisPlot","thesisPlot",10,10,1200,1200);
//    theCanvas->SetLeftMargin(.11);
//    theCanvas->SetRightMargin(.11);
//    theCanvas->SetTopMargin(.11);
//    theCanvas->SetBottomMargin(.11);
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//
////    const int NRGBs = 7, NCont = 50;
////    gStyle->SetNumberContours(NCont);
//////    Double_t stops[NRGBs] = { 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
//////    Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
//////    Double_t green[NRGBs] = { 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
//////    Double_t blue[NRGBs]  = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
////    Double_t stops[NRGBs] = { 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
////    Double_t blue[NRGBs]   = { 1.00, 0.4, 0.40, 0.20, 0.20, 0.20, 0.00 };
////    Double_t green[NRGBs] = { 1.00, 1.0, 0.40, 0.40, 0.40, 0.20, 0.00 };
////    Double_t red[NRGBs]  = { 1.00, 1.0, 1.00, 1.00, 0.40, 0.40, 0.00 };
////    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
////    theCanvas->SetGrayscale();
//
//    flux2DHist.Draw("CONTZ");
//    theCanvas->SaveAs(GraphsDir+"neutronBeamProfileSim140619_pretty.pdf");

///Make BGO attenuation graph
//makeBGOSingleEnergyGraphs();

///Make PLO Graph
//    double kk0[12]={0.9994634,0.9996557,0.999513,0.9963012,0.9978678,0.9997649,0.998882,0.9991077,0.9932819,0.9982696,0.9984449,0.9987223};
//    double kk1[12]={3.986864e-05,0.0001689317,5.064656e-05,7.003159e-05,-4.994234e-05,-0.0002110803,-3.110464e-05,-0.000195551,-0.0002297314,-3.916408e-05,2.321166e-05,-1.501302e-05};
//    double kk2[12]={2.370452e-06,3.815403e-07,2.070204e-06,1.066584e-05,6.392543e-06,1.425570e-06,4.000117e-06,3.397849e-06,2.053844e-05,5.875048e-06,5.764383e-06,5.524947e-06};
//    double kk3[12]={1.751450e-08,-1.849636e-08,-1.176445e-09,2.041891e-08,6.756355e-08,1.650474e-09,3.531353e-08,5.803961e-08,3.642482e-08,6.059123e-08,5.107648e-08,9.335619e-08};
//    double kk4[12]={-5.952597e-12,9.001259e-10,1.595719e-10,-6.976799e-10,5.757778e-11,7.425284e-10,1.138028e-10,-2.687829e-10,-4.223121e-09,-6.766136e-10,-1.236748e-09,-1.018134e-09};
//    double kk5[12]={-2.620916e-12,2.281255e-12,-1.274924e-12,-4.08466e-12,-8.450335e-12,3.617150e-12,-4.012398e-12,-5.354198e-12,-7.169351e-13,-8.380383e-12,-3.582307e-12,-1.27932e-11};
//    double kk6[12]={6.048114e-14,-3.857624e-14,1.131760e-14,7.401177e-14,3.295636e-14,-1.26734e-14,2.109971e-14,3.401312e-14,3.269036e-13,6.981891e-14,1.798018e-13,1.564707e-13};
//    double norm_func[12]={1.015887,1.013419,1.011222,1.028474,1.025036,1.017557,1.017506,1.009917,1.023982,1.014295,1.018611,1.019129};
//
//    TH1D* ploHist=new TH1D("ploHist",";Distance from Center of Crystal (mm);Relative Response to Crystal Center",101,-100,100);
//    double gain=0;
//    double zTrans;
//    for (int i = 0;i< 101;i++)
//    {
//        zTrans=-100+i*2;
//        for (int DetNum = 0;DetNum< 12;DetNum++)
//        {
//            gain +=(1./12.)*(kk0[DetNum] + kk1[DetNum]*zTrans + kk2[DetNum]*pow(zTrans,2) + kk3[DetNum]*pow(zTrans,3) + kk4[DetNum]*pow(zTrans,4) + kk5[DetNum]*pow(zTrans,5) + kk6[DetNum]*pow(zTrans,6));
//        }
//
//        ploHist->SetBinContent(i+1,gain);
//
//        gain=0;
//    }
//
//    convertTH1ToTXT(ploHist,HistsDir+"averagePLOHist140610.txt");
//
//    TCanvas* theCanvas = getThesisPlotCanvas();
//    ploHist->Draw("L");
//    theCanvas->SaveAs(GraphsDir+"averagePLOHist140610.pdf");

///Plot Cs137Co57
//    TH1D* calHist= getTabSeperatedHist(HistsDir+"Cs137Co57_S502_Det1.txt",2018,-.5,2017.5);
//    calHist->SetTitle(";Uncalibrated BGO Peak Height Channel;Counts");
//    calHist->Rebin(4);
//    calHist->SetMarkerColor(kRed);
//    calHist->SetLineColor(kRed);
//    TCanvas* theCanvas = getThesisPlotCanvas();
//    calHist->Draw("L");
//    theCanvas->SaveAs(GraphsDir+"Cs137Co57140610.pdf");

///Create central light response curve with correct region by region experiment weighting with Averkiev
//    const int numGraphs=4;
//
//
//
//    const int numDewarPoints=8;
//    double dewarEnergies[numDewarPoints]={23.3,31.7,59.5,81.0,88.0,99.8,123.7,185.7};
//    double dewarLOs[numDewarPoints]={.734,.812,.888,.916,.907,.865,.892,.933};
//    double dewarLOErrors[numDewarPoints]={.03,.025,.01,.015,.01,.015,.015,.02};
//    TGraphErrors dewarGraph(numDewarPoints,dewarEnergies,dewarLOs,NULL,dewarLOErrors);
//    dewarGraph.SetName("DewarGraph");
//    TString dewarLegendString="RDK II BGO in Dewar [MCA]";
//    dewarGraph.SetLineColor(kRed);
//    dewarGraph.SetMarkerColor(kRed);
//    dewarGraph.SetMarkerStyle(kFullCircle);
//    dewarGraph.SetMarkerSize(1.3);
//
//    const int numCalRDK2Points=6;
//    double calRDK2Energies[numCalRDK2Points]={32.3,42.7,78.7,79.7,121.9,356};
//    double calRDK2LOs[numCalRDK2Points]={0.81,0.87,0.94,0.93,0.92,0.986};
//    double calRDK2LOErrors[numCalRDK2Points]={0.02,0.02,0.01,0.01,0.01,0.002};
//    TGraphErrors calRDK2Graph(numCalRDK2Points,calRDK2Energies,calRDK2LOs,NULL,calRDK2LOErrors);
//    calRDK2Graph.SetName("CalRDK2Graph");
//    TString calRDK2LegendString="RDK II BGO Array [DAQ]";
//    calRDK2Graph.SetLineColor(kViolet);
//    calRDK2Graph.SetMarkerColor(kViolet);
//    calRDK2Graph.SetMarkerStyle(kFullSquare);
//    calRDK2Graph.SetMarkerSize(1.3);
//
//    //TGraph* khodyukGraph=getKhodyukBGOResponseGraphOfficial();
//    //Assumed error was similar to moszynski and was 2%
//    const int numKhodyukPoints=375;
//    double khodyukEnergies[numKhodyukPoints]={11,11,12,12,13,13.025,13.05,13.075,13.1,13.125,13.15,13.175,13.2,13.225,13.25,13.275,13.3,13.325,13.35,13.375,
//        13.4,13.425,13.45,13.475,13.5,13.525,13.55,13.575,13.6,13.625,13.65,13.675,13.7,13.725,13.75,13.775,13.8,13.825,13.85,13.875,13.9,13.925,13.95,13.975,14,
//        14.025,14.05,14.075,14.1,14.125,14.15,14.175,14.2,14.225,14.25,14.275,14.3,14.325,14.35,14.375,14.4,14.425,14.45,14.475,14.5,14.525,14.55,14.575,14.6,14.625,
//        14.65,14.675,14.7,14.725,14.75,14.775,14.8,14.825,14.85,14.875,14.9,14.925,14.95,14.975,15,15.025,15.05,15.075,15.1,15.125,15.15,15.175,15.2,15.225,15.25,
//        15.275,15.3,15.325,15.35,15.375,15.4,15.425,15.45,15.475,15.5,15.525,15.55,15.575,15.6,15.625,15.65,15.675,15.7,15.725,15.75,15.775,15.8,15.825,15.85,15.875,
//        15.9,15.925,15.95,15.975,16,16.025,16.05,16.075,16.1,16.125,16.15,16.175,16.2,16.225,16.25,16.275,16.3,16.325,16.35,16.375,16.4,16.425,16.45,16.475,16.5,16.525,
//        16.55,16.575,16.6,16.625,16.65,16.675,16.7,16.725,16.75,16.775,16.8,16.825,16.85,16.875,16.9,16.925,16.95,16.975,17,17.025,17.05,17.075,17.1,17.125,17.15,17.175,
//        17.2,17.225,17.25,17.275,17.3,17.325,17.35,17.375,17.4,17.425,17.45,17.475,17.5,17.525,17.55,17.575,17.6,17.625,17.65,17.675,17.7,17.725,17.75,17.775,17.8,17.825,
//        17.85,17.875,17.9,17.925,17.95,17.975,18,18.025,18.05,18.075,18.1,18.125,18.15,18.175,18.2,18.225,18.25,18.275,18.3,18.325,18.35,18.375,18.4,18.425,18.45,18.475,
//        18.5,18.525,18.55,18.575,18.6,18.625,18.65,18.675,18.7,18.725,18.75,18.775,18.8,18.825,18.85,18.875,18.9,18.925,18.95,18.975,19,19.025,19.05,19.075,19.1,19.125,
//        19.15,19.175,19.2,19.225,19.25,19.275,19.3,19.325,19.35,19.375,19.4,19.425,19.45,19.475,19.5,19.525,19.55,19.575,19.6,19.625,19.65,19.675,19.7,19.725,19.75,19.775,
//        19.8,19.825,19.85,19.875,19.9,19.925,19.95,19.975,20,20,25,25,30,30,30,30,35,35,40,40,45,45,45,45,50,50,55,55,60,60,65,65,70,70,75,75,80,80,85,85,89.5,89.6,89.7,
//        89.8,89.9,90,90.1,90.2,90.3,90.4,90.5,90.6,90.7,90.8,90.9,91,91.1,91.2,91.3,91.4,91.5,91.6,91.7,91.8,91.9,92,92.1,92.2,92.3,92.4,92.5,92.6,92.7,92.8,92.9,93,93.1,
//        93.2,93.3,93.4,93.5,93.6,93.7,93.8,93.9,94,94.1,94.2,94.3,94.4,94.5,94.6,94.7,94.8,94.9,95,95,100,100};
//    double khodyukLOs[numKhodyukPoints]={0.700,0.707,0.725,0.724,0.740,0.742,0.737,0.743,0.743,0.741,0.737,0.740,0.738,0.736,0.743,0.746,0.743,0.739,0.742,0.742,
//        0.736,0.722,0.703,0.704,0.705,0.695,0.703,0.701,0.698,0.699,0.702,0.702,0.704,0.695,0.701,0.696,0.705,0.696,0.706,0.701,0.696,0.693,0.697,0.697,0.694,0.700,0.699,0.698,
//        0.699,0.702,0.697,0.692,0.693,0.696,0.697,0.689,0.697,0.683,0.684,0.693,0.687,0.698,0.695,0.702,0.699,0.697,0.697,0.699,0.697,0.699,0.694,0.691,0.699,0.691,0.692,0.694,
//        0.694,0.702,0.697,0.699,0.701,0.695,0.698,0.699,0.702,0.698,0.697,0.700,0.702,0.698,0.697,0.698,0.701,0.708,0.703,0.697,0.699,0.698,0.701,0.705,0.704,0.706,0.702,0.700,
//        0.704,0.704,0.706,0.712,0.701,0.705,0.704,0.704,0.702,0.698,0.697,0.702,0.703,0.698,0.699,0.693,0.701,0.700,0.701,0.700,0.694,0.697,0.698,0.698,0.699,0.695,0.701,0.704,
//        0.698,0.701,0.697,0.701,0.697,0.700,0.694,0.697,0.691,0.687,0.684,0.690,0.692,0.685,0.689,0.687,0.693,0.686,0.694,0.697,0.690,0.691,0.690,0.688,0.695,0.690,0.690,0.690,
//        0.690,0.694,0.693,0.693,0.690,0.688,0.692,0.687,0.689,0.692,0.687,0.694,0.695,0.692,0.693,0.693,0.698,0.694,0.695,0.692,0.699,0.694,0.693,0.695,0.695,0.696,0.692,0.694,
//        0.690,0.694,0.697,0.698,0.700,0.700,0.700,0.696,0.701,0.696,0.696,0.698,0.704,0.705,0.696,0.694,0.698,0.701,0.702,0.703,0.697,0.701,0.702,0.701,0.698,0.708,0.702,0.702,
//        0.701,0.702,0.699,0.706,0.705,0.700,0.714,0.699,0.710,0.705,0.704,0.702,0.704,0.703,0.703,0.703,0.705,0.706,0.703,0.701,0.712,0.701,0.710,0.707,0.709,0.704,0.709,0.713,
//        0.706,0.711,0.710,0.709,0.719,0.709,0.705,0.711,0.712,0.711,0.708,0.709,0.709,0.710,0.710,0.716,0.716,0.712,0.709,0.717,0.713,0.715,0.711,0.716,0.716,0.717,0.714,0.715,
//        0.711,0.720,0.714,0.718,0.715,0.711,0.720,0.716,0.717,0.722,0.722,0.724,0.720,0.718,0.768,0.770,0.803,0.800,0.800,0.800,0.827,0.826,0.844,0.845,0.861,0.864,0.858,0.858,
//        0.875,0.874,0.880,0.880,0.885,0.885,0.890,0.889,0.896,0.893,0.897,0.898,0.901,0.895,0.906,0.909,0.911,0.911,0.912,0.910,0.911,0.910,0.909,0.910,0.908,0.908,0.899,0.885,
//        0.885,0.884,0.885,0.881,0.883,0.880,0.882,0.881,0.879,0.879,0.878,0.881,0.879,0.879,0.880,0.880,0.882,0.882,0.877,0.878,0.878,0.880,0.876,0.879,0.880,0.881,0.879,0.878,
//        0.878,0.879,0.878,0.880,0.880,0.879,0.878,0.878,0.877,0.878,0.880,0.877,0.878,0.878,0.878,0.881,0.880,0.882,0.883};
//    TGraph khodyukGraph(numKhodyukPoints,khodyukEnergies,khodyukLOs);
//    TString khodyukLegendString="I. Khodyuk and P. Dorenbos";
//    khodyukGraph.SetLineColor(kBlue);
//    khodyukGraph.SetMarkerColor(kBlue);
//    khodyukGraph.SetMarkerStyle(kFullCircle);
//
//
//    const int numMoszynskiPoints=30;
//    double moszynskiEnergies[numMoszynskiPoints]={14.4,16.6,22,31,32,45,52,59.6,74,81,84,88,122,244.1,302.9,320.1,335.8,356,411.8,435,493.7,511,569.7,586,662,758.8,834.8,1063.7,1173.2,1275};
//    double moszynskiLOs[numMoszynskiPoints]={0.711,0.719,0.775,0.829,0.8412,0.853,0.9,0.893,0.91,0.926,0.937,0.933,0.915,0.964,0.979,0.964,0.97,0.986,0.97,0.995,1,0.996,1,0.999,1,1,1,1.01,1.013,1.01};
//    double moszynskiLOErrors[numMoszynskiPoints]={0.01422,0.01438,0.0155,0.01658,0.016824,0.02559,0.018,0.01786,0.0182,0.01852,0.01874,0.01866,0.0183,0.02892,0.01958,0.01928,0.0291,0.01972,0.0194,0.02985,0.03,0.01992,0.02,0.02997,0.02,0.03,0.02,0.0202,0.02026,0.0202};
//    TGraphErrors moszynskiGraph(numMoszynskiPoints,moszynskiEnergies,moszynskiLOs,NULL,moszynskiLOErrors);
//    moszynskiGraph.SetName("MoszynskiGraph");
//    TString moszynskiLegendString="Moszynski et. al.";
//    moszynskiGraph.SetLineColor(kGreen);
//    moszynskiGraph.SetMarkerColor(kGreen);
//    moszynskiGraph.SetMarkerStyle(kFullTriangleUp);
//    moszynskiGraph.SetMarkerSize(1.3);
//
////    const int numAverkievPoints=20;
////    double averkievEnergies[numAverkievPoints]={2.1,6.4,7.5,8.0,8.8,10.0,11.2,16.6,17.2,23.5,33.1,58.3,71.6,78.3,120.1,162.3,277.4,384.8,506.5,657.9};
////    double averkievLOs[numAverkievPoints]={0.854,0.777,0.739,0.706,0.678,0.657,0.641,0.579,0.590,0.685,0.731,0.837,0.844,0.827,0.884,0.903,0.959,0.977,0.986,0.996};
////    double averkievLOErrors[numAverkievPoints]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.035,0.035,0.035,0.035,0.035,0.02,0.015,0.01,0.005,0.005};
////    TGraphErrors averkievGraph(numAverkievPoints,averkievEnergies,averkievLOs,NULL,averkievLOErrors);
////    averkievGraph.SetName("AverkievGraph");
////    TString averkievLegendString="Averkiev et. al.";
////    averkievGraph.SetLineColor(kMagenta);
////    averkievGraph.SetMarkerColor(kMagenta);
////    averkievGraph.SetMarkerStyle(kOpenCircle);
////    averkievGraph.SetMarkerSize(1.5);
//
//
//    double legendSize=.05*numGraphs;
//    TLegend* theLegend = new TLegend(0.6,.15,0.9,0.5);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
////    theLegend->AddEntry(&quadraticGraph,quadraticLegendString,"L");
//    theLegend->AddEntry(&khodyukGraph,khodyukLegendString,"L");
//    theLegend->AddEntry(&moszynskiGraph,moszynskiLegendString,"P");
////    theLegend->AddEntry(&averkievGraph,averkievLegendString,"P");
//    theLegend->AddEntry(&dewarGraph,dewarLegendString,"P");
//    theLegend->AddEntry(&calRDK2Graph,calRDK2LegendString,"P");
//
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
//
////    khodyukGraph.GetYaxis()->SetTitleOffset(1.3);
//    khodyukGraph.GetXaxis()->SetLimits(0,1000);
//    khodyukGraph.GetYaxis()->SetRangeUser(.65,1.05);
//    khodyukGraph.SetTitle(";Energy Absorbed in BGO (keV);Light Output Relative to 662 keV");
//
//
//    khodyukGraph.Draw("A L");
////    quadraticGraph.Draw("L same");
//    moszynskiGraph.Draw("P E1 same");
////    averkievGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//    theLegend->Draw("same");
//    theCanvas->SaveAs(GraphsDir+"LOExperimentsGraph_140619.pdf");
//    khodyukGraph.GetXaxis()->SetLimits(0,120);
//    theCanvas->SaveAs(GraphsDir+"LOExperimentsGraphZoom_140619.pdf");

///LAAPD Low energy response
//    TString rootFilePath="/home/mjbales/school/rdk/misc/BNLUC3BeamlineApril2011/BNL_UC3_Beamline_1104_hists.root";
//    TFile rootFile(rootFilePath,"READ");
//
//    TH1I* hist880=(TH1I*) rootFile.Get("apr11_9_14");
//    TH1I* hist600=(TH1I*) rootFile.Get("apr11_9_41");
//    TH1I* histBackground1=(TH1I*) rootFile.Get("apr11_9_34"); //taken for 3 times the time as others
//
//    TH1I* histBackground2=(TH1I*) rootFile.Get("apr11_10_11");
//    TH1I* hist400=(TH1I*) rootFile.Get("apr11_10_6");
//    TH1I* hist500=(TH1I*) rootFile.Get("apr11_10_15");
//    TH1I* hist1400=(TH1I*) rootFile.Get("apr11_10_23");
//
//    TH1I  hist600BS("hist600",";Uncalibrated LAAPD Peak Height Response;Normalized Counts",1024,-.5,1023.5);
//    TH1I  hist880BS("hist880",";Uncalibrated LAAPD Peak Height Response;Normalized Counts",1024,-.5,1023.5);
//    TH1I  hist400BS("hist400",";Uncalibrated LAAPD Peak Height Response;Normalized Counts",1024,-.5,1023.5);
//    TH1I  hist500BS("hist500",";Uncalibrated LAAPD Peak Height Response;Normalized Counts",1024,-.5,1023.5);
//    TH1I  hist1400BS("hist1400",";Uncalibrated LAAPD Peak Height Response;Normalized Counts",1024,-.5,1023.5);
//
//    hist880BS.Add(hist880); hist880BS.Add(histBackground1,-.333);
//    hist600BS.Add(hist600); hist600BS.Add(histBackground1,-.333);
//    hist1400BS.Add(hist1400); hist1400BS.Add(histBackground2,-1);
//    hist500BS.Add(hist500); hist500BS.Add(histBackground2,-1);
//    hist400BS.Add(hist400); hist400BS.Add(histBackground2,-1);
//
//    rootFile.Close();
//
//    hist400BS.SetLineColor(kBlack);
//    hist500BS.SetLineColor(kBlue);
//    hist600BS.SetLineColor(kRed);
//    hist880BS.SetLineColor(kGreen+2);
//    hist1400BS.SetLineColor(kViolet);
//
//    hist400BS.SetLineStyle(1);
//    hist500BS.SetLineStyle(2);
//    hist600BS.SetLineStyle(3);
//    hist880BS.SetLineStyle(4);
//    hist1400BS.SetLineStyle(6);
//
//    hist400BS.SetLineWidth(1);
//    hist500BS.SetLineWidth(2);
//    hist600BS.SetLineWidth(2);
//    hist880BS.SetLineWidth(2);
//    hist1400BS.SetLineWidth(2);
//
//
//    scaleHistogramsTogether(&hist400BS,&hist500BS);
//    scaleHistogramsTogether(&hist400BS,&hist1400BS);
//    scaleHistogramsTogether(&hist400BS,&hist600BS);
//    scaleHistogramsTogether(&hist400BS,&hist880BS);
//
//    hist400BS.GetSumw2()->Set(0);
//    hist500BS.GetSumw2()->Set(0);
//    hist600BS.GetSumw2()->Set(0);
//    hist880BS.GetSumw2()->Set(0);
//    hist1400BS.GetSumw2()->Set(0);
//
//    hist400BS.GetXaxis()->SetRangeUser(15,180);
//    hist400BS.GetYaxis()->SetRangeUser(0,3000);
//
//    TCanvas* theCanvas = getThesisPlotCanvas();
//
//    hist400BS.Draw("L");
//    hist500BS.Draw("L same");
//    hist600BS.Draw("L same");
//    hist880BS.Draw("L same");
//    hist1400BS.Draw("L same");
//
//     TLegend* theLegend = new TLegend(0.6,.55,0.9,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
//
//    theLegend->AddEntry(&hist400BS,"400 eV","L");
//    theLegend->AddEntry(&hist500BS,"500 eV","L");
//    theLegend->AddEntry(&hist600BS,"600 eV","L");
//    theLegend->AddEntry(&hist880BS,"880 eV","L");
//    theLegend->AddEntry(&hist1400BS,"1400 eV","L");
//    theLegend->Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"laapdMonochrometerExamples140619.pdf");

//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID420", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID421", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID422", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID423", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID424", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID425", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID426", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID427", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID428", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID429", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID430", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID431", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID432", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID433", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID434", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID435", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID436", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID437", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID438", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID439", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID440", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID441", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID442", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID443", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 1. 0. mm\n/MRK/det/setGammaDetOffset 0. 1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y+ gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID444", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID445", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID446", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID447", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID448", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
//makeGeantMacroFilesForClusterEventGenerator(10,38, "RID449", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract

///Make BAPD Depth model graph

//    TH1D ceHist("cehist",";Energy Absorbtion Depth (#mum);Signal Collection Efficiency",601,-0.1,60.1);
//
//    for (int i = 0;i< 601;i++)
//    {
//        ceHist.SetBinContent(i+1,BAPDCollectionEfficiencyModel(500,i*.01));
//    }
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
//
//    ceHist.GetXaxis()->SetRangeUser(0,60);
//    ceHist.GetYaxis()->SetRangeUser(0,1.1);
//
//    ceHist.Draw("L");
//    theCanvas->SaveAs(GraphsDir+"laapdCEHist140620.pdf");

///initial event plots
///Make hists
//    HistDim ep0Dim={100,0.,1};
//    HistDim bgoDim={790,10.,800};
//    int numEventFiles=48;
//    RDK2MCAnalysis theAnalysis("AID48",  "20 nm Si Dead Layer on SBD with P Resolution", numEventFiles,  "RID358", "RID359",  numEventFiles, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
////
////    TH1* ep03BHist=theAnalysis.MakeHist("ep0Hist",THREEBODY,"ep0","",ep0Dim);
////    ep03BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ep03BHist,HistsDir+"ep0Hist3B140610.txt");
////
////    TH1* ep04BHist=theAnalysis.MakeHist("ep0Hist2",FOURBODY,"ep0","eg0 > 14",ep0Dim);
////    ep04BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ep04BHist,HistsDir+"ep0Hist4B140610.txt");
////
////    TH1* ee03BHist=theAnalysis.MakeHist("ee0Hist",THREEBODY,"ee0","",EDEPE_DIM);
////    ee03BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ee03BHist,HistsDir+"ee0Hist3B140610.txt");
////
////    TH1* ee04BHist=theAnalysis.MakeHist("ee0Hist2",FOURBODY,"ee0","eg0 > 14",EDEPE_DIM);
////    ee04BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ee04BHist,HistsDir+"ee0Hist4B140610.txt");
////
////    TH1* eg04BHist=theAnalysis.MakeHist("eg0Hist",FOURBODY,"eg0","eg0 > 10",bgoDim);
////    eg04BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(eg04BHist,HistsDir+"eg0Hist4B140610.txt");
//
//    TH1D* ep03BHist=getTabSeperatedHist(HistsDir+"ep0Hist3B140610.txt",ep0Dim);
//    TH1D* ep04BHist=getTabSeperatedHist(HistsDir+"ep0Hist4B140610.txt",ep0Dim);
//    TH1D* ee03BHist=getTabSeperatedHist(HistsDir+"ee0Hist3B140610.txt",EDEPE_DIM);
//    TH1D* ee04BHist=getTabSeperatedHist(HistsDir+"ee0Hist4B140610.txt",EDEPE_DIM);
//    TH1D* eg04BHist=getTabSeperatedHist(HistsDir+"eg0Hist4B140610.txt",bgoDim);
//
//
//    ee03BHist->Rebin(10);
//    ee04BHist->Rebin(10);
//    eg04BHist->Rebin(10);
//    ee03BHist->Scale(10);
//    ep03BHist->Scale(10);
//
//    scaleHistogramsTogether(ep03BHist,ep04BHist);
//    scaleHistogramsTogether(ee03BHist,ee04BHist);
//
//    TCanvas* theCanvas = getThesisPlotCanvas();
////    ep0Hist->Smooth(10);
//    ep03BHist->SetLineColor(kBlue);
//    ep04BHist->SetLineColor(kRed);
//    ep04BHist->SetLineStyle(2);
//    ep04BHist->SetLineWidth(2);
//
//    ep03BHist->SetTitle(";Initial Kinetic Energy from Decay (keV);Emission Probability (AU)");
//    ep03BHist->GetSumw2()->Set(0);
//    ep04BHist->GetSumw2()->Set(0);
////    ep0Hist->GetXaxis()->SetRangeUser(-1,1);
//    ep03BHist->Draw("L");
//    ep04BHist->Draw("L same");
//
//    TLegend* theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(ep03BHist,"Three Body Decays","L");
//    theLegend->AddEntry(ep04BHist,"Four Body Decays","L");
//    theLegend->Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"ep0Hist140610.pdf");
//    delete theLegend;
//
//    ee03BHist->SetLineColor(kBlue);
//    ee04BHist->SetLineColor(kRed);
//    ee04BHist->SetLineStyle(2);
//    ee04BHist->SetLineWidth(2);
//
//    ee03BHist->SetTitle(";Initial Kinetic Energy from Decay (keV);Emission Probability (AU)");
//    ee03BHist->GetSumw2()->Set(0);
//    ee04BHist->GetSumw2()->Set(0);
//    ee03BHist->GetYaxis()->SetRangeUser(0,0.25);
//    ee03BHist->Draw("L");
//    ee04BHist->Draw("L same");
//
//    theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(ee03BHist,"Three Body Decays","L");
//    theLegend->AddEntry(ee04BHist,"Four Body Decays","L");
//    theLegend->Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"ee0Hist140610.pdf");
//    delete theLegend;
//
//    eg04BHist->SetLineColor(kRed);
//    eg04BHist->SetLineStyle(2);
//    eg04BHist->SetLineWidth(2);
//
//    eg04BHist->SetTitle(";Initial Kinetic Energy from Decay (keV);Emission Probability (AU)");
//    eg04BHist->GetSumw2()->Set(0);
//
//    eg04BHist->Draw("L");
//
//
//    theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(eg04BHist,"Four Body Decays","L");
//    theLegend->Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"eg0Hist140610.pdf");
//    delete theLegend;

///Make Exp Sim comparison plots

//    RDK2AnalysisPlotter* thePlotter;
//    TH1* expHist;
//    TH1* mcHist;
//
//
//    thePlotter=new RDK2AnalysisPlotter("Latest and Greatest","EID3","AID48");
//
//    CoDet detType;
//    PlotVarType varType;
//
////    TCanvas* theCanvas = getThesisPlotCanvas();
////    TCanvas* theCanvas = getHalfPresentationPlot2();
//    TCanvas* theCanvas = getHalfPresentationPlot();
//
//    for (int i = 0;i< 2;i++)
//    {
//        detType=static_cast<CoDet>(i+1);
//        for (int j = 0;j< 4;j++)
//        {
//            varType=static_cast<PlotVarType>(j);
//
//            if(detType==DET_EP && varType==PLOTVAR_GE)
//                continue;
//
//            expHist = thePlotter->getExpHist(0,detType, varType);
//            mcHist = thePlotter->getMCHist(0,detType, varType);
//
//            Color_t theColor;
//            TString yAxisString;
//            TString imageString="compPlot_";
//            if(detType==DET_EP)
//            {
//                theColor=kBlue;
//                yAxisString="Counts (AU)";
//                imageString+="EP_";
//            }
//            else if(detType==DET_EPG)
//            {
//                theColor=kRed;
//                yAxisString="#frac{BGO ep#gamma counts}{ep counts}";
////                mcHist->SetLineStyle(2);
//                mcHist->SetLineWidth(2);
//                imageString+="EPG_";
//            }
//            else if(detType==DET_EPBG)
//            {
//                theColor=kGreen+2;
//                yAxisString="#frac{LAAPD ep#gamma counts}{ep counts}";
////                mcHist->SetLineStyle(3);
//                mcHist->SetLineWidth(2);
//                imageString+="EPBG_";
//            }
//
//            TString xAxisString;
//            int rebin;
//            if(varType==PLOTVAR_EE)
//            {
//                xAxisString="Calibrated Electron Peak Height (781.5=781.5 keV)";
//                imageString+="EE";
//                rebin=10;
//                if(detType!=DET_EP)
//                {
//                    rebin=20;
//                }
//
//            }
//            else if(varType==PLOTVAR_PE)
//            {
//                xAxisString="Calibrated Proton Peak Height (781.5=781.5 keV)";
//                imageString+="PE";
//                rebin=2;
//            }
//            else if(varType==PLOTVAR_PT)
//            {
//                xAxisString="Proton Time-of-Flight (s)";
//                imageString+="PT";
//                mcHist->GetXaxis()->SetRangeUser(2e-6, 10e-6);
//                rebin=9;
//            }
//            else if(varType==PLOTVAR_GE)
//            {
//                imageString+="GE";
//                if(detType==DET_EPG)
//                {
//                   xAxisString="Calibrated Photon Peak Height (511=511 keV)";
//                }
//                else if(detType==DET_EPBG)
//                {
//                    xAxisString="Calibrated Photon Peak Height (5.6=5.6 keV)";
//                }
//                rebin=10;
//            }
//
//            expHist->SetMarkerColor(kBlack);
//            expHist->SetLineColor(kBlack);
//            mcHist->SetLineColor(theColor);
//            mcHist->GetSumw2()->Set(0);
//            expHist->SetMarkerStyle(4);
//            mcHist->GetYaxis()->SetTitleOffset(1.7);
////            theCanvas->SetLeftMargin(0.2);
//
//            mcHist->Rebin(rebin);
//            expHist->Rebin(rebin);
//
//            gStyle->SetErrorX(0.0001);
//
//            mcHist->SetTitle(";"+xAxisString+";"+yAxisString);
//
//            scaleHistogramsTogether(expHist,mcHist);
//
//            if(varType==PLOTVAR_PT)
//            {
//
//                mcHist->GetXaxis()->SetRangeUser(2e-6, 10e-6);
//
//            }
//
//            mcHist->Draw("L");
//            expHist->Draw("P E0 Same");
//
//            TLegend* theLegend = new TLegend(0.70,.55,0.90,0.9);
//            theLegend->SetBorderSize(1);
//            theLegend->SetFillColor(kWhite);
//            theLegend->AddEntry(expHist,"Experiment","PE");
//            theLegend->AddEntry(mcHist,"Simulation","L");
//            theLegend->Draw("same");
//
////            imageString+="_140621.pdf";
//            imageString+="_140621c.png";
//            theCanvas->SaveAs(GraphsDir+imageString);
//            delete mcHist;
//            delete expHist;
//            delete theLegend;
//        }
//    }

///Make BGO attenuation graph

//    double x[40];
//    double y[40];
//
//    for (int i = 0;i< 40;i++)
//    {
//        x[i]=20.*(i+1);
//        TFile rootFile(ResultsDir+"Results_DirectBGO_140619_" +int2str(20*(i+1))+"keV.root");
//        TTree* theTree=(TTree*) rootFile.Get("g");
//        y[i]=theTree->Draw("BGO1EDepTotal","BGO1EDepTotal > 14");
//        y[i]/=100000;
//        rootFile.Close();
//    }
//
//    TGraph theGraph(40,x,y);
//    theGraph.SetTitle(";Incident Photon Energy(keV);Probability of BGO Energy Deposit (>14 keV)");
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
//
//    theGraph.SetLineColor(kBlue);
//    theGraph.GetXaxis()->SetRangeUser(0,800);
//    theGraph.GetYaxis()->SetRangeUser(0,1.05);
//    theGraph.Draw("AL");
//
//    theCanvas->SaveAs(GraphsDir+"BGOAttenuation140625.pdf");

///Make BAPD Plots
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"EXPBen_Hist131009_BAPD_AverageIndiv_EPG.txt",BAPD_DIM);
//
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248.txt",BAPD_DIM);
//    expHist->SetBinContent(2,0);
//    mcHist->SetBinContent(2,0);
//    mcHist->SetBinContent(1,0);
////cout << expHist->Integral() << endl;
////cout << mcHist->Integral() << endl;
////
//    scaleHistogramsTogether(expHist,mcHist);
//
//
//
//
////    TCanvas* theCanvas=getThesisPlotCanvas();
//    TCanvas* theCanvas=getHalfPresentationPlot2();
//    expHist->SetMarkerColor(kBlack);
//    expHist->SetLineColor(kBlack);
//    mcHist->SetLineColor(kGreen+2);
//    mcHist->GetSumw2()->Set(0);
//    expHist->SetMarkerStyle(4);
////    mcHist->SetLineStyle(3);
//    mcHist->SetLineWidth(2);
//
//    mcHist->Rebin(2);
//    expHist->Rebin(2);
//
//    gStyle->SetErrorX(0.0001);
//
//    mcHist->SetTitle(";Calibrated Photon Peak Height (5.9=5.9 keV);#frac{LAAPD ep#gamma counts}{ep counts}");
//
//
//    mcHist->GetYaxis()->SetRangeUser(0,5.5e-6);
//
//    mcHist->Draw("L");
//    expHist->Draw("P E0 Same");
//
//    TLegend* theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(expHist,"Experiment","PE");
//    theLegend->AddEntry(mcHist,"Simulation","L");
//    theLegend->Draw("same");
//
////     TLine line(0,0,20,0);
////     line.Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"compPlot_EPBG_BGE_140701b.png");
//
//    TH1D* ratioHist = (TH1D*) expHist->Clone("ratioHist");
//    ratioHist->Reset();
//    ratioHist->Sumw2();
//    ratioHist->SetTitle(";Calibrated Photon Peak Height (5.9=5.9 keV);Normalized Residual #frac{Exp - Sim}{Sim}");
//
//    ratioHist->Add(expHist,mcHist,1,-1);
//    ratioHist->Divide(mcHist);
////    ratioHist->Divide(expHist,mcHist);
////    ratioHist->Add(new TF1("constantn1","-1",0,20)); //subtract one for (exp-mc)/mc
//
//    ratioHist->SetLineColor(kGreen+2);
//    ratioHist->SetMarkerColor(kGreen+2);
//    ratioHist->GetYaxis()->SetRangeUser(-1.5,1.5);
//
//    ratioHist->Draw("P E0");
//
////    line.Draw("same");
//    ratioHist->Draw("P E0 same");
//
////    theCanvas->SaveAs(GraphsDir+"compPlot_EPBG_BGE_NResid_140701.pdf");
////    theCanvas->SaveAs(GraphsDir+"compPlot_EPBG_BGE_NResid_140701b.png");

///BGO Var Hists
//    TH1D* protoHist=new TH1D("epg_gEAvgVarHist","epg_gEAvgHist",20,BGO_VARHIST_LOWEDGES);
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"VarHist_MC_AID48_EPG_EDepGAvg.txt",protoHist);
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"VarHist_Exp_EID3_EPG_EDepGAvg.txt",protoHist);
////    expHist->SetBinContent(2,0);
////    mcHist->SetBinContent(2,0);
//
//    scaleHistogramsTogether(expHist,mcHist);
//
//
//
//
////    TCanvas* theCanvas=getThesisPlotCanvas();
//    TCanvas* theCanvas=getHalfPresentationPlot2();
//    expHist->SetMarkerColor(kBlack);
//    expHist->SetLineColor(kBlack);
//    mcHist->SetLineColor(kGreen+2);
////    mcHist->GetSumw2()->Set(0);
//    expHist->SetMarkerStyle(4);
//
////    mcHist->Rebin(2);
////    expHist->Rebin(2);
//
//    gStyle->SetErrorX(0.0001);
//
//    mcHist->SetTitle(";Calibrated Photon Peak Height (5.6=5.6 keV);#frac{LAAPD ep#gamma counts}{ep counts}");
//
//
////    mcHist->GetYaxis()->SetRangeUser(0,5.5e-6);
//
//
//
//    TH1D* ratioHist = (TH1D*) expHist->Clone("ratioHist");
//    ratioHist->Reset();
//    ratioHist->Sumw2();
//    ratioHist->SetTitle(";Calibrated Photon Peak Height (511=511 keV);Normalized Residual #frac{Exp - Sim}{Sim}");
//
//    ratioHist->Add(expHist,mcHist,1,-1);
//    ratioHist->Divide(mcHist);
////    ratioHist->Divide(expHist,mcHist);
////    ratioHist->Add(new TF1("constantn1","-1",0,20)); //subtract one for (exp-mc)/mc
//
//    ratioHist->SetLineColor(kRed);
//    ratioHist->SetMarkerColor(kRed);
////    ratioHist->SetLineStyle(2);
//    ratioHist->SetLineWidth(2);
//    ratioHist->GetYaxis()->SetRangeUser(-0.3,.3);
//
//    ratioHist->Draw("P E0");
//    theCanvas->SetLogx();
//    TLine line(0,0,800,0);
//    line.Draw("same");
//    ratioHist->Draw("P E0 same");
//
////    theCanvas->SaveAs(GraphsDir+"compPlot_EPG_GE_NResid_140701.pdf");
//    theCanvas->SaveAs(GraphsDir+"compPlot_EPG_GE_NResid_140701b.png");

///BGO Timing Hist
//    RDK2ExpAnalysis expAnalysis( "EID3",  "Standard", "PID8", "G7_Cycle5_8", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//
////    TH1D* avgHist=makeTH1DFromDim("avgHist",";e-#gamma Timing (Channels);Counts",EGTIME_DIM);
////    for (int i = 0;i< 12;i++)
////    {
////    	avgHist->Add(expAnalysis.MakeHist("egtime"+int2str(i+1),"BGODeltaT"+int2str(i+1),STD_EPCutSetB.CreateExpEPCut(),EGTIME_DIM));
////    }
////    convertTH1ToTXT(avgHist,HistsDir+"ExpHist_EP_EGTime_140701.pdf");
//
//    TH1D* avgHist=getTabSeperatedHist(HistsDir+"ExpHist_EP_EGTime_140701.pdf",EGTIME_DIM);
//    avgHist->SetTitle(";e-#gamma Timing (Channels);Counts");
//
//    avgHist->SetLineColor(kRed);
//    avgHist->SetMarkerColor(kRed);
//    avgHist->SetMarkerStyle(4);
//    avgHist->Rebin(5);
//
//    avgHist->GetXaxis()->SetRangeUser(-500,500);
//
//    avgHist->Draw("P E0");
//
//    theCanvas->SaveAs(GraphsDir+"ExpPlot_EP_EGTime_140701.pdf");

///BGO EPG/EP rate per series
//    RDK2ExpAnalysis expAnalysis( "EID3",  "Standard", "PID8", "G7_Cycle5_8", STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//
//    double x[NUM_EXP_SERIES_G7_Cycle5_8];
//    double y[NUM_EXP_SERIES_G7_Cycle5_8];
//    double yError[NUM_EXP_SERIES_G7_Cycle5_8];
//
//    for (int i = 0;i< NUM_EXP_SERIES_G7_Cycle5_8;i++)
//    {
//        RDK2ExpAnalysis expAnalysis( "EID3",  "Standard", "PID8", "S"+int2str(EXP_SERIES_ID_LIST_G7_Cycle5_8[i]), STD_EPCutSetB, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//        x[i]=EXP_SERIES_ID_LIST_G7_Cycle5_8[i];
//        y[i]=expAnalysis.CalcEPGPerEP(DET_EPG,yError[i]);
//    }
//
//    TGraphErrors theGraph(NUM_EXP_SERIES_G7_Cycle5_8,x,y,NULL,yError);
//
//
//    theGraph.SetTitle(";Series Number;#frac{BGO ep#gamma counts}{ep counts}");
//
//    theGraph.SetLineColor(kRed);
//    theGraph.SetMarkerColor(kRed);
//    theGraph.SetMarkerStyle(4);
//
//    theGraph.GetXaxis()->SetRangeUser(-500,500);
//    theGraph.GetYaxis()->SetTitleOffset(1.3);
//    theCanvas->SetLeftMargin(.13);
//
//    theGraph.Draw("A P E0");
//
//    theCanvas->SaveAs(GraphsDir+"ExpPlot_EPGRatio_BySeries_140701.pdf");

///Make thesis plot of nonpropmodels
//    const int n=200;
//    double x[3][n];
//    double y[3][n];
//    EDepModel theModel[3]={ALLEXP_PARAM_MODEL,ALLEXP_PARAM_MODEL_HIGH,ALLEXP_PARAM_MODEL_LOW};
//    const int colorList[12]= {kBlack, kRed, kBlue, 600, 616, 432,800, 820, 840, 860, 870,880};
//    TString modelString[3]={"Fit","High (+3%)","Low (-3%)"};
//
//    TGraph theGraphs[3];
//    TCanvas* c2 = getThesisPlotCanvas();
//    TLegend *theLegend;
//    theLegend = new TLegend(0.8,0.2,0.93,0.6);
//
//    for (int i = 0;i< 3;i++)
//    {
//        for (int j = 0;j< n;j++)
//        {
//            x[i][j]=(double) j+1;
//            y[i][j]=evalMattAllExpParamModel2(x[i][j], theModel[i]);
//        }
//
//
//
//    	theGraphs[i]=TGraph(n,x[i],y[i]);
//
//    	theGraphs[i].SetLineColor(colorList[i]);
//    	theGraphs[i].SetLineStyle(i+1);
//    	if(i==0)
//    	{
//    	    theGraphs[i].Draw("AL");
//    	    theGraphs[i].GetYaxis()->SetRangeUser(0.3,1.0);
//    	    theGraphs[i].SetTitle(";Energy Absorbed in BGO (keV);Light output relative to 662 keV");
//    	}
//    	else
//    	{
//    	    theGraphs[i].SetLineWidth(2);
//    	    theGraphs[i].Draw("Lsame");
//    	}
//    	theLegend->AddEntry(&theGraphs[i],modelString[i],"L");
//    }
//
//    theLegend->Draw();
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
//    TString imagePath=GraphsDir+"ParamModelComparison140705.png";
//    c2->SaveAs(imagePath);
//    delete c2;

///BGO Timing Hist - double check
//    TCanvas* theCanvas = getThesisPlotCanvas();
//
//    RDK2ExpAnalysis expAnalysis( "EID3",  "Standard", "PID8", "G7_Cycle5_8", STD_EPCutSetA, STD_BGO_WINDOWS_A, STD_BAPD_WINDOWS_A);
//    TH1D* avgHist=makeTH1DFromDim("avgHist",";e-#gamma Timing (Channels);Counts",EGTIME_DIM2);
//    for (int i = 0;i< 12;i++)
//    {
//    	avgHist->Add(expAnalysis.MakeHist("egtime"+int2str(i+1),"BGODeltaT"+int2str(i+1),STD_EPCutSetA.CreateExpEPCut(),EGTIME_DIM2));
//    }
//    convertTH1ToTXT(avgHist,HistsDir+"ExpHist_EP_EGTime_140709.txt");
//
////    TH1D* avgHist=getTabSeperatedHist(HistsDir+"ExpHist_EP_EGTime_140709.txt",EGTIME_DIM2);
//
//    avgHist->SetTitle(";e-#gamma Timing (Channels);Counts");
//
//    avgHist->SetLineColor(kRed);
//    avgHist->SetMarkerColor(kRed);
////    avgHist->SetMarkerStyle(4);
////    avgHist->Rebin(5);
//
////    avgHist->GetXaxis()->SetRangeUser(-500,500);
//    avgHist->GetYaxis()->SetRangeUser(0,avgHist->GetMaximum()*1.1);
//
//    avgHist->Draw("P E1");
//
//    theCanvas->SaveAs(GraphsDir+"ExpPlot_EP_EGTime_140709.pdf");

///BGO Timing Hist - Make Jeff's latest
//    TCanvas* theCanvas = getThesisPlotCanvas();
//
//
//    TH1D* avgHist=getTabSeperatedHist(HistsDir+"ExpHist_Jeff_EGTime_140710.txt",1800,-650.5,1149.5);
//
//    avgHist->SetTitle(";e-#gamma Timing (Channels);Counts");
//
//    avgHist->SetLineColor(kRed);
//    avgHist->SetMarkerColor(kRed);
////    avgHist->SetMarkerStyle(4);
////    avgHist->Rebin(5);
//
//    avgHist->GetXaxis()->SetRangeUser(-600,600);
//    avgHist->GetYaxis()->SetRangeUser(0,2000);
//
//    avgHist->Draw("P E1");
//
//    TLine* line[4];
//    double xLine[4]={-475,-75,75,525};
//    for (int i =0;i< 4;i++)
//    {
//        line[i] = new TLine(xLine[i],0,xLine[i],2000);
//        line[i]->SetLineStyle(2);
//        line[i]->Draw("same");
//    }
//    avgHist->Draw("P same");
//
//    TPaveText textBox1(0.235786,0.7052632,0.3227425,0.9052632,"brNDC");
//    textBox1.SetTextSize(.03);
//    textBox1.AddText("Pre-peak (Background)");
//    textBox1.SetFillColor(kWhite);
//    textBox1.SetLineColor(kWhite);
//    textBox1.SetShadowColor(kWhite);
//    textBox1.SetLineWidth(0);
//    textBox1.Draw("same");
//
//
//    TPaveText textBox2(0.729097,0.8075188,0.8160535,0.8511278,"brNDC");
//    textBox2.SetTextSize(.03);
//    textBox2.AddText("Post-peak (Background)");
//    textBox2.SetFillColor(kWhite);
//    textBox2.SetLineColor(kWhite);
//    textBox2.SetShadowColor(kWhite);
//    textBox2.SetLineWidth(0);
//    textBox2.Draw("same");
//
//    TPaveText textBox3(0.4749164,0.1894737,0.5551839,0.2601504,"brNDC");
//    textBox3.SetTextSize(.03);
//    textBox3.AddText("On-peak");
//    textBox3.SetFillColor(kWhite);
//    textBox3.SetLineColor(kWhite);
//    textBox3.SetShadowColor(kWhite);
//    textBox3.SetLineWidth(0);
//    textBox3.Draw("same");
//
//
//
//
//    theCanvas->SaveAs(GraphsDir+"ExpPlot_EP_EGTime_140710.pdf");

///Create central light response curve with correct region by region experiment weighting with Averkiev
//    const int numGraphs=4;
//
//
//
//    const int numDewarPoints=8;
//    double dewarEnergies[numDewarPoints]={23.3,31.7,59.5,81.0,88.0,99.8,123.7,185.7};
//    double dewarLOs[numDewarPoints]={.734,.812,.888,.916,.907,.865,.892,.933};
//    double dewarLOErrors[numDewarPoints]={.03,.025,.01,.015,.01,.015,.015,.02};
//    TGraphErrors dewarGraph(numDewarPoints,dewarEnergies,dewarLOs,NULL,dewarLOErrors);
//    dewarGraph.SetName("DewarGraph");
//    TString dewarLegendString="RDK II BGO in Dewar [MCA]";
//    dewarGraph.SetLineColor(kRed);
//    dewarGraph.SetMarkerColor(kRed);
//    dewarGraph.SetMarkerStyle(kFullCircle);
//    dewarGraph.SetMarkerSize(1.3);
//
//    const int numCalRDK2Points=6;
//    double calRDK2Energies[numCalRDK2Points]={32.3,42.7,78.7,79.7,121.9,356};
//    double calRDK2LOs[numCalRDK2Points]={0.81,0.87,0.94,0.93,0.92,0.986};
//    double calRDK2LOErrors[numCalRDK2Points]={0.02,0.02,0.01,0.01,0.01,0.002};
//    TGraphErrors calRDK2Graph(numCalRDK2Points,calRDK2Energies,calRDK2LOs,NULL,calRDK2LOErrors);
//    calRDK2Graph.SetName("CalRDK2Graph");
//    TString calRDK2LegendString="RDK II BGO Array [DAQ]";
//    calRDK2Graph.SetLineColor(kViolet);
//    calRDK2Graph.SetMarkerColor(kViolet);
//    calRDK2Graph.SetMarkerStyle(kFullSquare);
//    calRDK2Graph.SetMarkerSize(1.3);
//
//    //TGraph* khodyukGraph=getKhodyukBGOResponseGraphOfficial();
//    //Assumed error was similar to moszynski and was 2%
//    const int numKhodyukPoints=375;
//    double khodyukEnergies[numKhodyukPoints]={11,11,12,12,13,13.025,13.05,13.075,13.1,13.125,13.15,13.175,13.2,13.225,13.25,13.275,13.3,13.325,13.35,13.375,
//        13.4,13.425,13.45,13.475,13.5,13.525,13.55,13.575,13.6,13.625,13.65,13.675,13.7,13.725,13.75,13.775,13.8,13.825,13.85,13.875,13.9,13.925,13.95,13.975,14,
//        14.025,14.05,14.075,14.1,14.125,14.15,14.175,14.2,14.225,14.25,14.275,14.3,14.325,14.35,14.375,14.4,14.425,14.45,14.475,14.5,14.525,14.55,14.575,14.6,14.625,
//        14.65,14.675,14.7,14.725,14.75,14.775,14.8,14.825,14.85,14.875,14.9,14.925,14.95,14.975,15,15.025,15.05,15.075,15.1,15.125,15.15,15.175,15.2,15.225,15.25,
//        15.275,15.3,15.325,15.35,15.375,15.4,15.425,15.45,15.475,15.5,15.525,15.55,15.575,15.6,15.625,15.65,15.675,15.7,15.725,15.75,15.775,15.8,15.825,15.85,15.875,
//        15.9,15.925,15.95,15.975,16,16.025,16.05,16.075,16.1,16.125,16.15,16.175,16.2,16.225,16.25,16.275,16.3,16.325,16.35,16.375,16.4,16.425,16.45,16.475,16.5,16.525,
//        16.55,16.575,16.6,16.625,16.65,16.675,16.7,16.725,16.75,16.775,16.8,16.825,16.85,16.875,16.9,16.925,16.95,16.975,17,17.025,17.05,17.075,17.1,17.125,17.15,17.175,
//        17.2,17.225,17.25,17.275,17.3,17.325,17.35,17.375,17.4,17.425,17.45,17.475,17.5,17.525,17.55,17.575,17.6,17.625,17.65,17.675,17.7,17.725,17.75,17.775,17.8,17.825,
//        17.85,17.875,17.9,17.925,17.95,17.975,18,18.025,18.05,18.075,18.1,18.125,18.15,18.175,18.2,18.225,18.25,18.275,18.3,18.325,18.35,18.375,18.4,18.425,18.45,18.475,
//        18.5,18.525,18.55,18.575,18.6,18.625,18.65,18.675,18.7,18.725,18.75,18.775,18.8,18.825,18.85,18.875,18.9,18.925,18.95,18.975,19,19.025,19.05,19.075,19.1,19.125,
//        19.15,19.175,19.2,19.225,19.25,19.275,19.3,19.325,19.35,19.375,19.4,19.425,19.45,19.475,19.5,19.525,19.55,19.575,19.6,19.625,19.65,19.675,19.7,19.725,19.75,19.775,
//        19.8,19.825,19.85,19.875,19.9,19.925,19.95,19.975,20,20,25,25,30,30,30,30,35,35,40,40,45,45,45,45,50,50,55,55,60,60,65,65,70,70,75,75,80,80,85,85,89.5,89.6,89.7,
//        89.8,89.9,90,90.1,90.2,90.3,90.4,90.5,90.6,90.7,90.8,90.9,91,91.1,91.2,91.3,91.4,91.5,91.6,91.7,91.8,91.9,92,92.1,92.2,92.3,92.4,92.5,92.6,92.7,92.8,92.9,93,93.1,
//        93.2,93.3,93.4,93.5,93.6,93.7,93.8,93.9,94,94.1,94.2,94.3,94.4,94.5,94.6,94.7,94.8,94.9,95,95,100,100};
//    double khodyukLOs[numKhodyukPoints]={0.700,0.707,0.725,0.724,0.740,0.742,0.737,0.743,0.743,0.741,0.737,0.740,0.738,0.736,0.743,0.746,0.743,0.739,0.742,0.742,
//        0.736,0.722,0.703,0.704,0.705,0.695,0.703,0.701,0.698,0.699,0.702,0.702,0.704,0.695,0.701,0.696,0.705,0.696,0.706,0.701,0.696,0.693,0.697,0.697,0.694,0.700,0.699,0.698,
//        0.699,0.702,0.697,0.692,0.693,0.696,0.697,0.689,0.697,0.683,0.684,0.693,0.687,0.698,0.695,0.702,0.699,0.697,0.697,0.699,0.697,0.699,0.694,0.691,0.699,0.691,0.692,0.694,
//        0.694,0.702,0.697,0.699,0.701,0.695,0.698,0.699,0.702,0.698,0.697,0.700,0.702,0.698,0.697,0.698,0.701,0.708,0.703,0.697,0.699,0.698,0.701,0.705,0.704,0.706,0.702,0.700,
//        0.704,0.704,0.706,0.712,0.701,0.705,0.704,0.704,0.702,0.698,0.697,0.702,0.703,0.698,0.699,0.693,0.701,0.700,0.701,0.700,0.694,0.697,0.698,0.698,0.699,0.695,0.701,0.704,
//        0.698,0.701,0.697,0.701,0.697,0.700,0.694,0.697,0.691,0.687,0.684,0.690,0.692,0.685,0.689,0.687,0.693,0.686,0.694,0.697,0.690,0.691,0.690,0.688,0.695,0.690,0.690,0.690,
//        0.690,0.694,0.693,0.693,0.690,0.688,0.692,0.687,0.689,0.692,0.687,0.694,0.695,0.692,0.693,0.693,0.698,0.694,0.695,0.692,0.699,0.694,0.693,0.695,0.695,0.696,0.692,0.694,
//        0.690,0.694,0.697,0.698,0.700,0.700,0.700,0.696,0.701,0.696,0.696,0.698,0.704,0.705,0.696,0.694,0.698,0.701,0.702,0.703,0.697,0.701,0.702,0.701,0.698,0.708,0.702,0.702,
//        0.701,0.702,0.699,0.706,0.705,0.700,0.714,0.699,0.710,0.705,0.704,0.702,0.704,0.703,0.703,0.703,0.705,0.706,0.703,0.701,0.712,0.701,0.710,0.707,0.709,0.704,0.709,0.713,
//        0.706,0.711,0.710,0.709,0.719,0.709,0.705,0.711,0.712,0.711,0.708,0.709,0.709,0.710,0.710,0.716,0.716,0.712,0.709,0.717,0.713,0.715,0.711,0.716,0.716,0.717,0.714,0.715,
//        0.711,0.720,0.714,0.718,0.715,0.711,0.720,0.716,0.717,0.722,0.722,0.724,0.720,0.718,0.768,0.770,0.803,0.800,0.800,0.800,0.827,0.826,0.844,0.845,0.861,0.864,0.858,0.858,
//        0.875,0.874,0.880,0.880,0.885,0.885,0.890,0.889,0.896,0.893,0.897,0.898,0.901,0.895,0.906,0.909,0.911,0.911,0.912,0.910,0.911,0.910,0.909,0.910,0.908,0.908,0.899,0.885,
//        0.885,0.884,0.885,0.881,0.883,0.880,0.882,0.881,0.879,0.879,0.878,0.881,0.879,0.879,0.880,0.880,0.882,0.882,0.877,0.878,0.878,0.880,0.876,0.879,0.880,0.881,0.879,0.878,
//        0.878,0.879,0.878,0.880,0.880,0.879,0.878,0.878,0.877,0.878,0.880,0.877,0.878,0.878,0.878,0.881,0.880,0.882,0.883};
//    TGraph khodyukGraph(numKhodyukPoints,khodyukEnergies,khodyukLOs);
//    TString khodyukLegendString="I. Khodyuk and P. Dorenbos";
//    khodyukGraph.SetLineColor(kBlue);
//    khodyukGraph.SetMarkerColor(kBlue);
//    khodyukGraph.SetMarkerStyle(kFullCircle);
//
//
//    const int numMoszynskiPoints=30;
//    double moszynskiEnergies[numMoszynskiPoints]={14.4,16.6,22,31,32,45,52,59.6,74,81,84,88,122,244.1,302.9,320.1,335.8,356,411.8,435,493.7,511,569.7,586,662,758.8,834.8,1063.7,1173.2,1275};
//    double moszynskiLOs[numMoszynskiPoints]={0.711,0.719,0.775,0.829,0.8412,0.853,0.9,0.893,0.91,0.926,0.937,0.933,0.915,0.964,0.979,0.964,0.97,0.986,0.97,0.995,1,0.996,1,0.999,1,1,1,1.01,1.013,1.01};
//    double moszynskiLOErrors[numMoszynskiPoints]={0.01422,0.01438,0.0155,0.01658,0.016824,0.02559,0.018,0.01786,0.0182,0.01852,0.01874,0.01866,0.0183,0.02892,0.01958,0.01928,0.0291,0.01972,0.0194,0.02985,0.03,0.01992,0.02,0.02997,0.02,0.03,0.02,0.0202,0.02026,0.0202};
//    TGraphErrors moszynskiGraph(numMoszynskiPoints,moszynskiEnergies,moszynskiLOs,NULL,moszynskiLOErrors);
//    moszynskiGraph.SetName("MoszynskiGraph");
//    TString moszynskiLegendString="Moszynski et. al.";
//    moszynskiGraph.SetLineColor(kGreen);
//    moszynskiGraph.SetMarkerColor(kGreen);
//    moszynskiGraph.SetMarkerStyle(kFullTriangleUp);
//    moszynskiGraph.SetMarkerSize(1.3);
//
////    const int numAverkievPoints=20;
////    double averkievEnergies[numAverkievPoints]={2.1,6.4,7.5,8.0,8.8,10.0,11.2,16.6,17.2,23.5,33.1,58.3,71.6,78.3,120.1,162.3,277.4,384.8,506.5,657.9};
////    double averkievLOs[numAverkievPoints]={0.854,0.777,0.739,0.706,0.678,0.657,0.641,0.579,0.590,0.685,0.731,0.837,0.844,0.827,0.884,0.903,0.959,0.977,0.986,0.996};
////    double averkievLOErrors[numAverkievPoints]={0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.035,0.035,0.035,0.035,0.035,0.02,0.015,0.01,0.005,0.005};
////    TGraphErrors averkievGraph(numAverkievPoints,averkievEnergies,averkievLOs,NULL,averkievLOErrors);
////    averkievGraph.SetName("AverkievGraph");
////    TString averkievLegendString="Averkiev et. al.";
////    averkievGraph.SetLineColor(kMagenta);
////    averkievGraph.SetMarkerColor(kMagenta);
////    averkievGraph.SetMarkerStyle(kOpenCircle);
////    averkievGraph.SetMarkerSize(1.5);
//
//
//    double legendSize=.05*numGraphs;
//    TLegend* theLegend = new TLegend(0.25,.15,0.9,0.5);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//
////    theLegend->AddEntry(&quadraticGraph,quadraticLegendString,"L");
//    theLegend->AddEntry(&khodyukGraph,khodyukLegendString,"L");
//    theLegend->AddEntry(&moszynskiGraph,moszynskiLegendString,"P");
////    theLegend->AddEntry(&averkievGraph,averkievLegendString,"P");
//    theLegend->AddEntry(&dewarGraph,dewarLegendString,"P");
//    theLegend->AddEntry(&calRDK2Graph,calRDK2LegendString,"P");
//
//
//    TCanvas* theCanvas=getHalfPresentationPlot();
//
////    khodyukGraph.GetYaxis()->SetTitleOffset(1.3);
//    khodyukGraph.GetXaxis()->SetLimits(0,120);
//    khodyukGraph.GetYaxis()->SetRangeUser(.65,1.05);
//    khodyukGraph.GetYaxis()->SetTitleOffset(1.8);
//    khodyukGraph.SetTitle(";Energy Absorbed in BGO (keV);Light Output Relative to 662 keV");
//
//    khodyukGraph.Draw("A L");
////    quadraticGraph.Draw("L same");
//    moszynskiGraph.Draw("P E1 same");
////    averkievGraph.Draw("P E1 same");
//    dewarGraph.Draw("P E1 same");
//    calRDK2Graph.Draw("P E1 same");
//
//    theCanvas->SaveAs(GraphsDir+"LOExperimentsGraphZoom_140713.png");
//    theLegend->Draw("same");
//    khodyukGraph.GetXaxis()->SetLimits(0,1000);
//    theCanvas->SaveAs(GraphsDir+"LOExperimentsGraph_140713.png");

//Cal Hist
//TH1D* calHist= getTabSeperatedHist(HistsDir+"Cs137Co57_S502_Det1.txt",2018,-.5,2017.5);
//    calHist->SetTitle(";Uncalibrated BGO Peak Height Channel;Counts");
//    calHist->Rebin(4);
//    calHist->SetMarkerColor(kRed);
//    calHist->SetLineColor(kRed);
//    TCanvas* theCanvas = getHalfPresentationPlot2();
//    calHist->Draw("P E1");
//    theCanvas->SaveAs(GraphsDir+"Cs137Co57140610.png");

///Average PLO Hist

//    TH1D* ploHist=getTabSeperatedHist(HistsDir+"averagePLOHist140610.txt",101,-100,100);
//    ploHist->SetTitle(";Distance from Center of Crystal (mm);Relative Response to Crystal Center");
//
//    TCanvas* theCanvas = getThesisPlotCanvas();
//    ploHist->SetLineColor(kRed);
//    ploHist->GetSumw2()->Set(0);
//
//    ploHist->Draw("L");
//    theCanvas->SaveAs(GraphsDir+"averagePLOHist140610.root");

///z0 Hist for thesis (check to see what no ToF cuts looks like)
//    HistDim z0Dim={280,-.4, .72};
//
//    RDK2MCAnalysis theAnalysis("AID36",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
//
//
//    TH1* z0EPHist=theAnalysis.MakeHist("zoHist",THREEBODY,"z0","pDChn.SBDEDepBlur > 7 && pDChn.SBDEDepBlur < 31 && eDChn.SBDEDepBlur > 10 && eDChn.SBDEDepBlur < 800",z0Dim);
//    convertTH1ToTXT(z0EPHist,HistsDir+"z0EPHist140723NoToF.txt");
//
////    TH1* z0EPHist=getTabSeperatedHist(HistsDir+"z0EPHist140723NoToF.txt",z0Dim);
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
//
//
//    for (int i = 1;i<=z0Dim.numBins ;i++)
//    {
//        z0EPHist->SetBinError(i,0);
//
//
//    }
//
//    z0EPHist->Rebin(2);
//
//
//    scaleHistTo(z0EPHist,.5*50);
//
//
//
//    z0EPHist->SetTitle(";Decay Position in Neutron Beam (m);Detection Probability (AU)");
//
//
//    z0EPHist->SetLineColor(kBlue);
//
//
//    z0EPHist->GetXaxis()->SetRangeUser(-.39,.45);
////    z0EPHist->GetXaxis()->SetLabelSize(.08);
////    z0EPHist->GetXaxis()->SetTitleSize(.08);
//
//    z0EPHist->GetYaxis()->SetRangeUser(0,.020*50);
//    z0EPHist->GetYaxis()->SetNdivisions(508);
//
//    z0EPHist->Draw("L");
//
//    TLegend legend(0.11,0.6,0.42,0.9,NULL,"brNDC");
////    legend.SetTextSize(1.4);
//    legend.SetFillColor(kWhite);
//    legend.AddEntry(z0EPHist,"EP");
//
//    theCanvas->SaveAs(GraphsDir+"z0Plot140723NoToF.pdf");
//
//    delete theCanvas;

///z0 Hist for thesis  (fix left scale and other improvements) 140623
//   HistDim z0Dim={280,-.4, .72};
//
////    RDK2MCAnalysis theAnalysis("AID36",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
////    TH1* z0EPHist=theAnalysis.MakeStandardHist("epHist",DET_EP,"z0",z0Dim);
////    TH1* z0EPGHist=theAnalysis.MakeStandardHist("epgHist",DET_EPG,"z0",z0Dim);
////    TH1* z0EPBGHist=theAnalysis.MakeStandardHist("epbgHist",DET_EPBG,"z0",z0Dim);
////
////    convertTH1ToTXT(z0EPHist,HistsDir+"z0EPHist140511.txt");
////    convertTH1ToTXT(z0EPGHist,HistsDir+"z0EPGHist140511.txt");
////    convertTH1ToTXT(z0EPBGHist,HistsDir+"z0EPBGHist140511.txt");
//
//    TCanvas* theCanvas=new TCanvas("thesisPlot","thesisPlot",10,10,2400,1000);
//    theCanvas->SetLeftMargin(.085);
//    theCanvas->SetRightMargin(.01);
//    theCanvas->SetTopMargin(.05);
//    theCanvas->SetBottomMargin(.18);
//    gPad->SetTickx(1);
//    gPad->SetTicky(1);
//    gPad->SetFillColor(kWhite);
//    TH1* z0EPHist=getTabSeperatedHist(HistsDir+"z0EPHist140518.txt",z0Dim);
//    TH1* z0EPGHist=getTabSeperatedHist(HistsDir+"z0EPGHist140518.txt",z0Dim);
//    TH1* z0EPBGHist=getTabSeperatedHist(HistsDir+"z0EPBGHist140518.txt",z0Dim);
//
//    for (int i = 1;i<=z0Dim.numBins ;i++)
//    {
//        z0EPHist->SetBinError(i,0);
//        z0EPGHist->SetBinError(i,0);
//        z0EPBGHist->SetBinError(i,0);
//
//    }
//
//    z0EPHist->Rebin(2);
//    z0EPGHist->Rebin(2);
//    z0EPBGHist->Rebin(2);
//
//    scaleHistTo(z0EPHist,1.*1000);
//    scaleHistTo(z0EPGHist,.3*1000);
//    scaleHistTo(z0EPBGHist,.15*1000);
//
//
//    z0EPHist->SetTitle(";Decay Position in Neutron Beam (m);Detection Prob. Per Unit Length (AU)");
////    z0EPHist->SetTitle(";;");
//
//    z0EPHist->SetLineColor(kBlue);
//    z0EPGHist->SetLineColor(kRed);
//    z0EPBGHist->SetLineColor(kGreen);
//
//    z0EPHist->SetLineWidth(2);
//    z0EPGHist->SetLineWidth(2);
//    z0EPBGHist->SetLineWidth(2);
//
//    z0EPGHist->SetLineStyle(2);
//    z0EPBGHist->SetLineStyle(3);
//
//    z0EPHist->GetXaxis()->SetRangeUser(-.39,.45);
//    z0EPHist->GetXaxis()->SetLabelSize(.08);
//    z0EPHist->GetXaxis()->SetTitleSize(.07);
//
//    z0EPHist->GetYaxis()->SetRangeUser(0,.020*1000);
//    z0EPHist->GetYaxis()->SetTitleOffset(.6);
//    z0EPHist->GetYaxis()->SetLabelSize(.08);
//    z0EPHist->GetYaxis()->SetTitleSize(.055);
//    z0EPHist->GetYaxis()->SetNdivisions(508);
//
//
//
//    z0EPHist->Draw("L");
//    z0EPGHist->Draw("L same");
//    z0EPBGHist->Draw("L same");
//    z0EPHist->Draw("L same");
//
//    TLegend legend(0.11,0.4,0.42,0.9,NULL,"brNDC");
////    legend.SetTextSize(1.4);
//    legend.SetFillColor(kWhite);
//    legend.AddEntry(z0EPHist,"ep");
//    legend.AddEntry(z0EPGHist,"ep#gamma (BGO)");
//    legend.AddEntry(z0EPBGHist,"ep#gamma (LAAPD)");
//    legend.Draw();
//
//    theCanvas->SaveAs(GraphsDir+"z0Plot140723.pdf");
//
//    delete theCanvas;

///Make BAPD Plots (fix things) 140623
//    const HistDim BAPD_DIMb={38,0.5,10};
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"EXPBen_Hist131009_BAPD_AverageIndiv_EPGb.txt",BAPD_DIMb);
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248b.txt",BAPD_DIMb);
//
//    scaleHistogramsTogether(expHist,mcHist);
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
////    TCanvas* theCanvas=getHalfPresentationPlot2();
//
//    TH1F* hr=theCanvas->DrawFrame(0,0,10,5.5e-6);
//    hr->SetXTitle("Photon Peak Height (5.9=5.9 keV)");
//    hr->SetYTitle("#frac{LAAPD ep#gamma counts}{ep counts}");
//    expHist->SetMarkerColor(kBlack);
//    expHist->SetLineColor(kBlack);
//    mcHist->SetLineColor(kGreen+2);
//    mcHist->GetSumw2()->Set(0);
//    expHist->SetMarkerStyle(4);
//    mcHist->SetLineStyle(3);
//    mcHist->SetLineWidth(2);
//
//    mcHist->Rebin(2);
//    expHist->Rebin(2);
//
//    gStyle->SetErrorX(0.0001);
//
//    mcHist->SetTitle(";Photon Peak Height (5.9=5.9 keV);#frac{LAAPD ep#gamma counts}{ep counts}");
//
//
//    mcHist->GetYaxis()->SetRangeUser(0,5.5e-6);
//    mcHist->GetXaxis()->SetRangeUser(0,10);
//
//    mcHist->Draw("L same");
//    expHist->Draw("P E0 Same");
//
//    TLegend* theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(expHist,"Experiment","PE");
//    theLegend->AddEntry(mcHist,"Simulation","L");
//    theLegend->Draw("same");
//
//     TLine line(0,0,10,0);
////     line.Draw("same");
//
////    theCanvas->SaveAs(GraphsDir+"compPlot_EPBG_BGE_140723b.png");
//    theCanvas->SaveAs(GraphsDir+"compPlot_EPBG_BGE_140723.pdf");
//
//    hr = theCanvas->DrawFrame(0,-1.5,10,1.5);
//    hr->SetXTitle("Photon Peak Height (5.9=5.9 keV)");
//    hr->SetYTitle("Normalized Residual #frac{Exp - Sim}{Sim}");
//    TH1D* ratioHist = (TH1D*) expHist->Clone("ratioHist");
//    ratioHist->Reset();
//    ratioHist->Sumw2();
//    ratioHist->SetTitle(";Photon Peak Height (5.9=5.9 keV);Normalized Residual #frac{Exp - Sim}{Sim}");
//
//    ratioHist->Add(expHist,mcHist,1,-1);
//    ratioHist->Divide(mcHist);
////    ratioHist->Divide(expHist,mcHist);
////    ratioHist->Add(new TF1("constantn1","-1",0,20)); //subtract one for (exp-mc)/mc
//
//    ratioHist->SetLineColor(kGreen+2);
//    ratioHist->SetMarkerColor(kGreen+2);
//    ratioHist->GetYaxis()->SetRangeUser(-1.5,1.5);
////    ratioHist->GetXaxis()->SetRangeUser(0,10);
//
//    ratioHist->Draw("P E0 same");
//
//    line.Draw("same");
//    ratioHist->Draw("P E0 same");
//
////    theCanvas->SaveAs(GraphsDir+"compPlot_EPBG_BGE_NResid_140723b.png");
////    theCanvas->SaveAs(GraphsDir+"compPlot_EPBG_BGE_NResid_140723.pdf");

///Make Exp/Sim Plots 140723
////    RDK2AnalysisPlotter* thePlotter;
////
////    thePlotter=new RDK2AnalysisPlotter("Latest and Greatest","EID3","AID48");
////    TH1* expHist = thePlotter->getExpHist(0,DET_EPG, PLOTVAR_GE);
////    TH1* mcHist = thePlotter->getMCHist(0,DET_EPG, PLOTVAR_GE);
////
////    convertTH1ToTXT(expHist,HistsDir+"GEHistExp.txt");
////    convertTH1ToTXT(mcHist,HistsDir+"GEHistMC.txt");
//
//    RDK2AnalysisPlotter* thePlotter;
//    TH1* expHist;
//    TH1* mcHist;
//
//    double xMin=0,xMax=1;
//
//
//    thePlotter=new RDK2AnalysisPlotter("Latest and Greatest","EID3","AID48");
//
//    CoDet detType;
//    PlotVarType varType;
//
//
//    for (int k = 0;k< 3;k++)
//    {
//
//        TCanvas* theCanvas;
//        if(k==0)
//        {
//            theCanvas=getThesisPlotCanvas();
//        }
//        else if(k==1)
//        {
//            theCanvas = getHalfPresentationPlot2();
//        }
//        else
//        {
//            theCanvas = getHalfPresentationPlot();
//        }
//
//        TH1* tempHist;
//        for (int i = 1;i< 2;i++)
//        {
//            detType=static_cast<CoDet>(i+1);
//            for (int j = 0;j< 1;j++)
//            {
//                varType=static_cast<PlotVarType>(j);
//
//                if(detType==DET_EP && varType==PLOTVAR_GE)
//                    continue;
//
//
//                if(j<3)
//                {
//                    expHist = thePlotter->getExpHist(0,detType, varType);
//                    mcHist = thePlotter->getMCHist(0,detType, varType);
//                }
//                else
//                {
//                    expHist=getTabSeperatedHist(HistsDir+"GEHistExp.txt",780,10,780);
//                    mcHist=getTabSeperatedHist(HistsDir+"GEHistMC.txt",780,10,780);
//                }
//
//                Color_t theColor;
//                TString yAxisString;
//                TString imageString="compPlot_";
//                if(detType==DET_EP)
//                {
//                    theColor=kBlue;
//                    yAxisString="Counts (AU)";
//                    imageString+="EP_";
//                }
//                else if(detType==DET_EPG)
//                {
//                    theColor=kRed;
//                    yAxisString="#frac{BGO ep#gamma counts}{ep counts}";
//                    mcHist->SetLineStyle(2);
//                    mcHist->SetLineWidth(2);
//                    imageString+="EPG_";
//                }
//                else if(detType==DET_EPBG)
//                {
//                    theColor=kGreen+2;
//                    yAxisString="#frac{LAAPD ep#gamma counts}{ep counts}";
//    //                mcHist->SetLineStyle(3);
//                    mcHist->SetLineWidth(2);
//                    imageString+="EPBG_";
//                }
//
//                TString xAxisString;
//                int rebin;
//                if(varType==PLOTVAR_EE)
//                {
//                    xAxisString="Electron Kinetic Energy Detected (keV)";
//                    imageString+="EE";
//                    rebin=10;
//
//                    if(detType!=DET_EP)
//                    {
//                        rebin=20;
//                    }
//                    xMax=800;
//                    xMin=100;
//
//
//                }
//                else if(varType==PLOTVAR_PE)
//                {
//                    xAxisString="Proton Kinetic Energy Detected (keV)";
//                    imageString+="PE";
//                    rebin=2;
//                    xMax=31;
//                    xMin=7;
//                }
//                else if(varType==PLOTVAR_PT)
//                {
//                    xAxisString="Proton Time-of-Flight (s)";
//                    imageString+="PT";
//                    xMax=10e-6;
//                    xMin=2e-6;
//
//                    rebin=9;
//                }
//                else if(varType==PLOTVAR_GE)
//                {
//                    imageString+="GE";
//                    if(detType==DET_EPG)
//                    {
//                       xAxisString="BGO Photon Peak Height (511=511 keV)";
//                       xMax=800;
//                        xMin=0;
//                    }
//                    else if(detType==DET_EPBG)
//                    {
//                        xAxisString="LAAPD Photon Peak Height (5.6=5.6 keV)";
//                    }
//                    rebin=10;
////                    tempHist=expHist;  expHist=trimHist(expHist,expHist->GetName()+TString("trim"),10); delete tempHist;
////                    tempHist=mcHist;  mcHist=trimHist(mcHist,mcHist->GetName()+TString("trim"),10); delete tempHist;
//
//                }
//
//                expHist->SetMarkerColor(kBlack);
//                expHist->SetLineColor(kBlack);
//                mcHist->SetLineColor(theColor);
//                mcHist->GetSumw2()->Set(0);
//                expHist->SetMarkerStyle(4);
//                mcHist->GetYaxis()->SetTitleOffset(1.7);
//    //            theCanvas->SetLeftMargin(0.2);
//
//                mcHist->Rebin(rebin);
//                expHist->Rebin(rebin);
//
//                gStyle->SetErrorX(0.0001);
//
//                mcHist->SetTitle(";"+xAxisString+";"+yAxisString);
//
//                scaleHistogramsTogether(expHist,mcHist);
//
//                TH1F* hr=theCanvas->DrawFrame(xMin,0,xMax,expHist->GetMaximum()*1.15);
//                hr->SetXTitle(xAxisString);
//                hr->SetYTitle(yAxisString);
//
//                mcHist->Draw("L Same");
//                expHist->Draw("P E0 Same");
//
//                TLegend* theLegend = new TLegend(0.70,.55,0.90,0.9);
//                theLegend->SetBorderSize(1);
//                theLegend->SetFillColor(kWhite);
//                theLegend->AddEntry(expHist,"Experiment","PE");
//                theLegend->AddEntry(mcHist,"Simulation","L");
//                theLegend->Draw("same");
//
//                imageString+="_140723";
//                if(k==0)
//                {
//                    imageString += ".pdf";
//                }
//                else if(k==1)
//                {
//                    imageString+="b.png";
//                }
//                else if(k==2)
//                {
//                    imageString+="c.png";
//                }
//
//                cout << GraphsDir+imageString << endl;
//                theCanvas->SaveAs(GraphsDir+imageString);
//                delete mcHist;
//                delete expHist;
//                delete theLegend;
//                delete hr;
//            }
//        }
//        delete theCanvas;
//    }

///BGO Var Hists
//    TH1D* protoHist=new TH1D("epg_gEAvgVarHist","epg_gEAvgHist",20,BGO_VARHIST_LOWEDGES);
//    TH1D* expHist=getTabSeperatedHist(HistsDir+"VarHist_MC_AID48_EPG_EDepGAvg.txt",protoHist);
//    TH1D* mcHist=getTabSeperatedHist(HistsDir+"VarHist_Exp_EID3_EPG_EDepGAvg.txt",protoHist);
////    expHist->SetBinContent(2,0);
////    mcHist->SetBinContent(2,0);
//
//    scaleHistogramsTogether(expHist,mcHist);
//
//
//
//
//    TCanvas* theCanvas=getThesisPlotCanvas();
////    TCanvas* theCanvas=getHalfPresentationPlot2();
//    expHist->SetMarkerColor(kBlack);
//    expHist->SetLineColor(kBlack);
//    mcHist->SetLineColor(kGreen+2);
////    mcHist->GetSumw2()->Set(0);
//    expHist->SetMarkerStyle(4);
//
////    mcHist->Rebin(2);
////    expHist->Rebin(2);
//
//    gStyle->SetErrorX(0.0001);
//
//    mcHist->SetTitle(";Calibrated Photon Peak Height (5.6=5.6 keV);#frac{LAAPD ep#gamma counts}{ep counts}");
//
//
////    mcHist->GetYaxis()->SetRangeUser(0,5.5e-6);
//
//
//
//    TH1D* ratioHist = (TH1D*) expHist->Clone("ratioHist");
//    ratioHist->Reset();
//    ratioHist->Sumw2();
//    ratioHist->SetTitle(";Photon Peak Height (511=511 keV);Normalized Residual #frac{Exp - Sim}{Sim}");
//
//    ratioHist->Add(expHist,mcHist,1,-1);
//    ratioHist->Divide(mcHist);
////    ratioHist->Divide(expHist,mcHist);
////    ratioHist->Add(new TF1("constantn1","-1",0,20)); //subtract one for (exp-mc)/mc
//
//    ratioHist->SetLineColor(kRed);
//    ratioHist->SetMarkerColor(kRed);
////    ratioHist->SetLineStyle(2);
//    ratioHist->SetLineWidth(2);
//    ratioHist->GetYaxis()->SetRangeUser(-0.3,.3);
//
//    ratioHist->Draw("P E0");
//    theCanvas->SetLogx();
//    TLine line(0,0,800,0);
//    line.Draw("same");
//    ratioHist->Draw("P E0 same");
//
//    theCanvas->SaveAs(GraphsDir+"compPlot_EPG_GE_NResid_140723.pdf");
////    theCanvas->SaveAs(GraphsDir+"compPlot_EPG_GE_NResid_140701b.png");

///initial event plots
///Make hists
//    HistDim ep0Dim={100,0.,1};
//    HistDim bgoDim={790,10.,800};
//    int numEventFiles=48;
//    RDK2MCAnalysis theAnalysis("AID48",  "20 nm Si Dead Layer on SBD with P Resolution", numEventFiles,  "RID358", "RID359",  numEventFiles, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
////
////    TH1* ep03BHist=theAnalysis.MakeHist("ep0Hist",THREEBODY,"ep0","",ep0Dim);
////    ep03BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ep03BHist,HistsDir+"ep0Hist3B140610.txt");
////
////    TH1* ep04BHist=theAnalysis.MakeHist("ep0Hist2",FOURBODY,"ep0","eg0 > 14",ep0Dim);
////    ep04BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ep04BHist,HistsDir+"ep0Hist4B140610.txt");
////
////    TH1* ee03BHist=theAnalysis.MakeHist("ee0Hist",THREEBODY,"ee0","",EDEPE_DIM);
////    ee03BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ee03BHist,HistsDir+"ee0Hist3B140610.txt");
////
////    TH1* ee04BHist=theAnalysis.MakeHist("ee0Hist2",FOURBODY,"ee0","eg0 > 14",EDEPE_DIM);
////    ee04BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(ee04BHist,HistsDir+"ee0Hist4B140610.txt");
////
////    TH1* eg04BHist=theAnalysis.MakeHist("eg0Hist",FOURBODY,"eg0","eg0 > 10",bgoDim);
////    eg04BHist->Scale(1./(numEventFiles*1.e6));
////    convertTH1ToTXT(eg04BHist,HistsDir+"eg0Hist4B140610.txt");
//
//    TH1D* ep03BHist=getTabSeperatedHist(HistsDir+"ep0Hist3B140610.txt",ep0Dim);
//    TH1D* ep04BHist=getTabSeperatedHist(HistsDir+"ep0Hist4B140610.txt",ep0Dim);
//    TH1D* ee03BHist=getTabSeperatedHist(HistsDir+"ee0Hist3B140610.txt",EDEPE_DIM);
//    TH1D* ee04BHist=getTabSeperatedHist(HistsDir+"ee0Hist4B140610.txt",EDEPE_DIM);
//    TH1D* eg04BHist=getTabSeperatedHist(HistsDir+"eg0Hist4B140610.txt",bgoDim);
//
//
//    ee03BHist->Rebin(10);
//    ee04BHist->Rebin(10);
//    eg04BHist->Rebin(10);
//    ee03BHist->Scale(10);
//    ep03BHist->Scale(10);
//
//    scaleHistogramsTogether(ep03BHist,ep04BHist);
//    scaleHistogramsTogether(ee03BHist,ee04BHist);
//
//    TCanvas* theCanvas = getThesisPlotCanvas();
////    ep0Hist->Smooth(10);
//    ep03BHist->SetLineColor(kBlue);
//    ep04BHist->SetLineColor(kRed);
//    ep04BHist->SetLineStyle(2);
//    ep04BHist->SetLineWidth(2);
//
//    ep03BHist->SetTitle(";Initial Kinetic Energy from Decay (keV);Emission Probability (AU)");
//    ep03BHist->GetSumw2()->Set(0);
//    ep04BHist->GetSumw2()->Set(0);
////    ep0Hist->GetXaxis()->SetRangeUser(-1,1);
//    ep03BHist->Draw("L");
//    ep04BHist->Draw("L same");
//
//    TLegend* theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(ep03BHist,"Three Body Decays","L");
//    theLegend->AddEntry(ep04BHist,"Four Body Decays","L");
//    theLegend->Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"ep0Hist140610.pdf");
//    delete theLegend;
//
//    ee03BHist->SetLineColor(kBlue);
//    ee04BHist->SetLineColor(kRed);
//    ee04BHist->SetLineStyle(2);
//    ee04BHist->SetLineWidth(2);
//
//    ee03BHist->SetTitle(";Initial Kinetic Energy from Decay (keV);Emission Probability (AU)");
//    ee03BHist->GetSumw2()->Set(0);
//    ee04BHist->GetSumw2()->Set(0);
//    ee03BHist->GetYaxis()->SetRangeUser(0,0.25);
//    ee03BHist->Draw("L");
//    ee04BHist->Draw("L same");
//
//    theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(ee03BHist,"Three Body Decays","L");
//    theLegend->AddEntry(ee04BHist,"Four Body Decays","L");
//    theLegend->Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"ee0Hist140610.pdf");
//    delete theLegend;
//
//    eg04BHist->SetLineColor(kRed);
//    eg04BHist->SetLineStyle(2);
//    eg04BHist->SetLineWidth(2);
//
//    eg04BHist->SetTitle(";Photon Energy from Decay (keV);Emission Probability (AU)");
//    eg04BHist->GetSumw2()->Set(0);
//
//    eg04BHist->Draw("L");
//
//
//    theLegend = new TLegend(0.70,.55,0.90,0.9);
//    theLegend->SetBorderSize(1);
//    theLegend->SetFillColor(kWhite);
//    theLegend->AddEntry(eg04BHist,"Four Body Decays","L");
//    theLegend->Draw("same");
//
//    theCanvas->SaveAs(GraphsDir+"eg0Hist140610.pdf");
//    delete theLegend;

///Co-plot new bAPD
////    createMCAnalysisAndPlot( "AID48",  "V1.036 Correct SBD Pos", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetB);
//    TH1* expHist[2]={getTabSeperatedHist(HistsDir+"ExpHist_Jeff_EPBG_EDepBGAvg.txt",BAPD_DIM),getTabSeperatedHist(HistsDir+"EXPBen_Hist131009_BAPD_AverageIndiv_EPG.txt",BAPD_DIM)};
//    expHist[0]->Sumw2();
//    expHist[0]->Scale(3);
//    TH1* mcHist[2]={getTabSeperatedHist(HistsDir+"Hist_MC_AID48_EPBG_EDepBGAvg.txt",BAPD_DIM),getTabSeperatedHist(HistsDir+"MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248.txt",BAPD_DIM)};;
//
//
//    for (int i = 0;i< 2;i++)
//    {
//        for (int j = 0;j< 2;j++)
//        {
//            expHist[j]->SetBinContent(i+1,0);
//            expHist[j]->SetBinError(i+1,0);
//            mcHist[j]->SetBinContent(i+1,0);
//            mcHist[j]->SetBinError(i+1,0);
//
//        }
//
//    }
//    double error;
//    cout << "Exp New: " << expHist[0]->IntegralAndError(1,80,error);cout << " +/- " << error  << endl;
//    cout << "Exp Old: " << expHist[1]->IntegralAndError(1,80,error);cout << " +/- " << error  << endl;
//    cout << "MC New: " << mcHist[0]->IntegralAndError(1,80,error);cout << " +/- " << error  << endl;
//    cout << "MC Old: " << mcHist[1]->IntegralAndError(1,80,error);cout << " +/- " << error  << endl;
//
//
//
//    expHist[0]->SetTitle("Jeff (New);bAPD Peak Height (5.9=5.9 keV);#frac{ep#gamma}{ep}");
//    expHist[1]->SetTitle("Ben (Old);bAPD Peak Height (5.9=5.9 keV);#frac{ep#gamma}{ep}");
//    mcHist[0]->SetTitle("MC (New);bAPD Peak Height (5.9=5.9 keV);#frac{ep#gamma}{ep}");
//    mcHist[1]->SetTitle("MC (Old);bAPD Peak Height (5.9=5.9 keV);#frac{ep#gamma}{ep}");
//
//
//    plotExpVersusMCToImage(1,expHist,1,mcHist,"Exp vs. Sim: bAPD","bAPDExpVsSim140730.png");
//    plotResidualsToImage(1,expHist,1,mcHist,"Exp vs. Sim: bAPD  Norm. Residuals","bAPDExpVsSim_Resid140730.png");

	/*Make Multiplicity Hist - 141005*/
//    RDK2MCAnalysis theAnalysis("AID48",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
//    TH1* multHist=theAnalysis.MakeMultiplicityHist();
//    convertTH1ToTXT(multHist,HistsDir+"MultHist_AID48_BGO.txt");
//    TH1* multHist = getTabSeperatedHist(HistsDir+"MultHist_AID48_BGO.txt",13,-.5,12.5);
//
//    TCanvas* theCanvas = getThesisPlotCanvas();
//    theCanvas->SetLogy();
//    multHist->Sumw2();
//    multHist->Scale(1./multHist->Integral(2,13));
//    multHist->SetTitle(";Multiplicity of Radiative Events in MC;Relative Number of Events");
//    multHist->GetXaxis()->SetRangeUser(.5,5.5);
//
//
//    multHist->Draw("colz");
//    theCanvas->SaveAs(GraphsDir+"MultHist_AID48_BGO.png");
//    delete multHist;
//    delete theCanvas;
	/*Make Multiplicity 2 Hist 141015*/
//    RDK2MCAnalysis theAnalysis("AID48",  "20 nm Si Dead Layer on SBD with P Resolution", 48,  "RID358", "RID359",  48, "RID360",  "RID361", "RID324",  STD_EPCutSetB);
//    TH1* mult2Hist=theAnalysis.MakeRadAvgHist("",DET_EPG,"",1,BGO_DIM);
//    convertTH1ToTXT(mult2Hist,HistsDir+"Mult2Hist_AID48_BGO141015.txt");
//    TH1* mcHist = getTabSeperatedHist(HistsDir+"Mult2Hist_AID48_BGO141015.txt",BGO_DIM);
//    TH1* expHist = getTabSeperatedHist(HistsDir+"Jeff_BGO_Mult2_141006.txt",200,0,1000);
//
//    for(int i=0;i<10;i++)
//    	mcHist->SetBinContent(i+1,0);
//
////    TCanvas* theCanvas = getThesisPlotCanvas();
//    mcHist->Rebin(5);
//    mcHist->SetTitle("MC;Peak Height (511=511 keV);Scaled Counts");
//    expHist->SetTitle("Exp;Peak Height (511=511 keV);Scaled Counts");
//
//    scaleHistogramsTogether(expHist,mcHist);
//
//    TH1* mcHists[1]={mcHist};
//    TH1* expHists[1]={expHist};
//
//    mcHist->Rebin(2);
//    expHist->Rebin(2);
//
//    plotExpVersusMCToImage(1, expHists,1,mcHists,"Multiplicity 2","Mult2Hist_AID48_BGO141015.png");
//
//    double more, less, ratio;
//    double moreError, lessError, ratioError;
//
//    more=expHist->IntegralAndError(11,100,moreError);
//    less=expHist->IntegralAndError(1,10,lessError);
//    ratio=more/less;
//    ratioError=getBasicErrorProgation(more,less,moreError,lessError,DIVISION);
//    cout << "Exp Ratio of >100 keV to <100 keV: " << ratio << " +/- " << ratioError<< endl;
//
//
//    more=mcHist->IntegralAndError(11,100,moreError);
//	less=mcHist->IntegralAndError(1,10,lessError);
//	ratio=more/less;
//	ratioError=getBasicErrorProgation(more,less,moreError,lessError,DIVISION);
//	cout << "MC Ratio of >100 keV to <100 keV: " << ratio << " +/- " << ratioError<< endl;
//
////    mcHist->SetLineColor(kRed);
////
////    mcHist->Draw("");
////    //expHist->Draw("P E same");
////	mcHist->SaveAs(GraphsDir+"Mult2Hist_AID48_BGO141015.png");
////	delete theCanvas;
//	delete mcHist;
//	delete expHist;
	/*Test Events ranges*/

//    cEvents theEvents;
//    double br100eV = theEvents.calcBR(0.1,EKEMAX,0,10000000,true,880.);
//    double br1keV = theEvents.calcBR(1,EKEMAX,0,10000000,true,880.);
//    double br10keV = theEvents.calcBR(10,EKEMAX,0,10000000,true,880.);
//
//    cout << "br100eV: " << br100eV << endl;
//    cout << "br1keV: " << br1keV << endl;
//    cout << "br10keV: " << br10keV << endl;
	/*Make cumulative ToF plot 141113*/

//
//       TH1D* mcHist=getTabSeperatedHist(HistsDir+"Hist_MC_AID48_EPG_ToFP.txt",TOFP_DIM2);
//       TH1D* mcHistEP=getTabSeperatedHist(HistsDir+"Hist_MC_AID48_EP_ToFP.txt",TOFP_DIM2);
//
//       TH1D* expHist=getTabSeperatedHist(HistsDir+"Hist_Exp_EID3_EPG_ToFP.txt",TOFP_DIM2);
//       TH1D* expHistEP=getTabSeperatedHist(HistsDir+"Hist_Exp_EID3_EP_ToFP.txt",TOFP_DIM2);
//
//
//
//       scaleHistTo(mcHistEP,1);
//       scaleHistTo(expHistEP,1);
//
//
//
//       TH1D* integralHistMC=makeTH1DFromDim("integralHist","Proton Time of Flight, Integral High to Low;Proton Time of Flight (s);EPG/EP",TOFP_DIM2);
//       TH1D* integralHistEPMC=makeTH1DFromDim("integralHistEP","Proton Time of Flight, Integral High to Low;Proton Time of Flight (s);EPG/EP",TOFP_DIM2);
//       TH1D* integralHistExp=makeTH1DFromDim("integralHistExp","Proton Time of Flight, Integral High to Low;Proton Time of Flight (s);EPG/EP",TOFP_DIM2);
//       TH1D* integralHistEPExp=makeTH1DFromDim("integralHistEPExp","Proton Time of Flight, Integral High to Low;Proton Time of Flight (s);EPG/EP",TOFP_DIM2);
//       integralHistMC->Sumw2();
//       integralHistEPMC->Sumw2();
//       integralHistExp->Sumw2();
//       integralHistEPExp->Sumw2();
//
//        double expIntegral,expError;
//        double mcIntegral, mcError;
//        for (int i = 0;i< TOFP_DIM2.numBins;i++)
//        {
//
//            //High to low
//
//            mcIntegral=mcHist->IntegralAndError(i+1,800,mcError);
//
//            if(mcIntegral != 0)
//            {
//                integralHistMC->SetBinContent(i+1,mcIntegral);
//                integralHistMC->SetBinError(i+1,mcError);
//            }
//
//            mcIntegral=mcHistEP->IntegralAndError(i+1,800,mcError);
//
//			if(mcIntegral != 0)
//			{
//				integralHistEPMC->SetBinContent(i+1,mcIntegral);
//				integralHistEPMC->SetBinError(i+1,mcError);
//			}
//
//			expIntegral=expHist->IntegralAndError(i+1,800,expError);
//
//			if(expIntegral != 0)
//			{
//				integralHistExp->SetBinContent(i+1,expIntegral);
//				integralHistExp->SetBinError(i+1,expError);
//			}
//
//			expIntegral=expHistEP->IntegralAndError(i+1,800,expError);
//
//			if(expIntegral != 0)
//			{
//				integralHistEPExp->SetBinContent(i+1,expIntegral);
//				integralHistEPExp->SetBinError(i+1,expError);
//			}
//
//
//
//
//        }
//
//    TH1* ratioHistMC=makeTH1DFromDim("ratioHistMC","Simulation;Proton Time of Flight (s);EPG/EP",TOFP_DIM2);
//    ratioHistMC->Sumw2();
//    ratioHistMC->Divide(integralHistMC,integralHistEPMC);
//    ratioHistMC->GetXaxis()->SetRangeUser(2e-6,10e-6);
//    convertTH1ToTXT(ratioHistMC,HistsDir+"CumulativeMC_TOFP_HighToLowIntegral_141113.txt");
//
//    TH1* ratioHistExp=makeTH1DFromDim("ratioHistExp","Experiment;Proton Time of Flight (s);EPG/EP",TOFP_DIM2);
//	ratioHistExp->Sumw2();
//	ratioHistExp->Divide(integralHistExp,integralHistEPExp);
//	ratioHistExp->GetXaxis()->SetRangeUser(2e-6,10e-6);
//	ratioHistExp->GetYaxis()->SetRangeUser(0.06e-3,0.1e-3);
//	convertTH1ToTXT(ratioHistExp,HistsDir+"CumulativeExp_TOFP_HighToLowIntegral_141113.txt");
//
//
//
//    TH1* mcHists[1]={ratioHistMC};
//    TH1* expHists[1]={ratioHistExp};
//    plotExpVersusMCToImage(1,expHists,1,mcHists,"Cumulative Proton Time of Flight, Integral High to Low","CumulativeMCExp_TOFP_HighToLowIntegral_141113.png");
	/*Remake RID449*/
	// makeGeantMacroFilesForClusterEventGenerator(0,48, "RID449", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. -1. 0. mm\n/MRK/det/setGammaDetOffset 0. -1. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift y- gamma det+electrostatic mirror +Correct SBD Length Contract
	/* Make RID470-509 for beam and magnetic field shifts*/
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID470", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID471", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID472", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID473", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID474", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID475", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID476", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID477", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID478", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID479", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID480", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID481", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID482", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID483", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID484", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID485", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID486", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID487", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID488", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID489", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID490", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID491", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID492", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID493", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID494", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID495", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID496", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID497", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID498", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID499", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID500", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID501", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID502", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID503", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID504", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID505", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID506", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID507", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID508", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID509", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
	/* Analyze higher stat detector shifts 150204*/
//   	createMCAnalysisAndPlot( "AID52",  "V1.036 Correct SBD Pos 1 mm shift z+ gamma det", 48,  "RID420", "RID421",  48, "RID422",  "RID423", "RID424",  STD_EPCutSetB);
//   	createMCAnalysisAndPlot( "AID53",  "V1.036 Correct SBD Pos 1 mm shift z- gamma det", 48,  "RID425", "RID426",  48, "RID427",  "RID428", "RID429",  STD_EPCutSetB);
//   	createMCAnalysisAndPlot( "AID54",  "V1.036 Correct SBD Pos 1 mm shift x+ gamma det", 48,  "RID430", "RID431",  48, "RID432",  "RID433", "RID434",  STD_EPCutSetB);
//   	createMCAnalysisAndPlot( "AID55",  "V1.036 Correct SBD Pos 1 mm shift x- gamma det", 48,  "RID435", "RID436",  48, "RID437",  "RID438", "RID439",  STD_EPCutSetB);
//   	createMCAnalysisAndPlot( "AID56",  "V1.036 Correct SBD Pos 1 mm shift y+ gamma det", 48,  "RID440", "RID441",  48, "RID442",  "RID443", "RID444",  STD_EPCutSetB);
//   	createMCAnalysisAndPlot( "AID57",  "V1.036 Correct SBD Pos 1 mm shift y- gamma det", 48,  "RID445", "RID446",  48, "RID447",  "RID448", "RID449",  STD_EPCutSetB);
//
//
//   	RDK2AnalysisPlotter* thePlotter;
//
//   	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm z","","AID48 AID52 AID53"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//   	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm z","","AID48 AID52 AID53"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//   	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm x","","AID48 AID54 AID55"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//   	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm x","","AID48 AID54 AID55"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//   	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm y","","AID48 AID56 AID57"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//   	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm y","","AID48 AID56 AID57"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
	/* Double checking False Triple coincidence (false epg) systematic */
	///Calc false EPG
//    RDK2MCAnalysis* theAnal=new RDK2MCAnalysis( "AID48",  "V1.036 Correct SBD Pos", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetB);
//    double error;
//    double answer=theAnal->CalcEPPer3Decay(error);
//    cout << "EP Rate: " << answer << " +/- " << error << endl;
//    answer=theAnal->CalcEPFalseGPerEP(DET_EPG,error);
//    cout << "False EPG: " << answer << " +/- " << error << endl;
//    delete theAnal;
	/* Make remaining events (11-48 million) for RID470-509 for beam and magnetic field shifts*/
	/*  Skipping three body events as statistics are not needed, per Herbert*/
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID472", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID473", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID474", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID477", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID478", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID479", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 0. 0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID482", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID483", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID484", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset 1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID487", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID488", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID489", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/setBeamOffset -1. 0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- neutron beam
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID492", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID493", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID494", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. 1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z+ B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID497", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID498", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID499", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 0.  0. -1. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift z- B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID502", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID503", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID504", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset 1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x+ B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID507", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID508", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
//	makeGeantMacroFilesForClusterEventGenerator(10,38, "RID509", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm\n/MRK/field/setCurrentFieldToModify 0\n/MRK/field/setFieldOffset -1.0  0. 0. mm"); //Length contraction corrected SBD (2 mm from bend) & 1 mm shift x- B Field
	/*From this point on using the alternative model of nonlinearity described in http://iron56.physics.umd.edu:8080/Detector/126 150210*/
////    createMCAnalysisAndPlot( "AID61",  "Standard using newer linearity model", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetB);
////    createMCAnalysisAndPlot( "AID62",  "1 mm shift z+ neutron beam", 10,  "RID470", "RID471",  48, "RID472",  "RID473", "RID474",  STD_EPCutSetB);
////    createMCAnalysisAndPlot( "AID63",  "1 mm shift z- neutron beam", 10,  "RID475", "RID476",  48, "RID477",  "RID478", "RID479",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID64",  "1 mm shift x+ neutron beam", 10,  "RID480", "RID481",  48, "RID482",  "RID483", "RID484",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID65",  "1 mm shift x- neutron beam", 10,  "RID485", "RID486",  48, "RID487",  "RID488", "RID489",  STD_EPCutSetB);
//      createMCAnalysisAndPlot( "AID66",  "1 mm shift z+ B Field", 10,  "RID490", "RID491",  48, "RID492",  "RID493", "RID494",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID67",  "1 mm shift z- B Field", 10,  "RID495", "RID496",  48, "RID497",  "RID498", "RID499",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID68",  "1 mm shift x+ B Field", 10,  "RID500", "RID501",  48, "RID502",  "RID503", "RID504",  STD_EPCutSetB);
//    createMCAnalysisAndPlot( "AID69",  "1 mm shift x- B Field", 10,  "RID505", "RID506",  48, "RID507",  "RID508", "RID509",  STD_EPCutSetB);
//
	//RDK2AnalysisPlotter* thePlotter;
//
////	thePlotter=new RDK2AnalysisPlotter("Neutron Beam +/- 1 mm z","","AID61 AID62 AID63"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Neutron Beam +/- 1 mm x","","AID61 AID64 AID65"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//thePlotter=new RDK2AnalysisPlotter("Magnetic Field +/- 1 mm z","","AID61 AID66 AID67"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Magnetic Field +/- 1 mm x","","AID61 AID68 AID69"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//
////	thePlotter=new RDK2AnalysisPlotter("Neutron Beam +/- 1 mm z","","AID61 AID62 AID63"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Neutron Beam +/- 1 mm x","","AID61 AID64 AID65"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//thePlotter=new RDK2AnalysisPlotter("Magnetic Field +/- 1 mm z","","AID61 AID66 AID67"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Magnetic Field +/- 1 mm x","","AID61 AID68 AID69"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
	/* Fix issue with incorrectly identified RID numbers 150204*/
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID510", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. 1. mm\n/MRK/det/setGammaDetOffset 0. 0. 1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z+ gamma det+electrostatic mirror +Correct SBD Length Contract
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID511", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 0. 0. -1. mm\n/MRK/det/setGammaDetOffset 0. 0. -1. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift z- gamma det+electrostatic mirror +Correct SBD Length Contract
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID512", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset 1. 0. 0. mm\n/MRK/det/setGammaDetOffset 1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x+ gamma det+electrostatic mirror +Correct SBD Length Contract
//	makeGeantMacroFilesForClusterEventGenerator(0,48, "RID513", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/field/setCurrentFieldToModify 1\n/MRK/field/setFieldOffset -1. 0. 0. mm\n/MRK/det/setGammaDetOffset -1. 0. 0. mm\n/MRK/det/setSBDDetectorBendOffset 0 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.31359045121 0 -1.87394264285 mm"); //1 mm shift x- gamma det+electrostatic mirror +Correct SBD Length Contract
	//createMCAnalysisAndPlot( "AID70",  "V1.037 Correct SBD Pos 1 mm shift z+ gamma det", 48,  "RID420", "RID421",  48, "RID422",  "RID423", "RID510",  STD_EPCutSetB);
//	createMCAnalysisAndPlot( "AID71",  "V1.036 Correct SBD Pos 1 mm shift z- gamma det", 48,  "RID425", "RID426",  48, "RID427",  "RID428", "RID511",  STD_EPCutSetB);
//	createMCAnalysisAndPlot( "AID72",  "V1.036 Correct SBD Pos 1 mm shift x+ gamma det", 48,  "RID430", "RID431",  48, "RID432",  "RID433", "RID512",  STD_EPCutSetB);
//	createMCAnalysisAndPlot( "AID73",  "V1.036 Correct SBD Pos 1 mm shift x- gamma det", 48,  "RID435", "RID436",  48, "RID437",  "RID438", "RID513",  STD_EPCutSetB);
//	RDK2AnalysisPlotter* thePlotter;
//
//	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm z","","AID61 AID70 AID71"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm z","","AID61 AID70 AID71"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm x","","AID61 AID72 AID73"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Corrected Gamma Det Position +/- 1 mm x","","AID61 AID72 AID73"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
	/* Plot 10 million SBD x+/- */
//	RDK2AnalysisPlotter* thePlotter;
//
//		thePlotter=new RDK2AnalysisPlotter("SBD +/- 1 mm x","","AID48 AID59 AID60"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//		thePlotter=new RDK2AnalysisPlotter("SBD +/- 1 mm x","","AID48 AID59 AID60"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
	/* Create additional data for SBD Shifts*/
//    makeGeantMacroFilesForClusterEventGenerator(10,38, "RID457", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.6726951503 0 -2.03899024875 mm"); //Length contraction corrected SBD (2 mm from bend)+1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(10,38, "RID458", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset 1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset 0.6726951503 0 -2.03899024875 mm"); //Length contraction corrected SBD (2 mm from bend)+1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(10,38, "RID462", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset -1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -1.2998760527 0 -1.70889503695 mm"); //Length contraction corrected SBD (2 mm from bend)-1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(10,38, "RID463", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT_WBAPD\n/MRK/det/setSBDDetectorBendOffset -1 0 -1.9 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -1.2998760527 0 -1.70889503695 mm"); //Length contraction corrected SBD (2 mm from bend)-1 mm x SBD
//    makeGeantMacroFilesForClusterEventGenerator(10,38, "RID467", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(10,38, "RID468", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT");
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID514", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -4 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift towards bend
//    makeGeantMacroFilesForClusterEventGenerator(0,10, "RID515", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -4 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift towards bend
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID516", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -4 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift towards bend
//    makeGeantMacroFilesForClusterEventGenerator(0,48, "RID517", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT\n/MRK/det/setSBDDetectorBendOffset 0 0 -4 mm\n/MRK/field/setCurrentFieldToModify 2\n/MRK/field/setFieldOffset -0.6601904236 0 -3.945142406 mm"); //After contraction, +2mm SBD shift towards bend
	/*Make EPG/EP plot based on proton/electron energy cut*/
//	const int n=5;
//	double xStart=90;
//	double xIncrement=5;
//	double x[n], y[n], yError[n];
//	TH1D* theHist =new TH1D("theHist","theHist",n,xStart-(0.5*xIncrement),xStart+xIncrement*(n+0.5));
//	for (int i = 0; i < n; ++i)
//	{
//		x[i] = xStart + i*xIncrement;
//		//RDK2CutSet theCut("theCut", "theCut", 100., 800., x[i], 31., 2.e-6, 10.e-6, 10., 800., 0.1, 20.);  //Proton Energy
//		RDK2CutSet theCut("theCut", "theCut", x[i], 800., 10., 31., 2.e-6, 10.e-6, 10., 800., 0.1, 20.);  //Electron Energy
//		RDK2MCAnalysis* theAnal = new RDK2MCAnalysis( "AID61",  "Standard using newer linearity model", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  theCut);
//		y[i] = theAnal->CalcEPGPerEP(DET_EPG, yError[i]);
//		cout << "EPG for lower cut of " << x[i] << "keV: " << y[i] << " +/- " << yError[i] << endl;
//		delete theAnal;
//		theHist->SetBinContent(i+1,y[i]);
//		theHist->SetBinError(i+1,yError[i]);
//	}
//
//	//convertTH1ToTXT(&theHist,HistsDir+"EPG_EP_VS_PROTON_LOWER_CUTOFF.txt");
//	convertTH1ToTXT(theHist,HistsDir+"EPG_EP_VS_ELECTRON_LOWER_CUTOFF.txt");
//	//TH1* theHist=getTabSeperatedHist(HistsDir+"EPG_EP_VS_PROTON_LOWER_CUTOFF.txt",15,2.5,28.5);
//	TH1* theHist=getTabSeperatedHist(HistsDir+"EPG_EP_VS_ELECTRON_LOWER_CUTOFF.txt",5,90-2.5,115+2.5);
//	TCanvas* theCanvas = new TCanvas("theCanvas", "theCanvas", 10, 10, 1200, 1024);
//	theCanvas->SetLeftMargin(.15);
//	theCanvas->SetRightMargin(.05);
//	theCanvas->SetTopMargin(.15);
//	theCanvas->SetBottomMargin(.12);
//	gPad->SetTickx(1);
//	gPad->SetTicky(1);
//	gPad->SetFillColor(kWhite);
//	//theHist->SetTitle("Effect of Proton Lower Energy Cutoff in Simulation;Lower Proton Energy Cutoff (keV);#frac{epg}{ep}");
//	theHist->SetTitle("Effect of Electron Lower Energy Cutoff in Simulation;Lower Electron Energy Cutoff (keV);#frac{epg}{ep}");
//	theHist->GetYaxis()->SetRangeUser(79.e-6,83.e-6);
//	theHist->Draw();
//
////	theCanvas->SaveAs(GraphsDir+"EPG_EP_VS_PROTON_LOWER_CUTOFF.pdf");
//	theCanvas->SaveAs(GraphsDir+"EPG_EP_VS_ELECTRON_LOWER_CUTOFF.pdf");
	/* Make Kevin graph*/
//	RDK2MCAnalysis* theAnal = new RDK2MCAnalysis( "KevinPlot",  "Standard using newer linearity model", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetC);
//	for (int i = 0; i < 12; ++i)
//	{
//		TH1* theHist=theAnal->MakeHist("BGO"+int2str(i+1),FOURBODY,"gChn.BGO"+int2str(i+1)+"EDepTotal",STD_EPCutSetC.CreateMCEPGCut(1,0,0,i+1,DET_EPG),BGOCAL_DIM);
//		convertTH1ToTXT(theHist,HistsDir+"Hist150313_BGO"+int2str(i+1)+"EDepTotal.txt");
//		delete theHist;
//	}
//	double error;
//	double answer=theAnal->CalcEPPer3Decay(error);
//	cout << answer << " +/- " << error << endl;
	/*Make SBD Shift higher stat plots*/
	//Reusing gamma results here as they should remain unaffected
//	createMCAnalysisAndPlot( "AID74",  "2mm z+ shift SBD from bend", 10,  "RID514", "RID515",  48, "RID516",  "RID517", "RID375",  STD_EPCutSetB);
//	createMCAnalysisAndPlot( "AID75",  "1.9mm z- shift SBD from bend", 10,  "RID465", "RID466",  48, "RID467",  "RID468", "RID375",  STD_EPCutSetB);
//	RDK2AnalysisPlotter* thePlotter;
//
//	thePlotter=new RDK2AnalysisPlotter("SBD Position +/- 2 mm z","","AID61 AID74 AID75"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("SBD Position +/- 2 mm z","","AID61 AID74 AID75"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
	/* New proton energy cut*/
	//createMCAnalysisAndPlot( "AID76",  "Standard using new proton energy cut", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetC);
	//createMCAnalysisAndPlot( "AID77",  "Standard using new proton energy cut, 10-25us", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetD);
	/* Check BR for Kevin*/
//	cEvents theEvents;
//	cout << "Fermi BR: " << theEvents.calcBR(10,EKEMAX,0,10000000,true) << endl;
	/* Determine EPG/EP for 2-10 and 10-25*/
//	RDK2MCAnalysis* theAnal1 = new RDK2MCAnalysis( "KevinPlot",  "Standard using newer linearity model", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetC);
//	RDK2MCAnalysis* theAnal2 = new RDK2MCAnalysis( "KevinPlot",  "Standard using newer linearity model", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetD);
//
//	double error;
//	double answer;
//
//	answer=theAnal1->CalcEPGPerEP(DET_EPG,error);
//	cout << "2-10 us: " << answer << " +/- " << error << endl;
//
//	answer=theAnal2->CalcEPGPerEP(DET_EPG,error);
//	cout << "10-25 us: " << answer << " +/- " << error << endl;
	/* Make Kevin graph without incorrect gamma cut*/
//	RDK2MCAnalysis* theAnal = new RDK2MCAnalysis( "KevinPlot",  "Standard using newer linearity model", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetC);
//	for (int i = 0; i < 12; ++i)
//	{
//		TH1* theHist=theAnal->MakeHist("BGO"+int2str(i+1),FOURBODY,"gChn.BGO"+int2str(i+1)+"EDepTotal",STD_EPCutSetC.CreateMCEPCut(1,0),BGOCAL_DIM);
//		convertTH1ToTXT(theHist,HistsDir+"Hist150323_BGO"+int2str(i+1)+"EDepTotal.txt");
//		delete theHist;
//	}
//	TCanvas theCanvas("theCanvas", "theCanvas", 1200, 800);
//	for (int i = 0; i < 12; ++i)
//	{
//		TH1D* theHist=getTabSeperatedHist(HistsDir+"Hist150323_BGO"+int2str(i+1)+"EDepTotal.txt",BGOCAL_DIM);
//		theHist->Rebin(10);
//		theHist->GetYaxis()->SetRangeUser(0,5000);
//		theHist->SetBinContent(1,0);
//		theHist->Draw("L");
//		theCanvas.SaveAs(GraphsDir+"Hist150323_BGO"+int2str(i+1)+"EDepTotal.png");
//
//		delete theHist;
//	}
	/*Checking plotter for individual detectors*/
//	RDK2AnalysisPlotter* thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Beam +/- 1 mm x","","AID61 AID64 AID65"); thePlotter->makeEPGPlot(DET_EPG); delete thePlotter;
//	thePlotter=new RDK2AnalysisPlotter("Beam +/- 1 mm x","","AID61 AID64 AID65"); thePlotter->makeEPGPlot(DET_EPBG); delete thePlotter;
	/*Make cumulative ToF plots*/
//	RDK2MCAnalysis* theAnal = new RDK2MCAnalysis( "NoToF",  "Standard using newer linearity model", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  CutSetNoTOF);
//	TH1* pTHist=theAnal->MakeStandardHist("EP_pTHist",DET_EP,"pChn.SBDTimeFirst",TOFP_DIM3);
//    convertTH1ToTXT(pTHist,HistsDir+"Hist_MC_NoTOF_EP_ToFP.txt");
//    TH1* bottomToTop=makeIntegratedHistogram((TH1D*) pTHist,false);
//    bottomToTop->SetName("bottomToTop");
//    TH1* topToBottom=makeIntegratedHistogram((TH1D*) pTHist,true);
//    convertTH1ToTXT(bottomToTop,HistsDir+"Hist_MC_NoTOF_EP_ToFP_IntegratedBT.txt");
//    convertTH1ToTXT(topToBottom,HistsDir+"Hist_MC_NoTOF_EP_ToFP_IntegratedTB.txt");
//    TH1* stdHist=getTabSeperatedHist(HistsDir+"Hist_MC_NoTOF_EP_ToFP.txt",TOFP_DIM3);
//    TH1* btHist=getTabSeperatedHist(HistsDir+"Hist_MC_NoTOF_EP_ToFP_IntegratedBT.txt",TOFP_DIM3);
//    TH1* tbHist=getTabSeperatedHist(HistsDir+"Hist_MC_NoTOF_EP_ToFP_IntegratedTB.txt",TOFP_DIM3);
//    TH1* integratedHists[2]={tbHist,btHist};
//    TH1* stdHists[1]={stdHist};
//
//
//    btHist->Sumw2();
//    tbHist->Sumw2();
//    btHist->Scale(1./7.466621000000000e+06);
//    tbHist->Scale(1./7.466621000000000e+06);
//
//    stdHist->SetTitle("ToF Histogram;Time of Flight (s);Counts (A.U.)");
//    btHist->SetTitle("Integrated Histogram: X to 30 #mus;Time of Flight (s);Ratio of counts");
//    tbHist->SetTitle("Integrated Histogram: 0 to X;Time of Flight (s);Ratio of counts");
//
//    plotExpVersusMCToImage(0,NULL,1,stdHists,"Standard ToF Hist", "ToFMC.pdf");
	//plotExpVersusMCToImage(0,NULL,2,integratedHists,"Integrated ToF Hists for Standard Settings", "ToFIntegratedHistogramsMC.pdf");
	//plotExpVersusMCToImage(0,NULL,2,integratedHists,"Integrated ToF Hists for Standard Settings", "ToFIntegratedHistogramsMC_Log.pdf",true);
	/* Producing text file of nonlineary functions for Kevin */
//	    const int n=10000;
//	    double x[3][n];
//	    double y[3][n];
//	    EDepModel theModel[3]={ALLEXP_PARAM_MODEL,ALLEXP_PARAM_MODEL_HIGH,ALLEXP_PARAM_MODEL_LOW};
//	    const int colorList[12]= {kBlack, 632, 416, 600, 616, 432,800, 820, 840, 860, 870,880};
//	    TString modelString[3]={"Norm","High","Low"};
//
//	    TGraphErrors* theGraphs[3];
//	    TCanvas* c2 = new TCanvas("canvas","canvas",10,10,1024,768);
//	    c2->SetFillColor(kGray);
//	    c2->cd();
//	    c2->SetLeftMargin(.12);
//	    c2->SetRightMargin(.06);
//	    gPad->SetTickx(1);
//	    gPad->SetTicky(1);
//	    gPad->SetFillColor(kWhite);
//	    TLegend *theLegend;
//	    theLegend = new TLegend(0.8,0.5,0.93,0.9);
//
//	    for (int i = 0;i< 3;i++)
//	    {
//	        for (int j = 0;j< n;j++)
//	        {
//	            x[i][j]=((double) (j+1))/10.;
//	            y[i][j]=evalMattAllExpParamModel2(x[i][j], theModel[i]);
//	        }
//
//
//	    	theGraphs[i]=new TGraphErrors(n,x[i],y[i],NULL,NULL);
//
//	    	theGraphs[i]->SetLineColor(colorList[i]);
//	    	convertTGraphErrorsToTXT(theGraphs[i], HistsDir+"ParamModel"+int2str(i)+".txt");
//	    	if(i==0)
//	    	{
//	    	    theGraphs[i]->Draw("AL");
//	    	    theGraphs[i]->GetYaxis()->SetRangeUser(0.3,1.0);
//	    	    theGraphs[i]->SetTitle("Nonlinearity Model (Alternate);BGO Energy EDep (keV);Energy Relative to 662 keV");
//	    	}
//	    	else
//	    	{
//	    	    theGraphs[i]->Draw("Lsame");
//	    	}
//	    	theLegend->AddEntry(theGraphs[i],modelString[i],"L");
//	    }
//
//	    theLegend->Draw();
//	    theLegend->SetBorderSize(1);
//	    theLegend->SetFillColor(kWhite);
//
//	    TString imagePath=GraphsDir+"ParamModelComparison140526alternate.png";
//	    c2->SaveAs(imagePath);
//	    delete c2;
	/* Make standard histograms for nearly final set */
	//    createMCAnalysisAndPlot( "AID79",  "Standard using new proton energy cut and 2-25 us", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetE);
	/* Compare to Jeff's experimental histogram 150611 */
//	TH1* mcHist[1];
//	TH1* expHist[1];
//	const int numHists=7;
//	string expHistPrefix="Hist_Exp_Jeff10d139d140_";
//	string mcHistPrefix="Hist_MC_AID78_";
//	string histList[numHists]={"EP_EDepE","EP_ToFP","EP_EDepP","EPG_EDepE","EPG_ToFP","EPG_EDepP","EPG_EDepGAvg"};
//	HistDim mcDimList[numHists]={EDEPE_DIM2,TOFP_DIM5,EDEPP_DIM2,EDEPE_DIM2,TOFP_DIM5,EDEPP_DIM2,BGO_DIM2};
//	HistDim expDimList[numHists]={EDEPE_DIM2,TOFP_DIM5,EDEPP_DIM2,EDEPE_DIM2,TOFP_DIM5,EDEPP_DIM2,BGO_DIM2};
//	int mcRebinList[numHists]={10,1,1,20,5,10,50};
//
//	TString titleStringList[numHists]={
//		"Electron Energy Response for EP Coincidence (MC scaled to Exp); Energy Response (keV); Experimental Counts",
//		"Proton Time of Flight for EP Coincidence (MC scaled to Exp); Time of Flight (#mus); Experimental Counts",
//		"Proton Energy Response for EP Coincidence (MC scaled to Exp); Energy Response (keV); Experimental Counts",
//		"Electron Energy Response for EPG Coincidence (MC scaled to Exp); Energy Response (keV); A.U.",
//		"Proton Time of Flight for EPG Coincidence (MC scaled to Exp); Time of Flight (#mus); A.U.",
//		"Proton Energy Response for EPG Coincidence (MC scaled to Exp); Energy Response (keV); A.U.",
//		"BGO Energy Response for EPG Coincidence (MC scaled to Exp); Energy Response (keV); A.U."
//	};
//
//	double maxRange[numHists]={800,30e-6,40,800,30e-6,40,800};
//
//	for (int i = 0; i < numHists; ++i)
//	{
//
//		mcHist[0]=getTabSeperatedHist(HistsDir+mcHistPrefix+histList[i]+".txt",mcDimList[i]);
//		expHist[0]=getTabSeperatedHist(HistsDir+expHistPrefix+histList[i]+".txt",expDimList[i]);
//
//		cout << HistsDir+mcHistPrefix+histList[i]+".txt" << ":  " << mcHist[0]->Integral() << endl;;
//		cout << HistsDir+expHistPrefix+histList[i]+".txt" << ":  " << expHist[0]->Integral() << endl;
//
//		mcHist[0]->Rebin(mcRebinList[i]);
//
//		expHist[0]->SetTitle(titleStringList[i]);//sets axis titles easily
//		TString bareTitle=expHist[0]->GetTitle();
//
//		expHist[0]->GetXaxis()->SetRangeUser(0,maxRange[i]);
//
//		expHist[0]->SetTitle("Exp");
//		mcHist[0]->SetTitle("MC");
//
//
//
//		scaleHistogramsTogether(mcHist[0],expHist[0]);
//
//
//
//		plotExpVersusMCToImage(1,expHist,1,mcHist,bareTitle,"compPlot_"+histList[i]+"150611.pdf");
//
//
//		if(i==6)
//		{
//
//			expHist[0]->Chi2Test(mcHist[0],"WWP");
//
//			plotResidualsToImage(1,expHist,1,mcHist,bareTitle,"compPlotResid_"+histList[i]+"150611.pdf");
//		}
//
//		delete mcHist[0];
//		delete expHist[0];
//	}
	/*Determine Correct SBD Layer*/

//	int deadLayerLength=22;
//	TFile theResultsFile(MRKRunDir+"DeadSiLayerTest"+int2str(deadLayerLength)+"nm.root");
//	TTree* resultsTree=(TTree*) theResultsFile.Get("g");
//	createDerivedResultsFile(0,resultsTree,MRKRunDir+"DeadSiLayerTest"+int2str(deadLayerLength)+"nm_Derv.root");
//	TFile theDervFile(MRKRunDir+"DeadSiLayerTest"+int2str(deadLayerLength)+"nm_Derv.root");
//	TTree* dervTree=(TTree*) theDervFile.Get("rd");
//	gROOT->cd();
//	TH1D* mcHist=makeTH1DFromDim("mcHist","mcHist",EDEPP_DIM2);
//	dervTree->AddFriend(resultsTree,"g");
//	dervTree->Draw("SBDEDepBlur >> mcHist","SBDEDepBlur > 13 && SBDEDepBlur < 31 && g.SBDTimeFirst > 2e-6 && g.SBDTimeFirst < 25e-6");
//	convertTH1ToTXT(mcHist,HistsDir+"DeadSiLayerTest"+int2str(deadLayerLength)+"nm.txt");
//	theResultsFile.Close();
//	theDervFile.Close();
//
//	TH1D* expHist=getTabSeperatedHist(HistsDir+"Hist_Exp_Jeff10d139d140_EP_EDepP.txt",500,0,50);
//
//	cout << "Exp: " << expHist->GetMean() << " +/-" << expHist->GetMeanError() << endl;
//	cout << "MC: " << mcHist->GetMean() << " +/-" << mcHist->GetMeanError() << endl;
//
//	mcHist->Rebin(10);
//	expHist->Rebin(10);
//	mcHist->GetXaxis()->SetRangeUser(10,30);
//	expHist->GetXaxis()->SetRangeUser(10,30);
//	scaleHistogramsTogether(expHist,mcHist);
//
//	TH1* expHists[1]={expHist};
//	TH1* mcHists[1]={mcHist};
//	plotExpVersusMCToImage(1,expHists,1,mcHists,int2str(deadLayerLength)+" nm Dead Layer Length Proton Energy Comparison","DeadSiLayerTest"+int2str(deadLayerLength)+"nm.png");
//	delete mcHist;
//	delete expHist;
	/*Make Nonlinearity function which has greater uncertainty at lower energies*/
//	const int n=10000;
//	double x[3][n];
//	double y[3][n];
//	EDepModel theModel[3]={ALLEXP_PARAM_MODEL,ALLEXP_PARAM_MODEL_HIGH,ALLEXP_PARAM_MODEL_LOW};
//	const int colorList[12]= {kBlack, 632, 416, 600, 616, 432,800, 820, 840, 860, 870,880};
//	TString modelString[3]={"Norm","High","Low"};
//
//	TGraphErrors* theGraphs[3];
//	TCanvas* c2 = new TCanvas("canvas","canvas",10,10,1024,768);
//	c2->SetFillColor(kGray);
//	c2->cd();
//	c2->SetLeftMargin(.12);
//	c2->SetRightMargin(.06);
//	gPad->SetTickx(1);
//	gPad->SetTicky(1);
//	gPad->SetFillColor(kWhite);
//	TLegend *theLegend;
//	theLegend = new TLegend(0.7,0.1,0.9,0.4);
//
//
//
//	for (int i = 0;i< 3;i++)
//	{
//		for (int j = 0;j< n;j++)
//		{
//			x[i][j]=((double) (j+1))/10.;
//			y[i][j]=evalMattAllExpParamModel5(x[i][j], theModel[i]);
//		}
//
//
//		theGraphs[i]=new TGraphErrors(n,x[i],y[i],NULL,NULL);
//
//		theGraphs[i]->SetLineColor(colorList[i]);
//		convertTGraphErrorsToTXT(theGraphs[i], HistsDir+"ParamModel"+int2str(i)+".txt");
//		if(i==0)
//		{
//			theGraphs[i]->Draw("AL");
//			theGraphs[i]->GetYaxis()->SetRangeUser(0.3,1.0);
//			theGraphs[i]->GetXaxis()->SetRangeUser(0,20);
//			theGraphs[i]->SetTitle("Nonproportionality Model 5;BGO Energy EDep (keV);Energy Relative to 662 keV");
//		}
//		else
//		{
//			theGraphs[i]->Draw("Lsame");
//		}
//		theLegend->AddEntry(theGraphs[i],modelString[i],"L");
//	}
//
//	double peakX[1]={3};
//	double peakY[1]={.50};
//	TGraph peakGraph(1,peakX,peakY);
//	peakGraph.SetName("peakGraph");
//	peakGraph.SetMarkerColor(kBlue);
//	peakGraph.SetMarkerStyle(20);
//
//	theLegend->AddEntry(&peakGraph,"Fe Peak","P");
//	peakGraph.Draw("P same");
//
//	double nimX[1]={3};
//	double nimY[1]={.69};
//	double nimYErrLow[1]={.08};
//	double nimYErrHigh[1]={.03};
//	TGraphAsymmErrors nimGraph(1,nimX,nimY,NULL,NULL,nimYErrLow, nimYErrHigh);
//	nimGraph.SetName("nimGraph");
//	nimGraph.SetMarkerColor(kGray);
//	nimGraph.SetLineColor(kGray);
//	nimGraph.SetMarkerStyle(20);
//	theLegend->AddEntry(&nimGraph,"Fe NIM","P");
//	nimGraph.Draw("P same");
//
//	theLegend->Draw();
//	theLegend->SetBorderSize(1);
//	theLegend->SetFillColor(kWhite);
//
//	TString imagePath=GraphsDir+"ParamModel5Comparison.png";
//	c2->SaveAs(imagePath);
//	delete c2;
	/*Test NonPropModels*/
//	for(int i=0;i<1000;i++)
//	{
//		cout << "Energy?:  ";
//		double energy;
//		cin >> energy;
//		cout << energy*evalMattAllExpParamModel(energy, ALLEXP_PARAM_MODEL_LOW) << endl;
//	}
	/* Analysis now with nonprop model*/
////	createMCAnalysisAndPlot( "AID80",  "Standard Updated NonProp", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetE);
////
//	 RDK2MCAnalysis mcAnalysis( "AID81",  "Standard Updated NonProp", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetF);
////	 TH1D histProto("histProto","histProto",20,BGO_VARHIST_LOWEDGES);
////
////	 TH1* normHist=mcAnalysis.MakeRadAvgHist("Normal",DET_EPG,"",0,&histProto);
////	 TH1* highHist=mcAnalysis.MakeRadAvgHist("High",DET_EPG,"",2,&histProto);
////	 TH1* lowHist=mcAnalysis.MakeRadAvgHist("Low",DET_EPG,"",3,&histProto);
////
////	 normHist->Print();
////	 highHist->Print();
////	 lowHist->Print();
////
////	 convertTH1ToTXT(normHist,HistsDir+"ResultsParamModel2Norm.txt");
////	 convertTH1ToTXT(highHist,HistsDir+"ResultsParamModel2High.txt");
////	 convertTH1ToTXT(lowHist,HistsDir+"ResultsParamModel2Low.txt");
////
////	 TH1* theHists[3]={normHist,highHist,lowHist};
////
////	 plotResidualsToImage(3,theHists,0,NULL,"Second model of non-proportionality","ResultsNonPropModel2.pdf");
////
////	 delete normHist;
////	 delete highHist;
////	 delete lowHist;
//
//	 TH1* normHist=mcAnalysis.MakeRadAvgHist("Normal",DET_EPG,"",0,BGO_DIM2);
//	 TH1* highHist=mcAnalysis.MakeRadAvgHist("High",DET_EPG,"",2,BGO_DIM2);
//	 TH1* lowHist=mcAnalysis.MakeRadAvgHist("Low",DET_EPG,"",3,BGO_DIM2);
//
//	 convertTH1ToTXT(normHist,HistsDir+"ResultsParamModel3Norm.txt");
//	 convertTH1ToTXT(highHist,HistsDir+"ResultsParamModel3High.txt");
//	 convertTH1ToTXT(lowHist,HistsDir+"ResultsParamModel3Low.txt");
//
////	 TH1* normHist=getTabSeperatedHist(HistsDir+"ResultsParamModel2NormFixed.txt",BGO_DIM2);
////	 TH1* highHist=getTabSeperatedHist(HistsDir+"ResultsParamModel2HighFixed.txt",BGO_DIM2);
////	 TH1* lowHist=getTabSeperatedHist(HistsDir+"ResultsParamModel2LowFixed.txt",BGO_DIM2);
//
//	 normHist->Rebin(50);
//	 highHist->Rebin(50);
//	 lowHist->Rebin(50);
//
//	 TH1* theHists2[3]={normHist,highHist,lowHist};
//
//	 double integrals[3]={normHist->Integral(),highHist->Integral(),lowHist->Integral()	 };
//
//	 theHists2[0]->SetTitle("Norm: "+d2str(integrals[0]));
//	 theHists2[1]->SetTitle("High: +"+d2str((integrals[1]-integrals[0])*100/integrals[0])+"%");
//	 theHists2[2]->SetTitle("Low: "+d2str((integrals[2]-integrals[0])*100/integrals[0])+"%");
//
//	 plotExpVersusMCToImage(3,theHists2,0,NULL,"Third model of non-proportionality: 20 keV Limit","ResultsNonPropModel320keVLimit.pdf");
	/* Divided Canvas plot of NonpropModels */
//	const int n = 10000;
//	const int numModels=5;
//	const int numExtents=3;
//	double x[numExtents][numModels][n];
//	double y[numExtents][numModels][n];
//	EDepModel theExtents[numExtents] ={ ALLEXP_PARAM_MODEL, ALLEXP_PARAM_MODEL_HIGH, ALLEXP_PARAM_MODEL_LOW };
//	const int colorList[12] ={ kBlack, 632, 416, 600, 616, 432, 800, 820, 840, 860, 870, 880 };
//	TString extentString[numExtents] ={ "Norm", "High", "Low" };
//	TString legendExtraText[numExtents][numModels]={{"","","","",""},{": +1.8%",": +1.1%",": +1.7%",": +1.1%",": +1.1%"},{":  -1.5%",":  -1.1%",":  -2.3%",":  -1.9%",":  -1.6%"}};
//
//	TGraph* cutOffGraphs[numExtents][numModels];
//
//	TGraphErrors* theGraphs[numExtents][numModels];
//	TLegend* theLegend[numModels];
//
//	double peakX[1] ={ 5.9 };
//	double peakY[1] ={ .50 };
//	TGraph peakGraph(1, peakX, peakY);
//	peakGraph.SetName("peakGraph");
//	peakGraph.SetMarkerColor(kBlue);
//	peakGraph.SetMarkerStyle(20);
//
//	double nimX[1] ={ 5.9 };
//	double nimY[1] ={ .69 };
//	double nimYErrLow[1] ={ .08 };
//	double nimYErrHigh[1] ={ .03 };
//	TGraphAsymmErrors nimGraph(1, nimX, nimY, NULL, NULL, nimYErrLow, nimYErrHigh);
//	nimGraph.SetName("nimGraph");
//	nimGraph.SetMarkerColor(kGray);
//	nimGraph.SetLineColor(kGray);
//	nimGraph.SetMarkerStyle(20);
//
//	TCanvas* c2 = new TCanvas("canvas", "canvas", 10, 10, 1000, 1296);
//	c2->Divide(2,3);
//	c2->SetFillColor(kGray);
//
//
//
//
//	for (int k = 0; k < numModels; k++)
//	{
//		c2->cd(k+1);
//		c2->SetLeftMargin(.12);
//		c2->SetRightMargin(.06);
//		gPad->SetTickx(1);
//		gPad->SetTicky(1);
//		gPad->SetFillColor(kWhite);
//
//		theLegend[k] = new TLegend(0.5, 0.1, 0.9, 0.4);
//		for (int i = 0; i < numExtents; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				x[i][k][j] = ((double) (j + 1)) / 10.;
//
//				if(k+1 == 1)
//				{
//					y[i][k][j] = evalMattAllExpParamModel(x[i][k][j], theExtents[i]);
//				}
//				else if(k+1 == 2)
//				{
//					y[i][k][j] = evalMattAllExpParamModel2(x[i][k][j], theExtents[i]);
//				}
//				else if(k+1 == 3)
//				{
//					y[i][k][j] = evalMattAllExpParamModel3(x[i][k][j], theExtents[i]);
//				}
//				else if(k+1 == 4)
//				{
//					y[i][k][j] = evalMattAllExpParamModel4(x[i][k][j], theExtents[i]);
//				}
//				else if(k+1 == 5)
//				{
//					y[i][k][j] = evalMattAllExpParamModel5(x[i][k][j], theExtents[i]);
//				}
//			}
//
//			theGraphs[i][k] = new TGraphErrors(n, x[i][k], y[i][k], NULL, NULL);
//
//			convertTGraphErrorsToTXT(theGraphs[i][k], HistsDir+"ParamModel"+int2str(k+1)+"_"+extentString[i]+".txt");
//
//
//			theGraphs[i][k]->SetLineColor(colorList[i]);
//			theGraphs[i][k]->SetMarkerColor(colorList[i]);
//			theGraphs[i][k]->SetMarkerStyle(20);
//			if(i == 0)
//			{
//				theGraphs[i][k]->Draw("AL");
//				theGraphs[i][k]->GetYaxis()->SetRangeUser(0.3, 1.0);
//				theGraphs[i][k]->GetXaxis()->SetRangeUser(0, 20);
//				theGraphs[i][k]->SetTitle("Nonproportionality Model "+int2str(k+1)+ ";BGO Energy EDep (keV);Energy Relative to 662 keV");
//			}
//			else
//			{
//				theGraphs[i][k]->Draw("Lsame");
//			}
//			theLegend[k]->AddEntry(theGraphs[i][k], extentString[i]+legendExtraText[i][k], "LP");
//			double xCut[3]={14.145,13.691,14.626};
//			if(k != 4 || i==0)
//			{
//				double* xC=new double[1];
//				xC[0]=xCut[i];
//				double* yC=new double[1];
//				yC[0]=10./xCut[i];
//				cutOffGraphs[i][k]=new TGraph(1, xC, yC);
//				cutOffGraphs[i][k]->SetMarkerColor(colorList[i]);
//				cutOffGraphs[i][k]->SetMarkerStyle(20);
//				cutOffGraphs[i][k]->Draw("P same");
//			}
//		}
//		theLegend[k]->AddEntry(&peakGraph, "Fe Peak", "P");
//		peakGraph.Draw("P same");
//
//
//		theLegend[k]->AddEntry(&nimGraph, "Fe NIM", "P");
//		nimGraph.Draw("P same");
//
//		theLegend[k]->Draw();
//		theLegend[k]->SetBorderSize(1);
//		theLegend[k]->SetFillColor(kWhite);
//	}
//
//	TString imagePath = GraphsDir + "ParamModelComparisonAll150709.pdf";
//	c2->SaveAs(imagePath);
//	delete c2;
	/*Playing with sinsodial variations*/
//	vector<double> detectorResults[12];
//	vector<double> totalResults;
//	vector<double> x;
//	vector<double> y;
//	vector<double> yError;
//
//	double sinMean=8.8362e-5;
//	double sinAmp=4.2e-6;
//	double gausSigma=1e-5;
//	for(int j=0;j<12;j++)
//	{
//
//		double phi=TMath::Pi()*j/6;
//		double value=sinAmp*sin(phi)+sinMean;
//		x.push_back(j+1);
//		for(int i=0;i<100000;i++)
//		{
//			double result=gRandom->Gaus(value,gausSigma);
//			totalResults.push_back(result);
//			detectorResults[j].push_back(result);
//		}
//		y.push_back(meanArray(detectorResults[j].size(),&(detectorResults[j][0])));
//		yError.push_back(stDevArray(detectorResults[j].size(),&(detectorResults[j][0])));
//
//	}
//
//	double overallMean=meanArray(totalResults.size(),&totalResults[0]);
//	double overallSTD=stDevArray(totalResults.size(),&totalResults[0]);
//
//	TF1 f1("f1","[0]*sin([1]*x+[3])+[2]",0,12);
//	f1.SetParameters(sinAmp,TMath::Pi()/6.,sinMean,0);
//
//	TCanvas theCanvas("c2","c2");
//	TGraphErrors sinSimData(12,&x[0],&y[0],NULL,&yError[0]);
//	sinSimData.SetTitle("#sigma = "+d2str(gausSigma));
//
//	sinSimData.SetMarkerStyle(11);
//	sinSimData.Fit(&f1,"VMR+");
//
//	double sinFitAmplitude = f1.GetParameter(2);
//	double sinFitAmplitudeError = f1.GetParError(2);
//	sinSimData.Draw("AP");
//	TLatex theLatex;
//	theLatex.SetTextSize(0.035);
//	theLatex.SetNDC(kTRUE);
//	theLatex.DrawLatex(0.55,.8,"Direct: " + d2str(overallMean)+" +/- "+d2str(overallSTD));
//	theLatex.DrawLatex(0.55,.6,"Fit: " + d2str(sinFitAmplitude)+" +/- "+d2str(sinFitAmplitudeError));
//
//
//	theCanvas.SaveAs(GraphsDir+"sinTest3.png");
//
//
//
//	cout << "Overall: " << overallMean << " +/- " << overallSTD << endl;
//	cout << "Fit: " << sinFitAmplitude << " +/- " << sinFitAmplitudeError << endl;
	/*playing with backgrounds*/
//	const int numPoints=350;
//	const int numExperiments=1000;
//	double x[numPoints];
//	double y[numPoints];
//	double averageBackground=13;
//	int numPointsInWindow;
//	for(int i=0;i<numPoints;i++)
//	{
//		x[i]=i+1;
//		numPointsInWindow=0;
//		for(int j=0;j<i+1;j++)
//		{
//			numPointsInWindow+=gRandom->Poisson(13);
//		}
//		y[i]=numPointsInWindow/(averageBackground*(i+1));
//	}
//	TCanvas theCanvas("theCanvas","theCanvas", 10, 10, 1000, 1296);
//	TGraph theGraph(numPoints,x,y);
//	theGraph->SetTitle(";Width of Onpeak Window;Counts in On Peak Window/"
//	theGraph.Draw("AL");
//	theCanvas.SaveAs(GraphsDir+"test.pdf");
//	const double N_D = 850;
//	const double L=1800;
//
//
//
//	vector<double> x;
//
//	vector<double> y1;
//	vector<double> y2;
//	vector<double> y3;
//
//	for (int i = 1; i < 61; i++)
//	{
//		double w_P = i;
//		//const double w_B = L-w_P;
//		const double w_B = 1800;
//		const double N_B = 13 * w_B;
//		const double N_P = N_D+13*w_P;
//		const double lambda_m = N_B / w_B;
//		const double E_lambda = sqrt(N_B) / w_B;
//		x.push_back(w_P);
//		double E_D = sqrt(N_D + w_P * w_P * E_lambda * E_lambda * (0.25 * lambda_m * lambda_m + lambda_m * (1. + 1. / (w_P * E_lambda * E_lambda)) + 1.));
//		y2.push_back(E_D / N_D);
//		E_D = sqrt(N_D);
//		y1.push_back(E_D / N_D);
//	}
//
//	for (int i = 1; i < 61; i++)
//		{
//			double w_P = i;
//			const double w_B = L-w_P;
//			//const double w_B = 1300;
//			const double N_B = 13 * w_B;
//			const double N_P = N_D+13*w_P;
//			const double lambda_m = N_B / w_B;
//			const double E_lambda = sqrt(N_B) / w_B;
//			double E_D = sqrt(N_D + w_P * w_P * E_lambda * E_lambda * (0.25 * lambda_m * lambda_m + lambda_m * (1. + 1. / (w_P * E_lambda * E_lambda)) + 1.));
//			y3.push_back(E_D / N_D);
//		}
//
//	TCanvas theCanvas("theCanvas", "theCanvas", 10, 10, 1000, 1296);
//	theCanvas.SetLeftMargin(.17);
//	theCanvas.SetRightMargin(.1);
//	theCanvas.SetGrid(1,1);
//	gPad->SetTickx(1);
//	gPad->SetTicky(1);
//	TGraph theGraph1(x.size(), &x[0], &y1[0]);
//	TGraph theGraph2(x.size(), &x[0], &y2[0]);
//	TGraph theGraph3(x.size(), &x[0], &y3[0]);
//	theGraph2.SetLineColor(kRed);
//	theGraph3.SetLineColor(kBlue);
//	theGraph1.SetTitle(";Width of Onpeak Window;Relative Uncertainty in Detected Counts");
//	theGraph1.GetYaxis()->SetTitleOffset(2.4);
//	theGraph1.GetYaxis()->SetLimits(0,0.11);
//	theGraph1.GetYaxis()->SetRangeUser(0,0.11);
//	theGraph1.Draw("AL");
//	theGraph2.Draw("L same");
//	theGraph3.Draw("L same");
//	TLegend theLegend(0.55,.7,0.87,0.9);
//	theLegend.AddEntry(&theGraph1,"Only #sqrt{N_{D}}","L");
//	theLegend.AddEntry(&theGraph2,"Eq. 8, Fixed w_{B}","L");
//	theLegend.AddEntry(&theGraph3,"Eq. 8, Fixed L","L");
//	theLegend.Draw("same");
//	theCanvas.SaveAs(GraphsDir + "test.pdf");
	/* 150817 Finish direct detector modification to back layer analysis */
//	HistDim FeDim = { 100, 0, 9.5675 }; //From histogram Tom sent
//	TFile backFile(MRKRunDir + "Results_Fe55_G4RadDecay_150615_Set0.root");
//	TFile noBackFile(MRKRunDir + "Results_Fe55_G4RadDecay_NoBackRegion_150615.root");
//	TTree* backTree = (TTree*) backFile.Get("g");
//	TTree* noBackTree = (TTree*) noBackFile.Get("g");
//
//
//	TH1D* expHists[3];
//	TH1D* mcBackHists[3];
//	TH1D* mcNoBackHists[3];
//
//	double peak, peakError;
//	double tail, tailError;
//
//	double expRatio[3];
//	double expRatioError[3];
//
//	double mcBackRatio[3];
//	double mcBackRatioError[3];
//	double mcNoBackRatio[3];
//	double mcNoBackRatioError[3];
//
//	double det[3] = { 1., 2., 3. };
//	TString mcHistString;
//	for (int i = 0; i < 3; i++)
//	{
//		expHists[i] = getTabSeperatedHist(HistsDir + Form("Hist_Fe_Nov11_Det%i.txt", i + 1), FeDim);
//
//		peak = expHists[i]->IntegralAndError(49, 101, peakError);
//		tail = expHists[i]->IntegralAndError(4, 48, tailError);
//		expRatio[i] = tail / peak;
//		expRatioError[i] = getBasicErrorProgation(tail, peak, tailError, peakError, DIVISION);
//
//		cout << "Det " << i + 1 << endl;
//		cout << "Exp Peak: " << peak << " +/- " << peakError << endl;
//		cout << "Exp Tail: " << tail << " +/- " << tailError << endl;
//		cout << "Exp Ratio " << expRatio[i] << " +/- " << expRatioError[i] << endl;
//
//		mcHistString = Form("mcBackHist%i", i + 1);
//		mcBackHists[i] = makeTH1DFromDim(mcHistString, mcHistString + ";Energy(keV);Counts", FeDim);
//		backTree->Draw(Form("BAPD%iEDepTotal >> ", i + 1) + mcHistString);
//
//		mcHistString = Form("mcNoBackHist%i", i + 1);
//		mcNoBackHists[i] = makeTH1DFromDim(mcHistString, mcHistString + ";Energy(keV);Counts", FeDim);
//		noBackTree->Draw(Form("BAPD%iEDepTotal >> ", i + 1) + mcHistString);
//
//		peak = mcBackHists[i]->IntegralAndError(49, 101, peakError);
//		tail = mcBackHists[i]->IntegralAndError(4, 48, tailError);
//		mcBackRatio[i] = tail / peak;
//		mcBackRatioError[i] = getBasicErrorProgation(tail, peak, tailError, peakError, DIVISION);
//
//		cout << "Det " << i + 1 << endl;
//		cout << "MC Back Peak: " << peak << " +/- " << peakError << endl;
//		cout << "MC Back Tail: " << tail << " +/- " << tailError << endl;
//		cout << "MC Back Ratio " << mcBackRatio[i] << " +/- " << mcBackRatioError[i] << endl;
//
//		peak = mcNoBackHists[i]->IntegralAndError(49, 101, peakError);
//		tail = mcNoBackHists[i]->IntegralAndError(4, 48, tailError);
//		mcNoBackRatio[i] = tail / peak;
//		mcNoBackRatioError[i] = getBasicErrorProgation(tail, peak, tailError, peakError, DIVISION);
//
//		cout << "Det " << i + 1 << endl;
//		cout << "MC No Back Peak: " << peak << " +/- " << peakError << endl;
//		cout << "MC No Back Tail: " << tail << " +/- " << tailError << endl;
//		cout << "MC No Back Ratio " << mcNoBackRatio[i] << " +/- " << mcNoBackRatioError[i] << endl;
//	}
//
//	TGraphErrors expGraph(3, det, expRatio, NULL, expRatioError);
//	TGraphErrors mcBackGraph(3, det, mcBackRatio, NULL, mcBackRatioError);
//	TGraphErrors mcNoBackGraph(3, det, mcNoBackRatio, NULL, mcNoBackRatioError);
//	expGraph.SetTitle("Fe Tail/Peak Ratio;BAPD Detector;Tail counts /Peak counts");
//	expGraph.GetYaxis()->SetTitleOffset(1.9);
//	expGraph.GetYaxis()->SetLimits(0.03,0.16);
//	expGraph.GetYaxis()->SetRangeUser(0.03,0.16);
//	expGraph.SetLineColor(kRed);
//	mcBackGraph.SetLineColor(kBlue);
//	mcNoBackGraph.SetLineColor(kBlack);
//	expGraph.SetMarkerColor(kRed);
//	mcBackGraph.SetMarkerColor(kBlue);
//	mcNoBackGraph.SetMarkerColor(kBlack);
//	expGraph.SetMarkerStyle(8);
//	mcBackGraph.SetMarkerStyle(8);
//	mcNoBackGraph.SetMarkerStyle(8);
//	expGraph.SetMarkerSize(1.5);
//	mcBackGraph.SetMarkerSize(1.5);
//	mcNoBackGraph.SetMarkerSize(1.5);
//
//	TCanvas theCanvas("theCanvas", "theCanvas", 10, 10, 1600, 1400);
//	theCanvas.SetLeftMargin(.17);
//	theCanvas.SetRightMargin(.1);
//	theCanvas.SetGrid(1, 1);
//	gPad->SetTickx(1);
//	gPad->SetTicky(1);
//
//	expGraph.Draw("AP");
//	mcNoBackGraph.Draw("P");
//	mcBackGraph.Draw("P");
//
//
//	TLegend theLegend(0.2,.2,0.5,0.4);
//	theLegend.AddEntry(&expGraph,"Exp","P");
//	theLegend.AddEntry(&mcNoBackGraph,"MC (No Partial at Rear)","P");
//	theLegend.AddEntry(&mcBackGraph,"MC (With Partial at Rear)","P");
//	theLegend.Draw("same");
//	theCanvas.SaveAs(GraphsDir+"FeRatioPlotBoth.pdf");
	/* Make dummy e-g timing plot*/
	///BGO Timing Hist - Make Jeff's latest
//	    TCanvas* theCanvas = getThesisPlotCanvas();
//	    double sf=40e-3;
//
//
//	    TH1D* avgChannelHist=getTabSeperatedHist(HistsDir+"ExpHist_Jeff_EGTime_140710.txt",1800,-650.5,1149.5);
//	    TH1D* avgHist = adjustXScaleHist(avgChannelHist,-650.5*sf,1149.5*sf); //now in micro-s
//
//	    avgHist->SetTitle(";Difference in Arrival Time of Electron and Photon #font[22]{(#mus)};Counts");
//
//	    avgHist->SetLineColor(kRed);
//	    avgHist->SetMarkerColor(kRed);
//	//    avgHist->SetMarkerStyle(4);
//	//    avgHist->Rebin(5);
//
//	    avgHist->GetXaxis()->SetRangeUser(-400*sf,600*sf);
//	    avgHist->GetYaxis()->SetRangeUser(0,2000);
//
//	    avgHist->GetYaxis()->SetTitleOffset(0.9);
//
//	    avgHist->GetXaxis()->SetTitleSize(.05);
//	    avgHist->GetYaxis()->SetTitleSize(.05);
//
//	    avgHist->Draw("P E1");
//
//	    TLine* line[4];
//	    double xLine[4]={-325*sf,-75*sf,75*sf,575*sf};
//	    for (int i =0;i< 4;i++)
//	    {
//	        line[i] = new TLine(xLine[i],0,xLine[i],2000);
//	        line[i]->SetLineStyle(2);
//	        line[i]->Draw("same");
//	    }
//	    avgHist->Draw("P same");
//
//	    TPaveText textBox1(0.21,0.7052632,0.3227425,0.9052632,"brNDC");
//	    textBox1.SetTextSize(.05);
//	    textBox1.AddText("Pre-peak");
//	    textBox1.AddText("Window");
//	    textBox1.AddText("(Background)");
//	    textBox1.SetFillColor(kWhite);
//	    textBox1.SetLineColor(kWhite);
//	    textBox1.SetShadowColor(kWhite);
//	    textBox1.SetLineWidth(0);
//	    textBox1.Draw("same");
//
//
//	    TPaveText textBox2(0.65,0.7052632,0.75,0.9052632,"brNDC");
//	    textBox2.SetTextSize(.05);
//	    textBox2.AddText("Post-peak");
//	    textBox2.AddText("Window");
//	    textBox2.AddText("(Background)");
//	    textBox2.SetFillColor(kWhite);
//	    textBox2.SetLineColor(kWhite);
//	    textBox2.SetShadowColor(kWhite);
//	    textBox2.SetLineWidth(0);
//	    textBox2.Draw("same");
//
//	    TPaveText textBox3(0.39,0.15,0.47,0.3,"brNDC");
//	    textBox3.SetTextSize(.05);
//	    textBox3.AddText("On-peak");
//	    textBox3.AddText("Window");
//	    textBox3.AddText("(Signal)");
//
//	    textBox3.SetFillColor(kWhite);
//	    textBox3.SetLineColor(kWhite);
//	    textBox3.SetShadowColor(kWhite);
//	    textBox3.SetLineWidth(0);
//	    textBox3.Draw("same");
//
//
//
//
//	    theCanvas->SaveAs(GraphsDir+"ExpPlot_EP_EGTime_Dummy.pdf");
	//Make Exp/Sim Plots 140723
	//    RDK2AnalysisPlotter* thePlotter;
	//
	//    thePlotter=new RDK2AnalysisPlotter("Latest and Greatest","EID3","AID48");
	//    TH1* expHist = thePlotter->getExpHist(0,DET_EPG, PLOTVAR_GE);
	//    TH1* mcHist = thePlotter->getMCHist(0,DET_EPG, PLOTVAR_GE);
	//
	//    convertTH1ToTXT(expHist,HistsDir+"GEHistExp.txt");
	//    convertTH1ToTXT(mcHist,HistsDir+"GEHistMC.txt");
//	RDK2AnalysisPlotter* thePlotter;
//	TH1D* expHist;
//	TH1D* mcHist;
//
//	double xMin = 0, xMax = 1;
//
//	thePlotter = new RDK2AnalysisPlotter("Latest and Greatest", "EID3", "AID48");
//
//	CoDet detType;
//	PlotVarType varType;
//
//	for (int k = 0; k < 1; k++)
//	{
//
//		TCanvas* theCanvas;
//		if(k == 0)
//		{
//			theCanvas = getThesisPlotCanvas();
//		}
//		else if(k == 1)
//		{
//			theCanvas = getHalfPresentationPlot2();
//		}
//		else
//		{
//			theCanvas = getHalfPresentationPlot();
//		}
//
//		TH1* tempHist;
//		for (int i = 1; i < 2; i++)
//		{
//			detType = static_cast<CoDet>(i + 1);
//			for (int j = 3; j < 4; j++)
//			{
//				varType = static_cast<PlotVarType>(j);
//
//				if(detType == DET_EP && varType == PLOTVAR_GE) continue;
//
//				if(j < 3)
//				{
//					expHist = thePlotter->getExpHist(0, detType, varType);
//					mcHist = thePlotter->getMCHist(0, detType, varType);
//				}
//				else
//				{
//					expHist = getTabSeperatedHist(HistsDir + "GEHistExp.txt", 780, 10, 780);
//					mcHist = getTabSeperatedHist(HistsDir + "GEHistMC.txt", 780, 10, 780);
//				}
//
//				Color_t theColor;
//				TString yAxisString;
//				TString imageString = "compPlot_";
//				if(detType == DET_EP)
//				{
//					theColor = kBlue;
//					yAxisString = "Counts (AU)";
//					imageString += "EP_";
//				}
//				else if(detType == DET_EPG)
//				{
//					theColor = kRed;
//					yAxisString = "#frac{BGO ep#gamma counts}{1 million ep counts}";
//					mcHist->SetLineStyle(2);
//					mcHist->SetLineWidth(2);
//					mcHist->Scale(1e6);
//					expHist->Scale(1e6);
//					imageString += "EPG_";
//				}
//				else if(detType == DET_EPBG)
//				{
//					theColor = kGreen + 2;
//					yAxisString = "#frac{LAAPD ep#gamma counts}{ep counts}";
//					//                mcHist->SetLineStyle(3);
//					mcHist->SetLineWidth(2);
//					imageString += "EPBG_";
//				}
//
//				TString xAxisString;
//				int rebin;
//				double yMin=-1;
//				if(varType == PLOTVAR_EE)
//				{
//					xAxisString = "Electron Kinetic Energy Detected (keV)";
//					imageString += "EE";
//					rebin = 10;
//
//					if(detType != DET_EP)
//					{
//						rebin = 20;
//					}
//					xMax = 800;
//					xMin = 100;
//
//				}
//				else if(varType == PLOTVAR_PE)
//				{
//					xAxisString = "Proton Kinetic Energy Detected (keV)";
//					imageString += "PE";
//					rebin = 2;
//					xMax = 31;
//					xMin = 7;
//				}
//				else if(varType == PLOTVAR_PT)
//				{
//					xAxisString = "Proton Time-of-Flight (s)";
//					imageString += "PT";
//					xMax = 10e-6;
//					xMin = 2e-6;
//
//					rebin = 9;
//				}
//				else if(varType == PLOTVAR_GE)
//				{
//					imageString += "GE";
//					if(detType == DET_EPG)
//					{
//						xAxisString = "BGO Photon Peak Height (511=511 keV)";
//						xMax = 800;
//						xMin = 10;
//					}
//					else if(detType == DET_EPBG)
//					{
//						xAxisString = "LAAPD Photon Peak Height (5.6=5.6 keV)";
//					}
//					rebin = 10;
//					theCanvas->SetLogx();
//					//yMin=0.0001;
//					//                    tempHist=expHist;  expHist=trimHist(expHist,expHist->GetName()+TString("trim"),10); delete tempHist;
//					//                    tempHist=mcHist;  mcHist=trimHist(mcHist,mcHist->GetName()+TString("trim"),10); delete tempHist;
//
//				}
//				else if(varType == PLOTVAR_GEVAR)
//				{
//					imageString += "GEV";
//					if(detType == DET_EPG)
//					{
//						xAxisString = "BGO Photon Peak Height (511=511 keV)";
//						xMax = 800;
//						xMin = 0;
//					}
//					else if(detType == DET_EPBG)
//					{
//						xAxisString = "LAAPD Photon Peak Height (5.6=5.6 keV)";
//					}
//					//	rebin = 10;
//					//                    tempHist=expHist;  expHist=trimHist(expHist,expHist->GetName()+TString("trim"),10); delete tempHist;
//					//                    tempHist=mcHist;  mcHist=trimHist(mcHist,mcHist->GetName()+TString("trim"),10); delete tempHist;
//
//				}
//
//				expHist->SetMarkerColor(kBlack);
//				expHist->SetLineColor(kBlack);
//				mcHist->SetLineColor(theColor);
//				mcHist->GetSumw2()->Set(0);
//				//expHist->SetMarkerStyle(4);
//				expHist->SetMarkerSize(1);
//				mcHist->GetYaxis()->SetTitleOffset(1.7);
//				//            theCanvas->SetLeftMargin(0.2);
//
//				if(j != 4)
//				{
//					mcHist->Rebin(rebin);
//					expHist->Rebin(rebin);
//				}
//
//				TH1D residualHist;
//				residualHist.Add(expHist,mcHist,1,-1);
//
//
//
//
//				gStyle->SetErrorX(0.0001);
//
//				mcHist->SetTitle(";" + xAxisString + ";" + yAxisString);
//
//				scaleHistogramsTogether(expHist, mcHist);
//				gStyle->SetErrorX(.5);
//				TH1F* hr = theCanvas->DrawFrame(xMin, yMin, xMax, expHist->GetMaximum() * 1.15);
//
//				hr->GetXaxis()->SetTitleSize(0.05);
//				hr->GetYaxis()->SetTitleSize(0.05);
//
//				hr->GetXaxis()->SetLabelSize(0.05);
//				hr->GetYaxis()->SetLabelSize(0.05);
//				hr->GetXaxis()->SetMoreLogLabels(true);
//				hr->GetXaxis()->SetNdivisions(10,10,10,true);
//				hr->GetXaxis()->SetNoExponent(kTRUE);
//
//
//				hr->SetXTitle(xAxisString);
//				hr->SetYTitle(yAxisString);
//
//
//
//				mcHist->Draw("L Same");
//				expHist->Draw("E1 Same");
//
//				TLegend* theLegend = new TLegend(0.60, .65, 0.90, 0.9);
//				theLegend->SetBorderSize(1);
//				theLegend->SetFillColor(kWhite);
//				theLegend->AddEntry(expHist, "Experiment", "ep");
//				theLegend->AddEntry(mcHist, "Theory/Simulation", "L");
//				theLegend->Draw("same");
//
//				imageString += "_Dummy";
//				if(k == 0)
//				{
//					imageString += ".pdf";
//				}
//				else if(k == 1)
//				{
//					imageString += "b.png";
//				}
//				else if(k == 2)
//				{
//					imageString += "c.png";
//				}
//
//				cout << GraphsDir + imageString << endl;
//				theCanvas->SaveAs(GraphsDir + imageString);
//				delete mcHist;
//				delete expHist;
//				delete theLegend;
//				delete hr;
//			}
//		}
//		delete theCanvas;
//	}

	/* Data and residuals Plotting */
//	double dataToResidRatio=3;
//	TString imageString = "GE_DATA_RESIDUAL_DUMMY_LogY.pdf";
//	gStyle->SetErrorX(1);
//	TCanvas theCanvas("largeCanvas", "largeCanvas", 10, 10, 1200, 1200);
//	gPad->SetFillColor(kWhite);
//	TPad pad1("pad1", "Top Pad", 0.0, 1.-1./(dataToResidRatio+1.), 1.0, 1.0, 21);
//	TPad pad2("pad2", "Bottom Pad", 0.0, 0.0, 1.0, 1.-1./(dataToResidRatio+1.), 22);
//
//
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
//	TH1* expHist = getTabSeperatedHist(HistsDir + "GEHistExp.txt", 780, 10, 780);
//	TH1* mcHist = getTabSeperatedHist(HistsDir + "GEHistMC.txt", 780, 10, 780);
//	scaleHistogramsTogether(expHist, mcHist);
//
//	mcHist->SetLineColor(kRed);
//	//mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(10);
//	expHist->Rebin(10);
//	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
//	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 800, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";BGO Photon Peak Height (511=511 keV);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.50, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4, 800, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/*Again now for BAPDS*/

//	double dataToResidRatio = 3;
//	TString imageString = "BGE_DATA_RESIDUAL_DUMMY.pdf";
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
//	const HistDim BAPD_DIMb={38,0.5,10};
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "EXPBen_Hist131009_BAPD_AverageIndiv_EPGb.txt", BAPD_DIMb);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "MC_Hist130904_EDepGBAPD_AverageIndiv_EPG_RID248b.txt", BAPD_DIMb);
//	scaleHistogramsTogether(expHist, mcHist);
//
//	mcHist->SetLineColor(kGreen+3);
//	//mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
////	mcHist->Rebin(10);
////	expHist->Rebin(10);
//
//	pad2.cd();
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, -0.2, 10, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";APD Photon Peak Height (5.9=5.9 keV);#frac{APD ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.50, .7, 0.90, 0.95);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4, 10, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05 * dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05 * dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
//	//	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//
//	theCanvas.SaveAs(GraphsDir + imageString);

	/* Make the final macros! */
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID518", PROTON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //“Final” Run (1)
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID519", ELECTRON, "STD_3B_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //“Final” Run (1)
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID520", PROTON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //“Final” Run (1)
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID521", ELECTRON, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //“Final” Run (1)
//	makeGeantMacroFilesForClusterEventGenerator(0,99, "RID522", GAMMA, "STD_4B_100eV_Endpoint_130429", "/MRK/field/addNamedFields NEW_RDK2_MATT"); //“Final” Run (1)

	/* False EPG Triple Coincidence (Bremmstralung) for BAPD */
//	RDK2MCAnalysis* theAnal = new RDK2MCAnalysis("AID81", "Not quite Final, Standard", 48, "RID406", "RID407", 48, "RID408", "RID409", "RID375", STD_EPCutSetG);
//	double error;
//	double answer = theAnal->CalcEPPer3Decay(error);
//	cout << "EP Rate: " << answer << " +/- " << error << endl;
//	answer = theAnal->CalcEPGPerEP(DET_EPBG, error);
//	cout << "EPG Rate: " << answer << " +/- " << error << endl;
//	answer = theAnal->CalcEPFalseGPerEP(DET_EPBG, error);
//	cout << "False EPG Rate: " << answer << " +/- " << error << endl;
//	delete theAnal;

	/* Final? analysis 150930 */
//	createMCAnalysisAndPlot( "AID82",  "Final?", 48,  "RID518", "RID519",  48, "RID520",  "RID521", "RID522",  STD_EPCutSetG);
//
//	RDK2AnalysisPlotter comPlot("Final?","","AID82");
//	comPlot.makeAll();

	/*Final?? analysis 151015 */
//	createMCAnalysisAndPlot( "AID83",  "Final?? (fixed BAPD cuts)", 48,  "RID518", "RID519",  48, "RID520",  "RID521", "RID522",  STD_EPCutSetH);
//	RDK2AnalysisPlotter comPlot("Final??","","AID83");
//	comPlot.makeAll();

	/*AID 83 for BAPDS*/

//		double dataToResidRatio = 3;
//		TString imageString = "AID83_BAPD.pdf";
//		gStyle->SetErrorX(0);
//		TCanvas theCanvas("largeCanvas", "largeCanvas", 10, 10, 1200, 1200);
//		gPad->SetFillColor(kWhite);
//		TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
//		TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
//		pad1.SetFillColor(kWhite);
//		pad2.SetFillColor(kWhite);
//		pad1.SetLeftMargin(.15);
//		pad2.SetLeftMargin(.15);
//		pad1.SetRightMargin(.05);
//		pad2.SetRightMargin(.05);
//		pad1.SetBottomMargin(0);
//		pad1.SetTopMargin(0);
//		pad2.SetTopMargin(0);
//		pad2.SetBottomMargin(.12);
//		pad1.Draw();
//		pad2.Draw();
//
//		const HistDim BAPD_DIMb={38,0.5,10};
//		TH1D* expHist = getTabSeperatedHist(HistsDir + "bAPD_d236_d237_epg_gEn.txt", BAPD_DIM3);
//		TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPBG_EDepBGAvg.txt", BAPD_DIM3);
//		expHist->Scale(1./2.51376e+07);
//		expHist->Rebin(4);
//		mcHist->Rebin(4);
//		fitMCtoExp(expHist, mcHist,1,expHist->GetNbinsX());
////		scaleHistogramsTogether(expHist, mcHist);
//
//		mcHist->SetLineColor(kGreen+3);
//		//mcHist->SetLineStyle(2);
//		mcHist->SetLineWidth(2);
//
//		expHist->SetMarkerStyle(8);
//		expHist->SetMarkerSize(2);
//		expHist->SetLineWidth(2);
//
//		mcHist->Print();
//		expHist->Print();
//
//		double integral, error;
//		integral=expHist->IntegralAndError(1,expHist->GetNbinsX(),error);
//		cout << integral << " +/- " << error << endl;
//		integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//		cout << integral << " +/- " << error << endl;
//
//		mcHist->Scale(1e6);
//		expHist->Scale(1e6);
//
//
//	//	mcHist->Rebin(10);
//	//	expHist->Rebin(10);
//
//		pad2.cd();
//		pad2.SetGrid();
//
//
//
//		TH1F* hr = pad2.DrawFrame(0, -0.2, 20, expHist->GetMaximum() * 1.15);
//
//		hr->GetXaxis()->SetTitleSize(0.05);
//		hr->GetYaxis()->SetTitleSize(0.05);
//
//		hr->GetXaxis()->SetLabelSize(0.05);
//		hr->GetYaxis()->SetLabelSize(0.05);
//		hr->GetYaxis()->SetTitleOffset(1.2);
//		//hr->GetYaxis()->CenterTitle();
//		hr->SetTitle(";APD Photon Peak Height (5.9=5.9 keV);#frac{APD ep#gamma counts}{1 million ep counts}");
//		pad2.cd();
//		mcHist->Draw("HIST same");
//		expHist->Draw("E1 same");
//
//		TLegend theLegend(0.50, .7, 0.90, 0.95);
//		theLegend.SetBorderSize(1);
//		theLegend.SetFillColor(kWhite);
//		theLegend.AddEntry(expHist, "Experiment", "ep");
//		theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//		theLegend.Draw("same");
//
//		TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//		residualHist->Sumw2();
//		residualHist->Add(mcHist, -1);
//
//		residualHist->SetMarkerColor(kGreen+3);
//		residualHist->SetLineColor(kGreen+3);
//		residualHist->GetXaxis()->SetLabelSize(0.2);
//		residualHist->GetYaxis()->SetLabelSize(0.2);
//		residualHist->GetXaxis()->SetTitleSize(0.2);
//		residualHist->GetYaxis()->SetTitleSize(0.2);
//		//residualHist->GetXaxis()->SetNdivisions(0,0,0);
//
//		pad1.cd();
//		pad1.SetGrid();
//		TH1F* hr1 = pad1.DrawFrame(0, -1.4, 20, 1.4);
//
//
//
//		hr1->GetXaxis()->SetTitleSize(0.05);
//		hr1->GetYaxis()->SetTitleSize(0.05 * dataToResidRatio);
//		hr1->SetTitleOffset(.32, "y");
//
//		hr1->GetXaxis()->SetLabelSize(0.05);
//		hr1->GetYaxis()->SetLabelSize(0.05 * dataToResidRatio);
//		hr1->GetYaxis()->CenterTitle();
//		//	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//
//		hr1->SetTitle("");
//		hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//		hr1->GetYaxis()->CenterTitle();
//
//		TLine lineToDraw(0,0,20,0);
//		lineToDraw.SetLineColor(kBlack);
//		lineToDraw.Draw("same");
//
//		residualHist->Draw("E1 same");
//
//		TLegend theLegend2(0.50, .1, 0.90, 0.3);
//		theLegend2.SetBorderSize(1);
//		theLegend2.SetFillColor(kWhite);
//		theLegend2.AddEntry(residualHist, "Residuals", "ep");
//		theLegend2.Draw("same");
//
//		theCanvas.SaveAs(GraphsDir + imageString);

/* Run old mc with new cuts and dims */
//	createMCAnalysisAndPlot( "AID84",  "Run old mc with new cuts and dims", 48,  "RID406", "RID407",  48, "RID408",  "RID409", "RID375",  STD_EPCutSetH);
//	RDK2AnalysisPlotter thePlotter1("SBD shift x direction","","AID84 AID83"); thePlotter1.makeEPGPlot(DET_EPBG);

/* Calc branching ratio for */
//	 cEvents theEvents;
//	 cout << theEvents.calcBR(0.4,14,0,10000000) << endl;

/*AID 83 for BGOS*/

//	double dataToResidRatio = 3;
//	TString imageString = "AID83_BGO_scaled.pdf";
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
//
////	TH1D* expHist = getTabSeperatedHist(HistsDir + "GEHistExp.txt", BGO_DIM2);
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_gEn.txt", BGO_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_EDepGAvg.txt", BGO_DIM2);
////	TH1D* mcHist = getTabSeperatedHist(HistsDir + "GEHistExp.txt", 780, 10, 780);
//	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(100);
//	expHist->Rebin(100);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 800, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";BGO Photon Peak Height (511=511 keV);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.50, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4, 800, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

/*Electron energy*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EDepE_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_eEn.txt", EDEPE_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_EDepE.txt", EDEPE_DIM2);
//
//	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(100);
//	expHist->Rebin(100);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 800, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Electron Energy Deposited (KeV);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4, 800, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);
/* Three Body 4B*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EP_EDepE_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_ep_eEn.txt", EDEPE_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EP_EDepE.txt", EDEPE_DIM2);
//
////	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1./1e6);
//	expHist->Scale(1./1e6);
//	mcHist->Rebin(100);
//	expHist->Rebin(100);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 800, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Electron Energy Deposited (KeV);Millions of Counts");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4e-2, 800, 1.4e-2);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .7, 0.90, 0.9);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/*Proton energy*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EPG_EDepP_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_pEn.txt", EDEPE_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_EDepP.txt", EDEPE_DIM2);
//
//	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(2);
//	expHist->Rebin(2);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(10, 1e-3, 35, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Proton Energy Deposited (KeV);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(10, -1.4, 35, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/* Three Body*/
//		double dataToResidRatio = 3;
//		TString imageString = "AID83_EP_EDepP_scaled.pdf";
//		gStyle->SetErrorX(0);
//		TCanvas theCanvas("largeCanvas", "largeCanvas", 10, 10, 1200, 1200);
//		gPad->SetFillColor(kWhite);
//		TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
//		TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
//		pad1.SetFillColor(kWhite);
//		pad2.SetFillColor(kWhite);
//		pad1.SetLeftMargin(.15);
//		pad2.SetLeftMargin(.15);
//		pad1.SetRightMargin(.05);
//		pad2.SetRightMargin(.05);
//		pad1.SetBottomMargin(0);
//		pad1.SetTopMargin(0);
//		pad2.SetTopMargin(0);
//		pad2.SetBottomMargin(.12);
//		pad1.Draw();
//		pad2.Draw();
//
//
//
//		TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_ep_pEn.txt", EDEPP_DIM2);
//		TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EP_EDepP.txt", EDEPP_DIM2);
//
//	//	expHist->Scale(1./(2.31112e+08));
//		mcHist->SetLineColor(kRed);
//		mcHist->SetLineStyle(2);
//		mcHist->SetLineWidth(2);
//
//		expHist->SetMarkerStyle(8);
//		expHist->SetMarkerSize(2);
//		expHist->SetLineWidth(2);
//
//		mcHist->Print();
//		double error;
//		double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//		cout << "MC: " << integral << " +/- " << error << endl;
//		expHist->Print();
//
//		mcHist->Scale(1./1e6);
//		expHist->Scale(1./1e6);
//		mcHist->Rebin(1);
//		expHist->Rebin(1);
//		scaleHistogramsTogether(expHist,mcHist);
//	//	makeHistNiceForLog(expHist,0.);
//
//		mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//		pad2.cd();
//	//	pad2.SetLogy(1);
//		pad2.SetGrid();
//
//		TH1F* hr = pad2.DrawFrame(10, 1e-3, 35, expHist->GetMaximum() * 1.15);
//
//		hr->GetXaxis()->SetTitleSize(0.05);
//		hr->GetYaxis()->SetTitleSize(0.05);
//
//		hr->GetXaxis()->SetLabelSize(0.05);
//		hr->GetYaxis()->SetLabelSize(0.05);
//		hr->GetYaxis()->SetTitleOffset(1.2);
//		//hr->GetYaxis()->CenterTitle();
//		hr->SetTitle(";Proton Energy Deposited (KeV);Millions of Counts");
//		pad2.cd();
//		mcHist->Draw("HIST same");
//		expHist->Draw("E1 same");
//
//		TLegend theLegend(0.60, .65, 0.90, 0.9);
//		theLegend.SetBorderSize(1);
//		theLegend.SetFillColor(kWhite);
//		theLegend.AddEntry(expHist, "Experiment", "ep");
//		theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//		theLegend.Draw("same");
//
//		TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//		residualHist->Sumw2();
//		residualHist->Add(mcHist, -1);
//
//		residualHist->SetMarkerColor(kRed);
//		residualHist->SetLineColor(kRed);
//		residualHist->GetXaxis()->SetLabelSize(0.2);
//		residualHist->GetYaxis()->SetLabelSize(0.2);
//		residualHist->GetXaxis()->SetTitleSize(0.2);
//		residualHist->GetYaxis()->SetTitleSize(0.2);
//		//residualHist->GetXaxis()->SetNdivisions(0,0,0);
//
//		pad1.cd();
//		pad1.SetLogy(0);
//		pad1.SetGrid();
//		TH1F* hr1 = pad1.DrawFrame(10, -.06, 35, .06);
//
//		hr1->GetXaxis()->SetTitleSize(0.05);
//		hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//		hr1->SetTitleOffset(.32, "y");
//
//		hr1->GetXaxis()->SetLabelSize(0.05);
//		hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//		hr1->GetYaxis()->CenterTitle();
//	//	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//		hr1->SetTitle("");
//		hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//		hr1->GetYaxis()->CenterTitle();
//		residualHist->Draw("E1 same");
//
//		TLine lineToDraw(0,0,20,0);
//		lineToDraw.SetLineColor(kBlack);
//		lineToDraw.Draw("same");
//
//		TLegend theLegend2(0.50, .1, 0.90, 0.3);
//		theLegend2.SetBorderSize(1);
//		theLegend2.SetFillColor(kWhite);
//		theLegend2.AddEntry(residualHist, "Residuals", "ep");
//		theLegend2.Draw("same");
//		theCanvas.SaveAs(GraphsDir + imageString);

	/*Proton Time of Flight*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EPG_ToFP_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_epTOF.txt", TOFP_DIM5);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_ToFP.txt", TOFP_DIM5);
//
//	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(5);
//	expHist->Rebin(5);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 30e-6, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Proton Time of Flight (#mu s);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4, 30e-6, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/* Three Body*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EP_ToFP_scaled.pdf";
//	gStyle->SetErrorX(0);
//	TCanvas theCanvas("largeCanvas", "largeCanvas", 10, 10, 1200, 1200);
//	gPad->SetFillColor(kWhite);
//	TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
//	TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
//	pad1.SetFillColor(kWhite);
//	pad2.SetFillColor(kWhite);
//	pad1.SetLeftMargin(.15);
//	pad2.SetLeftMargin(.15);
//	pad1.SetRightMargin(.10);
//	pad2.SetRightMargin(.10);
//	pad1.SetBottomMargin(0);
//	pad1.SetTopMargin(0);
//	pad2.SetTopMargin(0);
//	pad2.SetBottomMargin(.12);
//	pad1.Draw();
//	pad2.Draw();
//
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_ep_epTOF.txt", TOFP_DIM5);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EP_ToFP.txt", TOFP_DIM5);
//
////	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(1);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1./1.e6);
//	expHist->Scale(1./1.e6);
//	mcHist->Rebin(1);
//	expHist->Rebin(1);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 30e-6, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Proton Time of Flight (#mu s);Millions of Counts");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4e-2, 30e-6, 1.4e-2);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);


/* With Normalized Residuals */
/*AID 83 for BGOS*/


/*Electron energy*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EDepE_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_eEn.txt", EDEPE_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_EDepE.txt", EDEPE_DIM2);
//
//	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(100);
//	expHist->Rebin(100);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 800, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Electron Energy Deposited (KeV);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4, 800, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);
/* Three Body Electron*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EP_EDepE_scaled_norm.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_ep_eEn.txt", EDEPE_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EP_EDepE.txt", EDEPE_DIM2);
//
////	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1./1e6);
//	expHist->Scale(1./1e6);
//	mcHist->Rebin(100);
//	expHist->Rebin(100);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 800, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Electron Energy Deposited (KeV);Millions of Counts");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4e-1, 800, 1.4e-1);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.04 * dataToResidRatio);
//	hr1->SetTitleOffset(.55, "y");
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
//	TLegend theLegend2(0.2, .77, 0.60, 0.97);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Normalized Residuals", "ep");
//	theLegend2.Draw("same");
//
//	theCanvas.SaveAs(GraphsDir + imageString);

	/*Proton energy*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EPG_EDepP_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_pEn.txt", EDEPE_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_EDepP.txt", EDEPE_DIM2);
//
//	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(2);
//	expHist->Rebin(2);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(10, 1e-3, 35, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Proton Energy Deposited (KeV);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(10, -1.4, 35, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/* Three Body*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EP_EDepP_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_ep_pEn.txt", EDEPP_DIM2);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EP_EDepP.txt", EDEPP_DIM2);
//
////	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1./1e6);
//	expHist->Scale(1./1e6);
//	mcHist->Rebin(1);
//	expHist->Rebin(1);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(10, 1e-3, 35, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Proton Energy Deposited (KeV);Millions of Counts");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(10, -.06, 35, .06);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/*Proton Time of Flight*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EPG_ToFP_scaled.pdf";
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
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_epg_epTOF.txt", TOFP_DIM5);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EPG_ToFP.txt", TOFP_DIM5);
//
//	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(2);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1e6);
//	expHist->Scale(1e6);
//	mcHist->Rebin(5);
//	expHist->Rebin(5);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 30e-6, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Proton Time of Flight (#mu s);#frac{BGO ep#gamma counts}{1 million ep counts}");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
//	theLegend.Draw("same");
//
//	TH1D* residualHist = (TH1D*) expHist->Clone("residualHist");
//	residualHist->Sumw2();
//	residualHist->Add(mcHist, -1);residualHist->Divide(mcHist, 1);
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4, 30e-6, 1.4);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.05*dataToResidRatio);
//	hr1->SetTitleOffset(.32, "y");
//
//	hr1->GetXaxis()->SetLabelSize(0.05);
//	hr1->GetYaxis()->SetLabelSize(0.05*dataToResidRatio);
//	hr1->GetYaxis()->CenterTitle();
////	hr1->SetTitle(";BGO Photon Peak Height (511=511 keV);#splitline{    Exp. -}{Theory/Sim.}");
//	hr1->SetTitle("");
//	hr1->GetYaxis()->SetNdivisions(5, 0, 0);
//	hr1->GetYaxis()->CenterTitle();
//	residualHist->Draw("E1 same");
//
//	TLine lineToDraw(0,0,20,0);
//	lineToDraw.SetLineColor(kBlack);
//	lineToDraw.Draw("same");
//
//	TLegend theLegend2(0.50, .1, 0.90, 0.3);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Residuals", "ep");
//	theLegend2.Draw("same");
//	theCanvas.SaveAs(GraphsDir + imageString);

	/* Three Body*/
//	double dataToResidRatio = 3;
//	TString imageString = "AID83_EP_ToFP_scaled_norm.pdf";
//	gStyle->SetErrorX(0);
//	TCanvas theCanvas("largeCanvas", "largeCanvas", 10, 10, 1200, 1200);
//	gPad->SetFillColor(kWhite);
//	TPad pad1("pad1", "The pad 20% of the height", 0.0, 1. - 1. / (dataToResidRatio + 1.), 1.0, 1.0, 21);
//	TPad pad2("pad2", "The pad 80% of the height", 0.0, 0.0, 1.0, 1. - 1. / (dataToResidRatio + 1.), 22);
//	pad1.SetFillColor(kWhite);
//	pad2.SetFillColor(kWhite);
//	pad1.SetLeftMargin(.15);
//	pad2.SetLeftMargin(.15);
//	pad1.SetRightMargin(.10);
//	pad2.SetRightMargin(.10);
//	pad1.SetBottomMargin(0);
//	pad1.SetTopMargin(0);
//	pad2.SetTopMargin(0);
//	pad2.SetBottomMargin(.12);
//	pad1.Draw();
//	pad2.Draw();
//
//
//
//	TH1D* expHist = getTabSeperatedHist(HistsDir + "BGO_d239_d240_ep_epTOF.txt", TOFP_DIM5);
//	TH1D* mcHist = getTabSeperatedHist(HistsDir + "Hist_MC_AID83_EP_ToFP.txt", TOFP_DIM5);
//
////	expHist->Scale(1./(2.31112e+08));
//	mcHist->SetLineColor(kRed);
//	mcHist->SetLineStyle(2);
//	mcHist->SetLineWidth(2);
//
//	expHist->SetMarkerStyle(8);
//	expHist->SetMarkerSize(1);
//	expHist->SetLineWidth(2);
//
//	mcHist->Print();
//	double error;
//	double integral=mcHist->IntegralAndError(1,mcHist->GetNbinsX(),error);
//	cout << "MC: " << integral << " +/- " << error << endl;
//	expHist->Print();
//
//	mcHist->Scale(1./1.e6);
//	expHist->Scale(1./1.e6);
//	mcHist->Rebin(1);
//	expHist->Rebin(1);
//	scaleHistogramsTogether(expHist,mcHist);
////	makeHistNiceForLog(expHist,0.);
//
//	mcHist->GetYaxis()->SetRangeUser(1e-3, 26);
//
//	pad2.cd();
////	pad2.SetLogy(1);
//	pad2.SetGrid();
//
//	TH1F* hr = pad2.DrawFrame(0, 1e-3, 30e-6, expHist->GetMaximum() * 1.15);
//
//	hr->GetXaxis()->SetTitleSize(0.05);
//	hr->GetYaxis()->SetTitleSize(0.05);
//
//	hr->GetXaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetLabelSize(0.05);
//	hr->GetYaxis()->SetTitleOffset(1.2);
//	//hr->GetYaxis()->CenterTitle();
//	hr->SetTitle(";Proton Time of Flight (#mu s);Millions of Counts");
//	pad2.cd();
//	mcHist->Draw("HIST same");
//	expHist->Draw("E1 same");
//
//	TLegend theLegend(0.60, .65, 0.90, 0.9);
//	theLegend.SetBorderSize(1);
//	theLegend.SetFillColor(kWhite);
//	theLegend.AddEntry(expHist, "Experiment", "ep");
//	theLegend.AddEntry(mcHist, "Theory/Simulation", "L");
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
//	TH1F* hr1 = pad1.DrawFrame(0, -1.4e-1, 30e-6, 1.4e-1);
//
//	hr1->GetXaxis()->SetTitleSize(0.05);
//	hr1->GetYaxis()->SetTitleSize(0.04 * dataToResidRatio);
//	hr1->SetTitleOffset(.55, "y");
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
//	TLegend theLegend2(0.2, .77, 0.60, 0.97);
//	theLegend2.SetBorderSize(1);
//	theLegend2.SetFillColor(kWhite);
//	theLegend2.AddEntry(residualHist, "Normalized Residuals", "ep");
//	theLegend2.Draw("same");
//
//	theCanvas.SaveAs(GraphsDir + imageString);

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

	double bgoCorrection=-.014;
	double apdCorrection=-.003;
	TString bgoExpID="d239_d240";
	TString apdExpID="d236_d237";
	double bgoEPCounts=2.31112e+08;
	double apdEPCounts=2.51376e+07;
	int mcID=87;
	makeOtherParameterPlots(bgoExpID, bgoCorrection, apdExpID,  apdCorrection,  mcID, true);
	return 0;  //end int main()
}
