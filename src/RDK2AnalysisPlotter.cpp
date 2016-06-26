#include <RDK2AnalysisPlotter.h>
#include "MRText.h"
#include "TStyle.h"

RDK2AnalysisPlotter::RDK2AnalysisPlotter()
{
    exp=nullptr;
    mc=nullptr;
    titleString="";
    numExp=0;
    numMC=0;
    deleteAllAnalysis=false;
}

RDK2AnalysisPlotter::RDK2AnalysisPlotter(TString inpTitleString,int inpNumExp, RDK2ExpAnalyzer** inpExp,int inpNumMC, RDK2MCAnalyzer** inpMC)
{
    exp=inpExp;
    mc=inpMC;
    titleString=inpTitleString;
    numExp=inpNumExp;
    numMC=inpNumMC;
    deleteAllAnalysis=false;
}

RDK2AnalysisPlotter::RDK2AnalysisPlotter(TString inpTitleString, TString expIDList,TString mcIDList)
{
    numExp=numItemsInStringList(expIDList);
    numMC=numItemsInStringList(mcIDList);
    mc=nullptr;
    exp=nullptr;
    titleString=inpTitleString;
    deleteAllAnalysis=true;

    if(numExp!=0)
    {
        exp= new RDK2ExpAnalyzer*[numExp];
    }

    if(numMC!=0)
    {
        mc= new RDK2MCAnalyzer*[numMC];
    }


    TString analysisName;
    for (int i = 0;i< numExp;i++)
    {
        TString idString=stringFromList(expIDList,i);
        TString filePath=ANALYSIS_DIR;
        filePath+=idString+".root";
        TFile theFile(filePath,"read");
        RDK2ExpAnalyzer* fromFile=(RDK2ExpAnalyzer*) theFile.Get(idString);
        analysisName=fromFile->GetName();
        gROOT->cd();
        exp[i]=(RDK2ExpAnalyzer*) fromFile->Clone();
        exp[i]->SetName(analysisName);
        theFile.Close();
    }

    for (int i = 0;i< numMC;i++)
    {
        TString idString=stringFromList(mcIDList,i);
        TString filePath=ANALYSIS_DIR;
        filePath+=idString+".root";
        TFile theFile(filePath,"read");
        RDK2MCAnalyzer* fromFile=(RDK2MCAnalyzer*) theFile.Get(idString);
        analysisName=fromFile->GetName();
        gROOT->cd();
        mc[i]=(RDK2MCAnalyzer*) fromFile->Clone();
        mc[i]->SetName(analysisName);
        theFile.Close();
    }



}


RDK2AnalysisPlotter::~RDK2AnalysisPlotter()
{
    clearPlotHists();
    if(deleteAllAnalysis)
    {
        for (int i = 0;i< numExp;i++) delete exp[i];
        for (int i = 0;i< numMC;i++) delete mc[i];
        delete[] exp;
        delete[] mc;
    }
}

void RDK2AnalysisPlotter::clearPlotHists()
{
    for (unsigned int i = 0;i< plotHists.size();i++)
    {
    	delete plotHists[i];
    }
    plotHists.clear();
}

void RDK2AnalysisPlotter::makeAll()
{
    makeEPPlot(DET_EP);
    makeEPPlot(DET_EPG);
    makeEPPlot(DET_EPBG);
    makeEPGPlot(DET_EPG);
    makeEPGPlot(DET_EPBG);
}

void RDK2AnalysisPlotter::makeEPPlot(CoDet detType)
{
    gROOT->cd();
    TCanvas* theCanvas;
    TPad* mainPad;
    TPad* titlePad;
    TPad* pTPad;
    TPad* eEPad;
    TPad* pEPad;
    TPad* pTPadSubs[3];
    TPad* eEPadSubs[3];
    TPad* pEPadSubs[3];


    int numTitleLines=getTitleBoxLines( detType);
    int canvasNumPixelsYPlot=600;
    int canvasNumPixelsYPlots=3*canvasNumPixelsYPlot;
    int canvasNumPixelsYTitle=40*numTitleLines;
    int canvasNumPixelsY=canvasNumPixelsYPlots+canvasNumPixelsYTitle;

    theCanvas = new TCanvas("EPExpMCAnalysisComparisonPlot","EPExpMCAnalysisComparisonPlot",10,10,1200,canvasNumPixelsY);
     mainPad=(TPad*) theCanvas->GetPad(0);


    double ylow,yhigh;

    yhigh=1;
    ylow=1.-canvasNumPixelsYTitle/(double)canvasNumPixelsY;
    titlePad=new TPad("titlePad", "titlePad", 0., ylow, 1., yhigh, -1, 1, 1);

    yhigh=1.-canvasNumPixelsYTitle/(double)canvasNumPixelsY;
    ylow=1.-canvasNumPixelsYTitle/(double)canvasNumPixelsY-canvasNumPixelsYPlot/(double)canvasNumPixelsY;
    pTPad=new TPad("pTPad", "pTPad", 0., ylow, 1., yhigh, -1, 1, 1);

    yhigh=1.-canvasNumPixelsYTitle/(double)canvasNumPixelsY-canvasNumPixelsYPlot/(double)canvasNumPixelsY;
    ylow=1.-canvasNumPixelsYTitle/(double)canvasNumPixelsY- 2*canvasNumPixelsYPlot/(double)canvasNumPixelsY;
    eEPad=new TPad("eEPad", "eEPad", 0., ylow, 1., yhigh, -1, 1, 1);

    yhigh=1.-canvasNumPixelsYTitle/(double)canvasNumPixelsY-2*canvasNumPixelsYPlot/(double)canvasNumPixelsY;
    ylow=0;
    pEPad=new TPad("pEPad", "pEPad", 0., ylow, 1., yhigh, -1, 1, 1);

    titlePad->Draw();
    pTPad->Draw();
    eEPad->Draw();
    pEPad->Draw();



//    titlePad=(TPad*) theCanvas->GetPad(1);
//    pTPad=(TPad*) theCanvas->GetPad(2);
//    eEPad=(TPad*) theCanvas->GetPad(3);
//    pEPad=(TPad*) theCanvas->GetPad(4);

    titlePad->SetMargin(0.05,0.0,0.0,0.0);
    pTPad->SetMargin(0.05,0.0,0.0,0.05);
    pEPad->SetMargin(0.05,0.0,0.0,0.05);
    eEPad->SetMargin(0.05,0.0,0.0,0.05);


    if(numExp+numMC >1)
    {
        pTPad->Divide(2);
        eEPad->Divide(2);
        pEPad->Divide(2);

        pTPadSubs[0]=(TPad*) pTPad->GetPad(1);
        eEPadSubs[0]=(TPad*) eEPad->GetPad(1);
        pEPadSubs[0]=(TPad*) pEPad->GetPad(1);

        TPad* tempPad;
        tempPad=(TPad*) pTPad->GetPad(2);
        tempPad->Divide(1,2);
        pTPadSubs[1]=(TPad*) tempPad->GetPad(1);
        pTPadSubs[2]=(TPad*) tempPad->GetPad(2);

        tempPad=(TPad*) pEPad->GetPad(2);
        tempPad->Divide(1,2);
        pEPadSubs[1]=(TPad*) tempPad->GetPad(1);
        pEPadSubs[2]=(TPad*) tempPad->GetPad(2);

        tempPad=(TPad*) eEPad->GetPad(2);
        tempPad->Divide(1,2);
        eEPadSubs[1]=(TPad*) tempPad->GetPad(1);
        eEPadSubs[2]=(TPad*) tempPad->GetPad(2);

        TPad* allBasePads[10]={titlePad,pTPadSubs[0],pTPadSubs[1],pTPadSubs[2],eEPadSubs[0],eEPadSubs[1],eEPadSubs[2],pEPadSubs[0],pEPadSubs[1],pEPadSubs[2]};
        for (int i = 1;i< 10;i++)
        {
            allBasePads[i]->SetGrid(1,1);
            allBasePads[i]->SetTickx(1);
            allBasePads[i]->SetTicky(1);
        }
        pTPadSubs[0]->SetMargin(0.13,0.1,.1,0.1);
        pEPadSubs[0]->SetMargin(0.13,0.1,.1,0.1);
        eEPadSubs[0]->SetMargin(0.13,0.1,.1,0.1);

        pTPadSubs[1]->SetMargin(0.1,0.1,.1,0.1);
        pEPadSubs[1]->SetMargin(0.1,0.1,.1,0.1);
        eEPadSubs[1]->SetMargin(0.1,0.1,.1,0.1);

        pTPadSubs[2]->SetMargin(0.1,0.1,.1,0.1);
        pEPadSubs[2]->SetMargin(0.1,0.1,.1,0.1);
        eEPadSubs[2]->SetMargin(0.1,0.1,.1,0.1);
        //Plot pT

        pTPadSubs[0]->cd(); drawPlot(detType, PLOTVAR_PT,PLOT_COMP);
        pTPadSubs[1]->cd(); drawPlot(detType, PLOTVAR_PT,PLOT_RESID);
        pTPadSubs[2]->cd(); drawPlot(detType, PLOTVAR_PT,PLOT_NORMRESID);


        //Plot eE
        eEPadSubs[0]->cd(); drawPlot(detType, PLOTVAR_EE,PLOT_COMP);
        eEPadSubs[1]->cd(); drawPlot(detType, PLOTVAR_EE,PLOT_RESID);
        eEPadSubs[2]->cd(); drawPlot(detType, PLOTVAR_EE,PLOT_NORMRESID);

        //Plot pE
        pEPadSubs[0]->cd(); drawPlot(detType, PLOTVAR_PE,PLOT_COMP);
        pEPadSubs[1]->cd(); drawPlot(detType, PLOTVAR_PE,PLOT_RESID);
        pEPadSubs[2]->cd(); drawPlot(detType, PLOTVAR_PE,PLOT_NORMRESID);
    }
    else
    {
        TPad* allBasePads[4]={titlePad,pTPad,eEPad,pEPad};
        for (int i = 1;i< 4;i++)
        {
            allBasePads[i]->SetGrid(1,1);
            allBasePads[i]->SetTickx(1);
            allBasePads[i]->SetTicky(1);
        }
        pTPad->SetMargin(0.13,0.1,.1,0.1);
        pEPad->SetMargin(0.13,0.1,.1,0.1);
        eEPad->SetMargin(0.13,0.1,.1,0.1);
        pTPad->cd(); drawPlot(detType, PLOTVAR_PT,PLOT_COMP);
        eEPad->cd(); drawPlot(detType, PLOTVAR_EE,PLOT_COMP);
        pEPad->cd(); drawPlot(detType, PLOTVAR_PE,PLOT_COMP);

    }
    mainPad->SetFillColor(kGray);

    ///Make Title box
    titlePad->cd();
    titlePad->SetFillColor(kGray);
    TPaveText* titleBox = makeTitleBox( detType);
    titleBox->Draw();





    TString coTypeString;
    if(detType==DET_EP)
    {
        coTypeString="EP";
    }
    else if(detType==DET_EPG)
    {
        coTypeString="EPG";
    }
    else if(detType==DET_EPBG)
    {
        coTypeString="EPBG";
    }
     TString imagePath=GRAPHS_DIR;
    imagePath+="layouts/PlotLayout";
    if(numMC>0) imagePath+="_MC";
    for (int i = 0;i< numMC;i++)
    {
        imagePath+=TString("_")+mc[i]->GetName();
    }
    if(numExp>0) imagePath+="_Exp";
    for (int i = 0;i< numExp;i++)
    {
        imagePath+=TString("_")+exp[i]->GetName();
    }

    imagePath+="_"+coTypeString+"_EPPlots.pdf";
    theCanvas->SaveAs(imagePath);

    delete theCanvas;
    clearPlotHists();

}

void RDK2AnalysisPlotter::makeEPGPlot(CoDet detType)
{
     gROOT->cd();
    TCanvas* theCanvas;
    TPad* mainPad;
    TPad* titlePad;
    TPad* gEPad;
    TPad* detPad;
    TPad* gEPadSubs[3];


    int numTitleLines=getTitleBoxLines( detType);
    int canvasNumPixelsYPlot=600;
    int canvasNumPixelsYPlots=2*canvasNumPixelsYPlot;
    int canvasNumPixelsYTitle=40*numTitleLines;
    int canvasNumPixelsY=canvasNumPixelsYPlots+canvasNumPixelsYTitle;

    theCanvas = new TCanvas("EPGExpMCAnalysisComparisonPlot","EPGExpMCAnalysisComparisonPlot",10,10,1200,canvasNumPixelsY);
     mainPad=(TPad*) theCanvas->GetPad(0);


    double ylow,yhigh;

    yhigh=1;
    ylow=1.-canvasNumPixelsYTitle/(double)canvasNumPixelsY;
    titlePad=new TPad("titlePad", "titlePad", 0., ylow, 1., yhigh, -1, 1, 1);

    yhigh=ylow;
    ylow=ylow-canvasNumPixelsYPlot/(double)canvasNumPixelsY;
    gEPad=new TPad("gEPad", "gEPad", 0., ylow, 1., yhigh, -1, 1, 1);

    yhigh=ylow;
    ylow=0;
    detPad=new TPad("detPad", "detPad", 0., ylow, 1., yhigh, -1, 1, 1);

    titlePad->Draw();
    gEPad->Draw();
    detPad->Draw();



//    titlePad=(TPad*) theCanvas->GetPad(1);
//    gEPad=(TPad*) theCanvas->GetPad(2);
//    detPad=(TPad*) theCanvas->GetPad(3);
//    egTPad=(TPad*) theCanvas->GetPad(4);

    titlePad->SetMargin(0.05,0.0,0.0,0.0);
    gEPad->SetMargin(0.05,0.0,0.0,0.05);

    if(numExp>1 && numMC>1)
    {
        gEPad->Divide(2);

        gEPadSubs[0]=(TPad*) gEPad->GetPad(1);

        TPad* tempPad;
        tempPad=(TPad*) gEPad->GetPad(2);
        tempPad->Divide(1,2);
        gEPadSubs[1]=(TPad*) tempPad->GetPad(1);
        gEPadSubs[2]=(TPad*) tempPad->GetPad(2);


        TPad* allBasePads[5]={titlePad,gEPadSubs[0],gEPadSubs[1],gEPadSubs[2],detPad};
         for (int i = 1;i< 5;i++)
        {
            allBasePads[i]->SetGrid(1,1);
            allBasePads[i]->SetTickx(1);
            allBasePads[i]->SetTicky(1);
        }

        gEPadSubs[2]->SetLogx();

        gEPadSubs[0]->SetMargin(0.15,0.1,.1,0.1);
        gEPadSubs[1]->SetMargin(0.1,0.1,.1,0.1);
        gEPadSubs[2]->SetMargin(0.1,0.1,.1,0.1);

        //Plot gE

        gEPadSubs[0]->cd(); drawPlot(detType, PLOTVAR_GE,PLOT_COMP);
        gEPadSubs[1]->cd(); drawPlot(detType, PLOTVAR_GE,PLOT_RESID);
        gEPadSubs[2]->cd(); drawPlot(detType, PLOTVAR_GEVAR,PLOT_NORMRESID);

    }
    else
    {

        TPad* allBasePads[3]={titlePad,gEPad,detPad};
        for (int i = 1;i< 3;i++)
        {
            allBasePads[i]->SetGrid(1,1);
            allBasePads[i]->SetTickx(1);
            allBasePads[i]->SetTicky(1);
        }

        gEPad->SetMargin(0.13,0.1,.1,0.1);


        //Plot gE

        gEPad->cd(); drawPlot(detType, PLOTVAR_GE,PLOT_COMP);
    }
    mainPad->SetFillColor(kGray);

    ///Make Title box
    titlePad->cd();
    titlePad->SetFillColor(kGray);
    TPaveText* titleBox = makeTitleBox( detType);
    titleBox->Draw();

    detPad->SetMargin(0.1,0.05,.1,0.1);

    //Plot det
     detPad->cd(); drawPlot(detType, PLOTVAR_GE,PLOT_DETS);

    TString coTypeString;
    if(detType==DET_EP)
    {
        coTypeString="EP";
    }
    else if(detType==DET_EPG)
    {
        coTypeString="EPG";
    }
    else if(detType==DET_EPBG)
    {
        coTypeString="EPBG";
    }
    TString imagePath=GRAPHS_DIR;
    imagePath+="layouts/PlotLayout";
    if(numMC) imagePath+="_MC";
    for (int i = 0;i< numMC;i++)
    {
        imagePath+=TString("_")+mc[i]->GetName();
    }
    if(numExp) imagePath+="_Exp";
    for (int i = 0;i< numExp;i++)
    {
        imagePath+=TString("_")+exp[i]->GetName();
    }

    imagePath+="_"+coTypeString+"_EPGPlots.pdf";
    theCanvas->SaveAs(imagePath);

    delete theCanvas;
    clearPlotHists();

}

int RDK2AnalysisPlotter::getTitleBoxLines(CoDet detType)
{
     int numLinesPerExpAnalysis;
     int numLinesPerMCAnalysis;
     if(detType!=DET_EP)
     {
        numLinesPerExpAnalysis=4;
        numLinesPerMCAnalysis=3;
     }
     else
     {
        numLinesPerExpAnalysis=2;
        numLinesPerMCAnalysis=2;
     }
    return 2+numExp*numLinesPerExpAnalysis+numMC*numLinesPerMCAnalysis;
}

TPaveText* RDK2AnalysisPlotter::makeTitleBox(CoDet detType)
{
    TString line;
    TPaveText* titleBox = new TPaveText(.00,.0,1,1,"NDC");
    titleBox->SetFillColor(kWhite);
    titleBox->SetLineColor(kWhite);
    titleBox->SetTextAlign(11);
    titleBox->SetTextFont(102);

    //Make outline for title
//    TBox* theBox;
//    int numLines=getTitleBoxLines( detType);
//    int numLinesPerExpAnalysis;
//    int numLinesPerMCAnalysis;
//    int currentLine=0;
//    double y0, y1;
//    if(detType!=DET_EP)
//    {
//        numLinesPerExpAnalysis=4;
//        numLinesPerMCAnalysis=3;
//    }
//    else
//    {
//        numLinesPerExpAnalysis=2;
//        numLinesPerMCAnalysis=2;
//    }

//    y0=1.-1/(double)numLines;
//    y1=1;
//    theBox=titleBox->AddBox(0,y0,1.,y1);
//    theBox->SetFillStyle(0);

//    titleBox->AddText("");

    TString coTypeString;
    if(detType==DET_EP)
    {
        coTypeString="EP";
    }
    else if(detType==DET_EPG)
    {
        coTypeString="EPG (BGO)";
    }
    else if(detType==DET_EPBG)
    {
        coTypeString="EPBG (BAPD)";
    }
    double answer,error;
    titleBox->AddText("Title: "+titleString+"  Coincidence: "+coTypeString+"  Date made: "+getCurrentDateString());
    for (int i = 0;i< numExp;i++)
    {
//        y0=y0-numLinesPerExpAnalysis/(double) numLines;
//        y1=y0;
//        theBox=titleBox->AddBox(0,y0,1,y1);
//        theBox->SetFillStyle(0);

        int theColor=colorList[(i+4)%GRAPHICS_NUM_TYPES];
        line="#color["+int2str(theColor)+"]{";
        line+=TString("Exp: ") + exp[i]->GetName()+" \"" + exp[i]->GetTitle() + "\" " + exp[i]->GetPIDString()+" " + exp[i]->GetGroupString();
        line+=TString("  EP Counts: ") + int2str(exp[i]->GetEPCounts());
        line+="}";
    	titleBox->AddText(line);

    	if(detType==DET_EPG)
        {
            line="#color["+int2str(theColor)+"]{";
            answer=exp[i]->GetEPGCounts(error);
            line+=TString("     EPG Counts: ") + int2str((int) (answer+0.5))+" #pm " + int2str((int) (error+0.5));
            answer=exp[i]->GetEPGPerEP(error);
            line+=TString("  #frac{EPG}{EP}: ") + d2str(answer,2,2)+" #pm " + d2str(error,1,2) +" (#pm"+d2str(error*100./answer,1,1)+"%)";
            line +="}";
            titleBox->AddText(line);
        }
        else if(detType==DET_EPBG)
        {
            line="#color["+int2str(theColor)+"]{";
            answer=exp[i]->GetEPBGCounts(error);
            line+=TString("     EPBG Counts: ") + int2str((int) (answer+0.5))+" #pm " + int2str((int) (error+0.5));
            answer=exp[i]->GetEPBGPerEP(error);
            line+=TString("  #frac{EPBG}{EP}: ") + d2str(answer,2,2)+" #pm " + d2str(error,1,2)+" (#pm"+d2str(error*100./answer,1,1)+"%)";
            line +="}";
            titleBox->AddText(line);
        }
        //Third Line
        if(detType != DET_EP)
        {
            line="#color["+int2str(theColor)+"]{";
            line+="     Prepeak BG Window: ["+d2str(exp[i]->GetWindowVal(detType,0))+","+d2str(exp[i]->GetWindowVal(detType,1))+"]";
            line+=" Peak Window: ["+d2str(exp[i]->GetWindowVal(detType,2))+","+d2str(exp[i]->GetWindowVal(detType,3))+"] ";
            line+=" Postpeak BG Window: ["+d2str(exp[i]->GetWindowVal(detType,4))+","+d2str(exp[i]->GetWindowVal(detType,5))+"] ";
            line +="}";
            titleBox->AddText(line);
        }

        //Fourth Line
        line="#color["+int2str(theColor)+"]{";
        line+="     "+exp[i]->GetCutSet().GetDisplayString( detType);
        line +="}";
        titleBox->AddText(line);


    }

    //Plot MC
    for (int i = 0;i< numMC;i++)
    {
//        y0=y0-numLinesPerMCAnalysis/(double) numLines;
//        y1=y0;
//        theBox=titleBox->AddBox(0,y0,1,y1);
//        theBox->SetFillStyle(0);
//        currentLine+=numLinesPerMCAnalysis;

        int theColor=colorList[(i)%GRAPHICS_NUM_TYPES];
        //First MC Line
        line="#color["+int2str(theColor)+"]{";
        line+=TString("MC: ") + mc[i]->GetName()+" \"" + mc[i]->GetTitle() + "\"";
        line+=" P_{3}"+mc[i]->GetP3ID()(3,3)+" E_{3}"+mc[i]->GetE3ID()(3,3);
        if(detType==DET_EPG || detType==DET_EPBG)
        {
            line+=" P_{4}"+mc[i]->GetP4ID()(3,3)+" E_{4}"+mc[i]->GetE4ID()(3,3)+" G_{4}"+mc[i]->GetG4ID()(3,3);
        }
        answer=mc[i]->GetEPPer3Decay(error);
    	line+=TString("  EP per 3 Decay: ") + d2str(answer,5,1)+" #pm " + d2str(error,5,1);
        line +="}";
    	titleBox->AddText(line);

    	//Second MC Line

    	if(detType==DET_EPG)
        {
            line="#color["+int2str(theColor)+"]{";
            answer=mc[i]->GetEPGPer4Decay(error);
            line+=TString("    EPG per 4 Decay: ") + d2str(answer,6,1)+" #pm " + d2str(error,6,1);
            answer=mc[i]->GetEPGPerEP(error);
            line+=TString("  #frac{EPG}{EP}: ") + d2str(answer,3,2)+" #pm " + d2str(error,1,2)+" (#pm"+d2str(error*100./answer,1,1)+"%)";
            if(numExp>0)
            {
                double expAnswer,expError;
                double mcAnswer, mcError;
                expAnswer=exp[0]->GetEPGPerEP(expError); mcAnswer=mc[i]->GetEPGPerEP(mcError);
                answer=(expAnswer/mcAnswer)-1;
                error=getBasicErrorProgation(expAnswer,mcAnswer,expError,mcError,DIVISION);
                line+=TString("  #frac{Exp - MC}{MC}: ") +d2str(answer*100,1,1)+"% #pm " + d2str(error*100,1,1)+"%";
            }
            line +="}";
            titleBox->AddText(line);
        }
        else if(detType==DET_EPBG)
        {
            line="#color["+int2str(theColor)+"]{";
            answer=mc[i]->GetEPBGPer4Decay(error);
            line+=TString("    EPBG per 4 Decay: ") + d2str(answer,6,1)+" #pm " + d2str(error,6,1);
            answer=mc[i]->GetEPBGPerEP(error);
            line+=TString("  #frac{EPBG}{EP}: ") + d2str(answer,3,2)+" #pm " + d2str(error,1,2)+" (#pm"+d2str(error*100./answer,1,1)+"%)";
            if(numExp>0)
            {
                double expAnswer,expError;
                double mcAnswer, mcError;
                expAnswer=exp[0]->GetEPGPerEP(expError); mcAnswer=mc[i]->GetEPGPerEP(mcError);
                answer=(expAnswer/mcAnswer)-1;
                error=getBasicErrorProgation(expAnswer,mcAnswer,expError,mcError,DIVISION);
                line+=TString("  #frac{Exp - MC}{MC}: ") +d2str(answer*100,1,1)+"% #pm " + d2str(error*100,1,1)+"%";
            }
            line +="}";
            titleBox->AddText(line);
        }
        //Third Line
        line="#color["+int2str(theColor)+"]{";
        line+="    "+mc[i]->GetCutSet().GetDisplayString( detType);
        line +="}";
        titleBox->AddText(line);

    }



    return titleBox;
}

void RDK2AnalysisPlotter::drawHist(int inpNumExp, TH1** expHists,int inpNumMC, TH1** mcHists)
{
    gStyle->SetTitleX(0.3);
    gROOT->cd();
    int totalNumHists=inpNumExp+inpNumMC;
    TH1* totalHists[GRAPHICS_NUM_TYPES];
    TString drawString;
    for (int i = 0;i< totalNumHists;i++)
    {

        if(i<inpNumExp)
        {
            totalHists[i]=expHists[i];
            drawString="P E1 X0";
        }
        else
        {
            totalHists[i]=mcHists[i-inpNumExp];
            drawString="HIST";

        }




        if(i!=0)
        {
            drawString+=" SAME";
        }


        totalHists[i]->Draw(drawString);

    }
}


void RDK2AnalysisPlotter::drawPlot(CoDet detType, PlotVarType plotVarType,PlotType plotType)
{
    gROOT->cd();



    TString histName;

    if(plotType ==PLOT_COMP)
    {
        drawPlotComp( detType,  plotVarType);
    }
    else if(plotType == PLOT_RESID && (numExp+numMC>1 ))
    {
        drawPlotResid( detType,  plotVarType, false);
    }
    else if(plotType == PLOT_NORMRESID && (numExp+numMC>1 ))
    {
        drawPlotResid( detType,  plotVarType, true);
    }
    else if(plotType == PLOT_DETS)
    {
        drawPlotDets(detType,plotVarType);
    }

}



void RDK2AnalysisPlotter::drawPlotComp(CoDet detType, PlotVarType plotVarType)
{
    gROOT->cd();


    TString histName;

    TH1** expHists=new TH1*[numExp];
    TH1** mcHists=new TH1*[numMC];
    int start=plotHists.size();
    //We will use these hists straight off
    for (int i = 0;i< numExp;i++)
    {
        histName="hist"+int2str(plotHists.size());
        expHists[i]=(TH1*) (getExpHist( i, detType, plotVarType)->Clone(histName));
        plotHists.push_back(expHists[i]);
        setDisplayQuantitiesForHist(i,1,expHists[i],detType,PLOT_COMP,plotVarType);
    }
    for (int i = 0;i< numMC;i++)
    {
        histName="hist"+int2str(plotHists.size());
        mcHists[i]=(TH1*) (getMCHist( i, detType, plotVarType)->Clone(histName));
        plotHists.push_back(mcHists[i]);
        setDisplayQuantitiesForHist(-i-1,1,mcHists[i],detType,PLOT_COMP,plotVarType);
        if(detType==DET_EP && numExp>0)//Scale MC hists to first Exp Integral
        {
            scaleHistogramsTogether(expHists[0],mcHists[i],1,expHists[0]->GetNbinsX());
        }

    }

    double max=0;
    for (int i = start;i< start+numExp+numMC;i++)
    {
        if(plotHists[i]->GetMaximum() > max) max=plotHists[i]->GetMaximum()+plotHists[i]->GetBinError(plotHists[i]->GetMaximumBin());
    }

    for (int i = start;i< start+numExp+numMC;i++)
    {
    	plotHists[i]->GetYaxis()->SetRangeUser(0,1.01*max);
    }


    drawHist(numExp, expHists, numMC, mcHists);
    delete[] expHists;
    delete[] mcHists;

}

void RDK2AnalysisPlotter::drawPlotResid(CoDet detType, PlotVarType plotVarType,bool normalize)
{
    gROOT->cd();
    PlotType plotType;
    if(normalize)
    {
        plotType=PLOT_NORMRESID;
    }
    else
    {
        plotType=PLOT_RESID;
    }


    TString histName;

    TH1** expHists=new TH1*[numExp];
    TH1** mcHists=new TH1*[numMC];
    int start=plotHists.size();
    //We will use these hists straight off
    for (int i = 0;i< numExp;i++)
    {
        histName="hist"+int2str(plotHists.size());
        expHists[i]=(TH1*) (getExpHist( i, detType, plotVarType)->Clone(histName));
        plotHists.push_back(expHists[i]);
    }
    for (int i = 0;i< numMC;i++)
    {
        histName="hist"+int2str(plotHists.size());
        mcHists[i]=(TH1*) (getMCHist( i, detType, plotVarType)->Clone(histName));
        plotHists.push_back(mcHists[i]);

        if(detType==DET_EP && numExp>0)//Scale MC hists to first Exp Integral
        {
            scaleHistogramsTogether(expHists[0],mcHists[i],1,expHists[0]->GetNbinsX());
        }
    }


    double origMax=0;
    for (int i = start;i< start+numExp+numMC;i++)
    {
        if(plotHists[i]->GetMaximum() >origMax) origMax=plotHists[i]->GetMaximum();
    }
    for (int i = 0;i< numExp;i++)
    {
        setDisplayQuantitiesForHist(i,origMax,expHists[i],detType,plotType,plotVarType);
    }

    for (int i = 0;i< numMC;i++)
    {
        setDisplayQuantitiesForHist(-i-1,origMax,mcHists[i],detType,plotType,plotVarType);
    }


    TH1** residHists;
    double max=-REALLY_BIG_DBL;
    double min=REALLY_BIG_DBL;
    double check;
    if(numExp==1)
    {
        ///Multiple MC
        residHists=new TH1*[numMC];
        for (int i = 0;i< numMC;i++)
        {
            histName="hist"+int2str(plotHists.size());
            residHists[i]=(TH1*) mcHists[i]->Clone(histName);
            plotHists.push_back(residHists[i]);
            residHists[i]->Reset();

            residHists[i]->Add(expHists[0],mcHists[i],1,-1);


            if(normalize)
            {
                residHists[i]->Divide(mcHists[i]);
            }

            check=GetIdealMaxForPlotting(residHists[i]);
            if(check > max) max=check;
            check=GetIdealMinForPlotting(residHists[i]);
            if(check < min) min=check;

        }

        if(normalize)
        {
            residHists[0]->GetYaxis()->SetTitle("Normalized Residual");
            residHists[0]->SetTitle(TString(plotHists[start]->GetTitle())+" Residuals #frac{Exp - MC}{MC}");
            if(detType==DET_EP)
            {
                residHists[0]->GetYaxis()->SetRangeUser(-.4,.4);
            }
            else
            {
                residHists[0]->GetYaxis()->SetRangeUser(-.5,.5);
            }
        }
        else
        {
            residHists[0]->SetTitle(TString(plotHists[start]->GetTitle())+" Residuals Exp - MC");
        }
        residHists[0]->GetYaxis()->SetRangeUser(min,max);

        drawHist(numMC, residHists, 0, nullptr);
    }
    else
    {
        ///Multiple Exp
        residHists=new TH1*[numExp];
        for (int i = 0;i< numExp;i++)
        {
            histName="hist"+int2str(plotHists.size());
            residHists[i]=(TH1*) expHists[i]->Clone(histName);
            plotHists.push_back(residHists[i]);
            residHists[i]->Reset();
            residHists[i]->Add(expHists[i],mcHists[0],1,-1);
            if(normalize)
            {
                residHists[i]->Divide(mcHists[0]);
            }
            check=GetIdealMaxForPlotting(residHists[i]);
            if(check > max) max=check;
            check=GetIdealMinForPlotting(residHists[i]);
            if(check < min) min=check;

        }
        if(normalize)
        {
            residHists[0]->SetTitle(TString(plotHists[start]->GetTitle())+" Residuals #frac{Exp - MC}{MC}");
            residHists[0]->GetYaxis()->SetTitle("Normalized Residual");
            if(detType==DET_EP)
            {
                residHists[0]->GetYaxis()->SetRangeUser(-.4,.4);
            }
            else
            {
                residHists[0]->GetYaxis()->SetRangeUser(-.5,.5);
            }
        }
        else
        {
            residHists[0]->SetTitle(TString(plotHists[start]->GetTitle())+" Residuals Exp - MC");
        }
        residHists[0]->GetYaxis()->SetRangeUser(min,max);



        drawHist(numExp, residHists, 0, nullptr);
    }

    delete[] residHists;
    delete[] expHists;
    delete[] mcHists;

}

void RDK2AnalysisPlotter::drawPlotDets(CoDet detType, PlotVarType plotVarType)
{
    gROOT->cd();

    int numDets;

    TString histTitleString;
    if(detType == DET_EPG)
    {
        numDets=12;
        histTitleString="Individual BGO Detectors;Detector Number;% Difference From Mean";
    }
    else
    {
        numDets=3;
        histTitleString="Individual BAPD Detectors;Detector Number;% Difference From Mean";
    }
    const int numHists=numExp+numMC;
    TH1* tempHist;
    double integral;
    double error;

    TH1** detsHists=new TH1*[numHists];
    TString histName;
    double max=-REALLY_BIG_DBL;
    double min=REALLY_BIG_DBL;

    for (int i = 0;i< numHists;i++)  //Set loop (either exp or mc)
    {
        histName="hist"+int2str(plotHists.size());
        detsHists[i]=new TH1D(histName,histTitleString,numDets,.5,numDets+.5);
        detsHists[i]->Sumw2();
        plotHists.push_back(detsHists[i]);

        for (int j = 0;j< numDets;j++)
        {
            if(i<numExp)
            {
                tempHist=getExpHist( i, detType, plotVarType,j+1);
            }
            else
            {
                tempHist=getMCHist( i-numExp, detType, plotVarType,j+1);
            }
            if(tempHist != nullptr)
            {
                integral = tempHist->IntegralAndError(1,tempHist->GetNbinsX(),error);
                detsHists[i]->SetBinContent(j+1,integral);
                detsHists[i]->SetBinError(j+1,error);
            }
        }

        if(i<numExp)
            setDisplayQuantitiesForHist(i,1,detsHists[i],detType,PLOT_DETS,plotVarType);
        else
            setDisplayQuantitiesForHist(-(i-numExp)-1,1,detsHists[i],detType,PLOT_DETS,plotVarType);

        double mean=detsHists[i]->Integral()/numDets; // Mean for the set
        detsHists[i]->Scale(100./mean);
        for (int j = 0;j< numDets;j++)
        {
            detsHists[i]->SetBinContent(j+1,detsHists[i]->GetBinContent(j+1)-100. );
        }
        if(detsHists[i]->GetMaximum() > max)
            max=detsHists[i]->GetMaximum()+detsHists[i]->GetBinError(detsHists[i]->GetMaximumBin());
        if(detsHists[i]->GetMinimum() < min)
            min=detsHists[i]->GetMinimum()-detsHists[i]->GetBinError(detsHists[i]->GetMinimumBin());
    }

    detsHists[0]->GetYaxis()->SetRangeUser(min-.05*abs(min),max+.05*abs(max));
    drawHist(numHists, detsHists, 0, nullptr);
    delete[] detsHists;

}


TH1* RDK2AnalysisPlotter::getExpHist(int i,CoDet detType, PlotVarType plotVarType, int det)
{
    TH1* outHist =nullptr;


    if(plotVarType==PLOTVAR_PT)
    {
        if(detType==DET_EP)
        {
            outHist= exp[i]->GetEP_pTHist();
        }
        else if(detType==DET_EPG)
        {
            outHist= exp[i]->GetEPG_pTHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= exp[i]->GetEPBG_pTHist();
        }
    }
    else if(plotVarType==PLOTVAR_EE)
    {
        if(detType==DET_EP)
        {
            outHist= exp[i]->GetEP_eEHist();
        }
        else if(detType==DET_EPG)
        {
            outHist= exp[i]->GetEPG_eEHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= exp[i]->GetEPBG_eEHist();
        }
    }
    else if(plotVarType==PLOTVAR_PE)
    {
        if(detType==DET_EP)
        {
            outHist= exp[i]->GetEP_pEHist();
        }
        else if(detType==DET_EPG)
        {
            outHist= exp[i]->GetEPG_pEHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= exp[i]->GetEPBG_pEHist();
        }
    }
    else if(plotVarType==PLOTVAR_GE)
    {
        if(det==0)
        {
            if(detType==DET_EPG)
            {
                outHist= exp[i]->GetEPG_gEAvgHist();
            }
            else if(detType==DET_EPBG)
            {
                outHist= exp[i]->GetEPBG_bGEAvgHist();
            }
        }
        else
        {
            if(detType==DET_EPG)
            {
                outHist= exp[i]->GetEPG_gEDetHists(det);
            }
            else if(detType==DET_EPBG)
            {
                outHist= exp[i]->GetEPBG_bGEDetHists(det);
            }
        }
    }
    else if(plotVarType==PLOTVAR_GEVAR)
    {
        if(detType==DET_EPG)
        {
            outHist= exp[i]->GetEPG_gEAvgVarHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= exp[i]->GetEPBG_bGEAvgVarHist();
        }
    }
    return outHist;
}

TH1* RDK2AnalysisPlotter::getMCHist(int i,CoDet detType, PlotVarType plotVarType,int det)
{
    TH1* outHist=nullptr;
    if(plotVarType==PLOTVAR_PT)
    {
        if(detType==DET_EP)
        {
            outHist= mc[i]->GetEP_pTHist();
        }
        else if(detType==DET_EPG)
        {
            outHist= mc[i]->GetEPG_pTHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= mc[i]->GetEPBG_pTHist();
        }
    }
    else if(plotVarType==PLOTVAR_EE)
    {
        if(detType==DET_EP)
        {
            outHist= mc[i]->GetEP_eEHist();
        }
        else if(detType==DET_EPG)
        {
            outHist= mc[i]->GetEPG_eEHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= mc[i]->GetEPBG_eEHist();
        }
    }
    else if(plotVarType==PLOTVAR_PE)
    {
        if(detType==DET_EP)
        {
            outHist= mc[i]->GetEP_pEHist();
        }
        else if(detType==DET_EPG)
        {
            outHist= mc[i]->GetEPG_pEHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= mc[i]->GetEPBG_pEHist();
        }

    }
    else if(plotVarType==PLOTVAR_GE)
    {
        if(det==0)
        {
            if(detType==DET_EPG)
            {
                outHist= mc[i]->GetEPG_gEAvgHist();
            }
            else if(detType==DET_EPBG)
            {
                outHist= mc[i]->GetEPBG_bGEAvgHist();
            }
        }
        else
        {
            if(detType==DET_EPG)
            {
                outHist= mc[i]->GetEPG_gEDetHists(det);
            }
            else if(detType==DET_EPBG)
            {
                outHist= mc[i]->GetEPBG_bGEDetHists(det);
            }
        }
    }
    else if(plotVarType==PLOTVAR_GEVAR)
    {
        if(detType==DET_EPG)
        {
            outHist= mc[i]->GetEPG_gEAvgVarHist();
        }
        else if(detType==DET_EPBG)
        {
            outHist= mc[i]->GetEPBG_bGEAvgVarHist();
        }
    }
    return outHist;
}

void RDK2AnalysisPlotter::setDisplayQuantitiesForHist(int i,double factor,TH1* inpHist, CoDet detType,PlotType plotType, PlotVarType plotVarType)
{
    int j=i;
    if(i<0)
    {
        j=-(i+1);
    }
    if(i>=0) //Exp
    {
        inpHist->SetLineColor(colorList[(j+4)%GRAPHICS_NUM_TYPES]);
        inpHist->SetMarkerColor(colorList[(j+4)%GRAPHICS_NUM_TYPES]);
        inpHist->SetMarkerStyle(markerList[j%GRAPHICS_NUM_TYPES]);
    }
    else  //MC
    {
        if(plotType==PLOT_COMP)
            inpHist->SetLineStyle(j%5+1);
        inpHist->SetLineColor(colorList[j%GRAPHICS_NUM_TYPES]);
        inpHist->SetMarkerColor(colorList[j%GRAPHICS_NUM_TYPES]);
        inpHist->SetMarkerStyle(markerList[j%GRAPHICS_NUM_TYPES]);
    }
    int rebin=1;
    double xLow=0,xHigh=1,yLow=0,yHigh=20;

    if(plotVarType==PLOTVAR_PT)
    {
        xLow=2e-6; xHigh=25e-6;
        if(detType==DET_EP)
        {
            rebin=3;
        }
        else
        {
            rebin=9;
        }
    }
    else if(plotVarType==PLOTVAR_EE)
    {
        xLow=50; xHigh=800;
        if(detType==DET_EP)
        {
            rebin=100;
        }
        else
        {
            rebin=200;
        }
    }
    else if(plotVarType==PLOTVAR_PE)
    {
        rebin=10; xLow=10; xHigh=32;
    }
    else if(plotVarType==PLOTVAR_GE && plotType != PLOT_DETS)
    {
        if(detType == DET_EPG)
        {
            rebin=100; xLow=10; xHigh=800;
        }
        else
        {
            rebin=10; xLow=.1; xHigh=20;
        }
    }

    inpHist->Rebin(rebin);

    if(plotType==PLOT_COMP || plotType == PLOT_DETS)
    {
        inpHist->GetYaxis()->SetTitleOffset(1.7);
    }
    else if(plotType==PLOT_RESID)
    {
        inpHist->GetYaxis()->SetTitleOffset(1.2);
        if(detType==DET_EP)
        {
            yLow=-0.4*factor; yHigh=0.4*factor;
        }
        else
        {
            yLow=-0.5*factor; yHigh=0.5*factor;
        }
        inpHist->GetYaxis()->SetRangeUser(yLow,yHigh);
    }
    else if(plotType==PLOT_NORMRESID)
    {
        inpHist->GetYaxis()->SetTitleOffset(1.2);
        if(detType==DET_EP)
        {
            yLow=-0.4; yHigh=0.4;
        }
        else
        {
            yLow=-0.5; yHigh=0.5;
        }
//          inpHist->GetYaxis()->SetRangeUser(yLow,yHigh);
    }

    if(plotType != PLOT_DETS)
        inpHist->GetXaxis()->SetRangeUser(xLow,xHigh);
    else
    {
        if(detType==DET_EPG)
        {
            inpHist->GetXaxis()->SetRangeUser(-0.5,12.5);
        }
        else
        {
            inpHist->GetXaxis()->SetRangeUser(-0.5,3.5);
        }
        inpHist->GetYaxis()->SetTitleOffset(1.0);
    }

    if(plotType==PLOT_NORMRESID && detType != DET_EP)
    {
        inpHist->GetXaxis()->SetMoreLogLabels(true);
        inpHist->GetXaxis()->SetTitleOffset(1.2);
    }



    inpHist->SetLineWidth(1);
    inpHist->SetMarkerSize(1);



    return;
}

void createExpAnalysisAndPlot(TString inpName, TString inpTitle, TString inpPIDString,TString inpGroupString, RDK2CutSet inpCutSet,const double inpBGOWindows[6],const double inpBAPDWindows[6])
{
    TString expFilePath;
    TFile* expFile;
    RDK2ExpAnalyzer* expAnalysis;
    expFilePath=ANALYSIS_DIR;
    expFilePath+=inpName+".root";
    expFile=new TFile(expFilePath,"recreate");
    expFile->cd();
    expAnalysis=new RDK2ExpAnalyzer(inpName,inpTitle,inpPIDString,inpGroupString,inpCutSet,inpBGOWindows,inpBAPDWindows);
    expAnalysis->MakeAllStandardHists();
    expFile->cd();
    expAnalysis->Write();
    expFile->Close();
    delete expAnalysis;
    delete expFile;
    RDK2AnalysisPlotter comPlot(inpTitle,inpName,"");
    comPlot.makeAll();
}

void createMCAnalysisAndPlot(TString inpName, TString inpTitle, int num3, TString p3RID, TString e3RID, int num4,TString p4RID, TString e4RID,TString g4RID, RDK2CutSet inpCutSet)
{
    TString mcFilePath;
    TFile* mcFile;
    RDK2MCAnalyzer* mcAnalysis;
    mcFilePath=ANALYSIS_DIR;
    mcFilePath+=inpName+".root";
    mcFile=new TFile(mcFilePath,"recreate");
    mcFile->cd();
    mcAnalysis=new RDK2MCAnalyzer( inpName,  inpTitle,  num3,  p3RID,  e3RID,  num4, p4RID, e4RID, g4RID,  inpCutSet);
    mcAnalysis->MakeAllStandardHists();
    mcFile->cd();
    mcAnalysis->Write();
    mcFile->Close();
    delete mcAnalysis;
    delete mcFile;
    RDK2AnalysisPlotter comPlot(inpTitle,"",inpName);
    comPlot.makeAll();

}

//
//#include "mattexpanalysis.h"
//void createAndPlotPESections(TString nameString,RDK2Set* inpThreeBodySet, int numSections, double begin, double end,int pOption, bool useExp)
//{
//    TString histName,cutString,histFilePath;
//    TString HistsDir=HISTS_DIR;
//    TH1* tempHist;
//    double pESections[numSections+1];
//    for (int i = 0;i< numSections+1;i++)
//    {
//        pESections[i]=begin+i*(end-begin)/numSections;
//    }
//
//    HistDim peDim={60,begin,end};
//    TString pEString="pChn.SBDEDepTotal";
//    TString histKicker="";
//    TString analType="MC";
//    if(pOption==1 && !useExp)
//    {
//        pEString="pDChn.SBDEDepBlur";
//        histKicker="WRes";
//    }
//    else if(pOption==2 && !useExp)
//    {
//        pEString="pChn.SBDKEFirst";
//        histKicker="KEFirst";
//    }
//    else if(useExp)
//    {
//        pEString="SBDEDepP";
//        analType="Exp";
//    }
//    int numSeries=NUM_EXP_SERIES_G7_Cycle5_8;
//    for (int i = 0;i< numSections;i++)
//    {
//        double lowEn=pESections[i];
//        double highEn=pESections[i+1];
//        histName="Hist_p_En_"+int2str(i);
//        if(useExp)
//        {
//            cutString="SBDEDepE > 50 && SBDEDepE < 800 && SBDToFP > 2e-6 && SBDToFP < 10e-6 && SBDEDepP > "+d2str(lowEn)+" && SBDEDepP < "+d2str(highEn);
//        }
//        else
//        {
//            cutString="eDChn.SBDEDepBlur > 50. && eDChn.SBDEDepBlur < 800. && "+pEString+" > "+d2str(lowEn)+" && "+pEString+" < "+d2str(highEn)+" && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 10.e-6";
//        }
//        cout << cutString << endl;
//        TCut theCut(cutString);
//         //pE
//
//    	if(useExp) tempHist=makeExpEPHist("PID8", numSeries, EXP_SERIES_ID_LIST_G7_Cycle5_8,pEString, theCut, peDim);
//        else tempHist=inpThreeBodySet->makeHist("temp",pEString,theCut,peDim);
//
//        histFilePath=HistsDir+analType+histName+"_Ep_EP"+histKicker+"_"+inpThreeBodySet->GetPResultsID()+"_"+getCurrentDateString()+".txt";
//        convertTH1ToTXT(tempHist,histFilePath);
//        delete tempHist;
//
//        //eE
//    	if(useExp) tempHist=makeExpEPHist("PID8", numSeries, EXP_SERIES_ID_LIST_G7_Cycle5_8,"SBDEDepE", theCut, EDEPE_DIM);
//        else tempHist=inpThreeBodySet->makeHist("temp","eDChn.SBDEDepBlur",theCut,EDEPE_DIM);
//        histFilePath=HistsDir+analType+histName+"_Ee_EP"+histKicker+"_"+inpThreeBodySet->GetPResultsID()+"_"+getCurrentDateString()+".txt";
//        convertTH1ToTXT(tempHist,histFilePath);
//        delete tempHist;
//
//        //pToF
//    	if(useExp) tempHist=makeExpEPHist("PID8", numSeries, EXP_SERIES_ID_LIST_G7_Cycle5_8,"SBDToFP/(40.0e-9)", theCut, TOFPCCHAN_DIM);
//        else tempHist=inpThreeBodySet->makeHist("temp","pChn.SBDTimeFirst/(40.0e-9)",theCut,TOFPCCHAN_DIM);
//        histFilePath=HistsDir+analType+histName+"_pToF_EP"+histKicker+"_"+inpThreeBodySet->GetPResultsID()+"_"+getCurrentDateString()+".txt";
//        convertTH1ToTXT(tempHist,histFilePath);
//        delete tempHist;
//    }
//
//    TH1* mcPeHists[numSections];
//    TH1* mcEeHists[numSections];
//    TH1* mcTHists[numSections];
//    int pEColorList[6]= {kBlack, kRed, kGreen, kBlue, kViolet, kGray};
//    for (int i = 0;i< numSections;i++)
//    {
//        histName="Hist_p_En_"+int2str(i);
//
//        mcPeHists[i]=getTabSeperatedHist(HistsDir+analType+histName+"_Ep_EP"+histKicker+"_"+inpThreeBodySet->GetPResultsID()+"_"+getCurrentDateString()+".txt",peDim);
//        mcEeHists[i]=getTabSeperatedHist(HistsDir+analType+histName+"_Ee_EP"+histKicker+"_"+inpThreeBodySet->GetPResultsID()+"_"+getCurrentDateString()+".txt",EDEPE_DIM);
//        mcTHists[i]=getTabSeperatedHist(HistsDir+analType+histName+"_pToF_EP"+histKicker+"_"+inpThreeBodySet->GetPResultsID()+"_"+getCurrentDateString()+".txt",TOFPCCHAN_DIM);
//
//        mcPeHists[i]->SetLineColor(pEColorList[i]);
//        mcEeHists[i]->SetLineColor(pEColorList[i]);
//        mcTHists[i]->SetLineColor(pEColorList[i]);
//
//        mcEeHists[i]->Rebin(10);
//        mcTHists[i]->Rebin(4);
//
//        mcPeHists[i]->SetTitle(nameString+";Proton Energy (keV);Counts");
//        mcEeHists[i]->SetTitle(";Electron Energy (keV);Counts");
//        mcTHists[i]->SetTitle(";e-p TOF (1ch=40ns);Counts");
//
//        mcPeHists[i]->GetXaxis()->SetRangeUser(begin,end);
//    }
//
//    TCanvas* theCanvas = new TCanvas("c","c",1024,2000);
//    theCanvas->Divide(1,3);
//
//    TLegend* theLegends[3];
//    TH1** theHists[3]={mcPeHists,mcEeHists,mcTHists};
//    for (int j = 0;j< 3;j++) //graph in layout loop
//    {
//        double max=0;
//        theCanvas->cd(j+1);
//        theLegends[j] = new TLegend(0.8,0.5,0.93,0.9);
//        theLegends[j]->SetBorderSize(1);
//        theLegends[j]->SetFillColor(kWhite);
//        for (int i = 0;i< numSections;i++)  //PEnergy group loop
//        {
//            double lowEn=pESections[i];
//            double highEn=pESections[i+1];
//            histName="Hist_p_En"+d2str(lowEn,2,1)+"_"+d2str(highEn,2,1);
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
//    TString imagePath=TString(GRAPHS_DIR)+analType+"_"+nameString+"sections_normalized_"+getCurrentDateString()+".pdf";
//    theCanvas->Print(imagePath);
//}

//void RDK2AnalysisPlotter::determineExtremeCutSet()
//{
//    extremeCutSet.SetEELow(REALLY_BIG_DBL);
//    extremeCutSet.SetPELOW(REALLY_BIG_DBL);
//    extremeCutSet.SetPTLOW(REALLY_BIG_DBL);
//    extremeCutSet.SetGELOW(REALLY_BIG_DBL);
//    extremeCutSet.SetEELow(REALLY_BIG_DBL);
//    extremeCutSet.SetPELOW(REALLY_BIG_DBL);
//    extremeCutSet.SetPTLOW(REALLY_BIG_DBL);
//    extremeCutSet.SetGELOW(REALLY_BIG_DBL);
//    //Check each exp
//    for (int i = 0;i< numExp;i++)
//    {
//        if(e
//
//    }
//
//}
