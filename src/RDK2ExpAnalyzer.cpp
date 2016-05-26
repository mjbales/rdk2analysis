#include "RDK2ExpAnalyzer.h"

#include "MRMisc.h"
#include "MRIO.h"
#include "MRText.h"
#include <iostream>
#include <cmath>
#include <time.h>

#include "RDK2IO.h"

using namespace std;

ClassImp(RDK2ExpAnalyzer);

RDK2ExpAnalyzer::RDK2ExpAnalyzer()
{
    SetName("RDK2ExpAnalysis");
    epCounts=0;
    epgCounts=0;
    epgPerEP=0;
    epgPerEPError=0;
    epbgCounts=0;
    epbgPerEP=0;
    epbgPerEPError=0;
    numLiveBGO=0;
    numLiveBAPD=0;
    NULLHists();
    expChain=nullptr;
    chainLoaded=false;
    for (int i = 0;i< Exp_NUMBGO;i++)
    {
    	epCountsForLiveBGO[i]=0;
    }
    for (int i = 0;i< Exp_NUMBAPD;i++)
    {
    	epCountsForLiveBAPD[i]=0;
    }
    averageEPCountsForLiveBGO=0;
    averageEPCountsForLiveBAPD=0;
}

RDK2ExpAnalyzer::~RDK2ExpAnalyzer()
{
    ResetHists();
}

void RDK2ExpAnalyzer::NULLHists()
{
    ep_eEHist=nullptr;
    ep_pEHist=nullptr;
    ep_pTHist=nullptr;
    epg_eEHist=nullptr;
    epg_pEHist=nullptr;
    epg_pTHist=nullptr;
    epbg_eEHist=nullptr;
    epbg_pEHist=nullptr;
    epbg_pTHist=nullptr;
    epg_gEAvgHist=nullptr;
    epg_gEAvgVarHist=nullptr;
    for (int i = 0;i< Exp_NUMBGO;i++)
    {
    	epg_gEDetHists[i]=nullptr;
    	epg_gEDetVarHists[i]=nullptr;
    }

    for (int i = 0;i< Exp_NUMBAPD;i++)
    {
    	epbg_bGEDetHists[i]=nullptr;
    	epbg_bGEDetVarHists[i]=nullptr;
    }
    epbg_bGEAvgHist=nullptr;
    epbg_bGEAvgVarHist=nullptr;
}

void RDK2ExpAnalyzer::ResetHists()
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
    for (int i = 0;i< Exp_NUMBGO;i++)
    {
        delete epg_gEDetHists[i];
        delete epg_gEDetVarHists[i];
    }
    for (int i = 0;i< Exp_NUMBAPD;i++)
    {
        delete epbg_bGEDetHists[i];
        delete epbg_bGEDetVarHists[i];
    }


    delete epbg_bGEAvgHist;
    delete epbg_bGEAvgVarHist;
    NULLHists();
}

RDK2ExpAnalyzer::RDK2ExpAnalyzer(TString inpName, TString inpTitle, TString inpPIDString,TString inpGroupString, RDK2CutSet inpCutSet,const double inpBGOWindows[6],const double inpBAPDWindows[6])
{
    NULLHists();
    numLiveBGO=0;
    numLiveBAPD=0;
    epCounts=0;
    epgCounts=0;
    epgPerEP=0;
    epgPerEPError=0;
    epbgCounts=0;
    epbgPerEP=0;
    epbgPerEPError=0;
    averageEPCountsForLiveBGO=0;
    averageEPCountsForLiveBAPD=0;
    for (int i = 0;i< Exp_NUMBGO;i++)
    {
    	epCountsForLiveBGO[i]=0;
    }
    for (int i = 0;i< Exp_NUMBAPD;i++)
    {
    	epCountsForLiveBAPD[i]=0;
    }
    chainLoaded=false;
    expChain=nullptr;
    SetName(inpName);
    SetTitle(inpTitle);

    PIDString=inpPIDString;
    groupString=inpGroupString;

    for (int i = 0;i< 6;i++)
    {
    	bgoWindows[i]=inpBGOWindows[i];
    	bapdWindows[i]=inpBAPDWindows[i];
    }

    theCutSet=inpCutSet;
}

void RDK2ExpAnalyzer::loadEventChain()
{
    int* seriesArray=nullptr;
    int numSeries=getSeriesArray(groupString, &seriesArray);
    expChain=makeExpChain(PIDString,numSeries,seriesArray);
    chainLoaded=true;
    delete[] seriesArray;
}

TH1* RDK2ExpAnalyzer::MakeHist(TString histName,TString drawString,TCut inpCut,HistDim inpHistDimX, HistDim inpHistDimY,HistDim inpHistDimZ)
{
    gROOT->cd();
    TH1* histPrototype;
    TH1* outHist;

    if(inpHistDimZ.numBins > 0)
        histPrototype=makeTH3DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY,inpHistDimZ);
    else if(inpHistDimY.numBins > 0)
        histPrototype=makeTH2DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY);
    else
        histPrototype=makeTH1DFromDim("prototype"+histName,"",inpHistDimX);


    outHist=MakeHist(histName,drawString,inpCut,histPrototype);

    delete histPrototype;
    return outHist;

}

TH1* RDK2ExpAnalyzer::MakeHist(TString histName,TString drawString,TCut inpCut,TH1* histPrototype)
{
    if(!chainLoaded)
    {
        loadEventChain();
    }
    gROOT->cd();
    TH1* outHist=(TH1*) histPrototype->Clone(histName+GetName());
    outHist->Reset();
    outHist->Sumw2();
     cout << "Drawing: " << drawString << endl;
    cout << "using cuts: " << inpCut.GetTitle() << endl;
    expChain->Draw(drawString+" >> "+histName+GetName(),inpCut);
    cout << "Events drawn: " << outHist->Integral() << endl;
    cout << "Hist Mean: " << outHist->GetMean() << endl;
    return outHist;
}

TH1* RDK2ExpAnalyzer::MakeRadHist(TString histName,CoDet detType,TString drawString,int peakOption,int det,TH1* histPrototype,bool scaleByEPCounts)
{
    TH1* outHist=nullptr;
    int* countsToUseArray=nullptr;
    TString detString;
    if(detType==DET_EPG)
    {
        countsToUseArray=epCountsForLiveBGO;
        detString="BGO";
    }
    else if(detType==DET_EPBG)
    {
        countsToUseArray=epCountsForLiveBAPD;
        detString="BAPD";
    }

    if(countsToUseArray[det-1] == 0) //Remember det starts at 1 and array starts at 0
    {
        CalcEPCountsForLiveDets(detType);
    }

    //If detector is dead return a blank hist
    if(countsToUseArray[det-1] == -1)
    {
        cout << detString << det << " is dead, no hist made" << endl;
        outHist=(TH1*) histPrototype->Clone(histName);
        outHist->Reset();
        outHist->Sumw2();
        return outHist;
    }

    if(detType==DET_EPG)
    {
        if(drawString=="")
        {
            drawString+="BGOEDep"+int2str(det);
        }
        outHist=MakeHist(histName,drawString,theCutSet.CreateExpEPGCut(0,0,0,det,detType)&&CreateWindowCut(detType,peakOption,det), histPrototype);
        if(scaleByEPCounts)
            outHist->Scale(1./countsToUseArray[det-1]);
    }
    else if(detType==DET_EPBG)
    {
        if(drawString=="")
        {
            drawString+="BAPDEDep"+int2str(det);
        }

        outHist=MakeHist(histName,drawString,theCutSet.CreateExpEPGCut(0,0,0,det,detType)&&CreateWindowCut(detType,peakOption,det), histPrototype);
        if(scaleByEPCounts)
            outHist->Scale(1./countsToUseArray[det-1]);
    }


    return outHist;
}

TH1* RDK2ExpAnalyzer::MakeRadHist(TString histName,CoDet detType,TString drawString,int peakOption,int det,HistDim inpHistDimX, HistDim inpHistDimY,HistDim inpHistDimZ,bool scaleByEPCounts)
{
    gROOT->cd();
    TH1* histPrototype;
    TH1* outHist;

    if(inpHistDimZ.numBins > 0)
        histPrototype=makeTH3DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY,inpHistDimZ);
    else if(inpHistDimY.numBins > 0)
        histPrototype=makeTH2DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY);
    else
        histPrototype=makeTH1DFromDim("prototype"+histName,"",inpHistDimX);

    outHist=MakeRadHist(histName,detType,drawString,peakOption,det,histPrototype, scaleByEPCounts);
    delete histPrototype;
    return outHist;
}

TH1* RDK2ExpAnalyzer::MakeRadAvgHist(TString histName,CoDet detType,TString drawString,int peakOption,TH1* histPrototype,bool scaleByEPCounts)
{
    int numDets=0;
    int numLiveDets=0;
    if(detType==DET_EPG)
    {
        numDets=Exp_NUMBGO;
        if(numLiveBGO==0)
        {
            CalcEPCountsForLiveDets(detType);
        }
        numLiveDets=numLiveBGO;
    }
    else if(detType==DET_EPBG)
    {
        numDets=Exp_NUMBAPD;
        if(numLiveBAPD==0)
        {
            CalcEPCountsForLiveDets(detType);
        }
        numLiveDets=numLiveBAPD;
    }
    TH1* outHist=(TH1*) histPrototype->Clone(histName);
    outHist->Reset();
    outHist->Sumw2();
    for (int i = 0;i< numDets;i++)
    {
    	TH1* tempHist = MakeRadHist("radIndivTemp"+int2str(i+1),detType,drawString,peakOption,i+1,histPrototype, scaleByEPCounts);
    	outHist->Add(tempHist,1./numLiveDets);
        delete tempHist;
    }
    return outHist;
}

TH1* RDK2ExpAnalyzer::MakeRadAvgHist(TString histName,CoDet detType,TString drawString,int peakOption,HistDim inpHistDimX, HistDim inpHistDimY,HistDim inpHistDimZ,bool scaleByEPCounts)
{
    gROOT->cd();
    TH1* histPrototype;
    TH1* outHist;

    if(inpHistDimZ.numBins > 0)
        histPrototype=makeTH3DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY,inpHistDimZ);
    else if(inpHistDimY.numBins > 0)
        histPrototype=makeTH2DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY);
    else
        histPrototype=makeTH1DFromDim("prototype"+histName,"",inpHistDimX);

    outHist=MakeRadAvgHist(histName,detType,drawString,peakOption,histPrototype, scaleByEPCounts);
    delete histPrototype;
    return outHist;
}


TH1* RDK2ExpAnalyzer::MakeBGSubHist(TString histName,CoDet detType,TString drawString,int det,HistDim inpHistDimX, HistDim inpHistDimY,HistDim inpHistDimZ,bool scaleByEPCounts)
{
    gROOT->cd();
    TH1* histPrototype;
    TH1* outHist;

    if(inpHistDimZ.numBins > 0)
        histPrototype=makeTH3DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY,inpHistDimZ);
    else if(inpHistDimY.numBins > 0)
        histPrototype=makeTH2DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY);
    else
        histPrototype=makeTH1DFromDim("prototype"+histName,"",inpHistDimX);

    outHist=MakeBGSubHist(histName,detType,drawString,det,histPrototype, scaleByEPCounts);
    delete histPrototype;
    return outHist;
}

TH1* RDK2ExpAnalyzer::MakeBGSubHist(TString histName,CoDet detType,TString drawString, int det,TH1* histPrototype,bool scaleByEPCounts)
{

    TH1* prePeakHist=MakeRadHist("prePeakHist",detType,drawString,0,det,histPrototype, scaleByEPCounts);
    TH1* onPeakHist=MakeRadHist("onPeakHist",detType,drawString,1,det,histPrototype, scaleByEPCounts);
//    TH1* postPeakHist=MakeRadHist("postPeakHist",detType,drawString,2,det,histPrototype);
    TH1* outHist=(TH1*) histPrototype->Clone(histName);
    outHist->Reset();
    outHist->Sumw2();
    outHist->Add(onPeakHist);
    outHist->Add(prePeakHist,-GetPrePeakScaleFactor(detType));
    delete prePeakHist;
    delete onPeakHist;
    return outHist;

}

TH1* RDK2ExpAnalyzer::MakeBGSubAvgHist(TString histName,CoDet detType,TString drawString,TH1* histPrototype,bool scaleByEPCounts)
{
    int numDets=0;
    int numLiveDets=0;
    if(detType==DET_EPG)
    {
        numDets=Exp_NUMBGO;
        if(numLiveBGO==0)
        {
            CalcEPCountsForLiveDets(detType);
        }
        numLiveDets=numLiveBGO;
    }
    else if(detType==DET_EPBG)
    {
        numDets=Exp_NUMBAPD;
        if(numLiveBAPD==0)
        {
            CalcEPCountsForLiveDets(detType);
        }
        numLiveDets=numLiveBAPD;
    }

    TH1* outHist=(TH1*) histPrototype->Clone(histName);
    outHist->Reset();
    outHist->Sumw2();

    for (int i = 0;i< numDets;i++)
    {
    	TH1* tempHist = MakeBGSubHist("radIndivTemp"+int2str(i+1),detType,drawString,i+1,histPrototype, scaleByEPCounts);
    	outHist->Add(tempHist,1./numLiveDets);
        delete tempHist;
    }
    return outHist;
}

TH1* RDK2ExpAnalyzer::MakeBGSubAvgHist(TString histName,CoDet detType,TString drawString,HistDim inpHistDimX, HistDim inpHistDimY,HistDim inpHistDimZ,bool scaleByEPCounts)
{
    gROOT->cd();
    TH1* histPrototype;
    TH1* outHist;

    if(inpHistDimZ.numBins > 0)
        histPrototype=makeTH3DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY,inpHistDimZ);
    else if(inpHistDimY.numBins > 0)
        histPrototype=makeTH2DFromDim("prototype"+histName,"",inpHistDimX,inpHistDimY);
    else
        histPrototype=makeTH1DFromDim("prototype"+histName,"",inpHistDimX);

    outHist=MakeBGSubAvgHist(histName,detType,drawString,histPrototype, scaleByEPCounts);
    delete histPrototype;
    return outHist;
}

TH1* RDK2ExpAnalyzer::MakeStandardHist(TString histName,CoDet detType,TString drawString,HistDim inpHistDimX)
{
    TH1* outHist;
    if(detType==DET_EPBG || detType==DET_EPG)
    {
        outHist=MakeBGSubAvgHist(histName,detType,drawString,inpHistDimX);
    }
    else
    {
        outHist=MakeHist(histName,drawString,theCutSet.CreateExpEPCut(),inpHistDimX);

    }
    return outHist;
}

void RDK2ExpAnalyzer::MakeAllStandardHists()
{
    cout << "Making all standard histograms for: " << GetName() << endl;
    time_t startTime,endTime;
    time (&startTime);  //Start Time for Run

//    CalcRates();  //Calculate EP and EPG/EPBG rates

    TString HistsDir=HISTS_DIR;
    TString baseHistName=TString("Hist_Exp_")+GetName()+"_";
    TString endHistName=".txt";

    TH1** eEHists[3]={&ep_eEHist,&epg_eEHist,&epbg_eEHist};
    TH1** pEHists[3]={&ep_pEHist,&epg_pEHist,&epbg_pEHist};
    TH1** pTHists[3]={&ep_pTHist,&epg_pTHist,&epbg_pTHist};


    CoDet detType;
    TString coString,yAxisString;
    int numDets=0;
    for (int i = 0;i< 2;i++)  //No BAPDs till we can do BAPDs
    {
        if(i==0)
        {
            cout << "=========================================================" << endl;
            cout << "Making EP Hists for " << GetName() << endl;
            cout << "=========================================================" << endl;
            detType=DET_EP;
            coString="EP";
            yAxisString="Counts";
        }
        else if(i==1)
        {
            cout << "=========================================================" << endl;
            cout << "Making EPG Hists for " << GetName() << endl;
            cout << "=========================================================" << endl;
            detType=DET_EPG;
            coString="EPG";
            numDets=Exp_NUMBGO;
            yAxisString="#frac{epg}{ep}";
        }
        else
        {
            cout << "=========================================================" << endl;
            cout << "Making EPBG Hists for " << GetName() << endl;
            cout << "=========================================================" << endl;
            detType=DET_EPBG;
            coString="EPBG";
            numDets=Exp_NUMBAPD;
            yAxisString="#frac{epbg}{ep}";
        }

        *(pTHists[i])=MakeStandardHist(coString+"_pTHist",detType,"SBDToFP",TOFP_DIM2);
        (*(pTHists[i]))->SetTitle("Proton Time of Flight "+coString+";Proton Time of Flight (s);"+yAxisString);
        convertTH1ToTXT(*(pTHists[i]),HistsDir+baseHistName+coString+"_ToFP"+endHistName);
        if(i==0) epCounts=(*(pTHists[i]))->Integral();


        *(pEHists[i])=MakeStandardHist(coString+"_pEHist",detType,"SBDEDepP",EDEPP_DIM);
        (*(pEHists[i]))->SetTitle("Proton Energy "+coString+";SBD Energy Channel (782=782 keV);"+yAxisString);
        convertTH1ToTXT(*(pEHists[i]),HistsDir+baseHistName+coString+"_EDepP"+endHistName);


        *(eEHists[i])=MakeStandardHist(coString+"_eEHist",detType,"SBDEDepE",EDEPE_DIM);
        (*(eEHists[i]))->SetTitle("Electron Energy "+coString+";SBD Energy Channel (782=782 keV);"+yAxisString);
        convertTH1ToTXT(*(eEHists[i]),HistsDir+baseHistName+coString+"_EDepE"+endHistName);


        if(i==1)
        {
            TH1* unscaledHist=makeTH1DFromDim("unscaledHist","unscaledHist",BGO_DIM);
            unscaledHist->Sumw2();
            epg_gEAvgHist=makeTH1DFromDim("epg_gEAvgHist","epg_gEAvgHist",BGO_DIM);
            epg_gEAvgHist->Sumw2();
            epg_gEAvgHist->SetTitle("Average BGO Gamma Energy "+coString+";BGO Energy Channel (511=511 keV);"+yAxisString);
            epg_gEAvgVarHist=new TH1D("epg_gEAvgVarHist","epg_gEAvgHist",20,BGO_VARHIST_LOWEDGES);
            epg_gEAvgVarHist->SetTitle("Average BGO Gamma Energy "+coString+";BGO Energy Channel (511=511 keV);"+yAxisString);
            for (int j = 0;j< numDets;j++)
            {
                epg_gEDetHists[j]=MakeBGSubHist("epg_gEDet"+int2str(j+1)+"Hist",detType,"",j+1,BGO_DIM);
                epg_gEDetHists[j]->SetTitle("BGO Det "+int2str(j+1)+" Gamma Energy "+coString+";BGO Energy Channel (511=511 keV);"+yAxisString);
                epg_gEDetVarHists[j]=MakeBGSubHist("epg_gEDet"+int2str(j+1)+"VarHist",detType,"",j+1,epg_gEAvgVarHist);
                epg_gEDetVarHists[j]->SetTitle("BGO Det "+int2str(j+1)+" Gamma Energy "+coString+";BGO Energy Channel (511=511 keV);"+yAxisString);
                convertTH1ToTXT(epg_gEDetHists[j],HistsDir+baseHistName+coString+"_EDepGDet"+int2str(j+1)+endHistName);
                convertTH1ToTXT(epg_gEDetVarHists[j],HistsDir+TString("VarHist_Exp_")+GetName()+"_"+coString+"_EDepGDet"+int2str(j+1)+endHistName);
                epg_gEAvgHist->Add(epg_gEDetHists[j],1./numLiveBGO);
                epg_gEAvgVarHist->Add(epg_gEDetVarHists[j],1./numLiveBGO);
                unscaledHist->Add(epg_gEAvgHist,epCountsForLiveBGO[j]);
            }
            convertTH1ToTXT(epg_gEAvgHist,HistsDir+baseHistName+coString+"_EDepGAvg"+endHistName);
            convertTH1ToTXT(epg_gEAvgVarHist,HistsDir+TString("VarHist_Exp_")+GetName()+"_"+coString+"_EDepGAvg"+endHistName);
            epgPerEP=epg_gEAvgHist->IntegralAndError(1,epg_gEAvgHist->GetNbinsX(),epgPerEPError);
            epgCounts=unscaledHist->IntegralAndError(1,unscaledHist->GetNbinsX(),epgCountsError);
            delete unscaledHist;
        }
        else if(i==2)
        {
            TH1* unscaledHist=makeTH1DFromDim("unscaledHist","unscaledHist",BAPD_DIM);
            unscaledHist->Sumw2();
            epbg_bGEAvgHist=makeTH1DFromDim("epbg_bGEAvgHist","epbg_bGEAvgHist",BAPD_DIM);
            epbg_bGEAvgHist->Sumw2();
            epbg_bGEAvgHist->SetTitle("Average BAPD Gamma Energy "+coString+";BAPD Energy Channel (5.9=5.9 keV);"+yAxisString);
            epbg_bGEAvgVarHist=new TH1D("epbg_bGEAvgVarHist","epbg_bGEAvgHist",20,BAPD_VARHIST_LOWEDGES);
            epbg_bGEAvgVarHist->SetTitle("Average BAPD Gamma Energy "+coString+";BAPD Energy Channel (5.9=5.9 keV);"+yAxisString);
            for (int j = 0;j< numDets;j++)
            {
                epbg_bGEDetHists[j]=MakeBGSubHist("epbg_bGEDet"+int2str(j+1)+"Hist",detType,"",j+1,BAPD_DIM);
                epbg_bGEDetHists[j]->SetTitle("BGO Det "+int2str(j+1)+" Gamma Energy "+coString+";BAPD Energy Channel (5.9=5.9 keV);"+yAxisString);
                epbg_bGEDetVarHists[j]=MakeBGSubHist("epbg_bGEDet"+int2str(j+1)+"VarHist",detType,"",j+1,epbg_bGEAvgVarHist);
                epbg_bGEDetVarHists[j]->SetTitle("BGO Det "+int2str(j+1)+" Gamma Energy "+coString+";BAPD Energy Channel (5.9=5.9 keV);"+yAxisString);
                convertTH1ToTXT(epbg_bGEDetHists[j],HistsDir+baseHistName+coString+"_EDepBGDet"+int2str(j+1)+endHistName);
                convertTH1ToTXT(epbg_bGEDetVarHists[j],HistsDir+TString("VarHist_Exp_")+GetName()+"_"+coString+"_EDepBGDet"+int2str(j+1)+endHistName);
                epbg_bGEAvgHist->Add(epbg_bGEDetHists[j],1./numLiveBAPD);
                epbg_bGEAvgVarHist->Add(epbg_bGEDetVarHists[j],1./numLiveBAPD);
                unscaledHist->Add(epg_gEAvgHist,epCountsForLiveBGO[j]);
            }

            convertTH1ToTXT(epbg_bGEAvgHist,HistsDir+baseHistName+coString+"_EDepBGAvg"+endHistName);
            convertTH1ToTXT(epbg_bGEAvgVarHist,HistsDir+TString("VarHist_Exp_")+GetName()+"_"+coString+"_EDepBGAvg"+endHistName);
            epbgPerEP=epbg_bGEAvgVarHist->IntegralAndError(1,epbg_bGEAvgVarHist->GetNbinsX(),epbgPerEPError);
            epgCounts=unscaledHist->IntegralAndError(1,unscaledHist->GetNbinsX(),epgCountsError);
            delete unscaledHist;
        }

    }

    PrintRates();  //Print the rates at the end
    time (&endTime);  //End Time for Run

}


int RDK2ExpAnalyzer::CalcEPCounts()
{
    HistDim simpleDim={3,-10,800};
    TH1* epRateHist=MakeHist("epRateHist","SBDEDepP",theCutSet.CreateExpEPCut(),simpleDim);
    epCounts=epRateHist->Integral(1,3);
    delete epRateHist;
    return epCounts;
}

double RDK2ExpAnalyzer::CalcEPGCounts(CoDet detType,double& error)
{
    double answer;
    HistDim simpleDim={3,-10,800};
    TH1* epgRateHist=MakeBGSubAvgHist("epgRateHist",detType,"",simpleDim,BLANK_DIM,BLANK_DIM,false);
    epgCounts=answer=epgRateHist->IntegralAndError(1,3,error);
    epgCountsError=error;
    delete epgRateHist;
    return answer;
}


double RDK2ExpAnalyzer::CalcEPGPerEP(CoDet detType,double& error) //Presuming error on EP is small
{
    double answer;
    HistDim simpleDim={3,-10,800};
    TH1* epgRateHist=MakeBGSubAvgHist("epgRateHist",detType,"",simpleDim,BLANK_DIM,BLANK_DIM,true); //Can't simply divide as we need to deal with a variable quantity of live detectors
    answer=epgPerEP=epgRateHist->IntegralAndError(1,3,error);
    epgPerEPError=error;
    delete epgRateHist;
    return answer;
}

void RDK2ExpAnalyzer::CalcRates()
{
    double tempError;
    CalcEPCounts();
    CalcEPCountsForLiveDets(DET_EPG);
    CalcEPGCounts(DET_EPG,tempError);
    CalcEPGPerEP(DET_EPG,tempError);
//    CalcEPCountsForLiveDets(DET_EPBG);
//    CalcEPGCounts(DET_EPBG,tempError);
//    CalcEPGPerEP(DET_EPBG,tempError);
}

void RDK2ExpAnalyzer::PrintRates()
{
    cout << GetName() << " Rates:" << endl;
    cout << "EP counts: " << epCounts << endl;
    cout << "EPG counts: " << epgCounts << " +/- " << epgCountsError << endl;
    cout << "EPG Per EP: " << epgPerEP << " +/- " << epgPerEPError << endl;
    cout << "EPBG Per Four Body Decay: " << epbgCounts << " +/- " << epbgCountsError << endl;
    cout << "EPBG Per EP: " << epbgPerEP << " +/- " << epbgPerEPError << endl;
}

TCut RDK2ExpAnalyzer::CreateWindowCut(CoDet detType, int peakOption,int det)
{
    double* windows;
    TString detString;
    if(detType==DET_EPG)
    {
        windows=bgoWindows;
        detString="BGO";
    }
    else if(detType==DET_EPBG)
    {
        windows=bapdWindows;
        detString="BAPD";
    }
    else
    {
        cout << "No window cut needed for this type" << endl;
        return TCut();
    }
    TString cutString=detString+"DeltaT"+int2str(det)+" > "+d2str(windows[peakOption*2])+" && "+detString+"DeltaT"+int2str(det) +" < "+d2str(windows[peakOption*2+1]);
    return TCut(cutString.Data());
}

double RDK2ExpAnalyzer::GetPrePeakScaleFactor(CoDet detType)
{
    double* windows;
    if(detType==DET_EPG)
    {
        windows=bgoWindows;
    }
    else if(detType==DET_EPBG)
    {
        windows=bapdWindows;
    }
    else
    {
        cout << "No window cut needed for this type" << endl;
        return 0;
    }
    double scaleFactor=round(windows[3]-windows[2])/round(windows[1]-windows[0]);
    return scaleFactor;
}

void RDK2ExpAnalyzer::CalcEPCountsForLiveDets(CoDet detType)
{
    int* theDetCounts;
    int numDets;
    string detString;
    double* average;
    int* liveDetCount;
    if(detType==DET_EPG)
    {
        theDetCounts=epCountsForLiveBGO;
        numDets=Exp_NUMBGO;
        detString="BGO";
        average=&averageEPCountsForLiveBGO;
        liveDetCount=&numLiveBGO;
    }
    else if(detType==DET_EPBG)
    {
        theDetCounts=epCountsForLiveBAPD;
        numDets=Exp_NUMBAPD;
        detString="BAPD";
        average=&averageEPCountsForLiveBAPD;
        liveDetCount=&numLiveBAPD;
    }
    else
    {
        cout << "Invalid dettype for live dets" << endl;
        return;
    }

    *liveDetCount=0;
    *average=0;

    for (int i = 0;i< numDets;i++)
    {
        HistDim simpleDim={3,-10,800};
        TCut gCut(detString+"EDep"+int2str(i+1)+ ">= 0.0"); //Because dead detectors were set to negative.  If the detectors didn't detect anything they are just 0
        TH1* liveDetHist=MakeHist("liveDetHist","SBDEDepE",gCut&&theCutSet.CreateExpEPCut(),simpleDim);
        theDetCounts[i]=liveDetHist->Integral();
        if(theDetCounts[i]==0) //Dead Detector
        {
            theDetCounts[i]=-1;
        }
        else
        {
            (*liveDetCount)++;
        }
        delete liveDetHist;
    }

    for (int i = 0;i< numDets;i++)
    {
        if(theDetCounts[i] != -1)
        {
            (*average)+=((double)(theDetCounts[i])/(double) numDets);
        }
    }

    return;
}

double RDK2ExpAnalyzer::GetWindowVal(CoDet detType, int i)
{
    if(detType==DET_EPBG)
    {
        return bapdWindows[i];
    }

    return bgoWindows[i];
}

const int NUM_EXP_SERIES_G1=62;
const int EXP_SERIES_ID_LIST_G1[NUM_EXP_SERIES_G1]={72,76,77,78,79,87,89,91,93,94,95,104,105,106,108,119,120,126,127,151,153,154,155,156,160,161,162,165,166,175,176,185,187,213,214,215,223,227,228,253,256,258,260,261,262,263,265,267,270,289,290,291,292,304,305,311,313,314,317,318,319,321};
const int NUM_EXP_SERIES_G2=54;
const int EXP_SERIES_ID_LIST_G2[NUM_EXP_SERIES_G2]={72,76,77,78,79,87,89,91,93,94,95,104,105,106,108,119,120,126,127,151,153,154,155,156,160,161,162,165,166,175,176,185,187,213,214,215,223,227,228,253,270,289,290,291,292,304,305,311,313,314,317,318,319,321}; //without G3
const int NUM_EXP_SERIES_G3=8;
const int EXP_SERIES_ID_LIST_G3[NUM_EXP_SERIES_G3]={256,258,260,261,262,263,265,267};  //Strange time of flight profile
const int NUM_EXP_SERIES_G4=22;
const int EXP_SERIES_ID_LIST_G4[NUM_EXP_SERIES_G4]={256,258,260,261,262,263,265,267,270,289,290,291,292,304,305,311,313,314,317,318,319,321};  //BAPDs
const int NUM_EXP_SERIES_G5=14;
const int EXP_SERIES_ID_LIST_G5[NUM_EXP_SERIES_G5]={270,289,290,291,292,304,305,311,313,314,317,318,319,321};  //BAPDs minus G3
const int NUM_EXP_SERIES_G6=61;
const int EXP_SERIES_ID_LIST_G6[NUM_EXP_SERIES_G6]={76,77,78,79,87,89,91,93,94,95,104,105,106,108,119,120,126,127,151,153,154,155,156,160,161,162,165,166,175,176,185,187,213,214,215,223,227,228,253,256,258,260,261,262,263,265,267,270,289,290,291,292,304,305,311,313,314,317,318,319,321}; //Updated BGO series from Jeff 121022
const int NUM_EXP_SERIES_GTEST=1;
const int EXP_SERIES_ID_LIST_GTEST[NUM_EXP_SERIES_GTEST]={87}; //Updated BGO series from Jeff 121022
const int NUM_EXP_SERIES_G6_Cycle3_4=18;
const int EXP_SERIES_ID_LIST_G6_Cycle3_4[NUM_EXP_SERIES_G6_Cycle3_4]={76,77,78,79,87,89,91,93,94,95,104,105,106,108,119,120,126,127}; //Updated BGO series from Jeff 121022
const int NUM_EXP_SERIES_G6_Cycle5_8=43;
const int EXP_SERIES_ID_LIST_G6_Cycle5_8[NUM_EXP_SERIES_G6_Cycle5_8]={151,153,154,155,156,160,161,162,165,166,175,176,185,187,213,214,215,223,227,228,253,256,258,260,261,262,263,265,267,270,289,290,291,292,304,305,311,313,314,317,318,319,321}; //Updated BGO series from Jeff 121022
const int NUM_EXP_SERIES_G7=46;
const int EXP_SERIES_ID_LIST_G7[NUM_EXP_SERIES_G7]={76,77,78,79,87,89,91,93,94,95,104,105,106,119,120,126,127,151,153,155,156,160,161,162,165,166,175,176,185,187,214,215,223,227,228,253,270,289,290,291,292,304,305,311,313,317}; //Updated BGO series from 2013 collab
const int NUM_EXP_SERIES_G7_Cycle3_4=17;
const int EXP_SERIES_ID_LIST_G7_Cycle3_4[NUM_EXP_SERIES_G7_Cycle3_4]={76,77,78,79,87,89,91,93,94,95,104,105,106,119,120,126,127}; //Updated BGO series from 2013 collab
const int NUM_EXP_SERIES_G7_Cycle5_8=29;
const int EXP_SERIES_ID_LIST_G7_Cycle5_8[NUM_EXP_SERIES_G7_Cycle5_8]={151,153,155,156,160,161,162,165,166,175,176,185,187,214,215,223,227,228,253,270,289,290,291,292,304,305,311,313,317}; //Updated BGO series from 2013 collab

const int NUM_EXP_SERIES_G7_Cycle3=10;
const int EXP_SERIES_ID_LIST_G7_Cycle3[NUM_EXP_SERIES_G7_Cycle3]={76,77,78,79,87,89,91,93,94,95};
const int NUM_EXP_SERIES_G7_Cycle4=7;
const int EXP_SERIES_ID_LIST_G7_Cycle4[NUM_EXP_SERIES_G7_Cycle3]={104,105,106,119,120,126,127};
const int NUM_EXP_SERIES_G7_Cycle5=13;
const int EXP_SERIES_ID_LIST_G7_Cycle5[NUM_EXP_SERIES_G7_Cycle5]={151,153,155,156,160,161,162,165,166,175,176,185,187}; //Updated BGO series from 2013 collab
const int NUM_EXP_SERIES_G7_Cycle6=5;
const int EXP_SERIES_ID_LIST_G7_Cycle6[NUM_EXP_SERIES_G7_Cycle6]={214,215,223,227,228}; //Updated BGO series from 2013 collab
const int NUM_EXP_SERIES_G7_Cycle7=2;
const int EXP_SERIES_ID_LIST_G7_Cycle7[NUM_EXP_SERIES_G7_Cycle7]={253,270}; //Updated BGO series from 2013 collab
const int NUM_EXP_SERIES_G7_Cycle8=9;
const int EXP_SERIES_ID_LIST_G7_Cycle8[NUM_EXP_SERIES_G7_Cycle8]={289,290,291,292,304,305,311,313,317}; //Updated BGO series from 2013 collab

const int NUM_EXP_SERIES_G8=3;
const int EXP_SERIES_ID_LIST_G8[NUM_EXP_SERIES_G8]={233,288,312}; //140305 e trigger only
const int NUM_EXP_SERIES_G9=2;
const int EXP_SERIES_ID_LIST_G9[NUM_EXP_SERIES_G9]={288,312}; //140305 e trigger only removed 233

int RDK2ExpAnalyzer::getSeriesArray(TString inpGroup, int** outArray)
{
    int out;

    if(inpGroup[0]=='S')
    {
        out=1;
        *outArray=new int[1];
        TString sNumString=inpGroup(1,inpGroup.Length()-1);
        int sNum=sNumString.Atoi();
        (*outArray)[0]=sNum;
        return out;
    }

    const int* original;
    if(inpGroup == "G7_Cycle3_4")
    {
        out=NUM_EXP_SERIES_G7_Cycle3_4;
        original=EXP_SERIES_ID_LIST_G7_Cycle3_4;
    }
    else if(inpGroup == "G7_Cycle5_8")
    {
        out=NUM_EXP_SERIES_G7_Cycle5_8;
        original=EXP_SERIES_ID_LIST_G7_Cycle5_8;
    }
    else if(inpGroup == "G7")
    {
        out=NUM_EXP_SERIES_G7;
        original=EXP_SERIES_ID_LIST_G7;
    }
    else
    {
        cout << "Group: " << inpGroup << " not recognized." << endl;
        out=0;
    }
    *outArray=new int[out];
    for (int i = 0;i< out;i++)
    {
    	(*outArray)[i]=original[i];
    }
    return out;
}
