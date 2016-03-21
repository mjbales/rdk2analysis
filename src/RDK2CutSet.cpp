#include "RDK2CutSet.h"

#include "MRText.h"

#include <iostream>

using namespace std;

ClassImp(RDK2CutSet);

RDK2CutSet::RDK2CutSet()
{
    SetName("RDK2CutSet");

    eELow=0.; //Electron Energy Cut Low
    eEHigh=0.; //Electron Energy Cut High
    pELow=0.; //Proton Energy Cut Low
    pEHigh=0.; //Proton Energy Cut High
    pTLow=0.; //Proton Time Cut Low
    pTHigh=0.; //Proton Time Cut High
    gELow=0.; //Gamma Energy Cut Low
    gEHigh=0.; //Gamma Energy Cut High
    bGELow=0.; //bapd Gamma Energy Cut Low
    bGEHigh=0.; //bapd Gamma Energy Cut High
}

RDK2CutSet::RDK2CutSet(TString inpName, TString inpTitle,double inpEELow, double inpEEHigh, double inpPELow, double inpPEHigh, double inpPTLow, double inpPTHigh, double inpGELow, double inpGEHigh, double inpBGELow, double inpBGEHigh)
{
    SetName(inpName);
    SetTitle(inpTitle);
    eELow=inpEELow; //Electron Energy Cut Low
    eEHigh=inpEEHigh; //Electron Energy Cut High
    pELow=inpPELow; //Proton Energy Cut Low
    pEHigh=inpPEHigh; //Proton Energy Cut High
    pTLow=inpPTLow; //Proton Time Cut Low
    pTHigh=inpPTHigh; //Proton Time Cut High
    gELow=inpGELow; //Gamma Energy Cut Low
    gEHigh=inpGEHigh; //Gamma Energy Cut High
    bGELow=inpBGELow; //bapd Gamma Energy Cut Low
    bGEHigh=inpBGEHigh; //bapd Gamma Energy Cut High
}

RDK2CutSet::~RDK2CutSet()
{
}





const TCut RDK2CutSet::CreateMCECut(int cutOption)
{
    TString cutString;
    if(cutOption==0)
    {
        cutString="eDChn.SBDEDepBlur > " + d2str(eELow)+ " && eDChn.SBDEDepBlur < " + d2str(eEHigh);
    }
    else
    {
        cout << "Error MCECut Option not valid!" << endl;
    }
    TCut outCut(cutString.Data());
    return outCut;
}

const TCut RDK2CutSet::CreateMCPCut(int cutOption)
{
    TString cutString;
    if(cutOption==0) //No resolution on proton energy
    {
        cutString+="pChn.SBDEDepTotal > " + d2str(pELow)+ " && pChn.SBDTimeFirst-eChn.SBDTimeFirst < " + d2str(pEHigh);
        cutString+=" && pChn.SBDTimeFirst-eChn.SBDTimeFirst > " + d2str(pTLow)+ " && pChn.SBDTimeFirst-eChn.SBDTimeFirstt < " + d2str(pTHigh);
    }
    else if(cutOption==1) //Use resolution on proton energy
    {
        cutString+="pDChn.SBDEDepBlur > " + d2str(pELow)+ " && pDChn.SBDEDepBlur < " + d2str(pEHigh);
        cutString+=" && pChn.SBDTimeFirst-eChn.SBDTimeFirst > " + d2str(pTLow)+ " && pChn.SBDTimeFirst-eChn.SBDTimeFirst < " + d2str(pTHigh);
    }
    else
    {
        cout << "Error MCPCut Option not valid!" << endl;
    }
    TCut outCut(cutString.Data());
    return outCut;
}

const TCut RDK2CutSet::CreateMCGCut(int cutOption,int detector,CoDet detType)
{
    TString cutString;
    if(detType==DET_EPG || detType==DET_G)
    {
        if(cutOption==0)
        {
            cutString+="gDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParam > " + d2str(gELow)+ " && gDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParam < " + d2str(gEHigh);
        }
        else if(cutOption==2)
        {
        	cutString+="gDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParamHigh > " + d2str(gELow)+ " && gDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParamHigh < " + d2str(gEHigh);
        }
        else if(cutOption==3)
		{
			cutString+="gDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParamLow > " + d2str(gELow)+ " && gDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParamLow < " + d2str(gEHigh);
		}
        else
        {
            cout << "Error MCGCut Option not valid!" << endl;
        }
    }
    else if(detType==DET_EPBG || detType==DET_BG)
    {
        if(cutOption==0)
        {
            cutString+="gChn.BAPD"+int2str(detector)+"EDepTotal > " + d2str(bGELow)+ " && gChn.BAPD"+int2str(detector)+"EDepTotal < " + d2str(bGEHigh);
        }
        else
        {
            cout << "Error MCBGCut Option not valid!" << endl;
        }
    }
    TCut outCut(cutString.Data());
    return outCut;
}

const TCut RDK2CutSet::CreateMCFalseGCut(int cutOption,int detector,CoDet detType)
{
    TString cutString;
    if(detType==DET_EPG)
    {
        if(cutOption==0)
        {
            cutString+="eDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParam > " + d2str(gELow)+ " && eDChn.BGO"+int2str(detector)+"EDepBlurLOAllExpParam < " + d2str(gEHigh);
        }
        else
        {
            cout << "Error MCFalseGCut Option not valid!" << endl;
        }
    }
    else
    {
        if(cutOption==0)
        {
            cutString+="eChn.BAPD"+int2str(detector)+"EDepTotal > " + d2str(bGELow)+ " && eChn.BAPD"+int2str(detector)+"EDepTotal < " + d2str(bGEHigh);
        }
        else
        {
            cout << "Error MCFalseBGCut Option not valid!" << endl;
        }
    }
    TCut outCut(cutString.Data());
    return outCut;
}

const TCut RDK2CutSet::CreateMCEPCut(int pCutOption,int eCutOption)
{

    TCut theEPCut=CreateMCECut(eCutOption) && CreateMCPCut(pCutOption);
    return theEPCut;
}


const TCut RDK2CutSet::CreateMCEPGCut(int pCutOption,int eCutOption,int gCutOption,int detector,CoDet detType)
{
    TCut theEPGCut=CreateMCEPCut(pCutOption,eCutOption) && CreateMCGCut(gCutOption,detector,detType);
    return theEPGCut;
}

const TCut RDK2CutSet::CreateMCEPFalseGCut(int pCutOption,int eCutOption,int gCutOption,int detector,CoDet detType)
{
    TCut theEPGCut=CreateMCEPCut(pCutOption,eCutOption) && CreateMCFalseGCut(gCutOption,detector,detType);
    return theEPGCut;
}

const TCut RDK2CutSet::CreateExpECut(int cutOption)
{
    TString cutString;
    if(cutOption==0)
    {
        cutString="SBDEDepE > " + d2str(eELow)+ " && SBDEDepE < " + d2str(eEHigh);
    }
    else
    {
        cout << "Error ExpECut Option not valid!" << endl;
    }
    TCut outCut(cutString.Data());
    return outCut;
}

const TCut RDK2CutSet::CreateExpPCut(int cutOption)
{
    TString cutString;
    if(cutOption==0)
    {
        cutString+="SBDEDepP > " + d2str(pELow)+ " && SBDEDepP < " + d2str(pEHigh);
        cutString+=" && SBDToFP >= " + d2str(pTLow)+ " && SBDToFP <= " + d2str(pTHigh);
    }
    else
    {
        cout << "Error ExpPCut Option not valid!" << endl;
    }
    TCut outCut(cutString.Data());
    return outCut;
}

const TCut RDK2CutSet::CreateExpGCut(int cutOption,int detector,CoDet detType)
{
    TString cutString;
    if(detType==DET_EPG)
    {
        if(cutOption==0)
        {
            cutString+="BGOEDep"+int2str(detector)+" > " + d2str(gELow)+ " && BGOEDep"+int2str(detector)+" < " + d2str(gEHigh);
        }
        else
        {
            cout << "Error ExpGCut Option not valid!" << endl;
        }
    }
    else
    {
        if(cutOption==0)
        {
            cutString+="BAPDEDep"+int2str(detector)+" > " + d2str(bGELow)+ " && BAPDEDep"+int2str(detector)+" < " + d2str(bGEHigh);
        }
        else
        {
            cout << "Error ExpBGCut Option not valid!" << endl;
        }
    }
    TCut outCut(cutString.Data());
    return outCut;
}

const TCut RDK2CutSet::CreateExpEPCut(int pCutOption,int eCutOption)
{

    TCut theEPCut=CreateExpECut(eCutOption) && CreateExpPCut(pCutOption);
    return theEPCut;
}


const TCut RDK2CutSet::CreateExpEPGCut(int pCutOption,int eCutOption,int gCutOption,int detector,CoDet detType)
{
    TCut theEPGCut=CreateExpEPCut(pCutOption,eCutOption) && CreateExpGCut(gCutOption,detector,detType);
    return theEPGCut;
}

TString RDK2CutSet::GetDisplayString(CoDet detType)
{
    TString out;

    out=d2str(pTLow*1e6)+" #mus < pT < "+d2str(pTHigh*1e6)+ " #mus";
    out+="    "+d2str(eELow)+" keV < eE < "+d2str(eEHigh)+ " keV";
    out+="    "+d2str(pELow)+" keV < pE < "+d2str(pEHigh)+ " keV";
    if(detType == DET_EPG)
    {
        out+="    "+d2str(gELow)+" keV < gE < "+d2str(gEHigh)+ " keV";
    }
    if(detType == DET_EPBG)
    {
        out+="    "+d2str(bGELow)+" keV < bGE < "+d2str(bGEHigh)+ " keV";
    }
    return out;
}



