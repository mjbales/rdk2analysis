#ifndef RDK2EXPANALYSIS_H
#define RDK2EXPANALYSIS_H

#include "TNamed.h"
#include "TString.h"
#include "TH1.h"


#include "RDK2Constants.h"
#include "RDK2Set.h"
#include "RDK2Vars.h"
#include "RDK2CutSet.h"

#define Exp_NUMBGO 12
#define Exp_NUMBAPD 3

//#define Exp_NUMBGO 1
//#define Exp_NUMBAPD 1

class RDK2ExpAnalysis : public TNamed
{
    public:

        ///Constructers/Destructers
        RDK2ExpAnalysis();
        RDK2ExpAnalysis(TString inpName, TString inpTitle, TString inpPIDString,TString inpGroupString, RDK2CutSet inpCutSet,const double inpBGOWindows[6],const double inpBAPDWindows[6]);
        virtual ~RDK2ExpAnalysis();

        ///Setters
        inline void SetCutSet(RDK2CutSet inp){theCutSet=inp;}
        inline void SetEPCounts(double inp){epCounts=inp;}
        inline void SetEPGCounts(double inp,double error){epgCounts=inp;epgCountsError=error;}
        inline void SetEPBGCounts(double inp,double error){epbgCounts=inp;epbgCountsError=error;}
        inline void SetEPGPerEPCounts(double inp,double error){epgPerEP=inp;epgPerEPError=error;}


        ///Getters
        inline RDK2CutSet GetCutSet(){return theCutSet;}
        inline TH1* GetEP_eEHist(){return ep_eEHist;}
        inline TH1* GetEP_pEHist(){return ep_pEHist;}
        inline TH1* GetEP_pTHist(){return ep_pTHist;}
        inline TH1* GetEPG_eEHist(){return epg_eEHist;}
        inline TH1* GetEPG_pEHist(){return epg_pEHist;}
        inline TH1* GetEPG_pTHist(){return epg_pTHist;}
        inline TH1* GetEPBG_eEHist(){return epbg_eEHist;}
        inline TH1* GetEPBG_pEHist(){return epbg_pEHist;}
        inline TH1* GetEPBG_pTHist(){return epbg_pTHist;}
        inline TH1* GetEPG_gEAvgHist(){return epg_gEAvgHist;}
        inline TH1* GetEPG_gEAvgVarHist(){return epg_gEAvgVarHist;}
        inline TH1* GetEPG_gEDetHists(int inp){return epg_gEDetHists[inp-1];}
        inline TH1* GetEPG_gEDetVarHists(int inp){return epg_gEDetVarHists[inp-1];}
        inline TH1* GetEPBG_bGEAvgHist(){return epbg_bGEAvgHist;}
        inline TH1* GetEPBG_bGEAvgVarHist(){return epbg_bGEAvgVarHist;}
        inline TH1* GetEPBG_bGEDetHists(int inp){return epbg_bGEDetHists[inp-1];}
        inline TH1* GetEPBG_bGEDetVarHists(int inp){return epbg_bGEDetVarHists[inp-1];}
        inline TString GetGroupString(){return groupString;}
        inline TString GetPIDString(){return PIDString;}
        inline int GetEPCounts(){return epCounts;}
        inline double GetEPGCounts(double& error){error=epgCountsError; return epgCounts;}
        inline double GetEPBGCounts(double& error){error=epbgCountsError; return epbgCounts;}
        inline double GetEPGPerEP(double& error){error=epgPerEPError; return epgPerEP;}
        inline double GetEPBGPerEP(double& error){error=epbgPerEPError; return epbgPerEP;}

        double GetWindowVal(CoDet detType, int i);


        TH1* MakeHist(TString histName,TString drawString,TCut inpCut,TH1* histPrototype);
        TH1* MakeHist(TString histName,TString drawString,TCut inpCut,HistDim inpHistDimX, HistDim inpHistDimY=BLANK_DIM,HistDim inpHistDimZ=BLANK_DIM);

        TH1* MakeRadHist(TString histName,CoDet detType,TString drawString, int peakOption,int det,TH1* histPrototype,bool scaleByEPCounts=true);  //Make Rad Hist using cutset
        TH1* MakeRadHist(TString histName,CoDet detType,TString drawString, int peakOption,int det,HistDim inpHistDimX, HistDim inpHistDimY=BLANK_DIM,HistDim inpHistDimZ=BLANK_DIM,bool scaleByEPCounts=true); //Make Rad Hist using cutset


        TH1* MakeRadAvgHist(TString histName,CoDet detType,TString drawString,int peakOption,TH1* histPrototype,bool scaleByEPCounts=true);  //Make RadAverage Hist using cutset
        TH1* MakeRadAvgHist(TString histName,CoDet detType,TString drawString,int peakOption,HistDim inpHistDimX, HistDim inpHistDimY=BLANK_DIM,HistDim inpHistDimZ=BLANK_DIM,bool scaleByEPCounts=true); //Make Rad Average Hist using cutset

        TH1* MakeBGSubHist(TString histName,CoDet detType,TString drawString, int det,TH1* histPrototype,bool scaleByEPCounts=true);  //Make BGSub Hist using cutset
        TH1* MakeBGSubHist(TString histName,CoDet detType,TString drawString, int det,HistDim inpHistDimX, HistDim inpHistDimY=BLANK_DIM,HistDim inpHistDimZ=BLANK_DIM,bool scaleByEPCounts=true); //Make BGSub Hist using cutset

        TH1* MakeBGSubAvgHist(TString histName,CoDet detType,TString drawString,TH1* histPrototype,bool scaleByEPCounts=true);  //Make BGSub Hist using cutset
        TH1* MakeBGSubAvgHist(TString histName,CoDet detType,TString drawString, HistDim inpHistDimX, HistDim inpHistDimY=BLANK_DIM,HistDim inpHistDimZ=BLANK_DIM,bool scaleByEPCounts=true); //Make BGSub Hist using cutset


        TH1* MakeStandardHist(TString histName, CoDet detType,TString drawString,HistDim inpHistDimX);

        void MakeAllStandardHists();

        int CalcEPCounts(); //Calculates epPerDecay using theCutSet
        double CalcEPGCounts(CoDet detType,double& error); //Calculates epgPerDecay using theCutSet
        double CalcEPGPerEP(CoDet detType, double& error); //Calculates epgPerDecay using theCutSet
        void CalcRates(); //Calculates EPG and EP rates
        void CalcEPCountsForLiveDets(CoDet detType);
        void PrintRates(); //Prints rates

    protected:
        ///Variables
        RDK2Set threeBodySet; //Three body Exp RDK2Set
        RDK2Set fourBodySet;  //Four body Exp RDK2Set

        RDK2CutSet theCutSet;  //E, P, and G cut set

        int epCounts;   //ep counts
        double epgCounts;  //epg counts
        double epgCountsError;  //epg counts error
        double epbgCounts; //epbg counts
        double epbgCountsError; //epbg counts error

        double epgPerEP; //EPG/EP detection ratio
        double epgPerEPError;//EPG/EP detection ratio error

        double epbgPerEP; //EPG/EP detection ratio
        double epbgPerEPError;//EPG/EP detection ratio error

        int epCountsForLiveBGO[Exp_NUMBGO]; //ep counts while detector was live Negative = dead
        int epCountsForLiveBAPD[Exp_NUMBAPD]; //ep counts while detector was live Negative = dead

        double averageEPCountsForLiveBGO;
        double averageEPCountsForLiveBAPD;

        int numLiveBGO;  //num live dets. Only has to be live for some of the series
        int numLiveBAPD; //num live des. Only has to be live for some of the series



        ///Hists
        TH1* ep_eEHist;
        TH1* ep_pEHist;
        TH1* ep_pTHist;
        TH1* epg_eEHist;
        TH1* epg_pEHist;
        TH1* epg_pTHist;
        TH1* epbg_eEHist;
        TH1* epbg_pEHist;
        TH1* epbg_pTHist;
        TH1* epg_gEAvgHist;
        TH1* epg_gEAvgVarHist;
        TH1* epg_gEDetHists[Exp_NUMBGO];
        TH1* epg_gEDetVarHists[Exp_NUMBGO];
        TH1* epbg_bGEAvgHist;
        TH1* epbg_bGEAvgVarHist;
        TH1* epbg_bGEDetHists[Exp_NUMBAPD];
        TH1* epbg_bGEDetVarHists[Exp_NUMBAPD];

        TString PIDString;
        TString groupString;

        double bgoWindows[6];
        double bapdWindows[6];


        bool chainLoaded;

        TChain* expChain;


        ///Methods
        void ResetHists();
        void NULLHists();
        void loadEventChain();

        TCut CreateWindowCut(CoDet detType, int peakOption,int det);

        double GetPrePeakScaleFactor(CoDet detType);

        int getSeriesArray(TString inpGroup, int** outArray);



    private:

    public:
        ClassDef(RDK2ExpAnalysis,1); //Class to analyze RDK experimental data
};

#endif // RDK2EXPANALYSIS_H
