#ifndef RDK2CUTSET_H
#define RDK2CUTSET_H

#include "RDK2Constants.h"
#include "TNamed.h"
#include "TString.h"
#include "TCut.h"


class RDK2CutSet : public TNamed
{
    public:
        RDK2CutSet();
        RDK2CutSet(TString inpName, TString inpTitle,double inpEELow, double inpEEHigh, double inpPELow, double inpPEHigh,
                   double inpPTLow, double inpPTHigh, double inpGELow, double inpGEHigh, double inpBGELow, double inpBGEHigh);
        virtual ~RDK2CutSet();

        ///Public Methods
        const TCut CreateMCECut(int cutOption=0);
        const TCut CreateMCPCut(int cutOption=0);
        const TCut CreateMCGCut(int cutOption=0, int detector=1,CoDet detType=DET_EPG);
        const TCut CreateMCFalseGCut(int cutOption=0, int detector=1,CoDet detType=DET_EPG);
        const TCut CreateMCEPCut(int pCutOption=0,int eCutOption=0);
        const TCut CreateMCEPGCut(int pCutOption=0,int eCutOption=0,int gCutOption=0,int detector=1,CoDet detType=DET_EPG);
        const TCut CreateMCEPFalseGCut(int pCutOption=0,int eCutOption=0,int gCutOption=0,int detector=1,CoDet detType=DET_EPG);
        const TCut CreateExpECut(int cutOption=0);
        const TCut CreateExpPCut(int cutOption=0);
        const TCut CreateExpGCut(int cutOption=0, int detector=1,CoDet detType=DET_EPG);
        const TCut CreateExpEPCut(int pCutOption=0,int eCutOption=0);
        const TCut CreateExpEPGCut(int pCutOption=0,int eCutOption=0,int gCutOption=0,int detector=1,CoDet detType=DET_EPG);
        inline bool CheckEE(double inp){return inp>eELow && inp<eEHigh;}
        inline bool CheckPE(double inp){return inp>pELow && inp<pEHigh;}
        inline bool CheckPT(double inp){return inp>pTLow && inp<pTHigh;}
        inline bool CheckGE(double inp){return inp>gELow && inp<gEHigh;}
        inline bool CheckBGE(double inp){return inp>bGELow && inp<bGEHigh;}
        inline bool CheckEP(double eE,double pE,double pT){return CheckEE(eE) && CheckPE(pE) && CheckPT(pT);}
        inline bool CheckEPG(double eE,double pE,double pT,double gE){return CheckEP(eE,pE,pT) && CheckGE(gE);}
        inline bool CheckEPBG(double eE,double pE,double pT,double bGE){return CheckEP(eE,pE,pT) && CheckBGE(bGE);}
        TString GetDisplayString(CoDet detType);

        ///Setters
        inline void SetEELow(double inp){eELow=inp;}
        inline void SetEEHigh(double inp){eEHigh=inp;}
        inline void SetPELow(double inp){pELow=inp;}
        inline void SetPEHigh(double inp){pEHigh=inp;}
        inline void SetPTLow(double inp){pTLow=inp;}
        inline void SetPTHigh(double inp){pTHigh=inp;}
        inline void SetGELow(double inp){gELow=inp;}
        inline void SetGEHigh(double inp){gEHigh=inp;}
        inline void SetBGELow(double inp){bGELow=inp;}
        inline void SetBGEHigh(double inp){bGEHigh=inp;}


        ///Getters
        inline double GetEELow(){return eELow;}
        inline double GetEEHigh(){return eEHigh;}
        inline double GetPELow(){return pELow;}
        inline double GetPEHigh(){return pEHigh;}
        inline double GetPTLow(){return pTLow;}
        inline double GetPTHigh(){return pTHigh;}
        inline double GetGELow(){return gELow;}
        inline double GetGEHigh(){return gEHigh;}
        inline double GetBGELow(){return bGELow;}
        inline double GetBGEHigh(){return bGEHigh;}


    protected:
        double eELow; //Electron Energy Cut Low
        double eEHigh; //Electron Energy Cut High
        double pELow; //Proton Energy Cut Low
        double pEHigh; //Proton Energy Cut High
        double pTLow; //Proton Time Cut Low
        double pTHigh; //Proton Time Cut High
        double gELow; //Gamma Energy Cut Low
        double gEHigh; //Gamma Energy Cut High
        double bGELow; //bAPD Gamma Energy Cut Low
        double bGEHigh; //bAPD Gamma Energy Cut Low

    private:

    public:
        ClassDef(RDK2CutSet,1); //Class to analyze RDK
};

#endif // RDK2CUTSET_H
