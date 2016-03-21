#ifndef RDK2VARS_H_INCLUDED
#define RDK2VARS_H_INCLUDED

#include <algorithm>

#include "TROOT.h"
#include "TCut.h"

#include "RDK2Constants.h"
#include "RDK2CutSet.h"

//const TCut STD_Mult1Cut("((gDChn.BGO1EDepBlur > 10.)*1 + (gDChn.BGO2EDepBlur > 10.)*1 + (gDChn.BGO3EDepBlur > 10.)*1 + (gDChn.BGO4EDepBlur > 10.)*1 + (gDChn.BGO5EDepBlur > 10.)*1 + (gDChn.BGO6EDepBlur > 10.)*1 + (gDChn.BGO7EDepBlur > 10.)*1 + (gDChn.BGO8EDepBlur > 10.)*1 + (gDChn.BGO9EDepBlur > 10.)*1 + (gDChn.BGO10EDepBlur > 10.)*1 + (gDChn.BGO11EDepBlur > 10.)*1 + (gDChn.BGO12EDepBlur > 10.)*1 ) ==1");
const TCut STD_Mult1Cut("BGOmultiplicityEdge==1");
const TCut STD_PCut("pChn.SBDEDepTotal > 10. && pChn.SBDEDepTotal < 31. && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 10.e-6");  //modified 130617
const TCut STD_OLDPCut("pChn.SBDEDepTotal > 7. && pChn.SBDEDepTotal < 31. && pChn.SBDTimeFirst > 2.0e-6 && pChn.SBDTimeFirst < 25.e-6");
const TCut STD_ECut("eDChn.SBDEDepBlur > 50. && eDChn.SBDEDepBlur < 800.");
const TCut STD_ECutB("eDChn.SBDEDepBlur > 100. && eDChn.SBDEDepBlur < 800."); //Higher electron energy cut
const TCut STD_ECutNoBlur("eChn.SBDEDepTotal > 50. && eChn.SBDEDepTotal < 800.");
const TCut STD_EPCut=STD_PCut && STD_ECut;
const TCut STD_EPCutB=STD_PCut && STD_ECutB; //Higher electron energy cut
const TCut STD_OLDEPCut=STD_OLDPCut && STD_ECut;
const TCut STD_EPGMult1Cut=STD_EPCut && STD_Mult1Cut;
const TCut STD_bZeroCut("bChn.dwcutb20 == 0");
const TCut STD_EPHmgCut=STD_EPCut && STD_bZeroCut;
const TCut STD_BGOTotalCut="(gDChn.BGO1EDepBlur + gDChn.BGO2EDepBlur + gDChn.BGO3EDepBlur + gDChn.BGO4EDepBlur + gDChn.BGO5EDepBlur + gDChn.BGO6EDepBlur + gDChn.BGO7EDepBlur + gDChn.BGO8EDepBlur + gDChn.BGO9EDepBlur + gDChn.BGO10EDepBlur + gDChn.BGO11EDepBlur + gDChn.BGO12EDepBlur) > 10";
const TCut STD_EPGCut=STD_EPCut && STD_BGOTotalCut;
const TCut STD_BGOTotalHighEnergyCut="(gDChn.BGO1EDepBlur + gDChn.BGO2EDepBlur + gDChn.BGO3EDepBlur + gDChn.BGO4EDepBlur + gDChn.BGO5EDepBlur + gDChn.BGO6EDepBlur + gDChn.BGO7EDepBlur + gDChn.BGO8EDepBlur + gDChn.BGO9EDepBlur + gDChn.BGO10EDepBlur + gDChn.BGO11EDepBlur + gDChn.BGO12EDepBlur) > 200";
const TCut STD_BGOTotalLowEnergyCut="(gDChn.BGO1EDepBlur + gDChn.BGO2EDepBlur + gDChn.BGO3EDepBlur + gDChn.BGO4EDepBlur + gDChn.BGO5EDepBlur + gDChn.BGO6EDepBlur + gDChn.BGO7EDepBlur + gDChn.BGO8EDepBlur + gDChn.BGO9EDepBlur + gDChn.BGO10EDepBlur + gDChn.BGO11EDepBlur + gDChn.BGO12EDepBlur) < 200";

const TCut STD_Mult1CutUnblur("((gChn.BGO1EDepTotal > 10.)*1 + (gChn.BGO2EDepTotal > 10.)*1 + (gChn.BGO3EDepTotal > 10.)*1 + (gChn.BGO4EDepTotal > 10.)*1 + (gChn.BGO5EDepTotal > 10.)*1 + (gChn.BGO6EDepTotal > 10.)*1 + (gChn.BGO7EDepTotal > 10.)*1 + (gChn.BGO8EDepTotal > 10.)*1 + (gChn.BGO9EDepTotal > 10.)*1 + (gChn.BGO10EDepTotal > 10.)*1 + (gChn.BGO11EDepTotal > 10.)*1 + (gChn.BGO12EDepTotal > 10.)*1 ) ==1");
const TCut STD_ECutUnblur("eChn.SBDEDepTotal > 50. && eChn.SBDEDepTotal < 800.");
const TCut STD_EPCutUnblur=STD_PCut && STD_ECutUnblur;
const TCut STD_EPGCutUnblur=STD_EPCutUnblur && STD_Mult1CutUnblur;
const TCut STD_BAPDTotalCut="(gChn.BAPD1EDepTotal + gChn.BAPD2EDepTotal + gChn.BAPD3EDepTotal) >.3";
const TCut STD_BAPD12TotalCut="(gChn.BAPD1EDepTotal + gChn.BAPD2EDepTotal) >.3";


const TString R_BGOMAX="max(gDChn.BGO1EDepBlur,max(gDChn.BGO2EDepBlur,max(gDChn.BGO3EDepBlur,max(gDChn.BGO4EDepBlur,max(gDChn.BGO5EDepBlur,max(gDChn.BGO6EDepBlur,max(gDChn.BGO7EDepBlur,max(gDChn.BGO8EDepBlur,max(gDChn.BGO9EDepBlur,max(gDChn.BGO10EDepBlur,max(gDChn.BGO11EDepBlur,gDChn.BGO12EDepBlur)))))))))))";
const TString R_BGOMAXUNBLUR="max(gChn.BGO1EDepTotal,max(gChn.BGO2EDepTotal,max(gChn.BGO3EDepTotal,max(gChn.BGO4EDepTotal,max(gChn.BGO5EDepTotal,max(gChn.BGO6EDepTotal,max(gChn.BGO7EDepTotal,max(gChn.BGO8EDepTotal,max(gChn.BGO9EDepTotal,max(gChn.BGO10EDepTotal,max(gChn.BGO11EDepTotal,gChn.BGO12EDepTotal)))))))))))";
const TString R_BGOTOTALEP="pDChn.BGO1EDepBlur  + pDChn.BGO2EDepBlur  + pDChn.BGO3EDepBlur  + pDChn.BGO4EDepBlur  + pDChn.BGO5EDepBlur  + pDChn.BGO6EDepBlur  + pDChn.BGO7EDepBlur  + pDChn.BGO8EDepBlur  + pDChn.BGO9EDepBlur  + pDChn.BGO10EDepBlur  + pDChn.BGO11EDepBlur  + pDChn.BGO12EDepBlur + eDChn.BGO1EDepBlur  + eDChn.BGO2EDepBlur  + eDChn.BGO3EDepBlur  + eDChn.BGO4EDepBlur  + eDChn.BGO5EDepBlur  + eDChn.BGO6EDepBlur  + eDChn.BGO7EDepBlur  + eDChn.BGO8EDepBlur  + eDChn.BGO9EDepBlur  + eDChn.BGO10EDepBlur  + eDChn.BGO11EDepBlur  + eDChn.BGO12EDepBlur";
const TString R_BGOTOTAL="(gDChn.BGO1EDepBlur + gDChn.BGO2EDepBlur + gDChn.BGO3EDepBlur + gDChn.BGO4EDepBlur + gDChn.BGO5EDepBlur + gDChn.BGO6EDepBlur + gDChn.BGO7EDepBlur + gDChn.BGO8EDepBlur + gDChn.BGO9EDepBlur + gDChn.BGO10EDepBlur + gDChn.BGO11EDepBlur + gDChn.BGO12EDepBlur)";
const TString R_BGOTOTALEDGE="(gDChn.BGO1EDepBlurEdge + gDChn.BGO2EDepBlurEdge + gDChn.BGO3EDepBlurEdge + gDChn.BGO4EDepBlurEdge + gDChn.BGO5EDepBlurEdge + gDChn.BGO6EDepBlurEdge + gDChn.BGO7EDepBlurEdge + gDChn.BGO8EDepBlurEdge + gDChn.BGO9EDepBlurEdge + gDChn.BGO10EDepBlurEdge + gDChn.BGO11EDepBlurEdge + gDChn.BGO12EDepBlurEdge)";
const TString R_BGOTOTALEDGEDBLGAUSS="gDChn.BGO1EDepBlurEdgeDBLGAUSS + gDChn.BGO2EDepBlurEdgeDBLGAUSS + gDChn.BGO3EDepBlurEdgeDBLGAUSS + gDChn.BGO4EDepBlurEdgeDBLGAUSS + gDChn.BGO5EDepBlurEdgeDBLGAUSS + gDChn.BGO6EDepBlurEdgeDBLGAUSS + gDChn.BGO7EDepBlurEdgeDBLGAUSS + gDChn.BGO8EDepBlurEdgeDBLGAUSS + gDChn.BGO9EDepBlurEdgeDBLGAUSS + gDChn.BGO10EDepBlurEdgeDBLGAUSS + gDChn.BGO11EDepBlurEdgeDBLGAUSS + gDChn.BGO12EDepBlurEdgeDBLGAUSS";
const TString R_BGOTOTALEDGEDBLRES="gDChn.BGO1EDepBlurEdgeDBLRES + gDChn.BGO2EDepBlurEdgeDBLRES + gDChn.BGO3EDepBlurEdgeDBLRES + gDChn.BGO4EDepBlurEdgeDBLRES + gDChn.BGO5EDepBlurEdgeDBLRES + gDChn.BGO6EDepBlurEdgeDBLRES + gDChn.BGO7EDepBlurEdgeDBLRES + gDChn.BGO8EDepBlurEdgeDBLRES + gDChn.BGO9EDepBlurEdgeDBLRES + gDChn.BGO10EDepBlurEdgeDBLRES + gDChn.BGO11EDepBlurEdgeDBLRES + gDChn.BGO12EDepBlurEdgeDBLRES";
const TString R_BGOTOTALEDGENORES="gDChn.BGO1EDepBlurEdgeNORES + gDChn.BGO2EDepBlurEdgeNORES + gDChn.BGO3EDepBlurEdgeNORES + gDChn.BGO4EDepBlurEdgeNORES + gDChn.BGO5EDepBlurEdgeNORES + gDChn.BGO6EDepBlurEdgeNORES + gDChn.BGO7EDepBlurEdgeNORES + gDChn.BGO8EDepBlurEdgeNORES + gDChn.BGO9EDepBlurEdgeNORES + gDChn.BGO10EDepBlurEdgeNORES + gDChn.BGO11EDepBlurEdgeNORES + gDChn.BGO12EDepBlurEdgeNORES";
const TString R_AZANGLEPROTON  = "(atan2(sqrt( mxp0 * mxp0 + myp0 * myp0), mzp0 )*57.29577951-90.)";
const TString R_AZANGLEELECTRON  = "(atan2(sqrt( mxe0 * mxe0 + mye0 * mye0), mze0 )*57.29577951-90.)";
const TString R_AZANGLEGAMMA  = "(atan2(sqrt( mxg0 * mxg0 + myg0 * myg0), mzg0 )*57.29577951-90.)";
const TString R_AZANGLEGAMMA_DIVSINTHETA  = "(atan2(sqrt( mxg0 * mxg0 + myg0 * myg0), mzg0 )*57.29577951-90.)/sqrt(1-mzg0*mzg0)";
const TString R_BGOONPEAK="(BGOEDep1 >10 && BGODeltaT1 > -75.01 && BGODeltaT1 < 74.99)*BGOEDep1 + (BGOEDep2 >10 && BGODeltaT2 > -75.01 && BGODeltaT2 < 74.99)*BGOEDep2 + (BGOEDep3 >10 && BGODeltaT3 > -75.01 && BGODeltaT3 < 74.99)*BGOEDep3 + (BGOEDep4 >10 && BGODeltaT4 > -75.01 && BGODeltaT4 < 74.99)*BGOEDep4 + (BGOEDep5 >10 && BGODeltaT5 > -75.01 && BGODeltaT5 < 74.99)*BGOEDep5 + (BGOEDep6 >10 && BGODeltaT6 > -75.01 && BGODeltaT6 < 74.99)*BGOEDep6 + (BGOEDep7 >10 && BGODeltaT7 > -75.01 && BGODeltaT7 < 74.99)*BGOEDep7 + (BGOEDep8 >10 && BGODeltaT8 > -75.01 && BGODeltaT8 < 74.99)*BGOEDep8 + (BGOEDep9 >10 && BGODeltaT9 > -75.01 && BGODeltaT9 < 74.99)*BGOEDep9 + (BGOEDep10 >10 && BGODeltaT10 > -75.01 && BGODeltaT10 < 74.99)*BGOEDep10 + (BGOEDep11 >10 && BGODeltaT11 > -75.01 && BGODeltaT11 < 74.99)*BGOEDep11 + (BGOEDep12 >10 && BGODeltaT12 > -75.01 && BGODeltaT12 < 74.99)*BGOEDep12";
const TString R_BGOPREPEAK="(BGOEDep1 >10 && BGODeltaT1 > -450.01 && BGODeltaT1 < -300.99)*BGOEDep1 + (BGOEDep2 >10 && BGODeltaT2 > -450.01 && BGODeltaT2 < -300.99)*BGOEDep2 + (BGOEDep3 >10 && BGODeltaT3 > -450.01 && BGODeltaT3 < -300.99)*BGOEDep3 + (BGOEDep4 >10 && BGODeltaT4 > -450.01 && BGODeltaT4 < -300.99)*BGOEDep4 + (BGOEDep5 >10 && BGODeltaT5 > -450.01 && BGODeltaT5 < -300.99)*BGOEDep5 + (BGOEDep6 >10 && BGODeltaT6 > -450.01 && BGODeltaT6 < -300.99)*BGOEDep6 + (BGOEDep7 >10 && BGODeltaT7 > -450.01 && BGODeltaT7 < -300.99)*BGOEDep7 + (BGOEDep8 >10 && BGODeltaT8 > -450.01 && BGODeltaT8 < -300.99)*BGOEDep8 + (BGOEDep9 >10 && BGODeltaT9 > -450.01 && BGODeltaT9 < -300.99)*BGOEDep9 + (BGOEDep10 >10 && BGODeltaT10 > -450.01 && BGODeltaT10 < -300.99)*BGOEDep10 + (BGOEDep11 >10 && BGODeltaT11 > -450.01 && BGODeltaT11 < -300.99)*BGOEDep11 + (BGOEDep12 >10 && BGODeltaT12 > -450.01 && BGODeltaT12 < -300.99)*BGOEDep12";
const TString R_BAPDONPEAK="(BAPDEDep1 >.3 && BAPDDeltaT1 > -49.01 && BAPDDeltaT1 < -8.99)*BAPDEDep1 + (BAPDEDep2 >.3 && BAPDDeltaT2 > -49.01 && BAPDDeltaT2 < -8.99)*BAPDEDep2 +(BAPDEDep3 >.3 && BAPDDeltaT3 > -49.01 && BAPDDeltaT3 < -8.99)*BAPDEDep3";
const TString R_BAPDPREPEAK="(BAPDEDep1 >.3 && BAPDDeltaT1 > -449.01 && BAPDDeltaT1 < -408.99)*BAPDEDep1 + (BAPDEDep2 >.3 && BAPDDeltaT2 > -449.01 && BAPDDeltaT2 < -408.99)*BAPDEDep2 +(BAPDEDep3 >.3 && BAPDDeltaT3 > -449.01 && BAPDDeltaT3 < -408.99)*BAPDEDep3";
const TString R_BAPDTOTAL="(gChn.BAPD1EDepTotal >.3)*gChn.BAPD1EDepTotal + (gChn.BAPD2EDepTotal >.3)*gChn.BAPD2EDepTotal + (gChn.BAPD3EDepTotal >.3)*gChn.BAPD3EDepTotal";
const TString R_BAPD12TOTAL="(gChn.BAPD1EDepTotal >.3)*gChn.BAPD1EDepTotal + (gChn.BAPD2EDepTotal >.3)*gChn.BAPD2EDepTotal";

const TString R_BGOSUM="(gDChn.BGO1EDepBlurLOAllExpParam > 10)*gDChn.BGO1EDepBlurLOAllExpParam+(gDChn.BGO2EDepBlurLOAllExpParam > 10)*gDChn.BGO2EDepBlurLOAllExpParam+(gDChn.BGO3EDepBlurLOAllExpParam > 10)*gDChn.BGO3EDepBlurLOAllExpParam+(gDChn.BGO4EDepBlurLOAllExpParam > 10)*gDChn.BGO4EDepBlurLOAllExpParam+(gDChn.BGO5EDepBlurLOAllExpParam > 10)*gDChn.BGO5EDepBlurLOAllExpParam+(gDChn.BGO6EDepBlurLOAllExpParam > 10)*gDChn.BGO6EDepBlurLOAllExpParam+(gDChn.BGO7EDepBlurLOAllExpParam > 10)*gDChn.BGO7EDepBlurLOAllExpParam+(gDChn.BGO8EDepBlurLOAllExpParam > 10)*gDChn.BGO8EDepBlurLOAllExpParam+(gDChn.BGO9EDepBlurLOAllExpParam > 10)*gDChn.BGO9EDepBlurLOAllExpParam+(gDChn.BGO10EDepBlurLOAllExpParam > 10)*gDChn.BGO10EDepBlurLOAllExpParam+(gDChn.BGO11EDepBlurLOAllExpParam > 10)*gDChn.BGO11EDepBlurLOAllExpParam+(gDChn.BGO12EDepBlurLOAllExpParam > 10)*gDChn.BGO12EDepBlurLOAllExpParam";

#endif // RDK2VARS_H_INCLUDED
