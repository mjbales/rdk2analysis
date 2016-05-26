#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED
#pragma once

//Physical or other constants (2006 CODATA)
#define PCHARGE 1.602176487e-19   	//Elementary charge in Columbs (Also joules/eV)
#define ECHARGE -1.602176487e-19
#define PMASS 1.672621637e-27 		//Mass of Proton in kg
#define PMASSE 938272.013           //Mass Energy of Proton in keV/c^2
#define EMASS 9.10938215e-31 		//Mass of Electron in kg
#define EMASSE  510.998910          //Mass Energy of Electron in keV/c^2
#define NMASS 1.67492721110-27      //Mass of Neutron in kg
#define NMASSE 939565.346           //Mass Energy of Neutron in keV/c^2
#define EKEMAX 781.571           //Maximum kinetic energy of electron in neutron decay
#define ETEMAX 1292.57           //Maximum total energy of electron in neturon decay
#define LITTLEA -0.103              //little a from neutron decay
#define CSPEED      2.99792458e8 			//Speed of Light in m/s
#define INVCSPEED   3.335640952e-9      //Speed of light inverse
#define CSQUARED    8.98755178736818e16 //Speed of Light squared in m/s
#define INVCSQUARED 1.11265005605362e-17 //To avoid a divide
#define FSCONST     7.2973525376e-3     //Fine structure constant
#define PI 3.14159265358979323846       //Pi...ya know math pi
#define SQRT2 1.41421356237309515 //Square root of 2
#define DEGPERRAD 57.29577951  //Degrees per Radian
#define RADPERDEG .01745329252 //Radians per Degree
#define ADIAB_POS_INCREMENT 0.00001  //In meters
#define RADDECAY_LAMBDA 1.2695
#define NLIFE 880.3                 //Neutron lifetime in seconds [ K.A. Olive et al.(Particle Data Group), Chin. Phys. C, 38, 090001 (2014) and 2015 update. (URL: http://pdg.lbl.gov)]
//#define NLIFE 885.7                 //Neutron lifetime in seconds (old)
#define GVCONSTANT 1.14962e-17
#define PLANKCONSTANT 4.13567e-18 //Plank's Constant in keV * s
//Tracker specific constant
#define INITIAL_MAX_TRACK_SIZE  500000
#define ULTIMATE_MAX_TRACK_SIZE 500000	//Limiter incase of too long of a track
#define PROTON_MAX_TIME 30e-6 //No need to run protons longer
#define ELECTRON_MAX_TIME 25e-8 //No need to run electrons longer
#define SAVE_EVENTS_FREQUENCY 1000
#define PROTON_STEPSIZE_DEFAULT 3e-10  //300 ps or 0.3 ns
#define ELECTRON_STEPSIZE_DEFAULT 1e-14  //100fs or 0.1 ps or 0.0001 ns (temp changed to 1e-14 from 1e-13)
#define NUM_EFIELDS 3
#define NUM_BFIELDS 3
#define DEFAULT_RESULTS_TREENAME "g"
#define DEFAULT_EVENTS_TREENAME "t"
#define DEFAULT_EVENTSDERV_TREENAME "ed"
#define DEFAULT_RESULTSDERV_TREENAME "rd"
#define DEFAULT_EXP_TREENAME "expTree"
#define DEFAULT_LIL_B_TREENAME "lilbtree"
#define SIG_SAVEEXIT 1
#define EVENTGEN_Z_MAX .62
#define EVENTGEN_Z_MIN -.18
#define BEAM_DIVERGENCE 0.01
//#define SBD_ENERGY_SIGMA 5. //in keV from Sn width 10/20/10
#define SBD_ENERGY_SIGMA 2.3 //in keV from Jeff 14/03
#define FLUXMAP_POS_Z 45.68
#define REALLY_BIG_DBL 1e99
#define ABSORB_DATA_NUM 501
#define ABSORB_DATA_START 0.1
#define ABSORB_DATA_SPACING 0.0598
//Hit codes

enum HitCode {HIT_NOTHING=0,
			HIT_DETECTOR=1, //Hits Surface Barrier Detector
			HIT_VOID=2,
			HIT_WALL=3,
			HIT_MIRROR=4,
			HIT_DOWNBEAM=5,
			HIT_NEAR=6,
			HIT_TOOSLOW=7,
			HIT_TOOFAST=8,
			HIT_SBDPOTCUT=9,
			HIT_SBDMINCUT=10,
			HIT_BACKSCATTER=11,
			HIT_OSCILLATION=12,
			HIT_ERROR=13,
			HIT_NOTRUN=14,
			HIT_BARE1=15,
			HIT_BARE2=16,
			HIT_BARE3=17,
			HIT_BGO1=18,
			HIT_BGO2=19,
			HIT_BGO3=20,
			HIT_BGO4=21,
			HIT_BGO5=22,
			HIT_BGO6=23,
			HIT_BGO7=24,
			HIT_BGO8=25,
			HIT_SIO2ABSORB=26,
			HIT_SITRANS=27
			};
#define MAXHITCODE 27
//Enums
enum TransMethod {RKCONSTANT,RKDYNAMIC,ADIABATIC};
enum ParType {NOTYPE,PROTON,ELECTRON,PHOTON};
const ParType GAMMA=PHOTON;
enum EveType {THREEBODY,FOURBODY};
enum TreeType {EVENT,RESULT,EVENTDERIVED,RESULTDERIVED};
enum BSModel{BS_BRIAN,BS_GEANTPEN,BS_GEANTLIV,BS_MCNP_ROB};
enum CoDet{DET_NONE,DET_EP,DET_EPG,DET_EPBG,DET_G,DET_BG};
//Analysis
#define SBD_BIAS -25000  //in Volts
#define PROTON_ARRIVAL_MIN 2.02e-6 //in sec
#define PROTON_ARRIVAL_MAX 24.98e-6  //in sec
#define ELECTRON_COMMON_POT 1400.
#define ELECTRON_MIN_ENERGY 50	    //in keV
#define ELECTRON_MAX_ENERGY 800	    //in keV
#define PROTON_MIN_ENERGY 7 //in keV
#define PROTON_MAX_ENERGY 31 //in keV
#define BGO_MULTIPLICITY_CUT_ENERGY 10 //in keV
#define BEND_ANGLE 9.5 //Angle of bend in degrees
#define COS_BEND 0.98628560153723140782496568567775  //Angle is 9.5 degrees
#define SIN_BEND 0.16504760586067764838263370078844
#define DEAD_DETECTOR_ENERGY -100000



#define MJBLAPTOP

//Folders
#ifdef MJBLAPTOP
#define EVENTS_DIR "/media/mjbales/mjbexternal/rdk/events/"
#define LITTLEB_DIR "/media/mjbales/mjbexternal/rdk/littleb/"
#define RESULTS_DIR "/media/mjbales/mjbexternal/rdk/finalresults/"
#define GRAPHS_DIR "/home/mjbales/work/rdk/graphs/"
#define HISTS_DIR "/home/mjbales/work/rdk/hists/"
#define MACROS_DIR "/home/mjbales/work/rdk/macros/"
#define MISC_DIR "/home/mjbales/work/rdk/misc/"
#define MRKGEANTDATA_DIR "/home/mjbales/work/rdk/mrkgeantrun/"
#define EXPANALYSIS_DIR "/media/mjbales/mjbexternal/rdk/analysis/"
#define CALANALYSIS_DIR "/home/mjbales/work/rdk/calanalysis/"
#define ANALYSIS_DIR "/home/mjbales/work/rdk/analysis/"
#endif

#ifndef MJBLAPTOP
#define EVENTS_DIR "/home/MJBSpace/MRK/Run/"
#define LITTLEB_DIR "/home/MJBSpace/MRK/Run/"
#define RESULTS_DIR "/home/MJBSpace/MRK/Run/"
#define GRAPHS_DIR "/home/MJBSpace/MRK/Run/"
#define HISTS_DIR "/home/MJBSpace/MRK/Run/"
#define MACROS_DIR "/home/MJBSpace/MRK/Run/"
#define MISC_DIR "/home/MJBSpace/MRK/Run/"
#define MRKGEANTDATA_DIR "/home/MJBSpace/MRK/Run/"
#define EXPANALYSIS_DIR "/home/MJBSpace/MRK/expanalysis/"
#define CALANALYSIS_DIR "/home/MJBSpace/MRK/calanalysis/"
#endif

typedef struct {
   double EDepELow;
   double EDepEHigh;
   double ToFPLow;
   double ToFPHigh;
   double EDepGLow;
   double EDepGHigh;
} AnalysisCuts;


#endif // CONSTANTS_H_INCLUDED
