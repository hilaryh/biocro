/*
 *  BioCro/src/BioCro.c by Fernando Ezequiel Miguez  Copyright (C) 2007-2011
 *
 *
 */

#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "c3photo.h"
#include "AuxBioCro.h"
#include "Century.h"
#include "BioCro.h"
#include "AuxwillowGro.h"
#include "AuxcaneGro.h"
#include "crocent.h"
#include "c3canopy.h"

SEXP CropGro(SEXP LAT,                 /* Latitude                  1 */ 
      SEXP DOY,                 /* Day of the year           2 */
	    SEXP HR,                  /* Hour of the day           3 */
	    SEXP SOLAR,               /* Solar Radiation           4 */
	    SEXP TEMP,                /* Temperature               5 */
	    SEXP RH,                  /* Relative humidity         6 */ 
	    SEXP WINDSPEED,           /* Wind Speed                7 */ 
	    SEXP PRECIP,              /* Precipitation             8 */
	    SEXP KD,                  /* K D (ext coeff diff)      9 */
	    SEXP CHILHF,              /* Chi and Height factor    10 */
	    SEXP NLAYERS,             /* # Lay canop              11 */
	    SEXP RHIZOME,             /* Ini Rhiz                 12 */
      SEXP IRTL,                /* i rhiz to leaf           13 */
      SEXP SENCOEFS,            /* sene coefs               14 */
      SEXP TIMESTEP,            /* time step                15 */
      SEXP VECSIZE,             /* vector size              16 */
	    SEXP SPLEAF,              /* Spec Leaf Area           17 */
	    SEXP SPD,                 /* Spec Lefa Area Dec       18 */
	    SEXP DBPCOEFS,            /* Dry Bio Coefs            19 */
	    SEXP THERMALP,            /* Themal Periods           20 */
	    SEXP VMAX,                /* Vmax of photo            21 */
	    SEXP ALPHA,               /* Quantum yield            22 */
	    SEXP KPARM,               /* k parameter (photo)      23 */
	    SEXP THETA,               /* theta param (photo)      24 */
	    SEXP BETA,                /* beta param  (photo)      25 */
	    SEXP RD,                  /* Dark Resp   (photo)      26 */
	    SEXP CATM,                /* CO2 atmosph              27 */
	    SEXP B0,                  /* Int (Ball-Berry)         28 */
	    SEXP B1,                  /* Slope (Ball-Berry)       29 */
	    SEXP WS,                  /* Water stress flag        30 */
	    SEXP SOILCOEFS,           /* Soil Coefficients        31 */
	    SEXP ILEAFN,              /* Ini Leaf Nitrogen        32 */
	    SEXP KLN,                 /* Decline in Leaf Nitr     33 */
	    SEXP VMAXB1,              /* Effect of N on Vmax      34 */
	    SEXP ALPHAB1,             /* Effect of N on alpha     35 */
	    SEXP MRESP,               /* Maintenance resp         36 */
	    SEXP SOILTYPE,            /* Soil type                37 */
	    SEXP WSFUN,               /* Water Stress Func        38 */
	    SEXP CENTCOEFS,           /* Century coefficients     39 */
      SEXP CENTTIMESTEP,        /* Century timestep         40 */ 
	    SEXP CENTKS,              /* Century decomp rates     41 */
	    SEXP SOILLAYERS,          /* # soil layers            42 */
	    SEXP SOILDEPTHS,          /* Soil Depths              43 */
	    SEXP CWS,                 /* Current water status     44 */
	    SEXP HYDRDIST,            /* Hydraulic dist flag      45 */
	    SEXP SECS,                /* Soil empirical coefs     46 */
	    SEXP KPLN,                /* Leaf N decay             47 */
	    SEXP LNB0,                /* Leaf N Int               48 */
	    SEXP LNB1,                /* Leaf N slope             49 */
      SEXP LNFUN,               /* Leaf N func flag         50 */
      SEXP UPPERTEMP,           /* Temperature Limitations photoParms */
	    SEXP LOWERTEMP,
	    SEXP NNITROP)           /*temperature Limitation photoParms */
{
    int vecsize = INTEGER(VECSIZE)[0];
    int dailyvecsize = vecsize/24;
//    Rprintf("%i\n",vecsize);
   /*********** CROCENT VARIABLES***********************/
   struct cropcentlayer CROPCENT;
   int woody, Eflag;
   woody = 0 ; // No woody Material for now
   Eflag = 1; // For N simulations only
   double *fake;
   // Get Defaukt parameters for miscanthus
     assignParms(&CROPCENT, fake);
  // Get Initial Values oof Pool for C and N
     assignPools(&CROPCENT, fake);
  // More parameters   
     GetBioCroToCropcentParms(&CROPCENT.BcroTOCentParms,fake);
   // Timestep is alreadt set to 1440.0 minutes (1 day) in the assignParms. We need to change the parameters to daily time step
     CROPCENTTimescaling(&CROPCENT);
     
   struct InputToCropcent leaflitter,stemlitter,rootlitter,rhizomelitter;
   struct minerals leaflitterE,stemlitterE,rootlitterE,rhizomelitterE;
   // The below parameters aee RCESTR from fix.100 representing CE Ratio of structural material
   leaflitterE.CN=200.0;
   leaflitterE.CP=500.0;
   leaflitterE.CS=500.0;
   leaflitterE.CS=500.0;
   stemlitterE = leaflitterE;
   rootlitterE=  leaflitterE;
   rhizomelitterE=leaflitterE;
   /*************************************************/
   
   struct crop_phenology cropdbp;
   static struct miscanthus miscanthus, deltamiscanthus;
   createNULLmiscanthus(&miscanthus,vecsize);

//   miscanthus.leafvec[vecsize].newbiomass=(double)vecsize;
//   Rprintf("%f, %i\n",miscanthus.leafvec[vecsize].newbiomass, vecsize);
   double dailynetassim, CanopyAGross, dailyGPP;
   double mrespLeaf, mrespStem, mrespRoot, mrespRhizome;
   struct senthermaltemp senthermaltemp;
   struct canopyparms canopyparms;
   struct frostParms frostparms;
   dailynetassim=0.0;
   dailyGPP=0.0;
   mrespLeaf=0.0;
   mrespStem=0.0;
   mrespRoot=0.0;
   mrespRhizome=0.0;
   frostparms.leafT0=-8.0;
   frostparms.leafT100=-16.0;
   frostparms.stemT0=-16.0;
   frostparms.stemT100=-16.0;
   
   
   /***********Management Variables *************/
   struct management management;
   assignManagement(&management);
   /********************************************/
   
	/* This creates vectors which will collect the senesced plant
	   material. This is needed to calculate litter and therefore carbon
	   in the soil and then N in the soil. */

	 double upperT=REAL(UPPERTEMP)[0];
	 double lowerT=REAL(LOWERTEMP)[0];
/*Reading NitroP Variables */
	struct nitroParms nitroparms;
	double TEMPdoubletoint;
	nitroparms.ileafN=REAL(NNITROP)[0];
  nitroparms.kln=REAL(NNITROP)[1];
	nitroparms.Vmaxb1=REAL(NNITROP)[2];
	nitroparms.Vmaxb0=REAL(NNITROP)[3];
	nitroparms.alphab1=REAL(NNITROP)[4];
	nitroparms.alphab0=REAL(NNITROP)[5];
  nitroparms.Rdb1=REAL(NNITROP)[6];
	nitroparms.Rdb0=REAL(NNITROP)[7];
	nitroparms.kpLN=REAL(NNITROP)[8];
	nitroparms.lnb0=REAL(NNITROP)[9];
	nitroparms.lnb1=REAL(NNITROP)[10];
	TEMPdoubletoint=REAL(NNITROP)[11];
	nitroparms.lnFun=(int)TEMPdoubletoint;
	nitroparms.maxln=REAL(NNITROP)[12];
	nitroparms.minln=REAL(NNITROP)[13];
	nitroparms.daymaxln=REAL(NNITROP)[14];


///////////////////////////////////////////////////////////////////	
  double iSp, Sp , propLeaf;
	int i, i2, i3;

	double vmax1;
	double alpha1;
	double kparm1;
	double theta;
	double beta;
	double Rd1, Ca;
	double b01, b11;

	double Leaf, Stem, Root, LAI, Grain = 0.0;
	double TTc = 0.0;
	double kLeaf = 0.0, kStem = 0.0, kRoot = 0.0, kRhizome = 0.0, kGrain = 0.0;
	double newLeaf, newStem = 0.0, newRoot, newRhizome, newGrain = 0.0;

	/* Variables needed for collecting litter */
	double LeafLitter = REAL(CENTCOEFS)[20], StemLitter = REAL(CENTCOEFS)[21];
	double RootLitter = REAL(CENTCOEFS)[22], RhizomeLitter = REAL(CENTCOEFS)[23];
	double LeafLitter_d = 0.0, StemLitter_d = 0.0;
	double RootLitter_d = 0.0, RhizomeLitter_d = 0.0;
	double ALitter = 0.0, BLitter = 0.0;
	/* Maintenance respiration */

	double mrc1 = REAL(MRESP)[0];
	double mrc2 = REAL(MRESP)[1]; 

	double waterCont;
	double StomWS = 1, LeafWS = 1;
	int timestep;
	double CanopyA, CanopyT;

	double Rhizome;

	/* Soil Parameters*/
	double FieldC, WiltP, phi1, phi2, soilDepth;
	int soilType, wsFun;
	double LeafN, LeafN_0, kLN;
	double soilEvap, TotEvap;
	int soillayers = INTEGER(SOILLAYERS)[0];
	double cwsVec[soillayers];
	for(i2=0;i2<soillayers;i2++){
		cwsVec[i2] = REAL(CWS)[i2];
	}
	double cwsVecSum = 0.0;
	/* Some soil related empirical coefficients */
	double rfl = REAL(SECS)[0];  /* root factor lambda */
	double rsec = REAL(SECS)[1]; /* radiation soil evaporation coefficient */
	double rsdf = REAL(SECS)[2]; /* root soil depth factor */
	double scsf = REAL(SOILCOEFS)[6]; /* stomatal conductance sensitivity factor */ /* Rprintf("scsf %.2f",scsf); */
	double transpRes = REAL(SOILCOEFS)[7]; /* Resistance to transpiration from soil to leaf */
	double leafPotTh = REAL(SOILCOEFS)[8]; /* Leaf water potential threshold */

        /* Parameters for calculating leaf water potential */
	double LeafPsim = 0.0;

        /* Effect of Nitrogen */
	double kpLN = REAL(KPLN)[0];
	double lnb0 = REAL(LNB0)[0]; 
	double lnb1 = REAL(LNB1)[0];
	int lnfun = INTEGER(LNFUN)[0];

	/* Century */
	double MinNitro = REAL(CENTCOEFS)[19];
	int doyNfert = REAL(CENTCOEFS)[18];
	double Nfert;
	double SCCs[9];
	double Resp = 0.0;
	int centTimestep = INTEGER(CENTTIMESTEP)[0];

	double SeneLeaf, SeneStem, SeneRoot = 0.0, SeneRhizome = 0.0 ;
	double *sti , *sti2, *sti3, *sti4; 
	double Remob;
	int k = 0, q = 0, m = 0, n = 0;
	int ri = 0;

	struct Can_Str Canopy;
	struct ws_str WaterS;
	struct dbp_str dbpS;
	struct cenT_str centS; 
	struct soilML_str soilMLS;
	struct soilText_str soTexS; /* , *soTexSp = &soTexS; */
	soTexS = soilTchoose(INTEGER(SOILTYPE)[0]);

	centS.SCs[0] = 0.0;
	centS.SCs[1] = 0.0;
	centS.SCs[2] = 0.0;
	centS.SCs[3] = 0.0;
	centS.SCs[4] = 0.0;
	centS.SCs[5] = 0.0;
	centS.SCs[6] = 0.0;
	centS.SCs[7] = 0.0;
	centS.SCs[8] = 0.0;
	centS.Resp = 0.0;

	SEXP lists, names;

	SEXP DayofYear;
	SEXP Hour;
	SEXP CanopyAssim;
	SEXP CanopyTrans;
	SEXP Leafy;
	SEXP Stemy;
	SEXP Rooty;
	SEXP Rhizomey;
	SEXP Grainy;
	SEXP LAIc;
	SEXP TTTc;
	SEXP SoilWatCont;
	SEXP StomatalCondCoefs;
	SEXP LeafReductionCoefs;
	SEXP LeafNitrogen;
	SEXP AboveLitter;
	SEXP BelowLitter;
	SEXP VmaxVec;
	SEXP AlphaVec;
	SEXP SpVec;
	SEXP MinNitroVec;
	SEXP RespVec;
	SEXP SoilEvaporation;
	SEXP cwsMat;
	SEXP psimMat; /* Holds the soil water potential */
	SEXP rdMat;
	SEXP SCpools;
	SEXP SNpools;
	SEXP LeafPsimVec;
// From here, we have daily (instead of hourly) output vectors
  SEXP DayafterPlanting;
  SEXP GDD; // daily vector of growing degree day
  SEXP GPP; // Gross Primary Productivity
  SEXP NPP; // Net Primary Productivity
  SEXP autoRESP; // Autotrophic Respiration
  SEXP hetRESP; // Heterotrophic Respiration
  SEXP NER; // Net Ecosystem Respiration
  SEXP StemMResp;
  SEXP RootMResp;
  SEXP RhizomeMResp;
  SEXP LeafDarkResp;
  SEXP Stemd;
  SEXP Leafd;
  SEXP Rootd;
  SEXP Rhizomed;
  SEXP Stemlitterd;
  SEXP Leaflitterd;
  SEXP Rootlitterd;
  SEXP Rhizomelitterd;
  SEXP LAId;
  SEXP totalSOC;
  SEXP strucc1;
  SEXP strucc2;
  SEXP metabc1;
  SEXP metabc2;
  SEXP som1c1;
  SEXP som1c2;
  SEXP som2c1;
  SEXP som2c2;
  SEXP som3c;
  
// Declaring Daily variables
 double  accumulatedGDD=0.0;





//	vecsize = length(DOY);
	PROTECT(lists = allocVector(VECSXP,59));
	PROTECT(names = allocVector(STRSXP,59));

	PROTECT(DayofYear = allocVector(REALSXP,vecsize));
	PROTECT(Hour = allocVector(REALSXP,vecsize));
	PROTECT(CanopyAssim = allocVector(REALSXP,vecsize));
	PROTECT(CanopyTrans = allocVector(REALSXP,vecsize));
	PROTECT(Leafy = allocVector(REALSXP,vecsize));
	PROTECT(Stemy = allocVector(REALSXP,vecsize));
	PROTECT(Rooty = allocVector(REALSXP,vecsize));
	PROTECT(Rhizomey = allocVector(REALSXP,vecsize));
	PROTECT(Grainy = allocVector(REALSXP,vecsize));
	PROTECT(LAIc = allocVector(REALSXP,vecsize));
	PROTECT(TTTc = allocVector(REALSXP,vecsize));
	PROTECT(SoilWatCont = allocVector(REALSXP,vecsize));
	PROTECT(StomatalCondCoefs = allocVector(REALSXP,vecsize));
	PROTECT(LeafReductionCoefs = allocVector(REALSXP,vecsize));
	PROTECT(LeafNitrogen = allocVector(REALSXP,vecsize));
	PROTECT(AboveLitter = allocVector(REALSXP,vecsize));
	PROTECT(BelowLitter = allocVector(REALSXP,vecsize));
	PROTECT(VmaxVec = allocVector(REALSXP,vecsize));
	PROTECT(AlphaVec = allocVector(REALSXP,vecsize));
	PROTECT(SpVec = allocVector(REALSXP,vecsize));
	PROTECT(MinNitroVec = allocVector(REALSXP,vecsize));
	PROTECT(RespVec = allocVector(REALSXP,vecsize));
	PROTECT(SoilEvaporation = allocVector(REALSXP,vecsize));
	PROTECT(cwsMat = allocMatrix(REALSXP,soillayers,vecsize));
	PROTECT(psimMat = allocMatrix(REALSXP,soillayers,vecsize));
	PROTECT(rdMat = allocMatrix(REALSXP,soillayers,vecsize));
	PROTECT(SCpools = allocVector(REALSXP,9));
	PROTECT(SNpools = allocVector(REALSXP,9));
	PROTECT(LeafPsimVec = allocVector(REALSXP,vecsize));
  PROTECT(DayafterPlanting = allocVector(REALSXP,dailyvecsize));
  PROTECT(GDD = allocVector(REALSXP,dailyvecsize));
  PROTECT(GPP = allocVector(REALSXP,dailyvecsize));
  PROTECT(NPP = allocVector(REALSXP,dailyvecsize));
  PROTECT(autoRESP = allocVector(REALSXP,dailyvecsize));
  PROTECT(hetRESP = allocVector(REALSXP,dailyvecsize));
  PROTECT(NER = allocVector(REALSXP,dailyvecsize));
  PROTECT(StemMResp= allocVector(REALSXP,dailyvecsize));
  PROTECT(RootMResp = allocVector(REALSXP,dailyvecsize));
  PROTECT(RhizomeMResp = allocVector(REALSXP,dailyvecsize));
  PROTECT(LeafDarkResp = allocVector(REALSXP,dailyvecsize));
  PROTECT(Stemd = allocVector(REALSXP,dailyvecsize));
  PROTECT(Leafd = allocVector(REALSXP,dailyvecsize));
  PROTECT(Rootd = allocVector(REALSXP,dailyvecsize));
  PROTECT(Rhizomed = allocVector(REALSXP,dailyvecsize));
  PROTECT(Stemlitterd = allocVector(REALSXP,dailyvecsize));
  PROTECT(Leaflitterd = allocVector(REALSXP,dailyvecsize));
  PROTECT(Rootlitterd = allocVector(REALSXP,dailyvecsize));
  PROTECT(Rhizomelitterd = allocVector(REALSXP,dailyvecsize));
  PROTECT(LAId = allocVector(REALSXP,dailyvecsize));
  PROTECT(totalSOC = allocVector(REALSXP,dailyvecsize));
  PROTECT(strucc1 = allocVector(REALSXP,dailyvecsize));
   PROTECT(strucc2 = allocVector(REALSXP,dailyvecsize));
  PROTECT(metabc1 = allocVector(REALSXP,dailyvecsize));
  PROTECT(metabc2 = allocVector(REALSXP,dailyvecsize));
   PROTECT(som1c1 = allocVector(REALSXP,dailyvecsize));
  PROTECT(som1c2 = allocVector(REALSXP,dailyvecsize));
  PROTECT(som2c1 = allocVector(REALSXP,dailyvecsize));
   PROTECT(som2c2 = allocVector(REALSXP,dailyvecsize));
  PROTECT(som3c = allocVector(REALSXP,dailyvecsize));
	/* Picking vmax, alpha and kparm */
	vmax1 = REAL(VMAX)[0];
	alpha1 = REAL(ALPHA)[0];
	kparm1 = REAL(KPARM)[0];
	theta = REAL(THETA)[0];
	beta = REAL(BETA)[0];
	Rd1 = REAL(RD)[0];
	Ca = REAL(CATM)[0];
	b01 = REAL(B0)[0];
	b11 = REAL(B1)[0];

	LeafN_0 = REAL(ILEAFN)[0];
	LeafN = LeafN_0; /* Initial value of N in the leaf */
	kLN = REAL(KLN)[0];
	timestep = INTEGER(TIMESTEP)[0];

	Rhizome = REAL(RHIZOME)[0];
	Sp = REAL(SPLEAF)[0]; 
	SeneLeaf = REAL(SENCOEFS)[0];
	SeneStem = REAL(SENCOEFS)[1];
	SeneRoot = REAL(SENCOEFS)[2];
	SeneRhizome = REAL(SENCOEFS)[3];

	/* Soil Parameters */
	FieldC = REAL(SOILCOEFS)[0];
	WiltP = REAL(SOILCOEFS)[1];
	phi1 = REAL(SOILCOEFS)[2];
	phi2 = REAL(SOILCOEFS)[3];
	soilDepth = REAL(SOILCOEFS)[4];
	waterCont = REAL(SOILCOEFS)[5];
	wsFun = INTEGER(WSFUN)[0];
	soilType = INTEGER(SOILTYPE)[0];

	SCCs[0] = REAL(CENTCOEFS)[0];
	SCCs[1] = REAL(CENTCOEFS)[1];
	SCCs[2] = REAL(CENTCOEFS)[2];
	SCCs[3] = REAL(CENTCOEFS)[3];
	SCCs[4] = REAL(CENTCOEFS)[4];
	SCCs[5] = REAL(CENTCOEFS)[5];
	SCCs[6] = REAL(CENTCOEFS)[6];
	SCCs[7] = REAL(CENTCOEFS)[7];
	SCCs[8] = REAL(CENTCOEFS)[8];


/* Creating pointers to avoid calling functions REAL and INTEGER so much */
//  int *pt_doy=INTEGER(DOY);
	  int *pt_doy;
    pt_doy = malloc(vecsize*sizeof(pt_doy));
    pt_doy=INTEGER(DOY);
  
// 	int *pt_hr = INTEGER(HR);
    int *pt_hr;
    pt_hr = malloc(vecsize*sizeof(pt_hr));
    pt_hr=INTEGER(HR);
    
//	double *pt_solar = REAL(SOLAR);
   double *pt_solar;
    pt_solar = malloc(vecsize*sizeof(pt_solar));
    pt_solar=REAL(SOLAR);
  
//	double *pt_temp = REAL(TEMP);
    double *pt_temp;
    pt_temp = malloc(vecsize*sizeof(pt_temp));
    pt_temp=REAL(TEMP);
    
//	double *pt_rh = REAL(RH);
    double *pt_rh;
    pt_rh = malloc(vecsize*sizeof(pt_rh));
    pt_rh=REAL(RH);
//	double *pt_windspeed = REAL(WINDSPEED);
   double *pt_windspeed;
    pt_windspeed = malloc(vecsize*sizeof(pt_windspeed));
    pt_windspeed=REAL(WINDSPEED);
//	double *pt_precip = REAL(PRECIP);
   double *pt_precip;
    pt_precip = malloc(vecsize*sizeof(pt_precip));
    pt_precip=REAL(PRECIP);


  
   
	double lat = REAL(LAT)[0];
	int nlayers = INTEGER(NLAYERS)[0];
	int ws = INTEGER(WS)[0];
	double kd = REAL(KD)[0];
	double chil = REAL(CHILHF)[0];
	double hf = REAL(CHILHF)[1];

 
 /***Initialize Daily Variables *********/
  struct respirationParms RESP;
 RESP.growth.stem=0.3;
 RESP.growth.root=0.3;
 RESP.growth.rhizome=0.3;
 
 RESP.maint.Qstem=2.0;
 RESP.maint.mstem=0.004;
 
 RESP.maint.Qroot=2.0;
 RESP.maint.mroot=0.002;
 
 RESP.maint.Qrhizome=2.0;
 RESP.maint.mrhizome=0.002;
//  Resp.growth.stem=0.3;
  
  double LeafResp,StemResp,RootResp,RhizResp;  
  double gRespCoeff = 0.0;
  double dailydelTT = 0.0;
  double delTT;
  double Tbase=0.0;
  dailynetassim=0.0;
  senthermaltemp.leafcriticalT=REAL(SENCOEFS)[0];
  senthermaltemp.stemcriticalT=REAL(SENCOEFS)[1];
  senthermaltemp.rootcriticalT=REAL(SENCOEFS)[2];
  senthermaltemp.rhizomecriticalT=REAL(SENCOEFS)[3];
  canopyparms.kN=0.1;
  canopyparms.SLA=0.1;
  canopyparms.remobFac=0.1;
  canopyparms.leafNsen=40;
  frostparms.leafT0=-20.0; //REAL(FROSTP)[0];
  frostparms.leafT100=-20.0;//REAL(FROSTP)[1];
  frostparms.stemT0=-20.0;//REAL(FROSTP)[2];
  frostparms.stemT100=-20.0;//REAL(FROSTP)[3];
 
  propLeaf = REAL(IRTL)[0]; 
	/* It is useful to assume that there is a small amount of
	   leaf area at the begining of the growing season. */
//	Leaf = Rhizome * 0.001; 
	/* Initial proportion of the rhizome that is turned
	   into leaf the first hour */
//	Stem = Rhizome * 0.001;
//	Root = Rhizome * 0.001;
	  /**********Assining Canopy Parameters********************/
      canopyparms.remobFac=0.5;
      int dap=0;
      /*******************************************************/
  miscanthus.leaf.biomass=0.0;
  miscanthus.stem.biomass=0.0;
  miscanthus.root.biomass=0.0;
  miscanthus.rhizome.biomass=management.emergenceparms.plantingrate;
 
  int emergence=0;
  struct dailyclimate dailyclimate;
  TTc=0.0;
  REAL(TTTc)[0]=TTc;
  
  
 /**************************************/
  updateafteremergence(&miscanthus,&management);
  LAI = miscanthus.leaf.biomass*Sp;
  int phototype;
  phototype=1;
  // This is specific to willow to avoid harvesting based on day of year
  
	for(i=0;i<vecsize;i++)
//    for(i=0;i<3;i++)
	{
		/* First calculate the elapsed Thermal Time*/
		/* The idea is that here I need to divide by the time step
		   to calculate the thermal time. For example, a 3 hour time interval
		   would mean that the division would need to by 8 */
//       delTT=*(pt_temp+i) / (24/timestep);
        delTT=getThermaltime(*(pt_temp+i), Tbase);
        delTT=delTT/24;
//    dailydelTT+=delTT;
//    Rprintf("index=%i,temp=%f, delTT= %f\n", i,*(pt_temp+i),delTT);
//         LAI=6.0;
        if(emergence==0)
            {
            TTc +=delTT;
            REAL(TTTc)[i] =REAL(TTTc)[i-1]+delTT ;
            CanopyA = 0.0;
            CanopyAGross =0.0;
        		CanopyT = 0.0;
            miscanthus.autoresp.leafdarkresp=0;
            }
        else
            {
//         Rprintf("Before Canopy Function, Phototype = %i, i= %i, Assim=%f, Leaf=%f, LAI=%f, Specific Leaf Area = %f \n", phototype,i, Canopy.Assim, miscanthus.leaf.biomass, LAI,Sp);
	        	TTc +=delTT;
		        REAL(TTTc)[i] =REAL(TTTc)[i-1]+delTT ;
            
            if(phototype==1)
            {
  	        Canopy = CanAC(LAI, *(pt_doy+i), *(pt_hr+i),
			       *(pt_solar+i), *(pt_temp+i),
			       *(pt_rh+i), *(pt_windspeed+i),
			       lat, nlayers,
			       vmax1,alpha1,kparm1,
			       theta,beta,Rd1,Ca,b01,b11,StomWS,
			       ws, kd,
			       chil, hf,LeafN, kpLN, lnb0, lnb1, lnfun,upperT,lowerT,nitroparms);
            }
            if(phototype==2)
            {
             double jmax1,o2;
             vmax1=100.0;
             jmax1=180.0;
             o2=210.0;
             b01=0.08;
             b11=5.0;
             theta=0.7;
             Canopy = newc3CanAC(LAI, *(pt_doy+i), *(pt_hr+i),
  		       *(pt_solar+i), *(pt_temp+i),
			       *(pt_rh+i), *(pt_windspeed+i),
			       lat, nlayers,
			       vmax1,jmax1,Rd1,Ca,o2,b01,b11,theta,kd,hf,LeafN, kpLN, lnb0, lnb1, lnfun, StomWS, ws); 
            }
//             Rprintf("After Canopy Function,Phototype= %i, i= %i, Assim=%f, Leaf=%f, LAI=%f, Specific Leaf Area = %f \n", phototype,i, Canopy.Assim, miscanthus.leaf.biomass, LAI,Sp);
        		CanopyA = Canopy.Assim * timestep;
            CanopyAGross =Canopy.GrossAssim*timestep;
        		CanopyT = Canopy.Trans * timestep;
 //           Rprintf("NetA=%f, Gross A= %f, Trans=%f \n",CanopyA, CanopyAGross);
            }
		/* Inserting the multilayer model */
		  if(soillayers > 1)
            {
        		soilMLS = soilML(*(pt_precip+i), CanopyT, &cwsVec[0], soilDepth, REAL(SOILDEPTHS), FieldC, WiltP,
        	            phi1, phi2, soTexS, wsFun, INTEGER(SOILLAYERS)[0], miscanthus.root.biomass, 
        					    LAI, 0.68, *(pt_temp+i), *(pt_solar), *(pt_windspeed+i), *(pt_rh+i), 
        					    INTEGER(HYDRDIST)[0], rfl, rsec, rsdf);

            StomWS = soilMLS.rcoefPhoto;
            LeafWS = soilMLS.rcoefSpleaf;
            soilEvap = soilMLS.SoilEvapo;
              			for(i3=0;i3<soillayers;i3++)
                      {
              				cwsVec[i3] = soilMLS.cws[i3];
              				cwsVecSum += cwsVec[i3];
              				REAL(cwsMat)[i3 + i*soillayers] = soilMLS.cws[i3];
              				REAL(rdMat)[i3 + i*soillayers] = soilMLS.rootDist[i3];
              			  }

			      waterCont = cwsVecSum / soillayers;
			      cwsVecSum = 0.0;

		      }
      else
          {

      			soilEvap = SoilEvapo(LAI, 0.68, *(pt_temp+i), *(pt_solar+i), waterCont, FieldC, WiltP, 
                                                   *(pt_windspeed+i), *(pt_rh+i), rsec);
      			TotEvap = soilEvap + CanopyT;
      			WaterS = watstr(*(pt_precip+i),TotEvap,waterCont,soilDepth,FieldC,WiltP,phi1,phi2,soilType, wsFun);   
      			waterCont = WaterS.awc;
      			StomWS = WaterS.rcoefPhoto ; 
      			LeafWS = WaterS.rcoefSpleaf;
      			REAL(cwsMat)[i] = waterCont;
      			REAL(psimMat)[i] = WaterS.psim;
      		}

                if(wsFun == 4)
                  {
		        	    LeafPsim = WaterS.psim - (CanopyT * 1e3 * 1e-4 * 1.0/3600.0) * transpRes;
                  if(LeafPsim < leafPotTh){
				          StomWS = 1 - ((leafPotTh - LeafPsim)/1000 * scsf);
				          if(StomWS < 0.1) StomWS = 0.1;
		        	   }else{
			      	    StomWS = 1;
	          		  }
	              	}else{
	          		  LeafPsim = 0;
                }
                                    
/****************Evaluating Daily Maintenance R espiration and Gross canopy assimilation******************/

  StemResp=MRespiration(miscanthus.stem.biomass, RESP.maint.Qstem, RESP.maint.mstem, *(pt_temp+i), timestep);
  miscanthus.autoresp.stemmaint+=StemResp;
  RootResp=MRespiration(miscanthus.root.biomass, RESP.maint.Qroot, RESP.maint.mroot, *(pt_temp+i), timestep);
  miscanthus.autoresp.rootmaint+=RootResp;
  RhizResp=MRespiration(miscanthus.rhizome.biomass, RESP.maint.Qrhizome, RESP.maint.mrhizome, *(pt_temp+i), timestep);
  miscanthus.autoresp.rhizomemaint+=RhizResp;
  miscanthus.autoresp.leafdarkresp+=(CanopyAGross-CanopyA);
  dailynetassim+=CanopyA;//Net Canopy Assimilation
  miscanthus.GPP+=CanopyAGross;
 /**********************************************************************************/
//  Rprintf("%f,%f,%f,%f,%f\n",*(pt_solar+i),*(pt_temp+i),*(pt_precip+i),*(pt_rh+i),*(pt_windspeed+i));

   if(i % 24== 0)
   {
  // update daily climate structure
   getdailyclimate(&dailyclimate, pt_doy,pt_solar,pt_temp, pt_rh, pt_windspeed,pt_precip,i,vecsize);  
  // calculate GDDD for today only if crop is growing, otherwise set it to zero
   dailydelTT = (emergence ==1) ? getThermaltime(dailyclimate.temp, Tbase):0.0;
   // Update day after planting today
   REAL(DayafterPlanting)[dap]=dap;  
   // update GDD (accumulated thermal time) upto today
   accumulatedGDD+=((dap==0)? 0:dailydelTT);
//   REAL(GDD)[dap]=((dap==0)? 0:REAL(GDD)[dap-1])+dailydelTT;
     REAL(GDD)[dap]=accumulatedGDD;
     dailymiscanthus(&miscanthus, REAL(DBPCOEFS),REAL(THERMALP),accumulatedGDD, *(pt_temp+i), dailynetassim,&senthermaltemp, &canopyparms,&frostparms,i,dailydelTT,&RESP,emergence);     
        /***** Check if plant is already emerged, then growth needs to be simulated ***/ 
          if(emergence==1)
          {
            /*** Simulate miscanthus growth ***/
//              dailymiscanthus(&miscanthus, REAL(DBPCOEFS),REAL(THERMALP),accumulatedGDD, *(pt_temp+i), dailynetassim,&senthermaltemp, &canopyparms,&frostparms,i,dailydelTT,&RESP,emergence);
//              Rprintf("DPY= %i, GDD=%f\n",dailyclimate.doy,accumulatedGDD);
            /** Check if today is harvest day **/
                 if(dailyclimate.doy==management.harvestparms.harvestdoy)
                {
                  emergence=0;                        //Emergence is set back to zero
                  REAL(GDD)[dap]=0.0;                 //Set GDD back to zero to restart phenology from beginning
                  updateafterharvest(&miscanthus,&management); // Use harvest parameters to implement pracices such as removingor leaving residues 
                  Rprintf("in CropGRO, harvest day %i\n",i);
                }         
          }
          /** Here is calculation for situations when plan is not emerged ***/
          else
          {      
                  /** Check if today is emergence day, IF yes, it will set emergence =1 for next day simulation **/
                  emergence=CheckEmergence(&dailyclimate,management.emergenceparms.emergenceTemp); 
                  if((dailyclimate.doy==120)&&(phototype==2))emergence=1;
                  /** if today indeed is emergence day then we also need to initialze leaf biomass based on storage pool to initiate simulations **/
                  if(emergence==1)
                  {
//                    Rprintf("BEFORE leaf=%f, rhizome=%f \n",miscanthus.leaf.biomass,miscanthus.rhizome.biomass);
                  updateafteremergence(&miscanthus,&management);
//                         Rprintf("AFTER leaf=%f, rhizome=%f \n",miscanthus.leaf.biomass,miscanthus.rhizome.biomass);
                      Rprintf("in CropGRO, emergence day %i\n",i);
                      accumulatedGDD=0.0;
                      TTc=0.0;
//                  miscanthus.leaf.biomass=0.02*  miscanthus.rhizome.biomass;
                  LAI = miscanthus.leaf.biomass* Sp;
                  }
//                   dailymiscanthus(&miscanthus, REAL(DBPCOEFS),REAL(THERMALP),accumulatedGDD, *(pt_temp+i), dailynetassim,&senthermaltemp, &canopyparms,&frostparms,i,dailydelTT,&RESP,emergence);
                 // updatedormantstage(&miscanthus);
               
          }
                  LAI=miscanthus.leaf.biomass*Sp;
//                  Rprintf("%i,%i,%f,%f\n",emergence,dap, iSp, Sp);


/****************************************************************************/
// CROPCENT SIMULATION BEGINS HHERE    
//
// Perhaps it is better to insert a function ImplementManagement. Implement Management Functions can do following things
// Implement Harvest Management, removing biomass and leaving a fraction of litter, which can then be used to update structures for updating the input to soil biogeochemistry
// Implement N application, which will affect N available for plants and also N dynamics of soil
// Implement Irrigation, which can change soil water status, affecting plant growth and soil biogeochemistry
// Implement Tillage which will change decomposition rate of difference soil pools

//Now here We are simpyl taking litter vector and creating new structures to send to soil biogeochemistry simulations

//   BiocroToCrocent(&LeafLitter,leaf.fallrate,leaf.lignin, &leaf.E, isotoperatio, 1, 0,leaflitter);
   
    BiocroToCrocent(&miscanthus.leaf.litter,1.0,0.2, &leaflitterE, 1.0, 1, 0,&leaflitter);
    BiocroToCrocent(&miscanthus.stem.litter,1.0,0.2, &stemlitterE, 1.0, 1, 0,&stemlitter);
    BiocroToCrocent(&miscanthus.root.litter,1.0,0.2, &rootlitterE, 1.0, 0, 0,&leaflitter);
    BiocroToCrocent(&miscanthus.rhizome.litter,1.0,0.2, &rhizomelitterE, 1.0, 0, 0,&rhizomelitter);
    
    
//    Rprintf("Biomass in g C/m2=%f, CN =%f, Lignin = %f , surface=%i, woody=%i \n",leaflitter.C.totalC,leaflitter.E.CN,leaflitter.lignin,leaflitter.surface,leaflitter.woody);
     if(leaflitter.C.totalC >0.0) 
     {
      Rprintf("Before updating from litter strucc1 = %f,structCN=%f, metabcpool =%f \n", CROPCENT.strucc1.C.totalC,CROPCENT.strucc1.E.CN,CROPCENT.metabc1.C.totalC);
      UpdateCropcentPoolsFromBioCro(&CROPCENT, &leaflitter);
      Rprintf("After updating from litter strucc1 = %f, structCN=%f,metabcpool =%f \n", CROPCENT.strucc1.C.totalC, CROPCENT.strucc1.E.CN,CROPCENT.metabc1.C.totalC);
     }
      if(stemlitter.C.totalC >0.0) 
     {
      UpdateCropcentPoolsFromBioCro(&CROPCENT, &stemlitter);
     }
      if(rootlitter.C.totalC >0.0) 
     {
      UpdateCropcentPoolsFromBioCro(&CROPCENT, &rootlitter);
     }
      if(rhizomelitter.C.totalC >0.0) 
     {
      UpdateCropcentPoolsFromBioCro(&CROPCENT, &rhizomelitter);
     }
       assignENV(&CROPCENT,fake,fake,fake,fake,fake,fake,fake);
       assignFluxRatios(&CROPCENT);
       decomposeCROPCENT(&CROPCENT, woody,Eflag);
       updatecropcentpools(&CROPCENT);
       printcropcentout(CROPCENT,
                        &REAL(totalSOC)[dap],
                        &REAL(strucc1)[dap],
                        &REAL(strucc2)[dap],
                        &REAL(metabc1)[dap],
                        &REAL(metabc2)[dap],
                        &REAL(som1c1)[dap],
                        &REAL(som1c2)[dap],
                        &REAL(som2c1)[dap],
                        &REAL(som2c2)[dap],
                        &REAL(som3c)[dap]);

// 
// BiocroToCrocent(&RootLitter,root.fallrate,root.lignin, &root.E, isotoperatio, 0, 0,rootlitter);
// BiocroToCrocent(&RhizomeLitter,rhiz.fallrate,rhiz.lignin, &rhiz.E, isotoperatio, 0, 0,rhizomelitter);

// Also We need to update cropcent layer environment based on soil properties,Temperature, and soil moisture.This will influence actual decomposition rate

// Now need to use structures to change the Cadded to soils based on structures. It is like adding into soil pools based on characteristic of the litter 

// Now we can simulate decomposition of the SOC

// Now we need to revert back decomposition coefficients, which may have been modified by differen factors, such as tillage etc. so they do not keep recycled.

// Take all the outputs and send them to R for plotting purpose
/***************************************************************************/
   REAL(GPP)[dap]=miscanthus.GPP;
   REAL(LeafDarkResp)[dap]=miscanthus.autoresp.leafdarkresp;
   REAL(StemMResp)[dap]=miscanthus.autoresp.stemmaint;
   REAL(RootMResp)[dap]=miscanthus.autoresp.rootmaint;
   REAL(RhizomeMResp)[dap]=miscanthus.autoresp.rhizomemaint;
   REAL(autoRESP)[dap]= miscanthus.autoresp.total;
   miscanthus.NPP=miscanthus.GPP-miscanthus.autoresp.total;
   REAL(NPP)[dap]=miscanthus.NPP;
   REAL(Stemd)[dap]=miscanthus.stem.biomass;
   REAL(Leafd)[dap]=miscanthus.leaf.biomass;
   REAL(Rootd)[dap]=miscanthus.root.biomass;
   REAL(Rhizomed)[dap]=miscanthus.rhizome.biomass;
   REAL(Stemlitterd)[dap]=miscanthus.stem.litter;
   REAL(Leaflitterd)[dap]=miscanthus.leaf.litter;
   REAL(Rootlitterd)[dap]=miscanthus.root.litter;
   REAL(Rhizomelitterd)[dap]=miscanthus.rhizome.litter;
   REAL(LAId)[dap]=LAI;
   
    miscanthus.autoresp.leafdarkresp=0.0;
    miscanthus.autoresp.stemmaint=0.0;
    miscanthus.autoresp.rootmaint=0.0;
    miscanthus.autoresp.rhizomemaint=0.0;  
    miscanthus.GPP=0;
    miscanthus.NPP=0;
    dailynetassim=0.0;
    dailyGPP=0.0;
    dap+=1;
		}



		MinNitro = centS.MinN; /* These should be kg / m^2 per week? */
		Resp = centS.Resp;
		SCCs[0] = centS.SCs[0];
		SCCs[1] = centS.SCs[1];
		SCCs[2] = centS.SCs[2];
		SCCs[3] = centS.SCs[3];
		SCCs[4] = centS.SCs[4];
		SCCs[5] = centS.SCs[5];
		SCCs[6] = centS.SCs[6];
		SCCs[7] = centS.SCs[7];
		SCCs[8] = centS.SCs[8];



		ALitter = LeafLitter + StemLitter;
		BLitter = RootLitter + RhizomeLitter;
    
		/* Here I could add a soil and nitrogen carbon component. I have soil
		   moisture, I have temperature and root and rhizome biomass */

		REAL(DayofYear)[i] =  INTEGER(DOY)[i];
		REAL(Hour)[i] =  INTEGER(HR)[i];
		REAL(CanopyAssim)[i] =  CanopyA;
		REAL(CanopyTrans)[i] =  CanopyT; 
		REAL(Leafy)[i] = miscanthus.leaf.biomass; //Leaf;
		REAL(Stemy)[i] = miscanthus.stem.biomass; // Stem;
		REAL(Rooty)[i] =  miscanthus.root.biomass; //Root;
		REAL(Rhizomey)[i] = miscanthus.rhizome.biomass; //Rhizome;
		REAL(Grainy)[i] = Grain;
		REAL(LAIc)[i] = LAI;
		REAL(SoilWatCont)[i] = waterCont;
		REAL(StomatalCondCoefs)[i] = StomWS;
		REAL(LeafReductionCoefs)[i] = LeafWS;
		REAL(LeafNitrogen)[i] = LeafN;
		REAL(AboveLitter)[i] = ALitter;
		REAL(BelowLitter)[i] = BLitter;
		REAL(VmaxVec)[i] = vmax1;
		REAL(AlphaVec)[i] = alpha1;
		REAL(SpVec)[i] = Sp;
		REAL(MinNitroVec)[i] = MinNitro/ (24*centTimestep);
		REAL(RespVec)[i] = Resp / (24*centTimestep);
		REAL(SoilEvaporation)[i] = soilEvap;
		REAL(LeafPsimVec)[i] = LeafPsim;
    
}

/* Populating the results of the Century model */

		REAL(SCpools)[0] = centS.SCs[0];
		REAL(SCpools)[1] = centS.SCs[1];
		REAL(SCpools)[2] = centS.SCs[2];
		REAL(SCpools)[3] = centS.SCs[3];
		REAL(SCpools)[4] = centS.SCs[4];
		REAL(SCpools)[5] = centS.SCs[5];
		REAL(SCpools)[6] = centS.SCs[6];
		REAL(SCpools)[7] = centS.SCs[7];
		REAL(SCpools)[8] = centS.SCs[8];

		REAL(SNpools)[0] = centS.SNs[0];
		REAL(SNpools)[1] = centS.SNs[1];
		REAL(SNpools)[2] = centS.SNs[2];
		REAL(SNpools)[3] = centS.SNs[3];
		REAL(SNpools)[4] = centS.SNs[4];
		REAL(SNpools)[5] = centS.SNs[5];
		REAL(SNpools)[6] = centS.SNs[6];
		REAL(SNpools)[7] = centS.SNs[7];
		REAL(SNpools)[8] = centS.SNs[8];

	SET_VECTOR_ELT(lists,0,DayofYear);
	SET_VECTOR_ELT(lists,1,Hour);
	SET_VECTOR_ELT(lists,2,CanopyAssim);
	SET_VECTOR_ELT(lists,3,CanopyTrans);
	SET_VECTOR_ELT(lists,4,Leafy);
	SET_VECTOR_ELT(lists,5,Stemy);
	SET_VECTOR_ELT(lists,6,Rooty);
	SET_VECTOR_ELT(lists,7,Rhizomey);
	SET_VECTOR_ELT(lists,8,Grainy);
	SET_VECTOR_ELT(lists,9,LAIc);
	SET_VECTOR_ELT(lists,10,TTTc);
	SET_VECTOR_ELT(lists,11,SoilWatCont);
	SET_VECTOR_ELT(lists,12,StomatalCondCoefs);
	SET_VECTOR_ELT(lists,13,LeafReductionCoefs);
	SET_VECTOR_ELT(lists,14,LeafNitrogen);
	SET_VECTOR_ELT(lists,15,AboveLitter);
	SET_VECTOR_ELT(lists,16,BelowLitter);
	SET_VECTOR_ELT(lists,17,VmaxVec);
	SET_VECTOR_ELT(lists,18,AlphaVec);
	SET_VECTOR_ELT(lists,19,SpVec);
	SET_VECTOR_ELT(lists,20,MinNitroVec);
	SET_VECTOR_ELT(lists,21,RespVec);
	SET_VECTOR_ELT(lists,22,SoilEvaporation);
	SET_VECTOR_ELT(lists,23,cwsMat);
	SET_VECTOR_ELT(lists,24,psimMat);
	SET_VECTOR_ELT(lists,25,rdMat);
	SET_VECTOR_ELT(lists,26,SCpools);
	SET_VECTOR_ELT(lists,27,SNpools);
	SET_VECTOR_ELT(lists,28,LeafPsimVec);
  SET_VECTOR_ELT(lists,29,DayafterPlanting);
  SET_VECTOR_ELT(lists,30,GDD);
  SET_VECTOR_ELT(lists,31,GPP);
  SET_VECTOR_ELT(lists,32,NPP);
  SET_VECTOR_ELT(lists,33,autoRESP);  
  SET_VECTOR_ELT(lists,34,hetRESP);
  SET_VECTOR_ELT(lists,35,NER);
  SET_VECTOR_ELT(lists,36,StemMResp);
  SET_VECTOR_ELT(lists,37,RootMResp);
  SET_VECTOR_ELT(lists,38,RhizomeMResp);
  SET_VECTOR_ELT(lists,39,LeafDarkResp);
  SET_VECTOR_ELT(lists,40,Stemd);
  SET_VECTOR_ELT(lists,41,Leafd);
  SET_VECTOR_ELT(lists,42,Rootd);
  SET_VECTOR_ELT(lists,43,Rhizomed);
  SET_VECTOR_ELT(lists,44,Stemlitterd);
  SET_VECTOR_ELT(lists,45,Leaflitterd);
  SET_VECTOR_ELT(lists,46,Rootlitterd);
  SET_VECTOR_ELT(lists,47,Rhizomelitterd);
  SET_VECTOR_ELT(lists,48,LAId);
  SET_VECTOR_ELT(lists,49,totalSOC);
  SET_VECTOR_ELT(lists,50,strucc1);
  SET_VECTOR_ELT(lists,51,strucc2);
  SET_VECTOR_ELT(lists,52,metabc1);
  SET_VECTOR_ELT(lists,53,metabc1);
  SET_VECTOR_ELT(lists,54,som1c1);
  SET_VECTOR_ELT(lists,55,som1c2);
  SET_VECTOR_ELT(lists,56,som2c1);
  SET_VECTOR_ELT(lists,57,som2c2);
  SET_VECTOR_ELT(lists,58,som3c);

	SET_STRING_ELT(names,0,mkChar("DayofYear"));
	SET_STRING_ELT(names,1,mkChar("Hour"));
	SET_STRING_ELT(names,2,mkChar("CanopyAssim"));
	SET_STRING_ELT(names,3,mkChar("CanopyTrans"));
	SET_STRING_ELT(names,4,mkChar("Leaf"));
	SET_STRING_ELT(names,5,mkChar("Stem"));
	SET_STRING_ELT(names,6,mkChar("Root"));
	SET_STRING_ELT(names,7,mkChar("Rhizome"));
	SET_STRING_ELT(names,8,mkChar("Grain"));
	SET_STRING_ELT(names,9,mkChar("LAI"));
	SET_STRING_ELT(names,10,mkChar("ThermalT"));
	SET_STRING_ELT(names,11,mkChar("SoilWatCont"));
	SET_STRING_ELT(names,12,mkChar("StomatalCondCoefs"));
	SET_STRING_ELT(names,13,mkChar("LeafReductionCoefs"));
	SET_STRING_ELT(names,14,mkChar("LeafNitrogen"));
	SET_STRING_ELT(names,15,mkChar("AboveLitter"));
	SET_STRING_ELT(names,16,mkChar("BelowLitter"));
	SET_STRING_ELT(names,17,mkChar("VmaxVec"));
	SET_STRING_ELT(names,18,mkChar("AlphaVec"));
	SET_STRING_ELT(names,19,mkChar("SpVec"));
	SET_STRING_ELT(names,20,mkChar("MinNitroVec"));
	SET_STRING_ELT(names,21,mkChar("RespVec"));
	SET_STRING_ELT(names,22,mkChar("SoilEvaporation"));
	SET_STRING_ELT(names,23,mkChar("cwsMat"));
	SET_STRING_ELT(names,24,mkChar("psimMat"));
	SET_STRING_ELT(names,25,mkChar("rdMat"));
	SET_STRING_ELT(names,26,mkChar("SCpools"));
	SET_STRING_ELT(names,27,mkChar("SNpools"));
	SET_STRING_ELT(names,28,mkChar("LeafPsimVec"));
  SET_STRING_ELT(names,29,mkChar("DayafterPlanting"));
  SET_STRING_ELT(names,30,mkChar("GDD"));
  SET_STRING_ELT(names,31,mkChar("GPP"));
  SET_STRING_ELT(names,32,mkChar("NPP"));
  SET_STRING_ELT(names,33,mkChar("autoRESP"));
  SET_STRING_ELT(names,34,mkChar("hetRESP"));
  SET_STRING_ELT(names,35,mkChar("NER"));
  SET_STRING_ELT(names,36,mkChar("StemMResp"));
  SET_STRING_ELT(names,37,mkChar("RootMResp"));
  SET_STRING_ELT(names,38,mkChar("RhizomeMResp"));
  SET_STRING_ELT(names,39,mkChar("LeafDarkResp"));
  SET_STRING_ELT(names,40,mkChar("Stemd"));
  SET_STRING_ELT(names,41,mkChar("Leafd"));
  SET_STRING_ELT(names,42,mkChar("Rootd"));
  SET_STRING_ELT(names,43,mkChar("Rhizomed"));
  SET_STRING_ELT(names,44,mkChar("Stemlitterd"));
  SET_STRING_ELT(names,45,mkChar("Leaflitterd"));
  SET_STRING_ELT(names,46,mkChar("Rootlitterd"));
  SET_STRING_ELT(names,47,mkChar("Rhizomelitterd"));
  SET_STRING_ELT(names,48,mkChar("LAId"));
  SET_STRING_ELT(names,49,mkChar("totalSOC"));
  SET_STRING_ELT(names,50,mkChar("strucc1"));
  SET_STRING_ELT(names,51,mkChar("strucc2"));
  SET_STRING_ELT(names,52,mkChar("metabc1"));
  SET_STRING_ELT(names,53,mkChar("metabc1"));
  SET_STRING_ELT(names,54,mkChar("som1c1"));
  SET_STRING_ELT(names,55,mkChar("som1c2"));
  SET_STRING_ELT(names,56,mkChar("som2c1"));
  SET_STRING_ELT(names,57,mkChar("som2c2"));
  SET_STRING_ELT(names,58,mkChar("som3c"));
	setAttrib(lists,R_NamesSymbol,names);
	UNPROTECT(61);
	return(lists);
}
