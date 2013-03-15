/*
 *  BioCro/src/BioCro.h by Fernando Ezequiel Miguez  Copyright (C) 2007-2008
 *
 */

/* New Functions in caneGro*/

struct dbp_sugarcane_str{
	double kLeaf;
	double kStem;
	double kRoot;
	double kSeedcane;
	double kSugar;
	double kFiber;
	double kTassel;

};

struct soilML_str soilML_rootfront(double precipit, double transp, double *cws, double soildepth,
			 double *depths, double fieldc, double wiltp, double phi1, double phi2,
                         struct soilText_str soTexS, int wsFun, int layers, double rootDB,
			 double LAI, double k, double AirTemp, double IRad, double winds, double RelH,
			 int hydrDist, double rfl, double rsec, double rsdf,int optiontocalculaterootdepth, double rootfrontvelocity ,double dap);


struct dbp_sugarcane_str SUGARCANE_DBP_CUADRA(double TT, double TT0,double TTseed,double Tmaturity, double Rd, double Alm, double Arm, double Clstem, double Ilstem, double Cestem, double Iestem, double Clsuc, double Ilsuc, double Cesuc, double Iesuc,double Temperature);


double GrowthRespiration(double CanopyA, double fraction);


double MRespiration (double W,double Q, double m,double T, double deltime);


double seasonal (double maxLN, double minLN, double day, double daymaxLN, double dayinyear, double lat);