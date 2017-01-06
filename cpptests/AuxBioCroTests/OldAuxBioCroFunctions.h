#ifndef OLD_AUXBIOCRO_H
#define OLD_AUXBIOCRO_H

#include "../../src/BioCro.h"

Light_profile OldsunML(double Idir, double Idiff, double LAI, int nlayers,
                       double cosTheta, double kd, double chil, double heightf);

struct ET_Str OldEvapoTrans(double Rad, double Itot, double Airtemperature,
                            double RH, double WindSpeed, double LeafAreaIndex,
                            double CanopyHeight, double StomataWS, int ws,
                            double vmax2, double alpha2, double kparm,
                            double theta, double beta, double Rd2, double b02,
                            double b12, double upperT, double lowerT,
                            double Catm);

struct ET_Str OldEvapoTrans2(double Rad, double Iave, double Airtemperature,
                             double RH, double WindSpeed, double LeafAreaIndex,
                             double CanopyHeight, double stomatacond,
                             double leafw, int eteq);

struct ws_str Oldwatstr(double precipit, double evapo, double cws,
                        double soildepth, double fieldc, double wiltp,
                        double phi1, double phi2, int soiltype, int wsFun);

#endif
