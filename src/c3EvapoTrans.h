#ifndef C3EVAPOTRANS_H
#define C3EVAPOTRANS_H
/*
 *  /src/c3EvapoTrans.c by Fernando Ezequiel Miguez  Copyright (C) 2010
 *
 *  Part of the code here (sunML, EvapoTrans, SoilEvapo, TempTo and
 *  the *prof functions) are based on code in WIMOVAC. WIMOVAC is
 *  copyright of Stephen Long and Stephen Humphries.
 *  Documentation for WIMOVAC can be found at
 *  http://www.life.illinois.edu/plantbio/wimovac/ (checked 02-13-2010)
 *
 */

/* Decalring functions used in this function */

double TempToDdryA(double Temp);
double TempToLHV(double Temp);
double TempToSFS(double Temp);
double TempToSWVC(double Temp);

#endif

