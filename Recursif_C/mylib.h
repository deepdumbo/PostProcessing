#ifndef MAIN_H
#define MAIN_H

#include <math.h>

double RRG;
double RRGDur;
double EchoSpacing;
double StartPt;
double PVM_AcquisitionTime;

void Appel_Recursif(double, double, double, double, double, double, double*, double*, int);
double GetFirstEcho(double, double);
void UpdateEchoSpacing(void);
double fmodulus(double, double);

#endif
