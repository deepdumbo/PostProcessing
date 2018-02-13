#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mylib.h"



int main(void)
{
	printf("Bienvenue dans le main... \n\n");

	// Initialisation	
	double slew 	= 100.0;
  	double readDim 	= 128;
	int buff = 0;


	printf("Veuillez saisir le Read Dim : ");
	scanf("%d",&buff);
	readDim = (double)buff;
	printf("\n\n");

	double BP 	= 150000;
	double aqq  	= (double)(readDim/BP)*1000.0; // ms
	double TE 	= 2.0;
	double newES 	= 0;	
	double rt 	= 0.268; // ms
	//double d14  	= rt;
	double grad 	= 15.76;
	double gradDur 	= aqq + rt;
	double readMm 	= slew*(grad/100)*gradDur;
	double rewMm 	= 0;

	double durLim 	= 2*rt;
	double gradLim	= 90;
	double factor 	= 2;

	RRGDur 	= aqq + rt;
	RRG  	= grad;
	PVM_AcquisitionTime = aqq; 
	StartPt = 0.7; // ms
	
	struct output{
		double rewDur;
		double rewGrad;
	}out;

	out.rewDur  = RRGDur;
	out.rewGrad = RRG;

	// Print initial values
	printf(" ---------------------------------\n");
	printf(" Read Dimension   : %.3f\n Bandwidth        : %.3f Hz\n\n Acquisition Time : %.3f ms\n Rise Time        : %.3f ms\n\n Read Grad        : %.3f %%\n Read Grad Dur    : %.3f ms\n Read Momentum    : %.3f\n\n Rewind Grad      : %.3f %%\n Rewind Grad Dur  : %.3f ms\n", readDim,BP,aqq,rt,grad,gradDur,readMm,RRG,RRGDur);
	
	newES = rt - 0.000006 + 0.00002 + aqq + rt + 0.02 + RRGDur;
	printf(" Min EchoSpac. 	  : %.3f ms\n",newES);
	printf(" ---------------------------------\n\n");


	Appel_Recursif(readMm, durLim, gradLim, factor, slew, rewMm, &out.rewGrad, &out.rewDur, 1);
	
	RRG = out.rewGrad;
	RRGDur = out.rewDur;
	rewMm = slew*(RRG/100)*RRGDur;

	// Printf final values
	printf("\n ---------------------------------\n");
	printf(" Read Grad            : %.3f %%\n Read Grad Dur        : %.3f ms\n Read Momentum        : %.3f\n\n Rewind Grad          : %.3f %%\n Rewind Grad Dur      : %.3f ms (min = %.3f ms)\n Rewind Grad Momentum : %.3f\n", grad, gradDur, readMm,RRG,RRGDur, durLim, rewMm);
	
	newES = rt - 0.000006 + 0.00002 + aqq + rt + 0.02 + RRGDur;
	printf(" Min EchoSpac. 	      : %.3f ms\n",newES);
	printf(" ---------------------------------\n\n");

	// Update the EchoSpacing and get the first echo
	TE = GetFirstEcho(TE, 2.0);
	UpdateEchoSpacing();
	
	printf("\n ---------------------------------\n");
	printf(" First Echo Time : %.3f ms\n EchoSpacing : %.3f ms\n",TE, EchoSpacing);
	printf(" ---------------------------------\n\n");
	

	return 0;
}

void Appel_Recursif(double readMm, double durLim, double gradLim, double factor, double slew, double rewMm, double* rewGrad, double* rewDur, int i)
{
	printf("\n ---------- Tour = %d ----------\n", i);
	rewMm = slew*((*rewGrad)/100)*(*rewDur);

	if(rewMm == readMm){
		printf(" rewMm == readMm -> OK\n");
		if((*rewGrad)*factor <= gradLim && (*rewDur)/factor >= durLim){
				(*rewGrad) = (*rewGrad)*factor;
				(*rewDur) = (*rewDur)/factor;
				printf(" rewGrad*factor <= gradLim &&  (*rewDur)/factor >= durLim -> OK\n");

				i++;
				Appel_Recursif(readMm,durLim,gradLim,factor,slew,rewMm,rewGrad, rewDur, i);
				
		}
		else{
			printf(" rewGrad*factor <= gradLim &&  (*rewDur)/factor >= durLim -> NOK\n");

			if(gradLim/(*rewGrad) > 1 && (*rewDur)/durLim > 1){
				
				if(gradLim/(*rewGrad) > (*rewDur)/durLim){
					
					factor = (*rewDur)/durLim;
					printf(" New factor from Dur : %f\n", factor);
					i++;
					Appel_Recursif(readMm,durLim,gradLim,factor,slew,rewMm,rewGrad, rewDur, i);
				}
				
				else{
					
					factor = gradLim/(*rewGrad);
					printf(" New factor from Grad : %f\n", factor);
					i++;
					Appel_Recursif(readMm,durLim,gradLim,factor,slew,rewMm,rewGrad, rewDur, i);
				}

			}

			else
				printf(" New factor is lower than 1\n");
		}
	}
	else{
		printf(" rewMm == readMm -> NOK\n");
	}
}

double GetFirstEcho(double tmp, double minTE){
    
    double rest;
  
        // Get the closest Dixon first echo corresponding to user's desired TE
        rest = fmod(tmp,StartPt);
        
        while(rest != 0){
            if(fmod(tmp - rest,StartPt) == 0){
                tmp = tmp - rest;
                break;
            }
            else if(fmod(tmp + rest,StartPt) == 0){
                tmp = tmp + rest;
                break;
            }
            else{
                tmp = tmp - rest;
                rest=fmod(tmp,StartPt);
            }
        }
        
        while(tmp < minTE){
                    tmp += StartPt;
                }
        
        return tmp;
}

void UpdateEchoSpacing(void)
{
	//printf(" -->UpdateEchoSpacing\n\n");
	
	double minEchoSpacing, riseTime, denab, tmp, rounding;
	
	// Initialization
	riseTime = 0.268;
	denab = riseTime - 0.000006;
	minEchoSpacing = 
                        denab                   + 
			0.01			+ 
			PVM_AcquisitionTime 	+ 
                        riseTime		+
			0.02			+
			RRGDur;
        
        // Round minEchoSpacing value to .2f
        rounding = (int)(minEchoSpacing * 100 + .5);
        minEchoSpacing = (double)rounding / 100;
	
        // Re-initialize EchoSpacing before finding the minimum
        EchoSpacing = 0;
        tmp = 0;
        
	/* Calculate EchoSpacing */
	while( EchoSpacing < minEchoSpacing )
	{
		EchoSpacing += (StartPt / 2);  
        } 
        
        // fmodulus() is another way to calculate fmod() which is not working for unknown reasons
        while(fmodulus(EchoSpacing,StartPt) == 0)
        {
            tmp = EchoSpacing - (StartPt/2);
            if(tmp < minEchoSpacing){
                EchoSpacing += (StartPt / 2);
            }
            else{
                EchoSpacing = tmp;
            }
        }
     
	//printf("\n <--UpdateEchoSpacing\n");
}

double fmodulus(double a, double b)
{
    a = round(a*100);
    b = round(b*100);
    
    while(a > 0)
    {
        a = a - b;
        
        if(a < 0)
        {
            a = a + b;
            break;
        }
    }
	
	return a;
}












