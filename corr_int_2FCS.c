/*
 *  corr_int_2FCS.c
 *
 *    o
 *   /   Created by William Ridgeway 2007-2012.
 *  o
 *   \   
 *    o
 *   /   Copyright 2012 TSRI. All rights reserved.
 *  o
 *
 *	Part of the Triple Correlation Toolbox package of programs:
 *
 *	Ridgeway WK, Millar DP, Williamson JR, Vectorized data acquisition 
 *  and fast triple-correlation integrals for Fluorescence Triple 
 *  Correlation Spectroscopy, 2012
 *
 *  corr_int_2FCS.c
 *  FCS correlation integrals implementing Eq. 4 in the main manuscript.  
 *  Called by F3CS_2FCS.
 *
 */

#include "libF3CS.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>

//	***	Mini glossary of which arrays correspond to which data
//	datax		=	raw data, I_x(t)
//	gxbx		= 	correlation function as returned to the calling program
//	pxbx		=	temporary repository of products I(t)I(t+tau)
//	shatzelx	=	vector[n*(1+pmax)] which contains pseudo-logarithmically-binned intensity values. <<I(t+tau)>> (sumbox)
//	klausx		=	vector[pmax] which for a given time t, equals klausx[i] = sum_{i=0}^{2^i} I_{i+time} <<I(t)>> (tempa)
//					Roughly, the correlation function is the outer product of shatzelx (x) klausx.

void    *TripleCorr(void *threadarg)  {
	int i,j,k,l,jmin;	//Allocate memory
    struct triplecorr_struct *val;
    val = threadarg;
    
	float * g0b0 = val->fcsv->g0b0;
	float * g0b1 = val->fcsv->g0b1;
	float * g1b1 = val->fcsv->g1b1;
	float * g1b2 = val->fcsv->g1b2;
	float * g2b0 = val->fcsv->g2b0;
	float * g2b2 = val->fcsv->g2b2;
	float * klaus0 = val->fcsv->klaus0;
	float * klaus1 = val->fcsv->klaus1;
	float * klaus2 = val->fcsv->klaus2;
	float * shatzel0 = val->fcsv->shatzel0;
	float * shatzel1 = val->fcsv->shatzel1;
	float * shatzel2 = val->fcsv->shatzel2;
	float	previous_shatzel0[veclen];
	float	previous_shatzel1[veclen];
	float	previous_shatzel2[veclen];
	for		(i=0; i<veclen; i++)	previous_shatzel0[i]=shatzel0[i];//	Preserve the delayed intensity data that carries over from the previous run
	for		(i=0; i<veclen; i++)	previous_shatzel1[i]=shatzel1[i];
	for		(i=0; i<veclen; i++)	previous_shatzel2[i]=shatzel2[i];
	float * stemp0 = val->fcsv->stemp0;
	float * stemp1 = val->fcsv->stemp1;
	float * stemp2 = val->fcsv->stemp2;
	float * sc0 = val->fcsv->sc0;	// Should already be set to zero in the main routine
	float * sc1 = val->fcsv->sc1;
	float * sc2 = val->fcsv->sc2;
	
    
    unsigned short	*data0 =	val->array0;
    unsigned short	*data1 =	val->array1;
    unsigned short	*data2 =	val->array2;
	float						c0 = 0.0;
    float						c1 = 0.0;
    float						c2 = 0.0;
    float						p0b0[veclen];
    for(i=0; i<veclen; i++) p0b0[i] = 0.0;
	float						p0b1[veclen];
	for(i=0; i<veclen; i++) p0b1[i] = 0.0;
	float						p0b2[veclen];
	for(i=0; i<veclen; i++) p0b2[i] = 0.0;
	float						p1b0[veclen];
	for(i=0; i<veclen; i++) p1b0[i] = 0.0;
	float						p1b1[veclen];
	for(i=0; i<veclen; i++) p1b1[i] = 0.0;
	float						p1b2[veclen];
	for(i=0; i<veclen; i++) p1b2[i] = 0.0;
	float						p2b0[veclen];
	for(i=0; i<veclen; i++) p2b0[i] = 0.0;
	float						p2b1[veclen];
	for(i=0; i<veclen; i++) p2b1[i] = 0.0;
	float						p2b2[veclen];
	for(i=0; i<veclen; i++) p2b2[i] = 0.0;
    	
	float						deltaklaus0;
	float						deltaklaus1;
	float						deltaklaus2;
	float						deltashatzel0[veclen];
    for(i=0; i<veclen; i++) deltashatzel0[i] = 0.0;
	float						deltashatzel1[veclen];
    for(i=0; i<veclen; i++) deltashatzel1[i] = 0.0;
	float						deltashatzel2[veclen];
    for(i=0; i<veclen; i++) deltashatzel2[i] = 0.0;

	//	*****	Calculate <I(t+\tau)> = (1/T) \Sum_{i=1}^{T} [I(t_i +\tau)]	*****	//
    for(i=0; i < algorquantum; i++) {
        
        k = i;
        stemp0[0] = shatzel0[2*n-1];
        stemp1[0] = shatzel1[2*n-1];
        stemp2[0] = shatzel2[2*n-1];
        for(l=2*n-1; l>0; l--)  shatzel0[l] = shatzel0[l-1];
        for(l=2*n-1; l>0; l--)  shatzel1[l] = shatzel1[l-1];
        for(l=2*n-1; l>0; l--)  shatzel2[l] = shatzel2[l-1];
        shatzel0[0] = (float) data0[i];
        c0 		   += (float) data0[i];
        shatzel1[0] = (float) data1[i];
        c1 		   += (float) data1[i];        
        shatzel2[0] = (float) data2[i];
        c2 		   += (float) data2[i];
		
        for(j=0; j<2*n; j++)   {
			sc0[j] += shatzel0[j];
			sc1[j] += shatzel1[j];
			sc2[j] += shatzel2[j];
        }
        
        for(j=1; j<pmax; j++)  {
            if((k&1) == 1)  {
                jmin = j - 1;
                
                for(l=n; l<2*n; l++)	sc0[j*n+l] += shatzel0[j*n+l];
                for(l=n; l<2*n; l++)	sc1[j*n+l] += shatzel1[j*n+l];
                for(l=n; l<2*n; l++)	sc2[j*n+l] += shatzel2[j*n+l];
                
                stemp0[j] = shatzel0[n*(j+2)-1];
                stemp1[j] = shatzel1[n*(j+2)-1];
                stemp2[j] = shatzel2[n*(j+2)-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel0[l] = shatzel0[l-1];
                shatzel0[j*n+n] = shatzel0[j*n+n-1] + stemp0[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel1[l] = shatzel1[l-1];
                shatzel1[j*n+n] = shatzel1[j*n+n-1] + stemp1[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel2[l] = shatzel2[l-1];
                shatzel2[j*n+n] = shatzel2[j*n+n-1] + stemp2[j-1];
            }
		else    j = pmax;
		k = k >> 1;
        } 
    }
	
	klaus0[0] = 0.0;
	klaus1[0] = 0.0;
	klaus2[0] = 0.0;
		
	const float					c0avg = c0/(float)algorquantum;
	const float					c1avg = c1/(float)algorquantum;
	const float					c2avg = c2/(float)algorquantum;
	float						sc0avg[veclen];
	float						sc1avg[veclen];
	float						sc2avg[veclen];
	float						scaletemp = 1.0;
    for(i=0; i<2*n; i++)  	{
		sc0avg[i] = scaletemp*sc0[i]/(float)algorquantum;
		sc1avg[i] = scaletemp*sc1[i]/(float)algorquantum;
		sc2avg[i] = scaletemp*sc2[i]/(float)algorquantum;
	}
	for(i=2*n; i<veclen; i++)  	{
		if (i%n == 0)	scaletemp *= 2.0;
		sc0avg[i] = scaletemp*sc0[i]/(float)algorquantum;
		sc1avg[i] = scaletemp*sc1[i]/(float)algorquantum;
		sc2avg[i] = scaletemp*sc2[i]/(float)algorquantum;
	} 	

	//	Restore the delayed intensity data carried over from the previous run
	for(i=0; i<veclen; i++)	shatzel0[i] = previous_shatzel0[i];
	for(i=0; i<veclen; i++)	shatzel1[i] = previous_shatzel1[i];
	for(i=0; i<veclen; i++)	shatzel2[i] = previous_shatzel2[i];

	//	Load deltashatzel values before the fractal integral starts
	for(j=0; j<veclen; j++)   {
		deltashatzel0[j] = shatzel0[j] - sc0avg[j];	// \deltaI(t+\tau) = I(t+\tau) - \overbar{I+\tau}
		deltashatzel1[j] = shatzel1[j] - sc1avg[j];
		deltashatzel2[j] = shatzel2[j] - sc2avg[j];
	}
	
	//	*****	Calculate <I(t)I(t+\tau)> (T) = \Sum_{i=1}^{T} [I(t_i)I(t_i +\tau)]	*****	//
    for(i=0; i < algorquantum; i++) {
        
        k = i;
        stemp0[0] = shatzel0[2*n-1];
        stemp1[0] = shatzel1[2*n-1];
        stemp2[0] = shatzel2[2*n-1];
        for(l=2*n-1; l>0; l--)  shatzel0[l] = shatzel0[l-1];
        for(l=2*n-1; l>0; l--)  shatzel1[l] = shatzel1[l-1];
        for(l=2*n-1; l>0; l--)  shatzel2[l] = shatzel2[l-1];
        shatzel0[0] = (float) data0[i];
        klaus0[0]  += (float) data0[i];
        shatzel1[0] = (float) data1[i];
        klaus1[0]  += (float) data1[i];       
        shatzel2[0] = (float) data2[i];
        klaus2[0]  += (float) data2[i];
		for(j=0; j<2*n; j++)   {
			deltashatzel0[j] = shatzel0[j] - sc0avg[j];
			deltashatzel1[j] = shatzel1[j] - sc1avg[j];
			deltashatzel2[j] = shatzel2[j] - sc2avg[j];
		}

		deltaklaus0 = data0[i]-c0avg;	// \deltaI(t) = I(t) - \overbar{I}
		deltaklaus1 = data1[i]-c1avg;
		deltaklaus2 = data2[i]-c2avg;
		
        for(j=0; j<2*n; j++)   {
        	p0b0[j] = p0b0[j] + (float) deltaklaus0*deltashatzel0[j];
        	p0b1[j] = p0b1[j] + (float) deltaklaus0*deltashatzel1[j];
        	p0b2[j] = p0b2[j] + (float) deltaklaus0*deltashatzel2[j];
        }
		for(j=0; j<2*n; j++)   {
			p1b0[j] = p1b0[j] + (float) deltaklaus1*deltashatzel0[j];
        	p1b1[j] = p1b1[j] + (float) deltaklaus1*deltashatzel1[j];
        	p1b2[j] = p1b2[j] + (float) deltaklaus1*deltashatzel2[j];
        }
		for(j=0; j<2*n; j++)   {
			p2b0[j] = p2b0[j] + (float) deltaklaus2*deltashatzel0[j];
        	p2b1[j] = p2b1[j] + (float) deltaklaus2*deltashatzel1[j];
        	p2b2[j] = p2b2[j] + (float) deltaklaus2*deltashatzel2[j];
        } 
		
        for(j=1; j<pmax; j++)  {
            if((k&1) == 1)  {
                jmin = j - 1;

				deltaklaus0 = klaus0[jmin] - (c0avg*pow(2.0,(float)j));
				deltaklaus1 = klaus1[jmin] - (c1avg*pow(2.0,(float)j));
				deltaklaus2 = klaus2[jmin] - (c2avg*pow(2.0,(float)j));

                for(l=(j+1)*n; l<(j+2)*n; l++)	p0b0[l] = p0b0[l] + deltaklaus0*deltashatzel0[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p0b1[l] = p0b1[l] + deltaklaus0*deltashatzel1[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p0b2[l] = p0b2[l] + deltaklaus0*deltashatzel2[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p1b0[l] = p1b0[l] + deltaklaus1*deltashatzel0[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p1b1[l] = p1b1[l] + deltaklaus1*deltashatzel1[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p1b2[l] = p1b2[l] + deltaklaus1*deltashatzel2[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p2b0[l] = p2b0[l] + deltaklaus2*deltashatzel0[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p2b1[l] = p2b1[l] + deltaklaus2*deltashatzel1[l];
                for(l=(j+1)*n; l<(j+2)*n; l++)	p2b2[l] = p2b2[l] + deltaklaus2*deltashatzel2[l];
				
                klaus0[j] += klaus0[j-1];
                klaus0[j-1] = 0.0;
                klaus1[j] += klaus1[j-1];
                klaus1[j-1] = 0.0;
                klaus2[j] += klaus2[j-1];
                klaus2[j-1] = 0.0;
                stemp0[j] = shatzel0[n*(j+2)-1];
                stemp1[j] = shatzel1[n*(j+2)-1];
                stemp2[j] = shatzel2[n*(j+2)-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel0[l] = shatzel0[l-1];
                shatzel0[j*n+n] = shatzel0[j*n+n-1] + stemp0[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel1[l] = shatzel1[l-1];
                shatzel1[j*n+n] = shatzel1[j*n+n-1] + stemp1[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel2[l] = shatzel2[l-1];
                shatzel2[j*n+n] = shatzel2[j*n+n-1] + stemp2[j-1];
				
				for(l=(j+1)*n; l<(j+2)*n; l++)	{
					deltashatzel0[l] = shatzel0[l] - sc0avg[l];
					deltashatzel1[l] = shatzel1[l] - sc1avg[l];
					deltashatzel2[l] = shatzel2[l] - sc2avg[l];
				}
            }
		else    j = pmax;
		k = k >> 1;
        } 
    }
	val->fcsv->tc0 = val->fcsv->tc0 + c0;
    val->fcsv->tc0sq = val->fcsv->tc0sq + pow(c0,2.0f);
    val->fcsv->tc1 = val->fcsv->tc1 + c1;
    val->fcsv->tc1sq = val->fcsv->tc1sq + pow(c1,2.0f);
    val->fcsv->tc2 = val->fcsv->tc2 + c2;
    val->fcsv->tc2sq = val->fcsv->tc2sq + pow(c2,2.0f);
	
    for(i=0; i<veclen; i++)    g0b0[i] = g0b0[i] + (float)p0b0[i];
    for(i=0; i<veclen; i++)    g1b1[i] = g1b1[i] + (float)p1b1[i];
    for(i=0; i<veclen; i++)    g2b2[i] = g2b2[i] + (float)p2b2[i];
    
    for(i=0; i<veclen; i++)    g0b1[i] = g0b1[i] + 0.5f*((float)p0b1[i] + (float)p1b0[i]);
    for(i=0; i<veclen; i++)    g1b2[i] = g1b2[i] + 0.5f*((float)p1b2[i] + (float)p2b1[i]);
    for(i=0; i<veclen; i++)    g2b0[i] = g2b0[i] + 0.5f*((float)p2b0[i] + (float)p0b2[i]);
		
    pthread_exit(NULL);
}
