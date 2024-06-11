/*
 *  GlobalFit.c
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
 *  GlobalFit.c
 *  Self-contained code to generate F3CS_GlobalFit.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "libF3CS.h"

#undef veclen
#define	veclen n_integral*(1+pmax)
int	n_integral = n;
#undef  n

//	The definition of useplplot determines whether the program is complied as F3CS_GlobalFit or F3CS_GlobalFit_3D_Plot with the plplot libraries
// #define useplplot
#ifdef  useplplot
#include <plplot/plplot.h>
#endif

#define LEVELS  10				// Contour plot levels
#define num_curves_to_consider 12800
#define	fcs_model_params	19
#define nptw	41
#define	N	40

#ifdef useplplot
PLDLLIMPEXP void
plexit( const char *errormsg );
#endif


struct m_pass_struct    {
    float  N_t;
    float  theta;
    float  tauD_1;
    float  tauD_2;
    float  omega;
    float  g_inf;
    float  tauf_1;
    float  tauf_2;
    float  T0_1;
    float  T0_2;
    float  gamma_1;
    float  gamma_2;
    float  frac_gamma_1;
    float  frac_gamma_2;
    float  k2;
    float  k3;
    float  k4;
    float  int0;
	float  gamma_R;
    
    double * m;
	double * m_no_triplet;
	double * mGGG;
};

struct g_pass_struct    {
    double  NN[10];
    double  tauD_i[10];
    double  omega;
    double  g_infIxJ[7];
    double  tauf_i[10];
    double  T0_i[10];
    double  gamma_i[10];
    double  frac_gamma_i[10];
    double	q_i[3];
    double	eta[9];
    double  int0;
    double  int1;
    double  int2;
    double	isd0;
    double	isd1;
    double	isd2;
    double	intensityweight;
	double  gamma_R;
    
    double	eps[30];
    double	occ[30];
    double	bg[3];
    double	cc[3];
    double	k2[3];
    double	k3[3];
    double	k4[3];
	double  Veff[6];
    
    int	 p;
    int  fix_N_t;
    int  fix_theta;
    int  fix_tauD_1;
    int  fix_tauD_2;
    int  fix_omega;
    int  fix_g_inf;
    int  fix_tauf_1;
    int  fix_tauf_2;
    int  fix_T0_1;
    int  fix_T0_2;
    int  fix_gamma_1;
    int  fix_gamma_2;
    int  fix_frac_gamma_1;
    int  fix_frac_gamma_2;
    int  fit_start_triplecorr;
    int  fit_stop_triplecorr;
    int  fit_start;
    int  fit_stop;
    int	 Gfix[9];
    int	 Ion[3];
    int	 verbose;

    int	   * 	pmask;
    double * 	m;
	double *    m_no_triplet;
	double *	mGGG;
	double *	y; 
    double *	sigma;
    double *	fits;
    double * 	w_resids;
};

// ***	Calculate a pseudo-logarithmic series of tau values	***	//
double  tau_G( int i){
	if (i<0)		return 0.0;
	
    double t = 1.0 / 1250000.0;
    int	n_tau = n_integral;
    
    int j = (i-n_tau-1)/n_tau;
    if (i<2*n_tau)   	return t * i;

    double tempa = pow(2.0, j);
    
    t = t*(tempa*(i-j*n_tau)+1);
    
    return t;
}

//	***	Calculate a tau series for triple-correlation	**	//
double  tau_GGG( int i){

	if (i<0)		return 0.0;
	
    double t = 1.0 / 1250000.0;
    int	n_tau = n_GGG;
    
    int j = (i-n_tau-1)/n_tau;
    if (i<2*n_tau)   	return t * i;

    double tempa = pow(2.0, j);
    
    t = t*(tempa*(i-j*n_tau)+1);	// t = 1/freq * 2^((i-n-1)/n) * (i-((i-n-1)/n)*n - 1)
    
    return t;

}

//	****	Calculates double-correlation fit functions	****	//
void    m_tau(void * m_pass)    {
    struct m_pass_struct * a;
    a = m_pass;
    int i;

    double t = 0.0;
    double  tempa1 = 0.0;
    double  tempa2 = 0.0;
    double  tempb1 = 0.0;
    double  tempb2 = 0.0;
    double  tempc1 = 0.0;
    double  tempc2 = 0.0;
    double  tempd1 = 0.0;
    double  tempd2 = 0.0;
    double	tempg1 = 0.0;
    double	tempg2 = 0.0;
    double	tempe1 = 0.0;
    double	tempe2 = 0.0;

    for(i=0; i < veclen; i++)  {
        t = tau_G(i) * 1000000.0;
        tempa1 = 1.0 / (1.0 + t / a->tauD_1);
        tempa2 = 1.0 / (1.0 + t / a->tauD_2);
        tempb1 = pow(1.0 + t / (a->omega*a->omega*a->tauD_1), -0.5);
        tempb2 = pow(1.0 + t / (a->omega*a->omega*a->tauD_2), -0.5);
        tempc1 = 1.0 / (1.0 + t / (a->gamma_1*a->tauD_1));
        tempc2 = 1.0 / (1.0 + t / (a->gamma_2*a->tauD_2));
        tempd1 = pow(1.0 + t / (a->gamma_1*a->omega*a->omega*a->tauD_1), -0.5);
        tempd2 = pow(1.0 + t / (a->gamma_2*a->omega*a->omega*a->tauD_2), -0.5);
        tempg1 = a->T0_1*a->theta / (1.0 - a->T0_1*a->theta - a->T0_2*(1.0-a->theta));
        tempg2 = a->T0_2*(1.0-a->theta) / (1.0 - a->T0_1*a->theta - a->T0_2*(1.0-a->theta));
		tempe1 = -1.0*t/a->tauf_1;
        tempe1 = exp(tempe1);
        tempe2 = -1.0*t/a->tauf_2;
        tempe2 = exp(tempe2);
		a->m_no_triplet[i] = (1.0/a->N_t)*(a->theta*((1.0-a->frac_gamma_1)*tempa1*tempb1 + a->frac_gamma_1*tempc1*tempd1)
			+ (1.0-a->theta)*((1.0-a->frac_gamma_2)*tempa2*tempb2 + a->frac_gamma_2*tempc2*tempd2));
		a->m[i] = a->m_no_triplet[i] * (1+tempg1*tempe1+tempg2*tempe2) + a->g_inf;
		a->m_no_triplet[i] =a->m_no_triplet[i] + a->g_inf;
    }  

    return;
}

//	****	Calculates triple-correlation fit functions	****	//
void    m_tauGGG(void * m_pass)    {
    struct m_pass_struct * a;
    a = m_pass;
    int i,j;

    double t = 0.0;
    double  tempa1 = 0.0;
    double  tempa2 = 0.0;
    double  tempb1 = 0.0;
    double  tempb2 = 0.0;
    double  tempc1 = 0.0;
    double  tempc2 = 0.0;
    double  tempd1 = 0.0;
    double  tempd2 = 0.0;

    for(i=0; i < veclen_GGG; i++)  {	//	G(tau2-tau1) ==> m[tau1*veclen+tau2]
		for(j=i; j < veclen_GGG; j++)  {
			t = fabs(tau_GGG(i)-tau_GGG(j)) * 1000000.0;
			tempa1 = 1.0 / (1.0 + t / a->tauD_1);
			tempa2 = 1.0 / (1.0 + t / a->tauD_2);
			tempb1 = pow(1.0 + t / (a->omega*a->omega*a->tauD_1), -0.5);  // (1/(1+tau/(w^2 tau_D))^1/2
			tempb2 = pow(1.0 + t / (a->omega*a->omega*a->tauD_2), -0.5);
			tempc1 = 1.0 / (1.0 + t / (a->gamma_1*a->tauD_1));
			tempc2 = 1.0 / (1.0 + t / (a->gamma_2*a->tauD_2));
			tempd1 = pow(1.0 + t / (a->gamma_1*a->omega*a->omega*a->tauD_1), -0.5);
			tempd2 = pow(1.0 + t / (a->gamma_2*a->omega*a->omega*a->tauD_2), -0.5);
			a->mGGG[j*veclen_GGG+i] = (1.0/1.0)*(a->theta*((1.0-a->frac_gamma_1)*tempa1*tempb1 + a->frac_gamma_1*tempc1*tempd1)
				+ (1.0-a->theta)*((1.0-a->frac_gamma_2)*tempa2*tempb2 + a->frac_gamma_2*tempc2*tempd2));
		}
    }  
	for(i=1; i < veclen_GGG; i++)  {	// Fill in lower diag. of matrix; necessarily symmetric.
		for(j=0; j < i; j++)  {
			a->mGGG[j*veclen_GGG+i] = a->mGGG[i*veclen_GGG+j];	
		}
	}

    return;
}

//	****	Calculates triple-correlation fit functions	****	//
void    m_tau_tripleG(void * m_pass)    {
    struct m_pass_struct * a;
    a = m_pass;
    int i,j;

    double t1 = 0.0;
	double t2 = 0.0;
    double  tempa1 = 0.0;
    double  tempb1 = 0.0;
    double  tempc1 = 0.0;
    double  tempd1 = 0.0;

	double	td1	   = 1.0/(3.0*a->tauD_1*a->tauD_1);
	double	td1b   = 1.0/(3.0*a->tauD_1);
	double	td2    = 1.0/(3.0*a->omega*a->omega*a->omega*a->omega*a->tauD_1*a->tauD_1);
	double	td2b   = 1.0/(3.0*a->omega*a->omega*a->tauD_1);

	double	td3	   = 1.0/(3.0*a->gamma_1*a->gamma_1*a->tauD_1*a->tauD_1);
	double	td3b   = 1.0/(3.0*a->gamma_1*a->tauD_1);
	double	td4    = 1.0/(3.0*a->gamma_1*a->gamma_1*a->omega*a->omega*a->omega*a->omega*a->tauD_1*a->tauD_1);
	double	td4b   = 1.0/(3.0*a->gamma_1*a->omega*a->omega*a->tauD_1);
	double	prefactor = pow(4.0/3.0,3.0);	//	gamma_3/ gamma_2^2 ratio
	prefactor = a->gamma_R;		// Expt. optimized factor; comment out to keep gaussian result from previous line

    for(i=0; i < veclen_GGG; i++)  {	//	GGG(tau1,tau2) ==> m[tau1*veclen+tau2]
		for(j=i; j < veclen_GGG; j++)  {
			t1 = tau_GGG(i) * 1000000.0;
			t2 = tau_GGG(j) * 1000000.0;
			tempa1 = 1.0 / (1.0 + 4.0*t1*(t2-t1)*td1 + 4.0*t2*td1b);
			tempb1 = pow((1.0 + 4.0*t1*(t2-t1)*td2 + 4.0*t2*td2b), -0.5);
			tempc1 = 1.0 / (1.0 + 4.0*t1*(t2-t1)*td3 + 4.0*t2*td3b);
			tempd1 = pow((1.0 + 4.0*t1*(t2-t1)*td4 + 4.0*t2*td4b), -0.5);
			a->mGGG[j*veclen_GGG+i] = prefactor*((1.0-a->frac_gamma_1)*tempa1*tempb1 + a->frac_gamma_1*tempc1*tempd1);
		}
    }  
	for(i=1; i < veclen_GGG; i++)  {	// Fill in lower diag. of matrix.
		for(j=0; j < i; j++)  {
			a->mGGG[j*veclen_GGG+i] = a->mGGG[i*veclen_GGG+j];	
		}
	}

    return;
}


void    grp_m_tau(const gsl_vector * x, void * g_pass, gsl_vector * f)    {
    struct g_pass_struct * a;
    a = g_pass;
    struct m_pass_struct o;
    int i,j,k,w,p;
 
	//Set terms to trivial values to keep compatibility w/ gtauglobal m_tau routine
	o.N_t = 1.0;
	o.theta = 1.0;	
	o.tauD_2 = 1.0;
	o.g_inf = 0.0;
	o.tauf_2 = 1.0;
	o.T0_2 = 0.0;
	o.gamma_2 = 1.0;			
	o.frac_gamma_2 = 0.0;
	o.m = a->m;
	o.m_no_triplet = a->m_no_triplet;
	o.mGGG = a->mGGG;

    size_t	apd0[6]={0,1,2,0,1,2};
    size_t	apd1[6]={0,1,2,1,2,0};
    double tempa = 0.0;
	
	int memct=0;	
	memct = 6*veclen + 3*veclen_GGG*veclen_GGG;
	
	for(i=0; i<memct; i++)	a->fits[i] = 0.0;
	
	double	NN[10];
	double	tauD_i[10];
	double	omega;
	double  gamma_R;
	double	g_infIxJ[7];
	double	tauf_i[10];
	double	T0_i[10];
	double	gamma_i[10];
	double	frac_gamma_i[10];
	double	q_i[3];
	double	bg[3];
	double	cc[3];
	double	tempt  = 0.0;
	double	intensity[3];
	double	inv_intensity[6];


	j=0;
	for(i=0; i<10; i++)	NN[i] = a->NN[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			NN[a->pmask[i]-j] = gsl_vector_get(x, i);

//	Code Below fits to N / Theta, not N0,N1.  Comment out to change behaviour.	(currently out)
	tempt = NN[0];
//	NN[0] = tempt * NN[1];
//	NN[1] = tempt * (1.0-NN[1]);
	tempt = 0.0;

	j+= 10;
	for(i=0; i<10; i++)	tauD_i[i] = a->tauD_i[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			tauD_i[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	omega = a->omega;
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			omega = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<7; i++)	g_infIxJ[i] = a->g_infIxJ[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			g_infIxJ[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<10; i++)	tauf_i[i] = a->tauf_i[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			tauf_i[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<10; i++)	T0_i[i] = a->T0_i[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			T0_i[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<10; i++)	gamma_i[i] = a->gamma_i[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			gamma_i[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<10; i++)	frac_gamma_i[i] = a->frac_gamma_i[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			frac_gamma_i[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<3; i++)	q_i[i] = a->q_i[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			q_i[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<3; i++)	bg[i] = a->bg[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			bg[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	for(i=0; i<3; i++)	cc[i] = a->cc[i];
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			cc[a->pmask[i]-j] = gsl_vector_get(x, i);
	j+= 10;
	gamma_R = a->gamma_R;
	for(i=0; i<a->p; i++)
		if(a->pmask[i] > (j-1) && a->pmask[i] < (j+10))
			gamma_R = gsl_vector_get(x, i);
		
			

	double	epsilon[30];	//	The Brightness Matrix (Supp. Materials, 30S ribosomal subunit paper)
	for(i=0; i<10; i++)	{	// Species
 		for(j=0; j<3; j++)	{	// channel
			epsilon[10*j+i] = 0.0;
 			for(k=0; k<3; k++)	{	// dye
				epsilon[10*j+i] += 1000.0*a->eta[j*3+k]*cc[j]*q_i[k]*a->occ[10*k+i];
 	}	}	}
	 			
	for(j=0; j<3; j++)	{
		intensity[j] = bg[j]*1000.0;
		for(i=0; i<10; i++)		intensity[j] += NN[i]*epsilon[10*j+i];
	}
	
				
	for(j=0; j<10; j++)	{	//For each species, generate a correlation curve, then add it to the nascent correlation functions with appropriate weighting based on brightness
		if (NN[j] != 0)	{
			i=0;
			o.tauD_1 = tauD_i[j];
			o.omega = omega;
			o.tauf_1 = tauf_i[j];
			o.T0_1 = T0_i[j];
			o.gamma_1 = gamma_i[j];
			o.frac_gamma_1 = frac_gamma_i[j];
			o.gamma_R = gamma_R;	
			
			m_tau(&o);
			for(k=0; k<3; k++)	{
				tempa = epsilon[apd0[k]*10+j]*epsilon[apd1[k]*10+j]*NN[j]/a->Veff[k];
				for(i=0; i<veclen; i++)	a->fits[k*veclen+i] += tempa*o.m[i];
			}
			for(k=3; k<6; k++)	{
				tempa = epsilon[apd0[k]*10+j]*epsilon[apd1[k]*10+j]*NN[j]/a->Veff[k];
				for(i=0; i<veclen; i++)	a->fits[k*veclen+i] += tempa*o.m_no_triplet[i];
			}
						
	//	***	Triple Correlation.  oh yeah...	***	//
			m_tau_tripleG(&o);
			for(i=0; i<veclen_GGG; i++)	{	//tau1
				for(w=0; w<veclen_GGG; w++)	{	//tau2		
					for(p=0; p<3; p++)	{	//	{G0x1x2, G1x2x0, G2x0x1}
						tempa = epsilon[0*10+j]*epsilon[1*10+j]*epsilon[2*10+j]*NN[j];
						a->fits[6*veclen+(p*veclen_GGG*veclen_GGG+i*veclen_GGG+w)] += tempa*o.mGGG[i*veclen_GGG+w];	//	G(tau1,tau2) ==> m[tau1*veclen+tau2]
					}
				}
			}
	}    }	//End of Species loop
    
	
	double	triple_int_product = 1;
	for(j=0; j<3; j++)	triple_int_product = triple_int_product * intensity[j];	// Denominator for Triple correlation

	inv_intensity[0] = intensity[0]*intensity[0];
	inv_intensity[1] = intensity[1]*intensity[1];
	inv_intensity[2] = intensity[2]*intensity[2];
	inv_intensity[3] = intensity[0]*intensity[1];
	inv_intensity[4] = intensity[1]*intensity[2];
	inv_intensity[5] = intensity[2]*intensity[0];

	for(i=0; i<6; i++)	inv_intensity[i] = 1.0 / inv_intensity[i];

	for(i=0; i<veclen; i++)	{
		for(j=0; j<6; j++)	{
			a->fits[(i+j*veclen)] = a->fits[(i+j*veclen)] * inv_intensity[j]+g_infIxJ[j];
	}	}
	
	for(i=0; i<veclen_GGG; i++)	{
		for(w=0; w<veclen_GGG; w++)	{
			for(j=0; j<3; j++)	{
				a->fits[6*veclen+(j*veclen_GGG*veclen_GGG+i*veclen_GGG+w)] = (a->fits[6*veclen+(j*veclen_GGG*veclen_GGG+i*veclen_GGG+w)] / triple_int_product) + g_infIxJ[6];
	}	}	}
     
//Afterpulsing Corrections
	double	t;
    double	tempk_0  = 0.0;
    double	tempk2_0 = a->k2[0]/a->int0;
    double	tempk3_0 = a->k3[0]/a->int0;
    double	tempk4_0 = a->k4[0]/a->int0;
    double	tempk_1  = 0.0;
    double	tempk2_1 = a->k2[1]/a->int1;
    double	tempk3_1 = a->k3[1]/a->int1;
    double	tempk4_1 = a->k4[1]/a->int1;
	double	tempk_2  = 0.0;
    double	tempk2_2 = a->k2[2]/a->int2;
    double	tempk3_2 = a->k3[2]/a->int2;
    double	tempk4_2 = a->k4[2]/a->int2;
	
	i=1;	//Ignore i=0, as 1/0 = inf
	do	{
		t = tau_G(i) * 1000000.0;
		tempt = 1.0/tau_G(i);
		tempt = tempt / tau_G(i);
		tempk_0 =  tempk2_0 * tempt;
		tempk_1 =  tempk2_1 * tempt;
		tempk_2 =  tempk2_2 * tempt;
		tempt = tempt / tau_G(i);
		tempk_0 += tempk3_0 * tempt;
		tempk_1 += tempk3_1 * tempt;
		tempk_2 += tempk3_2 * tempt;
		tempt = tempt / tau_G(i);
		tempk_0 += tempk4_0 * tempt;
		tempk_1 += tempk4_1 * tempt;
		tempk_2 += tempk4_2 * tempt;

//	Comment out next three lines to disable Afterpulsing corrections
		a->fits[0*veclen+i] += tempk_0;
		a->fits[1*veclen+i] += tempk_1;
		a->fits[2*veclen+i] += tempk_2;
		i++;
	}	while (t < 100);
          
	for(i=0; i<veclen; i++)	{
		for(j=0; j<6; j++)	{
			tempa = a->fits[(i+j*veclen)] - a->y[(i+j*veclen)];
			a->w_resids[(i+j*veclen)] = tempa / a->sigma[(i+j*veclen)];
		}
	} 

	for(i=0; i<veclen_GGG; i++)	{
		for(w=0; w<veclen_GGG; w++)	{
			for(j=0; j<3; j++)	{
				tempa = a->fits[6*veclen+(j*veclen_GGG*veclen_GGG+i*veclen_GGG+w)] - a->y[6*veclen+(j*veclen_GGG*veclen_GGG+i*veclen_GGG+w)];
				a->w_resids[6*veclen+(j*veclen_GGG*veclen_GGG+i*veclen_GGG+w)] = tempa / a->sigma[6*veclen+(j*veclen_GGG*veclen_GGG+i*veclen_GGG+w)];
	}	}	}

	k=0;

	for(j=0; j<6; j++)	{
		if(a->Gfix[j]==1)	{
			for(i=0; i<(a->fit_stop-a->fit_start); i++)	{
				gsl_vector_set(f,i+k,a->w_resids[(i+(a->fit_start)+j*veclen)]);
			}
			k+= a->fit_stop-a->fit_start;
		}
	}

	int tempp = 0;
	tempp = 6 * veclen;	
	
	for(j=0; j<3; j++)	{
		if(a->Gfix[j+6]==1)	{
			for(i=0; i<(a->fit_stop_triplecorr-a->fit_start_triplecorr); i++)	{
				for(w=0; w<(a->fit_stop_triplecorr-a->fit_start_triplecorr); w++)	{
					gsl_vector_set(f,k,a->w_resids[tempp+(j*veclen_GGG*veclen_GGG+(i+a->fit_start_triplecorr)*veclen_GGG+w+a->fit_start_triplecorr)]);
					k++;
	}	}	}	}

	if (a->Ion[0] == 1)	{	gsl_vector_set(f,k, a->intensityweight*(intensity[0] - a->int0) / a->isd0);	k++;}
	if (a->Ion[1] == 1)	{	gsl_vector_set(f,k, a->intensityweight*(intensity[1] - a->int1) / a->isd1);	k++;}
	if (a->Ion[2] == 1)	{	gsl_vector_set(f,k, a->intensityweight*(intensity[2] - a->int2) / a->isd2);	k++;}
			
    return;
}

//	***	Routing functions for use with the gnuplot solver driver	***	//
int grp_tau_f (const gsl_vector * x, void * g_pass, gsl_vector * f)   {

	grp_m_tau(x, g_pass, f);

    return GSL_SUCCESS;
}

        /* Jacobian matrix J(i,j) = dfi / dxj, */
int grp_tau_df (const gsl_vector * x, void *g_pass, gsl_matrix * J)  {
	struct g_pass_struct * a = g_pass;
	int i,j,k,l,w;
	double	x_pl = 1.00001;
	double  x_mi = 0.99999;
	double	x_pl_add = 0.001;
	double	x_mi_add = -0.001;
	gsl_vector * plus  = gsl_vector_alloc (veclen*(6)+3+3*veclen_GGG*veclen_GGG);
	gsl_vector * minus = gsl_vector_alloc (veclen*(6)+3+3*veclen_GGG*veclen_GGG);
	gsl_vector * x_tweaked = gsl_vector_alloc (a->p);

	for(k=0; k<a->p; k++)	{	//k \in {parameters}
		for (i=0; i<a->p; i++)		gsl_vector_set (x_tweaked, i, gsl_vector_get(x, i));
		
		if(gsl_vector_get(x,k) != 0)	{
			gsl_vector_set (x_tweaked, k, gsl_vector_get(x, k)*x_pl);
			grp_m_tau(x_tweaked, g_pass, plus);
			
			gsl_vector_set (x_tweaked, k, gsl_vector_get(x, k)*x_mi);
			grp_m_tau(x_tweaked, g_pass, minus);

			l=0;
			for(j=0; j<6; j++)		{	//j \in {GIxJ curves}
				if(a->Gfix[j]==1)	{
					for(i=0; i<(a->fit_stop-a->fit_start); i++)		{	//i \in {taus}
						gsl_matrix_set(J,i+l,k,(gsl_vector_get(plus, i+l) - gsl_vector_get(minus,i+l))/((x_pl-x_mi)*gsl_vector_get(x,k)));
					}
					l+=a->fit_stop-a->fit_start;
				}
			}
			if (a->Ion[0] == 1)	{gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus, l))/((x_pl-x_mi)*gsl_vector_get(x,k))); l++;}
			if (a->Ion[1] == 1)	{gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus, l))/((x_pl-x_mi)*gsl_vector_get(x,k))); l++;}
			if (a->Ion[2] == 1)	{gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus, l))/((x_pl-x_mi)*gsl_vector_get(x,k))); l++;}
		
			for(i=0; i<(a->fit_stop_triplecorr-a->fit_start_triplecorr); i++)	{
				for(w=0; w<(a->fit_stop_triplecorr-a->fit_start_triplecorr); w++)	{
					for(j=0; j<3; j++)	{
						if(a->Gfix[j+6]==1)	{
							gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus,l))/((x_pl-x_mi)*gsl_vector_get(x,k)));
							l++;
			}	}	}	}
		}
		else	{	// if fit parameter == 0, do an additive perturbation,  not multiplicative
		gsl_vector_set (x_tweaked, k, gsl_vector_get(x, k)+x_pl_add);
			grp_m_tau(x_tweaked, g_pass, plus);
			
			gsl_vector_set (x_tweaked, k, gsl_vector_get(x, k)+x_mi_add);
			grp_m_tau(x_tweaked, g_pass, minus);

			l=0;
			for(j=0; j<6; j++)		{	//j \in {GIxJ curves}
				if(a->Gfix[j]==1)	{
					for(i=0; i<(a->fit_stop-a->fit_start); i++)		{	//i \in {taus}
						gsl_matrix_set(J,i+l,k,(gsl_vector_get(plus, i+l) - gsl_vector_get(minus,i+l))/(x_pl_add-x_mi_add));
					}
					l+=a->fit_stop-a->fit_start;
				}
			}
			if (a->Ion[0] == 1)	{gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus, l))/(x_pl_add-x_mi_add)); l++;}
			if (a->Ion[1] == 1)	{gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus, l))/(x_pl_add-x_mi_add)); l++;}
			if (a->Ion[2] == 1)	{gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus, l))/(x_pl_add-x_mi_add)); l++;}
		
			for(i=0; i<(a->fit_stop_triplecorr-a->fit_start_triplecorr); i++)	{
				for(w=0; w<(a->fit_stop_triplecorr-a->fit_start_triplecorr); w++)	{
					for(j=0; j<3; j++)	{
						if(a->Gfix[j+6]==1)	{
							gsl_matrix_set(J,l,k,(gsl_vector_get(plus, l) - gsl_vector_get(minus,l))/(x_pl_add-x_mi_add));
							l++;
			}	}	}	}

		}
	}
	
    return GSL_SUCCESS;
}

//	***	Routing functions for use with the gnuplot solver driver	***	//
int grp_tau_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
    grp_tau_f (x, data, f);
    grp_tau_df (x, data, J);

    return GSL_SUCCESS;
}

#ifdef useplplot
static PLOptionTable options[] = {{NULL, NULL, NULL, NULL, 0, NULL, NULL }};
#endif

//	***		Main Routine	***	//
int main(int argc, char *argv[])    {
	FILE * fp;
	int i,j,k,junk,w,i2,tau2,tau,num_selected_curves,comma_yet;
	double	G_max, G_min;
#ifdef useplplot	
	double tempa;
#endif
	struct g_pass_struct o;
	

#ifdef useplplot	
	if(argc != 5)  {
        printf("Usage: F3CS_GlobalFit_3DPlot inputfile.txt tag verbose?{0,1} fit_mode {plotdata | plotguess | fitdata | filegen}\ne.g.> F3CS_GlobalFit_3DPlot input.sA.txt sampleA 0 plotguess\n");
        return 1;
    }
#else
	if(argc != 5)  {
        printf("Usage: F3CS_GlobalFit inputfile.txt tag verbose?{0,1} fit_mode {plotdata | plotguess | fitdata | filegen}\ne.g.> F3CS_GlobalFit input.sA.txt sampleA 0 plotguess\n");
        return 1;
    }	
#endif
	
/************* Read in Fit Parameters from inputfile *************/
    
    fp = fopen(argv[1],"r");
    if (fp == NULL) {
        printf("Input file %s not found.  Exiting.\n", argv[1]);
        return 1;
    }
    int	verbose = 0;
    sscanf(argv[3],"%i", &i);
    if (i == 1)	{verbose = 1; printf("Verbose Mode\n");}
   	if (i != 1 && i != 0)   	{
        printf("Improper Verbose Flag (%s).  Print either 0 or 1.  Exiting.\n", argv[3]);
        return 1;
    }
    int	fit_mode = 0;
    char	fit_string[8] ="fitdata";
    char	plotdata_string[9] ="plotdata";
    char	plotguess_string[10] ="plotguess";
    char	filegen_string[11] ="filegen";
    if(strcmp(argv[4],fit_string) == 0){
    		fit_mode = 0; 
    		printf("Fit Mode\n");
    	}
    else if(strcmp(argv[4],plotdata_string)==0){
    		fit_mode = 1; 
    		printf("Plot Data Mode\n");
    	}
	else if(strcmp(argv[4],plotguess_string)==0){
    		fit_mode = 2; 
    		printf("Plot Guess Mode\n");
    	}
	else if(strcmp(argv[4],filegen_string)==0){
    		fit_mode = 4; 
    		printf("FileGen Mode\n");
    	}
	else	{	
		printf("Improper Fit Mode(%s).  fit_mode {plotdata | plotguess | fitdata | filegen}.  Exiting.\n",argv[4]);
		return 1;
    }
    
    //  Model Parameters
    for(i=0; i<30; i++)	o.eps[i] = -10.0;
	for(i=0; i<30; i++)	o.occ[i] = -10.0;
	for(i=0; i<3; i++)	o.q_i[i] = -10.0;
	for(i=0; i<9; i++)	o.eta[i] = -10.0;
    for(i=0; i<3; i++)	o.k2[i] = 0.0;
    for(i=0; i<3; i++)	o.k3[i] = 0.0;
    for(i=0; i<3; i++)	o.k4[i] = 0.0;
    for(i=0; i<3; i++)	o.bg[i] = 0.0;
    for(i=0; i<3; i++)	o.cc[i] = 0.0;
	for(i=0; i<6; i++)	o.Veff[i] = 0.0;
    
    int  fix_N_t[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int  fix_tauD_i[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int  fix_omega = -1;
	int	 fix_gamma_R = -1;
    int  fix_g_infIxJ[7] = {-1,-1,-1,-1,-1,-1,-1};
    int  fix_tauf_i[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int  fix_T0_i[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int  fix_gamma_i[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int  fix_frac_gamma_i[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int  fix_q_i[3]={-1,-1,-1};
    int  fix_bg[3]={-1,-1,-1};
    int	 fix_cc[3]={-1,-1,-1};
    
    char    line[512];
	char	* ignorereturnvalue;
	int linenum = 1;
    
    int ret = 0;
    ignorereturnvalue=fgets(line, 512, fp);linenum++;
    ignorereturnvalue=fgets(line, 512, fp);linenum++;
    ignorereturnvalue=fgets(line, 512, fp);linenum++;

	ret = fscanf(fp,"G0x0=%i\n", &o.Gfix[0]);linenum++;
	ret+= fscanf(fp,"G1x1=%i\n", &o.Gfix[1]);linenum++;
	ret+= fscanf(fp,"G2x2=%i\n", &o.Gfix[2]);linenum++;
	ret+= fscanf(fp,"G0x1=%i\n", &o.Gfix[3]);linenum++;
	ret+= fscanf(fp,"G1x2=%i\n", &o.Gfix[4]);linenum++;
	ret+= fscanf(fp,"G2x0=%i\n", &o.Gfix[5]);linenum++;
	ret+= fscanf(fp,"G0x1x2=%i\n", &o.Gfix[6]);linenum++;
	ret+= fscanf(fp,"G1x2x0=%i\n", &o.Gfix[7]);linenum++;
	ret+= fscanf(fp,"G2x0x1=%i\n", &o.Gfix[8]);linenum++;
    if(verbose)	printf("G0x0=%i\nG1x1=%i\nG2x2=%i\nG0x1=%i\nG1x2=%i\nG2x0=%i\n", o.Gfix[0], o.Gfix[1], o.Gfix[2], o.Gfix[3], o.Gfix[4], o.Gfix[5]);
    if (ret != 9 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    
    ignorereturnvalue=fgets(line, 512, fp);linenum++;
    ret = fscanf(fp,"Intensity_Weight=%lf\n", &o.intensityweight);
	if (ret != 1 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	if(verbose)	printf("Intensityweight=%lf\n", o.intensityweight);
    
    ignorereturnvalue=fgets(line, 512, fp);linenum++;
    ignorereturnvalue=fgets(line, 512, fp);linenum++;
    ignorereturnvalue=fgets(line, 512, fp);linenum++;

    junk = 0;
    for	(j=0; j<10; j++)	{
	    ret = fscanf(fp,"N_t_%i=%lf,%i\n", &junk, &o.NN[j], &fix_N_t[j]);
    	if (ret != 3 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
	for	(j=0; j<10; j++)	{
	    ret = fscanf(fp,"tauD_%i=%lf,%i\n", &junk, &o.tauD_i[j], &fix_tauD_i[j]);
    	if (ret != 3 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
    ret = fscanf(fp,"omega=%lf,%i\n", &o.omega, &fix_omega);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
	ret = fscanf(fp,"gamma_R=%lf,%i\n", &o.gamma_R, &fix_gamma_R);
    if(verbose) printf("gamma_R=%f, gamma_R_fix=%i\n", o.gamma_R, fix_gamma_R);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    j=0;
    ret = fscanf(fp,"g_inf0x0=%lf,%i\n", &o.g_infIxJ[j], &fix_g_infIxJ[j]);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    j++;
    ret = fscanf(fp,"g_inf1x1=%lf,%i\n", &o.g_infIxJ[j], &fix_g_infIxJ[j]);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    j++;
    ret = fscanf(fp,"g_inf2x2=%lf,%i\n", &o.g_infIxJ[j], &fix_g_infIxJ[j]);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    j++;
    ret = fscanf(fp,"g_inf0x1=%lf,%i\n", &o.g_infIxJ[j], &fix_g_infIxJ[j]);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    j++;
    ret = fscanf(fp,"g_inf1x2=%lf,%i\n", &o.g_infIxJ[j], &fix_g_infIxJ[j]);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    j++;
    ret = fscanf(fp,"g_inf2x0=%lf,%i\n", &o.g_infIxJ[j], &fix_g_infIxJ[j]);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    j++;
    ret = fscanf(fp,"g_inf0x1x2=%lf,%i\n", &o.g_infIxJ[j], &fix_g_infIxJ[j]);
    if(verbose) printf("ginf=%f, fix=%i, j=%i\n", o.g_infIxJ[6], fix_g_infIxJ[6],j);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    for	(j=0; j<10; j++)	{
	    ret = fscanf(fp,"tauf_%i=%lf,%i\n", &junk, &o.tauf_i[j], &fix_tauf_i[j]);
    	if (ret != 3 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
    for	(j=0; j<10; j++)	{
	    ret = fscanf(fp,"T0_%i=%lf,%i\n", &junk, &o.T0_i[j], &fix_T0_i[j]);
    	if (ret != 3 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
    for	(j=0; j<10; j++)	{
	    ret = fscanf(fp,"gamma_%i=%lf,%i\n", &junk, &o.gamma_i[j], &fix_gamma_i[j]);
    	if (ret != 3 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
        if (o.gamma_i[j] <= 0.0 )  {printf("gamma must be non-zero positive. line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
    for	(j=0; j<10; j++)	{
	    ret = fscanf(fp,"frac_gamma_%i=%lf,%i\n", &junk, &o.frac_gamma_i[j], &fix_frac_gamma_i[j]);
    	if (ret != 3 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
    for	(j=0; j<3; j++)	{
	    ret = fscanf(fp,"q_%i_%i=%lf,%i\n", &junk, &junk, &o.q_i[j], &fix_q_i[j]);
    	if (ret != 4 )  {printf("File format not valid, line %i q.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
    ret = fscanf(fp,"Xi={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}}\n", &o.eta[0], &o.eta[1], &o.eta[2], &o.eta[3], &o.eta[4], &o.eta[5], &o.eta[6], &o.eta[7], &o.eta[8]);
	if (ret != 9 )  {printf("File format not valid, line %i eta.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	ignorereturnvalue=fgets(line, 512, fp);
    linenum++;
    for	(j=0; j<3; j++)	{
	    ret = fscanf(fp,"bg_%i=%lf,%i\n", &junk, &o.bg[j], &fix_bg[j]);
    	if (ret != 3 )  {printf("File format not valid, line %i bg.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
	for	(j=0; j<3; j++)	{
	    ret = fscanf(fp,"Channel_Correction_%i=%lf,%i (Only Select N-1, where N=Num channels)\n", &junk, &o.cc[j], &fix_cc[j]);
	    if(verbose)	printf("jnk=%i, cc = %f, fix=%i\n", junk, o.cc[j], fix_cc[j]);
    	if (ret != 3 )  {printf("FFile format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }
    for	(j=0; j<3; j++)	{
	    ret = fscanf(fp,"occ_%i={%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf}\n", &junk, &o.occ[10*j+0], &o.occ[10*j+1], &o.occ[10*j+2], &o.occ[10*j+3], &o.occ[10*j+4], &o.occ[10*j+5], &o.occ[10*j+6], &o.occ[10*j+7], &o.occ[10*j+8], &o.occ[10*j+9]);
    	if (ret != 11 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    	linenum++;
    }

	
	ret = fscanf(fp,"Veff{0x0=%lf, 1x1=%lf, 2x2=%lf, 0x1=%lf, 1x2=%lf, 2x0=%lf}\n", &o.Veff[0],&o.Veff[1],&o.Veff[2],&o.Veff[3],&o.Veff[4],&o.Veff[5]);
	if (ret != 6) {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
    
    ignorereturnvalue=fgets(line, 512, fp);
    linenum++;
    ignorereturnvalue=fgets(line, 512, fp);
    linenum++;

    int num_entries = 0;
    ret = fscanf(fp,"num_entries=%i\n", &num_entries);
    if (ret != 1 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    if (num_entries > num_curves_to_consider)  {printf("recompile with num_curves_to_consider > 128.  Exiting.\n"); return 1; }

    char    filetag[128*num_curves_to_consider];
    float  f1[num_curves_to_consider];
    for(i=0; i<num_curves_to_consider; i++)	f1[i] = 0.0;
    for(i=0; i<num_entries; i++)    {
        ret = fscanf(fp,"tag=%s timepoint = %f\n", &filetag[128*i], &f1[i]);
		if(verbose) printf("Tag = %s, f1 = %f\n", &filetag[128*i], f1[i]);
        if (ret != 2 )  {printf("File format not valid, Data entry %i (line %i).  Exiting\n", i+1, linenum); fclose(fp); return 1;}    
        linenum++;
    }
    
    ignorereturnvalue=fgets(line, 512, fp);
    linenum++;
    ignorereturnvalue=fgets(line, 512, fp);
    linenum++;
    
	o.fit_start = 0;
	o.fit_stop = 0;
	ret = fscanf(fp,"fit_range_start = %i\nfit_range_stop = %i\n", &o.fit_start, &o.fit_stop);
	if (ret != 2 )  {printf("File format not valid, fit_start, fit_stop. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}    
	if(o.fit_stop > 0)	{printf("fit_stop should be negative.  Exiting\n"); return 1;}
	o.fit_stop = veclen + o.fit_stop;
	if(verbose)	printf("Data will be fit from i = (%i,%i]\n", o.fit_start, o.fit_stop);
	o.fit_start_triplecorr = 0;
	o.fit_stop_triplecorr = 0;
	ret = fscanf(fp,"fit_range_start_triplecorr = %i\nfit_range_stop_triplecorr = %i\n", &o.fit_start_triplecorr, &o.fit_stop_triplecorr);
	if (ret != 2 )  {printf("File format not valid, fit_start, fit_stop. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}    
	if(o.fit_stop_triplecorr > 0)	{printf("fit_stop_triplecorr should be negative.  Exiting\n"); return 1;}
	o.fit_stop_triplecorr = veclen_GGG + o.fit_stop_triplecorr;
	if(verbose)	printf("Data will be fit from i = (%i,%i]\n", o.fit_start_triplecorr, o.fit_stop_triplecorr);
	
	ignorereturnvalue=fgets(line, 512, fp);
    linenum++;

	ret = fscanf(fp,"k2=%lf,%lf,%lf\n", &o.k2[0], &o.k2[1], &o.k2[2]);
	if (ret != 3 )  {printf("File format not valid, k2. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	ret = fscanf(fp,"k3=%lf,%lf,%lf\n", &o.k3[0], &o.k3[1], &o.k3[2]);
	if (ret != 3 )  {printf("File format not valid, k3. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	ret = fscanf(fp,"k4=%lf,%lf,%lf\n", &o.k4[0], &o.k4[1], &o.k4[2]);
	if (ret != 3 )  {printf("File format not valid, k4. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	if(verbose)	printf("APD0 Correction Factors:\t%e\t%e\t%e\n", o.k2[0], o.k3[0], o.k4[0]);
	if(verbose)	printf("APD1 Correction Factors:\t%e\t%e\t%e\n", o.k2[1], o.k3[1], o.k4[1]);
	if(verbose)	printf("APD2 Correction Factors:\t%e\t%e\t%e\n", o.k2[2], o.k3[2], o.k4[2]);

	ignorereturnvalue=fgets(line, 512, fp);
	ignorereturnvalue = ignorereturnvalue;
    linenum++;

	double	Plot_Z_min = 0.0;
	double  Plot_Z_max = 0.0;
	int		Plot_Z_min_fit = 0;
	int		Plot_Z_max_fit = 0;
	ret = fscanf(fp,"Plot_Z_max = %lf,%i\n", &Plot_Z_max, &Plot_Z_max_fit);
	if (ret != 2 )  {printf("Plot_Z_max not valid, k3. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	ret = fscanf(fp,"Plot_Z_min = %lf,%i\n", &Plot_Z_min, &Plot_Z_min_fit);
	if (ret != 2 )  {printf("Plot_Z_max not valid, k3. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
    int verifystring = 0;
    ret = fscanf(fp, "Random string to verify file completeness:  %i\n", &verifystring);
    if  ( verifystring != 7392342)  {
        printf("Verification string (%i) not found - check file integrity\n", verifystring);
        return 1;
    }
    fclose(fp);
    
    //Verify fix parameters \in {0,1}
    for(j=0; j<6; j++)	{if(o.Gfix[j]<0 || o.Gfix[j]>1){printf("Gfix[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}}
    for(j=0; j<10; j++)	{if(fix_N_t[j]<0 || fix_N_t[j]>1){printf("fix_N_t[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}}
    for(j=0; j<10; j++)	{if(fix_tauD_i[j]<0 || fix_tauD_i[j]>1){printf("fix_tauD_i[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}}
    if(fix_omega<0 || fix_omega>1){printf("fix_omega needs to be 0 or 1.  Exiting\n"); fclose(fp); return 1;}
	if(fix_gamma_R<0 || fix_gamma_R>1){printf("fix_gamma_R needs to be 0 or 1.  Exiting\n"); fclose(fp); return 1;}
    for(j=0; j<7; j++)	{if(fix_g_infIxJ[j]<0 || fix_g_infIxJ[j]>1){printf("fix_g_infIxJ[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}} 
    for(j=0; j<10; j++)	{if(fix_tauf_i[j]<0 || fix_tauf_i[j]>1){printf("fix_tauf_i[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}}
    
    for(j=0; j<10; j++)	{if(fix_T0_i[j]<0 || fix_T0_i[j]>1){printf("fix_T0_i[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}}
    for(j=0; j<10; j++)	{if(fix_gamma_i[j]<0 || fix_gamma_i[j]>1){printf("fix_gamma_i[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}}
    for(j=0; j<10; j++)	{if(fix_frac_gamma_i[j]<0 || fix_frac_gamma_i[j]>1){printf("fix_frac_gamma_i[%i] needs to be 0 or 1.  Exiting\n", j); fclose(fp); return 1;}}



	//	Count Fit parameters, assign fit_param --> N map
	int	Gcount = 0;
	int	Icount = 0;
	int GGGcount = 0;
	for(j=0; j<3; j++)	o.Ion[j] = 0;
	for(j=0; j<6; j++)	if(o.Gfix[j]==1)	Gcount++;
	for(j=6; j<9; j++)	if(o.Gfix[j]==1)	GGGcount++;
	if(o.Gfix[0] == 1 || o.Gfix[3] == 1 || o.Gfix[5] == 1)	o.Ion[0] = 1;
	if(o.Gfix[1] == 1 || o.Gfix[3] == 1 || o.Gfix[4] == 1)	o.Ion[1] = 1;
	if(o.Gfix[2] == 1 || o.Gfix[4] == 1 || o.Gfix[5] == 1)	o.Ion[2] = 1;
	for(j=0; j<3; j++)	if(o.Ion[j]==1)	Icount++;
	if(o.Gfix[6] == 1) {
		Icount = 3;
		o.Ion[0] = 1;	o.Ion[1] = 1;	o.Ion[2] = 1;
	}
	if(verbose)	printf("Ifix[0] = %i\t Ifix[0] = %i\t Ifix[0] = %i\n", o.Ion[0], o.Ion[1], o.Ion[2]);

	int	pcount = 0;
	int	* pmask;
	double	*xtemp;
	pmask = malloc(10*100*sizeof(int));
	xtemp = malloc(10*100*sizeof(double));
	for(i=0; i<1000; i++)	pmask[i] = 0;
	for(i=0; i<1000; i++)	xtemp[i] = 0.0;
	j=0;
	for(i=0; i<10; i++)	{
		if(fix_N_t[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.NN[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<10; i++)	{
		if(fix_tauD_i[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.tauD_i[i];
			pcount++;
	}	}	j+= 10;
	if(fix_omega == 1)	{
		pmask[pcount] = j;
		xtemp[pcount] = o.omega;
		pcount++;
	}	j+= 10;
	for(i=0; i<7; i++)	{
		if(fix_g_infIxJ[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.g_infIxJ[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<10; i++)	{
		if(fix_tauf_i[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.tauf_i[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<10; i++)	{
		if(fix_T0_i[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.T0_i[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<10; i++)	{
		if(fix_gamma_i[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.gamma_i[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<10; i++)	{
		if(fix_frac_gamma_i[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.frac_gamma_i[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<3; i++)	{
		if(fix_q_i[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.q_i[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<3; i++)	{
		if(fix_bg[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.bg[i];
			pcount++;
	}	}	j+= 10;
	for(i=0; i<3; i++)	{
		if(fix_cc[i] == 1)	{
			pmask[pcount] = i+j;
			xtemp[pcount] = o.cc[i];
			pcount++;
	}	}	j+= 10;
	if(fix_gamma_R == 1)	{
		pmask[pcount] = j;
		xtemp[pcount] = o.gamma_R;
		pcount++;
	}	j+= 10;

	const int	p = pcount;
	o.p = p;
	gsl_vector * x = 0;
	x = gsl_vector_alloc (p);
	for (i=0; i<p; i++)		gsl_vector_set (x, i, xtemp[i]);
	o.pmask = pmask;

	o.verbose = verbose;

//	Size up first file, generate filelist
//  Read Data in From bin files
    char resultsfile01[80];	//Open file to write out data
    FILE * fpB;
    snprintf(resultsfile01, 50, "bin_clean%s.bin",&filetag[128*0]);
    fpB = fopen(resultsfile01, "rb");
    if (fpB == NULL)	{printf("***File %s Failed to Open!***\nQuiting\n", resultsfile01);	return 1;}
    int veclen_helper = 0;
    i=ret=fread(&veclen_helper, sizeof(int), 1, fpB);
    int	FiIDAQuant = 0;
    i=ret=fread(&FiIDAQuant, sizeof(int), 1, fpB);
    if(verbose)	printf("Veclen  (bin) = %i, Compiled Veclen = %i, FiIDAQuant = %i\n", veclen_helper, veclen, FiIDAQuant);
    if(veclen_helper != veclen)	{
    	printf("FCS data in file %s are not compatible with the way this program was compiled.  Recomplile this program with different n_integral and p_max values.\nCurrent pmax = %i\nCurrent n_integral = %i\n", resultsfile01, pmax, n_integral); 
    	return 1;
    }

    fseek (fpB , 0 , SEEK_END);
	long	file_size = ftell (fpB);
	if(verbose)	printf("File Size is %li\n", file_size);   
	long curves_per_file = file_size -  3*sizeof(int) - veclen_helper*sizeof(double);
	curves_per_file = curves_per_file / ((7+12*veclen_helper)*sizeof(double));
	if(verbose)	printf("Curves in file = %li (%li/%li)\n", curves_per_file, curves_per_file*((7+12*veclen_helper)*sizeof(double))+3*sizeof(int)+veclen_helper*sizeof(double),  file_size);
	if(fit_mode == 4)	{
		printf("***Start:To Assist in File Generation Only***\n");
		printf("num_entries=%li\n", curves_per_file);
		for(i=0; i<curves_per_file; i++)	printf("tag=%s timepoint = %.1f\n", &filetag[128*0], i*1.0);
		printf("***End:To Assist in File Generation Only***\n");
	}	
    fclose(fpB);
    
    char resultsfile01b[80];
    int veclen_a_helper_GGG = 0;
    if(o.Gfix[6] == 1 || o.Gfix[7] == 1 || o.Gfix[8] == 1)	{	//	Read in bin2 files for triple correlation
		FILE * fpB;
		snprintf(resultsfile01b, 50, "bin_clean%s.bin2",&filetag[128*0]);
		fpB = fopen(resultsfile01b, "rb");
		if (fpB == NULL)	{printf("***File %s Failed to Open!  Quitting***\nRun stream_cleaner_triple on triple-correlation data processed by stream_symm_triplecorr,\nor simply do not fit any of the G0x1x2 curves, i.e. G0x1x2=0, G1x2x0=0, G2x0x1=0\n", resultsfile01b); return 1;}
		
		i=ret=fread(&veclen_a_helper_GGG, sizeof(int), 1, fpB);
		FiIDAQuant = 0;
		i=ret=fread(&FiIDAQuant, sizeof(int), 1, fpB);
		if(verbose)	printf("Veclen (bin2) = %i, Compiled Veclen = %i, FiIDAQuant = %i\n", veclen_a_helper_GGG, veclen_GGG, FiIDAQuant);
		if(veclen_a_helper_GGG != veclen_GGG)	{
    	printf("Triple-Correlation FCS data in file %s are not compatible with the way this program was compiled.  Recomplile this program with different n_GGG and pmax_GGG values.\nCurrent pmax_GGG = %i\nCurrent n_GGG = %i\n", resultsfile01b, pmax_GGG, n_GGG); 
    	return 1;
    }
	
		fseek (fpB , 0 , SEEK_END);
		long	file_size = ftell (fpB);
		if(verbose)	printf("File Size is %li\n", file_size);   
		long curves_per_file = file_size -  3*sizeof(int) - veclen_GGG*sizeof(double);
		curves_per_file = curves_per_file / ((7+12*veclen_GGG+6*veclen_GGG*veclen_GGG)*sizeof(double));
		if(verbose)	printf("Curves in file = %li (%li/%li)\n", curves_per_file, curves_per_file*((7+12*veclen_GGG+6*veclen_GGG*veclen_GGG)*sizeof(double))+3*sizeof(int)+veclen_GGG*sizeof(double),  file_size);
		if(fit_mode == 4)	{
			printf("***Start:To Assist in File Generation Only***\n");
			printf("num_entries=%li\n", curves_per_file);
			for(i=0; i<curves_per_file; i++)	printf("tag=%s timepoint = %.1f\n", &filetag[128*0], i*1.0);
			printf("***End:To Assist in File Generation Only***\n");
		}	
		fclose(fpB);
    }
    
    
	//	***	Memory for data	***	//
    double  * tauarray;
    double	* time_array;
    double  * G_array;
    double  * G_stdev_array;
    double  * G_array_from_GGG;
    double	* G_stdev_array_from_GGG;
    double  * GGG_array;
    double  * GGG_stdev_array;
    double  * intensity0_array;
    double  * intensity0_stdev_array;
    double  * intensity1_array;
    double  * intensity1_stdev_array;
    double  * intensity2_array;
    double  * intensity2_stdev_array;
	double  * intensity0_array_GGG;
    double  * intensity1_array_GGG;
    double  * intensity2_array_GGG;
	double	* fits;
    
    
	tauarray = malloc(veclen_helper*sizeof(double));
	for(i=0; i<veclen_helper; i++) tauarray[i] = -6.6;
	
	time_array = malloc(num_entries*sizeof(double));
	intensity0_array = malloc(num_entries*sizeof(double));
	intensity0_stdev_array = malloc(num_entries*sizeof(double));
	intensity1_array = malloc(num_entries*sizeof(double));
	intensity1_stdev_array = malloc(num_entries*sizeof(double));
	intensity2_array = malloc(num_entries*sizeof(double));
	intensity2_stdev_array = malloc(num_entries*sizeof(double));
	intensity0_array_GGG = malloc(num_entries*sizeof(double));
	intensity1_array_GGG = malloc(num_entries*sizeof(double));
	intensity2_array_GGG = malloc(num_entries*sizeof(double));
	for(i=0; i<num_entries; i++) intensity0_array[i] = -6.6;
	for(i=0; i<num_entries; i++) intensity0_stdev_array[i] = -6.6;
	for(i=0; i<num_entries; i++) intensity1_array[i] = -6.6;
	for(i=0; i<num_entries; i++) intensity1_stdev_array[i] = -6.6;
	for(i=0; i<num_entries; i++) intensity2_array[i] = -6.6;
	for(i=0; i<num_entries; i++) intensity2_stdev_array[i] = -6.6;
	
	G_array 				= malloc(6*veclen_helper*sizeof(double));
	G_stdev_array 			= malloc(6*veclen_helper*sizeof(double));
	G_array_from_GGG 		= malloc(6*veclen_GGG*sizeof(double));
	G_stdev_array_from_GGG 	= malloc(6*veclen_GGG*sizeof(double));
	GGG_array		 		= malloc(3*veclen_GGG*veclen_GGG*sizeof(double));
	if(GGG_array==NULL)	printf("Malloc of GGG_array Failed\n");
	else	if(verbose) printf("Malloc of GGG_array Succeeded\n");
	GGG_stdev_array		 	= malloc(3*veclen_GGG*veclen_GGG*sizeof(double));
	for(i=0; i<6*veclen_helper; i++) G_array[i] = -6.6;
	for(i=0; i<6*veclen_helper; i++) G_stdev_array[i] = -6.6;
	for(i=0; i<6*veclen_GGG; i++) G_array_from_GGG[i] = -6.6;
	for(i=0; i<6*veclen_GGG; i++) G_stdev_array_from_GGG[i] = -6.6;
	for(i=0; i<3*veclen_GGG*veclen_GGG; i++) GGG_array[i] = -6.6;
	for(i=0; i<3*veclen_GGG*veclen_GGG; i++) GGG_stdev_array[i] = -6.6;
	fits 			= malloc((6*veclen_helper+3*veclen_GGG*veclen_GGG)*sizeof(double));
    o.w_resids 		= malloc((6*veclen_helper+3*veclen_GGG*veclen_GGG)*sizeof(double));
    for(i=0; i<6*veclen+3*veclen_GGG*veclen_GGG; i++)	fits[i] = 1.01;
	for(i=0; i<6*veclen+3*veclen_GGG*veclen_GGG; i++)	o.w_resids[i] = 1.01;
	
	//	***	Allocate memory before data-loops start	***	//
	
	o.m = malloc(veclen*sizeof(double));
	for(i=0; i<veclen; i++)	o.m[i] = 0.0;
	o.m_no_triplet = malloc(veclen*sizeof(double));
	for(i=0; i<veclen; i++)	o.m_no_triplet[i] = 0.0;
	o.mGGG = malloc(veclen_GGG*veclen_GGG*sizeof(double));
	for(i=0; i<veclen_GGG*veclen_GGG; i++)	o.mGGG[i] = 0.0;
	o.fits = malloc((6*veclen+3*veclen_GGG*veclen_GGG)*sizeof(double));
	for(i=0; i<6*veclen+3*veclen_GGG*veclen_GGG; i++)	o.fits[i] = 0.0;
    
    char resultsfile02[80];
    char datatag[80];
	char corr_mode[56] = "0x0  1x1  2x2  0x1  1x2  2x0  0x1x2";
	
    size_t nfull;
    nfull = veclen*6;		//Possibly +3 at later date to include intensity

    const gsl_multifit_fdfsolver_type *T;
	T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s = 0;
    gsl_matrix *covar = 0;
    int status;
    unsigned int iter = 0;
	size_t nFIT;
	nFIT = (o.fit_stop-o.fit_start)*Gcount+Icount+(o.fit_stop_triplecorr-o.fit_start_triplecorr)*(o.fit_stop_triplecorr-o.fit_start_triplecorr)*GGGcount;
    nfull = veclen*6 + veclen_GGG*veclen_GGG*3;

	gsl_multifit_function_fdf f;    
    f.f = &grp_tau_f;
    f.df = &grp_tau_df;
    f.fdf = &grp_tau_fdf;
    f.n = nFIT;
	f.params = &o; 

	double dof = 0.0;
	double chi = 0.0;
    double c = 0.0;
	

	gsl_vector * fff = 0;
	fff = gsl_vector_alloc (Gcount*(o.fit_stop-o.fit_start)+Icount+(o.fit_stop_triplecorr-o.fit_start_triplecorr)*(o.fit_stop_triplecorr-o.fit_start_triplecorr)*GGGcount);
	for (i=0; i<(Gcount*(o.fit_stop-o.fit_start)+Icount+(o.fit_stop_triplecorr-o.fit_start_triplecorr)*(o.fit_stop_triplecorr-o.fit_start_triplecorr)*GGGcount); i++)		gsl_vector_set (fff, i, 0.0);

    o.y = malloc(nfull*sizeof(double));
	o.sigma = malloc(nfull*sizeof(double));
	
	double	*	ptw;
	ptw = malloc(nptw*num_entries*sizeof(double));
	for(i=0; i<nptw*num_entries; i++)	ptw[i] = 0.0;
    
	int	current_entry = 0;
	for (current_entry=0; current_entry<num_entries; current_entry++)	{	//	Main loop: read, fit, write
		//	************************* Read in Data *************************	//
	
	    int numslices = 0;
        snprintf(resultsfile01, 50, "bin_clean%s.bin",&filetag[128*current_entry]);
        fpB = fopen(resultsfile01, "rb");
        if (fpB == NULL)	printf("***File %s Failed to Open!***\n", resultsfile01);
        else if(verbose) printf("File %s Opened\n", resultsfile01);
        fseek(fpB, 2*sizeof(int), SEEK_CUR);
        j = ret=fread(&numslices, sizeof(int), 1, fpB);
        if (j != 1)	printf("***File %s Read Problem! (1)  j = %i***\n", resultsfile01, j);
        if(current_entry==0)    {   //read in tau values
            j = ret=fread(tauarray, sizeof(double), veclen_helper, fpB);
            if (j != veclen_helper)	printf("***File %s Read Problem! (2)  j = %i***\n", resultsfile01, j);
        }
        else   fseek(fpB, veclen_helper*sizeof(double), SEEK_CUR); 
        //  ***	For Each timepoint, write Int0,Int1,Int2,SDInt0,SDInt1,SDInt2, 
        //      G0x0,G1x1,G2x2,G0x1,G1x2,G2x0,SDG0x0,SDG1x1,SDG2x2,SDG0x1,SDG1x2,SDG2x0	***	//
        
        fseek(fpB, (7+12*veclen_helper)*sizeof(double)*f1[current_entry], SEEK_CUR);
        //	***	read int's	***	//
        ret=fread(&time_array[current_entry], sizeof(double),1, fpB);
        ret=fread(&intensity0_array[current_entry], sizeof(double),1, fpB);
        ret=fread(&intensity1_array[current_entry], sizeof(double),1, fpB);
        ret=fread(&intensity2_array[current_entry], sizeof(double),1, fpB);
        ret=fread(&intensity0_stdev_array[current_entry], sizeof(double),1, fpB);
        ret=fread(&intensity1_stdev_array[current_entry], sizeof(double),1, fpB);
        ret=fread(&intensity2_stdev_array[current_entry], sizeof(double),1, fpB);
        
        for(j=0; j<6; j++)	ret=fread(&G_array[j*veclen_helper], sizeof(double),veclen_helper, fpB);
        for(j=0; j<6; j++)	ret=fread(&G_stdev_array[j*veclen_helper], sizeof(double),veclen_helper, fpB);

        fclose(fpB);
        
        if (o.Gfix[6]==1 || o.Gfix[7]==1 || o.Gfix[8]==1)	{

			snprintf(resultsfile01b, 50, "bin_clean%s.bin2",&filetag[128*current_entry]);
			fpB = fopen(resultsfile01b, "rb");
			if (fpB == NULL)	printf("***File %s Failed to Open!***\n", resultsfile01b);
			else if(verbose) printf("File %s Opened\n", resultsfile01b);
			fseek(fpB, 2*sizeof(int), SEEK_CUR);        //  Vestigial, Number of slices (timepoints)
			fseek(fpB, 1*sizeof(int), SEEK_CUR);		// numslices
			if(current_entry==0)    {   //read in tau values
				fseek(fpB, veclen_GGG*sizeof(double), SEEK_CUR);
			}
			else   fseek(fpB, veclen_GGG*sizeof(double), SEEK_CUR); 
			//	***	For Each timepoint, write Int0,Int1,Int2,SDInt0,SDInt1,SDInt2, 
			//      G0x0,G1x1,G2x2,G0x1,G1x2,G2x0,SDG0x0,SDG1x1,SDG2x2,SDG0x1,SDG1x2,SDG2x0	***	//
			
			fseek(fpB, (7+12*veclen_GGG+6*veclen_GGG*veclen_GGG)*sizeof(double)*f1[current_entry], SEEK_CUR);
			ret=fread(&time_array[current_entry], sizeof(double),1, fpB);
			ret=fread(&intensity0_array_GGG[current_entry], sizeof(double),1, fpB);
			ret=fread(&intensity1_array_GGG[current_entry], sizeof(double),1, fpB);
			ret=fread(&intensity2_array_GGG[current_entry], sizeof(double),1, fpB);
			fseek(fpB, 1*sizeof(double), SEEK_CUR);
			fseek(fpB, 1*sizeof(double), SEEK_CUR);
			fseek(fpB, 1*sizeof(double), SEEK_CUR);
			
			for(j=0; j<6; j++)	ret=fread(&G_array_from_GGG[j*veclen_GGG], sizeof(double),veclen_GGG, fpB);
			for(j=0; j<6; j++)	ret=fread(&G_stdev_array_from_GGG[j*veclen_GGG], sizeof(double),veclen_GGG, fpB);
			
			for(j=0; j<3; j++)	ret=fread(&GGG_array[j*veclen_GGG*veclen_GGG], sizeof(double),veclen_GGG*veclen_GGG, fpB);
			for(j=0; j<3; j++)	ret=fread(&GGG_stdev_array[j*veclen_GGG*veclen_GGG], sizeof(double),veclen_GGG*veclen_GGG, fpB);
				
			fclose(fpB);
		}
	
//	*************************** Fit Data ***************************	//

//	Transfer data to pass structure

	o.int0 = intensity0_array[current_entry];	// Currently no way to optimize Intensity...
	o.int1 = intensity1_array[current_entry];
	o.int2 = intensity2_array[current_entry];
	o.isd0 = intensity0_stdev_array[current_entry];
	o.isd1 = intensity1_stdev_array[current_entry];
	o.isd2 = intensity2_stdev_array[current_entry];
	
	if(verbose)	printf("Int0=%lf, Int1=%lf, Int2=%lf\n", o.int0,o.int1,o.int2);

	
	for (i=0; i<veclen*6; i++) {
        o.y[i] = G_array[i];
        o.sigma[i] = G_stdev_array[i];
    }

	for (i=0; i<veclen_GGG*veclen_GGG*3; i++) {
        o.y[i+veclen*6] = GGG_array[i];
        o.sigma[i+veclen*6] = GGG_stdev_array[i];
		if(o.sigma[i+veclen*6] <= 0) o.sigma[i+veclen*6] = 1000000.0;
		if(isnan(o.sigma[i+veclen*6])) o.sigma[i+veclen*6] = 2000000.0;
    }
	
	o.sigma[veclen_GGG*veclen_GGG*1+veclen*6] = 1000000.0;
	o.sigma[veclen_GGG*veclen_GGG*2+veclen*6] = 1000000.0;
	
	for(i=0; i<veclen_GGG*veclen_GGG; i+=veclen_GGG+1)	{
		if(i+veclen*6 > 0)				o.sigma[i+veclen*6] = 1000000.0;
	}
	for(i=veclen_GGG*veclen_GGG; i<2*veclen_GGG*veclen_GGG; i+=veclen_GGG+1)	{
		if(i+veclen*6 > 0)				o.sigma[i+veclen*6] = 1000000.0;
	}
	for(i=veclen_GGG*veclen_GGG*2; i<3*veclen_GGG*veclen_GGG; i+=veclen_GGG+1)	{
		if(i+veclen*6 > 0)				o.sigma[i+veclen*6] = 1000000.0;
	}
	
	
	grp_m_tau(x, &o, fff);

	if(verbose)	printf("Number of Parameters (2):\t%i\n", p);

	if(covar != 0)	{
		gsl_multifit_fdfsolver_free (s);
		gsl_matrix_free (covar);  
	}
    covar = gsl_matrix_alloc (p, p);
    for (i=0; i<p; i++)		gsl_vector_set (x, i, xtemp[i]);
    
    f.p = p;

    s = gsl_multifit_fdfsolver_alloc (T, nFIT, p);
    gsl_multifit_fdfsolver_set (s, &f, x);
  
  
	if(verbose)	{
		printf("x(strt)= {");
		for(i=0; i<p; i++)	printf(" %12.4le", gsl_vector_get(s->x,i));
		printf("}\n");
  	}
  	
  	
  	if(fit_mode == 0)	{
  	
		dof = nFIT - p;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);
			printf ("\t Status = %s, Red.Chi^2 = %.8f\n", gsl_strerror (status), pow(gsl_blas_dnrm2(s->f), 2.0)/dof);
                        if(verbose)	{
				printf("x(fit) = {");
				for(i=0; i<p; i++)	printf(" %12.4le", gsl_vector_get(s->x,i));
				printf("}\n");
			}
			if (status) break;
			status = gsl_multifit_test_delta (s->dx, s->x, 1e-9, 3e-5);
		} while (status == GSL_CONTINUE && iter < 2000);
		
		gsl_multifit_covar (s->J, 0.0, covar);
		chi = gsl_blas_dnrm2(s->f);	// Euclidean norm of f
		c = chi / sqrt(dof); 
		double	overallredchisq = pow(chi, 2.0) / dof;
		printf("chisq/dof = %g\n",  overallredchisq);
		
		//Fish out fitted parameters	
		j=0;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.NN[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.tauD_i[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.omega = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.g_infIxJ[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.tauf_i[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.T0_i[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.gamma_i[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.frac_gamma_i[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.q_i[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.bg[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.cc[o.pmask[i]-j] = gsl_vector_get(s->x, i);
		j+= 10;
		for(i=0; i<o.p; i++)
			if(o.pmask[i] > (j-1) && o.pmask[i] < (j+10))
				o.gamma_R = gsl_vector_get(s->x, i);
				
		
		grp_m_tau(s->x, &o, s->f);
		printf("N = {");
		for(i=0; i<10; i++)	printf(" %4.4le", o.NN[i]);
		printf("}\n");
		
		printf("Errors:");  
		for(i=0; i<p; i++)	printf("  %3.2e", c*sqrt(gsl_matrix_get(covar,i,i)));
		printf("\n");  
		printf("Fit Params:");
		for(i=0; i<p; i++)	printf("  %3.2e", gsl_vector_get(s->x,i));
		printf("\n");
	
		printf ("status = %i, %s\n", status, gsl_strerror (status));    
		printf("Overall Red. Chi Sq: %.4f\n", overallredchisq);
		
		for(i=0; i<10; i++)	{
			ptw[i+current_entry*nptw] = o.NN[i];
			ptw[i+current_entry*nptw+10] = o.tauD_i[i];
			ptw[i+current_entry*nptw+20] = -1.0;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == 0+i)	{
					ptw[i+current_entry*nptw+20] = c*sqrt(gsl_matrix_get(covar,k,k));
					k = 666;
			}	}
		}
	
		for(i=0; i<10; i++)	{
			if (isnan(ptw[i+current_entry*nptw+20]))  ptw[i+current_entry*nptw+20] = 1000000000;		
		}
	
		ptw[30+current_entry*nptw] = time_array[current_entry];	//time
		ptw[31+current_entry*nptw] = overallredchisq;	//ChiSq
		if(isnan(overallredchisq))	ptw[31+current_entry*nptw] = 1000000000;
		ptw[32+current_entry*nptw] = intensity0_array[current_entry];	//Int0
		ptw[33+current_entry*nptw] = intensity1_array[current_entry];	//Int1
		ptw[34+current_entry*nptw] = intensity2_array[current_entry];	//Int2
		ptw[35+current_entry*nptw] = intensity0_array[current_entry]+intensity1_array[current_entry]+intensity2_array[current_entry];	//IntT
		ptw[36+current_entry*nptw] = status;	//FitTerm
		ptw[37+current_entry*nptw] = o.q_i[0];	//q0
		ptw[38+current_entry*nptw] = o.q_i[1];	//q1
		ptw[39+current_entry*nptw] = o.q_i[2];	//q2
		ptw[40+current_entry*nptw] = o.gamma_R;	//q2
	
	
	// ******************** Print Fit Params to Screen ********************    //
		for(i=0; i<10; i++)	{
			printf("Species %i\n", i);
		
			j=0;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "  N=%10.3lf/%7.1le  ",o.NN[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(        "  N=%10.3lf          ",o.NN[i]);
			j= 40;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "tauf=%7.3lf/%7.1le  ",o.tauf_i[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "tauf=%7.3lf          ",o.tauf_i[i]);
			j=60;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "gamma=%7.3lf/%7.1le  ",o.gamma_i[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "gamma=%7.3lf          ",o.gamma_i[i]);
			j= 20;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j)	{
					printf(        "w=%7.3lf/%7.1le  ",o.omega,c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "w=%8.3lf          ",o.omega);
				
			j= 10;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "\n  tauD=%7.2lf/%7.1le  ",o.tauD_i[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "\n  tauD=%7.2lf          ",o.tauD_i[i]);
			j= 50;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "T0=%9.3lf/%7.1le  ",o.T0_i[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "T0=%9.3lf          ",o.T0_i[i]);
			j= 70;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "fracG=%7.3lf/%7.1le  \n",o.frac_gamma_i[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "fracG=%7.3lf          \n",o.frac_gamma_i[i]);
		
		}
		
		for(i=0; i<3; i++)	{
			j= 80;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "q_%i=%7.3lf/%7.1le  \t",i,o.q_i[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "q_%i=%7.3lf          \t",i,o.q_i[i]);
			j= 90;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "bg_%i=%7.3lf/%7.1le  ",i,o.bg[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "bg_%i=%7.3lf          ",i,o.bg[i]);
			j= 30;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i)	{
					printf(        "G%1ix%1i(inf)=%7.4lf/%7.1le  ",i,i,o.g_infIxJ[i],c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "G%1ix%1i(inf)=%7.4lf          ",i,i,o.g_infIxJ[i]);
			j= 30;
			for(k=0; k<o.p; k++)	{
				if(o.pmask[k] == j+i+3)	{
					printf(        "G%1ix%1i(inf)=%7.4lf/%7.1le  \n",i,(i+1)%3, o.g_infIxJ[i+3], c*sqrt(gsl_matrix_get(covar,k,k)));
					k = 666;
			}	}
			if(k<= 666)
				printf(   "G%1ix%1i(inf)=%7.4lf          \n",i,(i+1)%3, o.g_infIxJ[i+3]);
		}
		j= 100;
		for(k=0; k<o.p; k++)	{
			if(o.pmask[k] == j)	{
				printf(        "cc_0=%7.3lf/%7.1le    ",o.cc[0],c*sqrt(gsl_matrix_get(covar,k,k)));
				k = 666;
		}	}
		if(k<= 666)
			printf(   "cc_0=%7.3lf            ",o.cc[0]);
		j= 101;
		for(k=0; k<o.p; k++)	{
			if(o.pmask[k] == j)	{
				printf(        "cc_1=%7.3lf/%7.1le  ",o.cc[1],c*sqrt(gsl_matrix_get(covar,k,k)));
				k = 666;
		}	}
		if(k<= 666)
			printf(   "cc_1=%7.3lf          ",o.cc[1]);
		j= 102;
		for(k=0; k<o.p; k++)	{
			if(o.pmask[k] == j)	{
				printf(        "cc_2=%7.3lf/%7.1le\n",o.cc[2],c*sqrt(gsl_matrix_get(covar,k,k)));
				k = 666;
		}	}
		if(k<= 666)
			printf(   "cc_2=%7.3lf\n",o.cc[2]);
		j= 36;
		for(k=0; k<o.p; k++)	{
			if(o.pmask[k] == j)	{
				printf(        "G0x1x2(0,0)=%7.3le/%7.1le\t",o.g_infIxJ[6],c*sqrt(gsl_matrix_get(covar,k,k)));
				k = 666;
		}	}
		if(k<= 666)
			printf(   "G0x1x2(0,0)=%7.3lf\t",o.g_infIxJ[6]);
		j= 110;
		for(k=0; k<o.p; k++)	{
			if(o.pmask[k] == j)	{
				printf(        "gamma_R=%7.3lf/%7.1le\n",o.gamma_R,c*sqrt(gsl_matrix_get(covar,k,k)));
				k = 666;
		}	}
		if(k<= 666)
			printf(   "gamma_R=%7.3lf\n",o.gamma_R);
	

//	********************** Write Fits to File ********************	//
		
		double	G_max = 0;
		double	G_min = 10;
		
		if(verbose)	{
			for(j=0; j< 6; j++)	{
				if (o.Gfix[j]==1)	{
					for(i=o.fit_start+15; i<o.fit_stop; i++)	{
						if(G_array[j*veclen_helper+i] > G_max)	G_max = G_array[j*veclen_helper+i];
						if(G_array[j*veclen_helper+i] < G_min)	G_min = G_array[j*veclen_helper+i];
					}
				}
			}
			G_min = G_min - (G_max-G_min) * 0.025;
			G_max = G_max + (G_max-G_min) * 0.050;
			
			snprintf(resultsfile01, 50, "gtaul%s.%i.dat",&filetag[128*current_entry],current_entry);
			fpB = fopen(resultsfile01, "w");
			for(i=o.fit_start; i< o.fit_stop; i++){
				fprintf(fpB, "%e", tau_G(i));
				for(j=0; j< 6; j++)	{
					fprintf(fpB, "\t%e", G_array[j*veclen_helper+i]);
				}
				for(j=0; j< 6; j++)	{
					fprintf(fpB, "\t%e", G_stdev_array[j*veclen_helper+i]);
				}
				for(j=0; j< 6; j++)	{
					fprintf(fpB, "\t%e", o.fits[j*veclen_helper+i]);
				}
				for(j=0; j< 6; j++)	{
					fprintf(fpB, "\t%e", o.w_resids[j*veclen_helper+i]);
				}
				fprintf(fpB, "\n");
			}
			fclose(fpB);
				
			
			snprintf(resultsfile02, 50, "gnuplot%s.fits.%i.txt",&filetag[128*current_entry],current_entry);
			printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
			fpB = fopen(resultsfile02, "w");
			fprintf(fpB, "set logscale x\n");
			fprintf(fpB, "set xlabel \"tau (s)\"\n");
			fprintf(fpB, "set  ylabel \"G(tau)\"\n");
			fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
			fprintf(fpB, "unset colorbox\n");
			fprintf(fpB, "set term postscript eps enhanced color\n");
			snprintf(resultsfile02, 50, "gnuplot%s.fits.eps",&filetag[128*current_entry]);
			fprintf(fpB, "set output '%s'\n", resultsfile02);
			fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
			int	first_plot = 1;
			for(j=0; j< 6; j++)	{
				if (o.Gfix[j]==1)	{
					if(first_plot == 1)	first_plot = 0;
					else	fprintf(fpB, ", ");
					snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
					i=0;
					do {
						if(datatag[i] == '_') datatag[i] = '-';
						i++;
					}	while (datatag[i] != '\0');
					fprintf(fpB, "\"%s\" using 1:%i title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6.0);
				}
			}
			for(j=0; j< 6; j++)	{
				if (o.Gfix[j]==1)	{
					fprintf(fpB, ", ");
					snprintf(datatag, 80, "fit_%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
					i=0;
					do {
						if(datatag[i] == '_') datatag[i] = '-';
						i++;
					}	while (datatag[i] != '\0');
					fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*6, datatag, (float)j / (float) 6.0);
				}
			}
			fprintf(fpB, "\n");
			fprintf(fpB, "set term x11\n");
			fprintf(fpB, "replot\n");
		
			fclose(fpB);
			
			snprintf(resultsfile02, 50, "gnuplot%s.errors.%i.txt",&filetag[128*current_entry],current_entry);
			printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
			fpB = fopen(resultsfile02, "w");
			fprintf(fpB, "set logscale x\n");
			fprintf(fpB, "set xlabel \"tau (s)\"\n");
			fprintf(fpB, "set  ylabel \"G(tau)\"\n");
			fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
			fprintf(fpB, "unset colorbox\n");
			fprintf(fpB, "set term postscript eps enhanced color\n");
			snprintf(resultsfile02, 50, "gnuplot%s.fits.eps",&filetag[128*current_entry]);
			fprintf(fpB, "set output '%s'\n", resultsfile02);
			fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
			first_plot = 1;
			for(j=0; j< 6; j++)	{
				if (o.Gfix[j]==1)	{
					if(first_plot == 1)	first_plot = 0;
					else	fprintf(fpB, ", ");
					snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
					i=0;
					do {
						if(datatag[i] == '_') datatag[i] = '-';
						i++;
					}	while (datatag[i] != '\0');
					fprintf(fpB, "\"%s\" using 1:%i:%i with yerrorbars title \"%s\" lt palette frac %.2f", resultsfile01, j+2, j+2+6, datatag, (float)j / (float) 6.0);
				}
			}
			fprintf(fpB, "\n");
			fprintf(fpB, "set term x11\n");
			fprintf(fpB, "replot\n");
		
			fclose(fpB);
		
			snprintf(resultsfile02, 50, "gnuplot%s.resids.%i.txt",&filetag[128*current_entry],current_entry);
			printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
			fpB = fopen(resultsfile02, "w");
			fprintf(fpB, "reset\n");
			fprintf(fpB, "set logscale x\n");
			fprintf(fpB, "set xlabel \"tau (s)\"\n");
			fprintf(fpB, "set  ylabel \"G(tau)\"\n");
			fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
			fprintf(fpB, "unset colorbox\n");
			fprintf(fpB, "set term postscript eps enhanced color\n");
			snprintf(resultsfile02, 50, "gnuplot%s.fits.eps",&filetag[128*current_entry]);
			fprintf(fpB, "set output '%s'\n", resultsfile02);
			fprintf(fpB, "plot [%2.2e:%2.2e][] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1);
			first_plot = 1;
				
			for(j=0; j< 6; j++)	{
				if (o.Gfix[j]==1)	{
					if(first_plot == 1)	first_plot = 0;
					else	fprintf(fpB, ", ");
					snprintf(datatag, 80, "w_resid_%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
					i=0;
					do {
						if(datatag[i] == '_') datatag[i] = '-';
						i++;
					}	while (datatag[i] != '\0');
					fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+3*6, datatag, (float)j / (float) 6.0);
				}
			}
			fprintf(fpB, "\n");
			fprintf(fpB, "set term x11\n");
			fprintf(fpB, "replot\n");
		
			fclose(fpB);

	
		}	// 	End if(verbose)	write files
		}   // End if(fit_mode == 0)	{		
		
	//	*************** Plot Fit functions using guesses ***************	//
	if(fit_mode == 2)	{	//Plot Guess
		grp_m_tau(s->x, &o, s->f);
	} //End Plot Guess
	//	************************ Just Plot Data ************************	//
	if(fit_mode == 1)	{	//Plot Data
	
		G_max = 0;
		G_min = 10;
		
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				for(i=o.fit_start+n_integral/2; i<o.fit_stop; i++)	{
					if(G_array[j*veclen_helper+i] > G_max)	G_max = G_array[j*veclen_helper+i];
					if(G_array[j*veclen_helper+i] < G_min)	G_min = G_array[j*veclen_helper+i];
				}
			}
		}
		G_min = G_min - (G_max-G_min) * 0.025;
		G_max = G_max + (G_max-G_min) * 0.150;
		
		snprintf(resultsfile01, 50, "gtaul_Gs%s.%i.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		for(i=o.fit_start; i< o.fit_stop; i++){
			fprintf(fpB, "%e", tau_G(i));
			for(j=0; j< 6; j++)	{
				fprintf(fpB, "\t%e", G_array[j*veclen_helper+i]);
			}
			for(j=0; j< 6; j++)	{
				fprintf(fpB, "\t%e", G_stdev_array[j*veclen_helper+i]);
			}
			for(j=0; j< 6; j++)	{
				fprintf(fpB, "\t%e", o.fits[j*veclen_helper+i]);
			}
			for(j=0; j< 6; j++)	{
				fprintf(fpB, "\t%e", o.w_resids[j*veclen_helper+i]);
			}
			fprintf(fpB, "\n");
		}
		fclose(fpB);
		
		snprintf(resultsfile02, 50, "gnuplot%s.fits.%i.txt",&filetag[128*current_entry],current_entry);
		printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
		fpB = fopen(resultsfile02, "w");
		fprintf(fpB, "set logscale x\n");
		fprintf(fpB, "set xlabel \"tau (s)\"\n");
		fprintf(fpB, "set  ylabel \"G(tau)\"\n");
		fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
		fprintf(fpB, "unset colorbox\n");
		fprintf(fpB, "set term postscript eps enhanced color\n");
		snprintf(resultsfile02, 50, "gnuplot%s.fits.eps",&filetag[128*current_entry]);
		fprintf(fpB, "set output '%s'\n", resultsfile02);
		if(o.fit_start == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(1)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		else fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		int	first_plot = 1;
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				if(first_plot == 1)	first_plot = 0;
				else	fprintf(fpB, ", ");
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6.0);
			}
		}
		fprintf(fpB, "\n");
		fprintf(fpB, "set term x11\n");
		fprintf(fpB, "replot\n");
	
		fclose(fpB);

		snprintf(resultsfile01, 50, "gtaul_Gs_fromGGG%s.%i.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		for(i=o.fit_start_triplecorr; i< o.fit_stop_triplecorr; i++){
			fprintf(fpB, "%e", tau_GGG(i));
			for(j=0; j< 6; j++)	{
				fprintf(fpB, "\t%e", G_array_from_GGG[j*veclen_GGG+i]);
			}
			for(j=0; j< 6; j++)	{
				fprintf(fpB, "\t%e", G_stdev_array_from_GGG[j*veclen_GGG+i]);
			}
			fprintf(fpB, "\n");
		}
		fclose(fpB);
		
		if(verbose) printf("Veclen = %i\tveclen_GGG_helper = %i\n", veclen_helper, veclen_GGG);
		snprintf(resultsfile02, 50, "gnuplot%s_fromGGG.%i.txt",&filetag[128*current_entry],current_entry);
		printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
		fpB = fopen(resultsfile02, "w");
		fprintf(fpB, "set logscale x\n");
		fprintf(fpB, "set xlabel \"tau (s)\"\n");
		fprintf(fpB, "set  ylabel \"G(tau)\"\n");
		fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
		fprintf(fpB, "unset colorbox\n");
		fprintf(fpB, "set term postscript eps enhanced color\n");
		snprintf(resultsfile02, 50, "gnuplot%s.fits.eps",&filetag[128*current_entry]);
		fprintf(fpB, "set output '%s'\n", resultsfile02);
		if(o.fit_start_triplecorr == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, G_min, G_max);
		else	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, G_min, G_max);
		first_plot = 1;
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				if(first_plot == 1)	first_plot = 0;
				else	fprintf(fpB, ", ");
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6.0);
			}
		}
		fprintf(fpB, "\n");
		fprintf(fpB, "set term x11\n");
		fprintf(fpB, "replot\n");
	
		fclose(fpB);
		
		
		G_max = 0;
		G_min = 10;
		
		for(j=0; j<3; j++)	{
			if (o.Gfix[j+6]==1)	{
				for(i=o.fit_start_triplecorr+n_GGG/2; i<o.fit_stop_triplecorr; i++)	{
					for(w=o.fit_start_triplecorr+n_GGG/2; w<o.fit_stop_triplecorr; w++)	{
						if(GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i] > G_max)	G_max = GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i];
						if(GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i] < G_min)	G_min = GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i];
					}
				}
			}
		}
		G_min = G_min - (G_max-G_min) * 0.025;
		G_max = G_max + (G_max-G_min) * 0.150;

		
		snprintf(resultsfile01, 50, "gtaul_GGGs%s.%i.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		for(i=o.fit_start_triplecorr; i< o.fit_stop_triplecorr; i++){
			for(w=o.fit_start_triplecorr; w< o.fit_stop_triplecorr; w++){
				fprintf(fpB, "%9.9e\t%9.9e", tau_GGG(w),tau_GGG(i));
				for(j=0; j< 3; j++)	{
					fprintf(fpB, "\t%e", GGG_array[j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				}
				for(j=0; j< 3; j++)	{
					fprintf(fpB, "\t%e", GGG_stdev_array[j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				}
				for(j=0; j< 3; j++)	{
					fprintf(fpB, "\t%e", o.fits[6*veclen_helper+j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				}
				for(j=0; j< 3; j++)	{
					fprintf(fpB, "\t%e", o.w_resids[6*veclen_helper+j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				}
				fprintf(fpB, "\n");
			}
			fprintf(fpB, "\n");
		}
		fclose(fpB);
		
		// ***	Nscaled .dat Data File	*** //
		snprintf(resultsfile01, 50, "gtaug%s.%i.Nscaled.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		for(i=o.fit_start; i< o.fit_stop; i++){
			fprintf(fpB, "%e", tau_G(i));
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", (G_array[j*veclen_helper+i]-1.0)*1.0);
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", (G_stdev_array[j*veclen_helper+i])*1.0);
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", (o.fits[j*veclen_helper+i]-1.0)*1.0);
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", o.w_resids[j*veclen_helper+i]);
			fprintf(fpB, "\n");
		}
		fclose(fpB);
		
		//	***	Resid.s .dat Data File	***//
		snprintf(resultsfile01, 50, "gtaug%s.%i.resids.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		for(i=o.fit_start; i< o.fit_stop; i++){
			fprintf(fpB, "%e", tau_G(i));
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", o.w_resids[j*veclen_helper+i]);
			fprintf(fpB, "\n");
		}
		fclose(fpB);
		
		int i2,tau2,tau;
		int	num_selected_curves = 4;
		//	***	Write selected decays to a seperate file for plotting  *** //
		//	Array index order:	chirality, tau_2 slice, tau_1 vs tau_2 const, [tau_1,tau_2,G(tau_1,tau_2)] 
		snprintf(resultsfile01b, 50, "gtaul_select_GGGs%s.%i.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01b, "w");
		for(tau2 = 0; tau2 < veclen_GGG; tau2++)	{	//	tau_2 down vertical axis
			fprintf(fpB, "%i", tau2);
			tau = 0;	//	2x0x1 G(tau^0_1,X) | 1x2x0 G(X,tau^0_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[2*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[2*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			}
			tau = 0;	//	2x0x1 G(tau^0_1,X) | 1x2x0 G(X,tau^0_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[1*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[1*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			}
			tau = 0;	//	0x1x2 G(tau^1_1,X) | 2x0x1 G(X,tau^1_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[0*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[0*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			}
			tau = 0;	//	0x1x2 G(tau^1_1,X) | 2x0x1 G(X,tau^1_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[2*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[2*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			}
			tau = 0;	//	1x2x0 G(tau^0_2,X) | 0x1x2 G(X,tau^2_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[1*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[1*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			}
			tau = 0;	//	1x2x0 G(tau^2_1,X) | 0x1x2 G(X,tau^2_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[0*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[0*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			}
			fprintf(fpB, "\n");
		}
		fclose(fpB);

		snprintf(resultsfile01, 50, "gtaul_GGGs%s.%i.dat",&filetag[128*current_entry],current_entry);
		snprintf(resultsfile02, 50, "gnuplot%s_GGGs.%i.txt",&filetag[128*current_entry],current_entry);
		char	colorlib[90];
		char	datatag0[90];
		char	datatag1[90];
		printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
		fpB = fopen(resultsfile02, "w");
		fprintf(fpB, "Gmax = %f\n",G_max);
		fprintf(fpB, "Gmin = %f\n",G_min);
		fprintf(fpB, "set logscale x\n");
		fprintf(fpB, "set logscale y\n");
		fprintf(fpB, "set zlabel \"G^{0x1x2}(tau_1,tau_2)\" offset -1.75,1.5,0 rotate by 90 font \"Helvetica, 12\"\n");
		fprintf(fpB, "set ticslevel 0\n");
		fprintf(fpB, "set view 78,135\n");
		fprintf(fpB, "set title font \"Helvetica, 10\" \n");
		fprintf(fpB, "set palette defined (0 \"grey70\", 0.5 \"navy\", 1.0 \"orange\")\n");
		fprintf(fpB, "unset colorbox\n");
		fprintf(fpB, "set term postscript eps enhanced color\n");
		snprintf(resultsfile02, 50, "multiGGGt%s.%i.eps",&filetag[128*current_entry],current_entry);
		fprintf(fpB, "set output '%s'\n", resultsfile02);
		fprintf(fpB, "set multiplot\n");
		fprintf(fpB, "set size 0.7,1.0\n");
		fprintf(fpB, "set style line 10 linetype 1 linecolor rgb \"gray\"\n");
		fprintf(fpB, "set bmargin 0\n");
		fprintf(fpB, "set tmargin 2\n");
		fprintf(fpB, "set lmargin 10\n");
		fprintf(fpB, "set rmargin 2\n");
		fprintf(fpB, "set ztics font \"Helvetica, 10\"\n");
		fprintf(fpB, "unset hidden3d\n");
		for(j=0; j< 3; j++)	{
			if (o.Gfix[j+6]==1)	{
				fprintf(fpB, "set origin %2.2f,0.8\n",-0.1+j*0.5);
				fprintf(fpB, "set xtics offset 1.5,-0.5 font \"Helvetica, 10\"\n");
				fprintf(fpB, "set ytics offset -1.5,-0.5 font \"Helvetica, 10\"\n");
				fprintf(fpB, "set xlabel \"tau_1 (s)\" offset 1.5,-1.0 font \"Helvetica, 12\"\n");
				fprintf(fpB, "set ylabel \"tau_2 (s)\" offset -1.0,-1.0 font \"Helvetica, 12\"\n");
				if(j==1)	{
					fprintf(fpB, "set xtics offset 1.5,-1.75 font \"Helvetica, 10\"\n");
					fprintf(fpB, "set ytics offset -1.5,-1.75 font \"Helvetica, 10\"\n");
					fprintf(fpB, "set xlabel \"tau_1 (s)\" offset 1.5,-2.25 font \"Helvetica, 12\"\n");
					fprintf(fpB, "set ylabel \"tau_2 (s)\" offset -1.0,-2.25 font \"Helvetica, 12\"\n");
				}
				if(o.fit_start_triplecorr==0)	fprintf(fpB, "splot [%2.2e:%2.2e][%2.2e:%2.2e][Gmin:Gmax] ", tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				else fprintf(fpB, "splot [%2.2e:%2.2e][%2.2e:%2.2e][Gmin:Gmax] ", tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:2:%i with lines title \"\" lw 0 lt palette frac 0", resultsfile01, j+3);
				for(k=1; k<(num_selected_curves+1.0); k++)	{
					if(j==2)	snprintf(colorlib, 30, "blue");
					if(j==0)	snprintf(colorlib, 30, "orange");
					if(j==1)	snprintf(colorlib, 30, "red");
					fprintf(fpB, ", \"%s\" using %i:%i:%i with lines title \"\" lc rgb \"%s\" lw 2 lt 1", resultsfile01b, ((j+1)%3)*6*(num_selected_curves+1)+k*3+2,((j+1)%3)*6*(num_selected_curves+1)+k*3+3,((j+1)%3)*6*(num_selected_curves+1)+k*3+4,colorlib);
					if(j==1)	snprintf(colorlib, 30, "blue");
					if(j==2)	snprintf(colorlib, 30, "orange");
					if(j==0)	snprintf(colorlib, 30, "red");
					fprintf(fpB, ", \"%s\" using %i:%i:%i with lines title \"\" lc rgb \"%s\" lw 2 lt 2", resultsfile01b, (((j+2)%3)*6+3)*(num_selected_curves+1)+k*3+2,(((j+2)%3)*6+3)*(num_selected_curves+1)+k*3+3,(((j+2)%3)*6+3)*(num_selected_curves+1)+k*3+4,colorlib);
				}
				fprintf(fpB, "\n");
				if(j==0)	{fprintf(fpB, "set zlabel \"G^{1x2x0}(tau_1,tau_2)\" offset 3.0,1.5,0 rotate by 90 font \"Helvetica, 12\"\n"); fprintf(fpB, "set ztics format \"\"\n");}
				if(j==1)	 fprintf(fpB, "set zlabel \"G^{2x0x1}(tau_1,tau_2)\" offset 3.0,1.5,0 rotate by 90 font \"Helvetica, 12\"\n");
			}
		}
		fprintf(fpB, "set size 0.6,0.6\n");
		fprintf(fpB, "set ztics font \"Helvetica, 10\"\n");
		fprintf(fpB, "set palette defined ( 0 \"black\", 1 \"blue\", 2 \"green\", 3 \"red\" )\n");
		fprintf(fpB, "set ytics offset 0,0 font \"Helvetica, 10\"\n");
		fprintf(fpB, "set xtics offset 0,0 font \"Helvetica, 10\"\n");
		fprintf(fpB, "set ylabel \"G(tau)\" offset 1.75,0 font \"Helvetica, 12\"\n");
		fprintf(fpB, "unset xlabel\n");
		fprintf(fpB, "set xlabel \"tau (s)\" offset 0,0 font \"Helvetica, 12\"\n");
		fprintf(fpB, "unset logscale\n");
		fprintf(fpB, "set logscale x\n");
		for(j=0; j< 3; j++)	{
			if (o.Gfix[j+6]==1)	{
				fprintf(fpB, "set origin %2.2f,0.3\n",-0.05+j*0.5);
				if(o.fit_start_triplecorr == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][Gmin:Gmax] ", tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				else	fprintf(fpB, "plot [%2.2e:%2.2e][Gmin:Gmax] ", tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				snprintf(datatag0, 80, "%ix%ix%i(tau,0)s%s_%.2f", (j+2)%3,(j+3)%3,(j+4)%3, &filetag[128*current_entry], f1[current_entry]);
				snprintf(datatag1, 80, "%ix%ix%i(0,tau)s%s_%.2f", (j+1)%3,(j+2)%3,(j+3)%3, &filetag[128*current_entry], f1[current_entry]);
				i=0;
				do {
					if(datatag0[i] == '_') datatag0[i] = '-';
					i++;
				}	while (datatag0[i] != '\0');
				i=0;
				do {
					if(datatag1[i] == '_') datatag1[i] = '-';
					i++;
				}	while (datatag1[i] != '\0');
				k=0;
				snprintf(colorlib, 30, "black");
				fprintf(fpB, "\"%s\" using %i:%i with lines title \"%s\" lc rgb \"%s\" lw 2 lt 1, ", resultsfile01b, ((j)%3)*6*(num_selected_curves+1)+k*3+2,((j)%3)*6*(num_selected_curves+1)+k*3+4,datatag0,colorlib);
				fprintf(fpB, "\"%s\" using %i:%i with lines title \"%s\" lc rgb \"%s\" lw 2 lt 2, ", resultsfile01b, (((j)%3)*6+3)*(num_selected_curves+1)+k*3+3,(((j)%3)*6+3)*(num_selected_curves+1)+k*3+4,datatag1,colorlib);
				if(j==0)	snprintf(colorlib, 30, "blue");
				if(j==1)	snprintf(colorlib, 30, "orange");
				if(j==2)	snprintf(colorlib, 30, "red");
				snprintf(datatag0, 80, "%ix%ix%i(tau,X)s%s_%.2f", (j+2)%3,(j+3)%3,(j+4)%3, &filetag[128*current_entry], f1[current_entry]);
				snprintf(datatag1, 80, "%ix%ix%i(X,tau)s%s_%.2f", (j+1)%3,(j+2)%3,(j+3)%3, &filetag[128*current_entry], f1[current_entry]);
				i=0;
				do {
					if(datatag0[i] == '_') datatag0[i] = '-';
					i++;
				}	while (datatag0[i] != '\0');
				i=0;
				do {
					if(datatag1[i] == '_') datatag1[i] = '-';
					i++;
				}	while (datatag1[i] != '\0');
				for(k=1; k<(num_selected_curves+1.0); k++)	{
					if(k>1)	snprintf(datatag0, 80, " ");
					fprintf(fpB, "\"%s\" using %i:%i with lines title \"\" lc rgb \"%s\" lw 2 lt 1", resultsfile01b, ((j)%3)*6*(num_selected_curves+1)+k*3+2,((j)%3)*6*(num_selected_curves+1)+k*3+4,colorlib);
					fprintf(fpB, ", ");
				}
				for(k=1; k<(num_selected_curves+1.0); k++)	{
					if(k>1)	snprintf(datatag1, 80, " ");
					fprintf(fpB, "\"%s\" using %i:%i with lines title \"\" lc rgb \"%s\" lw 2 lt 2", resultsfile01b, (((j)%3)*6+3)*(num_selected_curves+1)+k*3+3,(((j)%3)*6+3)*(num_selected_curves+1)+k*3+4,colorlib);
					if(k==num_selected_curves)	fprintf(fpB, "\n");
					else	fprintf(fpB, ", ");
				}
				if(j==0)	{
					fprintf(fpB, "unset ylabel\n");
					fprintf(fpB, "set ytics format \"\"\n");
				}
			}
		}
		fprintf(fpB, "unset multiplot\n");
		fprintf(fpB, "set term x11\n");
		fprintf(fpB, "replot\n");
		fprintf(fpB, "reset\n");
	
		fclose(fpB);
	
		//	***	Plotting file for "all-in-one" multiplot style from gtauglobal	***	//
        G_max = 0;
		G_min = 10;
		
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				for(i=o.fit_start+n_integral/2; i<o.fit_stop; i++)	{
					if(G_array[j*veclen_helper+i] > G_max)	G_max = G_array[j*veclen_helper+i];
					if(G_array[j*veclen_helper+i] < G_min)	G_min = G_array[j*veclen_helper+i];
				}
			}
		}
		G_min = G_min - (G_max-G_min) * 0.025;
		G_max = G_max + (G_max-G_min) * 0.150;
		if (verbose) printf ("Gmax, Gmin = %f , %f\n", G_max, G_min);
		snprintf(resultsfile01, 50, "gtaul_Gs%s.%i.dat",&filetag[128*current_entry],current_entry);
		snprintf(resultsfile02, 50, "gnuplot%s_Gs.%i.txt",&filetag[128*current_entry],current_entry);
		printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
		fpB = fopen(resultsfile02, "w");
		fprintf(fpB, "reset\n");
		fprintf(fpB, "set term postscript eps enhanced color size 9, 6.3 solid\n");
		snprintf(resultsfile02, 50, "multiG.gnuplot%s.eps",&filetag[128*current_entry]);
		fprintf(fpB, "set output '%s'\n", resultsfile02);
		fprintf(fpB, "set multiplot\n");
		fprintf(fpB, "set size 0.5,0.5\n");
		fprintf(fpB, "set bmargin 0\n");
		fprintf(fpB, "set tmargin 2\n");
		fprintf(fpB, "set lmargin 10\n");
		fprintf(fpB, "set rmargin 2\n");
		fprintf(fpB, "set origin 0.0,0.5\n");
		fprintf(fpB, "set logscale x\n");
		fprintf(fpB, "set xlabel \"\"\n");
		fprintf(fpB, "set format x \"\"\n");
		fprintf(fpB, "set  ylabel \"G(tau)\"\n");
		fprintf(fpB, "set palette defined ( 0 \"blue\", 1 \"goldenrod\", 2 \"red\", 3 \"green\", 4 \"orange\", 5 \"violet\" )\n");
		fprintf(fpB, "unset colorbox\n");
		
		if(o.fit_start == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(1)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		else	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6);
				fprintf(fpB, ", ");
			}
		}
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				snprintf(datatag, 80, "fit_%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*6, datatag, (float)j / (float) 6);
				if (j != 6 -1) fprintf(fpB, ", ");
			}
		}
		fprintf(fpB, "\n");
	
		fprintf(fpB, "set origin 0.0,0.0\n");
		fprintf(fpB, "set key\n");
		fprintf(fpB, "set bmargin 3\n");
		fprintf(fpB, "set tmargin 0\n");
		fprintf(fpB, "set lmargin 10\n");
		fprintf(fpB, "set rmargin 2\n");
		fprintf(fpB, "set xlabel \"tau (s)\"\n");
		fprintf(fpB, "set format x\n");
		fprintf(fpB, "set  ylabel \"G(tau)\"\n");
		fprintf(fpB, "unset colorbox\n");
	
		if(o.fit_start == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(1)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		else fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:%i:%i with yerrorbars title \"%s\" lt palette frac %.2f", resultsfile01, j+2, j+2+6, datatag, (float)j / (float) 6);
				if (j != 6 -1) fprintf(fpB, ", ");
			}
		}
		fprintf(fpB, "\n");
	
		fprintf(fpB, "set origin 0.5,0.5\n");
		fprintf(fpB, "set nokey\n");
		fprintf(fpB, "set xlabel \"\"\n");
		fprintf(fpB, "set format x \"\"\n");
		fprintf(fpB, "set bmargin 0\n");
		fprintf(fpB, "set tmargin 2\n");
		fprintf(fpB, "set lmargin 8\n");
		fprintf(fpB, "set rmargin 2\n");
		fprintf(fpB, "set format x \"\"\n");
		fprintf(fpB, "set  ylabel \"G(tau) (Norm.)\"\n");
		fprintf(fpB, "unset colorbox\n");
		fprintf(fpB, "set ticslevel 0.0\n");
	
		snprintf(resultsfile01, 50, "gtaug%s.%i.Nscaled.dat",&filetag[128*current_entry],current_entry);
		if(o.fit_start == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][-0.05:1.2] ", tau_G(1)*0.9, tau_G(o.fit_stop-1)*1.1);
		else	fprintf(fpB, "plot [%2.2e:%2.2e][-0.05:1.2] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1);
		for(j=0; j< 6; j++)	{
			snprintf(datatag, 80, "scaled_%3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
			i=0;
			do {
				if(datatag[i] == '_') datatag[i] = '-';
				i++;
			}	while (datatag[i] != '\0');
			fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6);
			fprintf(fpB, ", ");
		}
		for(j=0; j< 6; j++)	{
			snprintf(datatag, 80, "scaledfit-%3s%s-%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
			i=0;
			do {
				if(datatag[i] == '_') datatag[i] = '-';
				i++;
			}	while (datatag[i] != '\0');
			fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*6, datatag, (float)j / (float) 6);
			if (j != 6 -1) fprintf(fpB, ", ");
		}
		fprintf(fpB, "\n");
	
		fprintf(fpB, "set origin 0.5,0.0\n");
		fprintf(fpB, "set nokey\n");
		fprintf(fpB, "set bmargin 3\n");
		fprintf(fpB, "set tmargin 0\n");
		fprintf(fpB, "set lmargin 8\n");
		fprintf(fpB, "set rmargin 2\n");
		fprintf(fpB, "set xlabel \"tau (s)\"\n");
		fprintf(fpB, "set format x\n");
		fprintf(fpB, "set  ylabel \"Weighted Residuals\"\n");
		fprintf(fpB, "unset colorbox\n");
		fprintf(fpB, "set ticslevel 0.0\n");
	
		snprintf(resultsfile01, 50, "gtaug%s.%i.resids.dat",&filetag[128*current_entry],current_entry);
	
		if(o.fit_start == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][] ", tau_G(1)*0.9, tau_G(o.fit_stop-1)*1.1);
		else	fprintf(fpB, "plot [%2.2e:%2.2e][] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1);
		for(j=0; j< 6; j++)	{
			snprintf(datatag, 80, "w_resid_%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
			i=0;
			do {
				if(datatag[i] == '_') datatag[i] = '-';
				i++;
			}	while (datatag[i] != '\0');
			fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6);
			if (j != 6 -1) fprintf(fpB, ", ");
		}
		fprintf(fpB, "\n");
		fprintf(fpB, "unset multiplot\n");
		fprintf(fpB, "set size 1, 1\n");
		fprintf(fpB, "set origin 0, 0\n");
		fprintf(fpB, "set term x11\n");
		fprintf(fpB, "replot\n");
		fclose(fpB);


	}	//  End if(fit_mode == 1 )	{		(plotdata)
	
	if(fit_mode == 2 || fit_mode ==0)	{	//Plot Guess
	
		//	Write Fits & Data to 3 .dat files
		G_max = 0;
		G_min = 10;
		
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				for(i=o.fit_start+n_integral/4; i<o.fit_stop; i++)	{
					if(G_array[j*veclen_helper+i] > G_max)	G_max = G_array[j*veclen_helper+i];
					if(G_array[j*veclen_helper+i] < G_min)	G_min = G_array[j*veclen_helper+i];
				}
			}
		}
		G_min = G_min - (G_max-G_min) * 0.025;
		G_max = G_max + (G_max-G_min) * 0.150;
	
		snprintf(resultsfile01, 80, "gtaul_Gs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
		if(fit_mode==0) snprintf(resultsfile01, 80, "gtaul_Gs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		for(i=0; i< veclen_helper; i++){
			fprintf(fpB, "%e", tau_G(i));
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", G_array[j*veclen_helper+i]);
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", G_stdev_array[j*veclen_helper+i]);
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", o.fits[j*veclen_helper+i]);
			for(j=0; j< 6; j++)	fprintf(fpB, "\t%e", o.w_resids[j*veclen_helper+i]);
			fprintf(fpB, "\n");
		}
		fclose(fpB);
		
		snprintf(resultsfile01, 80, "gtaul_GGGs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
		if(fit_mode==0)	snprintf(resultsfile01, 80, "gtaul_GGGs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		for(i=0; i< veclen_GGG; i++){
			for(w=0; w< veclen_GGG; w++){
				if(w>0 && i>0)	fprintf(fpB, "%9.9e\t%9.9e", tau_GGG(w), tau_GGG(i));
				if(w==0)		fprintf(fpB, "%9.9e\t%9.9e", tau_GGG(w)+0.000000001,tau_GGG(i));
				if(i==0)		fprintf(fpB, "%9.9e\t%9.9e", tau_GGG(w), tau_GGG(i)+0.000000001);
				for(j=0; j< 3; j++)	fprintf(fpB, "\t%e", GGG_array[j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				for(j=0; j< 3; j++)	fprintf(fpB, "\t%e", GGG_stdev_array[j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				for(j=0; j< 3; j++)	fprintf(fpB, "\t%e", o.fits[6*veclen_helper+j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				for(j=0; j< 3; j++)	fprintf(fpB, "\t%e", o.w_resids[6*veclen_helper+j*veclen_GGG*veclen_GGG+i*veclen_GGG+w]);
				fprintf(fpB, "\n");
			}
			fprintf(fpB, "\n");
		}
		fclose(fpB);

		num_selected_curves = 4;
		//	***	Write selected decays to a seperate file for plotting  *** //
		//	Array index order:	chirality, tau_2 slice, tau_1 vs tau_2 const, [tau_1,tau_2,G(tau_1,tau_2)] 
		snprintf(resultsfile01b, 80, "gtaul_select_GGGs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
		if(fit_mode==0)	snprintf(resultsfile01b, 80, "gtaul_select_GGGs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01b, "w");
		for(tau2 = 0; tau2 < veclen_GGG; tau2++)	{	//	tau_2 down vertical axis
			fprintf(fpB, "%i", tau2);
			tau = 0;	//	2x0x1 G(tau^0_1,X) | 1x2x0 G(X,tau^0_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[2*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[2*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			}
			tau = 0;	//	2x0x1 G(tau^0_1,X) | 1x2x0 G(X,tau^0_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[1*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[1*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			}
			tau = 0;	//	0x1x2 G(tau^1_1,X) | 2x0x1 G(X,tau^1_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[0*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[0*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			}
			tau = 0;	//	0x1x2 G(tau^1_1,X) | 2x0x1 G(X,tau^1_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[2*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[2*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			}
			tau = 0;	//	1x2x0 G(tau^0_2,X) | 0x1x2 G(X,tau^2_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[1*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau2), tau_GGG(tau), GGG_array[1*veclen_GGG*veclen_GGG+tau*veclen_GGG+tau2]);
			}
			tau = 0;	//	1x2x0 G(tau^2_1,X) | 0x1x2 G(X,tau^2_2)
			fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[0*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
				tau = o.fit_start_triplecorr+(o.fit_stop_triplecorr-o.fit_start_triplecorr-1)*i2/(num_selected_curves-1);
				fprintf(fpB, "\t%2.2e\t%2.2e\t%.6f", tau_GGG(tau), tau_GGG(tau2), GGG_array[0*veclen_GGG*veclen_GGG+tau2*veclen_GGG+tau]);
			}
			fprintf(fpB, "\n");
		}
		fclose(fpB);

		//	Gnuplot scripts to generate 4 graphs:	1x auto&cross, 3x triplecorr
		snprintf(resultsfile01, 80, "gtaul_Gs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
		snprintf(resultsfile02, 80, "gnuplot%s_Gs.%i.guess.txt",&filetag[128*current_entry],current_entry);
		if(fit_mode==0)	snprintf(resultsfile01, 80, "gtaul_Gs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
		if(fit_mode==0)	snprintf(resultsfile02, 80, "gnuplot%s_Gs.%i.fits.txt",&filetag[128*current_entry],current_entry);
		printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
		fpB = fopen(resultsfile02, "w");
		fprintf(fpB, "reset\n");
		fprintf(fpB, "set term postscript eps enhanced color size 9, 6.3 solid\n");
		snprintf(resultsfile02, 50, "gtaul_Gs%s.%i.guess.eps",&filetag[128*current_entry],current_entry);
		if(fit_mode==0)	snprintf(resultsfile02, 50, "gtaul_Gs%s.%i.fits.eps",&filetag[128*current_entry],current_entry);
		fprintf(fpB, "set output '%s'\n", resultsfile02);
		fprintf(fpB, "set logscale x\n");
		fprintf(fpB, "set xlabel \"\"\n");
		fprintf(fpB, "set format x \"\"\n");
		fprintf(fpB, "set  ylabel \"G(tau)\"\n");
		fprintf(fpB, "set palette defined ( 0 \"blue\", 1 \"goldenrod\", 2 \"red\", 3 \"green\", 4 \"orange\", 5 \"violet\" )\n");
		fprintf(fpB, "unset colorbox\n");
		
		if(o.fit_start == 0)	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(1)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		else	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(o.fit_start)*0.9, tau_G(o.fit_stop-1)*1.1, G_min, G_max);
		comma_yet = 0;
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				if(comma_yet==0)	{
					fprintf(fpB, "\"%s\" using 1:%i with points title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6);
					comma_yet = 1;
				}
				fprintf(fpB, ", \"%s\" using 1:%i with points title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) 6);
			}
		}
		for(j=0; j< 6; j++)	{
			if (o.Gfix[j]==1)	{
				snprintf(datatag, 80, "fit_%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[j]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, ", \"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*6, datatag, (float)j / (float) 6);
			}
		}
		fprintf(fpB, "\n");
		fprintf(fpB, "set term x11\n");
		fprintf(fpB, "replot\n");
		fprintf(fpB, "\n");
		fclose(fpB);
		
		G_max = 0;
		G_min = 10;
		
		for(j=0; j<3; j++)	{
			if (o.Gfix[j+6]==1)	{
				for(i=o.fit_start_triplecorr+n_GGG/2; i<o.fit_stop_triplecorr; i++)	{
					for(w=o.fit_start_triplecorr+n_GGG/2; w<o.fit_stop_triplecorr; w++)	{
						if(GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i] > G_max)	G_max = GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i];
						if(GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i] < G_min)	G_min = GGG_array[j*veclen_GGG*veclen_GGG+w*veclen_GGG+i];
					}
				}
			}
		}
		G_min = G_min - (G_max-G_min) * 0.025;
		G_max = G_max + (G_max-G_min) * 0.150;
		
		for(j=0; j< 3; j++)	{
			if (o.Gfix[j+6]==1)	{
				snprintf(resultsfile01, 80, "gtaul_GGGs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
				snprintf(resultsfile02, 80, "gnuplot%s_GGGs.%i.%i.guess.txt",&filetag[128*current_entry],current_entry,j);
				snprintf(resultsfile01b, 80, "gtaul_select_GGGs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
				if(fit_mode==0)	snprintf(resultsfile01, 80, "gtaul_GGGs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
				if(fit_mode==0)	snprintf(resultsfile02, 80, "gnuplot%s_GGGs.%i.%i.fits.txt",&filetag[128*current_entry],current_entry,j);
				if(fit_mode==0)	snprintf(resultsfile01b, 80, "gtaul_select_GGGs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
				char	colorlib[90];
				printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
				fpB = fopen(resultsfile02, "w");
				fprintf(fpB, "reset\n");
				fprintf(fpB, "Gmax = %f\n",G_max);
				fprintf(fpB, "Gmin = %f\n",G_min);
				fprintf(fpB, "set logscale x\n");
				fprintf(fpB, "set logscale y\n");
				fprintf(fpB, "set zlabel \"G^{0x1x2}(tau_1,tau_2)\" offset -1.75,1.5,0 rotate by 90 font \"Helvetica, 12\"\n");
				fprintf(fpB, "set ticslevel 0\n");
				fprintf(fpB, "set view 78,135\n");
				fprintf(fpB, "set title font \"Helvetica, 10\" \n");
				fprintf(fpB, "set palette defined (0 \"grey70\", 0.5 \"navy\", 1.0 \"orange\")\n");
				fprintf(fpB, "unset colorbox\n");
				fprintf(fpB, "set term postscript eps enhanced color\n");
				snprintf(resultsfile02, 50, "gtaul_GGGs%s.%i.%i.guess.eps",&filetag[128*current_entry],current_entry,j);
				fprintf(fpB, "set output '%s'\n", resultsfile02);
				fprintf(fpB, "set style line 10 linetype 1 linecolor rgb \"gray\"\n");
				fprintf(fpB, "set ztics font \"Helvetica, 10\"\n");
				fprintf(fpB, "unset hidden3d\n");
				fprintf(fpB, "set xtics offset 1.5,-0.5 font \"Helvetica, 10\"\n");
				fprintf(fpB, "set ytics offset -1.5,-0.5 font \"Helvetica, 10\"\n");
				fprintf(fpB, "set xlabel \"tau_1 (s)\" offset 1.5,-1.0 font \"Helvetica, 12\"\n");
				fprintf(fpB, "set ylabel \"tau_2 (s)\" offset -1.0,-1.0 font \"Helvetica, 12\"\n");
				if(j==1)	{fprintf(fpB, "set zlabel \"G^{1x2x0}(tau_1,tau_2)\" offset 3.0,1.5,0 rotate by 90 font \"Helvetica, 12\"\n"); fprintf(fpB, "set ztics format \"\"\n");}
				if(j==2)	 fprintf(fpB, "set zlabel \"G^{2x0x1}(tau_1,tau_2)\" offset 3.0,1.5,0 rotate by 90 font \"Helvetica, 12\"\n");
				if(o.fit_start_triplecorr==0)	fprintf(fpB, "splot [%2.2e:%2.2e][%2.2e:%2.2e][Gmin:Gmax] ", tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				else fprintf(fpB, "splot [%2.2e:%2.2e][%2.2e:%2.2e][Gmin:Gmax] ", tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:2:%i with lines title \"\" lw 0 lt palette frac 0", resultsfile01, j+3);
				fprintf(fpB, ",\"%s\" using 1:2:%i with lines title \"\" lw 0", resultsfile01, j+9);
				for(k=1; k<(num_selected_curves+1.0); k++)	{
					if(j==2)	snprintf(colorlib, 30, "blue");
					if(j==0)	snprintf(colorlib, 30, "orange");
					if(j==1)	snprintf(colorlib, 30, "red");
					fprintf(fpB, ", \"%s\" using %i:%i:%i with lines title \"\" lc rgb \"%s\" lw 2 lt 1", resultsfile01b, ((j+1)%3)*6*(num_selected_curves+1)+k*3+2,((j+1)%3)*6*(num_selected_curves+1)+k*3+3,((j+1)%3)*6*(num_selected_curves+1)+k*3+4,colorlib);
					if(j==1)	snprintf(colorlib, 30, "blue");
					if(j==2)	snprintf(colorlib, 30, "orange");
					if(j==0)	snprintf(colorlib, 30, "red");
					fprintf(fpB, ", \"%s\" using %i:%i:%i with lines title \"\" lc rgb \"%s\" lw 2 lt 2", resultsfile01b, (((j+2)%3)*6+3)*(num_selected_curves+1)+k*3+2,(((j+2)%3)*6+3)*(num_selected_curves+1)+k*3+3,(((j+2)%3)*6+3)*(num_selected_curves+1)+k*3+4,colorlib);
				}
				fprintf(fpB, "\n");
				fprintf(fpB, "set term x11\n");
				fprintf(fpB, "replot\n\n");
				fclose(fpB);
			}
		}

		for(j=0; j< 3; j++)	{
			if (o.Gfix[j+6]==1)	{
				snprintf(resultsfile01, 80, "gtaul_GGGs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
				snprintf(resultsfile02, 80, "gnuplot%s_GGGs.%i.%i.resids.txt",&filetag[128*current_entry],current_entry,j);
				snprintf(resultsfile01b, 80, "gtaul_select_GGGs%s.%i.guess.dat",&filetag[128*current_entry],current_entry);
				if(fit_mode==0)	snprintf(resultsfile01, 80, "gtaul_GGGs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
				if(fit_mode==0)	snprintf(resultsfile02, 80, "gnuplot%s_GGGs.%i.%i.resids.txt",&filetag[128*current_entry],current_entry,j);
				if(fit_mode==0)	snprintf(resultsfile01b, 80, "gtaul_select_GGGs%s.%i.fits.dat",&filetag[128*current_entry],current_entry);
				printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
				fpB = fopen(resultsfile02, "w");
				fprintf(fpB, "reset\n");
				fprintf(fpB, "Gmax = %f\n",G_max);
				fprintf(fpB, "Gmin = %f\n",G_min);
				fprintf(fpB, "set logscale x\n");
				fprintf(fpB, "set logscale y\n");
				fprintf(fpB, "set zlabel \"G^{0x1x2}(tau_1,tau_2)\" offset -1.75,1.5,0 rotate by 90 font \"Helvetica, 12\"\n");
				fprintf(fpB, "set ticslevel 0\n");
				fprintf(fpB, "set view 78,135\n");
				fprintf(fpB, "set title font \"Helvetica, 10\" \n");
				fprintf(fpB, "set palette defined (0 \"grey70\", 0.5 \"navy\", 1.0 \"orange\")\n");
				fprintf(fpB, "unset colorbox\n");
				fprintf(fpB, "set term postscript eps enhanced color\n");
				snprintf(resultsfile02, 50, "gtaul_GGGs%s.%i.%i.guess.eps",&filetag[128*current_entry],current_entry,j);
				fprintf(fpB, "set output '%s'\n", resultsfile02);
				fprintf(fpB, "set style line 10 linetype 1 linecolor rgb \"gray\"\n");
				fprintf(fpB, "set ztics font \"Helvetica, 10\"\n");
				fprintf(fpB, "unset hidden3d\n");
				fprintf(fpB, "set xtics offset 1.5,-0.5 font \"Helvetica, 10\"\n");
				fprintf(fpB, "set ytics offset -1.5,-0.5 font \"Helvetica, 10\"\n");
				fprintf(fpB, "set xlabel \"tau_1 (s)\" offset 1.5,-1.0 font \"Helvetica, 12\"\n");
				fprintf(fpB, "set ylabel \"tau_2 (s)\" offset -1.0,-1.0 font \"Helvetica, 12\"\n");
				if(o.fit_start_triplecorr==0)	fprintf(fpB, "splot [%2.2e:%2.2e][%2.2e:%2.2e][] ", tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, tau_GGG(1)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				else fprintf(fpB, "splot [%2.2e:%2.2e][%2.2e:%2.2e][] ", tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1, tau_GGG(o.fit_start_triplecorr)*0.9, tau_GGG(o.fit_stop_triplecorr-1)*1.1);
				snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
				i=0;
				do {
					if(datatag[i] == '_') datatag[i] = '-';
					i++;
				}	while (datatag[i] != '\0');
				fprintf(fpB, "\"%s\" using 1:2:%i with lines title \"\" lw 0 lt palette frac 0", resultsfile01, j+12);
				fprintf(fpB, "\n");
				fprintf(fpB, "set term x11\n");
				fprintf(fpB, "replot\n\n");
				fclose(fpB);
			}
		}

		float toggle = 0.1;
		snprintf(resultsfile01, 80, "gnuplot%s.%i.allresids.dat",&filetag[128*current_entry],current_entry);
		fpB = fopen(resultsfile01, "w");
		i=0;
		for(j=0; j< Gcount; j++)	{
			for(i2=0; i2< o.fit_stop-o.fit_start; i2++)	{
				fprintf(fpB, "%e\t%e\t%e\t%e\n", (float)i, pow(gsl_vector_get(s->f,i),2.0), 1.0, toggle);
				i++;
			}
			if(toggle < 100)	toggle = 100;
			else toggle = 0.1;
		}
		for(j=0; j< GGGcount; j++)	{
			for(i2=0; i2< (o.fit_stop_triplecorr-o.fit_start_triplecorr)*(o.fit_stop_triplecorr-o.fit_start_triplecorr); i2++)	{
				fprintf(fpB, "%e\t%e\t%e\t%e\n", (float)i, pow(gsl_vector_get(s->f,i),2.0), 1.0, toggle);
				i++;
			}
			if(toggle < 100)	toggle = 100;
			else toggle = 0.1;
		}
		for(j=0; j<Icount; j++)	{
			fprintf(fpB, "%e\t%e\t%e\t%e\n", (float)i+((j+1)*100.0), pow(gsl_vector_get(s->f,i),2.0), 1.0, toggle);
			i++;
		}
		fclose(fpB);
						
		snprintf(resultsfile02, 80, "gnuplot%s.%i.allresids.txt",&filetag[128*current_entry],current_entry);
		printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
		fpB = fopen(resultsfile02, "w");
		fprintf(fpB, "reset\n");
		fprintf(fpB, "Gmax = %f\n",G_max);
		fprintf(fpB, "Gmin = %f\n",G_min);
		fprintf(fpB, "set logscale y\n");
		fprintf(fpB, "set title font \"Helvetica, 10\" \n");
		fprintf(fpB, "unset colorbox\n");
		fprintf(fpB, "set term postscript eps enhanced color\n");
		snprintf(resultsfile02, 50, "gnuplot%s.%i.%i.allresids.eps",&filetag[128*current_entry],current_entry,j);
		fprintf(fpB, "set output '%s'\n", resultsfile02);
		fprintf(fpB, "set style line 10 linetype 1 linecolor rgb \"gray\"\n");
		fprintf(fpB, "plot [][0.1:100]");
		snprintf(datatag, 80, "%3.3s%s_%.2f", &corr_mode[5*j], &filetag[128*current_entry], f1[current_entry]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:2 with points title \"All Weighted Resids\" lw 0 lt palette frac 0, \"%s\" using 1:3 with lines title \"\" lc rgb \"blue\", \"%s\" using 1:4 with lines title \"\" lc rgb \"red\" ", resultsfile01, resultsfile01, resultsfile01);
		fprintf(fpB, "\n");
		fprintf(fpB, "set term x11\n");
		fprintf(fpB, "replot\n\n");
		fclose(fpB);


#ifdef useplplot
		 // ***	Make shaded plots with plplot.  Code largely adapted from plplot 3dplot example ***	//
		 snprintf(resultsfile02, 80, "Shaded%s.%i.ps",&filetag[128*current_entry],current_entry);
		 PLFLT    *x_axis, *y_axis, **z;
		 PLFLT    ii[2], h[2], l[2], s[2]; 
		 PLFLT    zmin, zmax;
		 if (o.Gfix[6]==1 || o.Gfix[7]==1 || o.Gfix[8]==1) {
		 
		 plMergeOpts( options, "x08c options", NULL );
		 
		 plsdev("psc");		// Output = ps, in colour
		 plsfnam(resultsfile02);	// Output file name
		 plscol0(15, 0, 0, 0);
		 plinit();
		 
		 // Data structures 
		 x_axis = (PLFLT *) calloc( (o.fit_stop_triplecorr-o.fit_start_triplecorr), sizeof ( PLFLT ) );
		 y_axis = (PLFLT *) calloc( (o.fit_stop_triplecorr-o.fit_start_triplecorr), sizeof ( PLFLT ) );
		 plAlloc2dGrid( &z, o.fit_stop_triplecorr-o.fit_start_triplecorr, o.fit_stop_triplecorr-o.fit_start_triplecorr );

         // Calculate the plot fxn 
		 printf("bound = %i\n",o.fit_stop_triplecorr-o.fit_start_triplecorr);
		 for ( i = o.fit_start_triplecorr; i < o.fit_stop_triplecorr; i++ )	x_axis[i-o.fit_start_triplecorr] = log10(tau_GGG(i)); 
		 for ( i = o.fit_start_triplecorr; i < o.fit_stop_triplecorr; i++ )	y_axis[i-o.fit_start_triplecorr] = log10(tau_GGG(i));
		 
		 
		 for(j=0; j<10; j++)	{ // G0x1x2, G1x2x0, G2x0x1...
		 for ( i = 0; i < o.fit_stop_triplecorr-o.fit_start_triplecorr; i++ )	{
		 for ( w = 0; w < o.fit_stop_triplecorr-o.fit_start_triplecorr; w++ )	{
		 if (j<3)	z[i][w] = GGG_array[(j%3)*veclen_GGG*veclen_GGG+(i+o.fit_start_triplecorr)*veclen_GGG+w+o.fit_start_triplecorr];
		 else if(j < 6) z[i][w] = o.fits[6*veclen_helper+(j%3)*veclen_GGG*veclen_GGG+(i+o.fit_start_triplecorr)*veclen_GGG+w+o.fit_start_triplecorr];
		 else if(j ==9) z[i][w] = (GGG_array[(0)*veclen_GGG*veclen_GGG+(i+o.fit_start_triplecorr)*veclen_GGG+w+o.fit_start_triplecorr]+GGG_array[(1)*veclen_GGG*veclen_GGG+(i+o.fit_start_triplecorr)*veclen_GGG+w+o.fit_start_triplecorr]+GGG_array[(2)*veclen_GGG*veclen_GGG+(i+o.fit_start_triplecorr)*veclen_GGG+w+o.fit_start_triplecorr])/3.0;
		 else z[i][w] = o.w_resids[6*veclen_helper+(j%3)*veclen_GGG*veclen_GGG+(i+o.fit_start_triplecorr)*veclen_GGG+w+o.fit_start_triplecorr];
		 if(isnan(z[i][w])) z[i][w] = -1.0;
		 }		}
		 
		 plMinMax2dGrid((const PLFLT **) z, o.fit_stop_triplecorr-o.fit_start_triplecorr, o.fit_stop_triplecorr-o.fit_start_triplecorr, &zmax, &zmin );
		 zmax = zmax+0.1*(zmax-zmin);
		 if(j<6 && Plot_Z_max_fit == 1)		zmax = Plot_Z_max;
		 if(j<6 && Plot_Z_min_fit == 1)		zmin = Plot_Z_min;
		 if(j==9 && Plot_Z_max_fit == 1)          zmax = Plot_Z_max;
		 if(j==9 && Plot_Z_min_fit == 1)          zmin = Plot_Z_min;
		 if(zmax==zmin) {printf("zmax=zmin, manually adjusting to inputfile values\n"); zmax = Plot_Z_max; zmin = Plot_Z_min;}
		 
		 if(j==9) { //ReNorm. average plot by 10^something so that the shading function looks better
		 tempa = pow(10.0,ceil(-1.0*log10(zmax)));
		 printf("PLPlot Scaling factor for av. G(tau,tau) = %f\n", tempa);
		 for ( i = 0; i < o.fit_stop_triplecorr-o.fit_start_triplecorr; i++ )    {
		 for ( w = 0; w < o.fit_stop_triplecorr-o.fit_start_triplecorr; w++ )    {
		 z[i][w] = z[i][w] * tempa;
		 }	}	
		 zmax = zmax * tempa;
		 zmin = zmin * tempa;
		 }
		 
		 pllightsource( -10.0, 10.0, 10.0 );
		 pladv( 0 );
		 plvpor( 0.0, 1.0, 0.0, 0.9 );
		 plwind( -1.0, 1.0, -0.9, 1.1 );
		 plcol0( 15 );
		 snprintf(resultsfile02, 80, "G%1ix%1ix%1i%s", j%3, (j+1)%3,(j+2)%3, &filetag[128*current_entry]);
		 if(j>2) snprintf(resultsfile02, 80, "Fits");
		 if(j>5) snprintf(resultsfile02, 80, "Weighted Residuals");
		 if(j>8) snprintf(resultsfile02, 80, "G_Averaged%s", &filetag[128*current_entry]);
		 plcol0( 15 );
		 plmtex( "t", 1.0, 0.5, 0.5, resultsfile02 );
		 printf("%f\t%f -- %f\t%f\t--%f\t%f\n",2*x_axis[0]-x_axis[1], 2.0*x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-1]-x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-2], 
		 2*x_axis[0]-x_axis[1], 2.0*x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-1]-x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-2], zmin, zmax);
		 plw3d( 1.0, 1.0, 1.0, 2*x_axis[0]-x_axis[1], 2.0*x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-1]-x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-2], 
		 2*x_axis[0]-x_axis[1], 2.0*x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-1]-x_axis[o.fit_stop_triplecorr-o.fit_start_triplecorr-2], zmin, zmax,	20.0, 240 );
		 
		 if(j<6 || j==9) plbox3( "bnsltu", "tau_2", 0.0, 0, "bnsltu", "tau_1", 0.0, 0, "bcdmnstuv", "G(tau1,tau2)", 0.0, 0 );
		 else    plbox3( "bnsltu", "tau_1", 0.0, 0, "bnsltu", "tau_1", 0.0, 0, "bcdmnstuv", "W. Resids", 0.0, 0 );
		 
		 ii[0] = 0.0;    // left boundary 
		 ii[1] = 1.0;
		 h[0] = 0.0;     // hue 
		 h[1] = 0.0;     
		 l[0] = 0.4;     // lightness
		 l[1] = 1.0;     
		 s[0] = 0.0;     // minimum saturation 
		 s[1] = 0.0;      
		 plscmap1n( 256 );
		 c_plscmap1l( 0, 2, ii, h, l, s, NULL );
		 plcol0( 15 );
		 plfsurf3d( x_axis, y_axis, plf2ops_c(), (PLPointer) z, o.fit_stop_triplecorr-o.fit_start_triplecorr, o.fit_stop_triplecorr-o.fit_start_triplecorr, FACETED, NULL, 0 );
		 }
		 
		 // Free memory
		 free( (void *) x_axis );
		 free( (void *) y_axis );
		 plFree2dGrid( z, o.fit_stop_triplecorr-o.fit_start_triplecorr, o.fit_stop_triplecorr-o.fit_start_triplecorr );
		 plend();
		 }	// End plplot 
#endif
		

// ***** Calculate signal / noise ratio for 10us - 100us  ***** //
		int	sntau_min = 0;
		int	sntau_max = 0;
		for(j=0; j<veclen_GGG; j++)	{
			if(tau_GGG(j) > 0.00001 && sntau_min == 0)	sntau_min = j;
			if(tau_GGG(j) < 0.0001)	sntau_max = j;
		}
		double	snr = 0.0;
		int		pt_counter = 0;
		for(j=0; j<3; j++)	{
			for ( i = sntau_min; i < sntau_max; i++ )	{
				for ( w = sntau_min; w < sntau_max; w++ )	{
					snr += GGG_array[(j%3)*veclen_GGG*veclen_GGG+(i)*veclen_GGG+w]/GGG_stdev_array[(j%3)*veclen_GGG*veclen_GGG+(i)*veclen_GGG+w];
					pt_counter++;
		}	}	}
		snr = snr / (float) pt_counter;
		
		// repeat for double-correlations
		sntau_min = 0;
		sntau_max = 0;
		for(j=0; j<veclen; j++)	{
			if(tau_G(j) > 0.00001 && sntau_min == 0)	sntau_min = j;
			if(tau_G(j) < 0.0001)						sntau_max = j;
		}
		double	snr_d[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		pt_counter = 0;
		for(j=0; j<6; j++)	{
			for ( i = sntau_min; i < sntau_max; i++ )	{
				snr_d[j] += (G_array[j*veclen_helper+i]-1.0)/G_stdev_array[j*veclen_helper+i];
				pt_counter++;
		}	}
		pt_counter = pt_counter / 6;
		for(j=0; j< 6; j++)		snr_d[j] = snr_d[j] / (float) pt_counter;
		
		printf("SNR%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", &filetag[128*current_entry], snr, snr_d[0], snr_d[1], snr_d[2], snr_d[3], snr_d[4], snr_d[5]);

	}		// End if(fit_mode == 1)
	}	//	End main loop: read, fit, write
	
	
// ****************** Write Parameters to Disk **************** //
	if(fit_mode == 0 )	{	//Fit Mode
		snprintf(resultsfile01, 50, "gtaul%s.fit_params.dat",&filetag[128*current_entry]);
		fpB = fopen(resultsfile01, "w");
		
		fprintf(fpB,"N0%s\tN1%s\tN2%s\tN3%s\tN4%s\tN5%s\tN6%s\tN7%s\tN8%s\tN9%s\ttauD0%s\ttauD1%s\ttauD2%s\ttauD3%s\ttauD4%s\ttauD5%s\ttauD6%s\ttauD7%s\ttauD%s\ttauD9%s\tSDN0%s\tSDN1%s\tSDN2%s\tSDN3%s\tSDN4%s\tSDN5%s\tSDN6%s\tSDN7%s\tSDN8%s\tSDN9%s\tTim%s\tChiSq%s\tInt0%s\tInt1%s\tInt2%s\tIntT%s\tFitTerm%s\tq0%s\tq1%s\tq2%s\tGammaR%s\n", argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2], argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2], argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2], argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2]);
	
		for (current_entry=0; current_entry<num_entries; current_entry++)	{
			fprintf(fpB, "%e", ptw[current_entry*nptw + 0]);
			for(j=1; j<nptw; j++)	
				fprintf(fpB, "\t%e", ptw[current_entry*nptw + j]);
			fprintf(fpB, "\n");
		}
		
	
		fclose(fpB);
	} // End if(fit_mode == 0 )	{


/******************** Free memory ********************/

    free(tauarray);
    free(time_array);
    free(o.m);
	free(o.m_no_triplet);
	free(o.mGGG);
    free(o.fits);
    gsl_vector_free (x);
    gsl_vector_free (fff);
	gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);  
    free(o.y);
    free(o.sigma);
    free(G_array);
    free(G_stdev_array);
    free(G_array_from_GGG);
    free(G_stdev_array_from_GGG);
    free(GGG_array);
    free(GGG_stdev_array);
    free(intensity0_array);
    free(intensity0_stdev_array);
    free(intensity1_array);
    free(intensity1_stdev_array);
    free(intensity2_array);
    free(intensity2_stdev_array);
	free(intensity0_array_GGG);
	free(intensity1_array_GGG);
	free(intensity2_array_GGG);
    free(pmask);
    free(xtemp);
    free(fits);
    free(o.w_resids);
	free(ptw);

    printf("Program completed.\n");
    return 0;
}
