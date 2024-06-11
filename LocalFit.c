/*
 *  LocalFit.c
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
 *  LocalFit.c
 *  Self-contained code to generate F3CS_LocalFit
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

#define num_curves_to_consider 128 
#define	fcs_model_params	19

#define	N	40
#define verbose 0

int nn = n;
#undef n
int n;

int	global_one = 1;

//  F3CS_LocalFit - a program to globally fit a series of FCS curves
//
//  2/14/09- WKR


//	****	Data structure to pass to the functions which generate fit functions for FCS data	*** //
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
    
    double * m;
};


struct data { 
    size_t nFIT;
    size_t nfull;
    size_t p;
    double * y;
    double * sigma;
    int	* p_mask;
    double * x_array;
    int	num_entries;
    int	fit_start;
    int	fit_stop;
};

//	***		Calculate a pseudo-logarithmic series of tau values		*** //
double  tau_G( int i){
	if (i<0)		return 0.0;
	
    double t = 1.0 / 1250000.0;
    int	n_tau = n;
    
    int j = (i-n_tau-1)/n_tau;
    if (i<2*n_tau)   	return t * i;

    double tempa = pow(2.0, j);
    
    t = t*(tempa*(i-j*n_tau)+1);
    
    return t;
}

//	***	Fit function for double correlation FCS data (both auto- and cross-correlations)	***	//
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
    double	tempt  = 0.0;
    double	tempk  = 0.0;
    double	tempk2 = a->k2/a->int0;
    double	tempk3 = a->k3/a->int0;
    double	tempk4 = a->k4/a->int0;

    for(i=0; i < veclen; i++)  {
        t = tau_G(i) * 1000000.0;
        tempt = 1.0/tau_G(i);
        tempt = 1.0 * tempt / tau_G(i);
        tempk =  tempk2 * tempt;
        tempt = 1.0 * tempt / tau_G(i);
        tempk += tempk3 * tempt;
        tempt = 1.0 * tempt / tau_G(i);
        tempk += tempk4 * tempt;
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
        a->m[i] = (1.0/a->N_t)*(a->theta*((1.0-a->frac_gamma_1)*tempa1*tempb1 + a->frac_gamma_1*tempc1*tempd1)
			+ (1.0-a->theta)*((1.0-a->frac_gamma_2)*tempa2*tempb2 + a->frac_gamma_2*tempc2*tempd2))
			* (1+tempg1*tempe1+tempg2*tempe2)
			+ a->g_inf;
		if(tau_G(i) < 0.0001)	a->m[i] += tempk;
    }  

    return;
}

//	***	Intermediate function called by the solver routine.  This function collects the current fit parameter guesses together, then passes them to m_tau() function to calculate the correlation fit functions ***	//
int gtau_f (const gsl_vector * x, void *data, gsl_vector * f)   {
	struct data * a = data;
	struct m_pass_struct o;
	o.m = malloc(veclen*sizeof(double));
	
	size_t i,j;
	double	tempa;
	
	//	Vary the fit parameters requested by the user
	for(i=0; i<fcs_model_params*a->num_entries; i++)	
		if (a->p_mask[i] != -1) a->x_array[i] = gsl_vector_get(x, a->p_mask[i]);
	
    for(j=0; j<a->num_entries; j++)	{
		i=0;
		o.N_t = a->x_array[(i++)+j*fcs_model_params];
		o.theta = a->x_array[(i++)+j*fcs_model_params];
		o.tauD_1 = a->x_array[(i++)+j*fcs_model_params];
		o.tauD_2 = a->x_array[(i++)+j*fcs_model_params];
		o.omega = a->x_array[(i++)+j*fcs_model_params];
		o.g_inf = a->x_array[(i++)+j*fcs_model_params];
		o.tauf_1 = a->x_array[(i++)+j*fcs_model_params];
		o.tauf_2 = a->x_array[(i++)+j*fcs_model_params];
		o.T0_1 = a->x_array[(i++)+j*fcs_model_params];
		o.T0_2 = a->x_array[(i++)+j*fcs_model_params];
		o.gamma_1 = a->x_array[(i++)+j*fcs_model_params];
		o.gamma_2 = a->x_array[(i++)+j*fcs_model_params];
		o.frac_gamma_1 = a->x_array[(i++)+j*fcs_model_params];
		o.frac_gamma_2 = a->x_array[(i++)+j*fcs_model_params];
		o.k2 = a->x_array[(i++)+j*fcs_model_params];
		o.k3 = a->x_array[(i++)+j*fcs_model_params];
		o.k4 = a->x_array[(i++)+j*fcs_model_params];
		o.int0 = a->x_array[(i++)+j*fcs_model_params];

		m_tau(&o);
    
		for(i=0; i<(a->fit_stop-a->fit_start); i++)	{
			tempa = o.m[i+(a->fit_start)] - a->y[(i+(a->fit_start)+j*veclen)];
			tempa = tempa / a->sigma[(i+(a->fit_start)+j*veclen)];
			gsl_vector_set(f,i+j*(a->fit_stop-a->fit_start),tempa);
		}
    } 
        
    free(o.m);
    return GSL_SUCCESS;
}

//	**** Calcualates the Jacobian matrix J(i,j) = dfi / dxj, ****	//
int gtau_df (const gsl_vector * x, void *data, gsl_matrix * J)  {
	struct data * a = data;
	struct m_pass_struct o;
	o.m = malloc(veclen*sizeof(double));
	double	x_pl = 1.00001;
	double  x_mi = 0.99999;
	
	size_t p = a->p;
	size_t i,j,k;
	
	gsl_vector * plus = gsl_vector_alloc (veclen*(a->num_entries));
	gsl_vector * minus = gsl_vector_alloc (veclen*(a->num_entries));
	gsl_vector * x_tweaked = gsl_vector_alloc (p);
	
	for(k=0; k<p; k++)	{
		for (i=0; i<p; i++)		gsl_vector_set (x_tweaked, i, gsl_vector_get(x, i));
		gsl_vector_set (x_tweaked, k, gsl_vector_get(x, k)*x_pl);
		for(i=0; i<fcs_model_params*a->num_entries; i++)	
			if (a->p_mask[i] != -1) a->x_array[i] = gsl_vector_get(x_tweaked, a->p_mask[i]);
		
		for(j=0; j<a->num_entries; j++)	{
			i=0;
			o.N_t = a->x_array[(i++)+j*fcs_model_params];
			o.theta = a->x_array[(i++)+j*fcs_model_params];
			o.tauD_1 = a->x_array[(i++)+j*fcs_model_params];
			o.tauD_2 = a->x_array[(i++)+j*fcs_model_params];
			o.omega = a->x_array[(i++)+j*fcs_model_params];
			o.g_inf = a->x_array[(i++)+j*fcs_model_params];
			o.tauf_1 = a->x_array[(i++)+j*fcs_model_params];
			o.tauf_2 = a->x_array[(i++)+j*fcs_model_params];
			o.T0_1 = a->x_array[(i++)+j*fcs_model_params];
			o.T0_2 = a->x_array[(i++)+j*fcs_model_params];
			o.gamma_1 = a->x_array[(i++)+j*fcs_model_params];
			o.gamma_2 = a->x_array[(i++)+j*fcs_model_params];
			o.frac_gamma_1 = a->x_array[(i++)+j*fcs_model_params];
			o.frac_gamma_2 = a->x_array[(i++)+j*fcs_model_params];
			o.k2 = a->x_array[(i++)+j*fcs_model_params];
			o.k3 = a->x_array[(i++)+j*fcs_model_params];
			o.k4 = a->x_array[(i++)+j*fcs_model_params];
			o.int0 = a->x_array[(i++)+j*fcs_model_params];

			m_tau(&o);
			for (i=0; i<(a->fit_stop-a->fit_start); i++)	gsl_vector_set (plus, i+j*(a->fit_stop-a->fit_start), o.m[i+a->fit_start]);
		
		} 
		
		for (i=0; i<p; i++)		gsl_vector_set (x_tweaked, i, gsl_vector_get(x, i));
		gsl_vector_set (x_tweaked, k, gsl_vector_get(x, k)*x_mi);
		for(i=0; i<fcs_model_params*a->num_entries; i++)	
			if (a->p_mask[i] != -1) a->x_array[i] = gsl_vector_get(x_tweaked, a->p_mask[i]);
	
		for(j=0; j<a->num_entries; j++)	{
			i=0;
			o.N_t = a->x_array[(i++)+j*fcs_model_params];
			o.theta = a->x_array[(i++)+j*fcs_model_params];
			o.tauD_1 = a->x_array[(i++)+j*fcs_model_params];
			o.tauD_2 = a->x_array[(i++)+j*fcs_model_params];
			o.omega = a->x_array[(i++)+j*fcs_model_params];
			o.g_inf = a->x_array[(i++)+j*fcs_model_params];
			o.tauf_1 = a->x_array[(i++)+j*fcs_model_params];
			o.tauf_2 = a->x_array[(i++)+j*fcs_model_params];
			o.T0_1 = a->x_array[(i++)+j*fcs_model_params];
			o.T0_2 = a->x_array[(i++)+j*fcs_model_params];
			o.gamma_1 = a->x_array[(i++)+j*fcs_model_params];
			o.gamma_2 = a->x_array[(i++)+j*fcs_model_params];
			o.frac_gamma_1 = a->x_array[(i++)+j*fcs_model_params];
			o.frac_gamma_2 = a->x_array[(i++)+j*fcs_model_params];
			o.k2 = a->x_array[(i++)+j*fcs_model_params];
			o.k3 = a->x_array[(i++)+j*fcs_model_params];
			o.k4 = a->x_array[(i++)+j*fcs_model_params];
			o.int0 = a->x_array[(i++)+j*fcs_model_params];

			m_tau(&o);
			for (i=0; i<(a->fit_stop-a->fit_start); i++)	gsl_vector_set (minus, i+j*(a->fit_stop-a->fit_start), o.m[i+a->fit_start]);
		}
		for(j=0; j<a->num_entries; j++)
			for(i=0; i<(a->fit_stop-a->fit_start); i++)	gsl_matrix_set(J,i+j*(a->fit_stop-a->fit_start),k,(gsl_vector_get(plus, i+j*(a->fit_stop-a->fit_start)) - gsl_vector_get(minus,i+j*(a->fit_stop-a->fit_start)))/(a->sigma[i+j*(a->fit_stop-a->fit_start)+a->fit_start]*(x_pl-x_mi)*gsl_vector_get(x,k)));
	}

    free(o.m);
    gsl_vector_free(plus);
    gsl_vector_free(minus);
    gsl_vector_free(x_tweaked);
    return GSL_SUCCESS;

}

//	***	Intermediate function for compatibility w/ the GSL solver	***	//
int gtau_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
    gtau_f (x, data, f);
    gtau_df (x, data, J);

    return GSL_SUCCESS;
}

//	***		Main routine	***	//
int main(int argc, char *argv[])    {
	FILE * fp;
	int i,j,ff;
	n = nn;
	
	//	Usage
	if(argc != 3)  {
        printf("Usage: F3CS_LocalFit inputfile.txt tag\ne.g.> F3CS_LocalFit localfit1.txt sampleA\n");
        return 1;
    }
    
    fp = fopen(argv[1],"r");
    if (fp == NULL) {
        printf("Input file %s not found.  Exiting.\n", argv[1]);
        return 1;
    }
    
    
    //  Read the input file to obtain initial Model Parameters and establish which fit parameters to vary
    float  N_t = -10.0;
    float  theta = -10.0;
    float  tauD_1 = -10.0;
    float  tauD_2 = -10.0;
    float  omega = -10.0;
    float  g_inf = -10.0;
    float  tauf_1 = -10.0;
    float  tauf_2 = -10.0;
    float  T0_1 = -10.0;
    float  T0_2 = -10.0;
    float  gamma_1 = -10.0;
    float  gamma_2 = -10.0;
    float  frac_gamma_1 = -10.0;
    float  frac_gamma_2 = -10.0;
    float  k2[3]={0.0,0.0,0.0};
    float  k3[3]={0.0,0.0,0.0};
    float  k4[3]={0.0,0.0,0.0};
    
	//	fix_X	is a flag for whether to fit the model parameter X
    int  fix_N_t = -1;
    int  fix_theta = -1;
    int  fix_tauD_1 = -1;
    int  fix_tauD_2 = -1;
    int  fix_omega = -1;
    int  fix_g_inf = -1;
    int  fix_tauf_1 = -1;
    int  fix_tauf_2 = -1;
    int  fix_T0_1 = -1;
    int  fix_T0_2 = -1;
    int  fix_gamma_1 = -1;
    int  fix_gamma_2 = -1;
    int  fix_frac_gamma_1 = -1;
    int  fix_frac_gamma_2 = -1;
    
    char    line[512];
    
    int ret = 0;
    if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
	if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
	if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
	if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};

    int linenum = 5;
    ret = fscanf(fp,"N_t=%f,%i\n", &N_t, &fix_N_t);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    
    ret = fscanf(fp,"theta=%f,%i\n", &theta, &fix_theta);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;

    ret = fscanf(fp,"tauD_1=%f,%i\n", &tauD_1, &fix_tauD_1);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    ret = fscanf(fp,"tauD_2=%f,%i\n", &tauD_2, &fix_tauD_2);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    
    ret = fscanf(fp,"omega=%f,%i\n", &omega, &fix_omega);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    
    ret = fscanf(fp,"g_inf=%f,%i\n", &g_inf, &fix_g_inf);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    
    ret = fscanf(fp,"tauf_1=%f,%i\n", &tauf_1, &fix_tauf_1);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    ret = fscanf(fp,"tauf_2=%f,%i\n", &tauf_2, &fix_tauf_2);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    
    ret = fscanf(fp,"T0_1=%f,%i\n", &T0_1, &fix_T0_1);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    ret = fscanf(fp,"T0_2=%f,%i\n", &T0_2, &fix_T0_2);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
        
    ret = fscanf(fp,"gamma_1=%f,%i\n", &gamma_1, &fix_gamma_1);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    ret = fscanf(fp,"gamma_2=%f,%i\n", &gamma_2, &fix_gamma_2);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    
    ret = fscanf(fp,"frac_gamma_1=%f,%i\n", &frac_gamma_1, &fix_frac_gamma_1);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    ret = fscanf(fp,"frac_gamma_2=%f,%i\n", &frac_gamma_2, &fix_frac_gamma_2);
    if (ret != 2 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    
	if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
    linenum++;
    if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
    linenum++;

    int num_entries = 0;
    ret = fscanf(fp,"num_entries=%i\n", &num_entries);
    if (ret != 1 )  {printf("File format not valid, line %i.  Exiting\n", linenum); fclose(fp); return 1;}
    linenum++;
    if (num_entries > num_curves_to_consider)  {printf("recompile with num_curves_to_consider > 128.  Exiting.\n"); return 1; }

    char    filetag[128*num_curves_to_consider];
    float  f1[num_curves_to_consider];
    char    corr_mode[128*num_curves_to_consider];
    int     corr_int[num_curves_to_consider];
    for(i=0; i< num_curves_to_consider; i++)    corr_int[i] = -9;
    for(i=0; i<num_entries; i++)    {
        ret = fscanf(fp,"tag=%s timepoint = %f corr=%s\n", &filetag[128*i], &f1[i], &corr_mode[128*i]);
        if (ret != 3 )  {printf("File format not valid, Data entry %i (line %i).  Exiting\n", i+1, linenum); fclose(fp); return 1;}    
        if  (corr_mode[128*i] == '0' && corr_mode[128*i+1] == 'x' && corr_mode[128*i+2] == '0') corr_int[i] = 0;
        if  (corr_mode[128*i] == '1' && corr_mode[128*i+1] == 'x' && corr_mode[128*i+2] == '1') corr_int[i] = 1;
        if  (corr_mode[128*i] == '2' && corr_mode[128*i+1] == 'x' && corr_mode[128*i+2] == '2') corr_int[i] = 2;
        if  (corr_mode[128*i] == '0' && corr_mode[128*i+1] == 'x' && corr_mode[128*i+2] == '1') corr_int[i] = 3;
        if  (corr_mode[128*i] == '1' && corr_mode[128*i+1] == 'x' && corr_mode[128*i+2] == '2') corr_int[i] = 4;
        if  (corr_mode[128*i] == '2' && corr_mode[128*i+1] == 'x' && corr_mode[128*i+2] == '0') corr_int[i] = 5;
        if  (corr_int[i] == -9) { printf("Correlation Mode for line (%i) is invalid.  Exiting\n", linenum); return 1;}
        if(verbose) printf("Tag = %s, f1 = %f, corr_mode = %s\n", &filetag[128*i], f1[i], &corr_mode[128*i]);
        linenum++;
    }
    
	if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
    linenum++;
    if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
    linenum++;
    
	int	fit_start = 0;
	int	fit_stop = 0;
	ret = fscanf(fp,"fit_range_start = %i\nfit_range_stop = %i\n", &fit_start, &fit_stop);
	if (ret != 2 )  {printf("File format not valid, fit_start, fit_stop. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}    
	if(fit_stop > 0)	{printf("fit_stop should be negative.  Exiting\n"); return 1;}
	fit_stop = veclen + fit_stop;
	if(verbose) printf("Data will be fit from i = (%i,%i]\n", fit_start, fit_stop);
	if(fit_stop - fit_start < num_entries)	{printf("Not enough data to fit.  Readjust fit_start, fit_stop.  Exiting\n"); return 1;}
	
	if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
    linenum++;

	ret = fscanf(fp,"k2=%f,%f,%f\n", &k2[0], &k2[1], &k2[2]);
	if (ret != 3 )  {printf("File format not valid, k2. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	ret = fscanf(fp,"k3=%f,%f,%f\n", &k3[0], &k3[1], &k3[2]);
	if (ret != 3 )  {printf("File format not valid, k3. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	ret = fscanf(fp,"k4=%f,%f,%f\n", &k4[0], &k4[1], &k4[2]);
	if (ret != 3 )  {printf("File format not valid, k4. Lines %i.  Exiting\n", linenum); fclose(fp); return 1;}
	linenum++;
	if(verbose) printf("APD0 Correction Factors:\t%e\t%e\t%e\n", k2[0], k3[0], k4[0]);
	if(verbose) printf("APD1 Correction Factors:\t%e\t%e\t%e\n", k2[1], k3[1], k4[1]);
	if(verbose) printf("APD2 Correction Factors:\t%e\t%e\t%e\n", k2[2], k3[2], k4[2]);

	if(fgets(line, 512, fp)==NULL)	{printf("Input file is prematurely short.  Exiting\n"); fclose(fp); return 1;};
    linenum++;

    int verifystring = 0;
    ret = fscanf(fp, "Random string to verify file completeness:  %i\n", &verifystring);
    if  ( verifystring != 7392342)  {
        printf("Verification string (%i) not found - check file integrity\n", verifystring);
        return 1;
    }
    fclose(fp);
        

    //  Read Data in From bin files
    char resultsfile01[80];	//Open file to write out data
    FILE * fpB;
    snprintf(resultsfile01, 50, "bin_clean%s.bin",&filetag[128*0]);
    fpB = fopen(resultsfile01, "rb");
    if (fpB == NULL)	{printf("***File %s Failed to Open, please verify this file exists***\nExiting program\n", resultsfile01); return 1;}
    int veclen_helper = 0;
    i=fread(&veclen_helper, sizeof(int), 1, fpB);
    int	FiIDAQuant = 0;
    i=fread(&FiIDAQuant, sizeof(int), 1, fpB);
    if(verbose) printf("Veclen = %i (%i), FiIDAQuant = %i\n", veclen_helper, i, FiIDAQuant);
    
    fseek (fpB , 0 , SEEK_END);
	long	file_size = ftell (fpB);
	if(verbose) printf("File Size is %li\n", file_size);   
	long curves_per_file = file_size -  3*sizeof(int) - veclen_helper*sizeof(double);
	curves_per_file = curves_per_file / ((7+12*veclen_helper)*sizeof(double));
	if(verbose) printf("Curves in file = %li (%li/%li)\n", curves_per_file, curves_per_file*((7+12*veclen_helper)*sizeof(double))+3*sizeof(int)+veclen_helper*sizeof(double),  file_size);	
	
    fclose(fpB);
    
    double  * tauarray;
    double	* time_array;
    double  * G_array;
    double  * G_stdev_array;
    double  * intensity0_array;
    double  * intensity0_stdev_array;
    double  * intensity1_array;
    double  * intensity1_stdev_array;
    double  * intensity2_array;
    double  * intensity2_stdev_array;
    tauarray = malloc(veclen_helper*sizeof(double));
    time_array = malloc(num_entries*sizeof(double));
    G_array = malloc(num_entries*veclen_helper*sizeof(double));
    G_stdev_array = malloc(num_entries*veclen_helper*sizeof(double));
    intensity0_array = malloc(num_entries*sizeof(double));
    intensity0_stdev_array = malloc(num_entries*sizeof(double));
    intensity1_array = malloc(num_entries*sizeof(double));
    intensity1_stdev_array = malloc(num_entries*sizeof(double));
    intensity2_array = malloc(num_entries*sizeof(double));
    intensity2_stdev_array = malloc(num_entries*sizeof(double));
    for(i=0; i<veclen_helper; i++) tauarray[i] = -6.6;
    for(i=0; i<num_entries*veclen_helper; i++) G_array[i] = -6.6;
    for(i=0; i<num_entries*veclen_helper; i++) G_stdev_array[i] = -6.6;
    for(i=0; i<num_entries; i++) intensity0_array[i] = -6.6;
    for(i=0; i<num_entries; i++) intensity0_stdev_array[i] = -6.6;
    for(i=0; i<num_entries; i++) intensity1_array[i] = -6.6;
    for(i=0; i<num_entries; i++) intensity1_stdev_array[i] = -6.6;
    for(i=0; i<num_entries; i++) intensity2_array[i] = -6.6;
    for(i=0; i<num_entries; i++) intensity2_stdev_array[i] = -6.6;
    
    
    int num_slices = 0;
    for (ff=0; ff<num_entries; ff++)   {
        snprintf(resultsfile01, 50, "bin_clean%s.bin",&filetag[128*ff]);
        fpB = fopen(resultsfile01, "rb");
        if (fpB == NULL)	printf("***File %s Failed to Open!***\n", resultsfile01);
        else if(verbose) printf("File %s Opened\n", resultsfile01);
        //  0   int = veclen
        fseek(fpB, 2*sizeof(int), SEEK_CUR);
        //  Vestigial, Number of slices (timepoints)
        j = fread(&num_slices, sizeof(int), 1, fpB);
        if (j != 1)	printf("***File %s Read Problem! (1)  j = %i***\n", resultsfile01, j);
        //  2-veclen+2  tau values
        if(ff==0)    {   //read in tau values
            j = fread(tauarray, sizeof(double), veclen_helper, fpB);
            if (j != veclen_helper)	printf("***File %s Read Problem! (2)  j = %i***\n", resultsfile01, j);
        }
        else   fseek(fpB, veclen_helper*sizeof(double), SEEK_CUR); 
        //  For Each timepoint, write Int0,Int1,Int2,SDInt0,SDInt1,SDInt2, 
        //      G0x0,G1x1,G2x2,G0x1,G1x2,G2x0,SDG0x0,SDG1x1,SDG2x2,SDG0x1,SDG1x2,SDG2x0
        
        //fseek for (6+12*veclen_helper)*sizeof(double)*timepoint
        fseek(fpB, (7+12*veclen_helper)*sizeof(double)*f1[ff], SEEK_CUR);
        //read int's
        if(fread(&time_array[ff], sizeof(double),1, fpB)!= 1) printf("Read Error: time\n");
        if(fread(&intensity0_array[ff], sizeof(double),1, fpB)!= 1) printf("Read Error: Int1\n");
        if(fread(&intensity1_array[ff], sizeof(double),1, fpB)!= 1) printf("Read Error: Int2\n");
        if(fread(&intensity2_array[ff], sizeof(double),1, fpB)!= 1) printf("Read Error: Int3\n");
        if(fread(&intensity0_stdev_array[ff], sizeof(double),1, fpB)!= 1) printf("Read Error: SDInt1\n");
        if(fread(&intensity1_stdev_array[ff], sizeof(double),1, fpB)!= 1) printf("Read Error: SDInt2\n");
        if(fread(&intensity2_stdev_array[ff], sizeof(double),1, fpB)!= 1) printf("Read Error: SDInt3\n");
        
		if(verbose) printf("Int original = %e\n", intensity1_array[ff]);
        
        //fseek for veclen*sizeof(double)*corr
        fseek(fpB, veclen_helper*sizeof(double)*corr_int[ff], SEEK_CUR);
        if(fread(&G_array[ff*veclen_helper], sizeof(double),veclen_helper, fpB)!= veclen_helper) printf("Read Error: Garray\n");
        //fseek for 6*veclen*sizeof(double)    
        fseek(fpB, veclen_helper*sizeof(double)*5, SEEK_CUR);
        if(fread(&G_stdev_array[ff*veclen_helper], sizeof(double),veclen_helper, fpB)!= veclen_helper) printf("Read Error: SDGarray\n");
            
        fclose(fpB);
        
    }
        
    //	Create an Array, populate w/ model params
    double	* x_array;
    double	* x_std_array;
    x_array = malloc(fcs_model_params*num_entries*sizeof(double));
    x_std_array = malloc(fcs_model_params*num_entries*sizeof(double));
    for(i=0; i<fcs_model_params*num_entries; i++)	x_array[i] = 0.0;
    for(i=0; i<fcs_model_params*num_entries; i++)	x_std_array[i] = 0.0;
	for(i=0; i<num_entries; i++)	{
		x_array[0+i*fcs_model_params] = N_t;
		x_array[1+i*fcs_model_params] = theta;
		x_array[2+i*fcs_model_params] = tauD_1;
		x_array[3+i*fcs_model_params] = tauD_2;
		x_array[4+i*fcs_model_params] = omega;
		x_array[5+i*fcs_model_params] = g_inf;
		x_array[6+i*fcs_model_params] = tauf_1;
		x_array[7+i*fcs_model_params] = tauf_2;
		x_array[8+i*fcs_model_params] = T0_1;
		x_array[9+i*fcs_model_params] = T0_2;
		x_array[10+i*fcs_model_params] = gamma_1;
		x_array[11+i*fcs_model_params] = gamma_2;
		x_array[12+i*fcs_model_params] = frac_gamma_1;
		x_array[13+i*fcs_model_params] = frac_gamma_2;
		x_array[14+i*fcs_model_params] = 0.0;
		x_array[15+i*fcs_model_params] = 0.0;
		x_array[16+i*fcs_model_params] = 0.0;
		x_array[17+i*fcs_model_params] = 0.0;
		x_array[18+i*fcs_model_params] = 0.0;
		if (corr_int[i] ==	0)	{
			x_array[14+i*fcs_model_params] = k2[0];
			x_array[15+i*fcs_model_params] = k3[0];
			x_array[16+i*fcs_model_params] = k4[0];
			x_array[17+i*fcs_model_params] = intensity0_array[i];
		}
		if (corr_int[i] ==	1)	{
			x_array[14+i*fcs_model_params] = k2[1];
			x_array[15+i*fcs_model_params] = k3[1];
			x_array[16+i*fcs_model_params] = k4[1];
			x_array[17+i*fcs_model_params] = intensity1_array[i];
		}
		if (corr_int[i] ==	2)	{
			x_array[14+i*fcs_model_params] = k2[2];
			x_array[15+i*fcs_model_params] = k3[2];
			x_array[16+i*fcs_model_params] = k4[2];
			x_array[17+i*fcs_model_params] = intensity2_array[i];
		}
		if (corr_int[i] ==	3)	{
			x_array[17+i*fcs_model_params] = intensity0_array[i];
			x_array[18+i*fcs_model_params] = intensity1_array[i];
		}
		if (corr_int[i] ==	4)	{
			x_array[17+i*fcs_model_params] = intensity1_array[i];
			x_array[18+i*fcs_model_params] = intensity2_array[i];
		}
		if (corr_int[i] ==	5)	{
			x_array[17+i*fcs_model_params] = intensity2_array[i];
			x_array[18+i*fcs_model_params] = intensity0_array[i];
		}
	}
		
	//	p_mask:	if (p_mask[i] > -1) x_array[i] = x[x_mask[i]]
	int	* p_mask;
	p_mask = malloc(fcs_model_params*num_entries*sizeof(int));
	double	* x_temp;
	x_temp = malloc(fcs_model_params*num_entries*sizeof(double));
	for(i=0; i<fcs_model_params*num_entries; i++)	p_mask[i] = -1;
	
	
	int	p_count = 0;
	int	p_placeholder = 0;
	if	(fix_N_t == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_N_t == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_N_t >4 || fix_N_t <0)	{printf("Input Error:  Invalid Fix Parameter for N_t (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_N_t); fix_N_t = 0;}
	
	p_placeholder++;
	
	if	(fix_theta == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_theta == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_theta >4 || fix_theta <0)	{printf("Input Error:  Invalid Fix Parameter for theta (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_theta); fix_theta = 0;}

	
	p_placeholder++;
	if	(fix_tauD_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauD_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_tauD_1 >4 || fix_tauD_1 <0)	{printf("Input Error:  Invalid Fix Parameter for tauD_1 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_tauD_1); fix_tauD_1 = 0;}

	
	p_placeholder++;
	if	(fix_tauD_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauD_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_tauD_2 >4 || fix_tauD_2 <0)	{printf("Input Error:  Invalid Fix Parameter for tauD_2 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_tauD_2); fix_tauD_2 = 0;}

	p_placeholder++;
	if	(fix_omega == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_omega == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_omega >4 || fix_omega <0)	{printf("Input Error:  Invalid Fix Parameter for omega (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_omega); fix_omega = 0;}

	p_placeholder++;
	if	(fix_g_inf == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_g_inf == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_g_inf >4 || fix_g_inf <0)	{printf("Input Error:  Invalid Fix Parameter for G(inf) (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_g_inf); fix_g_inf = 0;}

	p_placeholder++;
	if	(fix_tauf_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauf_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_tauf_1 >4 || fix_tauf_1 <0)	{printf("Input Error:  Invalid Fix Parameter for tauf_1 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_tauf_1); fix_tauf_1 = 0;}

	p_placeholder++;
	if	(fix_tauf_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauf_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_tauf_2 >4 || fix_tauf_2 <0)	{printf("Input Error:  Invalid Fix Parameter for tauf_2 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_tauf_2); fix_tauf_2 = 0;}

	p_placeholder++;
	if	(fix_T0_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_T0_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_T0_1 >4 || fix_T0_1 <0)	{printf("Input Error:  Invalid Fix Parameter for T0_1 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_T0_1); fix_T0_1 = 0;}

	p_placeholder++;
	if	(fix_T0_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_T0_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_T0_2 >4 || fix_T0_2 <0)	{printf("Input Error:  Invalid Fix Parameter for T0_2 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_T0_2); fix_T0_2 = 0;}

	p_placeholder++;
	if	(fix_gamma_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_gamma_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_gamma_1 >4 || fix_gamma_1 <0)	{printf("Input Error:  Invalid Fix Parameter for gamma_1 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_theta); fix_theta = 0;}

	p_placeholder++;
	if	(fix_gamma_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_gamma_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_gamma_2 >4 || fix_gamma_2 <0)	{printf("Input Error:  Invalid Fix Parameter for gamma_2 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_gamma_2); fix_gamma_2 = 0;}

	p_placeholder++;
	if	(fix_frac_gamma_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_frac_gamma_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_frac_gamma_1 >4 || fix_frac_gamma_1 <0)	{printf("Input Error:  Invalid Fix Parameter for frac_gamma_1 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_frac_gamma_1); fix_frac_gamma_1 = 0;}

	p_placeholder++;
	if	(fix_frac_gamma_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_frac_gamma_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	else	if	(fix_frac_gamma_2 >4 || fix_frac_gamma_2 <0)	{printf("Input Error:  Invalid Fix Parameter for frac_gamma_2 (user input = %i).  \nAppropriate values: \n\t0 (fixed), \n\t1 (fit with one value for all species on the second pass only), \n\t2 (fit with one value for per species on the second pass only),\n\t3 (fit with one value for all species on the first and second passes),\n\t4 (fit with one value for per species on the first and second passes).\n---> Value changed to 0.\n\n",fix_frac_gamma_2); fix_frac_gamma_2 = 0;}

	p_placeholder++;

    printf("Fitting round 1:\tNumber of Parameters:\t%i\n", p_count);

	
/******************** Initiate GSL Levenberg-Marquardt Solver ********************/

    const gsl_multifit_fdfsolver_type *T;
	T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s = 0;
    gsl_matrix *covar = 0;
    int status;
    unsigned int iter = 0;
	double *	y; 
    double *	sigma;
	size_t nFIT;
    size_t nfull;
	nFIT = (fit_stop-fit_start)*num_entries;
    nfull = veclen*num_entries;
	y = malloc(nfull*sizeof(double));
	sigma = malloc(nfull*sizeof(double));
    size_t p;
    struct data d;
	d.nFIT = nFIT;
	d.nfull = nfull;
	d.y = y;
	d.sigma = sigma;
	d.num_entries = num_entries;
	d.fit_start = fit_start;
	d.fit_stop = fit_stop;
	d.x_array = x_array;
	gsl_multifit_function_fdf f;    
    f.f = &gtau_f;
    f.df = &gtau_df;
    f.fdf = &gtau_fdf;
    f.n = nFIT;
	f.params = &d; 
	for (i=0; i<nfull; i++) {
        y[i] = G_array[i];
        sigma[i] = G_stdev_array[i];
    }
	gsl_vector * x = 0;
	double dof = 0.0;
	double chi = 0.0;
    double c = 0.0; 

	//	***	First Fitting w/ N, G(inf)	***	//
	if(p_count > 0)	{
		p = p_count;
		d.p = p;
		d.p_mask = p_mask;
		covar = gsl_matrix_alloc (p, p);
		x = gsl_vector_alloc (p);
		for (i=0; i<p; i++)		gsl_vector_set (x, i, x_temp[i]);
		
		f.p = p;
	
		s = gsl_multifit_fdfsolver_alloc (T, nFIT, p);
		gsl_multifit_fdfsolver_set (s, &f, x);
	  
		dof = nFIT - p;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);
			printf ("\tStatus = %s, Red.Chi^2 = %.8f\n", gsl_strerror (status), pow(gsl_blas_dnrm2(s->f), 2.0)/dof);
			if (status) break;
			status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
		} while (status == GSL_CONTINUE && iter < 500);
		
		gsl_multifit_covar (s->J, 0.0, covar);
		chi = gsl_blas_dnrm2(s->f);	// Euclidean norm of f
		c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
		printf("Fitting Round 1 complete:  chisq/dof = %g\nFit parameters:\n",  pow(chi, 2.0) / dof);
		
		for(i=0; i<fcs_model_params*num_entries; i++)	
			if (p_mask[i] != -1) x_array[i] = gsl_vector_get(s->x, p_mask[i]);
		for(i=0; i<fcs_model_params*num_entries; i++)	
			if (p_mask[i] != -1) x_std_array[i] = sqrt(gsl_matrix_get(covar,p_mask[i],p_mask[i]));
					
	}	//End clause (if p_count > 0)
		
	p_count = 0;
	p_placeholder = 0;
	if	(fix_N_t == 1 || fix_N_t == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_N_t == 2 || fix_N_t == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_theta == 1 || fix_theta == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_theta == 2 || fix_theta == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_tauD_1 == 1 || fix_tauD_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauD_1 == 2 || fix_tauD_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_tauD_2 == 1 || fix_tauD_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauD_2 == 2 || fix_tauD_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_omega == 1 || fix_omega == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_omega == 2 || fix_omega == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_g_inf == 1 || fix_g_inf == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_g_inf == 2 || fix_g_inf == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_tauf_1 == 1 || fix_tauf_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauf_1 == 2 || fix_tauf_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_tauf_2 == 1 || fix_tauf_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_tauf_2 == 2 || fix_tauf_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_T0_1 == 1 || fix_T0_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_T0_1 == 2 || fix_T0_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_T0_2 == 1 || fix_T0_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_T0_2 == 2 || fix_T0_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_gamma_1 == 1 || fix_gamma_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_gamma_1 == 2 || fix_gamma_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_gamma_2 == 1 || fix_gamma_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_gamma_2 == 2 || fix_gamma_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_frac_gamma_1 == 1 || fix_frac_gamma_1 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_frac_gamma_1 == 2 || fix_frac_gamma_1 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;
	if	(fix_frac_gamma_2 == 1 || fix_frac_gamma_2 == 3)	{
		for(j=0; j< num_entries; j++)	p_mask[p_placeholder+j*fcs_model_params] = p_count;
		x_temp[p_count] = x_array[p_placeholder];
		p_count++;
	}
	else	if	(fix_frac_gamma_2 == 2 || fix_frac_gamma_2 == 4)	{
		for(j=0; j< num_entries; j++)	{
			p_mask[p_placeholder+j*fcs_model_params] = p_count;
			x_temp[p_count] = x_array[p_placeholder+j*fcs_model_params];
			p_count++;
		}
	}
	p_placeholder++;

		
	for(j=0; j< num_entries; j++)	{
		printf("N: %.4f/%.4f Th: %.4f/%.4f TD1: %.4f/%.4f TD2 :%.4f/%.4f \n\tw: %.4f/%.4f G(inf):%.4f/%.4f \n\tTf1: %.4f/%.4f Tf2: %.4f/%.4f T0,1: %.4f/%.4f T0,2: %.4f/%.4f \n\tGm1: %.4f/%.4f Gm2: %.4f/%.4f fGm1: %.4f/%.4f fGm2: %.4f/%.4f\n", x_array[0+j*fcs_model_params], x_std_array[0+j*fcs_model_params], x_array[1+j*fcs_model_params], x_std_array[1+j*fcs_model_params], x_array[2+j*fcs_model_params], x_std_array[2+j*fcs_model_params], x_array[3+j*fcs_model_params], x_std_array[3+j*fcs_model_params], x_array[4+j*fcs_model_params], x_std_array[4+j*fcs_model_params], x_array[5+j*fcs_model_params], x_std_array[5+j*fcs_model_params], x_array[6+j*fcs_model_params], x_std_array[6+j*fcs_model_params], x_array[7+j*fcs_model_params], x_std_array[7+j*fcs_model_params], x_array[8+j*fcs_model_params], x_std_array[8+j*fcs_model_params], x_array[9+j*fcs_model_params], x_std_array[9+j*fcs_model_params], x_array[10+j*fcs_model_params], x_std_array[10+j*fcs_model_params], x_array[11+j*fcs_model_params], x_std_array[11+j*fcs_model_params], x_array[12+j*fcs_model_params], x_std_array[12+j*fcs_model_params], x_array[13+j*fcs_model_params], x_std_array[13+j*fcs_model_params]);
		printf("\n");
	}


	if(covar != 0)	{
		gsl_vector_free (x);
		gsl_multifit_fdfsolver_free (s);
		gsl_matrix_free (covar);  

	}

	/***	Second Fitting w/ All Parameters	***/
	printf("\nFitting round 2:\tNumber of Parameters:\t%i\n", p_count);
	double	overallredchisq=0.0;
	p=0;
	if(p_count > 0)	{
		p = p_count;
		d.p = p;
		d.p_mask = p_mask;
		covar = gsl_matrix_alloc (p, p);
		x = gsl_vector_alloc (p);
		for (i=0; i<p; i++)		gsl_vector_set (x, i, x_temp[i]);
		
		f.p = p;

		s = gsl_multifit_fdfsolver_alloc (T, nFIT, p);
		gsl_multifit_fdfsolver_set (s, &f, x);
	  
		dof = nFIT - p;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);
			printf ("\tStatus = %s, Red.Chi^2 = %.8f\n", gsl_strerror (status), pow(gsl_blas_dnrm2(s->f), 2.0)/dof);
			if (status) break;
			status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
		} while (status == GSL_CONTINUE && iter < 500);
		
		gsl_multifit_covar (s->J, 0.0, covar);
		chi = gsl_blas_dnrm2(s->f);	// Euclidean norm of f
		c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
		printf("Fitting Round 2 complete:  chisq/dof = %g\nFit parameters:\n",  pow(chi, 2.0) / dof);    
		overallredchisq = pow(chi, 2.0) / dof;
	}
	else	printf("*********************** Input Error *****************\nNo parameters were set to fit in the input file.  \nSet at least one, e.g.  N_t=0.556,2\n*****************************************************\n");

/******************** Print Curve-Specific Chi Sq's to Screen ********************/    
	//	Generate Fit Lines    
	for(i=0; i<fcs_model_params*num_entries; i++)	
		if (p_mask[i] != -1) x_array[i] = gsl_vector_get(s->x, p_mask[i]);
	for(i=0; i<fcs_model_params*num_entries; i++)	
		if (p_mask[i] != -1) x_std_array[i] = sqrt(gsl_matrix_get(covar,p_mask[i],p_mask[i]));
    
    double	* fits;
    double	* w_resids;
    fits = malloc(num_entries*veclen*sizeof(double));
    w_resids = malloc(num_entries*veclen*sizeof(double));
    for(i=0; i<num_entries*veclen; i++)	fits[i] = 0.0;
	for(i=0; i<num_entries*veclen; i++)	w_resids[i] = 0.0;
    
	struct m_pass_struct o;
	o.m = malloc(veclen*sizeof(double));

    for(j=0; j<num_entries; j++)	{
		i=0;
		o.N_t = x_array[(i++)+j*fcs_model_params];
		o.theta = x_array[(i++)+j*fcs_model_params];
		o.tauD_1 = x_array[(i++)+j*fcs_model_params];
		o.tauD_2 = x_array[(i++)+j*fcs_model_params];
		o.omega = x_array[(i++)+j*fcs_model_params];
		o.g_inf = x_array[(i++)+j*fcs_model_params];
		o.tauf_1 = x_array[(i++)+j*fcs_model_params];
		o.tauf_2 = x_array[(i++)+j*fcs_model_params];
		o.T0_1 = x_array[(i++)+j*fcs_model_params];
		o.T0_2 = x_array[(i++)+j*fcs_model_params];
		o.gamma_1 = x_array[(i++)+j*fcs_model_params];
		o.gamma_2 = x_array[(i++)+j*fcs_model_params];
		o.frac_gamma_1 = x_array[(i++)+j*fcs_model_params];
		o.frac_gamma_2 = x_array[(i++)+j*fcs_model_params];
		o.k2 = x_array[(i++)+j*fcs_model_params];
		o.k3 = x_array[(i++)+j*fcs_model_params];
		o.k4 = x_array[(i++)+j*fcs_model_params];
		o.int0 = x_array[(i++)+j*fcs_model_params];

		m_tau(&o);
		
		for(i=0; i<veclen; i++)	fits[i+j*veclen] = o.m[i];
		for(i=0; i<veclen; i++)	w_resids[i+j*veclen] = (o.m[i] - G_array[i+j*veclen])/(G_stdev_array[i+j*veclen]);
    } 
    free(o.m);

	double * chi_vec;
	chi_vec = malloc(num_entries*sizeof(double));
	for(j=0; i< num_entries; i++)	chi_vec[j] = 0.0;
	gsl_vector	* ch = gsl_vector_alloc (fit_stop-fit_start);
	dof = dof / (float) num_entries;	//Not Strictly true...
    for(j=0; j< num_entries; j++)	{
    	for(i=0; i<fit_stop-fit_start; i++)	gsl_vector_set(ch, i, w_resids[i+j*veclen+fit_start]);
    	chi_vec[j] = gsl_blas_dnrm2(ch);
	}
    gsl_vector_free(ch);
    
/***********	Print Final Parameters To Screen	***********/		
	char datatag[80] = "";
	snprintf(resultsfile01, 50, "F3CS_LocalFit%s.fit_params.txt",&filetag[128*0]);
	fpB = fopen(resultsfile01, "w");
	fprintf(fpB, "Overall Red. Chi Sq: %.4f\n", overallredchisq);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		printf("%s\tRed.Chi^2: %.4f N: %.4f/%.4f Th: %.4f/%.4f\n\tTD1: %.4f/%.4f TD2 :%.4f/%.4f w: %.4f/%.4f G(inf):%.4f/%.4f \n\tTf1: %.4f/%.4f Tf2: %.4f/%.4f T0,1: %.4f/%.4f T0,2: %.4f/%.4f \n\tGm1: %.4f/%.4f Gm2: %.4f/%.4f fGm1: %.4f/%.4f fGm2: %.4f/%.4f\n", datatag, pow(chi_vec[j], 2.0) / dof, x_array[0+j*fcs_model_params], x_std_array[0+j*fcs_model_params], x_array[1+j*fcs_model_params], x_std_array[1+j*fcs_model_params], x_array[2+j*fcs_model_params], x_std_array[2+j*fcs_model_params], x_array[3+j*fcs_model_params], x_std_array[3+j*fcs_model_params], x_array[4+j*fcs_model_params], x_std_array[4+j*fcs_model_params], x_array[5+j*fcs_model_params], x_std_array[5+j*fcs_model_params], x_array[6+j*fcs_model_params], x_std_array[6+j*fcs_model_params], x_array[7+j*fcs_model_params], x_std_array[7+j*fcs_model_params], x_array[8+j*fcs_model_params], x_std_array[8+j*fcs_model_params], x_array[9+j*fcs_model_params], x_std_array[9+j*fcs_model_params], x_array[10+j*fcs_model_params], x_std_array[10+j*fcs_model_params], x_array[11+j*fcs_model_params], x_std_array[11+j*fcs_model_params], x_array[12+j*fcs_model_params], x_std_array[12+j*fcs_model_params], x_array[13+j*fcs_model_params], x_std_array[13+j*fcs_model_params]);
		printf("\n");
		fprintf(fpB, "%s\tRed.Chi^2: %.4f N: %.4f/%.4f Th: %.4f/%.4f\n\tTD1: %.4f/%.4f TD2 :%.4f/%.4f w: %.4f/%.4f G(inf):%.4f/%.4f \n\tTf1: %.4f/%.4f Tf2: %.4f/%.4f T0,1: %.4f/%.4f T0,2: %.4f/%.4f \n\tGm1: %.4f/%.4f Gm2: %.4f/%.4f fGm1: %.4f/%.4f fGm2: %.4f/%.4f\n", datatag, pow(chi_vec[j], 2.0) / dof, x_array[0+j*fcs_model_params], x_std_array[0+j*fcs_model_params], x_array[1+j*fcs_model_params], x_std_array[1+j*fcs_model_params], x_array[2+j*fcs_model_params], x_std_array[2+j*fcs_model_params], x_array[3+j*fcs_model_params], x_std_array[3+j*fcs_model_params], x_array[4+j*fcs_model_params], x_std_array[4+j*fcs_model_params], x_array[5+j*fcs_model_params], x_std_array[5+j*fcs_model_params], x_array[6+j*fcs_model_params], x_std_array[6+j*fcs_model_params], x_array[7+j*fcs_model_params], x_std_array[7+j*fcs_model_params], x_array[8+j*fcs_model_params], x_std_array[8+j*fcs_model_params], x_array[9+j*fcs_model_params], x_std_array[9+j*fcs_model_params], x_array[10+j*fcs_model_params], x_std_array[10+j*fcs_model_params], x_array[11+j*fcs_model_params], x_std_array[11+j*fcs_model_params], x_array[12+j*fcs_model_params], x_std_array[12+j*fcs_model_params], x_array[13+j*fcs_model_params], x_std_array[13+j*fcs_model_params]);
		fprintf(fpB, "\n");
	}
	fclose(fpB);
	
	printf("Errors:");  
	for(i=0; i<p; i++)	printf("  %3.2e", c*sqrt(gsl_matrix_get(covar,i,i)));
	printf("\n");  

	printf("\nOverall Red. Chi Sq: %.4f\n", overallredchisq);
    
/******************** Write Data to Disk ********************/
    
    double	G_max = 0;
    double	G_min = 10;
    for(j=0; j< num_entries; j++)	{
    	for(i=fit_start; i<fit_stop; i++)	{
    		if(G_array[j*veclen_helper+i] > G_max)	G_max = G_array[j*veclen_helper+i];
    		if(G_array[j*veclen_helper+i] < G_min)	G_min = G_array[j*veclen_helper+i];
		}
    }
    G_min = G_min - (G_max-G_min) * 0.025;
    G_max = G_max + (G_max-G_min) * 0.025;
    
    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.dat",&filetag[128*0]);
	fpB = fopen(resultsfile01, "w");
	for(i=fit_start; i< fit_stop; i++){
        fprintf(fpB, "%e", tau_G(i));
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", G_array[j*veclen_helper+i]);
        }
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", G_stdev_array[j*veclen_helper+i]);
        }
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", fits[j*veclen_helper+i]);
        }
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", w_resids[j*veclen_helper+i]);
        }
    	fprintf(fpB, "\n");
	}
	fclose(fpB);
	
	char resultsfile02[80];
	snprintf(resultsfile02, 50, "gnuplot%s.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02);
	fpB = fopen(resultsfile02, "w");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s-%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
    	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");
	fclose(fpB);
	
	
	snprintf(resultsfile02, 50, "gnuplot%s.fits.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
	fpB = fopen(resultsfile02, "w");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.fits.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		fprintf(fpB, ", ");
	}
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "fit_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");

	fclose(fpB);

	snprintf(resultsfile02, 50, "gnuplot%s.errors.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02); 
	fpB = fopen(resultsfile02, "w");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.errors.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i:%i with yerrorbars title \"%s\" lt palette frac %.2f", resultsfile01, j+2, j+2+num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");
	fclose(fpB);

    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.resids.dat",&filetag[128*0]);
	fpB = fopen(resultsfile01, "w");
	for(i=fit_start; i< fit_stop; i++){
        fprintf(fpB, "%e", tau_G(i));
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", w_resids[j*veclen+i]);
        }
    	fprintf(fpB, "\n");
	}
	fclose(fpB);

	snprintf(resultsfile02, 50, "gnuplot%s.resids.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02);
	fpB = fopen(resultsfile02, "w");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.resids.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "w_resid_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");
	fclose(fpB);

    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.Nscaled.dat",&filetag[128*0]);
	fpB = fopen(resultsfile01, "w");
	for(i=fit_start; i< fit_stop; i++){
        fprintf(fpB, "%e", tau_G(i));
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", (G_array[j*veclen_helper+i] - x_array[5+j*fcs_model_params]) * x_array[0+j*fcs_model_params]);
        }
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", G_stdev_array[j*veclen_helper+i] * x_array[0+j*fcs_model_params]);
        }
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", (fits[j*veclen_helper+i] - x_array[5+j*fcs_model_params]) * x_array[0+j*fcs_model_params]);
        }
        for(j=0; j< num_entries; j++)	{
        	fprintf(fpB, "\t%f", w_resids[j*veclen_helper+i]);
        }
    	fprintf(fpB, "\n");
	}
	fclose(fpB);

	snprintf(resultsfile02, 50, "gnuplot%s.Nscaled.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02);
	fpB = fopen(resultsfile02, "w");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.Nscaled.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][-0.05:1.2] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaled_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");
	fclose(fpB);

	snprintf(resultsfile02, 50, "gnuplot%s.Nscaledfits.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02);
	fpB = fopen(resultsfile02, "w");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.Nscaledfits.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][-0.05:1.2] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaled_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		fprintf(fpB, ", ");
	}

	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaledfit-%3s%s-%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	

	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");
	fclose(fpB);
	
	snprintf(resultsfile01, 50, "F3CS_LocalFit%s.fit_params.%s.dat",&filetag[128*0], argv[2]);
	fpB = fopen(resultsfile01, "w");
	fprintf(fpB, "RxnTag_%s\tN_%s\tTheta_%s\tTauD_1_%s\tTauD_2_%s\tOmega_%s\tG_inf_%s\tTauf_1_%s\tTauf_2_%s\tT0_1_%s\tT0_2_%s\tGamma_1_%s\tGamma_2_%s\tFrac_gamma_1_%s\tFrac_gamma_2_%s\tStD_N_%s\tStD_Theta_%s\tStD_TauD_1_%s\tStD_TauD_2_%s\tStD_Omega_%s\tStD_G_inf_%s\tStD_Tauf_1_%s\tStD_Tauf_2_%s\tStD_T0_1_%s\tStD_T0_2_%s\tStD_Gamma_1_%s\tStD_Gamma_2_%s\tStD_Frac_gamma_1_%s\tStD_Frac_gamma_2_%s\tChiSq_%s\tInt_0_%s\tInt_1_%s\tInt_2_%s\n", argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2],argv[2]);
	for(j=0; j< num_entries; j++)	{
        fprintf(fpB, "%3s%s", &corr_mode[128*j], &filetag[128*j]);
		for(i=0; i< 14; i++)	fprintf(fpB, "\t%f", x_array[i+j*fcs_model_params]);
        for(i=0; i< 14; i++)	fprintf(fpB, "\t%f", x_std_array[i+j*fcs_model_params]);
    	fprintf(fpB, "\t%.9f", pow(chi_vec[j], 2.0) / dof);
        fprintf(fpB, "\t%f", intensity0_array[j]);
        fprintf(fpB, "\t%f", intensity1_array[j]);
        fprintf(fpB, "\t%f\n", intensity2_array[j]);
	}
	fclose(fpB);
	printf(".dat file: %s\n",resultsfile01);


/*******************  Prepare gnuplot scripts   ***********************/
	snprintf(resultsfile01, 50, "F3CS_LocalFit%s.dat",&filetag[128*0]);

	snprintf(resultsfile02, 50, "gnuplot%s.all.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02);
	fpB = fopen(resultsfile02, "w");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s-%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");

	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.fits.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		fprintf(fpB, ", ");
	}
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "fit_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");

	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set  ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.errors.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i:%i with yerrorbars title \"%s\" lt palette frac %.2f", resultsfile01, j+2, j+2+num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");

    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.dat", &filetag[128*0]);

	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set logscale y\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot ");
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != 8 -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");
	fprintf(fpB, "unset logscale\n");

    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.resids.dat",&filetag[128*0]);

	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.resids.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][] ",tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "w_resid_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");

    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.Nscaled.dat",&filetag[128*0]);

	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.Nscaled.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][-0.05:1.2] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaled_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");

	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"tau (s)\"\n");
	fprintf(fpB, "set ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");
	fprintf(fpB, "set term postscript eps enhanced color\n");
	snprintf(resultsfile02, 50, "gnuplot%s.Nscaledfits.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "plot [%2.2e:%2.2e][-0.05:1.2] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaled_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		fprintf(fpB, ", ");
	}
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaledfit-%3s%s-%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}

	fprintf(fpB, "\n");
	fprintf(fpB, "set term x11\n");
	fprintf(fpB, "replot\n");
	fclose(fpB);

/*********** Plot all plots together on one page, good for notebooks  ***********/
	snprintf(resultsfile01, 50, "F3CS_LocalFit%s.dat",&filetag[128*0]);
	snprintf(resultsfile02, 50, "gnuplot%s.multiploteps.txt",&filetag[128*0]);
	printf("For Plotting: gnuplot> load \'%s\'\n",resultsfile02);
	fpB = fopen(resultsfile02, "w");
	
	fprintf(fpB, "set term postscript eps enhanced color size 9, 6.3 solid\n");
	snprintf(resultsfile02, 50, "multieps.gnuplot%s.eps",&filetag[128*0]);
	fprintf(fpB, "set output '%s'\n", resultsfile02);
	fprintf(fpB, "set multiplot\n");
	fprintf(fpB, "set size 0.5,0.5\n");
	fprintf(fpB, "set bmargin 0\n");
	fprintf(fpB, "set tmargin 2\n");
	fprintf(fpB, "set lmargin 10\n");
	fprintf(fpB, "set rmargin 2\n");
	fprintf(fpB, "set origin 0.0,0.5\n");
	fprintf(fpB, "set nokey\n");
	fprintf(fpB, "set logscale x\n");
	fprintf(fpB, "set xlabel \"\"\n");
	fprintf(fpB, "set format x \"\"\n");
	fprintf(fpB, "set ylabel \"G(tau)\"\n");
	fprintf(fpB, "set palette defined ( 0 \"green\", 1 \"blue\", 2 \"red\", 3 \"orange\" )\n");
	fprintf(fpB, "unset colorbox\n");

	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		fprintf(fpB, ", ");
	}
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "fit_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
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
	fprintf(fpB, "set ylabel \"G(tau)\"\n");
	fprintf(fpB, "unset colorbox\n");

	fprintf(fpB, "plot [%2.2e:%2.2e][%f:%f] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1, G_min, G_max);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i:%i with yerrorbars title \"%s\" lt palette frac %.2f", resultsfile01, j+2, j+2+num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
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
	fprintf(fpB, "set ylabel \"G(tau) (Norm.)\"\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");

    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.Nscaled.dat",&filetag[128*0]);
	fprintf(fpB, "plot [%2.2e:%2.2e][-0.05:1.2] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaled_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		fprintf(fpB, ", ");
	}
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "scaledfit-%3s%s-%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2+2*num_entries, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
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
	fprintf(fpB, "set ylabel \"Weighted Residuals\"\n");
	fprintf(fpB, "unset colorbox\n");
	fprintf(fpB, "set ticslevel 0.0\n");

    snprintf(resultsfile01, 50, "F3CS_LocalFit%s.resids.dat",&filetag[128*0]);

	fprintf(fpB, "plot [%2.2e:%2.2e][] ", tau_G(fit_start)*0.9, tau_G(fit_stop-1)*1.1);
	for(j=0; j< num_entries; j++)	{
		snprintf(datatag, 80, "w_resid_%3s%s_%.2f", &corr_mode[128*j], &filetag[128*j], f1[j]);
		i=0;
		do {
			if(datatag[i] == '_') datatag[i] = '-';
			i++;
		}	while (datatag[i] != '\0');
		fprintf(fpB, "\"%s\" using 1:%i with lines title \"%s\" lt palette frac %.2f", resultsfile01, j+2, datatag, (float)j / (float) num_entries);
		if (j != num_entries -1) fprintf(fpB, ", ");
	}
	fprintf(fpB, "\n");
	fprintf(fpB, "unset multiplot\n");
	fprintf(fpB, "set size 1, 1\n");
	fprintf(fpB, "set origin 0, 0\n");
	fclose(fpB);

/********* Print Covariance Matrix To Screen **********/
    snprintf(resultsfile01, 50, "F3CS_LocalFit_covar%s.dat",&filetag[128*0]);
	printf("Covariance Matrix Written to: %s\n",resultsfile01);
	fpB = fopen(resultsfile01, "w");
    for(i=0; i<p; i++)
        for(j=0; j<p; j++)
            fprintf(fpB, "%i\t%i\t%3.8e\n", i, j, gsl_matrix_get(covar,i,j));
    fclose(fpB);
    
    int   gridsize = 6;
    snprintf(resultsfile01, 50, "F3CS_LocalFit_normcovar%s.dat",&filetag[128*0]);
	printf("Norm. Covariance Matrix Written to: %s\n",resultsfile01);
	fpB = fopen(resultsfile01, "w");
    for(i=0; i<p*gridsize; i++)  {
        for(j=0; j<p*gridsize; j++)  {
            fprintf(fpB, "%i\t%i\t%3.8e\n", i, j, gsl_matrix_get(covar,i/gridsize,j/gridsize)/(sqrt(gsl_matrix_get(covar,i/gridsize,i/gridsize)*gsl_matrix_get(covar,j/gridsize,j/gridsize))));
        }
        fprintf(fpB, "\n");
    }
    fclose(fpB);
    
/******************** Free memory ********************/
    free(tauarray);
    free(time_array);
    gsl_vector_free (x);
	gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);  
    free(y);
    free(sigma);
    free(G_array);
    free(G_stdev_array);
    free(intensity0_array);
    free(intensity0_stdev_array);
    free(intensity1_array);
    free(intensity1_stdev_array);
    free(intensity2_array);
    free(intensity2_stdev_array);
    free(x_array);
    free(x_std_array);
    free(p_mask);
    free(x_temp);
    free(fits);
    free(w_resids);
    free(chi_vec);

    printf("Program completed.\n");
    return 0;
}
