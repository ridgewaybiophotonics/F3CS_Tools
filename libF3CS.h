/*
 *  libF3CS.h
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
 *  libF3CS.h:
 *  Sets the prototypes, memory structures and constants used by the entire 
 *  suite.  Of particular note, the n and pmax values for FCS integrals and 
 *  n_GGG and pmax_GGG values for F3CS integrals are set here on lines 30-34.  
 *  (The programs must be recompiled (make clean; make) for these changes to 
 *  have effect.)
 */

#include <stdlib.h>
#include <math.h>

// Users can change the n and pmax values below in order to change the time-range and spacing of time-points in the correlation integrals calculated by F3CS_2FCS and F3CS_AxAxA, F3CS_AxAxB and F3CS_AxBxG.  Note discussion in Fig. 7 of the main manuscript.  (The programs must be recompilied (make clean; make) for changes to these values to have effect.)
#define pmax       12
#define n          32
//  Same as above, but for computing triple correlations:
#define pmax_GGG   15
#define n_GGG      4


int bits_to_words(unsigned short *one, unsigned short *two, unsigned short *thr, int len, unsigned int *remainder, unsigned char *data);
void    *stream_sipper_thread(void *threadarg);
void    *stream_sipper_thread_test(void *threadarg);

void    *TripleCorr(void *threadarg);
void    *zbz_zboCorr(void *threadarg);
void    *obo_Corr(void *threadarg);

#ifndef _NI_uInt8_DEFINED_
#define _NI_uInt8_DEFINED_ 
typedef unsigned char uInt8;
#endif

#define FCS_BUFFER  4096
#define BURST_BUFFER 1024
#define rambuffer  100
#define samplesize 131073
#define curvestobin 16
#define anticipatedtime 12
#define	veclen n*(1+pmax)
#define	veclen_GGG n_GGG*(1+pmax_GGG)
//#define num_colours 2

#define NUM_THREADS     32
#define FILE_STRING_LENGTH 256

#define	num_curves_per_thread	2

unsigned short  *array0_0;
unsigned short  *array0_1;
unsigned short	*array0_2;
unsigned short  *array1_0;
unsigned short  *array1_1;
unsigned short	*array1_2;
unsigned short  *interleaved;

float           g0b0[veclen];
float           g0b1[(veclen)];
float           g0b2[(veclen)];
float           g1b0[(veclen)];
float           g1b1[(veclen)];
float           g1b2[(veclen)];
float           g2b0[(veclen)];
float           g2b1[(veclen)];
float           g2b2[(veclen)];
int             klaus0[pmax];
int             klaus1[pmax];
int				klaus2[pmax];
long long int   shatzel0[(veclen)];
long long int   shatzel1[(veclen)];
long long int   shatzel2[(veclen)];
unsigned short int  h_shatzel0[2*n];
unsigned short int  h_shatzel1[2*n];
unsigned short int  h_shatzel2[2*n];
int             stemp0[pmax];
int             stemp1[pmax];
int             stemp2[pmax];

long long int   total_counter;
double   tc0;
double   tc1;
double	 tc2;
double   tc0sq;
double   tc1sq;
double	 tc2sq;        
int             algorquantum;
int             algorquantum0b0;
int             algorquantum1b1;

struct fcs_variables {
	float *		g0b0;
	float *		g0b1;
	float *		g1b1;
	float *		g1b2;
	float *		g2b0;
	float *		g2b2;
	float *		g0b1b2;
	float *		g1b2b0;
	float *		g2b0b1;
	float *		klaus0;
	float *		klaus1;
	float *		klaus2;
	float *		shatzel0;
	float *		shatzel1;
	float *		shatzel2;
	float *		stemp0;
	float *		stemp1;
	float *		stemp2;
	double		tc0;
	double		tc1;
	double		tc2;
	float *		sc0;
	float *		sc1;
	float *		sc2;
	double		tc0sq;
	double		tc1sq;
	double		tc2sq;
};

//struct triplecorrstruct {
//    unsigned short  *array0;
//    unsigned short  *array1;
//    unsigned short	*array2;
//};

struct triplecorr_struct {
    unsigned short  *array0;
    unsigned short  *array1;
    unsigned short	*array2;
	struct fcs_variables * fcsv;
};

struct autocorrdata {
    float *data;
    int    len;
    float *G;
};

struct readrubbish {
    float *array;
    int    length;
    char  *filename;
};

struct	stream_struct	{
	int	is_first_file;
	char *file_nametag;
   	int	first_file;
   	int last_file;    
    int curvestobin_rxn;
    double	time_per_file;
    int	core_id;
    int	u_shorts_per_file;
	int	debug_flag;
};

