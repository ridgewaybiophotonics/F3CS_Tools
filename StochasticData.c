/*
 *  StochasticData.c
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
 *  StochasticData.c
 *  Self-contained code to generate F3CS_StochasticData.  The rates 
 *  displayed in the tutorial figure 1 can be set by altering lines 
 *  120-125.  The fake optical artifact can be eliminated on lines 
 *  226-229.
 *
 */

#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define	algorquant		131072
#define TUNE_BUFFER		2048

#define NN 10

struct  paramAcq    {
    int     run;
    int     contReadFreq[4];
    int     activeStates;
    int     clock_frequency_array_value;
    int     correlation_mode;
    int     clock_freq_scalar_value;
};
struct  paramAcq NICard;

struct write_stuff {
	int	key;
	int	filenum;
	int	len;
	char *tag;
	unsigned short * array;
};


void    *getchar_fxn  (void * val)	{
	do	{
		printf("Press 'q' to halt data generation\n");
	}	while ((char) getchar() != 'q');
	NICard.run = 0;
	printf("* * * * * * * * * * * *\nData generation stopping...\n* * * * * * * * * * * *\n");
	pthread_exit(NULL);
}

void    *write_fxn  (void * pass)	{
	struct write_stuff * val;
	val = pass;
	char resultsfile1[80];
	char numb[5];
	snprintf(numb, 5, "%i",val->filenum);
	strcpy(resultsfile1, "rxn_");
	strcat(resultsfile1, numb);
	strcat(resultsfile1, val->tag);
	strcat(resultsfile1, ".w3c");
	
	FILE * fp;
    fp = fopen(resultsfile1, "wb");
	if (fp == NULL) {
		printf ("Error:  File %s Failed To Open!\nKilling program...\n", resultsfile1);
		NICard.run = 0;
		}
	int	data_written = 0;
	data_written = fwrite(val->array, sizeof(unsigned short), val->len, fp);
	if (data_written != val->len)	{
		printf ("Error:  Could not write all data to file %s!\nKilling program...\n", resultsfile1);
		NICard.run = 0;
	}
	else	printf("\t\t\tFile write complete: %s \n", resultsfile1);
	fclose(fp);
	val->key = 0;
	pthread_exit(NULL);
}

int main(int argc, char *argv[])    {

	if(argc != 3)  {
        printf("Usage: F3CS_StochasticData tag_string time(s) \ne.g.> F3CS_StochasticData _test1 100\n");
        return 1;
    }
	
	int     m;
	int		l = 0;
	int		h = 0;
    int     i = 0;
    int     clock_freq = 20000000;
	unsigned short *circ_buffer1;
    circ_buffer1= (unsigned short *) calloc(3*128*algorquant,sizeof(unsigned short));
	for(i=0; i<3*128*algorquant; i++)	circ_buffer1[i] = 0;
	i=0;
	double tim = 0;
	sscanf(argv[2],"%i", &i);
	double max_tim = clock_freq*(float)i;

    // Set rate matrix here: (Currently set to a futile cycle, k 1>2>3>1 faster than k 2>1>3>2
    double k12 = 10/1000.0;
    double k13 = 1/1000.0;
    double k21 = 1/1000.0;
    double k23 = 10/1000.0;
    double k31 = 10/1000.0;
    double k32 = 1/1000.0;
    
	printf("***********************************************\nGenerates data corresponding to a 3-state Stochastic process:\n***********************************************\n");

	/*	Launch the "press q to quit" thread   */
	int rc;
	pthread_t getchar_thread;
	pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int *status;
    status = 0;
    struct paramAcq pass_to_Int;    
   	NICard.run = 1;
    
    rc = pthread_create(&getchar_thread, &attr, getchar_fxn, (void *) &pass_to_Int);
	if (rc){
		printf("ERROR; return code from pthread_create() is %d\n", rc);
		exit(-1);
	}
    
	struct write_stuff val;
	val.filenum = 0;
	val.key = 0;
	val.array = circ_buffer1;
	char	temp_char[512];
	sprintf(temp_char, "%s", argv[1]);
	val.tag = temp_char;
	pthread_t write_thread;

	int	writes = 0;
	int	failed_writes = 0;
	
	const gsl_rng_type * T;
	gsl_rng * r;
	
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	
	NICard.run = 1;
		
	/*	loop structure (for 96MB files, 1.0/13.42s)
		While time < time_max
			Write to disk: 128*algorquantum
				Check APD:	1 algorquantum (~10Hz)
					Read from card:	TUNE_BUFFER/16 samples 1/(16*64) algorquantum
	*/
	int ii = 0;
    int jj = 0;

    double uu = 0;
    double up = 0;

    int Ni[NN];  // Vector specifying which molecule is in which state
    for(jj=0; jj<NN; jj++) Ni[jj] = (int) floor(gsl_rng_uniform(r) * 3.0); // Initialize from equal distribution
    double deltaT[NN];
    for(jj=0; jj<NN; jj++) deltaT[jj] = 0;
    double upsilon[9] = {1.0,0.2,0.0,0.02,1.0,0.2,0.01,0.02,1.0};
    
    while (NICard.run == 1 && tim<max_tim)     {
		l = 0;
		for(m=0; m<128; m++)  {
			for(h=0; h<16*64; h++)	{
				for(i=0; i< 128; i++) {
                    circ_buffer1[l] = 0;
                    circ_buffer1[l+1] = 0;
                    circ_buffer1[l+2] = 0;
                    
                    // Figure out which molecules need to jump, and jump them.
                    for(jj=0; jj<NN; jj++)   {
                        deltaT[jj]--;
                        if (deltaT[jj] <= 0 )   {   // If time is up, then jump species to one of the other two states
                            uu = gsl_rng_uniform(r);
                            up = gsl_rng_uniform(r);
                            if (Ni[jj] == 0)    {
                                Ni[jj] = 1;
                                if ((uu * (k12+k13)) < k13) Ni[jj] = 2; // k is 1-indexed, Ni is zero indexed
                                deltaT[jj] = log(1/up) / (k12+k23);
                            }
                            else if (Ni[jj] == 1)    {
                                Ni[jj] = 2;
                                if ((uu * (k21+k23)) < k21) Ni[jj] = 0; // k is 1-indexed, Ni is zero indexed
                                deltaT[jj] = log(1/up) / (k21+k23);
                            }
                            else {
                                Ni[jj] = 0;
                                if ((uu * (k31+k32)) < k32) Ni[jj] = 1; // k is 1-indexed, Ni is zero indexed
                                deltaT[jj] = log(1/up) / (k31+k32);
                            }
                        }
                    }
                    
                    // Calculate random photon emission from each molecule
                    
                    for(jj=0; jj<NN; jj++){
                        circ_buffer1[l] += gsl_ran_poisson(r,upsilon[3*Ni[jj]+0]); // i.e. photon emission of each particular state is a poisson process
                        circ_buffer1[l+1] += gsl_ran_poisson(r,upsilon[3*Ni[jj]+1]);
                        circ_buffer1[l+2] += gsl_ran_poisson(r,upsilon[3*Ni[jj]+2]);
                    }
                    
                    if(tim >= 10*clock_freq && tim < 10*clock_freq+8*32*algorquant)  {
//                        circ_buffer1[l] += 1.0*gsl_ran_poisson(r,200.0); // Fake optical artifact in channel alpha at time = 10s
                        circ_buffer1[l] += 1000000.0*(l%300);
                    }
                    
                    l = l+3;
                    ii++;
				} 
			}			
			tim = tim + 16*algorquant;
		}
		
		writes++;
		val.len = l;
		val.key = 1;
		val.filenum++;
		printf("\t\tFile write started.  Files written so far = %i, Overruns = %i\n", writes, failed_writes);
		rc = pthread_create(&write_thread, &attr, write_fxn, (void *) &val);
		if (rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
		pthread_join(write_thread, (void **)status);
	}
	
    while (val.key!=0)   {
        printf("Waiting for the file I/O to complete before closing...\n");
        sleep(1);
    }
	
    printf("Small section of calculated intensities:\n");
	for(i=0; i< 128; i++) {  
		printf("I_alpha(%i) = %3.3hu \t",i,circ_buffer1[i*3]);
		printf("I_beta (%i) = %3.3hu \t",i,circ_buffer1[i*3+1]);
		printf("I_gamma(%i) = %3.3hu \n",i,circ_buffer1[i*3+2]);
	} 
	printf("        ...          \t        ...          \t        ...          \n");
	
    /*	Clean up and Exit	*/
    pthread_cancel(getchar_thread);
    NICard.run = 1;
	free(circ_buffer1);
	pthread_attr_destroy(&attr);
	gsl_rng_free (r);
	printf("Program completed.\n");
	pthread_exit(NULL);
}
