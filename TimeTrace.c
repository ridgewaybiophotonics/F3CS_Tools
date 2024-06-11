/*
 *  TimeTrace.c
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
 *  TimeTrace.c
 *  Self-contained code to generate F3CS_TimeTrace.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define	num_colors 3

unsigned short  *interleaved;

int main(int argc, char *argv[])    {
    if(argc != 5)  {
        printf("Usage: F3CS_TimeTrace filename time_quantum_ns output_hz tag\ne.g. >F3CS_TimeTrace rxn_1_test1.w3c 800 2000 test\n");
        return 1;
    }
    int i, fp_counter, filereadblock;
    float float_tc0   = 0.0f;
    float float_tc0sq = 0.0f;
    float float_tc1   = 0.0f;
    float float_tc1sq = 0.0f;
    float float_tc2   = 0.0f;
    float float_tc2sq = 0.0f;
    
    FILE * fpw;
    FILE * fpr;
    fpr = fopen(argv[1], "rb");
    if(fpr == NULL) {
        printf("File (%s) could not be opened.  Exiting.\n", argv[1]);
        return 0;
    }
    else    printf("Reading bursts from file %s...\n", argv[1]);
//	char * resultsfile1="3ctrace";
	char resultsfile1[100];
	strcpy(resultsfile1, "TimeTrace_");
	strcat(resultsfile1, argv[4]);
	strcat(resultsfile1, ".dat");
	fpw = fopen(resultsfile1, "w");
	if(fpw == NULL) {
        printf("File (%s) could not be opened.  Exiting.\n", argv[1]);
        return 0;
    }
    fprintf(fpw, "time_%s\tI_0_%s\tI_1_%s\tI_2_%s\tstdev_0_%s\tstdev_1_%s\tstdev_2_%s\n", argv[4], argv[4], argv[4], argv[4], argv[4], argv[4], argv[4]);
	
    sscanf(argv[2],"%i", &i);
    float time_quantum = i*1.0e-9;

    sscanf(argv[3],"%i", &i);
    int		out_freq = (int)i;
    if(out_freq < 1) {
        printf("Minimum out_freq = 1 Hz. (input was  %s Hz)  Exiting.\n", argv[3]);
        return 0;
    }
	float	freq_multiplier = 1.0 / time_quantum;
    printf("Calculating trace for %i kHz data binned to %i Hz\n", (int)(0.001f/time_quantum), out_freq);
    
    filereadblock = 1 / (time_quantum*out_freq);
    printf ("filereadblock = %i\n", filereadblock);
    
    interleaved = (unsigned short *) malloc (num_colors*filereadblock*sizeof(unsigned short));
    for(i=0; i<num_colors*filereadblock; i++)    interleaved[i] = 0;
    
    double  m = 0; 
    int     l = 0;
    float ftmp0, ftmp1, ftmp2;
    while (!feof (fpr))   {
        float_tc0   = 0.0f;
		float_tc0sq = 0.0f;
		float_tc1   = 0.0f;
		float_tc1sq = 0.0f;
		float_tc2   = 0.0f;
		float_tc2sq = 0.0f;
        fp_counter = fread(interleaved, sizeof(unsigned short), num_colors*filereadblock, fpr);
        if(fp_counter == num_colors*filereadblock)	{
			m = m + fp_counter;
			for(i=0; i<filereadblock; i++)  {
				ftmp0 = (float) interleaved[num_colors*i];;
				ftmp1 = (float) interleaved[num_colors*i+1];;
				ftmp2 = (float) interleaved[num_colors*i+2];;
				float_tc0 = float_tc0 + ftmp0;
				float_tc0sq = float_tc0sq + (ftmp0*ftmp0);
				float_tc1 = float_tc1 + ftmp1;
				float_tc1sq = float_tc1sq + (ftmp1*ftmp1);
				float_tc2 = float_tc2 + ftmp2;
				float_tc2sq = float_tc2sq + (ftmp2*ftmp2);
			}
			
			float_tc0 = float_tc0/(float)filereadblock;
			float_tc1 = float_tc1/(float)filereadblock;
			float_tc2 = float_tc2/(float)filereadblock;
			float_tc0sq = float_tc0sq/(float)filereadblock;
			float_tc1sq = float_tc1sq/(float)filereadblock;
			float_tc2sq = float_tc2sq/(float)filereadblock;
			float_tc0sq = powf((float_tc0sq - (float_tc0*float_tc0)),0.5f);
			float_tc1sq = powf((float_tc1sq - (float_tc1*float_tc1)),0.5f);
			float_tc2sq = powf((float_tc2sq - (float_tc2*float_tc2)),0.5f);
			float_tc0 = float_tc0 * freq_multiplier;
			float_tc1 = float_tc1 * freq_multiplier;
			float_tc2 = float_tc2 * freq_multiplier;
			float_tc0sq = float_tc0sq * freq_multiplier;
			float_tc1sq = float_tc1sq * freq_multiplier;
			float_tc2sq = float_tc2sq * freq_multiplier;
			
			printf("<0>=%2.6e +/- %2.2e, <1>=%2.6e, +/- %2.2e, <2>=%2.6e, +/- %2.2e\n", float_tc0, float_tc0sq, float_tc1, float_tc1sq, float_tc2, float_tc2sq);	
			fprintf(fpw, "%9.9f\t%9.9f\t%9.9f\t%9.9f\t%9.9f\t%9.9f\t%9.9f\n", (float)l/(float)out_freq, float_tc0, float_tc1, float_tc2, float_tc0sq, float_tc1sq, float_tc2sq);
			l++;
		}
    }
    
    fclose(fpw);
    fclose(fpr);
    free(interleaved);
    
    printf("Program completed.\n");
    return 0;
}
