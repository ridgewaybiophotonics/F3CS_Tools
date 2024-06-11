/*
 *  Reverser.c
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
 *  Reverser.c
 *  Self-contained code to generate F3CS_Reverser.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "libF3CS.h"


//	Takes a series of reaction data and reverses the time-order.
//	The output is renamed, and can be processed using normal correlation programs
//  

int main(int argc, char *argv[])    {
    if(argc != 4)  {
        printf("Usage: F3CS_Reverser filenametag #first_file #last_file \ne.g. >F3CS_Reverser _s7 1 44 \n");
        return 1;
    }
    int i,j;
    
	FILE * fp;
	FILE * fp_r;
	char resultsfile01[80];
	char resultsfile02[80];
	
	char	temp_char[512];
	sprintf(temp_char, "%s", argv[1]);
	char *file_nametag;
	file_nametag = temp_char;
    sscanf(argv[2],"%i", &i);
   	int	first_file = i;
   	sscanf(argv[3],"%i", &i);
   	int last_file = i;
	
	for(j=first_file; j<=last_file; j++)	{		
		snprintf(resultsfile01, 50, "rxn_%i%s.w3c",j,file_nametag);
		fp = fopen(resultsfile01, "rb");
		if(fp == NULL) {
			printf("File (%s) could not be opened.  Exiting.\n", resultsfile01);
			exit(-1);
		}	
		else{
			printf("File (%s) Exists.  Good.\n", resultsfile01);
			fclose(fp);
		} 
	}
	
	
	//Begin variables from 3_xheader
    unsigned short  *reverse_interleaved;
	unsigned short  *interleaved;
	
    algorquantum = n * (int) powl(2,pmax);
    reverse_interleaved = malloc (sizeof(unsigned short)*algorquantum*4*128);
    interleaved = malloc (sizeof(unsigned short)*algorquantum*4*128);
    for(i=0; i<algorquantum*4*128; i++)	reverse_interleaved[i] = 0;
    
	int	fp_counter = 0;
	int eof = 0;
	for(j=first_file; j<=last_file; j++)	{		
		
		snprintf(resultsfile01, 50, "rxn_%i%s.w3c",j,file_nametag);
		snprintf(resultsfile02, 50, "rxn_%i%s.r.w3c",last_file-j+1,file_nametag);
		fp = fopen(resultsfile01, "rb");
		if(fp == NULL)	printf("File (%s) could not be opened.  ERROR!.\n", resultsfile01);
		fp_r = fopen(resultsfile02, "wb");
		if(fp_r == NULL)	printf("File (%s) could not be opened.  ERROR!.\n", resultsfile02);
		else{
			fp_counter = fread(interleaved, sizeof(unsigned short), 4*algorquantum*128, fp); // larger than file size to trip eof
			eof = feof(fp);
			if (fp_counter != 3*algorquantum*128 || eof == 0)	printf("(Previous) Files must be exactly 128 algorquantum in size.  ERROR!\nCurrent file size: %i, EOF? = %i\n",fp_counter/algorquantum, eof);
//			printf("Current file size: %i algorquanta (384 expected), EOF? = %i\n",fp_counter/algorquantum, eof);
			fp_counter = 3*algorquantum*128;
			for(i=0; i< fp_counter; i += 3)  {
				reverse_interleaved[fp_counter-i-3] = interleaved[i];
				reverse_interleaved[fp_counter-i-2] = interleaved[i+1];
				reverse_interleaved[fp_counter-i-1] = interleaved[i+2];
			}
			fp_counter = fwrite(reverse_interleaved, sizeof(unsigned short), 3*algorquantum*128, fp_r);
			if (fp_counter != 3*algorquantum*128)	printf("Error:  Not all data written to file (check free space on HD)\n");
			else printf("%s -> %s\n", resultsfile01,resultsfile02);
			fclose(fp);
			fclose(fp_r);
		} 
	}
	
	free(reverse_interleaved);
	free(interleaved);
    printf("Program completed.\n");
	pthread_exit(NULL);
}


