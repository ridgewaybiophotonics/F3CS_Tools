/*
 *  Difference.c
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
 *  Difference.c
 *  Self-contained code for computing Delta G(tau_1,tau_2) as described 
 *  in the manuscript section "Time reversal asymmetry detection for 
 *  analysis of irreversible processes."
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "libF3CS.h"

#define	initial_read_length 	128
#define	mask_vector_l 1024*4

#ifdef n
#undef n
#define n n_GGG
#endif

#ifdef pmax
#undef pmax
#define pmax pmax_GGG
#endif

#ifdef veclen
#undef veclen
#define veclen veclen_GGG
#endif

int	max_masked(int *mask, int length, double *data)	{
	int i = 0;
	int j;
	double test = -1.0/0.0;
	for(j=0; j<length; j++)	{
		if(data[j]>test && mask[j]!=0)	{
			test = data[j];
			i = j;	
		}
	}
	return i;
}

int main(int argc, char *argv[])    {
	FILE * fp;
	int i,j,j2,k,num_curves;
		
    double *	time_array = NULL;
	double * 	intensity0_array = NULL;
	double * 	intensity1_array = NULL;
	double * 	intensity2_array = NULL;
	double * 	intensity0_stdev_array = NULL;
	double * 	intensity1_stdev_array = NULL;
	double * 	intensity2_stdev_array = NULL;
	double *	G_array = NULL;
	double * 	G_stdev_array = NULL;
	double *    G3_array = NULL;
	double *    G3_stdev_array = NULL;
	double *	GInt_array = NULL;
	double *	GInt_av = NULL;
	double *	GInt_stdev = NULL;
    int		valid_passes = 120;
    float   *garray0b0 = NULL;
    float   *garray0b1 = NULL;
    float   *garray1b1 = NULL;
    float   *garray1b2 = NULL;
    float   *garray2b0 = NULL;
    float   *garray2b2 = NULL;
    float   *garray0b1b2 = NULL;
    float   *garray1b2b0 = NULL;
    float   *garray2b0b1 = NULL;
    int     garraysize = valid_passes*veclen;
    double	tauarray[veclen];
	double 	dump[veclen*6];
	double  * triple_dump = NULL;
	double 	dump2[veclen*6];
	double  * triple_dump2 = NULL;
	double	curveheader[7];
	double 	* G0x0 = NULL;
	double 	* G1x1 = NULL;
	double 	* G2x2 = NULL;
	double 	* G0x1 = NULL;
	double 	* G1x2 = NULL;
	double 	* G2x0 = NULL;
	double 	* G0x0_prebin = NULL;
	double 	* G1x1_prebin = NULL;
	double 	* G2x2_prebin = NULL;
	double 	* G0x1_prebin = NULL;
	double 	* G1x2_prebin = NULL;
	double 	* G2x0_prebin = NULL;
	double 	* G0x1x2 = NULL;
	double 	* G1x2x0 = NULL;
	double 	* G2x0x1 = NULL;
	double 	* G0x1x2_prebin = NULL;
	double 	* G1x2x0_prebin = NULL;
	double 	* G2x0x1_prebin = NULL;
	double  * tt = NULL;
	double	* curve_time = NULL;
	double	* curve_int0 = NULL;
	double	* curve_int1 = NULL;
	double	* curve_int2 = NULL;
    double * tempa = NULL;
    double * av_array = NULL;

	if(argc != 5)  {
        printf("Usage: F3CS_Difference string first_file sigmas binwidth\ne.g. >F3CS_Difference _test1 1 3.2 8\n");
        return 1;
    }
    
    double sigmas;     
    char *endp;

    sigmas = strtod(argv[3], &endp);
    if(sigmas == 0)	{
    	printf("Input a floatingpoint number for sigmas.  Program terminated early.\n");
    	return 1;
    }
    
 	char	temp_char[512];
	sprintf(temp_char, "%s", argv[1]);
 	char *file_nametag;
	file_nametag = temp_char;
    sscanf(argv[2],"%i", &i);
   	int	first_file = i;
   	int	file_delta = num_curves_per_thread;
    sscanf(argv[4],"%i", &i);
   	int	bw = i;
    if (mask_vector_l < bw)	{
    	printf("Recompile with defined mask_vector_l > %i\nProgram terminating\n", bw);
    	return 1;
    }
        
	char Garrayfile[80];
	snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, first_file);
	int	file_length = 0;
	int	number_of_scans_total = 0;
	int	number_of_scans_prev = 0;
	int	*scans_per_file;
	scans_per_file = (int*) malloc(mask_vector_l*512*2*sizeof(int)); // small database of where data are located, which file and scannumber
	for(i=0; i<mask_vector_l*512*2; i++)	scans_per_file[i] = 0;
	int	data_quantum = 7+6*veclen+3*veclen*veclen;
	
    int file_index = first_file;
    //Check existance of files...
	do	{		
		snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, file_index);
		fp = fopen(Garrayfile, "rb");
		if(fp != NULL) {
			fseek(fp, 0L, SEEK_END);
			file_length = ftell(fp) / sizeof(double);
			number_of_scans_prev = number_of_scans_total;
			number_of_scans_total += (file_length - 1 -veclen) / data_quantum;
			for(i=number_of_scans_prev; i<number_of_scans_total; i++)	{
				scans_per_file[i*2] = file_index;
				scans_per_file[i*2+1] = i-number_of_scans_prev;
			}
			fclose(fp);
			file_index += file_delta;
		}
	}
	while (fp != NULL);
    int high_file = file_index - file_delta;
	
	// Curve N --> File X
    garray0b0 = (float *)malloc(garraysize*sizeof(float));
    garray0b1 = (float *)malloc(garraysize*sizeof(float));
    garray1b1 = (float *)malloc(garraysize*sizeof(float));
    garray1b2 = (float *)malloc(garraysize*sizeof(float));
    garray2b0 = (float *)malloc(garraysize*sizeof(float));
    garray2b2 = (float *)malloc(garraysize*sizeof(float));
    garray0b1b2 = (float *)malloc(garraysize*veclen*sizeof(float));
    garray1b2b0 = (float *)malloc(garraysize*veclen*sizeof(float));
    garray2b0b1 = (float *)malloc(garraysize*veclen*sizeof(float));
    
    for(i=0; i<garraysize; i++) garray0b0[i] = 1.0*i;
    for(i=0; i<garraysize; i++) garray1b1[i] = 2.0*i;
    for(i=0; i<garraysize; i++) garray2b2[i] = 3.0*i;
    for(i=0; i<garraysize; i++) garray0b1[i] = 4.0*i;
    for(i=0; i<garraysize; i++) garray1b2[i] = 5.0*i;
    for(i=0; i<garraysize; i++) garray2b0[i] = 6.0*i;
    for(i=0; i<garraysize*veclen; i++) garray0b1b2[i] = 4.0*i;
    for(i=0; i<garraysize*veclen; i++) garray1b2b0[i] = 4.0*i;
    for(i=0; i<garraysize*veclen; i++) garray2b0b1[i] = 4.0*i;
    
    double *d;
	const int	d_length = 7+6*veclen;
	d = malloc(sizeof(double)*d_length);
	
	//Open file, populate array	
	triple_dump = malloc(sizeof(double)*3*veclen*veclen);
	triple_dump2 = malloc(sizeof(double)*3*veclen*veclen);
	for(i=0; i< 6*veclen; i++)		dump[i] = -1.0;
	for(i=0; i< 6*veclen; i++)		dump2[i] = -3.0;
	for(i=0; i< 3*veclen*veclen; i++)		triple_dump[i] = -2.0;
	for(i=0; i< 3*veclen*veclen; i++)		triple_dump2[i] = -4.0;
	G0x0 = malloc(sizeof(double)*veclen*initial_read_length);
	for(i=0; i< veclen*initial_read_length; i++)	G0x0[i] = 0.0;
	G1x1 = malloc(sizeof(double)*veclen*initial_read_length);
	for(i=0; i< veclen*initial_read_length; i++)	G1x1[i] = 0.0;
	G2x2 = malloc(sizeof(double)*veclen*initial_read_length);
	for(i=0; i< veclen*initial_read_length; i++)	G2x2[i] = 0.0;
	G0x1 = malloc(sizeof(double)*veclen*initial_read_length);
	for(i=0; i< veclen*initial_read_length; i++)	G0x1[i] = 0.0;
	G1x2 = malloc(sizeof(double)*veclen*initial_read_length);
	for(i=0; i< veclen*initial_read_length; i++)	G1x2[i] = 0.0;
	G2x0 = malloc(sizeof(double)*veclen*initial_read_length);
	for(i=0; i< veclen*initial_read_length; i++)	G2x0[i] = 0.0;
	
	G0x0_prebin = malloc(sizeof(double)*veclen);
	for(i=0; i< veclen; i++)	G0x0_prebin[i] = 0.0;
	G1x1_prebin = malloc(sizeof(double)*veclen);
	for(i=0; i< veclen; i++)	G1x1_prebin[i] = 0.0;
	G2x2_prebin = malloc(sizeof(double)*veclen);
	for(i=0; i< veclen; i++)	G2x2_prebin[i] = 0.0;
	G0x1_prebin = malloc(sizeof(double)*veclen);
	for(i=0; i< veclen; i++)	G0x1_prebin[i] = 0.0;
	G1x2_prebin = malloc(sizeof(double)*veclen);
	for(i=0; i< veclen; i++)	G1x2_prebin[i] = 0.0;
	G2x0_prebin = malloc(sizeof(double)*veclen);
	for(i=0; i< veclen; i++)	G2x0_prebin[i] = 0.0;
	
	G0x1x2 = malloc(sizeof(double)*veclen*veclen*initial_read_length);
    if(G0x1x2 == NULL)  {printf("Malloc failed, likely not enough memory.  Exiting\n"); goto cleanupandexit;}
	for(i=0; i< veclen*veclen*initial_read_length; i++)	G0x1x2[i] = 0.0;
	G1x2x0 = malloc(sizeof(double)*veclen*veclen*initial_read_length);
    if(G1x2x0 == NULL)  {printf("Malloc failed, likely not enough memory.  Exiting\n"); goto cleanupandexit;}
	for(i=0; i< veclen*veclen*initial_read_length; i++)	G1x2x0[i] = 0.0;
	G2x0x1 = malloc(sizeof(double)*veclen*veclen*initial_read_length);
    if(G2x0x1 == NULL)  {printf("Malloc failed, likely not enough memory.  Exiting\n"); goto cleanupandexit;}
	for(i=0; i< veclen*veclen*initial_read_length; i++)	G2x0x1[i] = 0.0;
	
	G0x1x2_prebin = malloc(sizeof(double)*veclen*veclen);
    if(G0x1x2_prebin == NULL)  {printf("Malloc failed, likely not enough memory.  Exiting\n"); goto cleanupandexit;}
	for(i=0; i< veclen*veclen; i++)	G0x1x2_prebin[i] = 0.0;
	G1x2x0_prebin = malloc(sizeof(double)*veclen*veclen);
    if(G1x2x0_prebin == NULL)  {printf("Malloc failed, likely not enough memory.  Exiting\n"); goto cleanupandexit;}
	for(i=0; i< veclen*veclen; i++)	G1x2x0_prebin[i] = 0.0;
	G2x0x1_prebin = malloc(sizeof(double)*veclen*veclen);
    if(G2x0x1_prebin == NULL)  {printf("Malloc failed, likely not enough memory.  Exiting\n"); goto cleanupandexit;}
	for(i=0; i< veclen*veclen; i++)	G2x0x1_prebin[i] = 0.0;
	
	tt   = malloc(sizeof(double)*6*initial_read_length*2);
	curve_time = (double *)malloc(sizeof(double)*initial_read_length*2);
	for(i=0; i<initial_read_length*2; i++)	curve_time[i] = 0;
	curve_int0 = (double *)malloc(sizeof(double)*initial_read_length*2);
	for(i=0; i<initial_read_length*2; i++)	curve_int0[i] = 0;
	curve_int1 = (double *)malloc(sizeof(double)*initial_read_length*2);
	for(i=0; i<initial_read_length*2; i++)	curve_int1[i] = 0;
	curve_int2 = (double *)malloc(sizeof(double)*initial_read_length*2);
	for(i=0; i<initial_read_length*2; i++)	curve_int2[i] = 0;

	double	time_offset;

	for(i=0; i<veclen; i++) tauarray[i] = 0.0f;
	for(i=0; i<veclen*6; i++)	d[i] = 0.0f;
	
	FILE * fp2;	
    file_index = first_file;
    snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, file_index);
    fp = NULL;
	fp = fopen(Garrayfile, "rb");
    if (fp==NULL)	{   printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1); }
    else                printf("Opening file (%s)\n", Garrayfile);

	snprintf(Garrayfile, 50, "G_array%s.r_%i.bin2",file_nametag, scans_per_file[2*(number_of_scans_total-1-0)]);
    int fp2_filenumber = scans_per_file[2*(number_of_scans_total-1-0)];
    fp2 = NULL;
	fp2 = fopen(Garrayfile, "rb");
	if (fp2==NULL)	{printf(".r file %s not found.  Exiting.\n", Garrayfile); exit(-1);}
    else			printf("Opening .r file (%s)\n", Garrayfile);

	j = fread(&time_offset, sizeof(double), 1, fp);
	j2 = fread(&time_offset, sizeof(double), 1, fp2);
	j = fread(tauarray, sizeof(double), veclen, fp);
	j2 = fread(tauarray, sizeof(double), veclen, fp2);
	
	int prebin_counter=0;
	file_index = first_file;
	num_curves=0;
	int num_individual_curves = 0;
	do {
		//	Time	Int_0	Int_1	Int_2	IntSd_0	IntSD_1	IntSD_2				(7)
        
        // if the .r file needed isn't the one that is open, then close the old file and open the new file
        if(fp2_filenumber != scans_per_file[2*(number_of_scans_total-1-num_individual_curves)]) {
            fclose(fp2); fp2 = NULL;
            snprintf(Garrayfile, 50, "G_array%s.r_%i.bin2",file_nametag, scans_per_file[2*(number_of_scans_total-1-num_individual_curves)]);
			printf("Opening .r File %s\n", Garrayfile);
			fp2 = fopen(Garrayfile, "rb");
			if (fp2==NULL)	{printf(".r file %s not found.  Exiting.\n", Garrayfile); goto cleanupandexit;}

            fp2_filenumber = scans_per_file[2*(number_of_scans_total-1-num_individual_curves)];
        }
        
		if(fseek(fp2, (veclen+1)*sizeof(double), SEEK_SET)!=0)	printf("Fseek1 Error!\n");
		
		j = fread(curveheader, sizeof(double), 7, fp);
		if (j!=7 && j != 0)	printf("****Something strange going on - file sizes are not as expected****\n");
		curve_time[num_curves] = curveheader[0];
		curve_int0[num_curves] = curveheader[1];
		curve_int1[num_curves] = curveheader[2];
		curve_int2[num_curves] = curveheader[3];
		
        if(fseek(fp2, data_quantum*scans_per_file[2*(number_of_scans_total-1-num_individual_curves)+1]*sizeof(double), SEEK_CUR)!=0 && j== (6*veclen+3*veclen*veclen))	printf("Fseek2.1 Error!\n");
		j2 = fread(curveheader, sizeof(double), 7, fp2);
		if (j2!=7 && j2 != 0)	printf("****Something strange is going on - .r file sizes are not as expected****\n");
		curve_time[num_curves+initial_read_length] = curveheader[0];
		curve_int0[num_curves+initial_read_length] = curveheader[1];
		curve_int1[num_curves+initial_read_length] = curveheader[2];
		curve_int2[num_curves+initial_read_length] = curveheader[3];
		
		j =  fread(dump, sizeof(double), 6*veclen, fp);
		j += fread(triple_dump, sizeof(double), 3*veclen*veclen, fp);
		
		j2 =  fread(dump2, sizeof(double), 6*veclen, fp2);
		j2 += fread(triple_dump2, sizeof(double), 3*veclen*veclen, fp2);
		
		if (j!= (6*veclen+3*veclen*veclen) && file_index < high_file)	{
			if(j != 0)	printf("Left-over bytes = %i\n", j);
			fclose(fp);
			
			file_index += file_delta;
			snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, file_index);
			fp = fopen(Garrayfile, "rb");
			if (fp==NULL)	{printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1);}
            else 	printf("Opening file (%s)\n", Garrayfile);
			
			j = fread(&time_offset, sizeof(double), 1, fp);
			fseek(fp, (veclen)*sizeof(double),SEEK_CUR);
			
			j = fread(curveheader, sizeof(double), 7, fp);			
			if (j!=7 && j != 0)	printf("****Something strange going on - file sizes are not as expected****\n");
			curve_time[num_curves] = curveheader[0];
			curve_int0[num_curves] = curveheader[1];
			curve_int1[num_curves] = curveheader[2];
			curve_int2[num_curves] = curveheader[3];
						
			j =  fread(dump, sizeof(double), 6*veclen, fp);
			j += fread(triple_dump, sizeof(double), 3*veclen*veclen, fp);
		}	
		if(j== (6*veclen+3*veclen*veclen))	{	// Write data to 8X prebinning arrays.  Average them.
			for(i=0; i< veclen; i++)	{
				G0x0_prebin[i] += dump2[i+0*veclen]-dump[i+0*veclen];
				G1x1_prebin[i] += dump2[i+1*veclen]-dump[i+1*veclen];
				G2x2_prebin[i] += dump2[i+2*veclen]-dump[i+2*veclen];
				G0x1_prebin[i] += dump2[i+3*veclen]-dump[i+3*veclen];
				G1x2_prebin[i] += dump2[i+4*veclen]-dump[i+4*veclen];
				G2x0_prebin[i] += dump2[i+5*veclen]-dump[i+5*veclen];
			}
			
			for(i=0; i< veclen*veclen; i++)	{
				G0x1x2_prebin[i] += triple_dump2[i+0*veclen*veclen]-triple_dump[i+0*veclen*veclen];
				G1x2x0_prebin[i] += triple_dump2[i+1*veclen*veclen]-triple_dump[i+1*veclen*veclen];
				G2x0x1_prebin[i] += triple_dump2[i+2*veclen*veclen]-triple_dump[i+2*veclen*veclen];
			}
			prebin_counter++;
		}
		if(prebin_counter == 8){
			for(i=0; i< veclen; i++)	{
				G0x0[i+num_curves*veclen] = G0x0_prebin[i]/8.0;
				G1x1[i+num_curves*veclen] = G1x1_prebin[i]/8.0;
				G2x2[i+num_curves*veclen] = G2x2_prebin[i]/8.0;
				G0x1[i+num_curves*veclen] = G0x1_prebin[i]/8.0;
				G1x2[i+num_curves*veclen] = G1x2_prebin[i]/8.0;
				G2x0[i+num_curves*veclen] = G2x0_prebin[i]/8.0;
			}
			
			for(i=0; i< veclen*veclen; i++)	{
				G0x1x2[i+num_curves*veclen*veclen] = G0x1x2_prebin[i]/8.0;
				G1x2x0[i+num_curves*veclen*veclen] = G1x2x0_prebin[i]/8.0;
				G2x0x1[i+num_curves*veclen*veclen] = G2x0x1_prebin[i]/8.0;
			}
			
			// Set prebin arrays to 0.
			for(i=0; i< veclen; i++)	{
				G0x0_prebin[i] = 0.0;
				G1x1_prebin[i] = 0.0;
				G2x2_prebin[i] = 0.0;
				G0x1_prebin[i] = 0.0;
				G1x2_prebin[i] = 0.0;
				G2x0_prebin[i] = 0.0;
			}
			
			for(i=0; i< veclen*veclen; i++)	{
				G0x1x2_prebin[i] = 0.0;
				G1x2x0_prebin[i] = 0.0;
				G2x0x1_prebin[i] = 0.0;
			}

			prebin_counter=0;
			num_curves++;
		}
		num_individual_curves++;
	}	while (j==(6*veclen+3*veclen*veclen) && num_curves<initial_read_length && number_of_scans_total >num_individual_curves);
	
	if(fp!= NULL)	{fclose(fp); fp = NULL;}
	if(fp2!=NULL)	{fclose(fp2); fp2 = NULL;}
	
	if(num_curves==initial_read_length)	printf("*******\nDataset has more than initial_read_length=%i curves.\nRecompile this program with a larger initial_read_length value\n********\n",initial_read_length);
	printf("Total number of curves in the data = %i\n", num_curves);
	
    if(bw<0)	{
		bw = num_curves;
		if(bw <= mask_vector_l)		printf("Binning matched to number of curves (%i)\n", bw);
		else {
			bw = mask_vector_l;
			printf("Binning set to maximum (%i).\n If a larger bin size is desired, recompile with a larger mask_vector_l\n", bw);
		}
	}
	if(bw>num_curves) {
		bw = num_curves;
		printf("Binning decreased to number of curves (%i)\n", bw);
	}
	
	float	bw_mask[mask_vector_l];
	for(i=0; i<mask_vector_l; i++)	bw_mask[i] = 1.0;
	
	int	num_slices = num_curves / bw;

	time_array = (double *) malloc(num_slices*sizeof(double));
	intensity0_array = (double *) malloc(num_slices*sizeof(double));
	intensity1_array = (double *) malloc(num_slices*sizeof(double));
	intensity2_array = (double *) malloc(num_slices*sizeof(double));
	intensity0_stdev_array = (double *) malloc(num_slices*sizeof(double));
	intensity1_stdev_array = (double *) malloc(num_slices*sizeof(double));
	intensity2_stdev_array = (double *) malloc(num_slices*sizeof(double));
	G_array = (double *) malloc(num_slices*veclen*6*sizeof(double));
	if(G_array == NULL)	printf("Malloc Failed!\n");
	G_stdev_array = (double *) malloc(num_slices*veclen*6*sizeof(double));
	G3_array = (double *) malloc(num_slices*veclen*veclen*3*sizeof(double));
	if(G3_array == NULL)	printf("Malloc Failed!\n");
	G3_stdev_array = (double *) malloc(num_slices*veclen*veclen*3*sizeof(double));
	if(G3_stdev_array == NULL)	printf("Malloc Failed!\n");
	GInt_array = (double *) malloc(bw*6*sizeof(double));
	GInt_av = (double *) malloc(6*sizeof(double));
	GInt_stdev = (double *) malloc(6*sizeof(double));
	
	for(i=0; i< num_slices; i++)   time_array[i] = 0.0;
	for(i=0; i< num_slices; i++)   intensity0_array[i] = 0.0;
	for(i=0; i< num_slices; i++)   intensity1_array[i] = 0.0;
	for(i=0; i< num_slices; i++)   intensity2_array[i] = 0.0;
	for(i=0; i< num_slices; i++)   intensity0_stdev_array[i] = 0.0;
	for(i=0; i< num_slices; i++)   intensity1_stdev_array[i] = 0.0;
	for(i=0; i< num_slices; i++)   intensity2_stdev_array[i] = 0.0;
	for(i=0; i< num_slices*veclen*6; i++)   G_array[i] = 0.0;
	for(i=0; i< num_slices*veclen*6; i++)   G_stdev_array[i] = 0.0;
	for(i=0; i< num_slices*veclen*veclen*3; i++)   G3_array[i] = 0.0;
	for(i=0; i< num_slices*veclen*veclen*3; i++)   G3_stdev_array[i] = 0.0;
	for(i=0; i< bw*6; i++)   	GInt_array[i] = 0.0;
	for(i=0; i< 6; i++)   		GInt_av[i] = 0.0;
	for(i=0; i< 6; i++)   		GInt_stdev[i] = 0.0;

	int	slices = 0;
	int	worst_curve = 0;
	double worst_value = 0.0;
	double worst_temp = 0.0;
	double remaining_bw_curves = 0.0;
	const double rootbw = 1.0 / pow((double) bw, 0.5);
	tempa = (double *) malloc(6*mask_vector_l*sizeof(double));
	av_array = (double *) malloc((6+1)*bw*sizeof(double));
	for(i=0; i< (6+1)*bw; i++)   av_array[i] = 0.0;
	for(slices = 0; slices < num_slices; slices++)	{	//Slice up dataset into bw-sized chunks, do math
        for(i=0; i<mask_vector_l; i++)	bw_mask[i] = 1.0;
        //For each curve, if any G[tau] values are NaN, discard the entire curve
        for(i=0; i<bw; i++)	{
            for(k=0; k<veclen; k++)	{		//Integrate each curve | tau = [10us-100ms]
                if(isnan(G0x0[(slices*bw+i)*veclen+k])) bw_mask[i] = 0;
                if(isnan(G1x1[(slices*bw+i)*veclen+k])) bw_mask[i] = 0;
                if(isnan(G2x2[(slices*bw+i)*veclen+k])) bw_mask[i] = 0;
                if(isnan(G0x1[(slices*bw+i)*veclen+k])) bw_mask[i] = 0;
                if(isnan(G1x2[(slices*bw+i)*veclen+k])) bw_mask[i] = 0;
                if(isnan(G2x0[(slices*bw+i)*veclen+k])) bw_mask[i] = 0;
            }
            if(bw_mask[i] == 0) printf("Curve %i contains NaN values and was discarded\n", slices*bw+i);
        }
        
		for(i=0; i<bw; i++) remaining_bw_curves += bw_mask[i]*1.0;		//Find Outliers
		
		//Find Outliers
		do{
			remaining_bw_curves = 0.0;
			for(i=0; i<bw; i++) remaining_bw_curves += bw_mask[i]*1.0;
			for(i=0; i< 6*bw; i++)	GInt_array[i] = 0.0;
			for(i=0; i<bw; i++)	{
				for(k=12; k<veclen; k++)	{		//Integrate each curve | tau = [10us-100ms]
					GInt_array[i+0*bw] += bw_mask[i]*G0x0[(slices*bw+i)*veclen+k];
					GInt_array[i+1*bw] += bw_mask[i]*G1x1[(slices*bw+i)*veclen+k];
					GInt_array[i+2*bw] += bw_mask[i]*G2x2[(slices*bw+i)*veclen+k];
					GInt_array[i+3*bw] += bw_mask[i]*G0x1[(slices*bw+i)*veclen+k];
					GInt_array[i+4*bw] += bw_mask[i]*G1x2[(slices*bw+i)*veclen+k];
					GInt_array[i+5*bw] += bw_mask[i]*G2x0[(slices*bw+i)*veclen+k];
				}
			}
			for(i=0; i< 6*bw; i++)	GInt_array[i] = GInt_array[i] / (veclen-12.0);
			
			//  Average each GaXb
			for(j=0; j< 6; j++)   {	 
				GInt_av[j] = 0.0;
				for(i=0; i<bw; i++)	{
					GInt_av[j] += bw_mask[j]*GInt_array[i+j*bw];
				}
				GInt_av[j] = GInt_av[j] / (double) remaining_bw_curves;
			}
			
			//	Calc. stdev for each GaXb	
			for(j=0; j< 6; j++)   {	 
				GInt_stdev[j] = 0.0;
				for(i=0; i<bw; i++)	{
					GInt_stdev[j] += bw_mask[i]*pow(GInt_av[j] - GInt_array[i+j*bw], 2.0);
				}
				GInt_stdev[j] = GInt_stdev[j] / (remaining_bw_curves - 1.0);
				GInt_stdev[j] = pow(GInt_stdev[j], 0.5);
			}
			
			//	Find largest difference^2
			worst_curve = 0;
			worst_value = 0.0;
			for(j=0; j< 6; j++)   {	
				for(i=0; i<bw; i++)	{
					worst_temp = pow(GInt_av[j] - GInt_array[i+j*bw], 2.0) / pow(GInt_stdev[j], 2.0);
					if (bw_mask[i])	{
						if (worst_temp > worst_value)	{
							worst_value = worst_temp;
							worst_curve = i;
						}
					}	
				}
			}
			
			//	If Difference is > threshold, omit, repeat loop
			worst_value = pow(worst_value, 0.5);
			if (worst_value > sigmas) {
                printf("Outlier found with error %f x sigma\n", worst_value);
                bw_mask[worst_curve] = 0;
            }
		}	while (worst_value > sigmas);
		
		
		//  Do Math on the Remaining Curves
		time_array[slices] = 0.0;
		intensity0_array[slices] = 0.0;
		intensity1_array[slices] = 0.0;
		intensity2_array[slices] = 0.0;
		intensity0_stdev_array[slices] = 0.0;
		intensity1_stdev_array[slices] = 0.0;
		intensity2_stdev_array[slices] = 0.0;
		
		for(j=0; j<bw; j++)	time_array[slices] += curve_time[slices*bw+j];
		time_array[slices] = time_array[slices] / bw;
		
        if(bw-remaining_bw_curves==1) printf("%.0f curve was discarded from this bin\n", bw-remaining_bw_curves);
        else  printf("%.0f curves were discarded from this bin\n", bw-remaining_bw_curves);
		
		for(j=0; j<bw; j++)	intensity0_array[slices] += bw_mask[j]*curve_int0[slices*bw+j];
		intensity0_array[slices] = intensity0_array[slices] / remaining_bw_curves;
		for(j=0; j<bw; j++)	intensity1_array[slices] += bw_mask[j]*curve_int1[slices*bw+j];
		intensity1_array[slices] = intensity1_array[slices] / remaining_bw_curves;
		for(j=0; j<bw; j++)	intensity2_array[slices] += bw_mask[j]*curve_int2[slices*bw+j];
		intensity2_array[slices] = intensity2_array[slices] / remaining_bw_curves;
		
		for(j=0; j<bw; j++)	intensity0_stdev_array[slices] += bw_mask[j]*pow(curve_int0[slices*bw+j] - intensity0_array[slices], 2.0);
		for(j=0; j<bw; j++)	intensity1_stdev_array[slices] += bw_mask[j]*pow(curve_int1[slices*bw+j] - intensity1_array[slices], 2.0);
		for(j=0; j<bw; j++)	intensity2_stdev_array[slices] += bw_mask[j]*pow(curve_int2[slices*bw+j] - intensity2_array[slices], 2.0);
		
		intensity0_stdev_array[slices] = intensity0_stdev_array[slices] / (remaining_bw_curves-1.0);
		intensity1_stdev_array[slices] = intensity1_stdev_array[slices] / (remaining_bw_curves-1.0);
		intensity2_stdev_array[slices] = intensity2_stdev_array[slices] / (remaining_bw_curves-1.0);
		
		intensity0_stdev_array[slices] = rootbw * pow(intensity0_stdev_array[slices], 0.5);
		intensity1_stdev_array[slices] = rootbw * pow(intensity1_stdev_array[slices], 0.5);
		intensity2_stdev_array[slices] = rootbw * pow(intensity2_stdev_array[slices], 0.5);
				
		for(k=0; k<veclen; k++)	{
			G_array[slices*veclen*6 + 0*veclen + k] = 0.0;
			G_array[slices*veclen*6 + 1*veclen + k] = 0.0;
			G_array[slices*veclen*6 + 2*veclen + k] = 0.0;
			G_array[slices*veclen*6 + 3*veclen + k] = 0.0;
			G_array[slices*veclen*6 + 4*veclen + k] = 0.0;
			G_array[slices*veclen*6 + 5*veclen + k] = 0.0;
			G_stdev_array[slices*veclen*6 + 0*veclen + k] = 0.0;
			G_stdev_array[slices*veclen*6 + 1*veclen + k] = 0.0;
			G_stdev_array[slices*veclen*6 + 2*veclen + k] = 0.0;
			G_stdev_array[slices*veclen*6 + 3*veclen + k] = 0.0;
			G_stdev_array[slices*veclen*6 + 4*veclen + k] = 0.0;
			G_stdev_array[slices*veclen*6 + 5*veclen + k] = 0.0;
			for(j=0; j<bw; j++) G_array[slices*veclen*6 + 0*veclen + k] += bw_mask[j]*G0x0[(slices*bw+j)*veclen+k];
			for(j=0; j<bw; j++) G_array[slices*veclen*6 + 1*veclen + k] += bw_mask[j]*G1x1[(slices*bw+j)*veclen+k];
			for(j=0; j<bw; j++) G_array[slices*veclen*6 + 2*veclen + k] += bw_mask[j]*G2x2[(slices*bw+j)*veclen+k];
			for(j=0; j<bw; j++) G_array[slices*veclen*6 + 3*veclen + k] += bw_mask[j]*G0x1[(slices*bw+j)*veclen+k];
			for(j=0; j<bw; j++) G_array[slices*veclen*6 + 4*veclen + k] += bw_mask[j]*G1x2[(slices*bw+j)*veclen+k];
			for(j=0; j<bw; j++) G_array[slices*veclen*6 + 5*veclen + k] += bw_mask[j]*G2x0[(slices*bw+j)*veclen+k];
			G_array[slices*veclen*6 + 0*veclen + k] = G_array[slices*veclen*6 + 0*veclen + k] / remaining_bw_curves;
			G_array[slices*veclen*6 + 1*veclen + k] = G_array[slices*veclen*6 + 1*veclen + k] / remaining_bw_curves;
			G_array[slices*veclen*6 + 2*veclen + k] = G_array[slices*veclen*6 + 2*veclen + k] / remaining_bw_curves;
			G_array[slices*veclen*6 + 3*veclen + k] = G_array[slices*veclen*6 + 3*veclen + k] / remaining_bw_curves;
			G_array[slices*veclen*6 + 4*veclen + k] = G_array[slices*veclen*6 + 4*veclen + k] / remaining_bw_curves;
			G_array[slices*veclen*6 + 5*veclen + k] = G_array[slices*veclen*6 + 5*veclen + k] / remaining_bw_curves;
			
			for(j=0; j<bw; j++) G_stdev_array[slices*veclen*6 + 0*veclen + k] += bw_mask[j]*pow(G0x0[(slices*bw+j)*veclen+k]-G_array[slices*veclen*6 + 0*veclen + k], 2.0);
			for(j=0; j<bw; j++) G_stdev_array[slices*veclen*6 + 1*veclen + k] += bw_mask[j]*pow(G1x1[(slices*bw+j)*veclen+k]-G_array[slices*veclen*6 + 1*veclen + k], 2.0);
			for(j=0; j<bw; j++) G_stdev_array[slices*veclen*6 + 2*veclen + k] += bw_mask[j]*pow(G2x2[(slices*bw+j)*veclen+k]-G_array[slices*veclen*6 + 2*veclen + k], 2.0);
			for(j=0; j<bw; j++) G_stdev_array[slices*veclen*6 + 3*veclen + k] += bw_mask[j]*pow(G0x1[(slices*bw+j)*veclen+k]-G_array[slices*veclen*6 + 3*veclen + k], 2.0);
			for(j=0; j<bw; j++) G_stdev_array[slices*veclen*6 + 4*veclen + k] += bw_mask[j]*pow(G1x2[(slices*bw+j)*veclen+k]-G_array[slices*veclen*6 + 4*veclen + k], 2.0);
			for(j=0; j<bw; j++) G_stdev_array[slices*veclen*6 + 5*veclen + k] += bw_mask[j]*pow(G2x0[(slices*bw+j)*veclen+k]-G_array[slices*veclen*6 + 5*veclen + k], 2.0);
			G_stdev_array[slices*veclen*6 + 0*veclen + k] = G_stdev_array[slices*veclen*6 + 0*veclen + k] / (remaining_bw_curves-1.0);
			G_stdev_array[slices*veclen*6 + 1*veclen + k] = G_stdev_array[slices*veclen*6 + 1*veclen + k] / (remaining_bw_curves-1.0);
			G_stdev_array[slices*veclen*6 + 2*veclen + k] = G_stdev_array[slices*veclen*6 + 2*veclen + k] / (remaining_bw_curves-1.0);
			G_stdev_array[slices*veclen*6 + 3*veclen + k] = G_stdev_array[slices*veclen*6 + 3*veclen + k] / (remaining_bw_curves-1.0);
			G_stdev_array[slices*veclen*6 + 4*veclen + k] = G_stdev_array[slices*veclen*6 + 4*veclen + k] / (remaining_bw_curves-1.0);
			G_stdev_array[slices*veclen*6 + 5*veclen + k] = G_stdev_array[slices*veclen*6 + 5*veclen + k] / (remaining_bw_curves-1.0);
			G_stdev_array[slices*veclen*6 + 0*veclen + k] = rootbw * pow(G_stdev_array[slices*veclen*6 + 0*veclen + k], 0.5);
			G_stdev_array[slices*veclen*6 + 1*veclen + k] = rootbw * pow(G_stdev_array[slices*veclen*6 + 1*veclen + k], 0.5);
			G_stdev_array[slices*veclen*6 + 2*veclen + k] = rootbw * pow(G_stdev_array[slices*veclen*6 + 2*veclen + k], 0.5);
			G_stdev_array[slices*veclen*6 + 3*veclen + k] = rootbw * pow(G_stdev_array[slices*veclen*6 + 3*veclen + k], 0.5);
			G_stdev_array[slices*veclen*6 + 4*veclen + k] = rootbw * pow(G_stdev_array[slices*veclen*6 + 4*veclen + k], 0.5);
			G_stdev_array[slices*veclen*6 + 5*veclen + k] = rootbw * pow(G_stdev_array[slices*veclen*6 + 5*veclen + k], 0.5);
		}
		
		for(k=0; k<veclen*veclen; k++)	{
			G3_array[3*slices*veclen*veclen + 0*veclen*veclen + k] = 0.0;
			G3_array[3*slices*veclen*veclen + 1*veclen*veclen + k] = 0.0;
			G3_array[3*slices*veclen*veclen + 2*veclen*veclen + k] = 0.0;
			G3_stdev_array[3*slices*veclen*veclen + 0*veclen*veclen + k] = 0.0;
			G3_stdev_array[3*slices*veclen*veclen + 1*veclen*veclen + k] = 0.0;
			G3_stdev_array[3*slices*veclen*veclen + 2*veclen*veclen + k] = 0.0;
			
			for(j=0; j<bw; j++) G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k] += bw_mask[j]*G0x1x2[(slices*bw+j)*veclen*veclen+k];
			for(j=0; j<bw; j++) G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k] += bw_mask[j]*G1x2x0[(slices*bw+j)*veclen*veclen+k];
			for(j=0; j<bw; j++) G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k] += bw_mask[j]*G2x0x1[(slices*bw+j)*veclen*veclen+k];
			
			G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k] / remaining_bw_curves;
			G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k] / remaining_bw_curves;
			G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k] / remaining_bw_curves;
			
            for(j=0; j<bw; j++) G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] += bw_mask[j]*pow(G0x1x2[(slices*bw+j)*veclen*veclen+k]-G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k], 2.0);
            for(j=0; j<bw; j++) G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] += bw_mask[j]*pow(G1x2x0[(slices*bw+j)*veclen*veclen+k]-G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k], 2.0);
            for(j=0; j<bw; j++) G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] += bw_mask[j]*pow(G2x0x1[(slices*bw+j)*veclen*veclen+k]-G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k], 2.0);
			
			G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] / (remaining_bw_curves-1.0);
			G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] / (remaining_bw_curves-1.0);
			G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] / (remaining_bw_curves-1.0);

			G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = rootbw * pow(G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k], 0.5);
			G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = rootbw * pow(G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k], 0.5);
			G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = rootbw * pow(G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k], 0.5);
		}
		
		for(i=0; i<bw; i++) {
			tempa[i+0*bw] = 0.0;
			for(k=37; k<143; k++) tempa[i+0*bw] += G0x0[(slices*bw+i)*veclen+k];
			tempa[i+1*bw] = 0.0;
			for(k=37; k<143; k++) tempa[i+1*bw] += G1x1[(slices*bw+i)*veclen+k];
			tempa[i+2*bw] = 0.0;
			for(k=37; k<143; k++) tempa[i+2*bw] += G2x2[(slices*bw+i)*veclen+k];
			tempa[i+3*bw] = 0.0;
			for(k=37; k<143; k++) tempa[i+3*bw] += G0x1[(slices*bw+i)*veclen+k];
			tempa[i+4*bw] = 0.0;
			for(k=37; k<143; k++) tempa[i+4*bw] += G1x2[(slices*bw+i)*veclen+k];
			tempa[i+5*bw] = 0.0;
			for(k=37; k<143; k++) tempa[i+5*bw] += G2x0[(slices*bw+i)*veclen+k];
		}
    }
	
	char resultsfile90[80];	//Open file to write out data
	snprintf(resultsfile90, 50, "clean_Gs%s.dat",file_nametag);
	FILE * fpG;
    fpG = fopen(resultsfile90, "w");
	
	char resultsfile100[80];	//Open file to write out data
	snprintf(resultsfile100, 50, "clean_GGGs%s.dat",file_nametag);
	FILE * fpGGG;
    fpGGG = fopen(resultsfile100, "w");
    
    char resultsfile110[80];	//Open file to write out data
	snprintf(resultsfile110, 50, "samples_GGGs%s.dat",file_nametag);
	FILE * fpsGGG;
    fpsGGG = fopen(resultsfile110, "w");

    fprintf(fpG,"tau%s", argv[1]);
    for(slices = 0; slices < num_slices; slices++)	{
		fprintf(fpG,"\tG0x0_%.2f%s\tG1x1_%.2f%s\tG2x2_%.2f%s\tG0x1_%.2f%s\tG1x2_%.2f%s\tG2x0_%.2f%s\tGSD0x0_%.2f%s\tGSD1x1_%.2f%s\tGSD2x2_%.2f%s\tGSD0x1_%.2f%s\tGSD1x2_%.2f%s\tGSD2x0_%.2f%s", time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1]);
	}
    fprintf(fpG,"\n");
    
    int tau = 0;
	for(tau = 0; tau < veclen; tau++)	{	//  Write G arrays to file
    	
		fprintf(fpG, "%2.2e\t", tauarray[tau]);
		for(j=0; j<num_slices-1; j++)	{
			for(k=0; k<6; k++)	{
				fprintf(fpG, "%.6f\t", G_array[6*veclen*j+k*veclen+tau]);
			}
			for(k=0; k<6; k++)	{
				fprintf(fpG, "%.6f\t", G_stdev_array[6*veclen*j+k*veclen+tau]);
			}
		}
		
		for(j=num_slices-1; j<num_slices; j++)	{
			for(k=0; k<6; k++)	{
				fprintf(fpG, "%.6f\t", G_array[6*veclen*j+k*veclen+tau]);
			}
			for(k=0; k<5; k++)	{
				fprintf(fpG, "%.6f\t", G_stdev_array[6*veclen*j+k*veclen+tau]);
			}
			fprintf(fpG, "%.6f\n", G_stdev_array[6*veclen*j+5*veclen+tau]);
		}
		
	}
    fclose(fpG);
    
    fprintf(fpGGG,"tau_1%s\ttau_2%s", argv[1], argv[1]);
    for(slices = 0; slices < num_slices; slices++)	{
		fprintf(fpGGG,"\tG0x1x2_%.2f%s\tG1x2x0_%.2f%s\tG2x0x1_%.2f%s\tG_stdev_0x1x2_%.2f%s\tG_stdev_1x2x0_%.2f%s\tG_stdev_2x0x1_%.2f%s\tGallthree_%.2f%s", time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1]);
	}
    fprintf(fpGGG,"\n");
    
    int tau2 = 0;
	for(tau = 1; tau < veclen; tau++)	{	//  Write Triple Corr arrays to file.  tau1  tau2    G_0x1x2_0(tau1,tau2)  G_1x2x0_0(tau1,tau2)    G_2x0x1_0(tau1,tau2)  G_0x1x2_1(tau1,tau2)  G_1x2x0_1(tau1,tau2)    G_2x0x1_1(tau1,tau2)  G_stdev_0x1x2_0(tau1,tau2)  G_stdev_1x2x0_0(tau1,tau2)    G_stdev_2x0x1_0(tau1,tau2)  G_stdev_0x1x2_1(tau1,tau2)  G_stdev_1x2x0_1(tau1,tau2)    G_stdev_2x0x1_1(tau1,tau2)...
	   for(tau2=1; tau2 < veclen; tau2++)  {
    	
            fprintf(fpGGG, "%2.2e\t%2.2e\t", tauarray[tau2], tauarray[tau]);
            for(j=0; j<num_slices-1; j++)	{
                for(k=0; k<3; k++)	{
                    fprintf(fpGGG, "%.6f\t", G3_array[3*j*veclen*veclen+k*veclen*veclen+tau*veclen+tau2]);
                }
                for(k=0; k<3; k++)	{
                    fprintf(fpGGG, "%.6f\t", G3_stdev_array[3*j*veclen*veclen+k*veclen*veclen+tau*veclen+tau2]);
                }
				fprintf(fpGGG, "%.6f\t", (G3_array[3*j*veclen*veclen+0*veclen*veclen+tau*veclen+tau2]+G3_array[3*j*veclen*veclen+1*veclen*veclen+tau*veclen+tau2]+G3_array[3*j*veclen*veclen+2*veclen*veclen+tau*veclen+tau2])/3.0);
            }
            
            for(j=num_slices-1; j<num_slices; j++)	{
                for(k=0; k<3; k++)	{
                    fprintf(fpGGG, "%.6f\t", G3_array[3*j*veclen*veclen+k*veclen*veclen+tau*veclen+tau2]);
                }
                for(k=0; k<2; k++)	{
                    fprintf(fpGGG, "%.6f\t", G3_stdev_array[3*j*veclen*veclen+k*veclen*veclen+tau*veclen+tau2]);
                }
                fprintf(fpGGG, "%.6f\t", G3_stdev_array[3*j*veclen*veclen+2*veclen*veclen+tau*veclen+tau2]);
				fprintf(fpGGG, "%.6f\n", (G3_array[3*j*veclen*veclen+0*veclen*veclen+tau*veclen+tau2]+G3_array[3*j*veclen*veclen+1*veclen*veclen+tau*veclen+tau2]+G3_array[3*j*veclen*veclen+2*veclen*veclen+tau*veclen+tau2])/3.0);
            }
        }    
        fprintf(fpGGG, "\n");        
	}
    fclose(fpGGG);

//	***	Write selected decays to a seperate file for plotting  *** //
//	Array index order:	Slices, chirality, tau_2 slice, tau_1 vs tau_2 const, [tau_1,tau_2,G(tau_1,tau_2)] 
    int i2;
    int	num_selected_curves = 4;

	for(tau2 = 0; tau2 < veclen; tau2++)	{	//	tau_2 down vertical axis
		fprintf(fpsGGG, "%i", tau2);
		for(j=0; j<num_slices; j++)	{		// 	blocks of curves according to slice number
			for(k=0; k<3; k++)	{				// 	one of each chirality (0x1x2, 1x2x0...)
				for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
					tau = veclen*i2/(num_selected_curves);
					fprintf(fpsGGG, "\t%2.2e\t%2.2e\t%.6f", tauarray[tau2], tauarray[tau], G3_array[3*j*veclen*veclen+k*veclen*veclen+tau*veclen+tau2]);
					fprintf(fpsGGG, "\t%2.2e\t%2.2e\t%.6f", tauarray[tau], tauarray[tau2], G3_array[3*j*veclen*veclen+k*veclen*veclen+tau2*veclen+tau]);
                }
			}
			for(k=0; k<3; k++)	{				// 	one diagonal of each chirality (0x1x2, 1x2x0...)
				fprintf(fpsGGG, "\t%2.2e\t%2.2e\t%.6f", tauarray[tau2], tauarray[tau2], G3_array[3*j*veclen*veclen+k*veclen*veclen+tau2*veclen+tau2]);
			}
		}
		fprintf(fpsGGG, "\n");
	}
            
	fclose(fpsGGG);

	// ***	Script to generate gnuplot commands to display triple correlations *** //
	char scriptfile01[80];	//Open file to write out data
	snprintf(scriptfile01, 50, "samples_GGGs%s.script",file_nametag);
	FILE * fpScript;
    fpScript = fopen(scriptfile01, "w");
	fprintf(fpScript, "set logscale x\nset logscale y\nset ticslevel 0.8\nset view 70, 135\nunset colorbox\nset ylabel \"tau1 (s)\"\nset xlabel \"tau2 (s)\"\nunset key\n");

	fprintf(fpScript, "splot [][][]");
	fprintf(fpScript, " \"%s\" u 1:2:3 with pm3d at b", resultsfile100);
	fprintf(fpScript, ", \"%s\" u 1:2:3 with dots", resultsfile100);
	for(k=0; k<3; k++)	{				// 	one of each chirality (0x1x2, 1x2x0...)
		for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
			fprintf(fpScript, ", \"%s\" u %i:%i:%i with lines", resultsfile110, 2+3*i2+2*2*3*num_selected_curves*k, 3+2*3*i2*2+3*num_selected_curves*k,4+3*i2+3*num_selected_curves*k);
		}
	}
	fprintf(fpScript, "\n");
    
    fclose(fpScript);
    //	***
    
    char resultsfile01[80];	//Open file to write out data
	snprintf(resultsfile01, 50, "bin_clean%s.bin2",file_nametag);
	FILE * fpB;
    fpB = fopen(resultsfile01, "wb");
    //  0   int = veclen
    int veclen_helper = veclen;
    fwrite(&veclen_helper, sizeof(int), 1, fpB);
    veclen_helper = 1;
    fwrite(&veclen_helper, sizeof(int), 1, fpB);
    veclen_helper = veclen;
    //  2 Number of slices (timepoints)
    fwrite(&slices, sizeof(int), 1, fpB);
    //  tau values
    fwrite(tauarray, sizeof(double), veclen_helper, fpB);
    //  For Each timepoint, write Time, Int0,Int1,Int2,SDInt0,SDInt1,SDInt2, 
    //      G0x0,G1x1,G2x2,G0x1,G1x2,G2x0,SDG0x0,SDG1x1,SDG2x2,SDG0x1,SDG1x2,SDG2x0
    for(j=0; j<num_slices; j++)	{
		fwrite(&time_array[j], sizeof(double), 1, fpB);
        fwrite(&intensity0_array[j], sizeof(double), 1, fpB);
        fwrite(&intensity1_array[j], sizeof(double), 1, fpB);
        fwrite(&intensity2_array[j], sizeof(double), 1, fpB);
        fwrite(&intensity0_stdev_array[j], sizeof(double), 1, fpB);
        fwrite(&intensity1_stdev_array[j], sizeof(double), 1, fpB);
        fwrite(&intensity2_stdev_array[j], sizeof(double), 1, fpB);
        for(k=0; k<6; k++)	fwrite(&G_array[6*veclen*j+k*veclen], sizeof(double), veclen_helper, fpB);
        for(k=0; k<6; k++)	fwrite(&G_stdev_array[6*veclen*j+k*veclen], sizeof(double),veclen_helper, fpB);
		for(k=0; k<3; k++)	fwrite(&G3_array[3*veclen*veclen*j+k*veclen*veclen], sizeof(double), veclen*veclen, fpB);
        for(k=0; k<3; k++)	fwrite(&G3_stdev_array[3*veclen*veclen*j+k*veclen*veclen], sizeof(double),veclen*veclen, fpB);
    }
    fclose(fpB);
    printf("Binary traces written to file = %s\n", resultsfile01);
    	
	//for(i=0; i<8; i++)	printf("*^*^*^ %9.9f\n", F1_array[i]);
    
    
cleanupandexit:
	free(time_array);
	free(triple_dump);
	free(triple_dump2);
	free(intensity0_array);
	free(intensity1_array);
	free(intensity2_array);
	free(intensity0_stdev_array);
	free(intensity1_stdev_array);
	free(intensity2_stdev_array);
	free(G_array);
	free(G_stdev_array);
	free(G3_array);
	if(G3_stdev_array!=NULL)        free(G3_stdev_array);
	if(d!=NULL)                     free(d);
	if(tt!=NULL)                    free(tt);
	if(G0x0!=NULL)                  free(G0x0);
	if(G1x1!=NULL)                  free(G1x1);
	if(G2x2!=NULL)                  free(G2x2);
	if(G0x1!=NULL)                  free(G0x1);
	if(G1x2!=NULL)                  free(G1x2);
	if(G2x0!=NULL)                  free(G2x0);
	if(G0x1x2!=NULL)                free(G0x1x2);
	if(G1x2x0!=NULL)                free(G1x2x0);
	if(G2x0x1!=NULL)                free(G2x0x1);
	if(G0x0_prebin!=NULL)           free(G0x0_prebin);
	if(G1x1_prebin!=NULL)           free(G1x1_prebin);
	if(G2x2_prebin!=NULL)           free(G2x2_prebin);
	if(G0x1_prebin!=NULL)           free(G0x1_prebin);
	if(G1x2_prebin!=NULL)           free(G1x2_prebin);
	if(G2x0_prebin!=NULL)           free(G2x0_prebin);
	if(G0x1x2_prebin!=NULL)         free(G0x1x2_prebin);
	if(G1x2x0_prebin!=NULL)         free(G1x2x0_prebin);
	if(G2x0x1_prebin!=NULL)         free(G2x0x1_prebin);
	if(garray0b0!=NULL)             free(garray0b0);
	if(garray0b1!=NULL)             free(garray0b1);
	if(garray1b1!=NULL)             free(garray1b1);
	if(garray1b2!=NULL)             free(garray1b2);
	if(garray2b0!=NULL)             free(garray2b0);
	if(garray2b2!=NULL)             free(garray2b2);
	if(garray0b1b2!=NULL)           free(garray0b1b2);
	if(garray1b2b0!=NULL)           free(garray1b2b0);
	if(garray2b0b1!=NULL)           free(garray2b0b1);
	if(curve_time!=NULL)            free(curve_time);
	if(curve_int0!=NULL)            free(curve_int0);
	if(curve_int1!=NULL)            free(curve_int1);
	if(curve_int2!=NULL)            free(curve_int2);
	free(tempa);
	free(av_array);
	if(GInt_array!=NULL)            free(GInt_array);
	if(GInt_av!=NULL)               free(GInt_av);
	if(GInt_stdev!=NULL)            free(GInt_stdev);
	if(scans_per_file!=NULL)        free(scans_per_file);
    
    printf("Program completed.\n");
    return 0;
}

