/*
 *  Outlier2.c
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
 *  Outlier2.c
 *  Self-contained code to generate F3CS_Outlier2.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "libF3CS.h"

#define	initial_read_length 	25000
#define	mask_vector_l 23423 

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
	int i,j,k,num_curves;
		
	if(argc != 5)  {
        printf("Usage: F3CS_Outlier2 string first_file sigmas binwidth\ne.g. >F3CS_Outlier2 _test1 1 3.2 8\n");
        return 1;
    }
    printf("F3CS_Outlier2 V.3.1\n");
    
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
	snprintf(Garrayfile, 50, "G_array%s_%i.bin",file_nametag, first_file);

    int file_index = first_file;
    //Check existance of files...
	do	{		
		snprintf(Garrayfile, 50, "G_array%s_%i.bin",file_nametag, file_index);
		fp = fopen(Garrayfile, "rb");
		if(fp != NULL) {
			printf("Opening file (%s)\n", Garrayfile);
			fclose(fp);
			file_index += file_delta;
		}
	}
	while (fp != NULL);
    if (file_index == first_file){
    	printf("File (%s) not found, program will quit\n", Garrayfile);
    	return 1;
    }
    int high_file = file_index - file_delta;

    int		valid_passes = 120;
    float   *garray0b0;
    float   *garray0b1;
    float   *garray1b1;
    float   *garray1b2;
    float   *garray2b0;
    float   *garray2b2;
    int garraysize = valid_passes*veclen;
    
    garray0b0 = (float *)malloc(garraysize*sizeof(float));
    garray0b1 = (float *)malloc(garraysize*sizeof(float));
    garray1b1 = (float *)malloc(garraysize*sizeof(float));
    garray1b2 = (float *)malloc(garraysize*sizeof(float));
    garray2b0 = (float *)malloc(garraysize*sizeof(float));
    garray2b2 = (float *)malloc(garraysize*sizeof(float));
    
    for(i=0; i<garraysize; i++) garray0b0[i] = 1.0*i;
    for(i=0; i<garraysize; i++) garray1b1[i] = 2.0*i;
    for(i=0; i<garraysize; i++) garray2b2[i] = 3.0*i;
    for(i=0; i<garraysize; i++) garray0b1[i] = 4.0*i;
    for(i=0; i<garraysize; i++) garray1b2[i] = 5.0*i;
    for(i=0; i<garraysize; i++) garray2b0[i] = 6.0*i;
    
    double *d;
	const int	d_length = 7+6*veclen;
	d = malloc(sizeof(double)*d_length);
	
	//Open file, populate array	
	double	tauarray[veclen];
	double 	dump[veclen*6];
	double	curveheader[7];
	double 	* G0x0;
	double 	* G1x1;
	double 	* G2x2;
	double 	* G0x1;
	double 	* G1x2;
	double 	* G2x0;
	double  * tt;
	double	* curve_time;
	double	* curve_int0;
	double	* curve_int1;
	double	* curve_int2;
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
	
	tt   = malloc(sizeof(double)*6*initial_read_length);
	curve_time = (double *)malloc(sizeof(double)*initial_read_length);
	for(i=0; i<initial_read_length; i++)	curve_time[i] = 0;
	curve_int0 = (double *)malloc(sizeof(double)*initial_read_length);
	for(i=0; i<initial_read_length; i++)	curve_int0[i] = 0;
	curve_int1 = (double *)malloc(sizeof(double)*initial_read_length);
	for(i=0; i<initial_read_length; i++)	curve_int1[i] = 0;
	curve_int2 = (double *)malloc(sizeof(double)*initial_read_length);
	for(i=0; i<initial_read_length; i++)	curve_int2[i] = 0;

	double	time_offset;

	for(i=0; i<veclen; i++) tauarray[i] = 0.0f;
	for(i=0; i<veclen*6; i++)	d[i] = 0.0f;

    file_index = first_file;
    snprintf(Garrayfile, 50, "G_array%s_%i.bin",file_nametag, file_index);
	fp = fopen(Garrayfile, "rb");
	if (fp==NULL)	{printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1);}

	j = fread(&time_offset, sizeof(double), 1, fp);
	j = fread(tauarray, sizeof(double), veclen, fp);
	
	file_index = first_file;
	num_curves=0;
	do {
		//	Time	Int_0	Int_1	Int_2	IntSd_0	IntSD_1	IntSD_2				(7)
		j = fread(curveheader, sizeof(double), 7, fp);
		if (j!=7 && j != 0)	printf("****Error (1), file is the wrong size****\n");
		curve_time[num_curves] = curveheader[0];
		curve_int0[num_curves] = curveheader[1];
		curve_int1[num_curves] = curveheader[2];
		curve_int2[num_curves] = curveheader[3];
		j = fread(dump, sizeof(double), 6*veclen, fp);
		if (j!= 6*veclen && file_index < high_file)	{
			if(j != 0)	printf("Left-over bytes = %i\n", j);
			fclose(fp);
			file_index += file_delta;
			snprintf(Garrayfile, 50, "G_array%s_%i.bin",file_nametag, file_index);
			fp = fopen(Garrayfile, "rb");
			if (fp==NULL)	{printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1);}
			j = fread(&time_offset, sizeof(double), 1, fp);
			fseek(fp, (veclen)*sizeof(double),SEEK_CUR);
			j = fread(curveheader, sizeof(double), 7, fp);
			if (j!=7 && j != 0)	printf("****Error (2), file is the wrong size****\n");
			curve_time[num_curves] = curveheader[0];
			curve_int0[num_curves] = curveheader[1];
			curve_int1[num_curves] = curveheader[2];
			curve_int2[num_curves] = curveheader[3];
            j = fread(dump, sizeof(double), 6*veclen, fp);
		}
		if(j== (6*veclen))	{
			for(i=0; i< veclen; i++)	{
				G0x0[i+num_curves*veclen] = dump[i+0*veclen];
				G1x1[i+num_curves*veclen] = dump[i+1*veclen];
				G2x2[i+num_curves*veclen] = dump[i+2*veclen];
				G0x1[i+num_curves*veclen] = dump[i+3*veclen];
				G1x2[i+num_curves*veclen] = dump[i+4*veclen];
				G2x0[i+num_curves*veclen] = dump[i+5*veclen];
			}
			num_curves++;
		}
	}	while (j==(6*veclen) && num_curves<initial_read_length);
	
	if(num_curves==initial_read_length)	printf("*******\nDataset has more than initial_read_length curves.\nRecompile this program to deal with the data\n********\n");
	if(num_curves > 1)  printf("Total number of curves in the data = %i, from %6.2f s data\n", num_curves, 2.0*curve_time[num_curves-1]-curve_time[num_curves-2]);
    else printf("Total number of curves in the data = %i\n", num_curves);
	
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
	double *	time_array;
	double * 	intensity0_array;
	double * 	intensity1_array;
	double * 	intensity2_array;
	double * 	intensity0_stdev_array;
	double * 	intensity1_stdev_array;
	double * 	intensity2_stdev_array;
	double *	G_array;
	double * 	G_stdev_array;
	double *	GInt_array;
	double *	GInt_av;
	double *	GInt_stdev;

	time_array = (double *) malloc(num_slices*sizeof(double));
	intensity0_array = (double *) malloc(num_slices*sizeof(double));
	intensity1_array = (double *) malloc(num_slices*sizeof(double));
	intensity2_array = (double *) malloc(num_slices*sizeof(double));
	intensity0_stdev_array = (double *) malloc(num_slices*sizeof(double));
	intensity1_stdev_array = (double *) malloc(num_slices*sizeof(double));
	intensity2_stdev_array = (double *) malloc(num_slices*sizeof(double));
	G_array = (double *) malloc(num_slices*veclen*6*sizeof(double));
	if(G_array == NULL)	printf("Malloc Failed!  Likely not enough RAM\n");
	G_stdev_array = (double *) malloc(num_slices*veclen*6*sizeof(double));
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
	for(i=0; i< bw*6; i++)   	GInt_array[i] = 0.0;
	for(i=0; i< 6; i++)   		GInt_av[i] = 0.0;
	for(i=0; i< 6; i++)   		GInt_stdev[i] = 0.0;

	int	slices = 0;
	int	worst_curve = 0;
	double worst_value = 0.0;
	double worst_temp = 0.0;
	double remaining_bw_curves = 0.0;
	double * tempa;
	const double rootbw = 1.0 / pow((double) bw, 0.5);
	tempa = (double *) malloc(6*mask_vector_l*sizeof(double));
	double * av_array;
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
			
			//Average each GaXb
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
                printf("Outlier at time %5.2f s\t(curve %i/%i), \tError = %6.3f x sigmas\n", curve_time[slices*bw+worst_curve], slices*bw+worst_curve,num_curves,worst_value);
                bw_mask[worst_curve] = 0;
            }
			
		}	while (worst_value > sigmas);
		
		
		//Do Math on the Remaining Curves
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
	}
	    
	char resultsfile90[80];	//Open file to write out data
	snprintf(resultsfile90, 50, "clean_Gs%s.dat",file_nametag);
	FILE * fpG;
    fpG = fopen(resultsfile90, "w");
    
    fprintf(fpG,"tau%s", argv[1]);
    for(slices = 0; slices < num_slices; slices++)	{
		fprintf(fpG,"\tG0x0_%.2f%s\tG1x1_%.2f%s\tG2x2_%.2f%s\tG0x1_%.2f%s\tG1x2_%.2f%s\tG2x0_%.2f%s\tGSD0x0_%.2f%s\tGSD1x1_%.2f%s\tGSD2x2_%.2f%s\tGSD0x1_%.2f%s\tGSD1x2_%.2f%s\tGSD2x0_%.2f%s", time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1]);
	}
    fprintf(fpG,"\n");
    
    int tau = 0;
	for(tau = 0; tau < veclen; tau++)	{	//Write G arrays to file
    	
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
    
    
    char resultsfile01[80];	//Open file to write out data
	snprintf(resultsfile01, 50, "bin_clean%s.bin",file_nametag);
	FILE * fpB;
    fpB = fopen(resultsfile01, "wb");
    //  0   int = veclen
    int veclen_helper = veclen;
    fwrite(&veclen_helper, sizeof(int), 1, fpB);
    //  1   int = 0 (Vestigial)
    veclen_helper = 0;
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
    }
    fclose(fpB);
    printf("Binary traces written to file: %s\n", resultsfile01);
    
    if (fp!=NULL) fclose(fp);
    
	/* Free Memory and Exit */
	free(time_array);
	free(intensity0_array);
	free(intensity1_array);
	free(intensity2_array);
	free(intensity0_stdev_array);
	free(intensity1_stdev_array);
	free(intensity2_stdev_array);
	free(G_array);
	free(G_stdev_array);
	free(d);
	free(tt);
	free(G0x0);
	free(G1x1);
	free(G2x2);
	free(G0x1);
	free(G1x2);
	free(G2x0);
	free(garray0b0);
	free(garray0b1);
	free(garray1b1);
	free(garray1b2);
	free(garray2b0);
	free(garray2b2);
	free(curve_time);
	free(curve_int0);
	free(curve_int1);
	free(curve_int2);
	free(tempa);
	free(av_array);
	free(GInt_array);
	free(GInt_av);
	free(GInt_stdev);
    
    printf("Program completed.\n");
    return 0;
}

