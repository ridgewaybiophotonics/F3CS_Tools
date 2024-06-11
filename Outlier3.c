/*
 *  Outlier3.c
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
 *  Outlier3.c
 *  Self-contained code to generate F3CS_Outlier3.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "libF3CS.h"

#define	initial_read_length 	10000
#define	mask_vector_l 4096

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
	int i,j,jj,k,num_curves;
		
	if(argc != 5)  {
        printf("Usage: F3CS_Outlier3 string first_file sigmas binwidth\ne.g. >F3CS_Outlier3 _test1 1 3.2 8\n");
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

    int file_index = first_file;
    //Check existance of files...
	do	{		
		snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, file_index);
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
    float   *garray0b1b2;
    float   *garray1b2b0;
    float   *garray2b0b1;
    int garraysize = valid_passes*veclen;
    
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
	double	tauarray[veclen];
	double 	dump[veclen*6];
	double  * triple_dump;
	double	curveheader[7];
	double 	* G0x0;
	double 	* G1x1;
	double 	* G2x2;
	double 	* G0x1;
	double 	* G1x2;
	double 	* G2x0;
	double 	* G0x0_prebin;
	double 	* G1x1_prebin;
	double 	* G2x2_prebin;
	double 	* G0x1_prebin;
	double 	* G1x2_prebin;
	double 	* G2x0_prebin;
	double 	* G0x1x2_prebin;
	double 	* G1x2x0_prebin;
	double 	* G2x0x1_prebin;
	double  * tt;
	double	* curve_time;
	double	* curve_int0;
	double	* curve_int1;
	double	* curve_int2;
	triple_dump = malloc(sizeof(double)*3*veclen*veclen);
	for(i=0; i< 3*veclen*veclen; i++)	     triple_dump[i] = 0.0;
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
	
	G0x1x2_prebin = malloc(sizeof(double)*veclen*veclen);
	for(i=0; i< veclen*veclen; i++)	G0x1x2_prebin[i] = 0.0;
	G1x2x0_prebin = malloc(sizeof(double)*veclen*veclen);
	for(i=0; i< veclen*veclen; i++)	G1x2x0_prebin[i] = 0.0;
	G2x0x1_prebin = malloc(sizeof(double)*veclen*veclen);
	for(i=0; i< veclen*veclen; i++)	G2x0x1_prebin[i] = 0.0;
	
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
    snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, file_index);
	fp = fopen(Garrayfile, "rb");
	if (fp==NULL)	{printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1);}

	j = fread(&time_offset, sizeof(double), 1, fp);
	j = fread(tauarray, sizeof(double), veclen, fp);
	
	int prebin_counter=0;
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
		j =  fread(dump, sizeof(double), 6*veclen, fp);
		j += fread(triple_dump, sizeof(double), 3*veclen*veclen, fp);
		if (j!= 6*veclen+3*veclen*veclen && file_index < high_file)	{
			if(j != 0)	printf("Left-over bytes = %i\n", j);
			fclose(fp);
			file_index += file_delta;
			snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, file_index);
			fp = fopen(Garrayfile, "rb");
			if (fp==NULL)	{printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1);}
			j = fread(&time_offset, sizeof(double), 1, fp);
			printf("Time Offset = %f / (%i)\n", time_offset, j);
			fseek(fp, (veclen)*sizeof(double),SEEK_CUR);
			j = fread(curveheader, sizeof(double), 7, fp);
			if (j!=7 && j != 0)	printf("****Error (2), file is the wrong size****\n");
			curve_time[num_curves] = curveheader[0];
			curve_int0[num_curves] = curveheader[1];
			curve_int1[num_curves] = curveheader[2];
			curve_int2[num_curves] = curveheader[3];
			j = fread(dump, sizeof(double), 6*veclen, fp);
			if(fseek(fp, (3*veclen*veclen)*sizeof(double), SEEK_CUR) == 0) j+= 3*veclen*veclen;
		}
		if(j== (6*veclen+3*veclen*veclen))	{	// Write data to 8X prebinning arrays.  Average them.
			for(i=0; i< veclen; i++)	{
				G0x0_prebin[i] += dump[i+0*veclen];
				G1x1_prebin[i] += dump[i+1*veclen];
				G2x2_prebin[i] += dump[i+2*veclen];
				G0x1_prebin[i] += dump[i+3*veclen];
				G1x2_prebin[i] += dump[i+4*veclen];
				G2x0_prebin[i] += dump[i+5*veclen];
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
			
			for(i=0; i< veclen; i++)	{
				G0x0_prebin[i] = 0.0;
				G1x1_prebin[i] = 0.0;
				G2x2_prebin[i] = 0.0;
				G0x1_prebin[i] = 0.0;
				G1x2_prebin[i] = 0.0;
				G2x0_prebin[i] = 0.0;
			}

			prebin_counter=0;
			num_curves++;
		}
	}	while (j==(6*veclen+3*veclen*veclen) && num_curves<initial_read_length);
	
	if(fp!= NULL)  fclose(fp);
	
	if(num_curves==initial_read_length)	printf("*******\nDataset has more than initial_read_length=%i curves.\nRecompile this program to deal with the data\n********\n",initial_read_length);
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
	double *	time_array;
	double * 	intensity0_array;
	double * 	intensity1_array;
	double * 	intensity2_array;
	double * 	intensity0_stdev_array;
	double * 	intensity1_stdev_array;
	double * 	intensity2_stdev_array;
	double *	G_array;
	double * 	G_stdev_array;
	double *    G3_array;
	double *    G3_stdev_array;
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
	if(G_array == NULL)	printf("Malloc Failed!  Likely need more RAM\n");
	G_stdev_array = (double *) malloc(num_slices*veclen*6*sizeof(double));
	G3_array = (double *) malloc(num_slices*veclen*veclen*3*sizeof(double));
	if(G3_array == NULL)	printf("Malloc Failed!  Likely need more RAM\n");
	G3_stdev_array = (double *) malloc(num_slices*veclen*veclen*3*sizeof(double));
	if(G3_stdev_array == NULL)	printf("Malloc Failed!  Likely need more RAM\n");
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
	double * tempa;
	const double rootbw = 1.0 / pow((double) bw, 0.5);
	tempa = (double *) malloc(6*mask_vector_l*sizeof(double));
	double * av_array;
	av_array = (double *) malloc((6+1)*bw*sizeof(double));
	for(i=0; i< (6+1)*bw; i++)   av_array[i] = 0.0;
	
	// Start file I/O for reading in triple-correlation data
	file_index = first_file;
    snprintf(Garrayfile, 50, "G_array%s_%i.bin2",file_nametag, file_index);
	fp = fopen(Garrayfile, "rb");
	if (fp==NULL)	{printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1);}
	j = fread(&time_offset, sizeof(double), 1, fp);
	j = fread(tauarray, sizeof(double), veclen, fp);
	
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
        
		for(i=0; i<bw; i++) remaining_bw_curves += bw_mask[i]*1.0;
        		
		//Find Outliers
		do{
			remaining_bw_curves = 0.0;
			for(i=0; i<bw; i++) remaining_bw_curves += bw_mask[i]*1.0;
			for(i=0; i<6*bw; i++)	GInt_array[i] = 0.0;
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
                printf("Outlier found with error %f x sigma\n", worst_value);
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
		
		
        // To keep memory usage down, the triple-correlation curves are read in and binned on-the-fly
		for(k=0; k<veclen*veclen; k++)	{
			G3_array[3*slices*veclen*veclen + 0*veclen*veclen + k] = 0.0;
			G3_array[3*slices*veclen*veclen + 1*veclen*veclen + k] = 0.0;
			G3_array[3*slices*veclen*veclen + 2*veclen*veclen + k] = 0.0;
			G3_stdev_array[3*slices*veclen*veclen + 0*veclen*veclen + k] = 0.0;
			G3_stdev_array[3*slices*veclen*veclen + 1*veclen*veclen + k] = 0.0;
			G3_stdev_array[3*slices*veclen*veclen + 2*veclen*veclen + k] = 0.0;
        }

        for(jj=0; jj<bw; jj++) {
            //  Read in data
            for(i=0; i< veclen*veclen; i++)	{        G0x1x2_prebin[i] = 0.0;        G1x2x0_prebin[i] = 0.0;        G2x0x1_prebin[i] = 0.0;}
            for(prebin_counter=0; prebin_counter <8; prebin_counter++) {
                j = fread(curveheader, sizeof(double), 7, fp);
                if (j!=7 && j != 0)	printf("****Error (3), file is the wrong size****\n");
                j=0;
                j =  fread(dump, sizeof(double), 6*veclen, fp);
                j += fread(triple_dump, sizeof(double), 3*veclen*veclen, fp);
                if (j!= 6*veclen+3*veclen*veclen && file_index < high_file)	{
                    if(j != 0)	printf("Left-over bytes = %i\n", j);
                    fclose(fp);
                    file_index += file_delta;
                    snprintf(Garrayfile, 50, "G_array%s_%i.bin2", file_nametag, file_index);
                    fp = fopen(Garrayfile, "rb");
                    if (fp==NULL)	{printf("file %s not found.  Exiting.\n", Garrayfile); exit(-1);}
                    j = fread(&time_offset, sizeof(double), 1, fp);
                    fseek(fp, (veclen)*sizeof(double),SEEK_CUR);
                    j = fread(curveheader, sizeof(double), 7, fp);
                    if (j!=7 && j != 0)	printf("****Error (4), file is the wrong size****\n");
                    j=0;
                    j =  fread(dump, sizeof(double), 6*veclen, fp);
                    j += fread(triple_dump, sizeof(double), 3*veclen*veclen, fp);
                    if (j!= 6*veclen+3*veclen*veclen && file_index < high_file)    printf("XXX Problem reading in triple correlations XXX\n");
                }
                for(i=0; i< veclen*veclen; i++)	{
                    G0x1x2_prebin[i] += triple_dump[i+0*veclen*veclen]/8.0;
                    G1x2x0_prebin[i] += triple_dump[i+1*veclen*veclen]/8.0;
                    G2x0x1_prebin[i] += triple_dump[i+2*veclen*veclen]/8.0;
                }
            }
            
            if(bw_mask[jj]>0)    {    
                for(k=0; k<veclen*veclen; k++)	{
                    G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k] += G0x1x2_prebin[k];
                    G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k] += G1x2x0_prebin[k];
                    G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k] += G2x0x1_prebin[k];
                    G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] += pow( G0x1x2_prebin[k], 2.0);
                    G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] += pow( G1x2x0_prebin[k], 2.0);
                    G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] += pow( G2x0x1_prebin[k], 2.0);
                }   
            }
        }
			
        for(k=0; k<veclen*veclen; k++)	{  //			stdev = sqrt (N/(N-1) * [E(x^2) - (E(x))^2])
        
			G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k] / remaining_bw_curves;
			G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k] / remaining_bw_curves;
			G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k] / remaining_bw_curves;
			
			G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] / (remaining_bw_curves);
			G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] / (remaining_bw_curves);
			G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] / (remaining_bw_curves);

            G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] - pow( G3_array[slices*3*veclen*veclen + 0*veclen*veclen + k],2.0);
			G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] - pow( G3_array[slices*3*veclen*veclen + 1*veclen*veclen + k],2.0);
			G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] - pow( G3_array[slices*3*veclen*veclen + 2*veclen*veclen + k],2.0);

			G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = rootbw * pow(((remaining_bw_curves)/(remaining_bw_curves-1.0))*G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k], 0.5);
			G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = rootbw * pow(((remaining_bw_curves)/(remaining_bw_curves-1.0))*G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k], 0.5);
			G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = rootbw * pow(((remaining_bw_curves)/(remaining_bw_curves-1.0))*G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k], 0.5);

			G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 0*veclen*veclen + k] ;
			G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 1*veclen*veclen + k] ;
			G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] = G3_stdev_array[slices*3*veclen*veclen + 2*veclen*veclen + k] ;
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
	
    if(fp!= NULL)  fclose(fp);

	char resultsfile100[80];	//Open file to write out data
	snprintf(resultsfile100, 50, "clean_GGGs%s.dat",file_nametag);
	FILE * fpGGG;
    fpGGG = fopen(resultsfile100, "w");
    
    char resultsfile110[80];	//Open file to write out data
	snprintf(resultsfile110, 50, "samples_GGGs%s.dat",file_nametag);
	FILE * fpsGGG;
    fpsGGG = fopen(resultsfile110, "w");
    
    fprintf(fpGGG,"tau_1%s\ttau_2%s", argv[1], argv[1]);
    for(slices = 0; slices < num_slices; slices++)	{
		fprintf(fpGGG,"\tG0x1x2_%.2f%s\tG1x2x0_%.2f%s\tG2x0x1_%.2f%s\tG_stdev_0x1x2_%.2f%s\tG_stdev_1x2x0_%.2f%s\tG_stdev_2x0x1_%.2f%s\tGallthree_%.2f%s", time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1], time_array[slices], argv[1]);
	}
    fprintf(fpGGG,"\n");
    
    int tau  = 0;
    int tau2 = 0;
	for(tau = 1; tau < veclen; tau++)	{	//Write Triple Corr arrays to file.  tau1  tau2    G_0x1x2_0(tau1,tau2)  G_1x2x0_0(tau1,tau2)    G_2x0x1_0(tau1,tau2)  G_0x1x2_1(tau1,tau2)  G_1x2x0_1(tau1,tau2)    G_2x0x1_1(tau1,tau2)  G_stdev_0x1x2_0(tau1,tau2)  G_stdev_1x2x0_0(tau1,tau2)    G_stdev_2x0x1_0(tau1,tau2)  G_stdev_0x1x2_1(tau1,tau2)  G_stdev_1x2x0_1(tau1,tau2)    G_stdev_2x0x1_1(tau1,tau2)...
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
//	Array index order:	Slices, corr (0x1x2, 1x2x0...), tau_2 slice, tau_1 vs tau_2 const, [tau_1,tau_2,G(tau_1,tau_2)] 
    int i2;
    int	num_selected_curves = 4;

	for(tau2 = 0; tau2 < veclen; tau2++)	{	//	tau_2 down vertical axis
		fprintf(fpsGGG, "%i", tau2);
		for(j=0; j<num_slices; j++)	{		// 	blocks of curves according to slice number
			for(k=0; k<3; k++)	{				// 	one of each corr (0x1x2, 1x2x0...)
				for(i2=0; i2 < num_selected_curves; i2++)  {	// one of each decay group
					tau = veclen*i2/(num_selected_curves-1);
					if(tau >= veclen) tau = veclen-1;
					fprintf(fpsGGG, "\t%2.2e\t%2.2e\t%.6f", tauarray[tau2], tauarray[tau], G3_array[3*j*veclen*veclen+k*veclen*veclen+tau*veclen+tau2]);
					fprintf(fpsGGG, "\t%2.2e\t%2.2e\t%.6f", tauarray[tau], tauarray[tau2], G3_array[3*j*veclen*veclen+k*veclen*veclen+tau2*veclen+tau]);
                }
			}
			for(k=0; k<3; k++)	{				// 	one diagonal of each corr (0x1x2, 1x2x0...)
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

	//fprintf(fpScript, "splot [4.00e-07:0.102][4.00e-07:0.102][]");
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
    
    char resultsfile01[80];	//Open file to write out data
	snprintf(resultsfile01, 50, "bin_clean%s.bin2",file_nametag);
	FILE * fpB;
    fpB = fopen(resultsfile01, "wb");
    //  0   int = veclen
    int veclen_helper = veclen;
    fwrite(&veclen_helper, sizeof(int), 1, fpB);
    veclen_helper = 0;  // (Vestigial)
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
    	
	/*  Free up memory and exit */
	free(time_array);
	free(triple_dump);
	free(intensity0_array);
	free(intensity1_array);
	free(intensity2_array);
	free(intensity0_stdev_array);
	free(intensity1_stdev_array);
	free(intensity2_stdev_array);
	free(G_array);
	free(G_stdev_array);
	free(G3_array);
	free(G3_stdev_array);
	free(d);
	free(tt);
	free(G0x0);
	free(G1x1);
	free(G2x2);
	free(G0x1);
	free(G1x2);
	free(G2x0);
	free(G0x0_prebin);
	free(G1x1_prebin);
	free(G2x2_prebin);
	free(G0x1_prebin);
	free(G1x2_prebin);
	free(G2x0_prebin);
	free(G0x1x2_prebin);
	free(G1x2x0_prebin);
	free(G2x0x1_prebin);
	free(garray0b0);
	free(garray0b1);
	free(garray1b1);
	free(garray1b2);
	free(garray2b0);
	free(garray2b2);
	free(garray0b1b2);
	free(garray1b2b0);
	free(garray2b0b1);
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

