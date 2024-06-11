/*
 *  2FCS.c
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
 *  2FCS.c:
 *  This is the main function for the F3CS_2FCS program, and handles the 
 *  busy work by performing file i/o, launching multithreading and calling 
 *  the correlation integrals coded in corr_int_2FCS.c.
 *
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "libF3CS.h"

int main(int argc, char *argv[])    {	// The Main function divides the datasets up into packets and sends each packet off to be processed by an individual thread
	if(argc != 5)  {
		printf("Usage: F3CS_2FCS filenametag #first_file #last_file #cores_to_use	\ne.g. >F3CS_2FCS _s7 1 44 6\n");
	return 1;
    }
    int i,j,t;
	int	debug_flag = 1;
    
	FILE * fp;
	
	//	Parse command-line inputs
	char	temp_char[512];
	sprintf(temp_char, "%s", argv[1]);
	char *file_nametag;
	file_nametag = temp_char;
    sscanf(argv[2],"%i", &i);
   	int	first_file = i;
   	sscanf(argv[3],"%i", &i);
   	int last_file = i;
	int curvestobin_rxn = 1;	// Vestigial option, leave set to 1
    sscanf(argv[4],"%i", &i);
    int	cores_to_use = i;
    if (cores_to_use < 1)	{
		cores_to_use = 1;
    	printf("# Cores requested < 1; %i core will be used.\n",cores_to_use);
    }
	else	if(cores_to_use>NUM_THREADS)	{
		printf("Too many cores requested (%i).  Request %i or less, or recompile with larger NUM_THREADS definition in main_sym.c\n", cores_to_use, NUM_THREADS); 
		exit(-1);
	}
	else	printf("%i Cores will be used\n", cores_to_use);
    
    //	Obtain time per file to estimate memory usage
	char test_file[80];
	char numb[5];
	int	fp_counter = 0;
	snprintf(numb, 5, "%i", first_file);
	strcpy(test_file, "rxn_");
	strcat(test_file, numb);
	strcat(test_file, file_nametag);
	strcat(test_file, ".w3c");
	fp = fopen(test_file, "rb");
	if(fp == NULL) {
		printf("File (%s) could not be opened.  Exiting.\n", test_file);
		exit(-1);
	}
	unsigned short * tester_array;
	tester_array = (unsigned short*) malloc(1024*3*sizeof(unsigned short));
	for(t=0; t<3*1024; t++)	tester_array[t] = 1;
	if(!tester_array)  printf("tester  Malloc Failed.\n");
	do	{
		fp_counter += fread(tester_array, sizeof(unsigned short), 3*1024, fp);
	}	while	(!feof(fp));
	fclose(fp);
	free(tester_array);
	int	u_shorts_per_file = fp_counter;
	const double	time_per_file = 0.0000008 * (double) (fp_counter) / 3.0;

	//Prepare data structures to send to threads
	struct	stream_struct thread_arg[NUM_THREADS];
	int	num_files_per_thread = (last_file-first_file+1)/ cores_to_use;
	if ((last_file-first_file+1) % cores_to_use > 0)	num_files_per_thread++;
	for(t=0; t<NUM_THREADS; t++)	{
		thread_arg[t].is_first_file = 0;
		thread_arg[t].file_nametag = temp_char;
		thread_arg[t].first_file = first_file;
		thread_arg[t].last_file = last_file;    
		thread_arg[t].curvestobin_rxn = curvestobin_rxn;
		thread_arg[t].time_per_file = time_per_file;
		thread_arg[t].core_id = t;
		thread_arg[t].u_shorts_per_file = u_shorts_per_file;
		thread_arg[t].debug_flag = debug_flag;
	}


	for(j=first_file; j<=last_file; j++)	{		
		snprintf(numb, 5, "%i",j);
		strcpy(test_file, "rxn_");
		strcat(test_file, numb);
		strcat(test_file, file_nametag);
		strcat(test_file, ".w3c");
		fp = fopen(test_file, "rb");
		if(fp == NULL) {
			printf("File (%s) could not be opened.  Exiting.\n", test_file);
			exit(-1);
		}	
		else{
			printf("File (%s) Exists.  Good.\n", test_file);
			fclose(fp);
		} 
	}

	int rc;
	pthread_t threads[NUM_THREADS];
	pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int *status;
    status = 0;
    
    //	Multithreaded dataprocessing works by calling function stream_sipper_thread with the names of the files to process.
    int	file_index = first_file;
    int is_first_file = 1;
    int	files_to_process = 0;
    int t_used = 0;
    while(file_index <= last_file)	{
		t_used = 0;
		for(t=0; t<cores_to_use; t++)	{
			if(file_index <= last_file)		{		
				thread_arg[t].is_first_file = 0;
				if(is_first_file==1)	{
					is_first_file = 0;
					thread_arg[t].is_first_file = 1;
				}
				thread_arg[t].file_nametag = temp_char;
				thread_arg[t].first_file = file_index;
				files_to_process = last_file - file_index + 1;
				if(files_to_process > num_curves_per_thread) files_to_process = num_curves_per_thread;
				file_index += files_to_process;
				thread_arg[t].last_file = file_index - 1;    
				thread_arg[t].curvestobin_rxn = curvestobin_rxn;
				thread_arg[t].time_per_file = time_per_file;
			
				if(debug_flag==1)		printf("going to start thread\n");
				rc = pthread_create(&threads[t], &attr, stream_sipper_thread, (void *) &thread_arg[t]);
				if (rc){
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
				else    {t_used++; if(debug_flag==1) printf("thread %i launched\n", t_used);}
			}
		}
		//	Let Threads Execute, then join and terminate the program
		for(t=0; t<t_used; t++)	pthread_join(threads[t], (void **)status);		
	}

    printf("Program completed.\n");
    pthread_attr_destroy(&attr);
	pthread_exit(NULL);
}

void    *stream_sipper_thread(void *threadarg)	{
	struct	stream_struct * val;
	val = threadarg;
	if(val->debug_flag==1)	printf("debug 1 -- entered threads\n");
    
	//	Allocate work arrays
    unsigned short  *array0_0;
	unsigned short  *array0_1;
	unsigned short	*array0_2;
	unsigned short  *array1_0;
	unsigned short  *array1_1;
	unsigned short	*array1_2;
	unsigned short  *interleaved;
	
	long long int   total_counter;
	struct triplecorr_struct targ1;
    struct triplecorr_struct targ2;
    
    struct fcs_variables pass_to_fcs;
    targ1.fcsv = &pass_to_fcs;
    targ2.fcsv = &pass_to_fcs;
    
	char *file_nametag = val->file_nametag;
   	int	first_file = val->first_file;
   	int last_file = val->last_file;

   	printf("Thread(%i):  First file = %i\tLast file =%i\n", val->core_id, first_file, last_file);
	int i,j,l;
   	
   	float	time_offset = val->time_per_file*(val->first_file-1.0);
	if(val->debug_flag==1)	printf("Timeoffset = %f, time_per_file = %f, first_file -1.0 = %f\n", time_offset, val->time_per_file, (val->first_file-1.0));
   	
    int curvestobin_rxn = val->curvestobin_rxn;    
	FILE * fp;
	FILE * fpW;
    int corr1;
	char resultsfile02[FILE_STRING_LENGTH];
	snprintf(resultsfile02, FILE_STRING_LENGTH, "rxn_%i%s.w3c", first_file, file_nametag);
	printf("Filename:\t%s\n", resultsfile02);
    
    fp = fopen(resultsfile02, "rb");
	if(fp == NULL) {
        printf("File (%s) could not be opened.  Exiting.\n", resultsfile02);
        return 0;
    }
    else    if(val->debug_flag==1)	printf("Reading bursts from file %s...\n", resultsfile02);
    
    float norm_tau[veclen];
    float gfinal0b0[veclen];
    float gfinal0b1[veclen];
    float gfinal1b1[veclen];
    float gfinal1b2[veclen];
    float gfinal2b0[veclen];
    float gfinal2b2[veclen];
    
    for(i=0; i<veclen; i++)    norm_tau[i] = 0;
    algorquantum = n * (int) powl(2,pmax);
    array0_0 = malloc (sizeof(unsigned short)*algorquantum);
    array0_1 = malloc (sizeof(unsigned short)*algorquantum);
    array0_2 = malloc (sizeof(unsigned short)*algorquantum);
    array1_0 = malloc (sizeof(unsigned short)*algorquantum);
    array1_1 = malloc (sizeof(unsigned short)*algorquantum);
    array1_2 = malloc (sizeof(unsigned short)*algorquantum);
    interleaved = malloc (sizeof(unsigned short)*algorquantum*3);
    for(i=0; i< algorquantum; i++)   array0_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_2[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_2[i] = 0;
    
    for(i=0; i<veclen; i++)   gfinal0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)   gfinal0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)   gfinal1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)   gfinal1b2[i] = 0.0f;
	for(i=0; i<veclen; i++)   gfinal2b0[i] = 0.0f;
    for(i=0; i<veclen; i++)   gfinal2b2[i] = 0.0f;
    
    //	***	Malloc pass_to_fcs arrays	***	//
	pass_to_fcs.g0b0 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.g1b1 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.g2b2 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.g0b1 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.g1b2 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.g2b0 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.shatzel0 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.shatzel1 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.shatzel2 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.stemp0 = (float *) malloc(pmax*sizeof(float));
	pass_to_fcs.stemp1 = (float *) malloc(pmax*sizeof(float));
	pass_to_fcs.stemp2 = (float *) malloc(pmax*sizeof(float));
	pass_to_fcs.klaus0 = (float *) malloc(pmax*sizeof(float));
	pass_to_fcs.klaus1 = (float *) malloc(pmax*sizeof(float));
	pass_to_fcs.klaus2 = (float *) malloc(pmax*sizeof(float));
	pass_to_fcs.sc0 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.sc1 = (float *) malloc(veclen*sizeof(float));
	pass_to_fcs.sc2 = (float *) malloc(veclen*sizeof(float));
    
    for(i=0; i<veclen; i++)   	pass_to_fcs.shatzel0[i] = 0.0;
    for(i=0; i<veclen; i++)   	pass_to_fcs.shatzel1[i] = 0.0;
    for(i=0; i<veclen; i++)   	pass_to_fcs.shatzel2[i] = 0.0;
    for(i=0; i<veclen; i++)		pass_to_fcs.sc0[i] = 0.0;
    for(i=0; i<veclen; i++)		pass_to_fcs.sc1[i] = 0.0;
    for(i=0; i<veclen; i++)		pass_to_fcs.sc2[i] = 0.0;
    for(i=0; i<pmax; i++)		pass_to_fcs.stemp0[i] = 0.0;
    for(i=0; i<pmax; i++)		pass_to_fcs.stemp1[i] = 0.0;
    for(i=0; i<pmax; i++)		pass_to_fcs.stemp2[i] = 0.0;
    for(i=0; i<pmax; i++)		pass_to_fcs.klaus0[i] = 0.0;
    for(i=0; i<pmax; i++)		pass_to_fcs.klaus1[i] = 0.0;
    for(i=0; i<pmax; i++)		pass_to_fcs.klaus2[i] = 0.0;
   
    //Calculate anticipated data size, and number of curves to bin
    printf("curves to avg = %i\n", curvestobin_rxn);
    float time_quantum = 800.0e-9f;

    targ1.array0 = (unsigned short*) array0_0;
    targ1.array1 = (unsigned short*) array0_1;
    targ1.array2 = (unsigned short*) array0_2;
    targ2.array0 = (unsigned short*) array1_0;
    targ2.array1 = (unsigned short*) array1_1;
    targ2.array2 = (unsigned short*) array1_2;

    int     fp_counter  = 0;
    int     fp_counter2 = 0;
    pthread_attr_t attr;
    pthread_t Acq_thread1;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//  To avoid including partially-weighted values of gIbJ into the std. dev. calc.s, 
	//	first do a "dummy scan" which increments sumbox & tempxIbJ only, then clear 
	//	gIbJ(sq)s, then resume calculations.
	if(val->debug_flag==1)	printf("Going into threads\n");

	if(val->is_first_file != 1)	{	//Take dummy-scan data from the previous file if it exists.  Thus all data (save the first set) are analyzed, none are wasted.
		FILE * fp_previous;
		snprintf(resultsfile02, FILE_STRING_LENGTH, "rxn_%i%s.w3c", first_file-1,file_nametag);
		if(val->debug_flag==1)	printf("Previous filename:\t%s\n", resultsfile02);
		
		fp_previous = fopen(resultsfile02, "rb");
		if(fp_previous == NULL) {
			printf("(Previous) File (%s) could not be opened.  Exiting.\n", resultsfile02);
			pthread_exit(NULL);
		}
		else    if(val->debug_flag==1)	printf("Reading bursts from (Previous) file %s...\n", resultsfile02);
    	fseek(fp_previous,(val->u_shorts_per_file - 3*algorquantum)*sizeof(unsigned short),SEEK_SET);
		
		fp_counter = fread(interleaved, sizeof(unsigned short), 3*algorquantum, fp_previous);
		if (fp_counter != 3*algorquantum)	{printf("(Previous) Files must be larger than 2 algorquanta (%s).  Program quitting...\n",resultsfile02); exit(-1);}
		for(i=0; i< fp_counter; i += 3)  {
			array1_0[i/3] = interleaved[i];
			array1_1[i/3] = interleaved[i+1];
			array1_2[i/3] = interleaved[i+2];
		}
		fclose(fp_previous);
	}
	else	{
		fp_counter = fread(interleaved, sizeof(unsigned short), 3*algorquantum, fp);
		if (fp_counter != 3*algorquantum)	{printf("Files must be larger than 2 algorquanta.  Program quitting...\n"); pthread_exit(NULL);}
		for(i=0; i< fp_counter; i += 3)  {
			array1_0[i/3] = interleaved[i];
			array1_1[i/3] = interleaved[i+1];
			array1_2[i/3] = interleaved[i+2];
		}
	}

    corr1 = pthread_create(&Acq_thread1, NULL, TripleCorr, (void *) &targ2);
    pthread_join(Acq_thread1, NULL);
	
    fp_counter = fread(interleaved, sizeof(unsigned short), 3*algorquantum, fp);    
    if (fp_counter != 3*algorquantum)	{printf("Files must be larger than 2 algorquanta.  Program quitting...\n"); exit(-1);}
    for(i=0; i< fp_counter; i += 3)  {
        array0_0[i/3] = interleaved[i];
        array0_1[i/3] = interleaved[i+1];
        array0_2[i/3] = interleaved[i+2];
    }
 
 	int file_index = first_file;
 
	//	Open G_array_tag_*.bin results file to write data
	//	G-array structure, v.2:	1+veclen+valid_passes*(7+6*veclen)
	//	time_offset	tau														(1+veclen)
	//		For Ea. valid_pass
	//	Time	Int_0	Int_1	Int_2	IntSd_0	IntSD_1	IntSD_2				(7)
	//	G0x0, G1x1, G2x2, G0x1, G1x2, G2x0									(6*veclen)

	snprintf(resultsfile02, FILE_STRING_LENGTH, "G_array%s_%i.bin", file_nametag, first_file);
	fpW = fopen(resultsfile02, "wb");  
	
	time_offset = val->time_per_file*(val->first_file-1.0);
	printf("Timeoffset = %f\n", time_offset);
	
    double * d;
    int	d_length = 7+6*veclen;
	d = (double *) malloc(sizeof(double)*d_length);
	if(!d)	printf("Array D Failed to malloc\n");
    
    d[0] = time_offset*1.0;
	for(i=0; i<2*n; i++)   	d[i+1] = (i)*time_quantum;
	for(j=1; j<pmax; j++)	for(l=0; l<n; l++)   d[l+n*(j+1)+1] = time_quantum*((powf(2,j)*(n+l))+1);
        
	//write time + tau to file  
	fwrite(d, sizeof(double), veclen+1, fpW);
	int	pointer_0 = 0;
	
	total_counter = 0;
    int m = total_counter / (3*algorquantum);
    int q = 0;
	int valid_passes = 0;
    while (!feof (fp))   {	//	Main loop -- run until you run out of data (feof(fp))
		//	Re-zero working memory for sums
        for(i=0; i<veclen; i++)    pass_to_fcs.g0b0[i] = 0.0f;
        for(i=0; i<veclen; i++)    pass_to_fcs.g0b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    pass_to_fcs.g1b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    pass_to_fcs.g1b2[i] = 0.0f;
        for(i=0; i<veclen; i++)    pass_to_fcs.g2b0[i] = 0.0f;
        for(i=0; i<veclen; i++)    pass_to_fcs.g2b2[i] = 0.0f;
        
        pass_to_fcs.tc0 = 0.0f;
        pass_to_fcs.tc0sq = 0.0f;
        pass_to_fcs.tc1 = 0.0f;
        pass_to_fcs.tc1sq = 0.0f;
        pass_to_fcs.tc2 = 0.0f;
        pass_to_fcs.tc2sq = 0.0f;
		for(i=0; i<veclen; i++)		  pass_to_fcs.sc0[i] = 0;
		for(i=0; i<veclen; i++)		  pass_to_fcs.sc1[i] = 0;
		for(i=0; i<veclen; i++)		  pass_to_fcs.sc2[i] = 0;
        
        total_counter = 0;
		q++;
		if (fp_counter == 3*algorquantum)   {   
			total_counter += fp_counter;
			
			//	Launch the final data processing loop (the pseudo-fractal integrals):
			corr1 = pthread_create(&Acq_thread1, NULL, TripleCorr, (void *) &targ1);
			if (corr1)  {
				printf("ERROR 7262736, multi-threading problem; return code from pthread_create() is %d\n", corr1);
				exit(-1);
			}
			pthread_join(Acq_thread1, NULL);
		   
			fp_counter = fread(interleaved, sizeof(unsigned short), 3*algorquantum, fp);	//	Read in the next data, if it is there
			if(fp_counter != 3*algorquantum && file_index < last_file)	{	//	If this file has been used up and more files exist in the 
																			//	assigned range of this thread, use the next file.
				if (fp_counter%3 != 0)	{	//Check for unexpectely small files
					printf("Possible file corruption / error in data processing as files are too small: size != 3*n*algorquantum, n (in) Integers\n");
					exit(-1);
				}
				fclose(fp);
				file_index++;
				snprintf(resultsfile02, FILE_STRING_LENGTH, "rxn_%i%s.w3c", file_index,file_nametag);
				printf("Current file:\t%s\n", resultsfile02);
				fp = fopen(resultsfile02, "rb");
				if(fp == NULL) {
					printf("File (%s) could not be opened (but it existed when the program started).  Exiting.\n", resultsfile02);
					exit(-1);
				}
				fp_counter2 = fread(&interleaved[fp_counter], sizeof(unsigned short), 3*algorquantum - fp_counter, fp);
				if (fp_counter2 == 3*algorquantum - fp_counter)	{
					fp_counter = 3*algorquantum;	//Allow loops to procede
				}
				else	{
					printf("File (%s) is abnormally short.  Exiting.\n", resultsfile02);
					exit(-1);
				}
			}			
			//	Take interleaved data, and write to separate arrays
			for(i=0; i< fp_counter; i += 3)  {
				array0_0[i/3] = interleaved[i];	//	changed from array1_X[i/3]...
				array0_1[i/3] = interleaved[i+1];
				array0_2[i/3] = interleaved[i+2];
			}                
		}
        else    if(val->debug_flag==1)	printf("q = %i / %i (%i)\n", q, curvestobin_rxn, val->core_id);
   
		// *******	Normallize double correlation data by the number of summation events and total intensity ******* //
		
		for(i=0; i<2*n; i++)   {
			norm_tau[i]  = total_counter/3.0f;	//  = T 
			//	In our publications, this factor is 1/T, not T.  The reason is that pass_to_fcs.scX[i]/tcX 
			//	are not yet normallized by a factor of T each, and (1/T)/(1/T^2) = T
			gfinal0b0[i] = pass_to_fcs.g0b0[i] * norm_tau[i] / ((float)pass_to_fcs.sc0[i]*pass_to_fcs.tc0);
			gfinal0b1[i] = pass_to_fcs.g0b1[i] * norm_tau[i] / ((float)pass_to_fcs.sc0[i]*pass_to_fcs.tc1);
			gfinal1b1[i] = pass_to_fcs.g1b1[i] * norm_tau[i] / ((float)pass_to_fcs.sc1[i]*pass_to_fcs.tc1); //Symmetric normallization
			gfinal1b2[i] = pass_to_fcs.g1b2[i] * norm_tau[i] / ((float)pass_to_fcs.sc1[i]*pass_to_fcs.tc2);
			gfinal2b0[i] = pass_to_fcs.g2b0[i] * norm_tau[i] / ((float)pass_to_fcs.sc2[i]*pass_to_fcs.tc0);
			gfinal2b2[i] = pass_to_fcs.g2b2[i] * norm_tau[i] / ((float)pass_to_fcs.sc2[i]*pass_to_fcs.tc2);
			
		}

		for(j=1; j<pmax; j++)  {
			for(l=0; l<n; l++)   {    // i is loop index (incr. below)
				//	Norm = T^2 / [2^j 2^j floor(T / 2^j)]  (written to safeguard against short datasets where 
				//		T is not a multiple of 2^j and processing terminates early)
				//	If T is a multiple of 2^j, this comes out to Norm = T / 2^j.
				norm_tau[i]  = (float)(total_counter*total_counter/9.0f)/(powl(2,j)*(powl(2,j)*(floor(total_counter / (3.0f*powl(2,j))))));
				gfinal0b0[i] = pass_to_fcs.g0b0[i] * norm_tau[i] / ((float)pass_to_fcs.sc0[i]*pass_to_fcs.tc0);
				gfinal0b1[i] = pass_to_fcs.g0b1[i] * norm_tau[i] / ((float)pass_to_fcs.sc0[i]*pass_to_fcs.tc1);
				gfinal1b1[i] = pass_to_fcs.g1b1[i] * norm_tau[i] / ((float)pass_to_fcs.sc1[i]*pass_to_fcs.tc1);
				gfinal1b2[i] = pass_to_fcs.g1b2[i] * norm_tau[i] / ((float)pass_to_fcs.sc1[i]*pass_to_fcs.tc2);
				gfinal2b0[i] = pass_to_fcs.g2b0[i] * norm_tau[i] / ((float)pass_to_fcs.sc2[i]*pass_to_fcs.tc0);
				gfinal2b2[i] = pass_to_fcs.g2b2[i] * norm_tau[i] / ((float)pass_to_fcs.sc2[i]*pass_to_fcs.tc2);
				
				i++;
			}
		}
		
		// *******	Reset symmetric intensity counters ******* //
		for(i=0; i<veclen; i++)		  pass_to_fcs.sc0[i] = 0;
		for(i=0; i<veclen; i++)		  pass_to_fcs.sc1[i] = 0;
		for(i=0; i<veclen; i++)		  pass_to_fcs.sc2[i] = 0;
		
		m = total_counter / (3*algorquantum);
		pass_to_fcs.tc0 = pass_to_fcs.tc0/(time_quantum*(float)m);
		pass_to_fcs.tc1 = pass_to_fcs.tc1/(time_quantum*(float)m);
		pass_to_fcs.tc2 = pass_to_fcs.tc2/(time_quantum*(float)m);
		pass_to_fcs.tc0sq = pass_to_fcs.tc0sq/(time_quantum*(float)m);
		pass_to_fcs.tc1sq = pass_to_fcs.tc1sq/(time_quantum*(float)m);
		pass_to_fcs.tc2sq = pass_to_fcs.tc2sq/(time_quantum*(float)m);
		pass_to_fcs.tc0sq = powf((pass_to_fcs.tc0sq - (pass_to_fcs.tc0*pass_to_fcs.tc0)),0.5f);
		pass_to_fcs.tc1sq = powf((pass_to_fcs.tc1sq - (pass_to_fcs.tc1*pass_to_fcs.tc1)),0.5f);
		pass_to_fcs.tc2sq = powf((pass_to_fcs.tc2sq - (pass_to_fcs.tc2*pass_to_fcs.tc2)),0.5f);
		pass_to_fcs.tc0 = pass_to_fcs.tc0/(float)algorquantum;
		pass_to_fcs.tc1 = pass_to_fcs.tc1/(float)algorquantum;
		pass_to_fcs.tc2 = pass_to_fcs.tc2/(float)algorquantum;
		pass_to_fcs.tc0sq = pass_to_fcs.tc0sq/(float)algorquantum;
		pass_to_fcs.tc1sq = pass_to_fcs.tc1sq/(float)algorquantum;
		pass_to_fcs.tc2sq = pass_to_fcs.tc2sq/(float)algorquantum;
		
		// *******	Transfer correlation data to a large array for storage ******* //
		pointer_0 = 0;
		d[pointer_0] = (float)valid_passes*time_quantum*(float)algorquantum*(float)curvestobin_rxn+time_offset;
		pointer_0++;
		d[pointer_0] = (float)pass_to_fcs.tc0;	//Intensity
		pointer_0++;
		d[pointer_0] = (float)pass_to_fcs.tc1;
		pointer_0++;
		d[pointer_0] = (float)pass_to_fcs.tc2;
		pointer_0++;
		d[pointer_0] = (float)pass_to_fcs.tc0sq;	//IntSD
		pointer_0++;
		d[pointer_0] = (float)pass_to_fcs.tc1sq;
		pointer_0++;
		d[pointer_0] = (float)pass_to_fcs.tc2sq;
		pointer_0++;
		
		for(i=0; i<veclen; i++)   {
			d[pointer_0+i+0*veclen] = gfinal0b0[i];		//	G(tau)0x0
			d[pointer_0+i+1*veclen] = gfinal1b1[i];
			d[pointer_0+i+2*veclen] = gfinal2b2[i];
			d[pointer_0+i+3*veclen] = gfinal0b1[i];
			d[pointer_0+i+4*veclen] = gfinal1b2[i];
			d[pointer_0+i+5*veclen] = gfinal2b0[i];
		}
		fwrite(d, sizeof(double), d_length, fpW);
		
		valid_passes++;
		printf("2FCS Iteration %i\n", valid_passes);
	}

	fclose(fp);		// Close data file
	fclose(fpW);	// Close results G_array file

	// *******	Free Memory ******* //
	free(pass_to_fcs.shatzel0);
	free(pass_to_fcs.shatzel1);
	free(pass_to_fcs.shatzel2);
	free(pass_to_fcs.stemp0);
	free(pass_to_fcs.stemp1);
	free(pass_to_fcs.stemp2);
	free(pass_to_fcs.klaus0);
	free(pass_to_fcs.klaus1);
	free(pass_to_fcs.klaus2);
    free(pass_to_fcs.sc0);
    free(pass_to_fcs.sc1);
    free(pass_to_fcs.sc2);
    free(array0_0);
    free(array0_1);
    free(array0_2);
    free(array1_0);
    free(array1_1);
    free(array1_2);
    free(interleaved);
    free(d);
	free(pass_to_fcs.g0b0);
	free(pass_to_fcs.g1b1);
	free(pass_to_fcs.g2b2);
	free(pass_to_fcs.g0b1);
	free(pass_to_fcs.g1b2);
	free(pass_to_fcs.g2b0);

	pthread_exit(NULL);
}

