/*
 *  DAQ_ASM.c
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
 *  DAQ_ASM.c
 *  Identical to DAQ.c, but for the fact that this file calls the 
 *  assembly-code versions of the fast TTL/SSE loops in the 
 *  bits_2048_16_3.s (etc) series of files.
 *
 */

#include <stdio.h>
#include <pthread.h>
#include <NIDAQmx.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include "complib.h"

#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else
#define	algorquant		131072
#define TUNE_BUFFER		2048

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
	int	inputchar;
	do	{
		inputchar = getchar();
		if((char) inputchar != 'q')	{
			printf("Press 'q' to halt data acquisition (char=%s)\n", (char*) &inputchar);
		}
	}	while ((char) inputchar != 'q');
	NICard.run = 0;
	printf("* * * * * * * * * * * *\nAcquisition Stopping...\n* * * * * * * * * * * *\n");
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

	if(argc != 4)  {
        printf("Usage: F3CS_DAQ_ASM tag_string time(s) maxcounts(kHz, int)\ne.g. >rxn_3c expt1 1000 800\n");
        return 1;
    }
	
	int     error=0;
	int     h,l,m;
    int     i = 0;
    int     clock_freq = 20000000;
	TaskHandle	taskHandle=0;
	TaskHandle	taskHandle3=0;
	uInt8   shutters[2]     = {0, 1}; /* {closed, open} */
	int32       *written = 0;
	char        errBuff[2048]={'\0'};
	double tim = 0;
	sscanf(argv[2],"%i", &i);
	double max_tim = clock_freq*(float)i;
	sscanf(argv[3],"%i", &i);
	double max_apd_counts = (double)i;
	max_apd_counts = max_apd_counts / 1000.0;
	printf("Laser will be shuttered > %.6f MHz\n", max_apd_counts);

	/*	Launch the "press q to quit" thread   */
	int rc;
	pthread_t getchar_thread;
	pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    struct paramAcq pass_to_Int;    
   	NICard.run = 1;
	printf("Press 'q' to halt data acquisition...\n");
    rc = pthread_create(&getchar_thread, &attr, getchar_fxn, (void *) &pass_to_Int);
	if (rc){
		printf("ERROR; return code from pthread_create() is %d\n", rc);
		exit(-1);
	}

    uInt8   data[2048];
    unsigned char  *scratch_one;
    unsigned char  *scratch_two;
    unsigned char  *scratch_thr;
    unsigned short *result_one;
    unsigned short *result_two;
    unsigned short *result_thr;
    unsigned short burst_data_uh0;
    unsigned short burst_data_uh1;
    unsigned short burst_data_uh2;
    double burst_data_d0;
    double burst_data_d1;
    double burst_data_d2;
    unsigned short *circ_buffer1;
    unsigned int   *rem;
    scratch_one= (unsigned char  *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    scratch_two= (unsigned char  *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    scratch_thr= (unsigned char  *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    circ_buffer1= (unsigned short *) calloc(3*128*algorquant,sizeof(unsigned short));
    rem = (unsigned int *) calloc(2,sizeof(unsigned int));
    
	int32		sampsRead;
	DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
    DAQmxErrChk (DAQmxCreateDIChan(taskHandle,"Dev1/port0/line0:7", "" , DAQmx_Val_ChanForAllLines));
    DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"/Dev1/20MhzTimebase", clock_freq, DAQmx_Val_Falling, DAQmx_Val_ContSamps, TUNE_BUFFER));
    rem[0] = 0x0;
    rem[1] = 0x0;
    l = 0;
    h = 0;
    
	/*	Run Shutter of Counter 0 Src (37 hot, 4 ground)    */
	DAQmxErrChk (DAQmxCreateTask("Shutter Task",&taskHandle3));
    DAQmxErrChk (DAQmxCreateDOChan(taskHandle3,"Dev1/port2/line0", "", DAQmx_Val_ChanForAllLines));

	struct write_stuff val;
	val.filenum = 0;
	val.key = 0;
	val.array = circ_buffer1;
	char	temp_char[512];
	sprintf(temp_char, "%s", argv[1]);
	val.tag = temp_char;
	pthread_t write_thread;

	/*    Purge buffers    */
	for(i=0; i<TUNE_BUFFER; i++)	data[i] = 0x0;
	bits_2048_16_3(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);

	int	writes = 0;
	int	failed_writes = 0;
	
	int	shutter_timer = -1;
	int error_manual=0;
	NICard.run = 1;
	
	DAQmxErrChk (DAQmxStartTask(taskHandle3));
	
	/*	Open Shutter, wait 1s   */
	printf("~~~~~ Opening Shutters ~~~~~~~~~~~~~*\n");
	DAQmxErrChk (DAQmxWriteDigitalU8(taskHandle3, 1, 0, 10, DAQmx_Val_GroupByScanNumber, &shutters[1], written, NULL)); /* Open Shutter, APDs	*/
	sleep(1);
	
	/*	loop structure (for 96MB files, 1.0/13.42s)
		While time < time_max
			Write to disk: 128*algorquantum
				Check APD:	1 algorquantum (~10Hz)
					Read from card:	TUNE_BUFFER/16 samples 1/(16*64) algorquantum
	*/
	 
	DAQmxErrChk (DAQmxStartTask(taskHandle));
	while (NICard.run == 1 && tim<max_tim)     {
		l = 0;
		for(m=0; m<128; m++)  {
			shutter_timer--;
			burst_data_d0 = 0;
			burst_data_d1 = 0;
			burst_data_d2 = 0;
			for(h=0; h<16*64; h++)	{
				burst_data_uh0 = 0;
				burst_data_uh1 = 0;
				burst_data_uh2 = 0;
				error_manual = DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL);
				if(error_manual < 0)	{	/*    Error Handling Time    */
					printf("Card Overrun.  Pausing briefly.  \n If Overruns occur too frequently, quit background processes \n and avoid scrolling windows, etc.\n");
					failed_writes++;
					if(taskHandle!=0)	{
						DAQmxStopTask(taskHandle);
						DAQmxClearTask(taskHandle);
					}
					else printf("* * * * * TaskHandle Fried * * * * *\n");
					DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
					DAQmxErrChk (DAQmxCreateDIChan(taskHandle,"Dev1/port0/line0:7", "" , DAQmx_Val_ChanForAllLines));
					DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"/Dev1/20MhzTimebase", clock_freq, DAQmx_Val_Falling, DAQmx_Val_ContSamps, TUNE_BUFFER));
					
					DAQmxGetExtendedErrorInfo(errBuff,2048);
					printf("DAQmx Error: %s\n",errBuff);
					for(i=l; i<3*128*algorquant; i++)	circ_buffer1[i] = 0;  /*	Zero-pad data array, bail from loop    */
					l = 3*128*algorquant;
					h = 1024*1024;
					m = 1024*1024;
					sleep(2);
					DAQmxStartTask(taskHandle);
					printf("task restarted\n");
				}
				else	{
					bits_2048_16_3(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);     
					for(i=0; i< 128; i++) {  
						burst_data_uh0 += result_one[i];
						circ_buffer1[l] = result_one[i];
						l++;
						burst_data_uh1 += result_two[i];
						circ_buffer1[l] = result_two[i];
						l++;
						burst_data_uh2 += result_thr[i];	
						circ_buffer1[l] = result_thr[i];
						l++;
					} 
					burst_data_d0 = burst_data_d0 + burst_data_uh0;
					burst_data_d1 = burst_data_d1 + burst_data_uh1;
					burst_data_d2 = burst_data_d2 + burst_data_uh2;
				}
			}			
			/*	Check for APD Saturation (@10Hz)    */
			if((burst_data_d0*0.00000953674316406 > max_apd_counts || burst_data_d1*0.00000953674316406 > max_apd_counts || burst_data_d2*0.00000953674316406 > max_apd_counts || shutter_timer < -4990)&&shutter_timer<0)	{	/* Max counts exceeded, or every 500s */
				DAQmxErrChk (DAQmxWriteDigitalU8(taskHandle3, 1, 0, 10, DAQmx_Val_GroupByScanNumber, &shutters[0], written, NULL)); /* Close Shutter  */
				printf("~~~~~ Closing Shutters ~~~~~*|        \n  Ct (MHz): %7.6f\t%7.6f\t%7.6f\n", burst_data_d0*0.00000953674316406, burst_data_d1*0.00000953674316406, burst_data_d2*0.00000953674316406);
				if(shutter_timer > -10)	{	/*	If shutter closed less than ~1s after reopening   */
					printf("**** Consistent high signal.  Shuttered for 5s ****\n");
					shutter_timer = 50;
				}
				else if (shutter_timer < -4990)	{
					printf("**** Routine Periodic Shutter closure.  Shuttered for 5s ****\n");
					shutter_timer = 50;
				}
				else	{
					shutter_timer = 10;	/*	count down 10 cycles to next opening (~1.0s)	*/
				}
			}
			if(shutter_timer == 0)	{
				printf("~~~~~ Opening Shutters ~~~~~~~~~~~~~*\n");
				DAQmxErrChk (DAQmxWriteDigitalU8(taskHandle3, 1, 0, 10, DAQmx_Val_GroupByScanNumber, &shutters[1], written, NULL)); /* Open Shutter	*/
			}
			if(shutter_timer%10 == 0)	{
				printf("Time:\t%6.2fs   %i   Ct (MHz): %7.6f\t%7.6f\t%7.6f\n",tim/clock_freq, val.key, burst_data_d0*0.00000953674316406, burst_data_d1*0.00000953674316406, burst_data_d2*0.00000953674316406);
			}
			tim = tim + 16*algorquant;
		}

		/*	Decide whether to write or stall	*/
		if(val.key == 0) {
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
		}
		else    {
			DAQmxStopTask(taskHandle);
			do {	
				printf("Overflow... sleeping for 1s\n");
				sleep(1);
			} while (val.key != 0);
			DAQmxStartTask(taskHandle);
		}
	}

	printf("~~~~~ Closing Shutters ~~~~~*|        \n");
	DAQmxErrChk (DAQmxWriteDigitalU8(taskHandle3, 1, 0, 10, DAQmx_Val_GroupByScanNumber, &shutters[0], written, NULL)); /*	Close Shutters	*/

    Error:
	if(DAQmxFailed(error))   printf("DAQmx Error: %s\n",errBuff);
	if(taskHandle!=0)	{
		printf("Clearing task...\n");
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
		DAQmxStopTask(taskHandle3);
		DAQmxClearTask(taskHandle3);
	}
    while (val.key!=0)   {
        printf("Waiting for the file I/O to complete before closing...\n");
        sleep(1);
    }

    /*	Clean up and Exit	*/
    NICard.run = 1;
    printf("NI Card stopped.\n");
    free(rem);
    free(scratch_one);
    free(scratch_two);
    free(scratch_thr);
    free(result_one);
    free(result_two);
    free(result_thr);
	free(circ_buffer1);
	pthread_attr_destroy(&attr);
	    
    printf("Program completed.\n");
	return 0;
}
