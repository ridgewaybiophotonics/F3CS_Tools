/*
 *  complib.h
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
 *  complib.h:
 *  Defines the prototypes for the functions that perform the vectorized 
 *  data decoding described in the main manuscript (Section titled:  
 *  "Data acquisition accomplished by mixing TTL circuitry with 128-bit 
 *  SSE operations," and Fig. 2).  The functions that end with "_c" are 
 *  coded in C in the file fast_sse_ttl_loops.c.  The duplicate functions 
 *  that do not end in "_c" are coded in assembly, and can be found in 
 *  their respective files: bits_2048_16_3.s, bits_2048_32_3.s, 
 *  bits_2048_64_3.s, bits_2048_128_3.s, bits_2048_256_3.s, 
 *  bits_2048_512_3.s, bits_2048_1024_3.s.
 */

/*
All functions take the output from the first pin of a counter chip, and
determine whether a pulse has taken place by subtracting the value
at time t+1 from the value at time t.  For example: 

00001111111100000111111000010001110 (time ->)
...|.......|....|.....|...||..|..|.

( . = no pulse, | = pulse )

The function can correlate up to 8 inputs on lines 0.0-0.7 as it only
looks at one bit at a time.  The data are streamed in in 
    2048 x 8bits    (one word)
    = Fifo size x 8 inputs
 And the functions itterate through 2040 of the 2048 bytes.  The last
 byte can't be analyzed as 
    t(2047) = data(2048) - data(2047)
 and of course data(2048) is not present, so in
 practice the last 8 bytes are written to the 64-bit address rem
 (for "remainder") which are read in the next time the function is called.
 Effectively, this boosts the total number of itterations up to 2048, so
 if the function is called repeatedly, no loss in data will be experienced.
 
 The function is in general able to process 8 time points at a time using
 128-bit registers and the SSE2/3 processor instruction set.  The different
 functions differ in how they they bin the data (Assuming 20Mhz Acq.):
 
 bits_2048_1_3    No binning;         50ns;    output_X[2048]
 bits_2048_16_3   16X binning;       800ns;    output_X[ 128]
 bits_2048_32_3   32X binning;       1.6us;    output_X[  64]
 bits_2048_64_3   64X binning;       3.2us;    output_X[  32]
 bits_2048_128_3  128X binning;      6.4us;    output_X[  16]
 bits_2048_256_3  256X binning;     12.8us;    output_X[   8]
 bits_2048_512_3  512X binning;     25.6us;    output_X[   4]
 bits_2048_1024_3 1024X binning;    51.2us;    output_X[   2]
 bits_2048_2048_3 2048X binning;   102.4us;    output_X[   1]
*/

//	Define National Instruments datatypes
#ifndef _NI_uInt8_DEFINED_
#define _NI_uInt8_DEFINED_ 
    typedef unsigned char uInt8;
#endif

void bits_2048_16_3(unsigned short *one, unsigned short *two, unsigned short *thr, unsigned int *rem, unsigned char *data);
void bits_2048_32_3(unsigned char *one, unsigned char *two, unsigned char *thr, unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_64_3(unsigned char *one, unsigned char *two, unsigned char *thr, unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_128_3(unsigned char *one, unsigned char *two, unsigned char *thr, unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_256_3(unsigned char *one, unsigned char *two, unsigned char *thr, unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_512_3(unsigned char *one, unsigned char *two, unsigned char *thr, unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_1024_3(unsigned char *one, unsigned char *two, unsigned char *thr, unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_2048_3(unsigned char *one, unsigned char *two, unsigned char *thr, unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);

void bits_2048_2048_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_1024_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_512_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_256_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_128_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_64_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_32_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_16_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);
void bits_2048_1_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data);

