/*
 *  fast_sse_ttl_loops.c
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
 *  fast_sse_ttl_loops.c
 *  C versions of the assembly code files described below.  The functions 
 *  were written to be highly amenable to auto-vectorization using GCC, 
 *  and developers will likely want to take advantage of the commented 
 *  out instructions in the makefile that enable a verbose accounting 
 *  of which loops were and were not vectorized.
 *
 */

#include <stdint.h>
#include <stdlib.h>

#define	div1024	2
#define	div512	4
#define	div256	8
#define	div128	16
#define	div64	32
#define	div32	64
#define	div16	128
#define	div8	256
#define	div1	2048

void bits_2048_1_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	for(j=0; j<2048; j++)		output_one[j] =  onn[j]&0x1;
	for(j=0; j<2048; j++)		output_two[j] = (onn[j]&0x4) >> 2;
	for(j=0; j<2048; j++)		output_thr[j] = (onn[j]&0x10) >> 4;
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_1024_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	i,j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2;
	
	for(i=0; i<div1024; i++)	{
		total2 = 0;
		for(j=i*2048/div1024; j<(i+1)*2048/div1024; j++)		total2 += onn[j]&0x1;
		output_one[i] = total2;
	}
	
	for(i=0; i<div1024; i++)	{
		total2 = 0;
		for(j=i*2048/div1024; j<(i+1)*2048/div1024; j++)		total2 += onn[j]&0x4;
		output_two[i] = total2 >> 2;
	}
	
	for(i=0; i<div1024; i++)	{
		total2 = 0;
		for(j=i*2048/div1024; j<(i+1)*2048/div1024; j++)		total2 += onn[j]&0x10;
		output_thr[i] = total2 >> 4;
	}
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_512_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	i,j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];

	
	unsigned int total2;
	
	for(i=0; i<div512; i++)	{
		total2 = 0;
		for(j=i*2048/div512; j<(i+1)*2048/div512; j++)		total2 += onn[j]&0x1;
		output_one[i] = total2;
	}
	
	for(i=0; i<div512; i++)	{
		total2 = 0;
		for(j=i*2048/div512; j<(i+1)*2048/div512; j++)		total2 += onn[j]&0x4;
		output_two[i] = total2 >> 2;
	}
	
	for(i=0; i<div512; i++)	{
		total2 = 0;
		for(j=i*2048/div512; j<(i+1)*2048/div512; j++)		total2 += onn[j]&0x10;
		output_thr[i] = total2 >> 4;
	}
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_256_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	i,j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2;
	
	for(i=0; i<div256; i++)	{
		total2 = 0;
		for(j=i*2048/div256; j<(i+1)*2048/div256; j++)		total2 += onn[j]&0x1;
		output_one[i] = total2;
	}
	
	for(i=0; i<div256; i++)	{
		total2 = 0;
		for(j=i*2048/div256; j<(i+1)*2048/div256; j++)		total2 += onn[j]&0x4;
		output_two[i] = total2 >> 2;
	}
	
	for(i=0; i<div256; i++)	{
		total2 = 0;
		for(j=i*2048/div256; j<(i+1)*2048/div256; j++)		total2 += onn[j]&0x10;
		output_thr[i] = total2 >> 4;
	}
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_128_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	i,j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2;
	
	for(i=0; i<div128; i++)	{
		total2 = 0;
		for(j=i*2048/div128; j<(i+1)*2048/div128; j++)		total2 += onn[j]&0x1;
		output_one[i] = total2;
	}
	
	for(i=0; i<div128; i++)	{
		total2 = 0;
		for(j=i*2048/div128; j<(i+1)*2048/div128; j++)		total2 += onn[j]&0x4;
		output_two[i] = total2 >> 2;
	}
	
	for(i=0; i<div128; i++)	{
		total2 = 0;
		for(j=i*2048/div128; j<(i+1)*2048/div128; j++)		total2 += onn[j]&0x10;
		output_thr[i] = total2 >> 4;
	}
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_64_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	i,j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2;
	
	for(i=0; i<div64; i++)	{
		total2 = 0;
		for(j=i*2048/div64; j<(i+1)*2048/div64; j++)		total2 += onn[j]&0x1;
		output_one[i] = total2;
	}
	
	for(i=0; i<div64; i++)	{
		total2 = 0;
		for(j=i*2048/div64; j<(i+1)*2048/div64; j++)		total2 += onn[j]&0x4;
		output_two[i] = total2 >> 2;
	}
	
	for(i=0; i<div64; i++)	{
		total2 = 0;
		for(j=i*2048/div64; j<(i+1)*2048/div64; j++)		total2 += onn[j]&0x10;
		output_thr[i] = total2 >> 4;
	}
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_32_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	i,j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2;
	
	for(i=0; i<div32; i++)	{
		total2 = 0;
		for(j=i*2048/div32; j<(i+1)*2048/div32; j++)		total2 += onn[j]&0x1;
		output_one[i] = total2;
	}
	
	for(i=0; i<div32; i++)	{
		total2 = 0;
		for(j=i*2048/div32; j<(i+1)*2048/div32; j++)		total2 += onn[j]&0x4;
		output_two[i] = total2 >> 2;
	}
	
	for(i=0; i<div32; i++)	{
		total2 = 0;
		for(j=i*2048/div32; j<(i+1)*2048/div32; j++)		total2 += onn[j]&0x10;
		output_thr[i] = total2 >> 4;
	}
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_16_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	i,j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
								onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2;
	
	for(i=0; i<div16; i++)	{
		total2 = 0;
		for(j=i*2048/div16; j<(i+1)*2048/div16; j++)		total2 += onn[j]&0x1;
		output_one[i] = total2;
	}
	
	for(i=0; i<div16; i++)	{
		total2 = 0;
		for(j=i*2048/div16; j<(i+1)*2048/div16; j++)		total2 += onn[j]&0x4;
		output_two[i] = total2 >> 2;
	}
	
	for(i=0; i<div16; i++)	{
		total2 = 0;
		for(j=i*2048/div16; j<(i+1)*2048/div16; j++)		total2 += onn[j]&0x10;
		output_thr[i] = total2 >> 4;
	}
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_2048_3_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned int *rem, unsigned char *data) {
	int	j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x1;
	output_one[0] = total2;

	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x4;
	output_two[0] = total2 >> 2;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x10;
	output_thr[0] = total2 >> 4;
	
	free(onn);
	free(twn);

	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}

void bits_2048_2048_8_c(unsigned short *output_one, unsigned short *output_two, unsigned short *output_thr, unsigned short *output_fou, unsigned short *output_fiv, unsigned short *output_six, unsigned short *output_sev, unsigned short *output_eig,unsigned int *rem, unsigned char *data) {
	int	j;
	unsigned char * remuc = (unsigned char *) rem;
	
	unsigned char * twn;
	twn = malloc(2048*sizeof(unsigned char));	//Vestigial code, but including it somehow speeds up the main loop...
	unsigned char * onn;
	onn = malloc(2048*sizeof(unsigned char));
	
	for(j=0; j<7; j++)			onn[j]   = remuc[j+1]^remuc[j];
	onn[7]   = data[0]^remuc[7];
	for(j=0; j<2040; j++)		onn[j+8] = data[j+1]^data[j];
	
	unsigned int total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x1;
	output_one[0] = total2;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x2;
	output_two[0] = total2 >> 1;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x4;
	output_thr[0] = total2 >> 2;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x8;
	output_fou[0] = total2 >> 3;
	
	for(j=0; j<2048; j++)		onn[j]  = onn[j] >> 4;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x10;
	output_fiv[0] = total2;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x20;
	output_six[0] = total2 >> 1;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x40;
	output_sev[0] = total2 >> 2;
	
	total2 = 0;
	for(j=0; j<2048; j++)		total2 += onn[j]&0x80;
	output_eig[0] = total2 >> 3;
	
	free(onn);
	free(twn);
	
	for(j=0; j<8; j++)	remuc[j] = data[2040+j];
	
    return;
}
