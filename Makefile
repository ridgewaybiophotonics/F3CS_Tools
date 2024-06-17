# *********************************************************************
#   Makefile
#
#     o
#    /   Created by William Ridgeway 2007-2012.
#   o
#    \   
#     o
#    /   Copyright 2012 TSRI. All rights reserved.
#   o
#
#	Part of the Triple Correlation Toolbox Suite of programs:
#
#	Ridgeway WK, Millar DP, Williamson JR, Vectorized data acquisition 
#   and fast triple-correlation integrals for Fluorescence Triple 
#   Correlation Spectroscopy, 2012
# *********************************************************************

#CC = -/usr/local/bin/gcc
CC = -/usr/bin/gcc
CFLAGS = -O3 -msse3 -g -fcommon -funroll-all-loops -D_FILE_OFFSET_BITS=64
LFLAGS = -O3 -msse3 -fcommon -funroll-all-loops -D_FILE_OFFSET_BITS=64

#	*******************************************
#	Make categories, organized by the libraries 
#	that are required by each executable
#	*******************************************


all: no_libraries_reqd gsl_reqd nidaqmx_reqd plplot_reqd nidaqmx_AMD_reqd

no_libraries_reqd:	F3CS_Difference F3CS_Outlier2 F3CS_Outlier3 F3CS_Reverser F3CS_TimeTrace F3CS_AxAxA F3CS_AxAxB F3CS_AxBxG F3CS_2FCS 

gsl_reqd:  F3CS_LocalFit F3CS_GlobalFit F3CS_StochasticData

nidaqmx_reqd: F3CS_DAQ F3CS_AlignTool 

plplot_reqd:  F3CS_GlobalFit_3D_Plot

nidaqmx_AMD_reqd: F3CS_DAQ_ASM

docker_static: LFLAGS += -static
docker_static: no_libraries_reqd gsl_reqd



#	*******************************************
#		     Programs, linking
#	*******************************************

F3CS_Difference: Difference.o
	$(CC) Difference.o -o F3CS_Difference -lm $(LFLAGS)
	
F3CS_Outlier2: Outlier2.o
	$(CC) Outlier2.o -o F3CS_Outlier2 -lm $(LFLAGS)
	
F3CS_Outlier3: Outlier3.o
	$(CC) Outlier3.o -o F3CS_Outlier3 -lm $(LFLAGS)
	
F3CS_Reverser: Reverser.o
	$(CC) Reverser.o -o F3CS_Reverser -lm $(LFLAGS)

F3CS_StochasticData: StochasticData.o
	$(CC) StochasticData.o -o F3CS_StochasticData -lgsl -lgslcblas -lm -pthread $(LFLAGS)

F3CS_AlignTool: AlignTool.o fast_sse_ttl_loops.o DAQ_lib.h complib.h
	$(CC) AlignTool.o fast_sse_ttl_loops.o -o F3CS_AlignTool -ldl -lm -pthread -lnidaqmx $(LFLAGS) -lnidaqmx `pkg-config --cflags gtk+-2.0 gthread-2.0` `pkg-config --libs gtk+-2.0 gthread-2.0`

F3CS_DAQ: DAQ.o fast_sse_ttl_loops.o DAQ_lib.h complib.h
	$(CC) DAQ.o fast_sse_ttl_loops.o -o F3CS_DAQ -ldl -lm -pthread -lnidaqmx $(LFLAGS) 

F3CS_DAQ_ASM: DAQ_ASM.o bits_2048_16_3.o DAQ_lib.h complib.h
	$(CC) DAQ_ASM.o bits_2048_16_3.o -o F3CS_DAQ_ASM -ldl -lm -pthread -lnidaqmx $(LFLAGS) 

F3CS_TimeTrace:	TimeTrace.o
	$(CC) TimeTrace.o -o F3CS_TimeTrace -lm $(LFLAGS)

F3CS_AxAxA:	AxAxA.o corr_int_AxAxA.o
	$(CC) AxAxA.o corr_int_AxAxA.o -o F3CS_AxAxA -lm -pthread $(LFLAGS)

F3CS_AxAxB:	AxAxB.o corr_int_AxAxB.o
	$(CC) AxAxB.o corr_int_AxAxB.o -o F3CS_AxAxB -lm -pthread $(LFLAGS)

F3CS_AxBxG:	AxBxG.o corr_int_AxBxG.o
	$(CC) AxBxG.o corr_int_AxBxG.o -o F3CS_AxBxG -lm -pthread $(LFLAGS)

F3CS_2FCS:	2FCS.o corr_int_2FCS.o
	$(CC) 2FCS.o corr_int_2FCS.o -o F3CS_2FCS -lm -pthread $(LFLAGS)

F3CS_LocalFit:	LocalFit.o
	$(CC) LocalFit.o -o F3CS_LocalFit -lgsl -lgslcblas -lm $(LFLAGS)

F3CS_GlobalFit:	GlobalFit.o
	$(CC) GlobalFit.o -o F3CS_GlobalFit -lgsl -lgslcblas -lm $(LFLAGS)
	
F3CS_GlobalFit_3D_Plot:	GlobalFit_3D_Plot.o
	$(CC) GlobalFit_3D_Plot.o -o F3CS_GlobalFit_3D_Plot -lgsl -lgslcblas -lplplotd -lm $(LFLAGS)
	
#	*******************************************
#		     Objects, compilation
#	*******************************************

Difference.o:	Difference.c libF3CS.h
	$(CC) -c Difference.c -o Difference.o $(CFLAGS)
	
Outlier2.o:	Outlier2.c libF3CS.h
	$(CC) -c Outlier2.c -o Outlier2.o $(CFLAGS)
	
Outlier3.o:	Outlier3.c libF3CS.h
	$(CC) -c Outlier3.c -o Outlier3.o $(CFLAGS)
		
Reverser.o:	Reverser.c libF3CS.h
	$(CC) -c Reverser.c -o Reverser.o $(CFLAGS)
		
StochasticData.o:	StochasticData.c
	$(CC) -c StochasticData.c -o StochasticData.o $(CFLAGS)

fast_sse_ttl_loops.o: fast_sse_ttl_loops.c complib.h
	$(CC) -c fast_sse_ttl_loops.c -fprefetch-loop-arrays -funroll-all-loops -ftree-vectorize -msse3 -O3
#	$(CC) -c fast_sse_ttl_loops.c -fprefetch-loop-arrays -funroll-all-loops -ftree-vectorize -ftree-vectorizer-verbose=2 -msse3 -O3
	
AlignTool.o: AlignTool.c  DAQ_lib.h complib.h
	$(CC) -c AlignTool.c -o AlignTool.o -funroll-all-loops -msse3 -O2 -Wall -D_FILE_OFFSET_BITS=64 `pkg-config --cflags gtk+-2.0 gthread-2.0`

DAQ.o: DAQ.c  DAQ_lib.h complib.h
	$(CC) -c DAQ.c -o DAQ.o -funroll-all-loops -msse3 -O2 -Wall -D_FILE_OFFSET_BITS=64
	
DAQ_ASM.o: DAQ_ASM.c  DAQ_lib.h complib.h
	$(CC) -c DAQ_ASM.c -o DAQ_ASM.o -funroll-all-loops -msse3 -O2 -Wall -D_FILE_OFFSET_BITS=64

TimeTrace.o:	TimeTrace.c
	$(CC) -c TimeTrace.c -o TimeTrace.o $(CFLAGS)

AxAxA.o:	AxAxA.c libF3CS.h
	$(CC) -c AxAxA.c -o AxAxA.o -pthread $(CFLAGS)
	
AxAxB.o:	AxAxB.c libF3CS.h
	$(CC) -c AxAxB.c -o AxAxB.o -pthread $(CFLAGS)
	
AxBxG.o:	AxBxG.c libF3CS.h
	$(CC) -c AxBxG.c -o AxBxG.o -pthread $(CFLAGS)
	
corr_int_AxAxA.o:	corr_int_AxAxA.c libF3CS.h
	$(CC) -c corr_int_AxAxA.c -o corr_int_AxAxA.o -ftree-vectorize -pthread $(CFLAGS)
#	$(CC) -c corr_int_AxAxA.c -o corr_int_AxAxA.o -ftree-vectorize -ftree-vectorizer-verbose=2 -pthread $(CFLAGS)
	
corr_int_AxAxB.o:	corr_int_AxAxB.c libF3CS.h
	$(CC) -c corr_int_AxAxB.c -o corr_int_AxAxB.o -ftree-vectorize -pthread $(CFLAGS)
#	$(CC) -c corr_int_AxAxB.c -o corr_int_AxAxB.o -ftree-vectorize -ftree-vectorizer-verbose=2 -pthread $(CFLAGS)

corr_int_AxBxG.o:	corr_int_AxBxG.c libF3CS.h
	$(CC) -c corr_int_AxBxG.c -o corr_int_AxBxG.o -ftree-vectorize -pthread $(CFLAGS)
#	$(CC) -c corr_int_AxBxG.c -o corr_int_AxBxG.o -ftree-vectorize -ftree-vectorizer-verbose=2 -pthread $(CFLAGS)

2FCS.o:	2FCS.c libF3CS.h
	$(CC) -c 2FCS.c -o 2FCS.o -ftree-vectorize -pthread $(CFLAGS)
#	$(CC) -c 2FCS.c -o 2FCS.o -ftree-vectorize -ftree-vectorizer-verbose=2 -pthread $(CFLAGS)

corr_int_2FCS.o:	corr_int_2FCS.c libF3CS.h
	$(CC) -c corr_int_2FCS.c -o corr_int_2FCS.o -ftree-vectorize -pthread $(CFLAGS)
#	$(CC) -c corr_int_2FCS.c -o corr_int_2FCS.o -ftree-vectorize -ftree-vectorizer-verbose=2 -pthread $(CFLAGS)

LocalFit.o:	LocalFit.c libF3CS.h
	$(CC) -c LocalFit.c -o LocalFit.o $(CFLAGS)

GlobalFit.o:	GlobalFit.c libF3CS.h
	$(CC) -c GlobalFit.c -o GlobalFit.o $(CFLAGS)
	
GlobalFit_3D_Plot.o:	GlobalFit_3D_Plot.c libF3CS.h
	$(CC) -c GlobalFit_3D_Plot.c -o GlobalFit_3D_Plot.o $(CFLAGS)

bits_2048_16_3.o: bits_2048_16_3.s complib.h
	$(CC) -c bits_2048_16_3.s  -msse3

clean:
	-rm -f *.o F3CS_Difference F3CS_Outlier2 F3CS_Outlier3 F3CS_Reverser F3CS_TimeTrace F3CS_AxAxA F3CS_AxAxB F3CS_AxBxG F3CS_LocalFit F3CS_GlobalFit F3CS_GlobalFit_3D_Plot F3CS_2FCS F3CS_AlignTool F3CS_DAQ F3CS_DAQ_ASM F3CS_StochasticData
