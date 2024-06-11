/*
 *  README.txt
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
 */
 
/* ===== List of programs in the Triple Correlation Toolbox ====== */

16 executables comprise the Triple Correlation Toolbox package:

Data acquisition
	F3CS_AlignTool  (Real-time tuning of microscope)
	F3CS_DAQ	(Robust data acquisition)
	F3CS_DAQ_ASM	(Assembly code version of F3CS_DAQ)

Simulated data acquisition for tutorial
	F3CS_StochasticData

Primary data manipulation
	F3CS_Reverser	(Reverse time-order of data)
	F3CS_TimeTrace	(Downsample data to lower time resolution)

Correlation integrals
	F3CS_2FCS		(FCS / FCCS correlation integrals)
	F3CS_AxAxA		(F3CS correlation integrals)
	F3CS_AxAxB		(F3CS correlation integrals)
	F3CS_AxBxG		(F3CS correlation integrals)

Outlier rejection
	F3CS_Outlier2		(For FCS / FCCS data)
	F3CS_Outlier3		(For F3CS data)
	F3CS_Difference	(Computes Delta G(tau1,tau2))

Data fitting
	F3CS_LocalFit	(Fit individual FCS / FCCS curves)
	F3CS_GlobalFit	(Globally fit FCS, FCCS and F3CS data)
	F3CS_GlobalFit_3D_Plot	(Provides extra graphics output)

/* ======================== Installation ========================= */

The F3CS_Manual.pdf contains complete installation instructions and library
dependencies, but installation follows these lines:

(1)  cd F3CS_Tools

(2)  One of the following 6 commands, depending on library availability:
make no_libraries_reqd
make gsl_reqd				(Requires Gnu Scientific Library (gsl) )
make nidaqmx_reqd			(Requires NIDAQmx and gtk+-2.0)
make plplot_reqd	  		(Requires plplot)
make nidaqmx_AMD_reqd:	(Requires NIDAQmx and an AMD CPU)
make all					(Requires  gsl, NIDAQmx, GTK +-2.x, plplot)

(3)  sudo cp F3CS_* /usr/local/bin		(Or any other location in the user's path)

						-----------------------

Programs complied by each make method:

make no_libraries_reqd:	
	F3CS_Difference 
	F3CS_Outlier2 
	F3CS_Outlier3 
	F3CS_Reverser 
	F3CS_TimeTrace 
	F3CS_AxAxA 
	F3CS_AxAxB 
	F3CS_AxBxG 
	F3CS_2FCS 

make gsl_reqd:	
	F3CS_LocalFit 
	F3CS_GlobalFit 
	F3CS_StochasticData

make nidaqmx_reqd:
	F3CS_DAQ 
	F3CS_AlignTool 

make plplot_reqd:
	F3CS_GlobalFit_3D_Plot 

make nidaqmx_AMD_reqd:
	F3CS_DAQ_ASM

make all:
	--all programs--
    
Note:  On various Linux distributions, several programs (F3CS_StochasticData, F3CS_2FCS, F3CS_AxAxA, F3CS_AxAxB, F3CS_AxBxG) generate segmentation faults at the end of their runs due to known errors in the pthread library.  The faults occur at the end of the programs, do not effect data integrity, and leak about 100 bytes of memory per fault.

/* ==================== Test run & Tutorial ====================== */

A tutorial, a means to generate sample data and sample output are 
provided in the user manual:  F3CS_Manual.pdf

/* =============== List of compilation dependencies ===============*/

Each of the executables is built using the following files:

F3CS_Difference: 
    Difference.c 
    libF3CS.h
	
F3CS_Outlier2: 
    Outlier2.c 
    libF3CS.h
	
F3CS_Outlier3: 
    Outlier3.c 
    libF3CS.h
	
F3CS_Reverser: 
    Reverser.c 
    libF3CS.h

F3CS_StochasticData: 
    StochasticData.c 

F3CS_AlignTool: 
    AlignTool.c  
    DAQ_lib.h 
    complib.h 
    fast_sse_ttl_loops.c

F3CS_DAQ: 
    DAQ.c
    DAQ_lib.h 
    complib.h 
    fast_sse_ttl_loops.c 

F3CS_DAQ_ASM: 
    DAQ_ASM.c  
    DAQ_lib.h 
    complib.h
    bits_2048_16_3.s

F3CS_TimeTrace:	
    TimeTrace.c

F3CS_AxAxA:	
    AxAxA.c 
    libF3CS.h 
    corr_int_AxAxA.c

F3CS_AxAxB:	
    AxAxB.c 
    libF3CS.h
    corr_int_AxAxB.c

F3CS_AxBxG:	
    AxBxG.c 
    libF3CS.h
    corr_int_AxBxG.o

F3CS_2FCS:	
    2FCS.c 
    libF3CS.h
    corr_int_2FCS.c

F3CS_LocalFit:	
    LocalFit.c 
    libF3CS.h

F3CS_GlobalFit:	
    GlobalFit.c 
    libF3CS.h
	
F3CS_GlobalFit_3D_Plot:	
    GlobalFit_3D_Plot.c 
    libF3CS.h

*/ =========== Description of each file (for developers) ========= */

-- Header files --
complib.h
    Defines the prototypes for the functions that perform the vectorized data decoding described in the main manuscript (Section titled:  "Data acquisition accomplished by mixing TTL circuitry with 128-bit SSE operations," and Fig. 2).  The functions that end with "_c" are coded in C in the file fast_sse_ttl_loops.c.  The duplicate functions that do not end in "_c" are coded in assembly, and can be found in their respective files: bits_2048_16_3.s, bits_2048_32_3.s, bits_2048_64_3.s, bits_2048_128_3.s, bits_2048_256_3.s, bits_2048_512_3.s, bits_2048_1024_3.s.


DAQ_lib.h
	Defines prototypes and memory sizes for the various data acquisition routines in F3CS_AlignTool, F3CS_DAQ and F3CS_DAQ_ASM.  The window control parameters for F3CS_AlignTool are mainly set in the last half of this file.  The pmax and n values defined here only affect the correlation integrals calculated in the "FCS" tab of F3CS_AlignTool, and do not impact later data processing.

libF3CS.h
	Sets the prototypes, memory structures and constants used by the entire suite.  Of particular note, the n and pmax values for FCS integrals and n_GGG and pmax_GGG values for F3CS integrals are set here on lines 30-34.  (The programs must be recompiled (make clean; make) for these changes to have effect.)

-- C source code files --

2FCS.c
	This is the main function for the F3CS_2FCS program, and handles the busy work by performing file i/o, launching multithreading and calling the correlation integrals coded in corr_int_2FCS.c.

AlignTool.c
    This file contains code to generate the user interface of F3CS_AlignTool, as well as to perform data acquisition and data processing.  

{AxAxA.c, AxAxB.c, AxBxG.c}
	These three files produce F3CS_AxAxA, F3CS_AxAxB, and F3CS_AxBxG respectively, and differ only by lines 37-39, where only one of three variables is defined.  Compile-time #ifdef statements throughout the rest of the file then use this information to determine which of the three programs is produced.  Like 2FCS.c, these files conduct the busy work of the correlation integral process and make repeated calls to the correlation integrals located in corr_int_AxAxA.c, corr_int_AxAxB.c, and corr_int_AxBxG.c.

DAQ.c
    This code provides the commandline interface for the fast, robust data acquisition process.  Calls the TTL/SSE loops from fast_sse_ttl_loops.c.
    
DAQ_ASM.c
    Identical to DAQ.c, but for the fact that this file calls the assembly-code versions of the fast TTL/SSE loops in the bits_2048_16_3.s (etc) series of files.

Difference.c
	Self-contained code for computing Delta G(tau_1,tau_2) as described in the manuscript section "Time reversal asymmetry detection for analysis of irreversible processes."

GlobalFit.c
	Self-contained code to generate F3CS_GlobalFit.

GlobalFit_3D_Plot.c
	The same as GlobalFit.c, but with extra routines that generate 3D plots (at the end of this file)

LocalFit.c
	Self-contained code to generate F3CS_LocalFit

Outlier2.c
	Self-contained code to generate F3CS_Outlier2.

Outlier3.c
	Self-contained code to generate F3CS_Outlier3.

Reverser.c
	Self-contained code to generate F3CS_Reverser.

StochasticData.c
	Self-contained code to generate F3CS_StochasticData.  The rates displayed in the tutorial figure 1 can be set by altering lines 120-125.  The fake optical artifact can be eliminated on lines 226-229.
	
TimeTrace.c
	Self-contained code to generate F3CS_TimeTrace.

corr_int_2FCS.c
	FCS correlation integrals implementing Eq. 4 in the main manuscript.  Called by F3CS_2FCS.
    
{corr_int_AxAxA.c, corr_int_AxAxB.c, corr_int_AxBxG.c}
	F3CS correlation integrals implementing Eq. 6 in the main manuscript.  Called by F3CS_AxAxA, F3CS_AxAxB and F3CS_AxBxG respectively.

fast_sse_ttl_loops.c
	C versions of the assembly code files described below.  The functions were written to be highly amenable to auto-vectorization using GCC, and developers will likely want to take advantage of the commented out instructions in the makefile that enable a verbose accounting of which loops were and were not vectorized.

-- Assembly code files --

{bits_2048_16_3.s, bits_2048_32_3.s, bits_2048_64_3.s, bits_2048_128_3.s, bits_2048_256_3.s, bits_2048_512_3.s, bits_2048_1024_3.s}
	Each of these files code for a single function that decodes data to complete the TTL/SSE scheme described in the main manuscript (Section titled:  "Data acquisition accomplished by mixing TTL circuitry with 128-bit SSE operations," and Fig. 2).  Each function processes 2048 samples of data, but bins them in different ways described by the second number in the filename: i.e., bits_2048_2048_3.s bins all 2048 samples together to return a single timepoint, while bits_2048_16_3.s bins every 16 data points together and returns 2048/16 = 128 timepoints.  Slower, but equivalent C code is provided in fast_sse_ttl_loops.c.


-- Sample control scripts --

input.global.txt
	Prototype input file for the program F3CS_Global

input.local.txt
	Prototype input file for the program F3CS_Local

Makefile
	Directions for compiling using the make set of tools.


