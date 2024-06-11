/*
 *  DAQ_lib.h
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
 *  DAQ_lib.h
 *  Defines prototypes and memory sizes for the various data acquisition routines in
 *  F3CS_AlignTool, F3CS_DAQ and F3CS_DAQ_ASM.  The window control parameters for
 *  F3CS_AlignTool are mainly set in the last half of this file.  The pmax and n values
 *  defined here only affect the correlation integrals calculated in the "FCS" tab of 
 *  F3CS_AlignTool, and do not impact later data processing.
 *
 */

#include <pthread.h>

/*  Combine a list of protoypes for linking */

#ifndef _NI_uInt8_DEFINED_
#define _NI_uInt8_DEFINED_ 
    typedef unsigned char uInt8;
#endif

int bits_to_words(unsigned short *one, unsigned short *two, unsigned short *thr, int len, unsigned int *remainder, unsigned char *data);
//  Treat the 64-bit remainder as an unsigned array of two 32-bit integers.

#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else
#define NUM_COLOURS 2
#define FRAME_WIDTH 900
#define BURST_FRAME_WIDTH 1024
#define NUM_GRAPH_POINTS 301   //Should equal FRAME_WIDTH / integer +1.
#define TUNE_BUFFER 2048
#define TEST_BUFFER 2048
#define TUNE_BINNING 512
#define tune_points_to_average_initial_value 50
#define LONG_AVERAGE_TIME 100
#define FCS_BUFFER  2048
#define BURST_BUFFER 1024


// pmax and n in this file do not effect dataprocessing in F3CS_2FCS or F3CS_AxBxG (etc), they are for the sole purpose of the FCS tab in F3CS_AlignToolß
#define pmax       12
#define n          32


#define samplesize 131073  //2^13*16 + 1
#define rambuffer  100
#define tuneDataBufferSize 32
#define tuneDataDivisor 4
#ifndef _NI_uInt8_DEFINED_
#define _NI_uInt8_DEFINED_ 
typedef unsigned char uInt8;
#endif

#define rambuffer  100
#define samplesize 131073
#define curvestobin 64
#define anticipatedtime 12
#define	veclen (n*(1+pmax))
#define	circ_buffer_size 6


unsigned short  *array0_0;
unsigned short  *array0_1;
unsigned short	*array0_2;
unsigned short  *array1_0;
unsigned short  *array1_1;
unsigned short	*array1_2;

float           g0b0[veclen];
float           g0b1[(veclen)];
float           g1b0[(veclen)];
float           g1b1[(veclen)];
float           g2b2[(veclen)];
float			gtau_min = 0.0f;
float			gtau_max = 2.0f;
unsigned short 	klaus0_A[pmax];
unsigned short 	klaus1_A[pmax];
unsigned short 	klaus2_A[pmax];
unsigned short 	shatzel0_A[(veclen)];
unsigned short 	shatzel1_A[(veclen)];
unsigned short 	shatzel2_A[(veclen)];
unsigned short 	klaus0_B[pmax];
unsigned short 	klaus1_B[pmax];
unsigned short 	klaus2_B[pmax];
unsigned short 	shatzel0_B[(veclen)];
unsigned short 	shatzel1_B[(veclen)];
unsigned short 	shatzel2_B[(veclen)];
unsigned int 	stemp0_A[pmax];
unsigned int 	stemp1_A[pmax];
unsigned int 	stemp2_A[pmax];
unsigned int 	stemp0_B[pmax];
unsigned int 	stemp1_B[pmax];
unsigned int 	stemp2_B[pmax];

unsigned long long int   total_counter;
double   		tc0_A;
double   		tc1_A;
double			tc2_A;
double   		tc0sq_A;
double   		tc1sq_A;
double   		tc0_B;
double   		tc1_B;
double			tc2_B;
double   		tc0sq_B;
double   		tc1sq_B;
double   		tc0;
double   		tc1;
double			tc2;
double   		tc0sq;
double   		tc1sq;    
int             algorquantum;
int             algorquantum0b0;
int             algorquantum1b1;

struct triplecorrstruct {
    unsigned short  *array0;
    unsigned short  *array1;
    unsigned short	*array2;
};

struct autocorrdata {
    float *data;
    int    len;
    float *G;
};

struct readrubbish {
    float *array;
    int    length;
    char  *filename;
};


int     tuneDataBuffer[tuneDataBufferSize];
int     tuneDataBufferIndex;

//typedef unsigned long      uInt32;
typedef float   burst_array_data;
/*  Acquisition will bin each acquision into (tuneDataBufferSize/tuneDataDivisor) points; update will put [0-tuneDataBufferSize] data points up per cycle, depending on which values it updated last time.  Under severe processor load, data will get skipped. */



/*********************************************************************************
*   Try to get a continuously-updating display by using timeouts from the mainloop
*   to asynchronously access data from the aquisition routine.
**********************************************************************************/

struct  label_data 
{
    GtkLabel    *label;
    GtkLabel    *long_average; 
    GtkLabel    *ymax_label;
    GtkLabel    *ymin_label;
    GtkLabel    *ymax_label_the_second;
    GtkLabel    *ymin_label_the_second;
    double      data;
    double      data2;
    double      data3;
	float		fcs_update_data1[NUM_GRAPH_POINTS];
	float		fcs_update_data2[NUM_GRAPH_POINTS];
    float       tuneDataBuffer[tuneDataBufferSize];
    float       tuneDataBuffer2[tuneDataBufferSize];
    float       chan0[NUM_GRAPH_POINTS];
    float       chan1[NUM_GRAPH_POINTS];
    float       chan2[NUM_GRAPH_POINTS];
    double      time; 
    double      time_past;
    double      trace[NUM_GRAPH_POINTS];
    double      trace2[NUM_GRAPH_POINTS];
    double      trace3[NUM_GRAPH_POINTS];
    GdkPoint    graphthis[NUM_GRAPH_POINTS];
    GdkPoint    graphthis2[NUM_GRAPH_POINTS];
    GdkPoint    graphthis3[NUM_GRAPH_POINTS];
    GtkWidget   *widget;
    GtkWidget   *widget_the_second;
    GtkWidget   *widget1;
    GtkWidget   *widget2;
    GtkWidget   *widget3;
    GtkWidget   *widget4;
    int         graph_max_y;
    int         graph_min_y;
    float       graf_max_y;
    float       graf_min_y;
    int         graph2_max_y;
    int         graph2_min_y;
    int         write_bursts;
    int         fcs_acq_time_sec;
    int         fcs_acq_time_samples;
    gboolean    notyetplotted;
    GtkToggleButton   *button_0;
    GtkToggleButton   *button_1;
    GtkToggleButton   *button_2;
    GtkToggleButton   *button_3;
    GtkToggleButton   *button_4;
    GtkToggleButton   *button_5;
    GtkToggleButton   *button_6;
    GtkToggleButton   *button_7;
    GtkToggleButton   *button_8;
    GtkToggleButton   *button_9;
    int         test;
	float		fcs_tune_intensity[3];
	float		fcs_tune_n[3];
	float		fcs_tune_e[3];
	float		fcs_tune_stdev[3];
	GtkLabel   *fcs_stats_intensity;
	GtkLabel   *fcs_stats_n;
	GtkLabel   *fcs_stats_e;
	GtkLabel   *fcs_stats_stdev;
	GtkLabel   *fcs_stats_reserved;
};

struct  paramAcq    {
    int     runLock;
    int     run;
    int     contReadFreq[5];
    int     activeStates;
    int     clock_frequency_array_value;
    int     correlation_mode;
    int     clock_freq_scalar_value;
};

struct  paramAcq NICard;

struct  stream  {
    int *checkout;
    int *file_counter;
    FILE * fp;
    unsigned short *array;
    int len;
};

void    tim_label_set_text  (GtkLabel *label, const gchar *format, ...) G_GNUC_PRINTF (2,3);
void    launchAcq  (GtkWidget *button, struct label_data *val);
void    closeAcq   (GtkWidget *button, void *val);
void    startFCS   (GtkWidget *button, struct label_data *val);
void    haltFCS    (GtkWidget *button, void *val);
void    startG0   (GtkWidget *button, struct label_data *val);
void    haltG0    (GtkWidget *button, void *val);
void    startBURST   (GtkWidget *button, struct label_data *val);
void    haltBURST    (GtkWidget *button, void *val);
void    startBurstWrite   (GtkWidget *button, struct label_data *val);
void    haltBurstWrite    (GtkWidget *button, void *val);
void    startSTREAM  (GtkWidget *button, struct label_data *val);
//void    reset_y    (GtkWidget *button, void *val);
void    reset_yNew  (GtkWidget *button, void *ii);
void    tune_freq  (GtkWidget *button, void *ii);
void    fcs_channels   (GtkWidget *button, void *ii);
void    *Acq(void *threadarg);
void    *AcqNew(void *threadarg);
void    *AcqFCS0x1(void *threadarg);
void    *AcqFCS1x2(void *threadarg);
void    *AcqFCS2x0(void *threadarg);
void    *AcqG0_0x0(void *threadarg);
void    *AcqG0_0x1(void *threadarg);
void    *AcqG0_1x1(void *threadarg);
void    *AcqG0_2x2(void *threadarg);
void    *AcqBurstWrite(void *threadarg);
void    *AcqBURST(void *threadarg);
void    *AcqBURSTNew(void *threadarg);
void    *AcqRxn(void *threadarg);
void    *AutoCorr(void *threadarg);
void    *CrossCorr(void *threadarg);
void    *StreamIt(void *threadarg);
void    setup_graph ( GtkWidget *widget);
void    eta_entry_cb(GtkWidget *widget, GtkWidget *entry);
void    fracd1_entry_cb(GtkWidget *widget, GtkWidget *entry);
void    taud1_entry_cb(GtkWidget *widget, GtkWidget *entry);
void    taud2_entry_cb(GtkWidget *widget, GtkWidget *entry);
void    r0_entry_cb(GtkWidget *widget, GtkWidget *entry);
void    z0_entry_cb(GtkWidget *widget, GtkWidget *entry);
void    tauf_entry_cb(GtkWidget *widget, GtkWidget *entry);
void    fracf_entry_cb(GtkWidget *widget, GtkWidget *entry);
gint    cb_update(void *ii);
gint    cb_reset_y(void *ii);
gint    cb_updateNew(void *ii);
gint    cb_update_text(void *ii);
gint    cb_update_text_New(void *ii);
gint    fcs_update(void *ii);
gint    G0_update(void *ii);
gint    burst_update(void *ii);
gint    calc_fcs_fit_raw_resid (void *ii);
static  gboolean delete_event( GtkWidget *widget, GdkEvent  *event, gpointer   data );
static  gboolean expose_event( GtkWidget *widget, GdkEventExpose *event );
static  gboolean configure_event( GtkWidget *widget, GdkEventConfigure *event );
static  gboolean expose_event_fcs_G( GtkWidget *widget, GdkEventExpose *event );
static  gboolean configure_event_fcs_G( GtkWidget *widget, GdkEventConfigure *event );
static  gboolean expose_event_burst1( GtkWidget *widget, GdkEventExpose *event );
static  gboolean expose_event_burst2( GtkWidget *widget, GdkEventExpose *event );
static  gboolean expose_event_burst3( GtkWidget *widget, GdkEventExpose *event );
static  gboolean expose_event_burst4( GtkWidget *widget, GdkEventExpose *event );
static  gboolean configure_event_burst1( GtkWidget *widget, GdkEventConfigure *event );
static  gboolean configure_event_burst2( GtkWidget *widget, GdkEventConfigure *event );
static  gboolean configure_event_burst3( GtkWidget *widget, GdkEventConfigure *event );
static  gboolean configure_event_burst4( GtkWidget *widget, GdkEventConfigure *event );
static  void burst_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin);
static  void tune_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin);
static  void g_low_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin);
static  void g_high_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin);
int     imax (int a, int b);
int     imin (int a, int b);
int     write_bursts_to_file (char * filename, float * burst_array, int len_of_array);
int     write_shortbursts_to_file (char * filename, unsigned short * burst_array, int len_of_array);
int 	 bits_to_words(unsigned short *one, unsigned short *two, unsigned short *thr, int len, unsigned int *remainder, unsigned char *data);
void    *TripleCorr_A(void *threadarg);
void    *TripleCorr_B(void *threadarg);
int ShatzelDriver_A(unsigned short *data_0, unsigned short *data_1, unsigned short *data_2, int amt_data, float time_quantum, float *gav0b0, float *gav0b1, float *gav1b1, float *gav2b2, float *gstd0b0, float *gstd0b1, float *gstd1b1, float *gstd2b2, float *counts);
int ShatzelDriver_B(unsigned short *data_0, unsigned short *data_1, unsigned short *data_2, int amt_data, float time_quantum, float *gav0b0, float *gav0b1, float *gav1b1, float *gav2b2, float *gstd0b0, float *gstd0b1, float *gstd1b1, float *gstd2b2, float *counts);
int ShatzelDriver_C(unsigned short *data_0, unsigned short *data_1, unsigned short *data_2, int amt_data, float time_quantum, float *gav0b0, float *gav0b1, float *gav1b1, float *gav2b2, float *gstd0b0, float *gstd0b1, float *gstd1b1, float *gstd2b2, float *counts);

//  Pixmaps (graphs)
static  GdkPixmap *pixmap = NULL;
static  GdkPixmap *pixmap_fcs_khz = NULL;
static  GdkPixmap *pixmap_fcs_G = NULL;
static  GdkPixmap *pixmap_burst1 = NULL;
static  GdkPixmap *pixmap_burst2 = NULL;
static  GdkPixmap *pixmap_burst3 = NULL;
static  GdkPixmap *pixmap_burst4 = NULL;
int     burst_pixmap_index = 0;

//  Global Arrays, Scalars
unsigned short   *fcs_data;
unsigned short   *fcr_data;
float   *burst_data;
float   *burst_data2;
float   *burst_data3;
float   *long_burst_data;
float   *GofTau;
float	GofTau_a[(veclen)];
float	GofTau_b[(veclen)];
float	GofTau_c[(veclen)];
int     updateGofTauPlot = 0;
float   *long_average_array;
int     Adatum, Idatum;
float   burst_scale;
int     tune_points_to_average;
int     write_time;
float   fitfcs_eta = 1.0f;
float   fitfcs_r0 = 1.0f;
float   fitfcs_z0 = 1.0f;
float   fitfcs_fracd1 = 0.0f;
float   fitfcs_taud1 = 1.0f;
float   fitfcs_taud2 = 1.0f;
float   fitfcs_tauf = 0.0f;
float   fitfcs_fracf = 0.0f;
int     circ1_checkout = 0;
int     circ2_checkout = 0;
int		bursts_running = 0;
pthread_mutex_t mutex_file_key;


//  Colours
static GdkColor cyan;
static GdkGC *cyan_gc;
static GdkColor tmr;
static GdkGC *tmr_gc;
static GdkColor green;
static GdkGC *green_gc;
static GdkColor red;
static GdkGC *red_gc;
static GdkColor background_Gray;


