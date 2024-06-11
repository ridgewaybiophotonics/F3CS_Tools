/*
 *  AlignTool.c
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
 *  AlignTool.c
 *  This file contains code to generate the user interface of F3CS_AlignTool, 
 *  as well as to perform data acquisition and data processing.  
 *  
 */

#include <stdio.h>
#include <pthread.h>
#include <NIDAQmx.h>
#include <stdlib.h>
#include "gtk/gtk.h"
#include "DAQ_lib.h"
#include "complib.h"
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <gtk/gtklabel.h>

//  For black fonts, delete or comment out the line below
#define WhiteFonts 1

int main( int   argc, char *argv[])     {
    NICard.runLock = 0;
    NICard.run = 1;
    NICard.activeStates = 0;
    NICard.contReadFreq[0] = 1e5;
    NICard.contReadFreq[1] = 5e4;
    NICard.contReadFreq[2] = 2.5e4;
    NICard.contReadFreq[3] = 1.25e4;    
    NICard.contReadFreq[4] = 6250;
    NICard.clock_frequency_array_value = 1;
    NICard.correlation_mode = 0;
    int w;
    for(w=0; w<tuneDataBufferSize; w++) tuneDataBuffer[w] = 0;
    tuneDataBufferIndex = 0;
    long_average_array = (float*)malloc(LONG_AVERAGE_TIME*sizeof(float));
    for(w=0; w<LONG_AVERAGE_TIME; w++) long_average_array[w] = 0;
    GofTau   = (float*)malloc((n*(1+pmax)+1)*sizeof(float)); 

	GtkWidget *tab_label1, *tab_label2, *tab_label3;
	GtkWidget *notebook;
    GtkWidget *window;
    GtkWidget *button;
    GtkWidget *start_button;
    GtkWidget *stop_button;
    GtkWidget *reset_button;
    GtkWidget *ymax_label;
    GtkWidget *ymin_label;
    GtkWidget *readout1;
    GtkWidget *readout2;
	GtkWidget *timelabel;
	GtkWidget *intlabel;
    GtkWidget *Avglabel;
    GtkWidget *Inslabel;
    GtkWidget *vspacer;
    GtkWidget *v2spacer;
    GtkWidget *v3spacer;
    GtkWidget *hspacer;
    GtkWidget *h2spacer;
    GtkWidget *h3spacer;
    GtkWidget *quitbox;
    GtkWidget *BigLabelBox1;
    GtkWidget *BigLabelBox2;
    GtkWidget *BigLabelBox;
	GtkWidget *IntensityBox;
    GtkWidget *NIControlBox;
    GtkWidget *GraphLabelBox;
    GtkWidget *the_matrix, *fcs_matrix, *burst_matrix;
    GtkWidget *drawing_area;
    GdkColor  tune_color;
    gdk_color_parse ("red", &tune_color);
    GdkColor  gray_color;
    gdk_color_parse ("black", &gray_color);
    char      *markup;
    
    g_thread_init(NULL);
    gdk_threads_init();
    gdk_threads_enter();
    gtk_init (&argc, &argv);
  
    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    cyan.red=0x0000;
    cyan.green=0xffff;
    cyan.blue=0xffff;
    tmr.red=0xffff;
    tmr.green=0xffff;
    tmr.blue=0x0000;
    red.red=0xffff;
    red.green=0x0000;
    red.blue=0x0000f;
    green.red=0x0000;
    green.green=0xffff;
    cyan.blue=0x0000;
    background_Gray.red=0x4000;
    background_Gray.green=0x4000;
    background_Gray.blue=0x4000;

    
	gtk_widget_modify_bg(window,GTK_STATE_NORMAL, &gray_color);
    g_signal_connect (G_OBJECT (window), "delete_event", G_CALLBACK (delete_event), NULL);
    gtk_container_set_border_width (GTK_CONTAINER (window), 2);
    the_matrix = gtk_table_new ( 5, 5, FALSE);
    notebook = gtk_notebook_new ();
    gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP );
    gtk_container_add (GTK_CONTAINER (window), notebook);
    
    fcs_matrix = gtk_table_new ( 5, 5, FALSE);
    burst_matrix = gtk_table_new (5, 6, FALSE);
	
	tab_label1 = gtk_label_new ("Alignment");
    #ifdef WhiteFonts	
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 18\" foreground=\"#FFFFFF\" weight=\"bold\">Alignment</span>");
    #else 
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 18\" foreground=\"#000000\" weight=\"bold\">Alignment</span>");
    #endif
    gtk_label_set_markup (GTK_LABEL(tab_label1), markup);
	gtk_widget_show (tab_label1);
	gtk_notebook_append_page (GTK_NOTEBOOK (notebook), the_matrix, tab_label1);

    tab_label2 = gtk_label_new ("Bursts");
    #ifdef WhiteFonts
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 18\" foreground=\"#FFFFFF\" weight=\"bold\">Bursts</span>");
    #else 
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 18\" foreground=\"#000000\" weight=\"bold\">Bursts</span>");
    #endif
    gtk_label_set_markup (GTK_LABEL(tab_label2), markup);
	gtk_widget_show (tab_label2);
	gtk_notebook_append_page (GTK_NOTEBOOK (notebook), burst_matrix, tab_label2);
    gtk_notebook_set_show_tabs( GTK_NOTEBOOK(notebook), 1);
    gtk_notebook_set_show_border( GTK_NOTEBOOK(notebook), 1);

	tab_label3 = gtk_label_new ("FCS    ");
    gtk_widget_modify_bg(tab_label3,GTK_STATE_NORMAL, &tune_color);
	gtk_widget_modify_bg(notebook,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_fg(notebook,GTK_STATE_NORMAL, &background_Gray);
    #ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 18\" foreground=\"#FFFFFF\" weight=\"bold\">   FCS   </span>");
    #else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 18\" foreground=\"#000000\" weight=\"bold\">   FCS   </span>");
    #endif
    gtk_label_set_markup (GTK_LABEL(tab_label3), markup);
	gtk_widget_show (tab_label3);
	gtk_notebook_append_page (GTK_NOTEBOOK (notebook), fcs_matrix, tab_label3);
    
    
    gtk_widget_modify_bg(tab_label2,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_fg(tab_label3,GTK_STATE_NORMAL, &background_Gray);
    
	gtk_widget_show (notebook);

    vspacer  = gtk_vseparator_new ();
    v2spacer = gtk_vseparator_new ();
    v3spacer = gtk_vseparator_new ();  
    hspacer  = gtk_hseparator_new ();
    h2spacer = gtk_hseparator_new ();
    h3spacer = gtk_hseparator_new ();
    
    gtk_widget_modify_fg(vspacer,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_fg(v2spacer,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_fg(v3spacer,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_fg(hspacer,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_fg(h2spacer,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_fg(h3spacer,GTK_STATE_NORMAL, &background_Gray);
    
    gtk_widget_set_size_request (hspacer,  FRAME_WIDTH, 5);
    gtk_widget_set_size_request (h2spacer, 100, 5);
    gtk_widget_set_size_request (h3spacer, 100, 5);
    gtk_widget_set_size_request (vspacer,  5, 250);
    gtk_widget_set_size_request (v2spacer, 5, 300);
    gtk_widget_set_size_request (v3spacer, 5, 100);
    gtk_table_attach(GTK_TABLE(the_matrix), vspacer,  0,1,1,2,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(the_matrix), v2spacer, 0,1,2,3,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(the_matrix), v3spacer, 0,1,3,4,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(the_matrix), hspacer,  2,3,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(the_matrix), h2spacer, 3,4,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(the_matrix), h3spacer, 1,2,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_widget_show(vspacer);
	gtk_widget_show(v2spacer);
	gtk_widget_show(v3spacer);
	gtk_widget_show(hspacer);
	gtk_widget_show(h2spacer);
	gtk_widget_show(h3spacer);
	
    BigLabelBox = gtk_vbox_new(FALSE,5);
    BigLabelBox1 = gtk_hbox_new(FALSE,5);
    BigLabelBox2 = gtk_hbox_new(FALSE,5);
    
	readout1 = gtk_label_new ("");
	readout2 = gtk_label_new ("");
	gtk_label_set_markup (GTK_LABEL (readout1), "<small>Small text</small>");
	gtk_label_set_markup (GTK_LABEL (readout2), "<small>Small text</small>");
	tim_label_set_text (GTK_LABEL (readout1), "Level %i", 5);
#ifdef WhiteFonts
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 96\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2f Hz</span>", 0.0);
    gtk_label_set_markup (GTK_LABEL(readout1), markup);
	tim_label_set_text (GTK_LABEL (readout2), "Level %i", 5);
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 96\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2f Hz</span>", 0.0);
    gtk_label_set_markup (GTK_LABEL(readout2), markup);
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 96\" foreground=\"#000000\" weight=\"bold\">%6.2f Hz</span>", 0.0);
    gtk_label_set_markup (GTK_LABEL(readout1), markup);
	tim_label_set_text (GTK_LABEL (readout2), "Level %i", 5);
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 96\" foreground=\"#000000\" weight=\"bold\">%6.2f Hz</span>", 0.0);
    gtk_label_set_markup (GTK_LABEL(readout2), markup);
#endif
    
    
    Avglabel   = gtk_label_new ("Averaged");
    Inslabel   = gtk_label_new ("Instant");
#ifdef WhiteFonts
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">Averaged\n</span>");
    gtk_label_set_markup (GTK_LABEL(Avglabel), markup);
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">Instant\n</span>");
    gtk_label_set_markup (GTK_LABEL(Inslabel), markup);
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">Averaged\n</span>");
    gtk_label_set_markup (GTK_LABEL(Avglabel), markup);
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">Instant\n</span>");
    gtk_label_set_markup (GTK_LABEL(Inslabel), markup);
#endif
	gtk_label_set_angle(GTK_LABEL (Avglabel), 90);
    gtk_label_set_angle(GTK_LABEL (Inslabel), 90);


    gtk_box_pack_start (GTK_BOX (BigLabelBox1), Avglabel, FALSE, FALSE, 10);
	gtk_box_pack_start (GTK_BOX(BigLabelBox1), readout2, FALSE, TRUE, 0);
    
    gtk_box_pack_start (GTK_BOX (BigLabelBox2), Inslabel, FALSE, FALSE, 10);
	gtk_box_pack_start (GTK_BOX(BigLabelBox2), readout1, FALSE, TRUE, 0);
	
    gtk_box_pack_start (GTK_BOX (BigLabelBox), BigLabelBox1, FALSE, FALSE, 10);
    gtk_box_pack_start (GTK_BOX (BigLabelBox), BigLabelBox2, FALSE, FALSE, 10);
    
    //	L,R, Top, Bottom
    gtk_table_attach(GTK_TABLE(the_matrix), BigLabelBox,1,3,1,2,GTK_SHRINK,GTK_FILL, 4,4);

    gtk_widget_show (BigLabelBox);
    gtk_widget_show (BigLabelBox1);
    gtk_widget_show (BigLabelBox2);
	gtk_widget_show (Avglabel);
    gtk_widget_show (Inslabel);
    gtk_widget_show (readout1);
	gtk_widget_show (readout2);

	// Establish the intensity graph
    drawing_area = gtk_drawing_area_new();
    gtk_widget_set_size_request(GTK_WIDGET(drawing_area),FRAME_WIDTH,300);
	//Note:  the above line sets the # of pixels in the drawing area;
	
	timelabel = gtk_label_new("Reads to Avg.");
    #ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">  Time  \n</span>");
    #else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">  Time  \n</span>");
    #endif
    gtk_label_set_markup (GTK_LABEL(timelabel), markup);
	
	IntensityBox = gtk_vbox_new(FALSE,5);
	gtk_box_pack_start (GTK_BOX(IntensityBox), drawing_area, FALSE, TRUE, 0);
	gtk_box_pack_end (GTK_BOX(IntensityBox), timelabel, FALSE, TRUE, 0);
	gtk_table_attach(GTK_TABLE(the_matrix), IntensityBox, 2,3,2,3,GTK_SHRINK,GTK_FILL, 4, 4);

    GraphLabelBox = gtk_vbox_new(TRUE,0);
    ymax_label = gtk_label_new ("Y Max");
#ifdef WhiteFonts
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 10\" foreground=\"#FFFFFF\" weight=\"bold\">  Ymax  \n</span>");
    gtk_label_set_markup (GTK_LABEL(ymax_label), markup);
    ymin_label = gtk_label_new ("Y Min");
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 10\" foreground=\"#FFFFFF\" weight=\"bold\">  Ymin  \n</span>");
    gtk_label_set_markup (GTK_LABEL(ymin_label), markup);
	intlabel   = gtk_label_new ("Intensity");
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">  Intensity  \n</span>");
    gtk_label_set_markup (GTK_LABEL(intlabel), markup);
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">  Ymax  \n</span>");
    gtk_label_set_markup (GTK_LABEL(ymax_label), markup);
    ymin_label = gtk_label_new ("Y Min");
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">  Ymin  \n</span>");
    gtk_label_set_markup (GTK_LABEL(ymin_label), markup);
	intlabel   = gtk_label_new ("Intensity");
	markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">  Intensity  \n</span>");
    gtk_label_set_markup (GTK_LABEL(intlabel), markup);
#endif
	gtk_label_set_angle(GTK_LABEL (intlabel), 90);
	gtk_box_pack_start (GTK_BOX (GraphLabelBox), ymax_label, FALSE, FALSE, 10);
	gtk_box_pack_start (GTK_BOX (GraphLabelBox), intlabel, FALSE, FALSE, 10);
	gtk_box_pack_start (GTK_BOX (GraphLabelBox), ymin_label, FALSE, FALSE, 10);
    gtk_table_attach(GTK_TABLE(the_matrix), GraphLabelBox, 1,2,2,3,GTK_SHRINK,GTK_FILL, 4, 4);
    gtk_widget_show(ymax_label);
    gtk_widget_show(ymin_label);
	gtk_widget_show(intlabel);
    gtk_widget_show(GraphLabelBox);
    gtk_widget_show(drawing_area);
	gtk_widget_show(timelabel);
	gtk_widget_show(IntensityBox);
    g_signal_connect (G_OBJECT (drawing_area), "expose_event", G_CALLBACK (expose_event), NULL);
    g_signal_connect (G_OBJECT (drawing_area),"configure_event", G_CALLBACK (configure_event), NULL);
    
    quitbox = gtk_vbox_new (FALSE, 2);
    button = gtk_button_new_with_label ("  Exit  ");
    g_signal_connect_swapped(G_OBJECT(button),"clicked",G_CALLBACK(delete_event), G_OBJECT(window));
    gtk_box_pack_end (GTK_BOX (quitbox), button, FALSE, FALSE, 10);
	gtk_table_attach(GTK_TABLE(the_matrix), quitbox, 3,4,3,5,GTK_SHRINK,GTK_FILL, 4,4);
    gtk_widget_show (button);
    gtk_widget_show (quitbox);

    NIControlBox = gtk_vbox_new (FALSE, 2);
    start_button = gtk_button_new_with_label ("Start");
    stop_button = gtk_button_new_with_label ("Stop");
    reset_button = gtk_button_new_with_label ("Reset Y-Axis");

    GtkWidget   *tune_freq_frame;
    tune_freq_frame = gtk_frame_new (NULL);

    GtkWidget   *tune_freq_box;
    tune_freq_box = gtk_vbox_new(FALSE, 7);
    gtk_widget_modify_fg(tune_freq_frame,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_bg(tune_freq_frame,GTK_STATE_NORMAL, &background_Gray);
    gtk_table_attach(GTK_TABLE(the_matrix), tune_freq_frame, 3, 4, 2, 3, GTK_SHRINK, GTK_FILL, 4, 4);
    gtk_frame_set_label (GTK_FRAME(tune_freq_frame), "");
    gtk_frame_set_label_align (GTK_FRAME(tune_freq_frame), 0.0, 1.0);
    gtk_container_add (GTK_CONTAINER (tune_freq_frame), tune_freq_box);

	//  Create a way to allow the user to determine how many data points are averaged into the top display
    tune_points_to_average = tune_points_to_average_initial_value;
    GtkWidget       *tune_scale_button;
    GtkWidget       *tune_scale_button_box;
    GtkAdjustment   *tune_scale_button_adj;
    tune_scale_button_box = gtk_vbox_new (FALSE, 0);
    GtkWidget       *tune_scale_button_label;
    tune_scale_button_label = gtk_label_new("Reads to Avg.");
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">Reads to Avg.\n</span>");
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">Reads to Avg.\n</span>");
#endif
    gtk_label_set_markup (GTK_LABEL(tune_scale_button_label), markup);
    gtk_box_pack_start (GTK_BOX(tune_scale_button_box), tune_scale_button_label, FALSE, TRUE, 0);
    tune_scale_button_adj = (GtkAdjustment *) gtk_adjustment_new (tune_points_to_average_initial_value, 1, NUM_GRAPH_POINTS, 1.0, 10.0, 0.0);
    tune_scale_button = gtk_spin_button_new (tune_scale_button_adj, tune_points_to_average_initial_value, 2);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (tune_scale_button), FALSE);
    gtk_widget_set_size_request (tune_scale_button, 100, -1);
    gtk_box_pack_start (GTK_BOX(tune_scale_button_box), tune_scale_button, FALSE, TRUE, 0);
    g_signal_connect (G_OBJECT(tune_scale_button_adj), "value_changed", G_CALLBACK(tune_scale_button_callback), (gpointer) tune_scale_button);
    gtk_box_pack_start(GTK_BOX(tune_freq_box), tune_scale_button_box, FALSE, FALSE, 10);
    gtk_widget_show (tune_scale_button_label);
    gtk_widget_show (tune_scale_button);
    gtk_widget_show (tune_scale_button_box);
    
    gtk_widget_show(tune_freq_frame);
    gtk_widget_show(tune_freq_box);

    struct label_data ilab;
    struct label_data jlab;
    ilab.label = GTK_LABEL(readout1);
    ilab.long_average = GTK_LABEL(readout2);
    ilab.ymax_label = GTK_LABEL(ymax_label);
    ilab.ymin_label = GTK_LABEL(ymin_label);
    ilab.data = 0;
    for(w=0; w< NUM_GRAPH_POINTS; w++) ilab.trace[w] = 0;
    ilab.time = 0;
    ilab.notyetplotted = TRUE;
    ilab.widget = GTK_WIDGET(drawing_area);
    ilab.graph_max_y = 100;
    jlab.label = GTK_LABEL(readout2);
    jlab.data = 0;

    g_signal_connect(G_OBJECT(start_button),"clicked",G_CALLBACK(launchAcq), &ilab);
    g_signal_connect(G_OBJECT(stop_button),"clicked",G_CALLBACK(closeAcq), &jlab);
    g_signal_connect(G_OBJECT(reset_button),"clicked",G_CALLBACK(reset_yNew), &ilab);
    gtk_box_pack_start (GTK_BOX (NIControlBox), start_button, FALSE, FALSE, 10);
    gtk_box_pack_start (GTK_BOX (NIControlBox), stop_button, FALSE, FALSE, 10);
    gtk_box_pack_start (GTK_BOX (NIControlBox), reset_button, FALSE, FALSE, 10);
	gtk_table_attach(GTK_TABLE(the_matrix), NIControlBox, 3, 4, 1, 2,GTK_SHRINK,GTK_FILL, 4, 4);    
    
    gtk_timeout_add(32, cb_updateNew, &ilab);
    gtk_timeout_add(200, cb_update_text_New, &ilab);
    gtk_timeout_add(8000, cb_reset_y, &ilab);

    gtk_widget_show(start_button);
    gtk_widget_show(stop_button);
    gtk_widget_show(reset_button);
    gtk_widget_show(NIControlBox);
    
	/************Configure Burst Interface**********************************/
    burst_scale = 1;
    GtkWidget *burst_v10spacer;
    GtkWidget *burst_v12spacer;
    GtkWidget *burst_v13spacer;
    GtkWidget *burst_v14spacer;
    GtkWidget *burst_h10spacer;
    GtkWidget *burst_h12spacer;
    GtkWidget *burst_quitbox;
    GtkWidget *burst_quitbut;
    
    burst_v10spacer  = gtk_vseparator_new ();
    burst_v12spacer = gtk_vseparator_new ();
    burst_v13spacer = gtk_vseparator_new (); 
    burst_v14spacer = gtk_vseparator_new ();
    burst_h10spacer  = gtk_hseparator_new ();
    burst_h12spacer = gtk_hseparator_new ();
    gtk_widget_set_size_request (burst_h10spacer,  BURST_FRAME_WIDTH, 5);
    gtk_widget_set_size_request (burst_h12spacer, 100, 5);
    gtk_widget_set_size_request (burst_v10spacer,  5, 140);
    gtk_widget_set_size_request (burst_v12spacer, 5, 150);
    gtk_widget_set_size_request (burst_v13spacer, 5, 150);
    gtk_widget_set_size_request (burst_v14spacer, 5, 150);
    gtk_table_attach(GTK_TABLE(burst_matrix), burst_v10spacer,  0,1,1,2,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(burst_matrix), burst_v12spacer, 0,1,2,3,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(burst_matrix), burst_v13spacer, 0,1,3,4,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(burst_matrix), burst_v14spacer, 0,1,4,5,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(burst_matrix), burst_h10spacer,  1,3,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(burst_matrix), burst_h12spacer, 3,4,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_widget_show(burst_v10spacer);
	gtk_widget_show(burst_v12spacer);
	gtk_widget_show(burst_v13spacer);
	gtk_widget_show(burst_v14spacer);
	gtk_widget_show(burst_h10spacer);
	gtk_widget_show(burst_h12spacer);
	
    burst_quitbox = gtk_vbox_new (FALSE, 2);
    burst_quitbut = gtk_button_new_with_label ("  Exit  ");
    g_signal_connect_swapped(G_OBJECT(burst_quitbut),"clicked",G_CALLBACK(delete_event), G_OBJECT(window));
    gtk_box_pack_end (GTK_BOX (burst_quitbox), burst_quitbut, FALSE, FALSE, 10);
	gtk_table_attach(GTK_TABLE(burst_matrix), burst_quitbox, 3,4,5,6,GTK_SHRINK,GTK_FILL, 4,4);
    gtk_widget_show (burst_quitbut);
    gtk_widget_show (burst_quitbox);
    
    GtkWidget   *burst_plot1;
    GtkWidget   *burst_plot2;
    GtkWidget   *burst_plot3;
    GtkWidget   *burst_plot4;
    burst_plot1 = gtk_drawing_area_new();
    gtk_widget_set_size_request(GTK_WIDGET(burst_plot1),BURST_FRAME_WIDTH,150);
    gtk_table_attach(GTK_TABLE(burst_matrix), burst_plot1, 1,3,1,2, GTK_SHRINK,GTK_FILL, 4, 4);
    burst_plot2 = gtk_drawing_area_new();
    gtk_widget_set_size_request(GTK_WIDGET(burst_plot2),BURST_FRAME_WIDTH,150);
    gtk_table_attach(GTK_TABLE(burst_matrix), burst_plot2, 1,3,2,3, GTK_SHRINK,GTK_FILL, 4, 4);
    burst_plot3 = gtk_drawing_area_new();
    gtk_widget_set_size_request(GTK_WIDGET(burst_plot3),BURST_FRAME_WIDTH,150);
    gtk_table_attach(GTK_TABLE(burst_matrix), burst_plot3, 1,3,3,4, GTK_SHRINK,GTK_FILL, 4, 4);
    burst_plot4 = gtk_drawing_area_new();
    gtk_widget_set_size_request(GTK_WIDGET(burst_plot4),BURST_FRAME_WIDTH,150);
    gtk_table_attach(GTK_TABLE(burst_matrix), burst_plot4, 1,3,4,5, GTK_SHRINK,GTK_FILL, 4, 4);
    
    gtk_widget_show(burst_plot1);
    gtk_widget_show(burst_plot2);
    gtk_widget_show(burst_plot3);
    gtk_widget_show(burst_plot4);
    g_signal_connect (G_OBJECT (burst_plot1), "expose_event", G_CALLBACK (expose_event_burst1), NULL);
    g_signal_connect (G_OBJECT (burst_plot1),"configure_event", G_CALLBACK (configure_event_burst1), NULL);
    g_signal_connect (G_OBJECT (burst_plot2), "expose_event", G_CALLBACK (expose_event_burst2), NULL);
    g_signal_connect (G_OBJECT (burst_plot2),"configure_event", G_CALLBACK (configure_event_burst2), NULL);
    g_signal_connect (G_OBJECT (burst_plot3), "expose_event", G_CALLBACK (expose_event_burst3), NULL);
    g_signal_connect (G_OBJECT (burst_plot3),"configure_event", G_CALLBACK (configure_event_burst3), NULL);
    g_signal_connect (G_OBJECT (burst_plot4), "expose_event", G_CALLBACK (expose_event_burst4), NULL);
    g_signal_connect (G_OBJECT (burst_plot4),"configure_event", G_CALLBACK (configure_event_burst4), NULL);
    

    GSList      *burst_freq_group;
    GtkWidget   *burst_freq_box;
    GtkWidget   *burst_freq_box_box;
    GtkWidget   *burst_freq_but1;   
    GtkWidget   *burst_freq_but2;
    GtkWidget   *burst_freq_but3;
    GtkWidget   *burst_freq_but4;
    GtkWidget   *burst_freq_but5;
    GtkWidget   *burst_freq_but6;
    GtkWidget   *burst_freq_but7;
    GtkWidget   *burst_freq_but8;
    GtkWidget   *burst_freq_but9;
    GtkWidget   *burst_freq_but10;
    GtkWidget   *burst_freq_frame;
    GtkWidget   *bin_width_label;
    burst_freq_frame = gtk_frame_new (NULL);
    gtk_widget_modify_fg(burst_freq_frame,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_bg(burst_freq_frame,GTK_STATE_NORMAL, &background_Gray);
    gtk_table_attach(GTK_TABLE(burst_matrix), burst_freq_frame, 3, 4, 3, 5, GTK_SHRINK, GTK_FILL, 4, 4);
    gtk_frame_set_label (GTK_FRAME(burst_freq_frame), ""); // Bin Width
    bin_width_label = gtk_label_new("Bin Width");
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">  Bin Width  \n</span>");
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">  Bin Width  \n</span>");
#endif
    gtk_label_set_markup (GTK_LABEL(bin_width_label), markup);
    gtk_frame_set_label_align (GTK_FRAME(burst_freq_frame), 0.0, 1.0);
    burst_freq_box_box = gtk_vbox_new(FALSE, 7);
    burst_freq_box = gtk_vbox_new(FALSE, 7);

    gtk_box_pack_start(GTK_BOX(burst_freq_box_box), bin_width_label, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(burst_freq_box_box), burst_freq_box, FALSE, FALSE, 0);
    gtk_container_add (GTK_CONTAINER (burst_freq_frame), burst_freq_box_box);
    burst_freq_but1 = gtk_radio_button_new_with_label (NULL, "50ns  ");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but1, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but1);
    gtk_widget_show(bin_width_label);
    
    burst_freq_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (burst_freq_but1));
    burst_freq_but2 = gtk_radio_button_new_with_label (burst_freq_group, "800ns");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but2, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but2);
    burst_freq_but3 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but2), "1.6us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but3, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but3);
    burst_freq_but4 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but3), "3.2us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but4, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but4);
    burst_freq_but5 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but4), "6.4us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but5, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but5);
    burst_freq_but6 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but5), "12.8us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but6, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but6);
    burst_freq_but7 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but6), "25.6us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but7, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but7);
    burst_freq_but8 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but7), "51.2us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but8, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but8);
    burst_freq_but9 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but8), "102.4us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but9, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but9);
    burst_freq_but10 = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(burst_freq_but8), "512us");
    gtk_box_pack_start(GTK_BOX(burst_freq_box), burst_freq_but10, FALSE, FALSE, 0);
    gtk_widget_show(burst_freq_but10);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (burst_freq_but2), TRUE);

    struct label_data blab;
    blab.label = GTK_LABEL(readout1);
    blab.ymax_label = GTK_LABEL(ymax_label);
    blab.ymin_label = GTK_LABEL(ymin_label);
    blab.data = 0;
    blab.time = 0;
    blab.notyetplotted = TRUE;
    blab.widget  = GTK_WIDGET(burst_plot1);
    blab.widget1 = GTK_WIDGET(burst_plot1);
    blab.widget2 = GTK_WIDGET(burst_plot2);
    blab.widget3 = GTK_WIDGET(burst_plot3);
    blab.widget4 = GTK_WIDGET(burst_plot4);
    blab.graph_max_y = 100;
    blab.button_0 = GTK_TOGGLE_BUTTON(burst_freq_but1);
    blab.button_1 = GTK_TOGGLE_BUTTON(burst_freq_but2);
    blab.button_2 = GTK_TOGGLE_BUTTON(burst_freq_but3);
    blab.button_3 = GTK_TOGGLE_BUTTON(burst_freq_but4);
    blab.button_4 = GTK_TOGGLE_BUTTON(burst_freq_but5);
    blab.button_5 = GTK_TOGGLE_BUTTON(burst_freq_but6);
    blab.button_6 = GTK_TOGGLE_BUTTON(burst_freq_but7);
    blab.button_7 = GTK_TOGGLE_BUTTON(burst_freq_but8);
    blab.button_8 = GTK_TOGGLE_BUTTON(burst_freq_but9);
    blab.button_9 = GTK_TOGGLE_BUTTON(burst_freq_but10);
    
    GtkWidget       *burst_scale_button;
    GtkWidget       *burst_scale_button_box;
    GtkAdjustment   *burst_scale_button_adj;
    burst_scale_button_box = gtk_vbox_new (FALSE, 0);
    GtkWidget       *burst_scale_button_label;
    burst_scale_button_label = gtk_label_new("Scale Value");
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">Scale Value</span>");
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">Scale Value</span>");
#endif
    gtk_label_set_markup (GTK_LABEL(burst_scale_button_label), markup);
    burst_scale_button_adj = (GtkAdjustment *) gtk_adjustment_new (1.0, 0, 1000000, 0.5, 100.0, 0.0);
    burst_scale_button = gtk_spin_button_new (burst_scale_button_adj, burst_scale, 2);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (burst_scale_button), FALSE);
    gtk_widget_set_size_request (burst_scale_button, 100, -1);
    g_signal_connect (G_OBJECT(burst_scale_button_adj), "value_changed", G_CALLBACK(burst_scale_button_callback), (gpointer) burst_scale_button);
    gtk_box_pack_end (GTK_BOX(burst_scale_button_box), burst_scale_button, FALSE, TRUE, 0);
    gtk_box_pack_end (GTK_BOX(burst_scale_button_box), burst_scale_button_label, FALSE, TRUE, 0);
	gtk_table_attach(GTK_TABLE(burst_matrix), burst_scale_button_box, 3, 4, 2, 3,GTK_SHRINK,GTK_FILL, 4, 4);
    gtk_widget_show (burst_scale_button_label);
    gtk_widget_show (burst_scale_button);

	//  Allow the user to change the total data accumulation time
    tune_points_to_average = tune_points_to_average_initial_value;
    gtk_widget_show (burst_scale_button_box);

    GtkWidget   *burstControlBox;
    burstControlBox = gtk_vbox_new (FALSE, 3);
    GtkWidget   *burst_start;
    burst_start = gtk_button_new_with_label ("Watch Bursts");
    GtkWidget   *burst_stop;
    burst_stop = gtk_button_new_with_label  ("     Stop     ");
    g_signal_connect(G_OBJECT(burst_freq_but1),"clicked",G_CALLBACK(tune_freq), &blab);
    g_signal_connect(G_OBJECT(burst_freq_but2),"clicked",G_CALLBACK(tune_freq), &blab);
    g_signal_connect(G_OBJECT(burst_freq_but3),"clicked",G_CALLBACK(tune_freq), &blab);
    g_signal_connect(G_OBJECT(burst_freq_but4),"clicked",G_CALLBACK(tune_freq), &blab);
    g_signal_connect(G_OBJECT(burst_freq_but5),"clicked",G_CALLBACK(tune_freq), &blab); 
    g_signal_connect(G_OBJECT(burst_freq_but6),"clicked",G_CALLBACK(tune_freq), &blab);
    g_signal_connect(G_OBJECT(burst_freq_but7),"clicked",G_CALLBACK(tune_freq), &blab);
    g_signal_connect(G_OBJECT(burst_freq_but8),"clicked",G_CALLBACK(tune_freq), &blab);
    g_signal_connect(G_OBJECT(burst_freq_but9),"clicked",G_CALLBACK(tune_freq), &blab);   
    g_signal_connect(G_OBJECT(burst_freq_but10),"clicked",G_CALLBACK(tune_freq), &blab);    
    g_signal_connect(G_OBJECT(burst_start),"clicked",G_CALLBACK(startBURST), &blab);
    g_signal_connect(G_OBJECT(burst_stop),"clicked",G_CALLBACK(haltBURST), &blab);
    gtk_box_pack_start (GTK_BOX (burstControlBox), burst_start, FALSE, FALSE, 5);
    gtk_box_pack_start (GTK_BOX (burstControlBox), burst_stop, FALSE, FALSE, 5);
	gtk_table_attach(GTK_TABLE(burst_matrix), burstControlBox, 3, 4, 1, 2,GTK_SHRINK,GTK_FILL, 4, 4);    
    gtk_widget_show(burst_start);
    gtk_widget_show(burst_stop);
    gtk_widget_show(burstControlBox);
    gtk_widget_show(burst_freq_frame);
    gtk_widget_show(burst_freq_box);
    gtk_widget_show(burst_freq_box_box);
    burst_data = (float*)malloc(TUNE_BUFFER*sizeof(float));
    burst_data2 = (float*)malloc(TUNE_BUFFER*sizeof(float));
    burst_data3 = (float*)malloc(TUNE_BUFFER*sizeof(float));

    gtk_timeout_add(1000, burst_update, &blab);

	/**************Setup FCS Window********************/
    GtkWidget *v10spacer;
    GtkWidget *v12spacer;
    GtkWidget *v13spacer;
    GtkWidget *h10spacer;
    GtkWidget *h12spacer;
    GtkWidget *h13spacer;
    GtkWidget *fcs_quitbox;
    GtkWidget *fcs_quitbut;
    
    v10spacer  = gtk_vseparator_new ();
    v12spacer = gtk_vseparator_new ();
    v13spacer = gtk_vseparator_new ();  
    h10spacer  = gtk_hseparator_new ();
    h12spacer = gtk_hseparator_new ();
    h13spacer = gtk_hseparator_new ();
    gtk_widget_set_size_request (h10spacer,  FRAME_WIDTH, 5);
    gtk_widget_set_size_request (h12spacer, 100, 5);
    gtk_widget_set_size_request (h13spacer, 100, 5);
    gtk_widget_set_size_request (v10spacer,  5, 20);
    gtk_widget_set_size_request (v12spacer, 5, 450);
    gtk_widget_set_size_request (v13spacer, 5, 10);
    gtk_table_attach(GTK_TABLE(fcs_matrix), v10spacer,  0,1,1,2,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(fcs_matrix), v12spacer, 0,1,2,3,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(fcs_matrix), v13spacer, 0,1,3,4,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(fcs_matrix), h10spacer,  2,3,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(fcs_matrix), h12spacer, 3,4,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_table_attach(GTK_TABLE(fcs_matrix), h13spacer, 1,2,0,1,GTK_SHRINK,GTK_FILL, 4,4);
	gtk_widget_show(v10spacer);
	gtk_widget_show(v12spacer);
	gtk_widget_show(v13spacer);
	gtk_widget_show(h10spacer);
	gtk_widget_show(h12spacer);
	gtk_widget_show(h13spacer);
	
    fcs_quitbox = gtk_vbox_new (FALSE, 2);
    fcs_quitbut = gtk_button_new_with_label ("  Exit  ");
    g_signal_connect_swapped(G_OBJECT(fcs_quitbut),"clicked",G_CALLBACK(delete_event), G_OBJECT(window));
    gtk_box_pack_end (GTK_BOX (fcs_quitbox), fcs_quitbut, FALSE, FALSE, 10);
	gtk_table_attach(GTK_TABLE(fcs_matrix), fcs_quitbox, 3,4,3,5,GTK_SHRINK,GTK_FILL, 4,4);
    gtk_widget_show (fcs_quitbut);
    gtk_widget_show (fcs_quitbox);
    
    GtkWidget   *fcs_plot;
    GtkWidget   *fcs_plot2;
    fcs_plot = gtk_vbox_new(FALSE,2);
    GtkWidget   *fcs_stats_intensity;
    fcs_stats_intensity = gtk_label_new ("Intensity");
    gtk_widget_show(fcs_stats_intensity);
    gtk_box_pack_start (GTK_BOX(fcs_plot), fcs_stats_intensity, FALSE, FALSE, 10);
    GtkWidget    *four_box;
    four_box = gtk_hbox_new(FALSE,2);
    gtk_box_pack_start (GTK_BOX(fcs_plot), four_box, FALSE, FALSE, 10);
    GtkWidget   *n_e_box;
    GtkWidget   *stdev_box;
    n_e_box = gtk_vbox_new(FALSE,3);
    stdev_box = gtk_vbox_new(FALSE,3);
    gtk_box_pack_start (GTK_BOX(four_box), n_e_box, FALSE, FALSE, 10);
    gtk_box_pack_start (GTK_BOX(four_box), stdev_box, FALSE, FALSE, 10);
    
    GtkWidget   *fcs_stats_n;
    fcs_stats_n = gtk_label_new ("fcs_stats_n");
    gtk_widget_show(fcs_stats_n);
    gtk_box_pack_start (GTK_BOX(n_e_box), fcs_stats_n, FALSE, FALSE, 10);
    
    GtkWidget   *fcs_stats_e;
    fcs_stats_e = gtk_label_new ("fcs_stats_e");
    gtk_widget_show(fcs_stats_e);
    gtk_box_pack_start (GTK_BOX(n_e_box), fcs_stats_e, FALSE, FALSE, 10);
    
    GtkWidget   *fcs_stats_stdev;
    fcs_stats_stdev = gtk_label_new ("fcs_stats_stdev1");
    gtk_widget_show(fcs_stats_stdev);
    gtk_box_pack_start (GTK_BOX(stdev_box), fcs_stats_stdev, FALSE, FALSE, 10);
    
    GtkWidget   *fcs_stats_reserved;
    fcs_stats_reserved = gtk_label_new ("");
    gtk_widget_show(fcs_stats_reserved);
    gtk_box_pack_start (GTK_BOX(stdev_box), fcs_stats_reserved, FALSE, FALSE, 10);

    gtk_widget_show (n_e_box);
    gtk_widget_show (stdev_box);
    gtk_widget_show (four_box);
    
	GtkWidget * Gtaulabel;
	GtkWidget * GtauLabelhBox;
	Gtaulabel = gtk_label_new("Gtau");
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">G(tau)\n</span>");
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">G(tau)\n</span>");
#endif
    gtk_label_set_markup (GTK_LABEL(Gtaulabel), markup);
	gtk_label_set_angle(GTK_LABEL (Gtaulabel), 90);
	
    gtk_table_attach(GTK_TABLE(fcs_matrix), fcs_plot, 1,3,1,2, GTK_SHRINK,GTK_FILL, 4, 4);
    fcs_plot2 = gtk_drawing_area_new();
    gtk_widget_set_size_request(GTK_WIDGET(fcs_plot2),FRAME_WIDTH,300);
	
	GtauLabelhBox = gtk_hbox_new(FALSE,5);

	gtk_box_pack_start (GTK_BOX (GtauLabelhBox), Gtaulabel, FALSE, FALSE, 10);
	gtk_box_pack_start (GTK_BOX (GtauLabelhBox), fcs_plot2, FALSE, FALSE, 10);
	gtk_table_attach(GTK_TABLE(fcs_matrix), GtauLabelhBox, 1,3,2,4, GTK_SHRINK,GTK_FILL,4, 4);
	
	gtk_widget_show(Gtaulabel);
	gtk_widget_show(GtauLabelhBox);
    
    gtk_widget_show(fcs_plot);
    gtk_widget_show(fcs_plot2);
    g_signal_connect (G_OBJECT (fcs_plot2), "expose_event", G_CALLBACK (expose_event_fcs_G), NULL);
    g_signal_connect (G_OBJECT (fcs_plot2),"configure_event", G_CALLBACK (configure_event_fcs_G), NULL);
    
    GtkWidget   *fcs_freq_box;
    GtkWidget   *fcs_0x1_button;   
    GtkWidget   *fcs_1x2_button;
    GtkWidget   *fcs_2x0_button;
    GtkWidget   *fcs_freq_frame;
    GtkWidget   *corr_mode_label;
    fcs_freq_frame = gtk_frame_new (NULL);
    gtk_widget_modify_fg(fcs_freq_frame,GTK_STATE_NORMAL, &background_Gray);
    gtk_widget_modify_bg(fcs_freq_frame,GTK_STATE_NORMAL, &background_Gray);
    gtk_table_attach(GTK_TABLE(fcs_matrix), fcs_freq_frame, 3, 4, 2, 3, GTK_SHRINK, GTK_FILL, 4, 4);
//    gtk_frame_set_label (GTK_FRAME(fcs_freq_frame), "Cor.Mode");
    gtk_frame_set_label_align (GTK_FRAME(fcs_freq_frame), 0.0, 1.0);
    fcs_freq_box = gtk_vbox_new(FALSE, 7);
    gtk_container_add (GTK_CONTAINER (fcs_freq_frame), fcs_freq_box);
    
    corr_mode_label = gtk_label_new("Correlation Mode");
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">  Correlation \n  Mode\n</span>");
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">  Correlation \n  Mode\n</span>");
#endif
    gtk_label_set_markup (GTK_LABEL(corr_mode_label), markup);
    gtk_box_pack_start(GTK_BOX(fcs_freq_box), corr_mode_label, FALSE, FALSE, 0);
    
    fcs_0x1_button = gtk_radio_button_new_with_label (NULL, " 0, 1 ");
    gtk_box_pack_start(GTK_BOX(fcs_freq_box), fcs_0x1_button, FALSE, FALSE, 0);
    gtk_widget_show(fcs_0x1_button);
    
    fcs_1x2_button = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(fcs_0x1_button), " 1, 2" );
    gtk_box_pack_start(GTK_BOX(fcs_freq_box), fcs_1x2_button, FALSE, FALSE, 0);
    gtk_widget_show(fcs_1x2_button);
    fcs_2x0_button = gtk_radio_button_new_with_label_from_widget (GTK_RADIO_BUTTON(fcs_1x2_button), " 2, 0 ");
    gtk_box_pack_start(GTK_BOX(fcs_freq_box), fcs_2x0_button, FALSE, FALSE, 0);
    gtk_widget_show(fcs_2x0_button);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (fcs_0x1_button), TRUE);

    struct label_data flab;
    flab.label = GTK_LABEL(readout1);
    flab.data = 0;
    flab.time = 0;
    flab.notyetplotted = TRUE;
    flab.widget = GTK_WIDGET(fcs_plot);
    flab.widget_the_second = GTK_WIDGET(fcs_plot2);
    flab.graph_max_y = 100;
    flab.button_0 = GTK_TOGGLE_BUTTON(fcs_0x1_button);
    flab.button_1 = GTK_TOGGLE_BUTTON(fcs_1x2_button);
    flab.button_2 = GTK_TOGGLE_BUTTON(fcs_2x0_button);
	flab.fcs_stats_intensity = GTK_LABEL(fcs_stats_intensity);
	flab.fcs_stats_n = GTK_LABEL(fcs_stats_n);
    flab.fcs_stats_e = GTK_LABEL(fcs_stats_e);
    flab.fcs_stats_stdev = GTK_LABEL(fcs_stats_stdev);
    flab.fcs_stats_reserved = GTK_LABEL(fcs_stats_reserved);
    flab.fcs_tune_intensity[0] = 0.0;
    flab.fcs_tune_intensity[1] = 0.0;
    flab.fcs_tune_intensity[2] = 0.0;
    flab.fcs_tune_n[0] = 0.0;
    flab.fcs_tune_n[1] = 0.0;
    flab.fcs_tune_n[2] = 0.0;
    flab.fcs_tune_e[0] = 0.0;
    flab.fcs_tune_e[1] = 0.0;
    flab.fcs_tune_e[2] = 0.0;
    flab.fcs_tune_stdev[0] = 0.0;
    flab.fcs_tune_stdev[1] = 0.0;
    flab.fcs_tune_stdev[2] = 0.0;

    GtkWidget       *g_low_scale_button;
    GtkWidget       *g_low_scale_button_label;
    GtkWidget       *g_high_scale_button_label;
    GtkWidget       *g_low_scale_button_box;
    GtkAdjustment   *g_low_scale_button_adj;
    g_low_scale_button_box = gtk_vbox_new (FALSE, 0);
    
    g_low_scale_button_label = gtk_label_new("Gtau");
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">G(tau) Min</span>");
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">G(tau) Min</span>");
#endif
    gtk_label_set_markup (GTK_LABEL(g_low_scale_button_label), markup);
    g_high_scale_button_label = gtk_label_new("Gtau");
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">G(tau) Max</span>");
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">G(tau) Max</span>");
#endif
    gtk_label_set_markup (GTK_LABEL(g_high_scale_button_label), markup);
    
//    g_low_scale_button_label = gtk_label_new("G(tau)_min");
    g_low_scale_button_adj = (GtkAdjustment *) gtk_adjustment_new (gtau_min, -10, 10, 0.001, 10.0, 0.0);
    g_low_scale_button = gtk_spin_button_new (g_low_scale_button_adj, burst_scale, 3);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (g_low_scale_button), FALSE);
    gtk_widget_set_size_request (g_low_scale_button, 100, -1);
    g_signal_connect (G_OBJECT(g_low_scale_button_adj), "value_changed", G_CALLBACK(g_low_scale_button_callback), (gpointer) g_low_scale_button);
    gtk_box_pack_end (GTK_BOX(g_low_scale_button_box), g_low_scale_button, FALSE, TRUE, 0);
    
    gtk_box_pack_end (GTK_BOX(g_low_scale_button_box), g_low_scale_button_label, FALSE, TRUE, 0);
    gtk_widget_show (corr_mode_label);
    gtk_widget_show (g_low_scale_button);
    gtk_widget_show (g_high_scale_button_label);
    gtk_widget_show (g_low_scale_button_label);

    GtkWidget       *g_high_scale_button;
    GtkWidget       *g_high_scale_button_box;
    GtkAdjustment   *g_high_scale_button_adj;
    g_high_scale_button_box = gtk_vbox_new (FALSE, 0);
//    g_high_scale_button_label = gtk_label_new("G(tau)_min");
    g_high_scale_button_adj = (GtkAdjustment *) gtk_adjustment_new (gtau_max, 0, 100, 0.001, 10.0, 0.0);
    g_high_scale_button = gtk_spin_button_new (g_high_scale_button_adj, burst_scale, 3);
    gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (g_high_scale_button), FALSE);
    gtk_widget_set_size_request (g_high_scale_button, 100, -1);
    g_signal_connect (G_OBJECT(g_high_scale_button_adj), "value_changed", G_CALLBACK(g_high_scale_button_callback), (gpointer) g_high_scale_button);
    gtk_box_pack_end (GTK_BOX(g_high_scale_button_box), g_high_scale_button, FALSE, TRUE, 0);
    gtk_box_pack_end (GTK_BOX(g_high_scale_button_box), g_high_scale_button_label, FALSE, TRUE, 0);
    gtk_widget_show (g_high_scale_button_label);
    gtk_widget_show (g_high_scale_button);

    GtkWidget   *FCSControlBox;
    FCSControlBox = gtk_vbox_new (FALSE, 2);
    GtkWidget   *fcs_start;
    fcs_start = gtk_button_new_with_label ("Start FCS");
    GtkWidget   *fcs_stop;
    fcs_stop = gtk_button_new_with_label ("Stop FCS");
    g_signal_connect(G_OBJECT(fcs_0x1_button),"clicked",G_CALLBACK(fcs_channels), &flab);
    g_signal_connect(G_OBJECT(fcs_1x2_button),"clicked",G_CALLBACK(fcs_channels), &flab);
    g_signal_connect(G_OBJECT(fcs_2x0_button),"clicked",G_CALLBACK(fcs_channels), &flab);   
    g_signal_connect(G_OBJECT(fcs_start),"clicked",G_CALLBACK(startFCS), &flab);
    g_signal_connect(G_OBJECT(fcs_stop),"clicked",G_CALLBACK(haltFCS), &flab);
    gtk_box_pack_start (GTK_BOX (FCSControlBox), fcs_start, FALSE, FALSE, 10);
    gtk_box_pack_start (GTK_BOX (FCSControlBox), fcs_stop, FALSE, FALSE, 10);
    gtk_box_pack_start (GTK_BOX (FCSControlBox), g_high_scale_button_box, FALSE, FALSE, 10);
    gtk_box_pack_start (GTK_BOX (FCSControlBox), g_low_scale_button_box, FALSE, FALSE, 10);

	gtk_table_attach(GTK_TABLE(fcs_matrix), FCSControlBox, 3, 4, 1, 2,GTK_SHRINK,GTK_FILL, 4, 4);    
    gtk_widget_show(fcs_start);
    gtk_widget_show(fcs_stop);
    gtk_widget_show(FCSControlBox);
    gtk_widget_show(fcs_freq_frame);
    gtk_widget_show(fcs_freq_box);
    gtk_widget_show(g_high_scale_button_box);
    gtk_widget_show(g_low_scale_button_box);

	updateGofTauPlot = -78;
    gtk_timeout_add(43, fcs_update, &flab);

	/**********Show the notebook pages*************************************/

    gtk_widget_show (the_matrix);
    gtk_widget_show (fcs_matrix);
    gtk_widget_show (burst_matrix);
    gtk_widget_show (window);
    
	/************Done Configuring.  Transfer Control to gtk_main()**********/    
    gtk_main ();
    gdk_threads_leave();
    return 0;
}

void    tim_label_set_text (GtkLabel *label, const gchar *format, ...)	{
     va_list args;
     gchar   *s; 
     g_return_if_fail (GTK_IS_LABEL (label));
     g_return_if_fail (format != NULL);
     if (format == NULL)
     {
          gtk_label_set_text ((GtkLabel*) label, "");
          return;
     }
     va_start (args, format);
     s = g_strdup_vprintf (format, args);
     va_end (args);
     gtk_label_set_text ((GtkLabel*) label, s);
     g_free (s);
}

void    launchAcq  (GtkWidget *button, struct label_data *val)	{
    int i;
    int pix_spacing;
    if(NICard.runLock == 1) {  //Only Run if the card is not already running...
        printf("NI Card is currently open for reading.  No Acquisition started\n");
        return;
    }
    else
    {
        if(val->notyetplotted)  {
            pix_spacing = ceil(val->widget->allocation.width/(NUM_GRAPH_POINTS-1));
            for(i=0; i<NUM_GRAPH_POINTS; i++)	{
				//Write 0 to data array, "0" to display array, and the time-base to the x arrary
                val->graphthis[i].y = val->widget->allocation.height;
                val->graphthis[i].x = val->widget->allocation.width - i*pix_spacing;
                val->graphthis2[i].y = val->widget->allocation.height;
                val->graphthis2[i].x = val->widget->allocation.width - i*pix_spacing;
                val->graphthis3[i].y = val->widget->allocation.height;
                val->graphthis3[i].x = val->widget->allocation.width - i*pix_spacing;
            }
            val->notyetplotted = FALSE;
        }
        (val->time) = 0;
        int rc;
        pthread_t Acq_thread;
        printf("Starting data collection thread\n");
        NICard.runLock = 1;
        rc = pthread_create(&Acq_thread, NULL, AcqNew, (void *) val);
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
}

void    startFCS  (GtkWidget *button, struct label_data *val)
{
    int i;
    int pix_spacing;
	int rc;

    if(NICard.runLock == 1)    {	//Only Run if the card is not already running...
        printf("NI Card is currently open for reading.  No Acquisition started\n");
        return;
    }
    else	{
        pix_spacing = ceil(val->widget->allocation.width/(NUM_GRAPH_POINTS-1));
        val->graph_max_y = 0;
        val->graph_min_y = 1e6;
        for(i=0; i<NUM_GRAPH_POINTS; i++)   {
			//Write 0 to data array, "0" to display array, and the time-base to the x arrary
            val->graphthis[i].y = val->widget->allocation.height;
            val->trace[i] = 0;
            val->graphthis[i].x = i*pix_spacing;
        }
        for(i=0; i< (n*(1+pmax)+1); i++)    GofTau[i] = -100.0f;
        (val->time) = 0;

        pthread_t Acq_thread;
        printf("NICard.correlation_mode = %i\n",NICard.correlation_mode);
        switch( NICard.correlation_mode )   {
            case 0 :    printf("0,1\n");
                        NICard.runLock = 1;
                        printf("Starting data collection thread\n");
                        rc = pthread_create(&Acq_thread, NULL, AcqFCS0x1, (void *) val);
                        break;
            case 1 :    printf("1,2\n");
                        NICard.runLock = 1;
                        printf("Starting data collection thread\n");
                        rc = pthread_create(&Acq_thread, NULL, AcqFCS1x2, (void *) val);
                        break;
            case 2 :    printf("2,0\n");
                        NICard.runLock = 1;
                        printf("Starting data collection thread\n");
                        rc = pthread_create(&Acq_thread, NULL, AcqFCS2x0, (void *) val);
                        break;
            default  : printf( "Software error.  Search for gherkin87867" );
                        break;
        }
    }
	i=rc;
}


void    haltFCS  (GtkWidget *button, void *val)	{
    NICard.run = 0;
    if(NICard.runLock == 1)	{
        printf("NI Card stopping...\n");
    }
    else	printf("NI Card was not running.\n");
    return;
}

void    startBURST  (GtkWidget *button, struct label_data *val)	{
    int i;
    int pix_spacing;
    if(NICard.runLock == 1)	{	//Only Run if the card is not already running...
        printf("NI Card is currently open for reading.  No Acquisition started\n");
        return;
    }
    else	{
        if(val->notyetplotted)	{
            pix_spacing = ceil(val->widget->allocation.width/(NUM_GRAPH_POINTS-1));
            for(i=0; i<NUM_GRAPH_POINTS; i++)
            {
				//Write 0 to data array, "0" to display array, and the time-base to the x arrary
                val->graphthis[i].y = val->widget->allocation.height;
                val->trace[i] = 0;
                val->graphthis[i].x = val->widget->allocation.width - i*pix_spacing;
            }
            val->notyetplotted = FALSE;
        }
        (val->time) = 0;
        int rc;
        pthread_t Acq_thread;
        printf("Starting data collection thread\n");
        NICard.runLock = 1;
        rc = pthread_create(&Acq_thread, NULL, AcqBURSTNew, (void *) val);
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
}

void    haltBURST  (GtkWidget *button, void *val)	{
    NICard.run = 0;
    if(NICard.runLock == 1)	printf("NI Card stopping...\n");
    else	printf("NI Card was not running.\n");
    return;
}

void    closeAcq  (GtkWidget *button, void *val)	{
    NICard.run = 0;
    if(NICard.runLock == 1)	printf("NI Card stopping...\n");
    else	printf("NI Card was not running.\n");
    return;
}


void    reset_yNew  (GtkWidget *button, void *ii)  {
    struct label_data *val;
    val = ii;
    int max_y,min_y,i;
    char *markup;
    max_y = 0;
    min_y = 1e9;
    for(i=0; i<NUM_GRAPH_POINTS; i++)   {
        if(val->chan0[i] > max_y) max_y = val->chan0[i];
        if(val->chan1[i] > max_y) max_y = val->chan1[i];
        if(val->chan2[i] > max_y) max_y = val->chan2[i];
        if(val->chan0[i] < min_y) min_y = val->chan0[i];
        if(val->chan1[i] < min_y) min_y = val->chan1[i];
        if(val->chan2[i] < min_y) min_y = val->chan2[i];
    }
    val->graph_max_y = (int)(max_y*1.1);
    val->graph_max_y = val->graph_max_y + 1;
    val->graph_min_y = (int)(min_y*0.9);
    printf("New Y-Range: %i - %i\n", val->graph_max_y, val->graph_min_y);
#ifdef WhiteFonts
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2e</span>", (float)val->graph_max_y);
    gtk_label_set_markup (val->ymax_label, markup);
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2e</span>", (float)val->graph_min_y);
    gtk_label_set_markup (val->ymin_label, markup);
#else
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">%6.2e</span>", (float)val->graph_max_y);
    gtk_label_set_markup (val->ymax_label, markup);
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">%6.2e</span>", (float)val->graph_min_y);
    gtk_label_set_markup (val->ymin_label, markup);
#endif
    return;
}

gint    cb_reset_y(void *ii) {
    struct label_data *val;
    val = ii;
    int max_y,min_y,i;
    char *markup;
    max_y = 0;
    min_y = 1e9;
	if (NICard.runLock == 1 && bursts_running == 0)	{
		for(i=0; i<NUM_GRAPH_POINTS; i++)   {
			if(val->chan0[i] > max_y) max_y = val->chan0[i];
			if(val->chan1[i] > max_y) max_y = val->chan1[i];
			if(val->chan2[i] > max_y) max_y = val->chan2[i];
			if(val->chan0[i] < min_y) min_y = val->chan0[i];
			if(val->chan1[i] < min_y) min_y = val->chan1[i];
			if(val->chan2[i] < min_y) min_y = val->chan2[i];
		}
		val->graph_max_y = (int)(max_y*1.1);
		val->graph_max_y = val->graph_max_y + 1;
		val->graph_min_y = (int)(min_y*0.9);
		printf("New Y-Range: %i - %i\n", val->graph_max_y, val->graph_min_y);
#ifdef WhiteFonts
        markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2e</span>", (float)val->graph_max_y);
		gtk_label_set_markup (val->ymax_label, markup);
		markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2e</span>", (float)val->graph_min_y);
		gtk_label_set_markup (val->ymin_label, markup);
#else
        markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">%6.2e</span>", (float)val->graph_max_y);
		gtk_label_set_markup (val->ymax_label, markup);
		markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">%6.2e</span>", (float)val->graph_min_y);
		gtk_label_set_markup (val->ymin_label, markup);
#endif
	}
    return 1;
}

gint    cb_updateNew(void *ii)  {
    struct  label_data *val;
    val     = ii;
    int     i,max_y,min_y;
    
	//Check for need to reset graph ranges
    if(val->chan0[0] > val->graph_max_y || val->chan1[0] > val->graph_max_y || val->chan2[0] > val->graph_max_y || val->chan0[0] < val->graph_min_y || val->chan1[0] < val->graph_min_y || val->chan2[0] < val->graph_min_y)  {
        char *markup1;
        max_y = 0;
        min_y = 1e9;
        for(i=0; i<NUM_GRAPH_POINTS; i++)   {
            if(val->chan0[i] > max_y) max_y = val->chan0[i];
            if(val->chan1[i] > max_y) max_y = val->chan1[i];
            if(val->chan2[i] > max_y) max_y = val->chan2[i];
            if(val->chan0[i] < min_y) min_y = val->chan0[i];
            if(val->chan1[i] < min_y) min_y = val->chan1[i];
            if(val->chan2[i] < min_y) min_y = val->chan2[i];
        }
        val->graph_max_y = (int)(max_y*1.1);
        val->graph_min_y = (int)(min_y*0.9);
#ifdef WhiteFonts
        markup1 = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2e</span>", (float)val->graph_max_y);
        gtk_label_set_markup (val->ymax_label, markup1);
        markup1 = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#FFFFFF\" weight=\"bold\">%6.2e</span>", (float)val->graph_min_y);
        gtk_label_set_markup (val->ymin_label, markup1);
#else
        markup1 = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">%6.2e</span>", (float)val->graph_max_y);
        gtk_label_set_markup (val->ymax_label, markup1);
        markup1 = g_markup_printf_escaped ("<span font_desc=\"Helvetica 12\" foreground=\"#000000\" weight=\"bold\">%6.2e</span>", (float)val->graph_min_y);
        gtk_label_set_markup (val->ymin_label, markup1);
#endif
        g_free(markup1);
    }
    
    for(i=0; i<NUM_GRAPH_POINTS; i++)   {
        val->graphthis[i].y =  (val->widget->allocation.height-((val->widget->allocation.height*(val->chan0[i] - val->graph_min_y))/(val->graph_max_y - val->graph_min_y)));
        val->graphthis2[i].y = (val->widget->allocation.height-((val->widget->allocation.height*(val->chan1[i] - val->graph_min_y))/(val->graph_max_y - val->graph_min_y)));
        val->graphthis3[i].y = (val->widget->allocation.height-((val->widget->allocation.height*(val->chan2[i] - val->graph_min_y))/(val->graph_max_y - val->graph_min_y)));
    }

    gdk_draw_rectangle (pixmap, val->widget->style->black_gc, TRUE, 0, 0, val->widget->allocation.width, val->widget->allocation.height);
    gdk_draw_lines(pixmap,cyan_gc,val->graphthis, NUM_GRAPH_POINTS);
    gdk_draw_lines(pixmap,tmr_gc,val->graphthis2, NUM_GRAPH_POINTS);
    gdk_draw_lines(pixmap,red_gc,val->graphthis3, NUM_GRAPH_POINTS);
    gdk_draw_drawable (val->widget->window, val->widget->style->fg_gc[GTK_WIDGET_STATE (val->widget)], pixmap, 0, 0, 0, 0, val->widget->allocation.width, val->widget->allocation.height);    
    return 1;
}

gint    cb_update_text_New(void *ii)    {
    struct  label_data *val;
    val     = ii;
    int     i;
    float   runningsum0;
    float   runningsum1;
    float   runningsum2;
    char *markup;
    runningsum0 = 0.0f;
    runningsum1 = 0.0f;
    runningsum2 = 0.0f;
    for(i=0; i<tune_points_to_average; i++)   {
        runningsum0 += val->chan0[i];
        runningsum1 += val->chan1[i];
        runningsum2 += val->chan2[i];
    }
    runningsum0 = runningsum0 / (float) tune_points_to_average;
    runningsum1 = runningsum1 / (float) tune_points_to_average;
    runningsum2 = runningsum2 / (float) tune_points_to_average;
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 56\" foreground=\"#00C0FF\" weight=\"bold\">%6.6f </span><span font_desc=\"Helvetica 56\" foreground=\"#FFFF00\" weight=\"bold\">%6.6f  </span><span font_desc=\"Helvetica 56\" foreground=\"#FF0000\" weight=\"bold\">%6.6f</span>", runningsum0*0.000001f, runningsum1*0.000001f, runningsum2*0.000001f);
    gtk_label_set_markup (val->long_average, markup);
    markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 56\" foreground=\"#00C0FF\" weight=\"bold\">%6.6f </span><span font_desc=\"Helvetica 56\" foreground=\"#FFFF00\" weight=\"bold\">%6.6f  </span><span font_desc=\"Helvetica 56\" foreground=\"#FF0000\" weight=\"bold\">%6.6f</span>", val->chan0[0]*0.000001f, val->chan1[0]*0.000001f, val->chan2[0]*0.000001f);
    gtk_label_set_markup (val->label, markup);
    g_free (markup);
    return 1;
}

gint    fcs_update(void *ii)    {
    struct label_data *val;
    val = ii;
    int i;
    float   max_yfloat, min_yfloat;
    float   pix_per_decade = 0.0f;
    char    *markup;
    markup = g_markup_printf_escaped (" ");
    if(val->time_past < val->time)  {
		for(i=0; i<NUM_GRAPH_POINTS; i++) {
            if(val->fcs_update_data1[i] > val->graph_max_y)    val->graph_max_y = val->fcs_update_data1[i];
        	if(val->fcs_update_data1[i] < val->graph_min_y)    val->graph_min_y = val->fcs_update_data1[i];
			if(val->fcs_update_data2[i] > val->graph_max_y)    val->graph_max_y = val->fcs_update_data2[i];
        	if(val->fcs_update_data2[i] < val->graph_min_y)    val->graph_min_y = val->fcs_update_data2[i];
        }
        if(val->data > val->graph_max_y)    val->graph_max_y = val->data;
        if(val->data < val->graph_min_y)    val->graph_min_y = val->data;
	    for(i=0; i<NUM_GRAPH_POINTS; i++) {
            val->graphthis[i].y = (5+val->widget->allocation.height-((val->widget->allocation.height*(val->fcs_update_data1[i] - val->graph_min_y))/(2.0*(val->graph_max_y - val->graph_min_y))));
        }
        val->time_past = val->time;
        gdk_draw_rectangle (pixmap_fcs_khz, val->widget->style->black_gc, TRUE, 0, 0, val->widget->allocation.width, val->widget->allocation.height);
		gdk_draw_lines(pixmap_fcs_khz,val->widget->style->white_gc,val->graphthis, NUM_GRAPH_POINTS);
        gdk_draw_drawable (val->widget->window, val->widget->style->fg_gc[GTK_WIDGET_STATE (val->widget)], pixmap_fcs_khz, 0, 0, 0, 0, val->widget->allocation.width, val->widget->allocation.height);    
    }
	if (updateGofTauPlot > -1)  {
        int     *incmask;
        int     *decmask;
        int     p,m;
        int     y_axis_height = 20;
        incmask  = (int*  )malloc((n*(1+pmax)+1)*sizeof(int  ));
        decmask  = (int*  )malloc((n*(1+pmax)+1)*sizeof(int  ));
        GdkPoint    graphthis[(n*(1+pmax)+1)];
        GdkPoint    graphthis1[(n*(1+pmax)+1)];
        GdkPoint    graphthis2[(n*(1+pmax)+1)];
        gdk_draw_rectangle (pixmap_fcs_G, val->widget_the_second->style->black_gc, TRUE, 0, 0, val->widget_the_second->allocation.width, val->widget_the_second->allocation.height);

		/*figure out the x axes*/
        float time_min, time_max, freq_offset;
        time_min = 16.0f / 20000000.0f;
        time_max = time_min * (pow(2,pmax-1)*(2*n+1));
        time_min = powf(10.0f,floorf(log10f(time_min)));
        time_min = 0.000001f;
        freq_offset = time_min;
        time_max = powf(10.0f,ceilf(log10f(time_max)));
        time_max = 0.1f;
        pix_per_decade = (float)val->widget_the_second->allocation.width / log10f(time_max / time_min);
        freq_offset = (pix_per_decade * log10f(freq_offset/time_min));
        int j;
        for(i=0; i<(int)log10f(time_max / time_min); i++)   {
            graphthis[30*i].x = (int)(pix_per_decade * i);
            graphthis[30*i+1].x = graphthis[30*i].x;
            graphthis[30*i+2].x = graphthis[30*i].x;
            graphthis[30*i].y = val->widget_the_second->allocation.height-20;
            graphthis[30*i+1].y = val->widget_the_second->allocation.height-10;
            graphthis[30*i+2].y = val->widget_the_second->allocation.height-20;
            for(j=1; j<10; j++) {
                graphthis[30*i+3*j].x = (int)(pix_per_decade * (log10f(pow(10,i)*(1+j))));
                graphthis[30*i+3*j+1].x = graphthis[30*i+3*j].x;
                graphthis[30*i+3*j+2].x = graphthis[30*i+3*j].x;
                graphthis[30*i+3*j].y = val->widget_the_second->allocation.height-20;
                graphthis[30*i+3*j+1].y = val->widget_the_second->allocation.height-15;
                graphthis[30*i+3*j+2].y = val->widget_the_second->allocation.height-20;
            }
        }
        gdk_draw_lines(pixmap_fcs_G,val->widget_the_second->style->white_gc,graphthis, 30*i);
        i = 0;
        PangoLayout *layout;
        PangoContext *context;
        PangoFontDescription *desc;
        context = gdk_pango_context_get_for_screen (gdk_drawable_get_screen (pixmap_fcs_G));
        layout = pango_layout_new (context);
        desc = pango_font_description_from_string ("Helvetica 14");
        pango_layout_set_font_description (layout, desc);
        
        if(time_min < 2e-6) {
            pango_layout_set_text (layout, "1us", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        if(time_min < 2e-5) {
            pango_layout_set_text (layout, "10us", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        if(time_min < 2e-4) {
            pango_layout_set_text (layout, "100us", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        if(time_min < 2e-3) {
            pango_layout_set_text (layout, "1ms", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        if(time_min < 2e-2) {
            pango_layout_set_text (layout, "10ms", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        if(time_min < 2e-1) {
            pango_layout_set_text (layout, "100ms", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        if(time_min < 2e-0) {
            pango_layout_set_text (layout, "1s", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        if(time_min < 2e+1) {
            pango_layout_set_text (layout, "10s", -1);
            gdk_draw_layout(pixmap_fcs_G, val->widget_the_second->style->white_gc, i * pix_per_decade + 2, val->widget_the_second->allocation.height-15, layout);
            i++;
        }
        
        max_yfloat = gtau_max;
        min_yfloat = gtau_min;

        int np1     = n+1;
        for(i=0; i<np1; i++)    {
            incmask[i] = i+1;
            decmask[i] = i;
        }
        for(p=0; p<pmax; p++)   {
            for(m=1; m<np1; m++)    {
                incmask[i-1] = (pow(2,p)*(n+m));
                decmask[i] = (pow(2,p)*(n+m));
                i++;
            }
        }
        int temp;
        for(i=1; i<(n*(1+pmax)+1); i++)  {
            graphthis[i-1].y = (val->widget_the_second->allocation.height-y_axis_height-(((val->widget_the_second->allocation.height-y_axis_height)*(GofTau_a[i] - min_yfloat))/(max_yfloat - min_yfloat)));
            temp = (int)(pix_per_decade*log10f((incmask[i]+decmask[i]-1)*0.5f));
            graphthis[i-1].x = temp + freq_offset;
            graphthis1[i-1].y = (val->widget_the_second->allocation.height-y_axis_height-(((val->widget_the_second->allocation.height-y_axis_height)*(GofTau_b[i] - min_yfloat))/(max_yfloat - min_yfloat)));
            graphthis1[i-1].x = temp + freq_offset;
            graphthis2[i-1].y = (val->widget_the_second->allocation.height-y_axis_height-(((val->widget_the_second->allocation.height-y_axis_height)*(GofTau_c[i] - min_yfloat))/(max_yfloat - min_yfloat)));
            graphthis2[i-1].x = temp + freq_offset;
        }
        gdk_draw_lines(pixmap_fcs_G,cyan_gc,graphthis, (n*(1+pmax)-1));
        gdk_draw_lines(pixmap_fcs_G,tmr_gc,graphthis1, (n*(1+pmax)-1));
        gdk_draw_lines(pixmap_fcs_G,red_gc,graphthis2, (n*(1+pmax)-1));
        gdk_draw_drawable (val->widget_the_second->window, val->widget_the_second->style->fg_gc[GTK_WIDGET_STATE (val->widget_the_second)],pixmap_fcs_G, 0, 0, 0, 0, val->widget_the_second->allocation.width, val->widget_the_second->allocation.height);
        free(incmask);
        free(decmask);
    }
    
		//Update Int, N, E, stdev...
		markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 56\" foreground=\"#00C0FF\" weight=\"bold\">%6.6f </span><span font_desc=\"Helvetica 56\" foreground=\"#FFFF00\" weight=\"bold\">%6.6f  </span><span font_desc=\"Helvetica 56\" foreground=\"#FF0000\" weight=\"bold\">%6.6f</span>", val->fcs_tune_intensity[0], val->fcs_tune_intensity[1], val->fcs_tune_intensity[2]);
        gtk_label_set_markup (val->fcs_stats_intensity, markup);
		markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 24\" foreground=\"#00C0FF\" weight=\"bold\">N\t%6.2f </span><span font_desc=\"Helvetica 24\" foreground=\"#FFFF00\" weight=\"bold\">%6.2f  </span><span font_desc=\"Helvetica 24\" foreground=\"#FF0000\" weight=\"bold\">%6.2f</span>", val->fcs_tune_n[0], val->fcs_tune_n[1], val->fcs_tune_n[2]);
        gtk_label_set_markup (val->fcs_stats_n, markup);
		markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 24\" foreground=\"#00C0FF\" weight=\"bold\">Y\t%6.2f </span><span font_desc=\"Helvetica 24\" foreground=\"#FFFF00\" weight=\"bold\">%6.2f  </span><span font_desc=\"Helvetica 24\" foreground=\"#FF0000\" weight=\"bold\">%6.2f</span>", val->fcs_tune_e[0], val->fcs_tune_e[1], val->fcs_tune_e[2]);
        gtk_label_set_markup (val->fcs_stats_e, markup);
		markup = g_markup_printf_escaped ("<span font_desc=\"Helvetica 24\" foreground=\"#00C0FF\" weight=\"bold\">S/N\t%6.2f </span><span font_desc=\"Helvetica 24\" foreground=\"#FFFF00\" weight=\"bold\">%6.2f  </span><span font_desc=\"Helvetica 24\" foreground=\"#FF0000\" weight=\"bold\">%6.2f</span>", val->fcs_tune_stdev[0], val->fcs_tune_stdev[1], val->fcs_tune_stdev[2]);
        gtk_label_set_markup (val->fcs_stats_stdev, markup);

    g_free (markup);

    return 1;
}

gint    burst_update(void *ii)  {
    struct label_data *val;
    val = ii;
    int i;
	int vertical_offset = val->widget1->allocation.height;
	float yscale = vertical_offset / (float) burst_scale;
	GdkPoint    graphthis[3*BURST_FRAME_WIDTH];
	GdkPoint    graphthis2[3*BURST_FRAME_WIDTH];
	GdkPoint    graphthis3[3*BURST_FRAME_WIDTH];
	if(bursts_running == 1)	{
		for(i=0; i< BURST_FRAME_WIDTH; i++)        {
			graphthis[3*i].y = 0;
			graphthis[3*i].x = i;
			graphthis[3*i+1].y = burst_data[i] * yscale;
			graphthis[3*i+1].x = i;
			graphthis[3*i+2].y = 0;
			graphthis[3*i+2].x = i;
			
			graphthis2[3*i].y = 0;
			graphthis2[3*i].x = i;
			graphthis2[3*i+1].y = burst_data2[i] * yscale;
			graphthis2[3*i+1].x = i;
			graphthis2[3*i+2].y = 0;
			graphthis2[3*i+2].x = i;
			
			graphthis3[3*i].y = 0;
			graphthis3[3*i].x = i;
			graphthis3[3*i+1].y = burst_data3[i] * yscale;
			graphthis3[3*i+1].x = i;
			graphthis3[3*i+2].y = 0;
			graphthis3[3*i+2].x = i;
		}
		gdk_draw_rectangle (pixmap_burst1, val->widget1->style->black_gc, TRUE, 0, 0, val->widget1->allocation.width, val->widget1->allocation.height);
		gdk_draw_lines(pixmap_burst1,cyan_gc,graphthis, 3*BURST_FRAME_WIDTH);
		gdk_draw_drawable (val->widget1->window, val->widget1->style->fg_gc[GTK_WIDGET_STATE (val->widget1)], pixmap_burst1, 0, 0, 0, 0, val->widget1->allocation.width, val->widget1->allocation.height);    
	
		gdk_draw_rectangle (pixmap_burst2, val->widget2->style->black_gc, TRUE, 0, 0, val->widget2->allocation.width, val->widget2->allocation.height);
		gdk_draw_lines(pixmap_burst2,tmr_gc,graphthis2, 3*BURST_FRAME_WIDTH);
		gdk_draw_drawable (val->widget2->window, val->widget2->style->fg_gc[GTK_WIDGET_STATE (val->widget2)], pixmap_burst2, 0, 0, 0, 0, val->widget2->allocation.width, val->widget2->allocation.height);    
	
		gdk_draw_rectangle (pixmap_burst3, val->widget3->style->black_gc, TRUE, 0, 0, val->widget3->allocation.width, val->widget3->allocation.height);
		gdk_draw_lines(pixmap_burst3,red_gc,graphthis3, 3*BURST_FRAME_WIDTH);
		gdk_draw_drawable (val->widget3->window, val->widget3->style->fg_gc[GTK_WIDGET_STATE (val->widget3)], pixmap_burst3, 0, 0, 0, 0, val->widget3->allocation.width, val->widget3->allocation.height);   
	}
    return 1;
}

void *AcqNew(void *threadarg)   {
    struct  label_data *val;
    int     j;
    int     clock_freq = 20000000;
    val = threadarg;
    printf("Starting Acq. thread\n");
	int         error=0;
	TaskHandle  taskHandle=0;

    uInt8   data[TUNE_BUFFER];
    unsigned short *one;
    unsigned short *two;
    unsigned short *thr;
    unsigned short *result_one;
    unsigned short *result_two;
    unsigned short *result_thr;
    unsigned int *remainder;
    one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_one = (unsigned short *) calloc(TEST_BUFFER,sizeof(unsigned short));
    result_two = (unsigned short *) calloc(TEST_BUFFER,sizeof(unsigned short));
    result_thr = (unsigned short *) calloc(TEST_BUFFER,sizeof(unsigned short));
    remainder = (unsigned int *) calloc(2,sizeof(unsigned int));
	int32		sampsRead;
	char        errBuff[2048]={'\0'};
	float      sumit1 = 0.0f;
	float      sumit2 = 0.0f;
    float      sumit3 = 0.0f;
	//DAQmxResetDevice ("Dev1");
	DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
    DAQmxErrChk (DAQmxCreateDIChan(taskHandle,"Dev1/port0/line0:7", "" , DAQmx_Val_ChanForAllLines));
    DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"/Dev1/20MhzTimebase", clock_freq, DAQmx_Val_Falling, DAQmx_Val_ContSamps, TUNE_BUFFER));

    DAQmxErrChk (DAQmxStartTask(taskHandle));
    remainder[0] = 0x0;
    remainder[1] = 0x0;
    while (NICard.run == 1)     {
		//  Acquire 2048 samples of 50ns each per acquistion.  Average 512 acquisitions to give a humanly-readable ~19Hz
        for(j=0; j< 512; j++)    {
            DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
            bits_2048_2048_3_c(&result_one[0], &result_two[0], &result_thr[0], &remainder[0], &data[0]);
            sumit1 += (float)result_one[0];
            sumit2 += (float)result_two[0];
            sumit3 += (float)result_thr[0];                 
        }                                                

        for(j=NUM_GRAPH_POINTS-1; j>0; j--)   {
            val->chan0[j] = val->chan0[j-1];
            val->chan1[j] = val->chan1[j-1];
            val->chan2[j] = val->chan2[j-1];
        }
        
        val->chan0[0] = sumit1*19.07348f;
        val->chan1[0] = sumit2*19.07348f;
        val->chan2[0] = sumit3*19.07348f;
        sumit1 = 0.0f;
        sumit2 = 0.0f;
        sumit3 = 0.0f;
	}

    Error:
	if( DAQmxFailed(error) )   {
		DAQmxGetExtendedErrorInfo(errBuff,2048);
        if( taskHandle!=0 ) 
        {
            DAQmxStopTask(taskHandle);
            DAQmxClearTask(taskHandle);
        }
		printf("*** Card Error:  Attempting to reset card, if Seg. Fault follows, reboot computer ***\n");
		DAQmxResetDevice ("Dev1");
        NICard.runLock = 0;
        NICard.run = 1;
    }
	if( taskHandle!=0 )		{
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
	}
	if( DAQmxFailed(error) )	printf("DAQmx Error (11): %s\n",errBuff);
    free(one);
    free(two);    
    free(thr);
    free(remainder);
    NICard.runLock = 0;
    NICard.run = 1;
    printf("NI Card stopped..\n");
    pthread_exit(NULL);
}

void *AcqFCS0x1(void *threadarg)	{
    printf("Starting FCS, please wait ~15s for data to appear on the screen...\n");
    struct  label_data *val;
    int     j,q;
    int     i=0;
    val = threadarg;
    int     clock_freq = 20000000;
    float time_quantum = 8.0e-7f;

	NICard.clock_freq_scalar_value = clock_freq/16;
    updateGofTauPlot = 5;
	int         error=0;
	TaskHandle  taskHandle=0;

    uInt8   data[TUNE_BUFFER];
    unsigned short *one;
    unsigned short *two;
    unsigned short *thr;
    unsigned int *remainder;
    unsigned short *result_one;
    unsigned short *result_two;
    unsigned short *result_thr;
    one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    remainder = (unsigned int *) calloc(2,sizeof(unsigned int));
	int circle, full_circle;
	unsigned short*		data_0;
	unsigned short* 	data_1;
	unsigned short* 	data_2;
	int passed = 9;
	int amt_data	= passed * n * (int) powl(2,pmax);
	data_0 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	data_1 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	data_2 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	float gav0b0[circ_buffer_size*veclen];
    float gav0b1[circ_buffer_size*veclen];
    float gav1b1[circ_buffer_size*veclen];
    float gav2b2[circ_buffer_size*veclen];    
    float gstd0b0[circ_buffer_size*veclen];
    float gstd0b1[circ_buffer_size*veclen];
    float gstd1b1[circ_buffer_size*veclen];
    float gstd2b2[circ_buffer_size*veclen];
	float counts[3];
    for(i=0; i<circ_buffer_size*veclen; i++)	gav0b0[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav0b1[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav1b1[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav2b2[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd0b0[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd0b1[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd1b1[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd2b2[i] = 0;
	float cb_gav0b0[veclen];
    float cb_gav0b1[veclen];
    float cb_gav1b1[veclen];
    float cb_gav2b2[veclen];    
    float cb_gstd0b0[veclen];
    float cb_gstd0b1[veclen];
    float cb_gstd1b1[veclen];
    float cb_gstd2b2[veclen];
    float Navg_0b0;
    float Navg_1b1;
    float Navg_2b2;

	int32		sampsRead;
	char        errBuff[2048]={'\0'};
	DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
    DAQmxErrChk (DAQmxCreateDIChan(taskHandle,"Dev1/port0/line0:7", "" , DAQmx_Val_ChanForAllLines));
    DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"/Dev1/20MhzTimebase", clock_freq, DAQmx_Val_Falling, DAQmx_Val_ContSamps, TUNE_BUFFER));
    DAQmxErrChk (DAQmxStartTask(taskHandle));
    remainder[0] = 0x0;
    remainder[1] = 0x0;
        
    i=0;
	for(i=0; i<NUM_GRAPH_POINTS; i++)	val->fcs_update_data1[i] = 0;
	for(i=0; i<NUM_GRAPH_POINTS; i++)	val->fcs_update_data2[i] = 0;

	//Prime the buffers!
	DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
	bits_2048_16_3_c(&result_one[0], &result_two[0], &result_thr[0], &remainder[0], &data[0]); 

	DAQmxStopTask(taskHandle);	

	full_circle = 0;
    while (NICard.run == 1)	{
		full_circle++;
		circle = (full_circle%circ_buffer_size)*veclen;
		counts[0] = 0.0f;
		counts[1] = 0.0f;
		counts[2] = 0.0f;
		
		DAQmxErrChk (DAQmxStartTask(taskHandle));	//collect algorquanta data
		for(q=0; q<amt_data; q+= TUNE_BUFFER/16)	{	
			DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
			bits_2048_16_3_c(&data_0[q], &data_1[q], &data_2[q], &remainder[0], &data[0]); 
		}
		DAQmxStopTask(taskHandle);
		
		ShatzelDriver_A(data_0, data_1, data_2, amt_data, time_quantum, &gav0b0[circle], &gav0b1[circle], &gav1b1[circle], &gav2b2[circle], &gstd0b0[circle], &gstd0b1[circle], &gstd1b1[circle], &gstd2b2[circle], &counts[0]);
		
		for(i=0; i<veclen; i++)	cb_gav0b0[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav0b0[i%veclen] += gav0b0[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav0b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav0b1[i%veclen] += gav0b1[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav1b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav1b1[i%veclen] += gav1b1[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav2b2[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav2b2[i%veclen] += gav2b2[i] / (float)circ_buffer_size;
		
		for(i=0; i<veclen; i++)	cb_gstd0b0[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd0b0[i%veclen] += powf((gav0b0[i]-cb_gav0b0[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd0b0[i] = powf(cb_gstd0b0[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd0b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd0b1[i%veclen] += powf((gav0b1[i]-cb_gav0b1[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd0b1[i] = powf(cb_gstd0b1[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd1b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd1b1[i%veclen] += powf((gav1b1[i]-cb_gav1b1[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd1b1[i] = powf(cb_gstd1b1[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd2b2[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd2b2[i%veclen] += powf((gav2b2[i]-cb_gav2b2[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd2b2[i] = powf(cb_gstd2b2[i],0.5f);
		
		float	sbr[3];	//Average (std. dev(G(tau))) over all n tau's in the p=1 level
		float	sav_0[n];
		float	sav_1[n];
		float	sav_2[n];
		float	ssd_0[n];
		float	ssd_1[n];
		float	ssd_2[n];
		for(i=0; i<n; i++)	sav_0[i] = 0.0f;
		for(i=0; i<n; i++)	sav_1[i] = 0.0f;
		for(i=0; i<n; i++)	sav_2[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_0[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_1[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_2[i] = 0.0f;
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_0[i] += gav0b0[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_1[i] += gav1b1[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_2[i] += gav2b2[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_0[i] += powf(gav0b0[j*veclen+2*n+i]-sav_0[i], 2.0) / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_1[i] += powf(gav1b1[j*veclen+2*n+i]-sav_1[i], 2.0) / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_2[i] += powf(gav2b2[j*veclen+2*n+i]-sav_2[i], 2.0) / (float)circ_buffer_size; }
		sbr[0] = 0.0f;
		sbr[1] = 0.0f;
		sbr[2] = 0.0f;
		for(i=0; i<n; i++)	sbr[0] += (sav_0[i]-1.0f) / (powf(ssd_0[i], 0.5f) * (float)n);
		for(i=0; i<n; i++)	sbr[1] += (sav_1[i]-1.0f) / (powf(ssd_1[i], 0.5f) * (float)n);
		for(i=0; i<n; i++)	sbr[2] += (sav_2[i]-1.0f) / (powf(ssd_2[i], 0.5f) * (float)n);
			
		val->fcs_tune_stdev[0] = sbr[0];
		val->fcs_tune_stdev[1] = sbr[1];
		val->fcs_tune_stdev[2] = sbr[2];
		
		Navg_0b0 = 0.0f;
		Navg_1b1 = 0.0f;
		Navg_2b2 = 0.0f;
		for(i=2*n; i<3*n; i++)	Navg_0b0 += cb_gav0b0[i]/n;
		for(i=2*n; i<3*n; i++)	Navg_1b1 += cb_gav1b1[i]/n;
		for(i=2*n; i<3*n; i++)	Navg_2b2 += cb_gav2b2[i]/n;
		Navg_0b0 = 1.0f/(Navg_0b0-1.0f);
		Navg_1b1 = 1.0f/(Navg_1b1-1.0f);
		Navg_2b2 = 1.0f/(Navg_2b2-1.0f);
		val->fcs_tune_n[0] = Navg_0b0;
		val->fcs_tune_n[1] = Navg_1b1;
		val->fcs_tune_n[2] = Navg_2b2;

		counts[0] = counts[0]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		counts[1] = counts[1]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		counts[2] = counts[2]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		val->fcs_tune_intensity[0] = counts[0];
		val->fcs_tune_intensity[1] = counts[1];
		val->fcs_tune_intensity[2] = counts[2];
		
		val->fcs_tune_e[0] = 1000.0f * counts[0] / Navg_0b0;
		val->fcs_tune_e[1] = 1000.0f * counts[1] / Navg_1b1;
		val->fcs_tune_e[2] = 1000.0f * counts[2] / Navg_2b2;

		for(i=0; i<veclen; i++)	GofTau_a[i] = cb_gav0b0[i] - 1.0;
		for(i=0; i<veclen; i++)	GofTau_b[i] = cb_gav0b1[i] - 1.0;
		for(i=0; i<veclen; i++)	GofTau_c[i] = cb_gav1b1[i] - 1.0;
	}

    free(one);
    free(two);
    free(thr);
    free(remainder);  
    free(data_0);
	free(data_1);
	free(data_2); 

    Error:
	if( DAQmxFailed(error) )   DAQmxGetExtendedErrorInfo(errBuff,2048);
	if( taskHandle!=0 )    {
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
	}
	if( DAQmxFailed(error) )   printf("DAQmx Error (3): %s\n",errBuff);

    NICard.runLock = 0;
    NICard.run = 1;
    val->time = 0;
    val->time_past = 0;
    updateGofTauPlot = -5;
    printf("NI Card stopped.\n");
    pthread_exit(NULL);
}

void *AcqFCS1x2(void *threadarg)	{
    printf("Starting FCS, please wait ~15s for data to appear on the screen...\n");
    struct  label_data *val;
    int     j,q;
    int     i=0;
    val = threadarg;
    int     clock_freq = 20000000;
    float time_quantum = 8.0e-7f;

	NICard.clock_freq_scalar_value = clock_freq/16;
    updateGofTauPlot = 5;
	int         error=0;
	TaskHandle  taskHandle=0;

    uInt8   data[TUNE_BUFFER];
    unsigned short *one;
    unsigned short *two;
    unsigned short *thr;
    unsigned int *remainder;
    unsigned short *result_one;
    unsigned short *result_two;
    unsigned short *result_thr;
    one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    remainder = (unsigned int *) calloc(2,sizeof(unsigned int));
	int circle, full_circle;
	unsigned short*		data_0;
	unsigned short* 	data_1;
	unsigned short* 	data_2;
	int passed = 9;
	int amt_data	= passed * n * (int) powl(2,pmax);
	data_0 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	data_1 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	data_2 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	float gav0b0[circ_buffer_size*veclen];
    float gav0b1[circ_buffer_size*veclen];
    float gav1b1[circ_buffer_size*veclen];
    float gav2b2[circ_buffer_size*veclen];    
    float gstd0b0[circ_buffer_size*veclen];
    float gstd0b1[circ_buffer_size*veclen];
    float gstd1b1[circ_buffer_size*veclen];
    float gstd2b2[circ_buffer_size*veclen];
	float counts[3];
    for(i=0; i<circ_buffer_size*veclen; i++)	gav0b0[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav0b1[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav1b1[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav2b2[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd0b0[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd0b1[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd1b1[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd2b2[i] = 0;
	float cb_gav0b0[veclen];
    float cb_gav0b1[veclen];
    float cb_gav1b1[veclen];
    float cb_gav2b2[veclen];    
    float cb_gstd0b0[veclen];
    float cb_gstd0b1[veclen];
    float cb_gstd1b1[veclen];
    float cb_gstd2b2[veclen];
    float Navg_0b0;
    float Navg_1b1;
    float Navg_2b2;

	int32		sampsRead;
	char        errBuff[2048]={'\0'};
	DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
    DAQmxErrChk (DAQmxCreateDIChan(taskHandle,"Dev1/port0/line0:7", "" , DAQmx_Val_ChanForAllLines));
    DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"/Dev1/20MhzTimebase", clock_freq, DAQmx_Val_Falling, DAQmx_Val_ContSamps, TUNE_BUFFER));
    DAQmxErrChk (DAQmxStartTask(taskHandle));
    remainder[0] = 0x0;
    remainder[1] = 0x0;
        
    i=0;
	for(i=0; i<NUM_GRAPH_POINTS; i++)	val->fcs_update_data1[i] = 0;
	for(i=0; i<NUM_GRAPH_POINTS; i++)	val->fcs_update_data2[i] = 0;

	//Prime the buffers!
	DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
	bits_2048_16_3_c(&result_one[0], &result_two[0], &result_thr[0], &remainder[0], &data[0]); 

	DAQmxStopTask(taskHandle);

	full_circle = 0;
    while (NICard.run == 1)	{
		full_circle++;
		circle = (full_circle%circ_buffer_size)*veclen;
		counts[0] = 0.0f;
		counts[1] = 0.0f;
		counts[2] = 0.0f;
		
		DAQmxErrChk (DAQmxStartTask(taskHandle));	//collect algorquanta data
		for(q=0; q<amt_data; q+= TUNE_BUFFER/16)	{	
			DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
			bits_2048_16_3_c(&data_0[q], &data_1[q], &data_2[q], &remainder[0], &data[0]); 
		}
		DAQmxStopTask(taskHandle);

		ShatzelDriver_B(data_0, data_1, data_2, amt_data, time_quantum, &gav0b0[circle], &gav0b1[circle], &gav1b1[circle], &gav2b2[circle], &gstd0b0[circle], &gstd0b1[circle], &gstd1b1[circle], &gstd2b2[circle], &counts[0]);
		
		for(i=0; i<veclen; i++)	cb_gav0b0[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav0b0[i%veclen] += gav0b0[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav0b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav0b1[i%veclen] += gav0b1[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav1b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav1b1[i%veclen] += gav1b1[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav2b2[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav2b2[i%veclen] += gav2b2[i] / (float)circ_buffer_size;
		
		for(i=0; i<veclen; i++)	cb_gstd0b0[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd0b0[i%veclen] += powf((gav0b0[i]-cb_gav0b0[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd0b0[i] = powf(cb_gstd0b0[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd0b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd0b1[i%veclen] += powf((gav0b1[i]-cb_gav0b1[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd0b1[i] = powf(cb_gstd0b1[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd1b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd1b1[i%veclen] += powf((gav1b1[i]-cb_gav1b1[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd1b1[i] = powf(cb_gstd1b1[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd2b2[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd2b2[i%veclen] += powf((gav2b2[i]-cb_gav2b2[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd2b2[i] = powf(cb_gstd2b2[i],0.5f);
		
		float	sbr[3];	//Average (std. dev(G(tau))) over all n tau's in the p=1 level
		float	sav_0[n];
		float	sav_1[n];
		float	sav_2[n];
		float	ssd_0[n];
		float	ssd_1[n];
		float	ssd_2[n];
		for(i=0; i<n; i++)	sav_0[i] = 0.0f;
		for(i=0; i<n; i++)	sav_1[i] = 0.0f;
		for(i=0; i<n; i++)	sav_2[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_0[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_1[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_2[i] = 0.0f;
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_1[i] += gav0b0[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_2[i] += gav1b1[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_0[i] += gav2b2[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_1[i] += powf(gav0b0[j*veclen+2*n+i]-sav_1[i], 2.0) / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_2[i] += powf(gav1b1[j*veclen+2*n+i]-sav_2[i], 2.0) / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_0[i] += powf(gav2b2[j*veclen+2*n+i]-sav_0[i], 2.0) / (float)circ_buffer_size; }
		sbr[0] = 0.0f;
		sbr[1] = 0.0f;
		sbr[2] = 0.0f;
		for(i=0; i<n; i++)	sbr[0] += (sav_0[i]-1.0f) / (powf(ssd_0[i], 0.5f) * (float)n);
		for(i=0; i<n; i++)	sbr[1] += (sav_1[i]-1.0f) / (powf(ssd_1[i], 0.5f) * (float)n);
		for(i=0; i<n; i++)	sbr[2] += (sav_2[i]-1.0f) / (powf(ssd_2[i], 0.5f) * (float)n);
			
		val->fcs_tune_stdev[0] = sbr[0];
		val->fcs_tune_stdev[1] = sbr[1];
		val->fcs_tune_stdev[2] = sbr[2];
		
		Navg_0b0 = 0.0f;
		Navg_1b1 = 0.0f;
		Navg_2b2 = 0.0f;
		for(i=2*n; i<3*n; i++)	Navg_1b1 += cb_gav0b0[i]/n;
		for(i=2*n; i<3*n; i++)	Navg_2b2 += cb_gav1b1[i]/n;
		for(i=2*n; i<3*n; i++)	Navg_0b0 += cb_gav2b2[i]/n;
		Navg_0b0 = 1.0f/(Navg_0b0-1.0f);
		Navg_1b1 = 1.0f/(Navg_1b1-1.0f);
		Navg_2b2 = 1.0f/(Navg_2b2-1.0f);
		val->fcs_tune_n[0] = Navg_0b0;
		val->fcs_tune_n[1] = Navg_1b1;
		val->fcs_tune_n[2] = Navg_2b2;

		counts[0] = counts[0]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		counts[1] = counts[1]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		counts[2] = counts[2]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		val->fcs_tune_intensity[0] = counts[0];
		val->fcs_tune_intensity[1] = counts[1];
		val->fcs_tune_intensity[2] = counts[2];
		
		val->fcs_tune_e[0] = 1000.0f * counts[0] / Navg_0b0;
		val->fcs_tune_e[1] = 1000.0f * counts[1] / Navg_1b1;
		val->fcs_tune_e[2] = 1000.0f * counts[2] / Navg_2b2;

		for(i=0; i<veclen; i++)	GofTau_a[i] = cb_gav0b0[i] - 1.0;
		for(i=0; i<veclen; i++)	GofTau_b[i] = cb_gav0b1[i] - 1.0;
		for(i=0; i<veclen; i++)	GofTau_c[i] = cb_gav1b1[i] - 1.0;
	}

    free(one);
    free(two);
    free(thr);
    free(remainder);  
    free(data_0);
	free(data_1);
	free(data_2); 

    Error:
	if( DAQmxFailed(error) )   DAQmxGetExtendedErrorInfo(errBuff,2048);
	if( taskHandle!=0 )    {
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
	}
	if( DAQmxFailed(error) )   printf("DAQmx Error (4): %s\n",errBuff);

    NICard.runLock = 0;
    NICard.run = 1;
    val->time = 0;
    val->time_past = 0;
    updateGofTauPlot = -5;
    printf("NI Card stopped.\n");
    pthread_exit(NULL);
}

void *AcqFCS2x0(void *threadarg)	{
    printf("Starting FCS, please wait ~15s for data to appear on the screen...\n");
    struct  label_data *val;
    int     j,q;
    int     i=0;
    val = threadarg;
    int     clock_freq = 20000000;
    float time_quantum = 8.0e-7f;

	NICard.clock_freq_scalar_value = clock_freq/16;
    updateGofTauPlot = 5;
	int         error=0;
	TaskHandle  taskHandle=0;

    uInt8   data[TUNE_BUFFER];
    unsigned short *one;
    unsigned short *two;
    unsigned short *thr;
    unsigned int *remainder;
    unsigned short *result_one;
    unsigned short *result_two;
    unsigned short *result_thr;
    one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    remainder = (unsigned int *) calloc(2,sizeof(unsigned int));
	int circle, full_circle;
	unsigned short*		data_0;
	unsigned short* 	data_1;
	unsigned short* 	data_2;
	int passed = 9;
	int amt_data	= passed * n * (int) powl(2,pmax);
	data_0 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	data_1 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	data_2 = (unsigned short*)malloc(amt_data*sizeof(unsigned short));
	float gav0b0[circ_buffer_size*veclen];
    float gav0b1[circ_buffer_size*veclen];
    float gav1b1[circ_buffer_size*veclen];
    float gav2b2[circ_buffer_size*veclen];    
    float gstd0b0[circ_buffer_size*veclen];
    float gstd0b1[circ_buffer_size*veclen];
    float gstd1b1[circ_buffer_size*veclen];
    float gstd2b2[circ_buffer_size*veclen];
	float counts[3];
    for(i=0; i<circ_buffer_size*veclen; i++)	gav0b0[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav0b1[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav1b1[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gav2b2[i]  = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd0b0[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd0b1[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd1b1[i] = 0;
    for(i=0; i<circ_buffer_size*veclen; i++)	gstd2b2[i] = 0;
	float cb_gav0b0[veclen];
    float cb_gav0b1[veclen];
    float cb_gav1b1[veclen];
    float cb_gav2b2[veclen];    
    float cb_gstd0b0[veclen];
    float cb_gstd0b1[veclen];
    float cb_gstd1b1[veclen];
    float cb_gstd2b2[veclen];
    float Navg_0b0;
    float Navg_1b1;
    float Navg_2b2;

	int32		sampsRead;
	char        errBuff[2048]={'\0'};
	DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
    DAQmxErrChk (DAQmxCreateDIChan(taskHandle,"Dev1/port0/line0:7", "" , DAQmx_Val_ChanForAllLines));
    DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"/Dev1/20MhzTimebase", clock_freq, DAQmx_Val_Falling, DAQmx_Val_ContSamps, TUNE_BUFFER));
    DAQmxErrChk (DAQmxStartTask(taskHandle));
    remainder[0] = 0x0;
    remainder[1] = 0x0;
        
    i=0;
	for(i=0; i<NUM_GRAPH_POINTS; i++)	val->fcs_update_data1[i] = 0;
	for(i=0; i<NUM_GRAPH_POINTS; i++)	val->fcs_update_data2[i] = 0;

	//Prime the buffers!
	DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
	bits_2048_16_3_c(&result_one[0], &result_two[0], &result_thr[0], &remainder[0], &data[0]); 

	DAQmxStopTask(taskHandle);

	full_circle = 0;
    while (NICard.run == 1)	{
		full_circle++;
		circle = (full_circle%circ_buffer_size)*veclen;
		counts[0] = 0.0f;
		counts[1] = 0.0f;
		counts[2] = 0.0f;
		
		DAQmxErrChk (DAQmxStartTask(taskHandle));	//collect algorquanta data
		for(q=0; q<amt_data; q+= TUNE_BUFFER/16)	{	
			DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TUNE_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
			bits_2048_16_3_c(&data_0[q], &data_1[q], &data_2[q], &remainder[0], &data[0]); 
		}
		DAQmxStopTask(taskHandle);
		
		ShatzelDriver_C(data_0, data_1, data_2, amt_data, time_quantum, &gav0b0[circle], &gav0b1[circle], &gav1b1[circle], &gav2b2[circle], &gstd0b0[circle], &gstd0b1[circle], &gstd1b1[circle], &gstd2b2[circle], &counts[0]);
		
		for(i=0; i<veclen; i++)	cb_gav0b0[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav0b0[i%veclen] += gav0b0[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav0b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav0b1[i%veclen] += gav0b1[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav1b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav1b1[i%veclen] += gav1b1[i] / (float)circ_buffer_size;
		for(i=0; i<veclen; i++)	cb_gav2b2[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gav2b2[i%veclen] += gav2b2[i] / (float)circ_buffer_size;
		
		for(i=0; i<veclen; i++)	cb_gstd0b0[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd0b0[i%veclen] += powf((gav0b0[i]-cb_gav0b0[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd0b0[i] = powf(cb_gstd0b0[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd0b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd0b1[i%veclen] += powf((gav0b1[i]-cb_gav0b1[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd0b1[i] = powf(cb_gstd0b1[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd1b1[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd1b1[i%veclen] += powf((gav1b1[i]-cb_gav1b1[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd1b1[i] = powf(cb_gstd1b1[i],0.5f);
		for(i=0; i<veclen; i++)	cb_gstd2b2[i] = 0;
		for(i=0; i<circ_buffer_size*veclen; i++)	cb_gstd2b2[i%veclen] += powf((gav2b2[i]-cb_gav2b2[i%veclen]),2.0) / ((float)circ_buffer_size-1.0f);
		for(i=0; i<veclen; i++)	cb_gstd2b2[i] = powf(cb_gstd2b2[i],0.5f);
		
		float	sbr[3];	//Average (std. dev(G(tau))) over all n tau's in the p=1 level
		float	sav_0[n];
		float	sav_1[n];
		float	sav_2[n];
		float	ssd_0[n];
		float	ssd_1[n];
		float	ssd_2[n];
		for(i=0; i<n; i++)	sav_0[i] = 0.0f;
		for(i=0; i<n; i++)	sav_1[i] = 0.0f;
		for(i=0; i<n; i++)	sav_2[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_0[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_1[i] = 0.0f;
		for(i=0; i<n; i++)	ssd_2[i] = 0.0f;
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_2[i] += gav0b0[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_0[i] += gav1b1[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	sav_1[i] += gav2b2[j*veclen+2*n+i] / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_2[i] += powf(gav0b0[j*veclen+2*n+i]-sav_2[i], 2.0) / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_0[i] += powf(gav1b1[j*veclen+2*n+i]-sav_0[i], 2.0) / (float)circ_buffer_size; }
		for(i=0; i<n; i++)	{ for(j=0; j<circ_buffer_size; j++)	ssd_1[i] += powf(gav2b2[j*veclen+2*n+i]-sav_1[i], 2.0) / (float)circ_buffer_size; }
		sbr[0] = 0.0f;
		sbr[1] = 0.0f;
		sbr[2] = 0.0f;
		for(i=0; i<n; i++)	sbr[0] += (sav_0[i]-1.0f) / (powf(ssd_0[i], 0.5f) * (float)n);
		for(i=0; i<n; i++)	sbr[1] += (sav_1[i]-1.0f) / (powf(ssd_1[i], 0.5f) * (float)n);
		for(i=0; i<n; i++)	sbr[2] += (sav_2[i]-1.0f) / (powf(ssd_2[i], 0.5f) * (float)n);
			
		val->fcs_tune_stdev[0] = sbr[0];
		val->fcs_tune_stdev[1] = sbr[1];
		val->fcs_tune_stdev[2] = sbr[2];
		
		Navg_0b0 = 0.0f;
		Navg_1b1 = 0.0f;
		Navg_2b2 = 0.0f;
		for(i=2*n; i<3*n; i++)	Navg_2b2 += cb_gav0b0[i]/n;
		for(i=2*n; i<3*n; i++)	Navg_0b0 += cb_gav1b1[i]/n;
		for(i=2*n; i<3*n; i++)	Navg_1b1 += cb_gav2b2[i]/n;
		Navg_0b0 = 1.0f/(Navg_0b0-1.0f);
		Navg_1b1 = 1.0f/(Navg_1b1-1.0f);
		Navg_2b2 = 1.0f/(Navg_2b2-1.0f);
		val->fcs_tune_n[0] = Navg_0b0;
		val->fcs_tune_n[1] = Navg_1b1;
		val->fcs_tune_n[2] = Navg_2b2;

		counts[0] = counts[0]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		counts[1] = counts[1]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		counts[2] = counts[2]*((float)passed*20.0f)/(((float)passed-1.0f)*16.0f*amt_data);
		val->fcs_tune_intensity[0] = counts[0];
		val->fcs_tune_intensity[1] = counts[1];
		val->fcs_tune_intensity[2] = counts[2];
		
		val->fcs_tune_e[0] = 1000.0f * counts[0] / Navg_0b0;
		val->fcs_tune_e[1] = 1000.0f * counts[1] / Navg_1b1;
		val->fcs_tune_e[2] = 1000.0f * counts[2] / Navg_2b2;

		for(i=0; i<veclen; i++)	GofTau_a[i] = cb_gav0b0[i] - 1.0;
		for(i=0; i<veclen; i++)	GofTau_b[i] = cb_gav0b1[i] - 1.0;
		for(i=0; i<veclen; i++)	GofTau_c[i] = cb_gav1b1[i] - 1.0;
	}

    free(one);
    free(two);
    free(thr);
    free(remainder);  
    free(data_0);
	free(data_1);
	free(data_2); 

    Error:
	if( DAQmxFailed(error) )   DAQmxGetExtendedErrorInfo(errBuff,2048);
	if( taskHandle!=0 )    {
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
	}
	if( DAQmxFailed(error) )   printf("DAQmx Error (5): %s\n",errBuff);

    NICard.runLock = 0;
    NICard.run = 1;
    val->time = 0;
    val->time_past = 0;
    updateGofTauPlot = -5;
    printf("NI Card stopped.\n");
    pthread_exit(NULL);
}

void *AcqBURSTNew(void *threadarg)  {
    struct  label_data *val;
    int     j,k,l,bin_num;
    int     i=0;
    val = threadarg;
    int     clock_freq;
    clock_freq = 20000000;
	int         error=0;
	float Sum_one = 0.0f;
	float Sum_two = 0.0f;
	float Sum_thr = 0.0f;
	TaskHandle  taskHandle=0;
    uInt8   data[TUNE_BUFFER];
    unsigned short *one;
    unsigned short *two;
    unsigned short *thr;
    unsigned short *result_one;
    unsigned short *result_two;
    unsigned short *result_thr;
    unsigned int *rem;
    one = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    two = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    thr = (unsigned short *) calloc(TUNE_BUFFER,sizeof(unsigned short));
    result_one = (unsigned short *) calloc(TEST_BUFFER,sizeof(unsigned short));
    result_two = (unsigned short *) calloc(TEST_BUFFER,sizeof(unsigned short));
    result_thr = (unsigned short *) calloc(TEST_BUFFER,sizeof(unsigned short));

    rem = (unsigned int *) calloc(2,sizeof(unsigned int));
	int32		sampsRead;
	char        errBuff[2048]={'\0'};

	bursts_running = 1;
	DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
    DAQmxErrChk (DAQmxCreateDIChan(taskHandle,"Dev1/port0/line0:7", "" , DAQmx_Val_ChanForAllLines));
    DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"/Dev1/20MhzTimebase", clock_freq, DAQmx_Val_Falling, DAQmx_Val_ContSamps, TUNE_BUFFER));
    DAQmxErrChk (DAQmxStartTask(taskHandle));
    rem[0] = 0x0;
    rem[1] = 0x0;
    k = 0;

    switch( NICard.clock_frequency_array_value )   {
            case 0 :    printf("50ns bins\n");
                        while (NICard.run == 1)     {
                            DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TEST_BUFFER, &sampsRead, NULL));
//                            bits_2048_2048_3_c(&one[0], &two[0], &thr[0], &rem[0], &data[0]);
                            bits_2048_1_3_c(&one[0], &two[0], &thr[0], &rem[0], &data[0]);
                            for (j=0; j<TEST_BUFFER; j++)     {
                                burst_data[j] = (float)one[j];
                                burst_data2[j] = (float)two[j];
                                burst_data3[j] = (float)thr[j];
                            }
                            (val->time)++;
                        }          
                        break;
            case 1 :    printf("800ns bins\n");
                        while (NICard.run == 1)     {
                            for(j=0; j< 16; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_16_3_c(&one[0], &two[0], &thr[0], &rem[0], &data[0]);                    
                                for(i=0; i< TEST_BUFFER/16; i++) {                            
                                    burst_data[k] = (float)one[i];
                                    burst_data2[k] = (float)two[i];
                                    burst_data3[k] = (float)thr[i];
                                    k++;
                                }
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 2 :    printf("1.6us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 32;
                            for(j=0; j< bin_num; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_32_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);                    
                                for(i=0; i< TEST_BUFFER/bin_num; i++) {                            
                                    burst_data[k] = (float)result_one[i];
                                    burst_data2[k] = (float)result_two[i];
                                    burst_data3[k] = (float)result_thr[i];
                                    k++;
                                }
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 3 :    printf("3.2us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 64;
                            for(j=0; j< bin_num; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_64_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);                    
                                for(i=0; i< TEST_BUFFER/bin_num; i++) {                            
                                    burst_data[k] = (float)result_one[i];
                                    burst_data2[k] = (float)result_two[i];
                                    burst_data3[k] = (float)result_thr[i];
                                    k++;
                                }
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 4 :    printf("6.4us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 128;
                            for(j=0; j< bin_num; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_128_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);                    
                                for(i=0; i< TEST_BUFFER/bin_num; i++) {                            
                                    burst_data[k] = (float)result_one[i];
                                    burst_data2[k] = (float)result_two[i];
                                    burst_data3[k] = (float)result_thr[i];
                                    k++;
                                }
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 5 :    printf("12.8us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 256;
                            for(j=0; j< bin_num; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_256_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);                    
                                for(i=0; i< TEST_BUFFER/bin_num; i++) {                            
                                    burst_data[k] = (float)result_one[i];
                                    burst_data2[k] = (float)result_two[i];
                                    burst_data3[k] = (float)result_thr[i];
                                    k++;
                                }
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 6 :    printf("25.6us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 512;
                            for(j=0; j< bin_num; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_512_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);                    
                                for(i=0; i< TEST_BUFFER/bin_num; i++) {                            
                                    burst_data[k] = (float)result_one[i];
                                    burst_data2[k] = (float)result_two[i];
                                    burst_data3[k] = (float)result_thr[i];
                                    k++;
                                }
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 7 :    printf("51.2us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 1024;
                            for(j=0; j< bin_num; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_1024_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);                    
                                for(i=0; i< TEST_BUFFER/bin_num; i++) {                            
                                    burst_data[k] = (float)result_one[i];
                                    burst_data2[k] = (float)result_two[i];
                                    burst_data3[k] = (float)result_thr[i];
                                    k++;
                                }
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 8 :    printf("102.4us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 2048;
                            for(j=0; j< bin_num; j++)    {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_2048_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);                                                
                                burst_data[k] = (float)result_one[i];
                                burst_data2[k] = (float)result_two[i];
                                burst_data3[k] = (float)result_thr[i];
                                k++;
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            case 9 :    printf("512us bins\n");
                        while (NICard.run == 1)     {
                            bin_num = 2048;
                            for(j=0; j< bin_num; j++)    {
                                for(l=0; l<5; l++)  {
                                DAQmxErrChk (DAQmxReadDigitalU8(taskHandle, TEST_BUFFER, 1.0, DAQmx_Val_GroupByScanNumber, data, TUNE_BUFFER, &sampsRead, NULL));
                                bits_2048_2048_3_c(&result_one[0], &result_two[0], &result_thr[0], &rem[0], &data[0]);
                                    Sum_one += (float)result_one[i];
                                    Sum_two += (float)result_two[i];
                                    Sum_thr += (float)result_thr[i];                 
                                }                                                
                                burst_data[k] = Sum_one;
                                burst_data2[k] = Sum_two;
                                burst_data3[k] = Sum_thr;
                                k++;
                                Sum_one = 0.0f;
                                Sum_two = 0.0f;
                                Sum_thr = 0.0f;
                            }
                            (val->time)++;
                            k = 0;
                        }
                        break;
            default  : printf( "Software error.  Search for gherkin972" );
                        break;
        }      

    Error:
	if(DAQmxFailed(error))   DAQmxGetExtendedErrorInfo(errBuff,2048);
	if(taskHandle!=0)	{
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
	}
	if(DAQmxFailed(error))   printf("DAQmx Error (8): %s\n",errBuff);

	bursts_running = 0;
    NICard.runLock = 0;
    NICard.run = 1;
    val->time = 1;
    val->time_past = 0;
    printf("NI Card stopped.\n");
    free(one);
    free(two);
    free(thr);
    free(rem);
    free(result_one);
    free(result_two);
    free(result_thr);
    pthread_exit(NULL);
}


void *StreamIt(void *threadarg)   {
    struct stream *val;
    val = threadarg;
    int i;
    for(i=0; i<1e6; i++);
    pthread_mutex_lock (&mutex_file_key);
    *(val->checkout) = 0;
    fwrite(val->array, sizeof(unsigned short), val->len, val->fp);
    pthread_mutex_unlock (&mutex_file_key);
    printf("...done writing\n");
    pthread_exit(NULL);
}

/* Create a new backing pixmap of the appropriate size */
static gboolean configure_event( GtkWidget *widget, GdkEventConfigure *event )
{
    if (pixmap) g_object_unref (pixmap);
    pixmap = gdk_pixmap_new (widget->window, widget->allocation.width, widget->allocation.height, -1);
    gdk_draw_rectangle (pixmap, widget->style->black_gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);
    gdk_color_alloc(gdk_colormap_get_system (), &cyan);
    gdk_color_alloc(gdk_colormap_get_system (), &tmr);
    gdk_color_alloc(gdk_colormap_get_system (), &red);
    cyan_gc = gdk_gc_new( widget->window );
    tmr_gc = gdk_gc_new( widget->window );
    red_gc = gdk_gc_new( widget->window );
    gdk_gc_set_foreground( cyan_gc, &cyan );
    gdk_gc_set_foreground( tmr_gc, &tmr );
    gdk_gc_set_foreground( red_gc, &red );
    return TRUE;
}

/* Redraw the screen from the backing pixmap */
static gboolean expose_event( GtkWidget *widget, GdkEventExpose *event )
{
  gdk_draw_drawable (widget->window,
		     widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		     pixmap,
		     event->area.x, event->area.y,
		     event->area.x, event->area.y,
		     event->area.width, event->area.height);
  return FALSE;
}

static  gboolean configure_event_fcs_G( GtkWidget *widget, GdkEventConfigure *event )
{
  if (pixmap_fcs_G) g_object_unref (pixmap_fcs_G);
  pixmap_fcs_G = gdk_pixmap_new (widget->window, widget->allocation.width, widget->allocation.height, -1);
  gdk_draw_rectangle (pixmap_fcs_G, widget->style->black_gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);
  return TRUE;
}

static  gboolean expose_event_fcs_G( GtkWidget *widget, GdkEventExpose *event )
{
  updateGofTauPlot = 1;
  gdk_draw_drawable (widget->window,
		     widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		     pixmap_fcs_G,
		     event->area.x, event->area.y,
		     event->area.x, event->area.y,
		     event->area.width, event->area.height);
  return FALSE;
}

static  gboolean expose_event_burst1( GtkWidget *widget, GdkEventExpose *event ) {
  gdk_draw_drawable (widget->window, widget->style->fg_gc[GTK_WIDGET_STATE (widget)], pixmap_burst1, event->area.x, event->area.y, event->area.x, event->area.y, event->area.width, event->area.height);
  return FALSE;
}
static  gboolean expose_event_burst2( GtkWidget *widget, GdkEventExpose *event ) {
  gdk_draw_drawable (widget->window, widget->style->fg_gc[GTK_WIDGET_STATE (widget)], pixmap_burst2, event->area.x, event->area.y, event->area.x, event->area.y, event->area.width, event->area.height);
  return FALSE;
}
static  gboolean expose_event_burst3( GtkWidget *widget, GdkEventExpose *event ) {
  gdk_draw_drawable (widget->window, widget->style->fg_gc[GTK_WIDGET_STATE (widget)], pixmap_burst3, event->area.x, event->area.y, event->area.x, event->area.y, event->area.width, event->area.height);
  return FALSE;
}
static  gboolean expose_event_burst4( GtkWidget *widget, GdkEventExpose *event ) {
  gdk_draw_drawable (widget->window, widget->style->fg_gc[GTK_WIDGET_STATE (widget)], pixmap_burst4, event->area.x, event->area.y, event->area.x, event->area.y, event->area.width, event->area.height);
  return FALSE;
}

static  gboolean configure_event_burst1( GtkWidget *widget, GdkEventConfigure *event )   {
    if (pixmap_burst1) g_object_unref (pixmap_burst1);
    pixmap_burst1 = gdk_pixmap_new (widget->window, widget->allocation.width, widget->allocation.height, -1);
    gdk_draw_rectangle (pixmap_burst1, widget->style->black_gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);
    gdk_color_alloc(gdk_colormap_get_system (), &cyan);
    gdk_color_alloc(gdk_colormap_get_system (), &tmr);
    cyan_gc = gdk_gc_new( widget->window );
    tmr_gc = gdk_gc_new( widget->window );
    gdk_gc_set_foreground( cyan_gc, &cyan );
    gdk_gc_set_foreground( tmr_gc, &tmr );
    return TRUE;
}
static  gboolean configure_event_burst2( GtkWidget *widget, GdkEventConfigure *event )   {
    if (pixmap_burst2) g_object_unref (pixmap_burst2);
    pixmap_burst2 = gdk_pixmap_new (widget->window, widget->allocation.width, widget->allocation.height, -1);
    gdk_draw_rectangle (pixmap_burst2, widget->style->black_gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);
    gdk_color_alloc(gdk_colormap_get_system (), &cyan);
    gdk_color_alloc(gdk_colormap_get_system (), &tmr);
    cyan_gc = gdk_gc_new( widget->window );
    tmr_gc = gdk_gc_new( widget->window );
    gdk_gc_set_foreground( cyan_gc, &cyan );
    gdk_gc_set_foreground( tmr_gc, &tmr );
    return TRUE;
}
static  gboolean configure_event_burst3( GtkWidget *widget, GdkEventConfigure *event )   {
    if (pixmap_burst3) g_object_unref (pixmap_burst3);
    pixmap_burst3 = gdk_pixmap_new (widget->window, widget->allocation.width, widget->allocation.height, -1);
    gdk_draw_rectangle (pixmap_burst3, widget->style->black_gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);
    gdk_color_alloc(gdk_colormap_get_system (), &cyan);
    gdk_color_alloc(gdk_colormap_get_system (), &tmr);
    cyan_gc = gdk_gc_new( widget->window );
    tmr_gc = gdk_gc_new( widget->window );
    gdk_gc_set_foreground( cyan_gc, &cyan );
    gdk_gc_set_foreground( tmr_gc, &tmr );
    return TRUE;
}
static  gboolean configure_event_burst4( GtkWidget *widget, GdkEventConfigure *event )   {
    if (pixmap_burst4) g_object_unref (pixmap_burst4);
    pixmap_burst4 = gdk_pixmap_new (widget->window, widget->allocation.width, widget->allocation.height, -1);
    gdk_draw_rectangle (pixmap_burst4, widget->style->black_gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);
    gdk_color_alloc(gdk_colormap_get_system (), &cyan);
    gdk_color_alloc(gdk_colormap_get_system (), &tmr);
    cyan_gc = gdk_gc_new( widget->window );
    tmr_gc = gdk_gc_new( widget->window );
    green_gc = gdk_gc_new( widget->window );
    gdk_gc_set_foreground( cyan_gc, &cyan );
    gdk_gc_set_foreground( tmr_gc, &tmr );
    return TRUE;
}

static gboolean delete_event( GtkWidget *widget, GdkEvent  *event, gpointer   data )    {
    int i;
    NICard.run = 0;
    while(NICard.runLock == 1) {printf("Waiting for NICard to close\n");    for(i=0; i<1e8; i++);}
    gtk_main_quit ();
    return FALSE;
}

void    tune_freq  (GtkWidget *button, void *ii)    {
    struct label_data *val;
    val = ii;
    if(gtk_toggle_button_get_active(val->button_0)) NICard.clock_frequency_array_value = 0;
    if(gtk_toggle_button_get_active(val->button_1)) NICard.clock_frequency_array_value = 1;
    if(gtk_toggle_button_get_active(val->button_2)) NICard.clock_frequency_array_value = 2;
    if(gtk_toggle_button_get_active(val->button_3)) NICard.clock_frequency_array_value = 3;
    if(gtk_toggle_button_get_active(val->button_4)) NICard.clock_frequency_array_value = 4;
    if(gtk_toggle_button_get_active(val->button_5)) NICard.clock_frequency_array_value = 5;
    if(gtk_toggle_button_get_active(val->button_6)) NICard.clock_frequency_array_value = 6;
    if(gtk_toggle_button_get_active(val->button_7)) NICard.clock_frequency_array_value = 7;
    if(gtk_toggle_button_get_active(val->button_8)) NICard.clock_frequency_array_value = 8;
    if(gtk_toggle_button_get_active(val->button_9)) NICard.clock_frequency_array_value = 9;
    return;
}

void    fcs_channels  (GtkWidget *button, void *ii)     {
    struct label_data *val;
    val = ii;
    if(gtk_toggle_button_get_active(val->button_0)) NICard.correlation_mode = 0;
    if(gtk_toggle_button_get_active(val->button_1)) NICard.correlation_mode = 1;
    if(gtk_toggle_button_get_active(val->button_2)) NICard.correlation_mode = 2;
    return;
}

static  void burst_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin)    {
    burst_scale = gtk_spin_button_get_value (spin);
}

static  void g_low_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin)    {
    gtau_min = gtk_spin_button_get_value (spin);
}

static  void g_high_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin)    {
    gtau_max = gtk_spin_button_get_value (spin);
}

static  void tune_scale_button_callback(GtkWidget *widget, GtkSpinButton *spin)    {
    tune_points_to_average = gtk_spin_button_get_value_as_int (spin);
}

int write_shortbursts_to_file (char * filename, unsigned short * burst_array, int len_of_array) {
    //Array is organized as colour0(t), colour1(t), colour0(t+1), colour1(t+1),...
    printf("Writing bursts to file...");
    int complete;
    FILE * fp;
    fp = fopen(filename, "wb");
    if (fp == NULL) {
        printf ("Error:  File Failed To Open!n");
        complete = -666;
    }
    else    {
        complete = fwrite(burst_array, sizeof(unsigned short), len_of_array, fp);
        if (fclose(fp) != 0) {
            printf("Error:  File Failed to Close");
            complete = -666;   
        }
        if (complete == len_of_array) {
            printf("bursts written successfully\n");
            complete = 0;
        }
        else printf("Failed to write\n");
    }
    return complete;
}

int imax (int a, int b)    {
    if (a > b) return a;
    else return b;
}

int imin (int a, int b)    {
    if (a < b) return a;
    else return b;
}

void    *TripleCorr_A(void *threadarg)  {
    int i,j,k,l;
    struct triplecorrstruct *val;
    val = threadarg;
	unsigned int    c0 = 0;
    unsigned int    c1 = 0;
    unsigned int	c2 = 0;
    unsigned short	*data0  = val->array0;
    unsigned short	*data1  = val->array1;
    unsigned short  *data2  = val->array2;
    unsigned int     			p0b0[veclen];
    for(i=0; i<veclen; i++)  	p0b0[i] = 0;
	unsigned int     			p0b1[veclen];
	for(i=0; i<veclen; i++)  	p0b1[i] = 0;
	unsigned int     			p1b0[veclen];
	for(i=0; i<veclen; i++)  	p1b0[i] = 0;
	unsigned int     			p1b1[veclen];
	for(i=0; i<veclen; i++)  	p1b1[i] = 0;
	unsigned int     			p2b2[veclen];
	for(i=0; i<veclen; i++)  	p2b2[i] = 0;
    for(i=0; i < algorquantum; i++) {
        k = i;
        stemp0_A[0] = shatzel0_A[2*n-1];
        stemp1_A[0] = shatzel1_A[2*n-1];
        stemp2_A[0] = shatzel2_A[2*n-1];
        for(l=2*n-1; l>0; l--)  shatzel0_A[l] = shatzel0_A[l-1];
        for(l=2*n-1; l>0; l--)  shatzel1_A[l] = shatzel1_A[l-1];
        for(l=2*n-1; l>0; l--)  shatzel2_A[l] = shatzel2_A[l-1];
        shatzel0_A[0] = data0[i];
        klaus0_A[0]  += data0[i];
        c0 		     += data0[i];
        shatzel1_A[0] = data1[i];
        klaus1_A[0]  += data1[i];
        c1 		     += data1[i]; 
        shatzel2_A[0] = data2[i];
        klaus2_A[0]  += data2[i];
        c2 		     += data2[i];        
        for(j=0; j<2*n; j++)   {
        	p0b0[j] = p0b0[j] + data0[i]*shatzel0_A[j];
        	p0b1[j] = p0b1[j] + data0[i]*shatzel1_A[j];
			p1b0[j] = p1b0[j] + data1[i]*shatzel0_A[j];
        	p1b1[j] = p1b1[j] + data1[i]*shatzel1_A[j];
        	p2b2[j] = p2b2[j] + data2[i]*shatzel2_A[j];
        }
        
        for(j=1; j<pmax; j++)  {
            if((k&1) == 1)  {
                for(l=n; l<2*n; l++)	p0b0[j*n+l] = p0b0[j*n+l] + klaus0_A[j-1]*shatzel0_A[j*n+l];
                for(l=n; l<2*n; l++)	p0b1[j*n+l] = p0b1[j*n+l] + klaus0_A[j-1]*shatzel1_A[j*n+l];
                for(l=n; l<2*n; l++)	p1b0[j*n+l] = p1b0[j*n+l] + klaus1_A[j-1]*shatzel0_A[j*n+l];
                for(l=n; l<2*n; l++)	p1b1[j*n+l] = p1b1[j*n+l] + klaus1_A[j-1]*shatzel1_A[j*n+l];
                for(l=n; l<2*n; l++)	p2b2[j*n+l] = p2b2[j*n+l] + klaus2_A[j-1]*shatzel2_A[j*n+l];
                klaus0_A[j] += klaus0_A[j-1];
                klaus0_A[j-1] = 0;
                klaus1_A[j] += klaus1_A[j-1];
                klaus1_A[j-1] = 0;
                klaus2_A[j] += klaus2_A[j-1];
                klaus2_A[j-1] = 0;
                stemp0_A[j] = shatzel0_A[n*(j+2)-1];
                stemp1_A[j] = shatzel1_A[n*(j+2)-1];
                stemp2_A[j] = shatzel2_A[n*(j+2)-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel0_A[l] = shatzel0_A[l-1];
                shatzel0_A[j*n+n] = shatzel0_A[j*n+n-1] + stemp0_A[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel1_A[l] = shatzel1_A[l-1];
                shatzel1_A[j*n+n] = shatzel1_A[j*n+n-1] + stemp1_A[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel2_A[l] = shatzel2_A[l-1];
                shatzel2_A[j*n+n] = shatzel2_A[j*n+n-1] + stemp2_A[j-1];
            }
		else    j = pmax;
		k = k >> 1;
        } 
    }
	tc0_A = tc0_A + c0;
    tc1_A = tc1_A + c1;
    tc2_A = tc2_A + c2;

    for(i=0; i<veclen; i++)    g0b0[i] = g0b0[i] + (float)p0b0[i];
    for(i=0; i<veclen; i++)    g1b1[i] = g1b1[i] + (float)p1b1[i];
    for(i=0; i<veclen; i++)    g2b2[i] = g2b2[i] + (float)p2b2[i];    
    for(i=0; i<veclen; i++)    g0b1[i] = g0b1[i] + 0.5f*((float)p0b1[i] + (float)p1b0[i]);
    pthread_exit(NULL);
}

void    *TripleCorr_B(void *threadarg)  {
    int i,j,k,l;
    struct triplecorrstruct *val;
    val = threadarg;
	unsigned int    c0 = 0;
    unsigned int    c1 = 0;
    unsigned int	c2 = 0;
    unsigned short	*data0  = val->array0;
    unsigned short	*data1  = val->array1;
    unsigned short  *data2  = val->array2;
    unsigned int     			p0b0[veclen];
    for(i=0; i<veclen; i++)  	p0b0[i] = 0;
	unsigned int     			p0b1[veclen];
	for(i=0; i<veclen; i++)  	p0b1[i] = 0;
	unsigned int     			p1b0[veclen];
	for(i=0; i<veclen; i++)  	p1b0[i] = 0;
	unsigned int     			p1b1[veclen];
	for(i=0; i<veclen; i++)  	p1b1[i] = 0;
	unsigned int     			p2b2[veclen];
	for(i=0; i<veclen; i++)  	p2b2[i] = 0;
    for(i=0; i < algorquantum; i++) {
        k = i;
        stemp0_B[0] = shatzel0_B[2*n-1];
        stemp1_B[0] = shatzel1_B[2*n-1];
        stemp2_B[0] = shatzel2_B[2*n-1];
        for(l=2*n-1; l>0; l--)  shatzel0_B[l] = shatzel0_B[l-1];
        for(l=2*n-1; l>0; l--)  shatzel1_B[l] = shatzel1_B[l-1];
        for(l=2*n-1; l>0; l--)  shatzel2_B[l] = shatzel2_B[l-1];
        shatzel0_B[0] = data0[i];
        klaus0_B[0]  += data0[i];
        c0 		     += data0[i];
        shatzel1_B[0] = data1[i];
        klaus1_B[0]  += data1[i];
        c1 		     += data1[i]; 
        shatzel2_B[0] = data2[i];
        klaus2_B[0]  += data2[i];
        c2 		     += data2[i];        
        for(j=0; j<2*n; j++)   {
        	p0b0[j] = p0b0[j] + data0[i]*shatzel0_B[j];
        	p0b1[j] = p0b1[j] + data0[i]*shatzel1_B[j];
			p1b0[j] = p1b0[j] + data1[i]*shatzel0_B[j];
        	p1b1[j] = p1b1[j] + data1[i]*shatzel1_B[j];
        	p2b2[j] = p2b2[j] + data2[i]*shatzel2_B[j];
        }
        
        for(j=1; j<pmax; j++)  {
            if((k&1) == 1)  {
                for(l=n; l<2*n; l++)	p0b0[j*n+l] = p0b0[j*n+l] + klaus0_B[j-1]*shatzel0_B[j*n+l];
                for(l=n; l<2*n; l++)	p0b1[j*n+l] = p0b1[j*n+l] + klaus0_B[j-1]*shatzel1_B[j*n+l];
                for(l=n; l<2*n; l++)	p1b0[j*n+l] = p1b0[j*n+l] + klaus1_B[j-1]*shatzel0_B[j*n+l];
                for(l=n; l<2*n; l++)	p1b1[j*n+l] = p1b1[j*n+l] + klaus1_B[j-1]*shatzel1_B[j*n+l];
                for(l=n; l<2*n; l++)	p2b2[j*n+l] = p2b2[j*n+l] + klaus2_B[j-1]*shatzel2_B[j*n+l];
                klaus0_B[j] += klaus0_B[j-1];
                klaus0_B[j-1] = 0;
                klaus1_B[j] += klaus1_B[j-1];
                klaus1_B[j-1] = 0;
                klaus2_B[j] += klaus2_B[j-1];
                klaus2_B[j-1] = 0;
                stemp0_B[j] = shatzel0_B[n*(j+2)-1];
                stemp1_B[j] = shatzel1_B[n*(j+2)-1];
                stemp2_B[j] = shatzel2_B[n*(j+2)-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel0_B[l] = shatzel0_B[l-1];
                shatzel0_B[j*n+n] = shatzel0_B[j*n+n-1] + stemp0_B[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel1_B[l] = shatzel1_B[l-1];
                shatzel1_B[j*n+n] = shatzel1_B[j*n+n-1] + stemp1_B[j-1];
                for(l=n*(j+2)-1; l>n*(j+1); l--)   shatzel2_B[l] = shatzel2_B[l-1];
                shatzel2_B[j*n+n] = shatzel2_B[j*n+n-1] + stemp2_B[j-1];
            }
		else    j = pmax;
		k = k >> 1;
        } 
    }
	tc0_B = tc0_B + c0;
    tc1_B = tc1_B + c1;
    tc2_B = tc2_B + c2;

    for(i=0; i<veclen; i++)    g0b0[i] = g0b0[i] + (float)p0b0[i];
    for(i=0; i<veclen; i++)    g1b1[i] = g1b1[i] + (float)p1b1[i];
    for(i=0; i<veclen; i++)    g2b2[i] = g2b2[i] + (float)p2b2[i];    
    for(i=0; i<veclen; i++)    g0b1[i] = g0b1[i] + 0.5f*((float)p0b1[i] + (float)p1b0[i]);
    pthread_exit(NULL);
}

int ShatzelDriver_A(unsigned short *data_0, unsigned short *data_1, unsigned short *data_2, int amt_data, float time_quantum, float *gav0b0, float *gav0b1, float *gav1b1, float *gav2b2, float *gstd0b0, float *gstd0b1, float *gstd1b1, float *gstd2b2, float *counts) {
    int i,j,l;
    
    int corr1,corr2;
    float norm_tau[veclen];
    float gfinal0b0[veclen];
    float gfinal0b1[veclen];
    float gfinal1b1[veclen];
    float gfinal2b2[veclen];

    for(i=0; i<veclen; i++)    norm_tau[i] = 0;
    algorquantum = n * (int) powl(2,pmax);	//2^18 (int = 2^31, ui = 2^32)
    //	ui would overflow after 215s of 20MHz signal.
    int	passes = (floor(amt_data/algorquantum) -1)/2;
    int B_offset = passes*algorquantum;
    array0_0 = malloc (sizeof(unsigned short)*algorquantum);
    array0_1 = malloc (sizeof(unsigned short)*algorquantum);
    array0_2 = malloc (sizeof(unsigned short)*algorquantum);
    array1_0 = malloc (sizeof(unsigned short)*algorquantum);
    array1_1 = malloc (sizeof(unsigned short)*algorquantum);
    array1_2 = malloc (sizeof(unsigned short)*algorquantum);

    for(i=0; i< algorquantum; i++)   array0_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_2[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_2[i] = 0;
    
    for(i=0; i<veclen; i++)    		gfinal0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal2b2[i] = 0.0f;
    
    for(i=0; i<veclen; i++)   		shatzel0_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel1_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel2_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp0_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp1_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp2_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus0_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus1_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus2_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel0_B[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel1_B[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel2_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp0_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp1_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp2_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus0_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus1_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus2_B[i] = 0;
    
    float   *garray0b0;
    float   *garray0b1;
    float   *garray1b1;
    float   *garray2b2;
    int garraysize = passes*veclen;
    garray0b0 = (float *)malloc(garraysize*sizeof(float));
    garray0b1 = (float *)malloc(garraysize*sizeof(float));
    garray1b1 = (float *)malloc(garraysize*sizeof(float));
    garray2b2 = (float *)malloc(garraysize*sizeof(float));
    for(i=0; i<garraysize; i++)    garray0b0[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray0b1[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray1b1[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray2b2[i] = 0.0f;
    struct triplecorrstruct targ1;
    struct triplecorrstruct targ2;
    targ1.array0 = (unsigned short*) array0_0;
    targ1.array1 = (unsigned short*) array0_1;
    targ1.array2 = (unsigned short*) array0_2;
    targ2.array0 = (unsigned short*) array1_0;
    targ2.array1 = (unsigned short*) array1_1;
    targ2.array2 = (unsigned short*) array1_2;

    pthread_attr_t attr;
    pthread_t Acq_thread1;
    pthread_t Acq_thread2;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//  To avoid including partially-weighted values of gIbJ into the std. dev. calc.s, first do a "dummy scan" which increments sumbox & tempxIbJ only, then clear gIbJ(sq)s, then resume calculations.

    for(i=0; i< algorquantum; i++)  {
        array0_0[i] = data_0[i];
        array0_1[i] = data_1[i];
        array0_2[i] = data_2[i];
        array1_0[i] = data_0[B_offset+i];
        array1_1[i] = data_1[B_offset+i];
        array1_2[i] = data_2[B_offset+i];
    }

	corr1 = pthread_create(&Acq_thread1, &attr, TripleCorr_A, (void *) &targ1);
	if (corr1)  { printf("ERROR; return code from pthread_create() is %d\n", corr1); exit(-1); }

	corr2 = pthread_create(&Acq_thread2, &attr, TripleCorr_B, (void *) &targ2);
	if (corr2)  { printf("ERROR; return code from pthread_create() is %d\n", corr2); exit(-1);}

	pthread_join(Acq_thread1, NULL);
	pthread_join(Acq_thread2, NULL);

    int valid_passes;
    for(valid_passes=1; valid_passes<=passes; valid_passes++)   {
        for(i=0; i<veclen; i++)    g0b0[i] = 0.0f;
        for(i=0; i<veclen; i++)    g0b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    g1b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    g2b2[i] = 0.0f;
        
        tc0_A = 0.0f;
        tc1_A = 0.0f;
        tc2_A = 0.0f;
        tc0_B = 0.0f;
        tc1_B = 0.0f;
        tc2_B = 0.0f;
        total_counter = 2*algorquantum;
		for(i=0; i<algorquantum; i++)  {
			array0_0[i] = data_0[valid_passes*algorquantum+i];
			array0_1[i] = data_1[valid_passes*algorquantum+i];
			array0_2[i] = data_2[valid_passes*algorquantum+i];
		}
		corr1 = pthread_create(&Acq_thread1, &attr, TripleCorr_A, (void *) &targ1);
		if (corr1)  { printf("ERROR; return code from pthread_create() is %d\n", corr1); exit(-1); }
		for(i=0; i<algorquantum; i++)  {
			array1_0[i] = data_0[B_offset+valid_passes*algorquantum+i];
			array1_1[i] = data_1[B_offset+valid_passes*algorquantum+i];
			array1_2[i] = data_2[B_offset+valid_passes*algorquantum+i];
		}
		corr2 = pthread_create(&Acq_thread2, &attr, TripleCorr_B, (void *) &targ2);
		if (corr2)  { printf("ERROR; return code from pthread_create() is %d\n", corr2); exit(-1); }
		pthread_join(Acq_thread1, NULL);
		pthread_join(Acq_thread2, NULL);

        tc0 = tc0_A+tc0_B;
        tc1 = tc1_A+tc1_B;
        tc2 = tc2_A+tc2_B;
        counts[0] += tc0;
        counts[1] += tc1;
        counts[2] += tc2;
		for(i=0; i<2*n; i++)   {
			norm_tau[i]  = total_counter;
			gfinal0b0[i] = g0b0[i] * norm_tau[i] / (tc0*tc0);
			gfinal0b1[i] = g0b1[i] * norm_tau[i] / (tc0*tc1);
			gfinal1b1[i] = g1b1[i] * norm_tau[i] / (tc1*tc1);
			gfinal2b2[i] = g2b2[i] * norm_tau[i] / (tc2*tc2);
		}
		for(j=1; j<pmax; j++)  {
			for(l=0; l<n; l++)   {    
				norm_tau[i]  = (float)(total_counter)/(powl(2,j));
				gfinal0b0[i] = g0b0[i] * norm_tau[i] / (tc0*tc0);
				gfinal0b1[i] = g0b1[i] * norm_tau[i] / (tc0*tc1);
				gfinal1b1[i] = g1b1[i] * norm_tau[i] / (tc1*tc1);
				gfinal2b2[i] = g2b2[i] * norm_tau[i] / (tc2*tc2);
				i++;
			}
		}

		for(i=0; i<veclen; i++)    garray0b0[i+(valid_passes-1)*veclen] = gfinal0b0[i];
		for(i=0; i<veclen; i++)    garray0b1[i+(valid_passes-1)*veclen] = gfinal0b1[i];
		for(i=0; i<veclen; i++)    garray1b1[i+(valid_passes-1)*veclen] = gfinal1b1[i];
		for(i=0; i<veclen; i++)    garray2b2[i+(valid_passes-1)*veclen] = gfinal2b2[i];
    }
       
    for(i=0; i<veclen; i++)    gav0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav2b2[i] = 0.0f;

    for(i=0; i<veclen; i++)    gstd0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd2b2[i] = 0.0f;
    
    for(j=0; j<passes; j++)   {
        for(i=0; i<veclen; i++)    gav0b0[i] += garray0b0[i+j*veclen];
        for(i=0; i<veclen; i++)    gav0b1[i] += garray0b1[i+j*veclen];
        for(i=0; i<veclen; i++)    gav1b1[i] += garray1b1[i+j*veclen];
        for(i=0; i<veclen; i++)    gav2b2[i] += garray2b2[i+j*veclen];
    }
    for(i=0; i<veclen; i++)    gav0b0[i] = gav0b0[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav0b1[i] = gav0b1[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav1b1[i] = gav1b1[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav2b2[i] = gav2b2[i] / (float)(passes);
    
    for(j=0; j<passes; j++)   {
        for(i=0; i<veclen; i++)    gstd0b0[i] += pow((garray0b0[i+j*veclen] - gav0b0[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd0b1[i] += pow((garray0b1[i+j*veclen] - gav0b1[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd1b1[i] += pow((garray1b1[i+j*veclen] - gav1b1[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd2b2[i] += pow((garray2b2[i+j*veclen] - gav2b2[i]), 2.0f);
    }
    
    for(i=0; i<veclen; i++)    gstd0b0[i] = pow((gstd0b0[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd0b1[i] = pow((gstd0b1[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd1b1[i] = pow((gstd1b1[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd2b2[i] = pow((gstd2b2[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));

    free(array0_0);
    free(array0_1);
    free(array1_0);
    free(array1_1);
    free(garray0b0);
    free(garray0b1);
    free(garray1b1);
    free(garray2b2);
    return 0;
}

int ShatzelDriver_B(unsigned short *data_0, unsigned short *data_1, unsigned short *data_2, int amt_data, float time_quantum, float *gav0b0, float *gav0b1, float *gav1b1, float *gav2b2, float *gstd0b0, float *gstd0b1, float *gstd1b1, float *gstd2b2, float *counts) {
    int i,j,l;
    
    int corr1,corr2;
    float norm_tau[veclen];
    float gfinal0b0[veclen];
    float gfinal0b1[veclen];
    float gfinal1b1[veclen];
    float gfinal2b2[veclen];

    for(i=0; i<veclen; i++)    norm_tau[i] = 0;
    algorquantum = n * (int) powl(2,pmax);	//2^18 (int = 2^31, ui = 2^32)
    //	ui would overflow after 215s of 20MHz signal.
    int	passes = (floor(amt_data/algorquantum) -1)/2;
    int B_offset = passes*algorquantum;
    array0_0 = malloc (sizeof(unsigned short)*algorquantum);
    array0_1 = malloc (sizeof(unsigned short)*algorquantum);
    array0_2 = malloc (sizeof(unsigned short)*algorquantum);
    array1_0 = malloc (sizeof(unsigned short)*algorquantum);
    array1_1 = malloc (sizeof(unsigned short)*algorquantum);
    array1_2 = malloc (sizeof(unsigned short)*algorquantum);

    for(i=0; i< algorquantum; i++)   array0_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_2[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_2[i] = 0;
    
    for(i=0; i<veclen; i++)    		gfinal0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal2b2[i] = 0.0f;
    
    for(i=0; i<veclen; i++)   		shatzel0_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel1_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel2_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp0_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp1_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp2_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus0_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus1_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus2_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel0_B[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel1_B[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel2_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp0_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp1_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp2_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus0_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus1_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus2_B[i] = 0;
    
    float   *garray0b0;
    float   *garray0b1;
    float   *garray1b1;
    float   *garray2b2;
    int garraysize = passes*veclen;
    garray0b0 = (float *)malloc(garraysize*sizeof(float));
    garray0b1 = (float *)malloc(garraysize*sizeof(float));
    garray1b1 = (float *)malloc(garraysize*sizeof(float));
    garray2b2 = (float *)malloc(garraysize*sizeof(float));
    for(i=0; i<garraysize; i++)    garray0b0[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray0b1[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray1b1[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray2b2[i] = 0.0f;
    struct triplecorrstruct targ1;
    struct triplecorrstruct targ2;
    targ1.array0 = (unsigned short*) array0_0;
    targ1.array1 = (unsigned short*) array0_1;
    targ1.array2 = (unsigned short*) array0_2;
    targ2.array0 = (unsigned short*) array1_0;
    targ2.array1 = (unsigned short*) array1_1;
    targ2.array2 = (unsigned short*) array1_2;

    pthread_attr_t attr;
    pthread_t Acq_thread1;
    pthread_t Acq_thread2;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//  To avoid including partially-weighted values of gIbJ into the std. dev. calc.s, first do a "dummy scan" which increments sumbox & tempxIbJ only, then clear gIbJ(sq)s, then resume calculations.

    for(i=0; i< algorquantum; i++)  {
        array0_0[i] = data_1[i];
        array0_1[i] = data_2[i];
        array0_2[i] = data_0[i];
        array1_0[i] = data_1[B_offset+i];
        array1_1[i] = data_2[B_offset+i];
        array1_2[i] = data_0[B_offset+i];
    }

	corr1 = pthread_create(&Acq_thread1, &attr, TripleCorr_A, (void *) &targ1);
	if (corr1)  { printf("ERROR; return code from pthread_create() is %d\n", corr1); exit(-1); }

	corr2 = pthread_create(&Acq_thread2, &attr, TripleCorr_B, (void *) &targ2);
	if (corr2)  { printf("ERROR; return code from pthread_create() is %d\n", corr2); exit(-1);}

	pthread_join(Acq_thread1, NULL);
	pthread_join(Acq_thread2, NULL);

    int valid_passes;
    for(valid_passes=1; valid_passes<=passes; valid_passes++)   {
        for(i=0; i<veclen; i++)    g0b0[i] = 0.0f;
        for(i=0; i<veclen; i++)    g0b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    g1b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    g2b2[i] = 0.0f;
        
        tc0_A = 0.0f;
        tc1_A = 0.0f;
        tc2_A = 0.0f;
        tc0_B = 0.0f;
        tc1_B = 0.0f;
        tc2_B = 0.0f;
        total_counter = 2*algorquantum;
		for(i=0; i<algorquantum; i++)  {
			array0_0[i] = data_1[valid_passes*algorquantum+i];
			array0_1[i] = data_2[valid_passes*algorquantum+i];
			array0_2[i] = data_0[valid_passes*algorquantum+i];
		}
		corr1 = pthread_create(&Acq_thread1, &attr, TripleCorr_A, (void *) &targ1);
		if (corr1)  { printf("ERROR; return code from pthread_create() is %d\n", corr1); exit(-1); }
		for(i=0; i<algorquantum; i++)  {
			array1_0[i] = data_1[B_offset+valid_passes*algorquantum+i];
			array1_1[i] = data_2[B_offset+valid_passes*algorquantum+i];
			array1_2[i] = data_0[B_offset+valid_passes*algorquantum+i];
		}
		corr2 = pthread_create(&Acq_thread2, &attr, TripleCorr_B, (void *) &targ2);
		if (corr2)  { printf("ERROR; return code from pthread_create() is %d\n", corr2); exit(-1); }
		pthread_join(Acq_thread1, NULL);
		pthread_join(Acq_thread2, NULL);

        tc0 = tc0_A+tc0_B;
        tc1 = tc1_A+tc1_B;
        tc2 = tc2_A+tc2_B;
        counts[1] += tc0;
        counts[2] += tc1;
        counts[0] += tc2;
		for(i=0; i<2*n; i++)   {
			norm_tau[i]  = total_counter;
			gfinal0b0[i] = g0b0[i] * norm_tau[i] / (tc0*tc0);
			gfinal0b1[i] = g0b1[i] * norm_tau[i] / (tc0*tc1);
			gfinal1b1[i] = g1b1[i] * norm_tau[i] / (tc1*tc1);
			gfinal2b2[i] = g2b2[i] * norm_tau[i] / (tc2*tc2);
		}
		for(j=1; j<pmax; j++)  {
			for(l=0; l<n; l++)   {    
				norm_tau[i]  = (float)(total_counter)/(powl(2,j));
				gfinal0b0[i] = g0b0[i] * norm_tau[i] / (tc0*tc0);
				gfinal0b1[i] = g0b1[i] * norm_tau[i] / (tc0*tc1);
				gfinal1b1[i] = g1b1[i] * norm_tau[i] / (tc1*tc1);
				gfinal2b2[i] = g2b2[i] * norm_tau[i] / (tc2*tc2);
				i++;
			}
		}

		for(i=0; i<veclen; i++)    garray0b0[i+(valid_passes-1)*veclen] = gfinal0b0[i];
		for(i=0; i<veclen; i++)    garray0b1[i+(valid_passes-1)*veclen] = gfinal0b1[i];
		for(i=0; i<veclen; i++)    garray1b1[i+(valid_passes-1)*veclen] = gfinal1b1[i];
		for(i=0; i<veclen; i++)    garray2b2[i+(valid_passes-1)*veclen] = gfinal2b2[i];
    }
       
    for(i=0; i<veclen; i++)    gav0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav2b2[i] = 0.0f;

    for(i=0; i<veclen; i++)    gstd0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd2b2[i] = 0.0f;
    
    for(j=0; j<passes; j++)   {
        for(i=0; i<veclen; i++)    gav0b0[i] += garray0b0[i+j*veclen];
        for(i=0; i<veclen; i++)    gav0b1[i] += garray0b1[i+j*veclen];
        for(i=0; i<veclen; i++)    gav1b1[i] += garray1b1[i+j*veclen];
        for(i=0; i<veclen; i++)    gav2b2[i] += garray2b2[i+j*veclen];
    }
    for(i=0; i<veclen; i++)    gav0b0[i] = gav0b0[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav0b1[i] = gav0b1[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav1b1[i] = gav1b1[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav2b2[i] = gav2b2[i] / (float)(passes);
    
    for(j=0; j<passes; j++)   {
        for(i=0; i<veclen; i++)    gstd0b0[i] += pow((garray0b0[i+j*veclen] - gav0b0[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd0b1[i] += pow((garray0b1[i+j*veclen] - gav0b1[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd1b1[i] += pow((garray1b1[i+j*veclen] - gav1b1[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd2b2[i] += pow((garray2b2[i+j*veclen] - gav2b2[i]), 2.0f);
    }
    
    for(i=0; i<veclen; i++)    gstd0b0[i] = pow((gstd0b0[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd0b1[i] = pow((gstd0b1[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd1b1[i] = pow((gstd1b1[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd2b2[i] = pow((gstd2b2[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));

    free(array0_0);
    free(array0_1);
    free(array1_0);
    free(array1_1);
    free(garray0b0);
    free(garray0b1);
    free(garray1b1);
    free(garray2b2);
    return 0;
}
int ShatzelDriver_C(unsigned short *data_0, unsigned short *data_1, unsigned short *data_2, int amt_data, float time_quantum, float *gav0b0, float *gav0b1, float *gav1b1, float *gav2b2, float *gstd0b0, float *gstd0b1, float *gstd1b1, float *gstd2b2, float *counts) {
    int i,j,l;
    
    int corr1,corr2;
    float norm_tau[veclen];
    float gfinal0b0[veclen];
    float gfinal0b1[veclen];
    float gfinal1b1[veclen];
    float gfinal2b2[veclen];

    for(i=0; i<veclen; i++)    norm_tau[i] = 0;
    algorquantum = n * (int) powl(2,pmax);	//2^18 (int = 2^31, ui = 2^32)
    //	ui would overflow after 215s of 20MHz signal.
    int	passes = (floor(amt_data/algorquantum) -1)/2;
    int B_offset = passes*algorquantum;
    array0_0 = malloc (sizeof(unsigned short)*algorquantum);
    array0_1 = malloc (sizeof(unsigned short)*algorquantum);
    array0_2 = malloc (sizeof(unsigned short)*algorquantum);
    array1_0 = malloc (sizeof(unsigned short)*algorquantum);
    array1_1 = malloc (sizeof(unsigned short)*algorquantum);
    array1_2 = malloc (sizeof(unsigned short)*algorquantum);

    for(i=0; i< algorquantum; i++)   array0_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array0_2[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_0[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_1[i] = 0;
    for(i=0; i< algorquantum; i++)   array1_2[i] = 0;
    
    for(i=0; i<veclen; i++)    		gfinal0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    		gfinal2b2[i] = 0.0f;
    
    for(i=0; i<veclen; i++)   		shatzel0_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel1_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel2_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp0_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp1_A[i] = 0;
    for(i=0; i<pmax; i++)           stemp2_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus0_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus1_A[i] = 0;
    for(i=0; i<pmax; i++)           klaus2_A[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel0_B[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel1_B[i] = 0;
    for(i=0; i<veclen; i++)   		shatzel2_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp0_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp1_B[i] = 0;
    for(i=0; i<pmax; i++)           stemp2_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus0_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus1_B[i] = 0;
    for(i=0; i<pmax; i++)           klaus2_B[i] = 0;
    
    float   *garray0b0;
    float   *garray0b1;
    float   *garray1b1;
    float   *garray2b2;
    int garraysize = passes*veclen;
    garray0b0 = (float *)malloc(garraysize*sizeof(float));
    garray0b1 = (float *)malloc(garraysize*sizeof(float));
    garray1b1 = (float *)malloc(garraysize*sizeof(float));
    garray2b2 = (float *)malloc(garraysize*sizeof(float));
    for(i=0; i<garraysize; i++)    garray0b0[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray0b1[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray1b1[i] = 0.0f;
    for(i=0; i<garraysize; i++)    garray2b2[i] = 0.0f;
    struct triplecorrstruct targ1;
    struct triplecorrstruct targ2;
    targ1.array0 = (unsigned short*) array0_0;
    targ1.array1 = (unsigned short*) array0_1;
    targ1.array2 = (unsigned short*) array0_2;
    targ2.array0 = (unsigned short*) array1_0;
    targ2.array1 = (unsigned short*) array1_1;
    targ2.array2 = (unsigned short*) array1_2;

    pthread_attr_t attr;
    pthread_t Acq_thread1;
    pthread_t Acq_thread2;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//  To avoid including partially-weighted values of gIbJ into the std. dev. calc.s, first do a "dummy scan" which increments sumbox & tempxIbJ only, then clear gIbJ(sq)s, then resume calculations.

    for(i=0; i< algorquantum; i++)  {
        array0_0[i] = data_2[i];
        array0_1[i] = data_0[i];
        array0_2[i] = data_1[i];
        array1_0[i] = data_2[B_offset+i];
        array1_1[i] = data_0[B_offset+i];
        array1_2[i] = data_1[B_offset+i];
    }

	corr1 = pthread_create(&Acq_thread1, &attr, TripleCorr_A, (void *) &targ1);
	if (corr1)  { printf("ERROR; return code from pthread_create() is %d\n", corr1); exit(-1); }

	corr2 = pthread_create(&Acq_thread2, &attr, TripleCorr_B, (void *) &targ2);
	if (corr2)  { printf("ERROR; return code from pthread_create() is %d\n", corr2); exit(-1);}

	pthread_join(Acq_thread1, NULL);
	pthread_join(Acq_thread2, NULL);

    int valid_passes;
    for(valid_passes=1; valid_passes<=passes; valid_passes++)   {
        for(i=0; i<veclen; i++)    g0b0[i] = 0.0f;
        for(i=0; i<veclen; i++)    g0b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    g1b1[i] = 0.0f;
        for(i=0; i<veclen; i++)    g2b2[i] = 0.0f;
        
        tc0_A = 0.0f;
        tc1_A = 0.0f;
        tc2_A = 0.0f;
        tc0_B = 0.0f;
        tc1_B = 0.0f;
        tc2_B = 0.0f;
        total_counter = 2*algorquantum;
		for(i=0; i<algorquantum; i++)  {
			array0_0[i] = data_2[valid_passes*algorquantum+i];
			array0_1[i] = data_0[valid_passes*algorquantum+i];
			array0_2[i] = data_1[valid_passes*algorquantum+i];
		}
		corr1 = pthread_create(&Acq_thread1, &attr, TripleCorr_A, (void *) &targ1);
		if (corr1)  { printf("ERROR; return code from pthread_create() is %d\n", corr1); exit(-1); }
		for(i=0; i<algorquantum; i++)  {
			array1_0[i] = data_2[B_offset+valid_passes*algorquantum+i];
			array1_1[i] = data_0[B_offset+valid_passes*algorquantum+i];
			array1_2[i] = data_1[B_offset+valid_passes*algorquantum+i];
		}
		corr2 = pthread_create(&Acq_thread2, &attr, TripleCorr_B, (void *) &targ2);
		if (corr2)  { printf("ERROR; return code from pthread_create() is %d\n", corr2); exit(-1); }
		pthread_join(Acq_thread1, NULL);
		pthread_join(Acq_thread2, NULL);

        tc0 = tc0_A+tc0_B;
        tc1 = tc1_A+tc1_B;
        tc2 = tc2_A+tc2_B;
        counts[2] += tc0;
        counts[0] += tc1;
        counts[1] += tc2;
		for(i=0; i<2*n; i++)   {
			norm_tau[i]  = total_counter;
			gfinal0b0[i] = g0b0[i] * norm_tau[i] / (tc0*tc0);
			gfinal0b1[i] = g0b1[i] * norm_tau[i] / (tc0*tc1);
			gfinal1b1[i] = g1b1[i] * norm_tau[i] / (tc1*tc1);
			gfinal2b2[i] = g2b2[i] * norm_tau[i] / (tc2*tc2);
		}
		for(j=1; j<pmax; j++)  {
			for(l=0; l<n; l++)   {    
				norm_tau[i]  = (float)(total_counter)/(powl(2,j));
				gfinal0b0[i] = g0b0[i] * norm_tau[i] / (tc0*tc0);
				gfinal0b1[i] = g0b1[i] * norm_tau[i] / (tc0*tc1);
				gfinal1b1[i] = g1b1[i] * norm_tau[i] / (tc1*tc1);
				gfinal2b2[i] = g2b2[i] * norm_tau[i] / (tc2*tc2);
				i++;
			}
		}

		for(i=0; i<veclen; i++)    garray0b0[i+(valid_passes-1)*veclen] = gfinal0b0[i];
		for(i=0; i<veclen; i++)    garray0b1[i+(valid_passes-1)*veclen] = gfinal0b1[i];
		for(i=0; i<veclen; i++)    garray1b1[i+(valid_passes-1)*veclen] = gfinal1b1[i];
		for(i=0; i<veclen; i++)    garray2b2[i+(valid_passes-1)*veclen] = gfinal2b2[i];
    }
       
    for(i=0; i<veclen; i++)    gav0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gav2b2[i] = 0.0f;

    for(i=0; i<veclen; i++)    gstd0b0[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd0b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd1b1[i] = 0.0f;
    for(i=0; i<veclen; i++)    gstd2b2[i] = 0.0f;
    
    for(j=0; j<passes; j++)   {
        for(i=0; i<veclen; i++)    gav0b0[i] += garray0b0[i+j*veclen];
        for(i=0; i<veclen; i++)    gav0b1[i] += garray0b1[i+j*veclen];
        for(i=0; i<veclen; i++)    gav1b1[i] += garray1b1[i+j*veclen];
        for(i=0; i<veclen; i++)    gav2b2[i] += garray2b2[i+j*veclen];
    }
    for(i=0; i<veclen; i++)    gav0b0[i] = gav0b0[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav0b1[i] = gav0b1[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav1b1[i] = gav1b1[i] / (float)(passes);
    for(i=0; i<veclen; i++)    gav2b2[i] = gav2b2[i] / (float)(passes);
    
    for(j=0; j<passes; j++)   {
        for(i=0; i<veclen; i++)    gstd0b0[i] += pow((garray0b0[i+j*veclen] - gav0b0[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd0b1[i] += pow((garray0b1[i+j*veclen] - gav0b1[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd1b1[i] += pow((garray1b1[i+j*veclen] - gav1b1[i]), 2.0f);
        for(i=0; i<veclen; i++)    gstd2b2[i] += pow((garray2b2[i+j*veclen] - gav2b2[i]), 2.0f);
    }
    
    for(i=0; i<veclen; i++)    gstd0b0[i] = pow((gstd0b0[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd0b1[i] = pow((gstd0b1[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd1b1[i] = pow((gstd1b1[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));
    for(i=0; i<veclen; i++)    gstd2b2[i] = pow((gstd2b2[i]/(float)(passes)), 0.5f)/(pow((float)(passes),0.5));

    free(array0_0);
    free(array0_1);
    free(array1_0);
    free(array1_1);
    free(garray0b0);
    free(garray0b1);
    free(garray1b1);
    free(garray2b2);
    printf("Program completed.\n");
    return 0;
}
