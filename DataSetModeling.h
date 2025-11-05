//  DataSetModeling.h
//  Automated_CSV_Data_Analysis
//  DavidRichardson02


#ifndef DataSetModeling_h
#define DataSetModeling_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "Integrators.h"
#include "StatisticalMethods.h"








/* ========================================================================
 * Unified MATLAB-Script Generation
 * ===================================================================== */
typedef enum MatlabScriptFlavor {
	MATLAB_SCRIPTS_INDIVIDUAL      = 1,  // one <field>_plot.m per field
	MATLAB_SCRIPTS_COMPREHENSIVE   = 2,  // one comprehensive_plots.m
	MATLAB_SCRIPTS_BOTH            = 3   // both of the above
} MatlabScriptFlavor;

/**
 * generate_matlab_scripts_unified
 * analysisDir : directory where *_histogram.txt and *_stats/_analysis files live
 * fieldNames  : array of field-name strings (length = numFields)
 * numFields   : number of fields
 * flavor      : INDIVIDUAL / COMPREHENSIVE / BOTH
 * scriptsSubdir (optional): if non-NULL/non-empty, scripts go to analysisDir/<scriptsSubdir>/ ;
 *                           otherwise they go directly in analysisDir.
 */
void generate_matlab_scripts_unified(const char *analysisDir, char **fieldNames, int numFields, MatlabScriptFlavor flavor, const char *scriptsSubdir);




// ------------- Helper Functions for Creating MATLAB(.m) Scripts/Files to Plot Some Model(various statistical models of a file's(.csv) data set, including but not limited to: histogram, theoretical distribution, standard normal distribution, etc.) of Data -------------
/// \{
void generate_matlab_model_scripts_unified(const char *analysisDir, char **fieldNames, int numFields, int overwrite);
/// \}



















#endif /* DataSetModeling_h */
