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




void gaussian_erf_plot(const char* filePathName);
// ------------- Helper Functions for Creating MATLAB(.m) Scripts/Files to Plot Some Model(various statistical models of a file's(.csv) data set, including but not limited to: histogram, theoretical distribution, standard normal distribution, etc.) of Data -------------
/// \{

void generate_individual_matlab_scripts(const char *resultsDir, char **fieldNames, int numFields);
void generate_comprehensive_matlab_script(const char *analysisDir, char **fieldNames, int numFields);

/// \}








void gaussian_error_function_plot(const char* filePathName);
void generate_fieldwise_modeling_scripts(const char *analysisDir, char **fieldNames, int numFields);
void generate_comprehensive_modeling_script(const char *analysisDir, char **fieldNames, int numFields);









void generate_modeling_matlab_scripts(const char *modelingDir, char **fieldNames, int numFields);
void generate_modeling_matlab_script_all_fields(const char *modelingDir, char **fieldNames, int numFields);
/**
 // Suppose 'resultsDirectory' is your top-level results from analysis
 char *modelingDir = combine_strings(resultsDirectory, "/Modeling");
 create_directory(modelingDir, "");
 
 // We assume you have a list of field names & the number of fields:
 generate_modeling_matlab_scripts(modelingDir, fieldNames, numFields);
 generate_modeling_matlab_script_all_fields(modelingDir, fieldNames, numFields);
 
 free(modelingDir);
 */
















#endif /* DataSetModeling_h */
