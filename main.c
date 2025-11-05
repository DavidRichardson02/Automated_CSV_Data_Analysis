//  main.c
//  Automated_CSV_Data_Analysis
//  DavidRichardson02

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "CommonDefinitions.h"
#include "GeneralUtilities.h"
#include "StringUtilities.h"
#include "FileUtilities.h"
#include "DataExtraction.h"
#include "Integrators.h"
#include "StatisticalMethods.h"
#include "DataSetModeling.h"
#include "DataAnalysis.h"
#include "DebuggingUtilities.h"




/**
 * run_directory_of_multiple_data_sets(const char *dataSetDirectoryPathName)
 *
 * PURPOSE:
 *   • Accepts a path to a directory containing multiple data sets.
 *   • Iterates over each data file to execute data extraction, preprocessing,
 *     statistical analysis, and modeling.
 *   • Returns an array of strings indicating the analysis results' directory paths.
 *
 * KEY IDEAS:
 *   • Each data file path is stored and processed through a high-level analysis routine.
 *   • The resulting directories from the analysis phase are collected into
 *     a return array, which allows for supplementary post-processing.
 *   • Memory is allocated dynamically to store file paths; the function
 *     expects the external utilities to handle safe memory usage and release.
 */
const char **run_directory_of_multiple_data_sets(const char *dataSetDirectoryPathName)
{
	
	/// Step 1: Extract directory properties from user-supplied path
	DirectoryProperties dataSetDirectoryProperties = extract_directory_properties(dataSetDirectoryPathName); 	// Extracts the file count and associated property details of each file
	
	
	/// Step 2: Allocate arrays to store file paths for data sets and analysis results
	const char **dataSetFilePathNames = allocate_memory_char_ptr_ptr(dataSetDirectoryProperties.fileCount + 1, MAX_STRING_SIZE); // +1 for safety/null termination
	const char **analysisDirectoryPathNames = allocate_memory_char_ptr_ptr(dataSetDirectoryProperties.fileCount + 1, MAX_STRING_SIZE); // +1 for safety/null termination
	
	
	/// Step 3: Loop through each data file to perform comprehensive analysis
	// Iterates over the file count, extracting each file path and performing analysis
	for(int i = 0; i <= dataSetDirectoryProperties.fileCount; i++)
	{
		dataSetFilePathNames[i] = dataSetDirectoryProperties.fileProperties[i].filePathName; // Capture and store the pathname for the current data file
		
		
		
		/// Step 4: Process and analyze the current data set
		DataSetAnalysis dataAnalysis = process_data_set_for_analysis(dataSetFilePathNames[i]); // Produces a DataSetAnalysis object, which holds references to preprocessed data
		const char *preprocessedDataDirectory = dataAnalysis.plottableDataDirectory.directoryPathName; // Obtain the directory path for data prepared for plotting
		print_directory_propertiees(dataAnalysis.plottableDataDirectory);
		
		
		
		/// Step 5: Execute detailed data analysis and modeling on the preprocessed data
		const char *analysisDir = perform_full_analysis_and_modeling(preprocessedDataDirectory);
		DirectoryProperties analysisDirectoryProperties = extract_directory_properties(analysisDir); // Fetch properties of the analysis output directory and display them
		print_directory_propertiees(analysisDirectoryProperties);
		analysisDirectoryPathNames[i] = analysisDir; // Store the analysis directory path for future reference
	}
	

	/// Step 6: Return an array of analysis directory paths for further use
	return analysisDirectoryPathNames;
}












/**
 * main
 *
 * The 'main' function provides a demonstration and entry point for processing
 * a single data set, including cleaning, preprocessing, and performing a
 * comprehensive analysis to produce results for further usage and interpretation.
 *
 * KEY POINTS:
 * 1) The user can specify the data set file path (currently hardcoded for demonstration).
 * 2) The function orchestrates the execution of data analysis routines:
 *    - 'process_data_set_for_analysis' for data cleaning/preprocessing,
 *    - 'perform_full_analysis_and_modeling' for deeper statistical analysis
 *      and modeling.
 * 3) The function outputs directory properties relevant to intermediate and
 *    final results for easy verification or logging.
 */
int main(int argc, const char * argv[])
{
	/*-----------   Step 1: Single Data Set Execution   -----------*/
	/// User-Provided pathname to a specific data set for demonstration
	const char *dataSetFilePathName = "/Users/98dav/Desktop/Interview Preparation/Datasets/multiemployerlist.csv";//"/Users/98dav/Desktop/Xcode/C-Programs/Automated_CSV_Data_Analysis/weather_measurements.txt";
	//"/Users/98dav/Desktop/Xcode/C-Programs/Automated_CSV_Data_Analysis/beach-weather-stations-automated-sensors-1.csv"; // The pathname to data set is hardcoded as 'dataSetFilePathName'
	
	
	
	
	/*-----------   Step 2: Run Data Set to Clean, Preprocess, and Identify+Capture Plottable Data  -----------*/
	/// Perform data set processing, which involves cleaning and identifying plottable data
	/// 'process_data_set_for_analysis' returns an object containing references
	/// to the preprocessed data for further use
	DataSetAnalysis dataAnalysis = process_data_set_for_analysis(dataSetFilePathName);
	const char *preprocessedDataDirectory = dataAnalysis.plottableDataDirectory.directoryPathName;
	print_directory_propertiees(dataAnalysis.plottableDataDirectory);
	
	
	
	
	/*-----------   Step 3: Execute statistical analysis and modeling procedures on the preprocessed data   -----------*/
	/// 'perform_full_analysis_and_modeling' applies advanced operations to derive insights and generate output directories for results
	const char *analysisDir = perform_full_analysis_and_modeling(preprocessedDataDirectory);
	DirectoryProperties analysisDirectoryProperties = extract_directory_properties(analysisDir);
	print_directory_propertiees(analysisDirectoryProperties);
	return 0;
}


