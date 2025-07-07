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





/*
 * main
 *
 * This function serves as the entry point for the program. It orchestrates the execution of data analysis on multiple data sets.
 * The function is responsible for setting up the data set paths, determining the number of data sets, and invoking the analysis
 * functions to process and analyze each data set. It also contains commented-out sections for single data set execution and
 * testing unit extraction, which can be used for further development or testing purposes.
 */
int main(int argc, const char * argv[])
{
	///
	/*-----------   Single Data Set Execution   -----------* /
	 ///
	 ///
	 ////*
	 // /*-----------   User-Provide Pathname(Hardcoded for now)   -----------*/
	const char *dataSetFilePathName = "/Users/98dav/Desktop/Xcode/C-Programs/Automated_CSV_Data_Analysis/physics_particles.txt";//"/Users/98dav/Desktop/Xcode/C-Programs/CSV_Data_Set_Analysis/physics_particles.txt"; // Pathname to data set
	


	
	//  /*-----------   Run Data Set   -----------* /
	DataSetAnalysis dataAnalysis = process_data_set_for_analysis(dataSetFilePathName);
	const char *preprocessedDataDirectory = dataAnalysis.plottableDataDirectory.directoryPathName;
	
	
	// Perform statistical analysis on the preprocessed data
	//print_written_data_set(preprocessedDataDirectory);
	print_directory_propertiees(dataAnalysis.plottableDataDirectory);
	
	
	
	const char *analysisDir = perform_full_analysis_and_modeling(preprocessedDataDirectory);
	DirectoryProperties analysisDirectoryProperties = extract_directory_properties(analysisDir);
	print_directory_propertiees(analysisDirectoryProperties);
	//print_written_data_set(analysisDir);
	//*/
	
	
	
	
	
	
	
	/// TESTING UNIT EXTRACTION
	/*
	 /*-----------   Begin Preprocessing File Contents to Standardize the Format and Achieve/Maintain Compatibility of the Contents   -----------* /
	 int fieldCount = count_data_fields(fileContents[0]);
	 char **formattedFileContents = fileContents;
	 char *headerLine = prune_and_trim_problematic_characters_from_string(fileContents[1], delimiter, fieldCount);
	 for(int i = 1; i < lineCount; i++)
	 {
	 formattedFileContents[i] = prune_and_trim_problematic_characters_from_string(fileContents[i], delimiter, fieldCount);
	 }
	 /// Individually preprocess, format, and capture each line of the data set in order to achieve a fully preprocessed/formatted data set
	 for(int i = 1; i < lineCount; i++)
	 {
	 /// Examine each data entry and filter out problematic characters, omit or replace disruptive aspects(such as repeated delimiters, whitespaces,
	 /// unexpected date/time formats, etc.), and standardize the arrangement of characters representing information for each data entry
	 char *dataEntry =  prune_and_trim_problematic_characters_from_string(fileContents[i], delimiter, fieldCount);
	 /// Pass the now preprocessed 'dataEntry' to the utility function, 'format_data_entry_for_plotting', where it will be formatted to ensure each one of it's
	 /// field's is of the correct data type, and correctly handling the case(s) when it they are not(i.e., if the field missing, mismatched data type, etc.), in order to
	 /// ensure compatibility with the rest of the program before being captured as an element of the 'formattedFileContents'
	 formattedFileContents[i] = format_data_entry_for_plotting(headerLine, dataEntry, fieldCount, delimiter);
	 }
	 char **extractedUnits = extract_units_from_fields(formattedFileContents[1], delimiter, fieldCount);
	 //char *extractedUnitsString = concatenate_string_array(extractedUnits, fieldCount, delimiter);
	 //printf("\n Extracted Units: %s", extractedUnitsString);
	 //print_string_array(extractedUnits, fieldCount, "extract_units_from_fields");
	 //*/
	
	return 0;
}


