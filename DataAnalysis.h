//  DataAnalysis.h
//  Automated_CSV_Data_Analysis
//  DavidRichardson02


#ifndef DataAnalysis_h
#define DataAnalysis_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "FileUtilities.h"




/**
 * DataSetProperties Structure: Encapsulates the key static properties of a data set for analysis and processing.
 *
 * The DataSetProperties structure provides a convenient way to store and pass around the key properties of a data set
 * once they have already been rigorously determined and compatibility with plotting tools has been ensured.
 *
 * Struct for data set properties members:
 *      - int entryCount: The number of data entries in the data set.
 *      - int fieldCount: The number of fields per data entry.
 *      - const char *delimiter: The delimiter character used to separate fields in the data set.
 *      - char *dataSetHeader: Header line of the data set.
 *      - char **fieldNameTypePairs: Array of strings storing pairs of field names and their corresponding types.
 *      - const char* dataSetFilePathName: Path to the data set file.
 *      - int *missingDataCount: Array of integers storing the count of missing values for each data entry.
 *      - int *plottabilityStatus: Array of integers indicating the plottability status of each field(1 for plottable and 0 for unplottable).
 *      - char** commonDataTypes: Array of strings where each string represents the most common data type of a field.
 *      - char **dataSetFileContents: Array of strings containing the entirety of contents of the data set file, the first pass data capture(everything in the original data set).
 */
typedef struct
{
	int entryCount;
	int fieldCount;
	
	
	const char *delimiter;
	char *dataSetFieldNames;
	char **fieldNameTypePairs;
	const char* dataSetFilePathName;
	
	
	int *missingDataCount;
	int *plottabilityStatus;
	char** commonDataTypes;
	
	char **dataSetFileContents;
} DataSetProperties;
DataSetProperties analyze_data_set_properties(const char *filePathName); // Function to analyze the properties of a data set for processing and analysis


/**
 * DataSetAnalysis Structure: Represents the properties of the data set that are required for data
 * analysis and operations, but not necessarily important for the actual meaning of data. For example, the
 * delimiter a data set uses has basically nothing to do with the contents of the data or their meaning, but
 * analyzing the data in any meaningful way would be impossible without a way to differentiate between various
 * fields and entries.
 *
 * This structure encapsulates the various attributes of a data set's properties that are relevant for
 * the capturing, formatting, and analysis of data, without which, there would be no way to perform
 * operations on data(NOTE: the properties are what's required, this encapsulation structure is just a convenient/efficient way to group them together).
 *
 * Each member variable is associated with a specific property of data set analysis, where each
 * property facilitates operations on the data, such as reading, writing, extraction of values, parsing, processing,
 * statistical operations, interpreting, etc.
 *
 *
 * Struct for data set analysis properties:
 *
 * 		- DataSetProperties dataSetProperties: The key properties of the data set required for analysis and processing, comprised of:
 *      		- int entryCount: The number of data entries in the data set.
 *      		- int fieldCount: The number of fields per data entry.
 *     		 	- const char *delimiter: The delimiter character used to separate fields in the data set.
 *     			- char *dataSetHeader: Header line of the data set.
 *      		- char **fieldNameTypePairs: Array of strings storing pairs of field names and their corresponding types.
 *      		- const char* dataSetFilePathName: Path to the data set file.
 *      		- int *missingDataCount: Count of missing values per field.
 *      		- int *plottabilityStatus: Array of integers indicating the plottability status of each field(1 for plottable and 0 for unplottable).
 *      		- char** commonDataTypes: Array of strings where each string represents the most common data type of a field.
 *      		- char **dataSetFileContents: Array of strings containing the entirety of contents of the data set file, the first pass data capture(everything in the original data set).
 *
 *      - double **dataSetRadixSortedContents: Radix sorted contents of the data set.
 *      - FileDirectory plottableDataDirectory: Directory for plottable data fields.
 *      - FileDirectory processedDataDirectory: Directory for processed data.
 *      - FileDirectory analysisDataDirectory: Directory for analysis data.
 */
typedef struct
{
	DataSetProperties *dataSetProperties;
	double **radixSortedData;
	
	
	DirectoryProperties plottableDataDirectory;
	DirectoryProperties analysisDataDirectory;
} DataSetAnalysis;






DataSetAnalysis configure_data_set_analysis(DataSetProperties dataSetProperties, const char *preprocessedDataDirectory); // Function to configure the analysis of a data set for processing and operations
char **blindly_extract_data_set(const char* dataSetFilePathName, int lineCount); // Function to extract the entire contents of a data set without any formatting or processing
char **extract_and_format_data_set(char **fileContents, int lineCount, int fieldCount, const char *delimiter); // Function to extract and format the contents of a data set for processing and analysis
DataSetAnalysis process_data_set_for_analysis(const char* dataSetFilePathName); // Function to process a data set for analysis, including extraction, formatting, and preparation for statistical operations
DataSetAnalysis analyze_preprocessed_data_set(DataSetAnalysis ); // Function to analyze a preprocessed data set for further operations and statistical analysis
const char *perform_full_analysis_and_modeling(const char *preprocessedDataDirectoryPath); // Function to perform full analysis and modeling on the preprocessed data set, generating results and output directories
















#endif /* DataAnalysis_h */
