//  DataAnalysis.c
//  Automated_CSV_Data_Analysis
//  DavidRichardson02


#include "DataAnalysis.h"
#include "GeneralUtilities.h"
#include "CommonDefinitions.h"
#include "StringUtilities.h"
#include "FileUtilities.h"
#include "Integrators.h"
#include "StatisticalMethods.h"
#include "DataExtraction.h"
#include "DataSetModeling.h"
#include "DebuggingUtilities.h"








/**
 The same basic steps apply regardless of the strategy of implementation, which even still, is only relevant for data preprocessing
 and extraction/parsing, this is because once the data has been sent to the files the remaining steps are relatively straightforward
 in comparison, and then the actual first step, just blindly capturing the data from the file, is trivial.
 
 The only strategy that is even somewhat tennative to change is the one for how data preprocessing is implemented, the current favorite
 
 
 Even then, the strategy of implementation is only relevant u the data has been captured, preprocessed, and extracted/parsed, then the
 I. CAPTURE DATA SET CONTENTS
 
 
 standardize data preprocessing to ensure compatibility throughout the entire program
 
 II. PREPROCESS DATA SET CONTENTS TO ENSURE COMPATIBILITY THROUGHOUT THE ENTIRE PROGRAM
 1. Capture the contents of the file in an array of strings
 2. Format the string array contents to be compatible with my CSV file operations to help standardize file's contents for further use and prune the array of strings to remove unwelcome things(things like: double commas,, lone - or + signs, etc.), ensure the file lines are organized in the order of header line with names of data fields followed by all data entries on subsequent lines. The header line is commonly a required prerequisite(depends on file contents/data set) and while data can still be plotted without any provided fields of data entries, it will likely lack key information, such as units, name, magnitude, categorization, etc.
 3. Parse the file contents to independently capture all entries for each of the data set's fields, so the name of each field, followed by all of it's data entries, will be captured. For example a data set on particles with the charge, mass, and name will have a main array of strings to capture all contents, which will then be used to create new files to capture the contents  and string/double/int/float/etc. array for each(along with the name of the field, which will be the first element in the field array)
 4. Categorize the data paramaters into plottable and unplottable values(ex: string is unplottable and a double value is plottable) and create appropriately named and located directories(based on og file pathname) for these new and place them into. The primary purpose of this step is to isolate the numerical data fields to allow the application of various analytical methods for plotting things like histograms and theoretical distributions. Additionally, these categorizations will be recorded and indexed to be used later to match things like the title of a graph to the data plotted on it.
 
 
 The entries of the data set will all be defined by the same number of fields(both numerical and character based fields are could possibly be NULL)
 The number of fields defining any given data entry, as read from the header line of the csv data set file
 
 Capture and categorize the contents of .csv data sets, continually refining, streamlining, and formatting data to prepare for plotting. Methodically
 
 
 
 
 II. DATA ANALYSIS
 */







/**
 * analyze_data_set_properties
 *
 *   This function orchestrates the crucial steps needed to convert a raw CSV-like file
 *   into a cohesive internal representation, capturing essential properties about the data set.
 *   It tackles file loading, determining how the data is delimited, formatting the content,
 *   classifying which fields are missing or plottable, and packaging these findings into
 *   a DataSetProperties structure for convenient downstream usage.
 *
 * Purpose & Role in Implementation:
 *   1) Read the file and determine the delimiter to standardize field parsing.
 *   2) Count how many lines and fields exist, then format the data for consistency.
 *   3) Detect each field’s plottability (is it numeric?) and track missing entries.
 *   4) Combine, store, and return all relevant metadata (field count, missing data counts,
 *      data set file path, etc.) in a unified DataSetProperties struct.
 *
 * Reasoning & Key Points:
 *   - By centralizing file reading, data cleaning, and field classification here, subsequent
 *     analysis steps can operate on well-structured data without re-checking basic assumptions.
 *   - The function returns DataSetProperties, a structure that captures both global file
 *     properties (e.g., total lines, delimiter) and field-specific details (e.g., numeric vs.
 *     nonnumeric, number of missing entries).
 */
DataSetProperties analyze_data_set_properties(const char *filePathName)
{
	int lineCount = count_file_lines(filePathName, MAX_NUM_FILE_LINES); // lineCount holds the total number of lines in the file (upper bounded by MAX_NUM_FILE_LINES).
	char **fileContents = parse_file_contents(filePathName, lineCount); // fileContents is a char** array that stores each line from the file.
	const char *delimiter = identify_delimiter(fileContents, lineCount); // delimiter is determined using a heuristic that checks the file contents.
	int fieldCount = count_data_fields(fileContents[0]); // fieldCount is derived from the file's header line (e.g., how many fields it contains).
	char **formattedContents = extract_and_format_data_set(fileContents, lineCount, fieldCount, delimiter); // Ensure all rows/fieds are standardized
	char **nameTypePairs = capture_data_set_header_for_plotting(formattedContents[0], formattedContents, delimiter); // Each field's name is paired with an inferred type, forming a crucial metadata for later steps.
	
	
	int *missingDataCount = count_missing_values(formattedContents, lineCount, fieldCount, delimiter, formattedContents[1]);; // Identify missing values, missingDataCount tracks the number of missing entries per field.
	int *plottabilityStatus = identify_plottable_fields(nameTypePairs, fieldCount, ";"); // A simple scheme that sets an entry to 1 for numeric fields (plottable) and 0 for non-numeric (unplottable).
	//char** commonDataTypes = determine_common_data_entry_types(formattedContents, lineCount, fieldCount, delimiter);
	char **fieldNamesTokenized = split_tokenized_string(formattedContents[0], delimiter, fieldCount); // Tokenize field names and accumulate plottable field names, a helper array to separate each field name
	char *plottableFieldNames;
	for(int i = 0; i < fieldCount; i++)
	{
		char *fieldNames;
		if(plottabilityStatus[i] == 0)
		{
			fieldCount--; // If a field is unplottable, reduce the fieldCount accordingly (though be mindful of side effects).
			
		}
		else if(plottabilityStatus[i] == 1)
		{
			fieldNames = combine_strings(fieldNames, fieldNamesTokenized[i]); // Combine each plottable field name into a single string.
		}
		else
		{
			perror("\n\nError: Unrecognized plottability status in 'analyze_data_set_properties'.\n");
			exit(1);
		}
		
		if(i == fieldCount)
		{
			plottableFieldNames = fieldNames; // The final iteration captures the entire concatenated field names.
		}
		
	}

	
	
	
	/// Populate the DataSetProperties structure and return it
	DataSetProperties dataSetProperties;
	dataSetProperties.entryCount = lineCount;
	dataSetProperties.fieldCount = fieldCount;
	dataSetProperties.delimiter = delimiter;
	dataSetProperties.dataSetFieldNames = formattedContents[0];
	dataSetProperties.fieldNameTypePairs =  nameTypePairs;
	dataSetProperties.dataSetFilePathName = filePathName;
	
	dataSetProperties.missingDataCount = missingDataCount;
	dataSetProperties.plottabilityStatus = plottabilityStatus;
	dataSetProperties.commonDataTypes = determine_common_data_entry_types(formattedContents, lineCount, fieldCount, delimiter);;//commonDataTypes;
	dataSetProperties.dataSetFileContents = formattedContents;
	

	
	return dataSetProperties; // Return the fully assembled structure containing all essential dataset properties.
}



DataSetAnalysis configure_data_set_analysis(DataSetProperties dataSetProperties, const char *preprocessedDataDirectory)
{
	double **radixSortedData = extract_radix_sorted_data(dataSetProperties.dataSetFileContents, dataSetProperties.entryCount, dataSetProperties.fieldCount, dataSetProperties.delimiter);
	
	
	DataSetAnalysis configuredDataSet;
	configuredDataSet.dataSetProperties = &dataSetProperties;
	configuredDataSet.radixSortedData = radixSortedData;
	configuredDataSet.plottableDataDirectory = extract_directory_properties(preprocessedDataDirectory);
	
	return configuredDataSet;
}




char **blindly_extract_data_set(const char* dataSetFilePathName, int lineCount)
{
	/*-----------   Capture File Contents in an Array of Strings   -----------*/
	char **fileContentsVerbatim = read_file_contents(dataSetFilePathName, lineCount);
	print_string_array(fileContentsVerbatim, lineCount, "Blind extraction of File Contents    -->    Formatting");
	return fileContentsVerbatim;
}




/**
 * extract_and_format_data_set
 *
 * Purpose:
 *   This function serves as a comprehensive step to refine raw file contents
 *   for compatibility within a larger CSV data analysis framework. It ensures
 *   each row of the dataset is free of problematic characters and correctly
 *   formatted for subsequent operations such as field parsing, data type
 *   standardization, and eventual plotting or analysis.
 *
 * Summary:
 *   1) Accepts raw file contents, total number of lines, field count, and the delimiter.
 *   2) Removes or transforms disruptive characters and irregular spacing.
 *   3) Standardizes each line to align every data field with expected formats
 *      (e.g., numeric fields are numeric, missing fields are marked).
 *   4) Returns a revised array of strings, preserving the original header and
 *      ensuring uniform data entries below it.
 */
char **extract_and_format_data_set(char **fileContents, int lineCount, int fieldCount, const char *delimiter)
{
	/*-----------   Begin Preprocessing File Contents to Standardize the Format and Achieve/Maintain Compatibility of the Contents   -----------*/
	char **formattedFileContents = fileContents; /// Prepare a reference to the same buffer for eventual formatted output
	
	
	// Extract and prune the header line (fileContents[1] is the actual header in current data)
	char *headerLine = prune_and_trim_problematic_characters_from_string(fileContents[1], delimiter, fieldCount);
	for(int i = 1; i < lineCount; i++) // Loop through each line to eliminate problematic characters, ensuring uniformity for all rows
	{
		formattedFileContents[i] = prune_and_trim_problematic_characters_from_string(fileContents[i], delimiter, fieldCount);
	}
	
	
	
	/// Individually refine each line by imposing stricter formatting rules and reconciling data fields
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
	//print_file_contents(formattedFileContents, lineCount);
	return formattedFileContents; // Return the fully standardized array of strings (header + formatted data)
}









/**
 * process_data_set_for_analysis
 *
 *   "process_data_set_for_analysis" is the main entry point for reading,
 *   standardizing, and storing CSV data within this analysis framework.
 *   It orchestrates the steps to parse raw file contents, identify the file
 *   delimiter, format the data into a consistent structure, and finally save
 *   it into organized directories for subsequent processing.
 *
 * Purpose & Rationale:
 *   • Automatically read the csv/text file and count its lines.
 *   • Parse the data, detect and apply the correct delimiter, and ensure that
 *     each field is standardized.
 *   • Produce a refined data set with consistent headers, plottability analysis,
 *     and missing-data tracking before storing it in directories.
 *   • Provide an easy-to-reference object (DataSetAnalysis) that captures all
 *     relevant properties of the data set, along with the location where
 *     processed results are written for further analysis.
 *
 * Means of Execution at a High Level:
 *   1. Read the file contents into memory and determine the file’s delimiter.
 *   2. Verify consistency between the recorded number of lines and actual
 *      content lines.
 *   3. Identify how many fields the data set contains and then extract &
 *      format the data to align with the expected format.
 *   4. Invoke specialized analysis (analyze_data_set_properties) to classify
 *      fields and detect numeric vs. textual data.
 *   5. Print and log the extracted data for debugging or user feedback.
 *   6. Write the organized output to files, categorizing data fields and
 *      storing them in a new directory path.
 *   7. Configure the final DataSetAnalysis structure that clients of this
 *      function can further utilize for more advanced analyses (plotting,
 *      statistics, etc.).
 */
DataSetAnalysis process_data_set_for_analysis(const char* dataSetFilePathName)
{
	/*-----------   Capture File Contents in an Array of Strings   -----------*/
	int lineCount = count_file_lines(dataSetFilePathName, MAX_NUM_FILE_LINES);
	char **fileContents = parse_file_contents(dataSetFilePathName, lineCount);
	const char *delimiter = identify_delimiter(fileContents, lineCount);
	if(count_array_strings(fileContents) != lineCount)
	{
		perror("\n\nError: count_array_strings != count_file_lines in 'process_data_set_for_analysis'.\n");
		//return dataSetFilePathName;
	}
	int fieldCount = count_data_fields(fileContents[0]); // Count the fields present in the header line (index 0) to assess the data's structure
	
	
	
	
	/*-----------   Extract and Format Data   -----------*/
	/// Extract and format all lines for consistent representation
	// 'extract_and_format_data_set' cleans, trims, and standardizes each row.
	/// This ensures each data row aligns with the expected CSV structure.
	char **formattedContents = extract_and_format_data_set(fileContents, lineCount, fieldCount, delimiter);

	
	
	/// Perform further examination of the data set, capturing missing-values, plottable fields, etc.
	// 'analyze_data_set_properties' returns a struct with many details about the data set (including field types).
	DataSetProperties dataSetProperties = analyze_data_set_properties(dataSetFilePathName);
	print_string_array(formattedContents, lineCount, "Formatted Extraction of File Contents    -->    Preprocessing");
	
	
	
	
	
	/*-----------   Process and Store Data   -----------*/
	/// Write the standardized contents out to new files for each field, then build a DataSetAnalysis object
	// 'write_data_set' writes out formatted data. The path to these new data files is stored in 'preprocessedDataDirectory'.
	const char *preprocessedDataDirectory = write_data_set(formattedContents, dataSetFilePathName, lineCount, fieldCount, delimiter);
	DataSetAnalysis dataAnalysis = configure_data_set_analysis(dataSetProperties, preprocessedDataDirectory); // sets up 'DataSetAnalysis' by linking the refined DataSetProperties with new directory paths, etc.
	
	
	
	
	
	/// Logging/feedback that describes the successful completion of data preprocessing
	print_string_array(dataAnalysis.dataSetProperties->fieldNameTypePairs, fieldCount, "Pairs of field names and their corresponding types in 'process_data_set_for_analysis'");
	printf("\n\n\n\n\n\n\n\n***********************************************************************");
	printf("\n Preliminary extraction and preprocessing of data set complete.");
	printf("\n		Summary: Preprocessed the data set to ensure operational compatibility throughout the program. \nExtracting and formatting the contents of the data set achieved by having: ");
	printf("\n		1. Extracted the data set contents into an array of strings.");
	printf("\n		2. Pruned and trimmed the data set contents to remove problematic characters.");
	printf("\n		3. Organized entries of the file corresponding to the expected line order of header and data entries.");
	printf("\n		4. Parsed the file contents to independently capture all entries for each of the data set's fields.");
	printf("\n		5. Categorized the data parameters into plottable and unplottable values.");
	printf("\n		6. Created directories for plottable fields of data entries, named & located with respect to the original file.");
	printf("\n		7. Data conditionally written into files for each field, named & placed in accordance with the aforementioned directory and nomenclature.");
	printf("\n\n***********************************************************************");
	return dataAnalysis; 	// Return the DataSetAnalysis instance for further usage by other parts of the program
}








/**
 * perform_full_analysis_and_modeling
 *
 * This function orchestrates the comprehensive analysis and modeling of a preprocessed data set.
 * It creates a dedicated results directory, identifies numeric data files, performs statistical evaluation,
 * writes out results summarizing the data’s distribution and statistical properties, and finally generates
 * MATLAB scripts for further visualization. The steps within this function ensure a streamlined approach
 * to examining the numeric fields in the provided directory while building a robust analytical framework.
 *
 * PURPOSE & RATIONALE:
 * 1. Automate comprehensive analysis steps for all plottable data fields.
 * 2. Generate per-field statistics, histograms, and other metrics in an
 *    easily accessible format (.txt files).
 * 3. Provide a convenient entry point for higher-level modeling tasks by
 *    producing MATLAB scripts for both individual fields and an
 *    all-encompassing script for further data exploration.
 *
 * EXECUTION STEPS (HIGH-LEVEL):
 * 1. Establish a dedicated directory for storing analysis results.
 * 2. Gather a list of .txt data files (each representing a single plottable field).
 * 3. For each data file, read in the numeric values, compute statistics (mean,
 *    standard deviation, skewness, etc.), and record to results files.
 * 4. Create MATLAB scripts for comprehensive data visualization and exploration.
 * 5. Return the path to the newly created analysis-results directory.
 */
const char *perform_full_analysis_and_modeling(const char *preprocessedDataDirectoryPath)
{
	/// Create new results directory for storing analysis outcomes.
	char *resultsDirectory = combine_strings(preprocessedDataDirectoryPath, "_Full_Analysis_Results");
	create_directory(resultsDirectory, ""); // create_directory ensures directory creation
	

	
	/// Acquire the full list of plottable/potential data files from the preprocessed directory.
	int fileCount = count_files_in_directory(preprocessedDataDirectoryPath);
	char **filePaths = get_file_pathnames_in_directory(preprocessedDataDirectoryPath);
	char **fieldNames = (char**)malloc(fileCount * sizeof(char*)); 	// Initialize an array to store the extracted field names.
	
	
	/// Loop over each file to verify file extension, parse contents, then perform statistical computations.
	for (int i = 0; i < fileCount; i++)
	{
		const char *fieldFilePath = filePaths[i];
		
		// Check extension to ensure it's a .txt data field file
		const char *ext = identify_file_extension(fieldFilePath);
		if (strcmp(ext, ".txt") != 0) continue;
		
		int lineCount = count_file_lines(fieldFilePath, MAX_NUM_FILE_LINES);
		//if (lineCount < 2) continue;  // Expect at least one header line + data lines
		
		char **contents = read_file_contents(fieldFilePath, lineCount);
		if (!contents) continue;
		
		
		
		// The first line is the field name
		char *fieldName = strdup(contents[0]);
		
		fieldNames[i] = (char*)malloc(strlen(fieldName) * sizeof(char));
		fieldNames[i] = duplicate_string(fieldName);
		
		
		// Remaining lines are numeric data
		double *data = allocate_memory_double_ptr(lineCount - 1);
		for (int j = 1; j < lineCount; j++)
		{
			// Shift indexing so that 'data' array starts at index 0.
			// Depending on the content at the final line, handle indexing carefully.
			if(j == (lineCount-2))
			{
				data[j] = atof(contents[j]);
			}
			else
			{
				data[j-1] = atof(contents[j]);
			}
			
		}
		
		//print_array(lineCount, data, fieldName);
		
		
		int n = lineCount - 1;
		deallocate_memory_char_ptr_ptr(contents, lineCount);
		//if (n < 2)
		//{
		//	free(data);
		//	free(fieldName);
		//	continue;
		//}
		
		
		/// Run statistical analysis on the extracted numeric data
		StatisticalReport results = analyze_numeric_data(data, n, resultsDirectory, fieldName);
	
		
		// Write the computed statistics to an output file in 'resultsDirectory'
		char *analysisFilePath = combine_strings(resultsDirectory, "/");
		analysisFilePath = combine_strings(analysisFilePath, fieldName);
		analysisFilePath = combine_strings(analysisFilePath, "_full_analysis.txt");
		
		FILE *analysisFile = fopen(analysisFilePath, "w+");
		if (analysisFile)
		{
			fprintf(analysisFile, "Field: %s\n", fieldName);
			fprintf(analysisFile, "Count: %d\n", n);
			fprintf(analysisFile, "Mean: %.17g\n", results.mean);
			fprintf(analysisFile, "Std Dev: %.17g\n", results.std_dev);
			fprintf(analysisFile, "Skewness: %.17g\n", results.skewness);
			fprintf(analysisFile, "Anderson-Darling Statistic (AD): %.17g\n", results.ad_stat);
			
			// Histogram data
			fprintf(analysisFile, "\nHistogram:\n");
			fprintf(analysisFile, "NumBins: %d\n", results.histogram.num_bins);
			fprintf(analysisFile, "BinWidth: %.17g\n", results.histogram.bin_width);
			fprintf(analysisFile, "Min: %.17g\nMax: %.17g\n", results.histogram.min_value, results.histogram.max_value);
			for (int b = 0; b < results.histogram.num_bins; b++)
			{
				double bin_start = results.histogram.min_value + b * results.histogram.bin_width;
				double bin_end = bin_start + results.histogram.bin_width;
				fprintf(analysisFile, "Bin %d: [%.17g, %.17g): %d\n", b, bin_start, bin_end, results.histogram.bins[b]);
			}
			
			fclose(analysisFile);
		}
		
	}
	/// Generate supplemental MATLAB scripts to visualize the data analysis.
	//generate_matlab_model_scripts_unified(resultsDirectory, fieldNames, fileCount, 1);  // Put per-field scripts under analysisDir/MATLAB_Scripts
	generate_matlab_scripts_unified(resultsDirectory, fieldNames, fileCount, MATLAB_SCRIPTS_BOTH, "MATLAB_Scripts");  // Put both per-field scripts under analysisDir/MATLAB_Scripts AND one big script at analysisDir:

	
	return resultsDirectory; // Return the path where the final results (analysis reports and scripts) are stored.
}




























