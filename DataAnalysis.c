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



DataSetProperties analyze_data_set_properties(const char *filePathName)
{
	int lineCount = count_file_lines(filePathName, MAX_NUM_FILE_LINES);
	char **fileContents = parse_file_contents(filePathName, lineCount);
	const char *delimiter = identify_delimiter(fileContents, lineCount);
	int fieldCount = count_data_fields(fileContents[0]);
	char **formattedContents = extract_and_format_data_set(fileContents, lineCount, fieldCount, delimiter);
	char **nameTypePairs = capture_data_set_header_for_plotting(formattedContents[0], formattedContents, delimiter);
	
	
	int *missingDataCount = count_missing_values(formattedContents, lineCount, fieldCount, delimiter, formattedContents[1]);;
	int *plottabilityStatus = identify_plottable_fields(nameTypePairs, fieldCount, ";");
	//char** commonDataTypes = determine_common_data_entry_types(formattedContents, lineCount, fieldCount, delimiter);
	char **fieldNamesTokenized = split_tokenized_string(formattedContents[0], delimiter, fieldCount);
	char *plottableFieldNames;
	for(int i = 0; i < fieldCount; i++)
	{
		char *fieldNames;
		if(plottabilityStatus[i] == 0)
		{
			fieldCount--;
			
		}
		else if(plottabilityStatus[i] == 1)
		{
			fieldNames = combine_strings(fieldNames, fieldNamesTokenized[i]);
		}
		else
		{
			perror("\n\nError: Unrecognized plottability status in 'analyze_data_set_properties'.\n");
			exit(1);
		}
		
		if(i == fieldCount)
		{
			plottableFieldNames = fieldNames;
		}
		
	}
	//printf("\n\nPlottable Field Names: %s", plottableFieldNames);
	
	
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
	
	
	
	//free(delimiter);
	//free(missingDataCount);
	//free(plottabilityStatus);
	//deallocate_memory_char_ptr_ptr(nameTypePairs, fieldCount);
	//deallocate_memory_char_ptr_ptr(commonDataTypes, fieldCount);
	//deallocate_memory_char_ptr_ptr(formattedContents, lineCount);
	//deallocate_memory_char_ptr_ptr(fileContents, lineCount);
	return dataSetProperties;
}


DataSetProperties capture_data_set_configurations(const char *filePathName, char **formattedContents, int lineCount, int fieldCount, const char *delimiter)
{
	char **nameTypePairs = capture_data_set_header_for_plotting(formattedContents[0], formattedContents, delimiter);
	int *missingDataCount = count_missing_values(formattedContents, lineCount, fieldCount, delimiter, formattedContents[1]);;
	int *plottabilityStatus = identify_plottable_fields(nameTypePairs, fieldCount, ";");
	
	
	char **commonDataTypes = determine_common_data_entry_types(formattedContents, lineCount, fieldCount, delimiter);
	
	
	char **fieldNamesTokenized = split_tokenized_string(formattedContents[0], delimiter, fieldCount);
	//printf("\n\nformattedContents fieldCount %d: %s", fieldCount, formattedContents[0]);
	//print_string_array(fieldNamesTokenized, fieldCount, "fieldNamesTokenized");
	
	char *plottableFieldNames;
	int plottableFieldCount = fieldCount;
	char *fieldNames;
	for(int i = 0; i < fieldCount; i++)
	{
		if(plottabilityStatus[i] == 0)
		{
			plottableFieldCount--;
		}
		else if(plottabilityStatus[i] == 1)
		{
			fieldNames = combine_strings(fieldNames,fieldNamesTokenized[i]);
			if(i < (fieldCount-1))
			{
				fieldNames = combine_strings(fieldNames, ",");
			}
			
			//printf("\n\nField Name%d: %s", i, fieldNames);
		}
		else
		{
			perror("\n\nError: Unrecognized plottability status in 'analyze_data_set_properties'.\n");
			exit(1);
		}
		
		
		
		if(i == (fieldCount-1))
		{
			plottableFieldNames = fieldNames;
		}
	}
	//printf("\n\nPlottable Field Names %d: %s",plottableFieldCount, plottableFieldNames);
	
	
	
	DataSetProperties dataSetProperties;
	dataSetProperties.dataSetFilePathName = filePathName;
	dataSetProperties.entryCount = lineCount;
	dataSetProperties.fieldCount = fieldCount;
	dataSetProperties.delimiter = delimiter;
	dataSetProperties.dataSetFileContents = formattedContents;
	
	
	dataSetProperties.dataSetFieldNames = formattedContents[0];
	dataSetProperties.fieldNameTypePairs =  nameTypePairs;
	dataSetProperties.missingDataCount = missingDataCount;
	dataSetProperties.plottabilityStatus = plottabilityStatus;
	dataSetProperties.commonDataTypes = commonDataTypes;
	
	
	
	
	//free(delimiter);
	//free(missingDataCount);
	//free(plottabilityStatus);
	//deallocate_memory_char_ptr_ptr(nameTypePairs, fieldCount);
	//deallocate_memory_char_ptr_ptr(commonDataTypes, fieldCount);
	//deallocate_memory_char_ptr_ptr(formattedContents, lineCount);
	//deallocate_memory_char_ptr_ptr(fileContents, lineCount);
	return dataSetProperties;
	
	
}

DataSetAnalysis configure_data_set_analysis(DataSetProperties dataSetProperties, const char *preprocessedDataDirectory)
{
	double **radixSortedData = extract_radix_sorted_data(dataSetProperties.dataSetFileContents, dataSetProperties.entryCount, dataSetProperties.fieldCount, dataSetProperties.delimiter);
	
	
	
	DataSetAnalysis configuredDataSet;
	configuredDataSet.dataSetProperties = &dataSetProperties;
	
	
	configuredDataSet.radixSortedData = radixSortedData;
	configuredDataSet.plottableDataDirectory = extract_directory_properties(preprocessedDataDirectory);
	
	
	
	//deallocate_memory_double_ptr_ptr(radixSortedData, dataSetProperties.entryCount);
	return configuredDataSet;
}







char **extract_and_format_data_set(char **fileContents, int lineCount, int fieldCount, const char *delimiter)
{
	/*-----------   Begin Preprocessing File Contents to Standardize the Format and Achieve/Maintain Compatibility of the Contents   -----------*/
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
	//print_file_contents(formattedFileContents, lineCount);
	return formattedFileContents;
}






DataSetAnalysis process_data_set_for_analysis(const char* dataSetFilePathName)
{
	/*-----------   Capture File Contents in an Array of Strings   -----------*/
	int lineCount = count_file_lines(dataSetFilePathName, MAX_NUM_FILE_LINES);
	char **fileContents = parse_file_contents(dataSetFilePathName, lineCount);
	const char *delimiter = identify_delimiter(fileContents, lineCount);
	
	
	if(count_array_strings(fileContents) != lineCount)
	{
		perror("\n\nError: count_array_strings != count_file_lines in 'preprocess_data_set_for_analysis'.\n");
		//return dataSetFilePathName;
	}
	int fieldCount = count_data_fields(fileContents[0]);
	
	
	/*-----------   Extract and Format Data   -----------*/ ///*footnote 1.1 in 'DataAnalysis.h'*
	char **formattedContents = extract_and_format_data_set(fileContents, lineCount, fieldCount, delimiter);
	//DataSetProperties dataSetProperties = capture_data_set_configurations(dataSetFilePathName, formattedContents, lineCount, fieldCount, delimiter);
	DataSetProperties dataSetProperties = analyze_data_set_properties(dataSetFilePathName);
	print_string_array(formattedContents, lineCount, "Formatted Extraction of File Contents    -->    Preprocessing");
	
	
	/*-----------   Process and Store Data   -----------*/ ///*footnote 1.2 in 'DataAnalysis.h'*
	const char *preprocessedDataDirectory = write_data_set(formattedContents, dataSetFilePathName, lineCount, fieldCount, delimiter);
	DataSetAnalysis dataAnalysis = configure_data_set_analysis(dataSetProperties, preprocessedDataDirectory);
	
	
	
	print_string_array(dataAnalysis.dataSetProperties->fieldNameTypePairs, fieldCount, "Pairs of field names and their corresponding types in 'process_data_set_for_analysis'");
	printf("\n\n\n\n\n\n\n\n***********************************************************************");
	printf("\n Preliminary extraction and preprocessing of data set complete.");
	// list the primary markers achieved in the preliminary phase of the data set analysis
	printf("\n		Summary: Preprocessed the data set to ensure operational compatibility throughout the program. \nExtracting and formatting the contents of the data set achieved by having: ");
	printf("\n		1. Extracted the data set contents into an array of strings.");
	printf("\n		2. Pruned and trimmed the data set contents to remove problematic characters.");
	printf("\n		3. Organized entries of the file corresponding to the expected line order of header and data entries.");
	printf("\n		4. Parsed the file contents to independently capture all entries for each of the data set's fields.");
	printf("\n		5. Categorized the data parameters into plottable and unplottable values.");
	printf("\n		6. Created directories for plottable fields of data entries, named & located with respect to the original file.");
	printf("\n		7. Data conditionally written into files for each field, named & placed in accordance with the aforementioned directory and nomenclature.");
	printf("\n\n***********************************************************************");
	return dataAnalysis;
}











const char *preprocess_data_set_for_analysis(const char* dataSetFilePathName)
{
	/*-----------   Capture File Contents in an Array of Strings   -----------*/
	int lineCount = count_file_lines(dataSetFilePathName, MAX_NUM_FILE_LINES);
	char **fileContents = parse_file_contents(dataSetFilePathName, lineCount);
	const char *delimiter = identify_delimiter(fileContents, lineCount);
	
	
	if(count_array_strings(fileContents) != lineCount)
	{
		perror("\n\nError: count_array_strings != count_file_lines in 'preprocess_data_set_for_analysis'.\n");
		return dataSetFilePathName;
	}
	int fieldCount = count_data_fields(fileContents[0]);
	
	
	/*-----------   Extract and Format Data   -----------*/ ///*footnote 1.1 in 'DataAnalysis.h'*
	char **formattedContents = extract_and_format_data_set(fileContents, lineCount, fieldCount, delimiter);
	
	
	print_file_contents(formattedContents, lineCount);
	
	/*-----------   Process and Store Data   -----------*/ ///*footnote 1.2 in 'DataAnalysis.h'*
	const char *preprocessedDataDirectory = write_data_set(formattedContents, dataSetFilePathName, lineCount, fieldCount, delimiter);
	
	
	
	
	
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
	return preprocessedDataDirectory;
}




char **blindly_extract_data_set(const char* dataSetFilePathName, int lineCount)
{
	/*-----------   Capture File Contents in an Array of Strings   -----------*/
	char **fileContentsVerbatim = read_file_contents(dataSetFilePathName, lineCount);
	print_string_array(fileContentsVerbatim, lineCount, "Blind extraction of File Contents    -->    Formatting");
	
	
	return fileContentsVerbatim;
}




const char *perform_full_analysis_and_modeling(const char *preprocessedDataDirectoryPath)
{
	// Create a results directory for analysis
	char *resultsDirectory = combine_strings(preprocessedDataDirectoryPath, "_Full_Analysis_Results");
	create_directory(resultsDirectory, ""); // create_directory ensures directory creation
	
	
	//
	/*
	 // Create a results directory for analysis
	 char *resultsDirectoryForAppend = (char*)malloc(string_length("_Full_Analysis_Results") + 1 * sizeof(char));
	 resultsDirectoryForAppend = "_Full_Analysis_Results";
	 
	 
	 char *resultsDirectory = combine_strings(preprocessedDataDirectoryPath, resultsDirectoryForAppend);
	 char *resultsDirectoryName = find_name_from_path(resultsDirectory);
	 create_directory(resultsDirectory, resultsDirectoryName); // create_directory ensures directory creation
	 
	 
	 */
	
	
	// Get list of plottable data files in preprocessed directory
	int fileCount = count_files_in_directory(preprocessedDataDirectoryPath);
	char **filePaths = get_file_pathnames_in_directory(preprocessedDataDirectoryPath);
	char **fieldNames = (char**)malloc(fileCount * sizeof(char*));
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
		
		// Analyze the numeric data
		StatisticalReport results = analyze_numeric_data(data, n, resultsDirectory, fieldName);
		
		// Write analysis results to a file
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
	//generational_individuality_matlabulous_scriptorion
	generate_individual_matlab_scripts(resultsDirectory, fieldNames, count_files_in_directory(resultsDirectory));
	generate_comprehensive_matlab_script(resultsDirectory, fieldNames, count_files_in_directory(resultsDirectory));
	
	
	
	return resultsDirectory;
}

































