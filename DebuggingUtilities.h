//  DebuggingUtilities.h
//  Automated_CSV_Data_Analysis
//  DavidRichardson02
/**
 * DebuggingUtilities code: Provides a number of utility functions for debugging
 * and testing purposes. These functions facilitate the printing of various data
 * types, such as strings, arrays, and files, to the console.
 */

#ifndef DebuggingUtilities_h
#define DebuggingUtilities_h

#include <stdio.h>

// ------------- Helper Functions for Printing Various Representations of String
// and Arrays for Debugging -------------
/// \{
void print_string(char *string, char *label); // Prints a single string.
void print_string_array(char **stringArray, int stringCount,
						char *label); // Prints an array of strings.
void print_string_array_array(
							  char ***stringArrayArray, int stringArraysCount, int stringSubArraysCount,
							  char *label); // NOTE: assumes the dimension of all sub-arrays is the same
void print_array(int n, double *data,
				 char *label); // Prints an array of doubles.
void print_array_array(double **data, int rows, int columns,
					   char *label); // Prints a 2D array of doubles.
void print_char_ptr_array(char *charPtrArr[], int stringCount,
						  char *label); // Prints an array of char pointers.
										/// \}

// ------------- Helper Functions for File Display -------------
/// \{
void print_data_file(
					 const char *filePathName); // Prints key information about a file as well as
												// its contents.
void print_file_contents(
						 char **fileContents,
						 int countFileLines); // Prints file contents from an array of strings.

void print_directory_data_files(
								const char *directoryPathName); // Prints the files in a directory.

// IMPORTANT: This function is only intended to be used with the returned
// pathname from: 'write_data_set', this is for two reasons
/// 1. The generic 'find_file_directory_path' includes the last '/' in the path,
/// so   =>  "/home/user/data_set_directory/"  <=  is returned, when this
/// expects the string to end in the directory name, not the directory name and
/// a trailing '/'
///
/// 2. This function relies heavily on the 'print_directory_data_files' which
/// also does not expect a trailing '/'
void print_written_data_set(
							const char *writtenDataSetDirectoryPathName); // Prints the directory where
																		  // 'write_data_set' wrote the
																		  // data set.
																		  /// \}

// print_data_set_properties(DataSetProperties dataSetProperties); // Prints the
// properties of a data set.
/**
 
 printf("\n\n\n\ndataSetProperties: ");
 printf("\n\ndataSetHeader: %s", dataSetProperties.dataSetFieldNames);
 printf("\ndelimiter: %s", dataSetProperties.delimiter);
 printf("\nentryCount: %d        fieldCount: %d      ",
 dataSetProperties.entryCount, dataSetProperties.fieldCount);
 
 print_string_array(dataSetProperties.fieldNameTypePairs,
 dataSetProperties.fieldCount, "fieldNameTypePairs");
 
 for (int i = 0; i < dataSetProperties.entryCount; i++)
 {
 //printf("\n%s", dataSetProperties.fieldNameTypePairs[i]);
 //printf("\n%s", dataSetProperties.fieldNameTypePairs[i]);
 
 }
 */

#endif /* DebuggingUtilities_h */
