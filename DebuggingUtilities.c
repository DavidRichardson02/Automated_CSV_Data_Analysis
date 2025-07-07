//  DebuggingUtilities.c
//  Automated_CSV_Data_Analysis
//  DavidRichardson02


#include "DebuggingUtilities.h"
#include "CommonDefinitions.h"
#include "GeneralUtilities.h"
#include "StringUtilities.h"
#include "FileUtilities.h"






/**
 * print_string
 * Prints a string with a label.
 *
 * @param string string to be printed.
 */
void print_string(char* string, char*label)
{
	// Check for NULL input and handle error.
	if (string == NULL){ perror("\n\nError: string was NULL in 'print_string'.\n");      exit(1); }
	printf("\n\n\n\n");
	printf("\n%s: ", label);
	//printf("\n\n-----------------------------------------------------------------------------------------\n");
	printf("'%s'", string);
	//printf("\n-----------------------------------------------------------------------------------------\n\n");
	printf("\n\n");
}




/**
 * print_string_array
 * Prints the elements of a string array with a label.
 *
 * @param stringArray Array of strings to be printed.
 * @param stringCount Number of strings in the array.
 * @param label Label to be printed before printing the array.
 */
void print_string_array(char** stringArray, int stringCount, char*label)
{
	// Check for NULL input and handle error.
	if(stringArray == NULL){ perror("\n\nError printing array of strings in 'print_string_array'.");      exit(1); }
	
	printf("\n\n\n\n\n\n----------------------------------------------------------------\n");
	printf(" %s:", label);
	printf("\n----------------------------------------------------------------");
	for(int i = 0; i < stringCount; i++)
	{
		printf("\n	%s", stringArray[i]);
	}
}




/**
 * print_string_array_array
 *
 * This function prints the contents of an array of arrays of strings. It is assumed that all sub-arrays have the same dimension.
 * The output is formatted with a label and each string is printed on a new line.
 * The function first checks for NULL input to prevent errors. It then prints a label for the array.
 * It then iterates over the main array and for each sub-array, it calculates the length of the label, allocates memory for it,
 * and then prints the sub-array using the 'print_string_array' function.
 *
 * @param stringArrayArray The 2D array of strings to be printed.
 * @param stringArraysCount The number of sub-arrays in the main array.
 * @param stringSubArraysCount The number of strings in each sub-array.
 * @param label The label to be printed before printing the array.
 */
void print_string_array_array(char*** stringArrayArray, int stringArraysCount, int stringSubArraysCount, char* label)
{
	// Check for NULL input and handle error.
	if (stringArrayArray == NULL){ perror("\n\nError: stringArrayArray was NULL in 'print_string_array_array'.\n");      exit(1); }
	
	
	printf("\n\n\n\n\n\n\n\nprint_string_array_array %s =========================================================================================", label);
	printf("\n\n\n\n%s: \n", label);
	// Iterate over the main array.
	for (int i = 0; i < stringArraysCount; i++)
	{
		// For each sub-array, calculate the length of the label and allocate memory for it.
		const char* firstElement = ARRAY_ELEMENT_TO_STRING(stringArrayArray[i], 0);
		if (firstElement == NULL)
		{
			perror("\n\nError: first element of a sub-array was NULL in 'print_string_array_array'.\n");
			continue;
		}
		
		size_t labelLength = string_length(firstElement);
		char* subArrayLabel = allocate_memory_char_ptr(string_length(label) + labelLength + 1);
		if (!subArrayLabel)
		{
			perror("\n\nError: Memory allocation failed for subArrayLabel in 'print_string_array_array'.\n");
			return;
		}
		
		// Initialize subArrayLabel and copy the label.
		subArrayLabel[0] = '\0';
		sprintf(subArrayLabel, "%s", firstElement);
		
		// Print the sub-array.
		print_string_array(stringArrayArray[i], stringSubArraysCount, subArrayLabel);
	}
	printf("\n\n\n=========================================================================================\n\n");
}




/**
 * print_array
 * Prints the elements of a double array with a label.
 *
 * @param n Number of elements in the array.
 * @param data Array of doubles to be printed.
 * @param label Label to be printed before printing the array.
 */
void print_array(int n, double *data, char*label)
{
	// Check for NULL input and handle error.
	if (data == NULL)
	{
		perror("\n\nError: data was NULL in 'print_array'.\n");
		//exit(1);
	}
	
	
	printf("\n\n%s: ", label);
	// Loop through the characters print
	for(int i = 0; i < n; i++)
	{
		printf("\n %-2d 	%.17g ", i, data[i]);
	}
	printf("\n");
}




/**
 * print_array_array
 *
 * This function prints the contents of a 2D array with each element
 * displayed to 17 decimal places. The output is formatted with a label and
 * surrounded by a visual border for clarity.
 *
 * @param data A pointer to a pointer of doubles representing a 2D array.
 * @param rows The number of rows in the 2D array.
 * @param columns The number of columns in each row.
 * @param label A string label for the array.
 */
void print_array_array(double **data, int rows, int columns, char*label)
{
	// Check for NULL input and handle error.
	if (data == NULL){ perror("\n\nError: data was NULL in 'print_array_array'.\n");      exit(1); }
	
	printf("\n\n\n\n\n\n\nprint_array_array %s =========================================================================================", label);
	printf("\n\n%s: \n", label);
	
	for (int i = 0; i < rows; i++)
	{
		if(data[i] == NULL)
		{
			perror("\n\nError: data[i] was NULL in 'print_array_array'.\n");
			continue;
		}
		
		// Make row label "Row i" where i is the row number.
		char *rowLabel = allocate_memory_char_ptr(strlen("Row ") + 1);
		snprintf(rowLabel, sizeof(rowLabel), "Row %d", i);
		print_array(columns, data[i], rowLabel);
		
		
		//for (int j = 0; j < columns; j++)
		//{
		//	printf("%.17g ", data[i][j]);
		//}
		
		
		printf("\n"); // Newline after each row
	}
	printf("\n\n=========================================================================================\n\n");
}




/**
 * print_char_ptr_array
 * Prints the elements of a string array.
 *
 * @param charPtrArr Array of strings to be printed.
 * @param stringCount Number of strings in the array.
 */
void print_char_ptr_array(char *charPtrArr[],  int stringCount, char* label)
{
	// Check for NULL input and handle error.
	if (charPtrArr == NULL){ perror("\n\nError: charPtrArr was NULL in 'print_char_ptr_array'.\n");      exit(1); }
	
	
	printf("\n\n\n\n\n\n\n\n%s: \n", label);
	// Loop through the elements array and print
	for (int i = 0; i < stringCount; i++)
	{
		printf("%s", charPtrArr[i]);
	}
}









/**
 * print_data_file
 * Prints the contents of a file to the console.
 * This function reads the contents of a file, including the file header, extension, directory path, name, and line count
 * and prints them to the console. It then prints the contents of the file line by line.
 *
 * @param filePathName The path to the file to be printed.
 */
void print_data_file(const char* filePathName)
{
	char *fileExtension = identify_file_extension(filePathName);
	char *fileName = find_name_from_path(filePathName);
	int countFileLines = count_file_lines(filePathName, MAX_NUM_FILE_LINES);
	char **fileContents = read_file_contents(filePathName, countFileLines);
	
	
	printf("\n\n\n\n\n\n--------------------------------------------------\n");
	printf(" %s", fileName);
	printf("\n--------------------------------------------------");
	printf("\n	Pathname: %s", filePathName);
	printf("\n	Extension: %s", fileExtension);
	printf("\n	Count of File Lines: %d", countFileLines);
	
	
	
	printf("\n\nFile Contents:");
	for(int i = 0; i < countFileLines; i++)
	{
		printf("\n    Line %-2d:       %-5s", i, fileContents[i]);
	}
	
	
	printf("\n***********************************************************************");
	
	//free(fileExtension);
	//free(fileName);
	//deallocate_memory_char_ptr_ptr(fileContents, countFileLines);
}




/**
 * print_file_contents
 * Prints the contents of a file from an array of strings to the console,
 * where each string represents a line in the file.
 *
 * @param fileContents The contents of the file to be printed.
 * @param countFileLines The number of lines in the file.
 */
void print_file_contents(char** fileContents, int countFileLines)
{
	if(fileContents == NULL)
	{
		perror("\n\nError printing array of strings in 'print_file_contents'.");
		exit(1);
	}
	
	printf("\n\n\n\n\n\nFile Contents: \n");
	for(int i = 0; i < countFileLines; i++)
	{
		printf("\n%-5s", fileContents[i]);
	}
}







void print_written_data_set(const char* writtenDataSetDirectoryPathName) // Prints the directory where 'write_data_set' wrote the data set.
{
	char *writtenDataSetDirectory = find_directory_name_from_path(writtenDataSetDirectoryPathName);
	char **writtenDataSetFilePaths = get_file_pathnames_in_directory(writtenDataSetDirectoryPathName);
	int countWrittenDataSetFiles = count_files_in_directory(writtenDataSetDirectoryPathName);
	
	
	printf("\n\n\n\n\n\n\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
	printf("\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	printf("\n Processed Dataset Directory %s\n", writtenDataSetDirectory);
	printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
	printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
	printf("\n|\n| Directory Pathname: %-5s", writtenDataSetDirectoryPathName);
	printf("\n| Count of Files: %-5d", countWrittenDataSetFiles);
	printf("\n|\n| Files in Directory: \n|");
	for(int i = 0; i < countWrittenDataSetFiles; i++)
	{
		//printf("\n| 	%d:       %-5s", i, writtenDataSetFilePaths[i]);
		print_data_file(writtenDataSetFilePaths[i]);
	}
	free(writtenDataSetDirectory);
	deallocate_memory_char_ptr_ptr(writtenDataSetFilePaths, countWrittenDataSetFiles);
}











// Remove later
void print_directory_data_files(const char* directoryPathName)
{
	//char *directoryPathName = find_file_directory_path(filePathName); // captures directory path with last '/' included?
	/** /// Trim the trailing '/' from the directory name
	 size_t directoryLength = string_length(directoryPathName);
	 if(directoryPathName[directoryLength - 1] == '/')
	 {
	 printf("\nNeed to trim the trailing '/' from the directory name.");
	 char *trimmedDirectoryPathName = allocate_memory_char_ptr(directoryLength-1); // Allocate memory for the trimmed directory path name
	 strncpy(trimmedDirectoryPathName, directoryPathName, directoryLength-1); // Copy the trimmed directory path name
	 printf("\n\nTrimmed Directory Pathname: %s", trimmedDirectoryPathName);
	 free(directoryPathName); // Free the original directory path name
	 directoryPathName = allocate_memory_char_ptr(directoryLength-1); // Allocate memory for the directory path name
	 strncpy(directoryPathName, directoryPathName, directoryLength-1); // Copy the trimmed directory path name
	 free(trimmedDirectoryPathName); // Free the trimmed directory path name
	 printf("\n\nDirectory Pathname: %s", directoryPathName);
	 }**/
	//char *directoryName = find_directory_name_from_path(directoryPathName);
	char **filePaths = get_file_pathnames_in_directory(directoryPathName);
	int countFiles = count_files_in_directory(directoryPathName);
	//printf("\n\n\n\n\n\n\n\n<<<<<<<<<<<<<<<<<<<<\n %s\n\n>>>>>>>>>>>>>>>>>>>>", directoryName);
	//printf("\n\n Directory Pathname: %s\n Count of Files: %d\n\n\n\nFiles in Directory: \n", directoryPathName, countFiles);
	//printf("\nCount of Files: %-5d", countFiles);
	for(int i = 0; i < countFiles; i++)
	{
		printf("\n    %-2d:       %-5s", i, filePaths[i]);
	}
	//free(directoryPathName); free(directoryName);
	deallocate_memory_char_ptr_ptr(filePaths, countFiles);
}

















