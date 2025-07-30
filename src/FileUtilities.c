//  FileUtilities.c
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02


#include "FileUtilities.h"
#include "CommonDefinitions.h"
#include "GeneralUtilities.h"
#include "StringUtilities.h"
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>






/**
 * determine_file_depth
 *
 * Calculates the depth of a file in a directory structure based on its path.
 * This function counts the number of '/' characters in the file path, assuming a Unix-like file system.
 * This function only works with string-based paths and assumes '/' as the directory separator.
 *
 *
 * @param filePathName A string representing the file path.
 * @return The depth of the file in the directory hierarchy.
 */
int determine_file_depth(const char* filePathName) //Only works for String-based path.
{
	//Simply count the number of '/'.
	int depthCount = 0;
	for(int i = 0; filePathName[i] != '\0'; i++)
	{
		if(filePathName[i] == '/')
		{
			depthCount++;
		}
	}
	return depthCount;
}


/**
 * identify_file_extension
 *
 * Identifies the file extension of a given file path.
 * This function uses the `strrchr` function to find the last occurrence of the '.' character,
 * which is assumed to be the start of the file extension. This function assumes that the file
 * path is a null-terminated string and that the extension is anything following the last '.'
 * in the path.
 *
 * @param filePathName A string representing the file path (i.e., "/home/user/file.txt").
 * @return A pointer to the file extension within the given file path, or NULL if no extension is found.
 */
char *identify_file_extension(const char* filePathName) //Only works for String-based path.
{
	if(count_character_occurrences(filePathName, '.') == 0) // If this pathname does not have a file extension(ex: a directory pathname)
	{
		return "\0";
	}
	
	
	//the strrchr function to get the position of the last occurrence of a character in a string.
	char *fileExtension = strrchr(filePathName, '.');
	
	
	return fileExtension;
}


/**
 * find_file_directory_path
 *
 * Extracts the directory path from a given full file path.
 * This function determines the directory path of a file by counting the characters up to the last '/'
 * character in the path. It dynamically allocates memory to store the directory path.
 * The function assumes Unix-like file paths and uses '/' as the directory separator.
 *
 * @param filePathName A string representing the full file path.
 * @return A dynamically allocated string containing the directory path.
 */
char* find_file_directory_path(const char* filePathName)
{
	int fileDepthLevel = determine_file_depth(filePathName);
	
	//Simply count the number of '/'.
	int directoryDepthLevel = 0;
	
	
	
	int directoryPathNameCharacterCount = 0; //before capturing the actual directory path name that this file is contained inside of, first need to determine how many
											 //characters the directory path name is from the full path name, i.e., need to count all characters except everything after the last '/'.
											 //Then, can use this number to allocate memory for a char* of an appropriate size and it will be straightforward to capture the all characters up
											 //until the cutoff point when the directory that this file is contained in is reached.
	
	
	for(int i = 0; filePathName[i] != '\0' && directoryDepthLevel != fileDepthLevel; i++) //for loop will terminate as soon as the last '/' is reached.
	{
		if(filePathName[i] == '/')
		{
			directoryDepthLevel++;
		}
		
		directoryPathNameCharacterCount++;
	}
	
	
	
	char* directoryPathName;
	directoryPathName = allocate_memory_char_ptr(directoryPathNameCharacterCount);
	for(int i = 0; i < directoryPathNameCharacterCount; i++)
	{
		directoryPathName[i] = filePathName[i];
	}
	
	
	
	return directoryPathName;
}




/**
 * find_directory_name_from_path
 *
 * Extracts and returns the directory name from a given directory path.
 * This function takes a directory path as input and determines the directory name
 * by analyzing the path string. It identifies the last '/' character and extracts
 * the subsequent string as the directory name. This is useful for isolating the
 * immediate directory name from a longer directory path.
 *
 * @param directoryPathName A string representing the directory path.
 * @return A dynamically allocated string containing the directory name.
 */
char *find_directory_name_from_path(const char* directoryPathName)
{
	//Simply count the number of '/'.
	int directoryDepthLevel = 0;
	
	//Count the number of characters in the directory path name.
	size_t directoryPathNameCharacterCount = string_length(directoryPathName);
	
	
	//Count the number of characters in the directory name.
	size_t directoryNameCharacterCount = 0;
	for(size_t i = directoryPathNameCharacterCount - 1; i >= 0; i--)
	{
		if(directoryPathName[i] == '/')
		{
			break; //Stop counting when the last '/' is reached/
		}
		directoryNameCharacterCount++;
	}
	
	
	char* directoryName;
	directoryName = allocate_memory_char_ptr(directoryNameCharacterCount);
	
	for(int i = 0; i < directoryNameCharacterCount; i++)
	{
		// The value in the bracket is the offset from the first character of the directory name to the first character of the directory path name.
		// i is the n_th character of the directory name.
		directoryName[i] = directoryPathName[directoryPathNameCharacterCount - directoryNameCharacterCount + i]; // Compute the offset from filepath[0] to the first character of the directory name.
	}
	
	
	
	
	return directoryName;
}




/**
 * find_name_from_path
 *
 * Extracts the filename from a given full file path.
 * This function isolates the filename from the full path by removing the directory path and the file extension.
 * It dynamically allocates memory for the filename. Assumes Unix-like string-based file paths.
 *
 * @param filePathName A string representing the full file path.
 * @return A dynamically allocated string containing the filename.
 */
char* find_name_from_path(const char* filePathName)
{
	char* fileExtension = identify_file_extension(filePathName);
	
	char* directoryPathName = find_file_directory_path(filePathName);
	
	
	
	/// full path ==
	/// 				directory(.../home/user/...[/filename])   +   file name   +   file extension(.txt, .csv, .dat, etc.)
	
	
	/// //+ file name + file extension(.txt, .csv, .dat, etc.)
	// Determine the length of the file name by subtracting the directory path length and the file extension length from the full path length.
	unsigned long filePathLength = strlen(filePathName);
	unsigned long fileExtensionLength = strlen(fileExtension);
	unsigned long directoryPathLength = strlen(directoryPathName);
	
	
	
	unsigned long fileNameLength = filePathLength - (directoryPathLength + fileExtensionLength);
	
	
	
	char *fileName;
	fileName = allocate_memory_char_ptr(fileNameLength);
	
	for(int i = 0; i < fileNameLength; i++)
	{
		fileName[i] = filePathName[i+directoryPathLength];
	}
	
	
	return fileName;
}




/**
 * count_data_fields
 *
 * Counts the number of data fields in a CSV-style line.
 * This function processes a line from a CSV file and counts the number of fields based on commas.
 * It also handles empty fields (i.e., two consecutive commas) by temporarily replacing them with spaces.
 * The function assumes that the input string is null-terminated and that the fields are
 * separated by commas. It modifies the input string temporarily but restores its original state.
 *
 * @param lineContents A string representing a line from a CSV file.
 * @return The number of data fields in the line.
 */
int count_data_fields(char* lineContents)
{
	for(int i = 0; i < strlen(lineContents)+1;i++)
	{
		if(lineContents[i] == ',' && lineContents[i+1] == ',')
		{
			lineContents[i] = ' ';
		}
	}
	
	
	int count = 0;
	// If the header is not empty, start with count 1 (the first field before any comma).
	if (lineContents && strlen(lineContents) > 0)
	{
		count = 1;
	}
	
	// Iterate over each character in the string.
	for (int i = 0; lineContents[i] != '\0'; i++)
	{
		// Increment count at each comma, which indicates a new field.
		if (lineContents[i] == ',')
		{
			count++;
		}
	}
	return count;
}




/**
 * count_plot_data_fields
 *
 * Counts the number of fields in a data entry for plotting.
 * This function tokenizes a CSV-style line using commas as delimiters and counts the number of tokens.
 * It is specifically intended for counting fields in data entries for plotting purposes. This function
 * creates a copy of the input line for tokenization and frees it after counting. It prints the token count
 * to standard output. Assumes a null-terminated input string and that the data is CSV formatted.
 *
 * @param lineContents A string representing a line from a CSV file.
 * @return The number of fields in the data entry.
 */
int count_plot_data_fields(char* lineContents, const char *delimiter)
{
	// Counts the number of fields in a data entry.
	int count = 0;
	char* dataCopy = strdup(lineContents);
	char* token = strtok(dataCopy, delimiter);
	while (token)
	{
		count++;
		token = strtok(NULL, delimiter);
	}
	free(dataCopy);
	
	printf("\n\ncount_plot_data_fields: %d", count);
	return count;
}




/**
 * count_file_lines
 *
 * Counts the number of lines in a file.
 * This function opens a file specified by its path (filePathName) and counts
 * the number of lines it contains, up to a maximum specified by maxLines.
 * It uses fgets to read each line into a buffer and increments a counter
 * for each line read.
 *
 * @param filePathName A string representing the path to the file.
 * @param maxLines An integer specifying the maximum line length buffer.
 * @return The total number of lines in the file.
 */
int count_file_lines(const char* filePathName, int maxLines)
{
	//printf("\nEntering: 'count_file_lines' function.\n");
	//Open the file at the specified path and ensure file is opened properly.
	FILE *tempFile = fopen(filePathName, "r");
	if (!tempFile)
	{
		perror("\n\nError: Unable to open file for 'count_file_lines'.\n");
		exit(1);
	}
	
	
	char line[maxLines]; // Buffer to hold each line from the file.
						 // Count the number of lines in the file from a starting point numLinesToSkip away from the first line.
	int count = 0;
	while((fgets(line, sizeof(line), tempFile)))
	{
		count++;
	}
	
	
	if (count >= maxLines)
	{
		printf("\n\nError: The number of lines in the file '%s' exceeded the expected number of maximum lines, %zu 'count_file_lines'.\n", filePathName, MAX_NUM_FILE_LINES);
		perror("\n\nError: File  'count_file_lines'.\n");
		exit(1);
	}
	fclose(tempFile);
	return count;
}




/**
 * count_file_lines_characters
 *
 * Counts the characters in each line of a file.
 * This function opens a file specified by filePathName and counts the number
 * of characters in each line, starting from a specified line (startLine) and
 * up to a total number of lines (lineCount). It uses getline for dynamic allocation
 * and handling of varying line lengths. The newline character is excluded from the count.
 *
 * @param filePathName A string representing the path to the file.
 * @param lineCount An integer specifying the number of lines to process.
 * @return A pointer to an array of integers, each representing the character count of a line.
 */
int* count_file_lines_characters(const char* filePathName, int lineCount)
{
	//Open the file at the specified path and ensure file is opened properly.
	FILE *tempFile = fopen(filePathName, "r");
	if (!tempFile)
	{
		perror("\n\nError: Unable to open file for 'count_file_lines_characters'.\n");
		exit(1);
	}
	
	
	int *fileCharCounts = allocate_memory_int_ptr(lineCount); //allocate memory to store the number of characters for each line.
															  // Initialize the array with zeros.
	for (int i = 0; i < lineCount; i++)
	{
		fileCharCounts[i] = 0;
	}
	
	
	// Declare variables for
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int currentLine = 0;
	
	// Read each line from the file.
	while ((read = getline(&line, &len, tempFile)) != -1 && currentLine < lineCount)
	{
		// Count characters in the line, excluding the newline character.
		for (int i = 0; i < read; i++)
		{
			if (line[i] != '\n')
			{
				fileCharCounts[currentLine]++;
			}
		}
		currentLine++;
	}
	
	free(line);
	fclose(tempFile);
	
	return fileCharCounts;
}



int determine_maximum_line_length(const char* filePathName, int lineCount) // Determines the maximum line length in a file using count_file_lines_characters
{
	int *charCounts = count_file_lines_characters(filePathName, lineCount);
	int maxLineLength = max_element(charCounts, lineCount);
	free(charCounts);
	return maxLineLength;
}

/**
 * count_characters_in_file_lines_range
 *
 * Counts the characters in each line of a file.
 * This function opens a file specified by filePathName and counts the number
 * of characters in each line, starting from a specified line (startLine) and
 * up to a total number of lines (lineCount). It uses getline for dynamic allocation
 * and handling of varying line lengths. The newline character is excluded from the count.
 *
 * @param filePathName A string representing the path to the file.
 * @param lineCount An integer specifying the number of lines to process.
 * @param startLine An integer specifying the line number to start counting from.
 * @return A pointer to an array of integers, each representing the character count of a line.
 */
int* count_characters_in_file_lines_range(const char* filePathName, int lineCount, int startLine)
{
	//Open the file at the specified path and ensure file is opened properly.
	FILE *tempFile = fopen(filePathName, "r");
	if (!tempFile)
	{
		perror("\n\nError: Unable to open file for 'count_file_lines_characters'.\n");
		exit(1);
	}
	
	
	
	int *fileCharCounts = allocate_memory_int_ptr(lineCount); //allocate memory to store the number of characters for each line.
	
	
	// Initialize the array with zeros.
	for (int i = 0; i < lineCount; i++)
	{
		fileCharCounts[i] = 0;
	}
	
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int currentLine = 0;
	
	// Read each line from the file.
	while ((read = getline(&line, &len, tempFile)) != -1 && currentLine < lineCount)
	{
		// If it's a line that needs to be skipped, just increment the line index.
		if (currentLine < startLine)
		{
			currentLine++;
			continue;
		}
		
		// Count characters in the line, excluding the newline character.
		for (int i = 0; i < read; i++)
		{
			if (line[i] != '\n')
			{
				fileCharCounts[currentLine]++;
			}
		}
		currentLine++;
	}
	
	free(line);
	fclose(tempFile);
	
	return fileCharCounts;
}








/**
 * read_file_contents
 *
 * Reads the contents of a file into an array of strings.
 * This function opens a file and reads its contents line by line, counting the characters in each line,
 * and then reading each line into a dynamically allocated array of strings (char pointers).
 * Each line is stored in the array up to the specified lineCount.
 *
 * @param filePathName A string representing the path of the file to be read.
 * @param lineCount An integer specifying the number of lines to read from the file.
 * @return A pointer to an array of strings, each string holding the content of a line.
 */
char** read_file_contents(const char* filePathName, int lineCount)
{
	// Checking if the filePathName is NULL or if lineCount is less than or equal to 0.
	// If so, return NULL as it's not possible to read the file contents under these conditions.
	if (filePathName == NULL || lineCount <= 0)
	{
		perror("\n\nError:filePathName is NULL or if lineCount is less than or equal to 0 in 'read_file_contents'.");
		return NULL;
	}
	
	
	//Open the file at the specified path and ensure file is opened properly.
	FILE *tempFile = fopen(filePathName, "r");
	if (!tempFile)
	{
		perror("\n\nError: Unable to open file for 'read_file_contents'.");
		exit(1);
	}
	
	
	// First, count the characters in each line.
	int *charCounts = count_file_lines_characters(filePathName, lineCount);
	size_t maxLineLength = max_element_int(charCounts, lineCount);
	
	
	// Allocate memory for the array of strings.
	char **fileContents = allocate_memory_char_ptr_ptr(maxLineLength+1, lineCount);
	
	size_t len = 0; // Buffer length for getline
	size_t read; //
	char *line = NULL;
	int currentLine = 0;
	
	
	
	// Read the file line by line and store each line in the corresponding string.
	while ((read = getline(&line, &len, tempFile)) != -1 && currentLine < lineCount)
	{
		// Allocate memory for the current line in the array.
		fileContents[currentLine] = allocate_memory_char_ptr((charCounts[currentLine] + 1));
		
		
		
		// Copy the contents of the line into the allocated memory
		strncpy(fileContents[currentLine], line, charCounts[currentLine]);
		fileContents[currentLine][charCounts[currentLine]] = '\0'; // Null-terminate the string
		
		currentLine++;
	}
	free(line);
	free(charCounts);
	fclose(tempFile);
	
	
	return fileContents;
}



/**
 * parse_file_contents
 *
 * Reads the contents of a file into an array of strings.
 * This function opens a file and reads its contents line by line, counting the characters in each line,
 * and then reading each line into a dynamically allocated array of strings (char pointers).
 * Each line is stored in the array up to the specified lineCount.
 *
 * @param filePathName A string representing the path of the file to be read.
 * @param lineCount An integer specifying the number of lines to read from the file.
 * @return A pointer to an array of strings, each string holding the content of a line.
 */
char** parse_file_contents(const char* filePathName, int lineCount) // Conditionally reads the contents of a file into a string array
{
	// Checking if the filePathName is NULL or if lineCount is less than or equal to 0.
	// If so, return NULL as it's not possible to read the file contents under these conditions.
	if (filePathName == NULL || lineCount <= 0)
	{
		perror("\n\nError:filePathName is NULL or if lineCount is less than or equal to 0 in 'read_file_contents'.");
		return NULL;
	}
	
	
	//Open the file at the specified path and ensure file is opened properly.
	FILE *tempFile = fopen(filePathName, "r");
	if (!tempFile)
	{
		perror("\n\nError: Unable to open file for 'read_file_contents'.");
		exit(1);
	}
	
	
	// First, count the characters in each line.
	int *charCounts = count_file_lines_characters(filePathName, lineCount);
	size_t maxLineLength = max_element_int(charCounts, lineCount);
	
	// Allocate memory for the array of strings.
	char **fileContents = allocate_memory_char_ptr_ptr(maxLineLength+1, lineCount);
	
	size_t len = 0; // Buffer length for getline
	size_t read; //
	char *line = NULL;
	int currentLine = 0;
	
	
	
	// Read the file line by line and store each line in the corresponding string.
	while ((read = getline(&line, &len, tempFile)) != -1 && currentLine < lineCount)
	{
		// Allocate memory for the current line in the array.
		fileContents[currentLine] = allocate_memory_char_ptr((charCounts[currentLine] + 1));
		
		
		
		// Copy the contents of the line into the allocated memory
		strncpy(fileContents[currentLine], line, charCounts[currentLine]);
		fileContents[currentLine][charCounts[currentLine]] = '\0'; // Null-terminate the string
		
		currentLine++;
	}
	free(line);
	free(charCounts);
	fclose(tempFile);
	
	
	
	// Replace double commas with a space to handle missing values.
	// This is meant to ensure that empty fields are accounted for.
	for (int i = 0; i < lineCount; i++)
	{
		char *tempData = fileContents[i];
		int extraSpaceNeeded = 0;
		
		// For each line, scan through the string and count how many times two commas appear consecutively. This will determine the additional space needed.
		for (int j = 0; tempData[j] != '\0'; j++)
		{
			if (tempData[j] == ',' && tempData[j + 1] == ',')
			{
				extraSpaceNeeded += strlen(",0.0,") - 1; // Subtract 1 because we replace two commas.
			}
		}
		
		
		// Allocate new memory for the modified string
		char *modifiedData = allocate_memory_char_ptr(strlen(tempData) + extraSpaceNeeded + 1);
		
		
		// Copy characters from the original string to the new string. When two consecutive commas are encountered, insert ", 0.0" and adjust the index accordingly.
		for (int j = 0, k = 0; tempData[j] != '\0'; j++, k++)
		{
			if (tempData[j] == ',' && tempData[j + 1] == ',')
			{
				strcpy(&modifiedData[k], ",0.0,");
				k += strlen(",0.0,") - 1; // Adjust index.
				j++; // Skip the next comma in the original string.
			}
			else
			{
				modifiedData[k] = tempData[j];
			}
		}
		modifiedData[strlen(tempData) + extraSpaceNeeded] = '\0'; // Null-terminate the new string.
		
		// Replace the original string with the modified string.
		free(fileContents[i]); // Free the memory allocated for the original string.
		fileContents[i] = modifiedData;
	}
	
	
	return fileContents;
}



/**
 * write_file_contents
 *
 * Writes an array of strings to a file.
 * This function opens a file for writing and writes each string from
 * the provided array (fileContents) to the file followed by a newline character.
 * If writing fails, an error is reported and the program exits.
 *
 * @note The function assumes that the array is NULL-terminated, and that the individual
 *       strings do not contain newline characters at the end.
 *
 * @param filename A string representing the name of the file to write to.
 * @param fileContents A pointer to an array of strings to be written.
 */
void write_file_contents(const char *filename, char **fileContents)
{
	// Open the file for append and update
	FILE *file = fopen(filename, "a+");
	
	if (file == NULL)
	{
		perror("\n\nError opening file for writing in 'write_file_contents'.");
		exit(1);
	}
	
	// Iterate through the array of strings
	for (int i = 0; fileContents[i] != NULL; i++)
	{
		// Write each string to the file
		if (fputs(fileContents[i], file) == EOF)
		{
			perror("\n\nError writing to file in 'write_file_contents'.");
			fclose(file);
			exit(1);
		}
		
		// Add a newline character to ensure the strings contain one
		if (fputc('\n', file) == EOF)
		{
			perror("\n\nError writing to file in 'write_file_contents'.");
			fclose(file);
			exit(1);
		}
	}
	
	fclose(file);
}


/**
 * write_file_numeric_data
 *
 * Writes an array of double values to a file.
 * This function opens a file in append mode and writes each double value from the provided
 * data array to the file, each on a new line. The precision of the double values
 * is maintained up to 17 significant digits. If writing fails, an error is reported
 * and the program exits.
 *
 * @param filename A string representing the name of the file to write to.
 * @param data A pointer to an array of double values.
 * @param countDataEntries An integer specifying the number of entries in the data array.
 */
void write_file_numeric_data(const char *filename, double *data, int countDataEntries, const char *dataFieldName)
{
	// Open the file for append and update
	FILE *file = fopen(filename, "a+");
	if (file == NULL)
	{
		perror("\n\nError opening file for writing in 'write_file_numeric_data'.");
		exit(1);
	}
	
	
	
	
	/// Write the name of the data field to the file... the two lines of code below are analogous, although fprintf is primarily used for writing formatted data to a file, whereas fputs is more generally used to write string(s) to a file.
	fprintf(file, "%s\n", dataFieldName);
	//fputs(dataFieldName, file); //fputc('\n', file);
	
	
	
	
	// Iterate through the array of strings
	for (int i = 1; i < countDataEntries; i++)
	{
		// Write each data entry to the file
		if (fprintf(file, "%.17g", data[i]) < 0)
		{
			perror("\n\nError writing to file in 'write_file_numeric_data'.");
			fclose(file);
			exit(1);
		}
		// Optionally add a newline character if the strings do not already contain one
		if (fputc('\n', file) == EOF)
		{
			perror("\n\nError writing to file in 'write_file_numeric_data'.");
			fclose(file);
			exit(1);
		}
		
		
		// Unnecessary, keeping for reference.
		//if(i == countDataEntries-1)
		//{
		//fputc('\n', file); //
		//}
	}
	
	
	fclose(file);
}








/**
 * generate_merged_filename
 *
 * Generates a merged filename based on two input file paths.
 * The merged filename is constructed by combining the base names of the input files,
 * appending "_merged" to the combined base names, and then prepending the directory path.
 *
 * @param filePath1 The path to the first file.
 * @param filePath2 The path to the second file.
 * @return A dynamically allocated string containing the new filename.
 */
char* generate_merged_filename(const char* filePath1, const char* filePath2)
{
	// Extract the base names of the files, without extension or path, to avoid directory paths in the merged filename.
	char* baseName1 = find_name_from_path(filePath1);
	char* baseName2 = find_name_from_path(filePath2);
	
	
	// Calculate the length needed for the new filename
	size_t mergedLength = string_length(baseName1) + string_length(baseName2) + string_length("_merged") + 2;  // +2 for underscore and null terminator
	char* mergedFilename = allocate_memory_char_ptr(mergedLength);
	if (mergedFilename == NULL)
	{
		perror("Memory allocation failed for merged filename");
		return NULL;
	}
	sprintf(mergedFilename, "%s_%s_merged", baseName1, baseName2); // Construct the new filename
	
	
	char *directoryPathName1 = find_file_directory_path(filePath1);
	char *directoryPathName2 = find_file_directory_path(filePath2);
	size_t characterCountMergedPathname = mergedLength;
	characterCountMergedPathname += (directoryPathName1 == directoryPathName2) ? string_length(directoryPathName1) : maximum(string_length(directoryPathName1), string_length(directoryPathName2)); // If they are contained within the same enclosing directory, append the first directory pathname, otherwise, use the longer pathname.
	
	
	char *mergedFilePathName = allocate_memory_char_ptr(characterCountMergedPathname); // Allocate memory for the merged file path name
	sprintf(mergedFilePathName, "%s%s", directoryPathName1, mergedFilename); // Construct the full path name for the merged file
	
	
	/// Cleanup memory
	free(baseName1);
	free(baseName2);
	free(mergedFilename);
	free(directoryPathName1);
	free(directoryPathName2);
	
	return mergedFilePathName; // Return the full path name of the merged file
}


/**
 * merge_two_files
 *
 * This function merges the contents of two files into a new file.
 * It reads the contents of the two input files, combines them, and writes the merged content to a new file.
 * The new file's name is generated based on the names of the input files.
 *
 * @param filePath1 The path to the first file.
 * @param filePath2 The path to the second file.
 * @return The path to the merged file.
 */
char* merge_two_files(const char* filePath1, const char* filePath2)
{
	// Count lines in each file
	int lineCount1 = count_file_lines(filePath1, MAX_NUM_FILE_LINES);
	int lineCount2 = count_file_lines(filePath2, MAX_NUM_FILE_LINES);
	
	char** contents1 = read_file_contents(filePath1, lineCount1);
	char** contents2 = read_file_contents(filePath2, lineCount2);
	
	
	// Create filename based on input filenames, "merged" string literal, and enclosing directory of the files
	char* mergedFilename = generate_merged_filename(filePath1, filePath2);  // Combines the base names of the files(extracted from full path) then appends this to the path to the enclosing directory of the files
	
	
	// Allocate memory for the merged contents array of strings
	char** mergedContents = allocate_memory_char_ptr_ptr(MAX_STRING_SIZE, lineCount1 + lineCount2 + 1);  // +1 for safety/null termination
	
	
	
	// Copy contents from both files into the merged array
	int i, j;
	for (i = 0; i < lineCount1; i++)
	{
		mergedContents[i] = duplicate_string(contents1[i]); // Copy contents from the first file
	}
	for (j = 0; j < lineCount2; j++, i++)
	{
		mergedContents[i] = duplicate_string(contents2[j]); // Copy contents from the second file
	}
	mergedContents[i] = NULL;  // Null-terminate the array
	
	
	
	// Write the merged contents to the new file
	write_file_contents(mergedFilename, mergedContents);
	
	// Cleanup memory
	deallocate_memory_char_ptr_ptr(contents1, lineCount1);
	deallocate_memory_char_ptr_ptr(contents2, lineCount2);
	deallocate_memory_char_ptr_ptr(mergedContents, lineCount1 + lineCount2);
	
	return mergedFilename; // Return the path to the merged file
}












/**
 * create_directory
 *
 * Creates a directory for storing dataset files. This function first determines
 * the directory path based on the provided file path, then checks if a directory
 * with the intended name already exists. If it does, the existing directory is
 * deleted to ensure a fresh start. Finally, a new directory is created.
 *
 * @param filePathName The path of the file for which the directory is to be created.
 * @return A pointer to a string containing the path of the newly created dataset directory.
 */
char *create_directory(const char *filePathName, const char *directoryName)
{
	// Extracts the directory path and dataset file name from the provided file path.
	char *directoryPathName = find_file_directory_path(filePathName);
	char *dataSetFileName = find_name_from_path(filePathName);
	dataSetFileName = combine_strings(dataSetFileName, directoryName);
	char *dataDirectoryPathName = combine_strings(directoryPathName, dataSetFileName);
	
	
	// Check if the directory already exists, and delete it if it does
	struct stat st;
	if (stat(dataDirectoryPathName, &st) == 0)
	{
		// If it exists, use delete_directory to remove it along with its contents.
		delete_directory(dataDirectoryPathName);
	}
	
	
	// Create the new directory with read, write, and execute permissions for the owner.
	int status = mkdir(dataDirectoryPathName, 0700);
	if (status == 0)
	{
		//printf("\n\nDirectory created in: '%s'.", dataDirectoryPathName); // Commented out the message for sucessful directory creation, annoying
	}
	else
	{
		perror("\n\nError: mkdir failed in 'directoryPathName'.");
	}
	
	
	free(directoryPathName);
	free(dataSetFileName);
	return dataDirectoryPathName;
}




/**
 * delete_directory
 *
 * Recursively deletes a directory and all its contents. This function opens the
 * specified directory, iterates through all its entries, and deletes each one.
 * If an entry is a subdirectory, the function calls itself recursively.
 * It is designed to handle normal files and directories but not special file types
 * like symbolic links.
 *
 * @param path The path of the directory to be deleted.
 * @return 0 on successful deletion, -1 on failure with errno set appropriately.
 */
int delete_directory(const char *path)
{
	// Open the directory and prepare for iteration over its contents.
	DIR *dir = opendir(path);
	size_t path_len = strlen(path);
	int r = -1;
	
	
	if (dir)
	{
		struct dirent *p;
		r = 0;
		
		
		// Iterate over each entry in the directory.
		while (!r && (p=readdir(dir)))
		{
			int r2 = -1;
			char *buf;
			size_t len;
			
			// Skip the names "." and ".." as we don't want to recurse on them.
			if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
			{
				continue;
			}
			
			
			// Prepare the full path of the entry.
			len = path_len + strlen(p->d_name) + 2;
			buf = malloc(len);
			
			if (buf)
			{
				struct stat statbuf;
				snprintf(buf, len, "%s/%s", path, p->d_name);
				
				// Delete the directory entry based on its type (file or directory).
				if (!stat(buf, &statbuf))
				{
					if (S_ISDIR(statbuf.st_mode))
					{
						r2 = delete_directory(buf);
					}
					else
					{
						r2 = remove(buf);
					}
				}
				free(buf);
			}
			r = r2;
		}
		closedir(dir);
	}
	
	
	// Remove the directory itself after its contents have been cleared.
	if (!r)
	{
		r = rmdir(path);
	}
	
	
	return r;
}











/**
 * count_files_in_directory
 *
 * This function counts the number of files within a specified directory.
 * It opens the directory, iterates over each entry, and increments a counter
 * for each valid file, excluding the current (".") and parent ("..") directory entries.
 * The function returns the total count of files found, or -1 if the directory cannot be opened.
 *
 * @param directoryPath: A string representing the path of the directory to scan.
 * @return An integer representing the number of files in the directory. Returns -1 if the directory cannot be opened.
 */
int count_files_in_directory(const char *directoryPath)
{
	DIR *dir;
	struct dirent *entry;
	int fileCount = 0;
	
	// Open the directory
	if ((dir = opendir(directoryPath)) == NULL)
	{
		perror("Unable to open directory");
		return -1; // Return -1 to indicate an error opening
	}
	
	// Iterate over each file in the directory
	while ((entry = readdir(dir)) != NULL)
	{
		// Skip "." and ".." entries
		if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
		{
			continue;
		}
		
		// Increment the file count
		fileCount++;
	}
	
	// Close the directory
	closedir(dir);
	
	return fileCount;
}







/**
 * get_file_pathnames_in_directory
 *
 * This function retrieves the full pathnames of files within a specified directory.
 * It opens the directory, iterates over each entry, and constructs the full path for each file,
 * excluding the current (".") and parent ("..") directory entries. It stores each valid file path
 * in a dynamically allocated array, which is then returned. The function handles memory allocation
 * for the array of file paths, and ensures directory resources are appropriately released upon completion.
 *
 * @param directoryPath: A string representing the path of the directory to scan.
 * @return An array of strings (char **) where each entry is the full pathname of a file in the directory.
 *         Returns NULL if the directory cannot be opened or memory allocation fails.
 */
char **get_file_pathnames_in_directory(const char *directoryPath)
{
	DIR *dir; // Pointer to a directory stream
	struct dirent *entry; // For storing information about each directory entry
	char **fileList = NULL; // Array to hold the paths of all files found in the directory
	int fileCount = 0; // Number of files found in the directory, as well as their indexing in the fileList array
	
	// Open the directory
	if ((dir = opendir(directoryPath)) == NULL)
	{
		perror("Unable to open directory");
		return NULL;
	}
	
	
	// Iterate over each file in the directory
	while ((entry = readdir(dir)) != NULL)
	{
		// Skip "." and ".." entries, which refer to the current and parent directories, respectively
		if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
		{
			continue;
		}
		
		
		// Construct the full path of the file
		char *fullPath = combine_strings(directoryPath, "/");
		fullPath = combine_strings(fullPath, entry->d_name);
		
		
		// Expand the fileList array to store the new file path
		fileList = realloc(fileList, sizeof(char *) * (fileCount + 1));
		if (fileList == NULL)
		{
			perror("Unable to allocate memory");
			closedir(dir);
			return NULL;
		}
		
		
		
		// Store the file path in the array
		fileList[fileCount] = strdup(fullPath);
		fileCount++;
		
		// Free the fullPath string
		free(fullPath);
	}
	
	
	// Close the directory
	closedir(dir);
	
	
	
	
	// print the fileList
	//printf("\n\n\n\n File List in Directory: %s\n", directoryPath);
	//for (int i = 0; i < fileCount; i++)
	//{
	//	printf("	File %d: %s\n", i, fileList[i]);
	//}
	//printf("\n\n\n\n");
	
	
	
	
	return fileList;
}


char **get_file_names_in_directory(const char *directoryPath)
{
	// Use get_file_pathnames_in_directory to get pathnames, then read names from pathnames
	char **filePaths = get_file_pathnames_in_directory(directoryPath);
	if (filePaths == NULL)
	{
		perror("Unable to retrieve file paths");
		return NULL;
	}
	
	int fileCount = count_files_in_directory(directoryPath);
	
	// Allocate memory for the array of file names
	char **fileNames = (char**)malloc(fileCount * sizeof(char*));// = allocate_memory_char_ptr_ptr(MAX_STRING_SIZE, count_files_in_directory(directoryPath) + 1); // +1 for safety/null termination
	
	
	// Extract the file names from the full pathnames
	for (int i = 0; i < fileCount && filePaths[i] != NULL; i++)
	{
		char* fName = find_name_from_path(filePaths[i]);
		size_t fNameLength = string_length(fName);
		fileNames[i] = allocate_memory_char_ptr(fNameLength);
		fileNames[i] = fName;
		free(fName);
	}
	
	// Free the file paths
	deallocate_memory_char_ptr_ptr(filePaths, fileCount);
	
	
	
	return fileNames;
}



FileProperties extract_file_properties(const char *filePathName)
{
	FileProperties fileProperties;
	
	fileProperties.fileName = find_name_from_path(filePathName);
	fileProperties.filePathName = filePathName;
	fileProperties.fileExtension = identify_file_extension(filePathName);
	fileProperties.fileLineCount = count_file_lines(filePathName, MAX_NUM_FILE_LINES);
	return fileProperties;
}


DirectoryProperties extract_directory_properties(const char *directoryPathName)
{
	DirectoryProperties fileDirectory;
	
	fileDirectory.directoryPathName = directoryPathName;
	fileDirectory.directoryName = find_directory_name_from_path(fileDirectory.directoryPathName);
	int fileCount = count_files_in_directory(fileDirectory.directoryPathName);
	fileDirectory.fileCount = fileCount;
	
	
	
	char **filePathNames = get_file_pathnames_in_directory(fileDirectory.directoryPathName);
	FileProperties *fileProerties = (FileProperties*)malloc(fileCount * sizeof(FileProperties));
	for (int i = 0; i < fileCount; i++)
	{
		fileProerties[i] = extract_file_properties(filePathNames[i]);
	}
	fileDirectory.fileProperties = fileProerties;
	
	
	return fileDirectory;
}


DirectoryProperties extract_directory_properties_from_file(const char *filePathName)
{
	DirectoryProperties fileDirectory;
	
	fileDirectory.directoryPathName = find_file_directory_path(filePathName);
	fileDirectory.directoryName = find_directory_name_from_path(fileDirectory.directoryPathName);
	int fileCount = count_files_in_directory(fileDirectory.directoryPathName);
	fileDirectory.fileCount = fileCount;
	
	
	
	
	char **filePathNames = get_file_pathnames_in_directory(fileDirectory.directoryPathName);
	FileProperties *fileProerties = (FileProperties*)malloc(fileCount * sizeof(FileProperties));
	for (int i = 0; i < fileCount; i++)
	{
		fileProerties[i] = extract_file_properties(filePathNames[i]);
	}
	fileDirectory.fileProperties = fileProerties;
	
	
	
	
	return fileDirectory;
}



void print_directory_propertiees(DirectoryProperties fileDirectory)
{
	printf("\n\n\n\n\n\n\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
	printf("\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	printf("\n Dataset Directory %s\n", fileDirectory.directoryName);
	printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
	printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
	printf("\n|\n| Directory Pathname: %-5s", fileDirectory.directoryPathName);
	printf("\n| Count of Files: %-5d", fileDirectory.fileCount);
	printf("\n|\n| Files in Directory: \n|");
	for(int i = 0; i < fileDirectory.fileCount; i++)
	{
		printf("\n| 	%d:       %-5s", i, fileDirectory.fileProperties[i].fileName);
		printf("\n| 	%d:       %-5s", i, fileDirectory.fileProperties[i].filePathName);
		printf("\n| 	%d:       %-5s", i, fileDirectory.fileProperties[i].fileExtension);
		printf("\n| 	%d:       %-5d\n|", i, fileDirectory.fileProperties[i].fileLineCount);
		
	}
	
	
	//free(writtenDataSetDirectory);
	//deallocate_memory_char_ptr_ptr(writtenDataSetFilePaths, countWrittenDataSetFiles);
	
}

