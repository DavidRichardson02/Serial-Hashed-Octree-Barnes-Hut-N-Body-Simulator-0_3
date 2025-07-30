//  DataExtraction.c
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02


#include "DataExtraction.h"
#include "CommonDefinitions.h"
#include "GeneralUtilities.h"
#include "StringUtilities.h"
#include "FileUtilities.h"
#include <ctype.h>






/**
 * parse_entire_file
 *
 * Function to parse the entire contents of a CSV file.
 * It categorizes and sorts all the values of each one of the dataset's fields and stores them in their own individual char** as a list of values.
 * The function takes in the file contents, line count, field count, and delimiter as parameters. It returns an array of arrays of strings,
 * where each of the ith subarrays of the highest level array is a representation of the ith data set fields and all the entries within it.
 *
 * @param fileContents Constant character pointer to the string containing the data.
 * @param lineCount Constant character pointer to the string containing the  data.
 * @param fieldCount Pointer to the integer value representing the number of paramaters a data set has.
 * @return An array of arrays of strings, where each of the ith subarrays of the highest level array is a representation of the ith data set fields and all the entries within it, while the
 */
char*** parse_entire_file(char **fileContents, int lineCount, int* fieldCount, const char *delimiter)
{
	// Capture the header of the data set for plotting
	char **fileHeader = capture_data_set_header_for_plotting(fileContents[0], fileContents, delimiter);
	int fieldCountCopy = count_data_fields(fileContents[0]);
	
	
	
	/// Capture the Data Field Names in a Single String Without the typeDelimiters
	const char *fieldNames[fieldCountCopy]; //The string of field names as they appear in the first element of the fileContents array of strings
	for (int i = 0; i < fieldCountCopy; i++)
	{
		const char *typeDelimiter = identify_delimiter(fileHeader, fieldCountCopy);
		
		// Split the token by the type delimiter to get the name
		char* typeDelimiterToken = strtok(fileHeader[i], typeDelimiter);
		if (typeDelimiterToken)
		{
			// Duplicate the token to the field names
			fieldNames[i] = strdup(typeDelimiterToken);
		}
	}
	
	
	// Allocate array of string arrays for each member
	char*** separatedData = malloc(fieldCountCopy * sizeof(char**));
	
	// Initialize arrays for each member
	for (int i = 0; i < fieldCountCopy; i++)
	{
		// Allocate memory for each member
		separatedData[i] = malloc((lineCount) * sizeof(char*)); // No need for +1, skipping header
		if (!separatedData[i])
		{
			// Free allocated memory in case of failure
			while (i--)
			{
				free(separatedData[i]);
			}
			free(separatedData);
			*fieldCount = 0;
			return NULL;
		}
		
		
		// Copy the key name
		separatedData[i][0] = strdup(fieldNames[i]);
		// Free fieldNames elements
		free((void*)fieldNames[i]);
	}
	
	
	// Parse the data
	for (int line = 1; line < lineCount; line++) // Start from line 1 to skip header
	{
		// Duplicate the data
		char* dataCopy = strdup(fileContents[line]);
		
		// Tokenize the data
		char* token = strtok(dataCopy, delimiter);
		for (int i = 0; i < fieldCountCopy && token; i++)
		{
			// Store the token in the separated data
			separatedData[i][line] = strdup(token);
			
			// Get the next token
			token = strtok(NULL, delimiter);
		}
		// Free the duplicated data
		free(dataCopy);
	}
	
	// Update the field count
	*fieldCount = fieldCountCopy;
	return separatedData;
}




/**
 * count_missing_values
 *
 * Counts the number of missing or incorrectly formatted values per field in the dataset.
 * Determines the data types of each field based on the second line of the file contents.
 * It iterates over each line in the dataset, tokenizes the line based on the delimiter, and checks each token.
 * If a token is a hyphen and the corresponding field is nonnumeric, it is considered a missing value and the count for that field is incremented.
 * If the determined data type of a token does not match the expected data type of the corresponding field, it is considered incorrectly formatted and the count for that field is incremented.
 * After all lines have been processed, it frees the memory allocated for the data types and returns the count of missing or incorrectly formatted values for each field.
 *
 *
 * @param fileContents Array of strings where each string is a line in the dataset.
 * @param lineCount Number of lines in the dataset.
 * @param fieldCount Number of fields in each line.
 * @param delimiter Delimiter used in the dataset.
 * @param headerLine Header line of the dataset to determine the data types of fields.
 * @return An array of integers where each integer represents the count of missing or incorrectly formatted values for a corresponding field.
 */
int* count_missing_values(char** fileContents, int lineCount, int fieldCount, const char *delimiter, const char *headerLine)
{
	if (!fileContents || lineCount <= 0 || fieldCount <= 0)
	{
		return NULL;
	}
	
	int *missingDataCount = calloc(fieldCount, sizeof(int)); // Initialize counts to zero
	if (!missingDataCount)
	{
		perror("Memory allocation failed for missingDataCount");
		return NULL;
	}
	
	// Analyze the data types based on the header line (assuming first line is the header)
	char **dataTypes = determine_data_entry_types(fileContents[1], fieldCount, delimiter);
	
	/// Iterate over each line in the dataset, tokenizing the line based on the delimiter, and checking each token.
	for (int i = 1; i < lineCount; i++) // Start from 1 to skip the header
	{
		char *dataCopy = strdup(fileContents[i]);
		char *token = strtok(dataCopy, delimiter);
		int fieldIndex = 0;
		
		
		/// Check each token
		while (token && fieldIndex < fieldCount)
		{
			while (isspace((unsigned char)*token)) token++; // Trim leading spaces from the token
			
			const char *typeDataEntry = determine_string_representation_type(token);
			
			
			
			/// Here the delimited value is nothing but a hyphen "-". This field is assumed to be a missing value(either N/A or not available)
			/// Furthermore, to ensure the value is not being used as a minus sign, a check is made to ensure the token is nonnumeric
			if (strcmp(token, "-") == 0 && strcmp(dataTypes[fieldIndex], "nonnumeric") == 0 )
			{
				missingDataCount[fieldIndex]++;  // This field is assumed to be a missing value(either N/A or not available)
			}
			else if (strcmp(typeDataEntry, dataTypes[fieldIndex]) != 0) //The determined data type does not match the expected data type
			{
				missingDataCount[fieldIndex]++;
			}
			
			
			token = strtok(NULL, delimiter);
			fieldIndex++;
		}
		
		free(dataCopy);
	}
	
	for (int i = 0; i < fieldCount; i++)
	{
		free(dataTypes[i]);
	}
	free(dataTypes);
	
	return missingDataCount; // Return the count of missing or incorrectly formatted values for each field
}




/**
 * capture_data_set_header_for_plotting
 *
 * Processes the header line of a dataset to pair each field with its corresponding data type
 * determined from the dataset's contents. Captures and formats the header information for
 * preparing data for plotting by knowing the type of each field.
 *
 *
 * @param headerContents A string containing the header line of the dataset.
 * @param fileContents An array of strings, each representing a line in the dataset.
 * @param delimiter A character used as a delimiter in the dataset.
 * @return An array of strings where each string is a field name from the header paired with its data type.
 */
char** capture_data_set_header_for_plotting(char* headerContents, char** fileContents, const char *delimiter)
{
	// Counts the number of fields in the header.
	int fieldCount = count_data_fields(headerContents);
	
	// Determines the data types of each field based on the second line of file contents.
	char **dataTypes = determine_data_entry_types(fileContents[1], fieldCount, delimiter);
	
	// Calculates total length needed for formatted data (field name + data type).
	int formatDataLength = 0;
	for (int i = 0; i < fieldCount; i++)
	{
		formatDataLength += strlen(dataTypes[i]) + 1;
	}
	
	
	// Allocates memory for storing formatted field names and types, indexing each field to capture all plotting relevant fields from the contents file
	char **fieldIndexing = allocate_memory_char_ptr_ptr(formatDataLength, fieldCount);
	
	
	
	// Duplicates the header string for tokenization.
	char* headerCopy = strdup(headerContents);
	
	
	//int primaryDelimiterCount = 1;
	//char **fileHeader = capture_data_set_header_for_plotting(headerCopy,  fileContents, delimiter);  // Capture header information for plotting.
	
	//const char *typeDelimiter = identify_delimiter_recursive(fileHeader, fieldCount, &delimiter, &primaryDelimiterCount);
	//printf("\n\ntypeDelimiter: %s\n", typeDelimiter);
	
	// Tokenize the header and pair each field with its data type.
	char* token = strtok(headerCopy, delimiter);
	int formatIndex = 0;
	while(token && formatIndex < fieldCount)
	{
		strcpy(fieldIndexing[formatIndex], token); // Copy the token into fieldIndexing
		strcat(fieldIndexing[formatIndex], ":");   // Append the type delimiter.
		strcat(fieldIndexing[formatIndex], dataTypes[formatIndex]); // Append the data type.
		
		
		// Get the next token
		formatIndex++;
		token = strtok(NULL, delimiter);
	}
	free(headerCopy);
	
	
	
	//print_string_array(fieldIndexing, formatIndex, "capture_data_set_header_for_plotting");
	
	return fieldIndexing; // Returns formatted field names with data types.
}






/**
 * determine_data_entry_types
 *
 * Determines the data types of fields in a data entry.
 * Analyzes a line from the dataset to infer the data type of each field.
 * It helps in understanding the kind of data each field in the dataset holds.
 *
 * @param dataEntry A string containing a line from the dataset.
 * @param fieldCount The number of fields in the data entry.
 * @param delimiter The character used as a delimiter in the dataset.
 * @return An array of strings where each string represents the inferred data type of a field.
 */
char** determine_data_entry_types(const char *dataEntry, int fieldCount, const char *delimiter)
{
	// Increase the buffer size to accommodate additional characters
	size_t characterCount = (strlen(dataEntry) + 1) * 3;
	
	// Allocate memory for storing data types of each field.
	char** dataFields = allocate_memory_char_ptr_ptr(characterCount, fieldCount);
	
	// Initialize the dataFields array with empty strings.
	for(int i = 0; i < fieldCount; i++)
	{
		strcpy(dataFields[i], "\0");
	}
	
	
	int fieldIndex = 0;
	char* dataCopy = strdup(dataEntry); // Create a copy of dataEntry to use with strtok (as it modifies the string)
	char* token = strtok(dataCopy, delimiter); // Tokenize and analyze each field
	while(token && fieldIndex < fieldCount)
	{
		// Check if the token is a double to determine the type of the field
		const char* type = determine_string_representation_type(token);
		strncat(dataFields[fieldIndex], type, characterCount - strlen(dataFields[fieldIndex]) - 1);
		
		// Get the next token
		fieldIndex++;
		token = strtok(NULL, delimiter);
		
	}
	
	
	free(dataCopy);
	
	
	// Ensure all allocated memory is freed if an early return or exception occurs
	if (token == NULL && fieldIndex < fieldCount)
	{
		for (int i = 0; i < fieldCount; i++)
		{
			free(dataFields[i]);
		}
		free(dataFields);
		return NULL;
	}
	
	//print_string_array(dataFields, fieldCount, "determine_data_entry_types");
	
	return dataFields;
}




/**
 * determine_common_data_entry_types
 *
 * Determines the most common data types (numeric or non-numeric) for each field across all data entries,
 * providing a comprehensive analysis of the entire dataset to establish an expected data type for each field.
 * It achieves this by examining the value found at a specified field, for each one of the data entries,
 * and establishing an expected data type for this field based on the most-common(highest occurences) data type found at said field.
 * (numeric or nonnumeric ==> plottable or unplottable)
 *
 * This function is essentially a more rigorous way of determining the data types of fields, one that takes into account
 * the value(s) found at a specified field across the entire dataset, as opposed to the 'determine_data_entry_types' function,
 * which considers only a single data entry(namely the value found at the entry's specified field).
 * The resulting data typing of fields can be used to perform conditional checks on the
 * data, where they are treated as the expected typing and compared to the data typing
 * found by examining the current data entry in isolation, meaning the expected typing
 * can be effectively used as insurance for cases of disruptive data entries, which are any data entries with missing,
 * blank, and/or unexpected typing of fields, where the examination of a single entry for a field would not be representative
 * of the typing of the field if examined through a broader-scope analysis of the entire dataset(this function).
 *
 * The functionality to perform conditional checks on the data found at a specified field by cross-referencing the expected
 * data typing with the data typing found by independently examining a data entry massively aids in the preprocessing phase of
 * data analysis, which is heavily dependent on both the accuracy of the expected data typings established in this function as
 * well as the appropriate grouping and formatting of types of data to distincly separate the plottable data from the unplottable data.
 *
 *
 * @param dataEntries An array of strings where each string is a data entry.
 * @param entryCount The number of data entries.
 * @param fieldCount The number of fields in each data entry.
 * @param delimiter The character used as a delimiter in the dataset.
 * @return An array of strings where each string represents the most common data type of a field.
 */
char** determine_common_data_entry_types(const char **dataEntries, int entryCount, int fieldCount, const char *delimiter)
{
	/** char*** dataTypesPerField: A 3D array to store the data types of each field in each data entry.
	 * There are three layers of arrays, the innermost, the middle, and the outermost.
	 * The elements of the innermost arrays are the individual data types of every data entry for a single field, then, these innermost arrays
	 * are themselves elements of the middle arrays, where each middle array represents the data typings for all the entries of a single field, and
	 * lastly, the middle arrays are elements of the outermost array, where each element represents the data typings for all the fields of the data set.
	 *
	 * innermost: dataTypesPerField[fieldIndex][entryIndex]		-		Data type of one single field value, extracted from a single entry
	 * middle: dataTypesPerField[fieldIndex]		-		Collection of data types found for a single field at each data entry
	 * outermost: dataTypesPerField			-		Collection of data types found for all fields at each data entry
	 *
	 *
	 * The structure of the dataTypesPerField array is as follows:
	 *
	 * 		Data Set Field Types
	 *				dataTypesPerField[0] 	// Field 0
	 *				|
	 *				|----> dataTypesPerField[0][0] 	// Field 0, Entry 0
	 *					 | dataTypesPerField[0][1] 	// Field 0, Entry 1
	 *					 ...
	 *					 | dataTypesPerField[0][entryCount] 	// Field 0, Entry entryCount
	 *				|
	 *				V
	 *				dataTypesPerField[1] 	// Field 1
	 *				|
	 *				|----> dataTypesPerField[1][0] 	// Field 1, Entry 0
	 *					 | dataTypesPerField[1][1] 	// Field 1, Entry 1
	 *					 ...
	 *					 | dataTypesPerField[1][entryCount] 	// Field 1, Entry entryCount
	 *				|
	 *				V
	 *				...
	 *				dataTypesPerField[fieldIndex]  	// Field fieldIndex
	 *				|
	 *				|----> dataTypesPerField[fieldIndex][0] 	// Field fieldIndex, Entry 0
	 *					 | dataTypesPerField[fieldIndex][1] 	// Field fieldIndex, Entry 1
	 *					 ...
	 *					 | dataTypesPerField[fieldIndex][entryCount] 	// Field fieldIndex, Entry entryCount
	 */
	char*** dataTypesPerField = (char***)malloc(fieldCount * sizeof(char**)); // Allocate memory for storing data types of each field.
	
	
	// Initialize the dataTypesPerField array with empty strings.
	for (int i = 0; i < fieldCount; i++)
	{
		dataTypesPerField[i] = NULL;  // Initialize to NULL to safely call realloc later
	}
	
	
	
	// Loop through each data entry.
	for(int i = 0; i < entryCount; i++)
	{
		char* dataCopy = strdup(dataEntries[i]); // Create a copy to use with strtok (as it modifies the string)
		char* token = strtok(dataCopy, delimiter);
		int fieldIndex = 0;
		
		// Tokenize and analyze each field in the current data entry.
		while(token && fieldIndex < fieldCount)
		{
			const char* type = determine_string_representation_type(token); // Determine the data type for the current token/field.
			dataTypesPerField[fieldIndex] = (char**)realloc(dataTypesPerField[fieldIndex], (i + 1) * sizeof(char*)); // Reallocate memory to store the new type.
			if (dataTypesPerField[fieldIndex] == NULL)
			{
				perror("\n\nError: Unable to reallocate memory in 'determine_common_data_entry_types'.\n");
				exit(1);
			}
			
			
			dataTypesPerField[fieldIndex][i] = strdup(type); // Store the determined type.
			if (dataTypesPerField[fieldIndex][i] == NULL)
			{
				perror("\n\nError: Unable to duplicate string in 'determine_common_data_entry_types'.\n");
				exit(1);
			}
			
			
			fieldIndex++;
			token = strtok(NULL, delimiter);
		}
		
		free(dataCopy);
	}
	
	// Prepare the output array for the most common data types.
	char** commonDataTypes = (char**)malloc(fieldCount * sizeof(char*));
	if (commonDataTypes == NULL)
	{
		perror("\n\nError: Unable to allocate memory for commonDataTypes in 'determine_common_data_entry_types'.\n");
		exit(1);
	}
	
	
	// Determine the most common data type for each field.
	for(int i = 0; i < fieldCount; i++)
	{
		commonDataTypes[i] = determine_most_common_string(dataTypesPerField[i], entryCount); // The
	}
	
	
	// Free the temporary storage.
	for(int i = 0; i < fieldCount; i++)
	{
		for(int j = 0; j < entryCount; j++)
		{
			free(dataTypesPerField[i][j]);
		}
		free(dataTypesPerField[i]);
	}
	free(dataTypesPerField);
	
	
	
	//print_string_array(commonDataTypes, fieldCount, "determine_common_data_entry_types");
	
	
	return commonDataTypes; // Return the array of most common data types.
}








/**
 * format_data_entry_for_plotting
 *
 * Formats a data entry for plotting by ensuring each field is of the correct data type.
 * Primarily used for preparing data lines for plotting where data types need to be consistent.
 * The function processes a line from the dataset, converting fields to a default value if they don't match
 * the expected data type (i.e., replacing non-numeric values with '0.0' in numeric fields).
 *
 * @param headerLine The header line of the dataset.
 * @param dataEntry A line from the dataset.
 * @param fieldCount Number of fields in the data entry.
 * @param delimiter The delimiter used in the dataset.
 * @return A new string with the formatted data entry.
 */
char *format_data_entry_for_plotting(char *headerLine, char* dataEntry, int fieldCount, const char *delimiter)
{
	if (!dataEntry)
	{
		perror("\n\nError: NULL argument provided for dataEntry in 'format_data_entry_for_plotting'.");
		return NULL;
	}
	
	
	
	// Create a buffer for the formatted data entry, making it larger to accommodate potential changes
	size_t bufferSize = strlen(dataEntry) + fieldCount * 4; // Extra space for "0.0" and commas
	char* formattedDataEntry = allocate_memory_char_ptr(bufferSize);
	
	strcpy(&formattedDataEntry[0], "\0"); // Initialize the string empty
	char* dataCopy = strdup(dataEntry);
	char* token;
	char* saveptr = dataCopy;
	token = strtok_r(dataCopy, delimiter, &saveptr);
	int formatIndex = 0;
	
	
	// Determine the data types for each field in the header.
	char** dataTypes = determine_data_entry_types(headerLine, fieldCount, delimiter);
	if (!dataTypes)
	{
		perror("\n\nError: NULL dataTypes in 'format_data_entry_for_plotting'.");
		free(formattedDataEntry);
		return NULL;
	}
	//print_string_array(dataTypes, fieldCount, "format_data_entry_for_plotting dataTypes");
	
	// Process each field based on its data type.
	while(token && formatIndex < fieldCount)
	{
		// Trim leading spaces from the token
		while (isspace((unsigned char)*token)) token++;
		
		// Determine the type of data in the token.
		const char* typeDataEntry = determine_string_representation_type(token);
		
		
		/// Handling incorrectly formatted fields.
		if (strcmp(token, "-") == 0) // Incorrectly formatted field, effectively a missing/blank value.
		{
			// If the field/token is comprised of only a single character(between delimiters within a data entry) which is found to be a hyphen, "-", then it is treated as a blank or missing space
			// As a result of the missing value for the field being represented as a nonnumeric character, the hyphen, then only the expected data type is relevant.
			if (strcmp(dataTypes[formatIndex], "numeric") == 0) // If the missing field was expected to be a numeric value, append a default value in place of the missing field.
			{
				strcat(formattedDataEntry, "0.0"); // Replace with default numeric value.
			}
			else if(strcmp(dataTypes[formatIndex], "numeric") != 0 ) // This "missing" field should be a nonnumeric value, append the hyphen as is or append some sort of arbitrary substitute(such as 0.0), so long as the functions for extracting numeric plottable data are not affected.
			{
				strcat(formattedDataEntry, "0.0"); // Replace with default numeric value... temporary, ideally the entire data set should adjust itself to ensure data is consistent.
			}
		}
		else if (strcmp(typeDataEntry, dataTypes[formatIndex]) != 0) //The determined data type does not match the expected data type
		{
			if (strcmp(dataTypes[formatIndex], "numeric") == 0) // Will be true when the field should have a numeric value but instead has a non-numeric value
			{
				strcat(formattedDataEntry, "0.0"); // Replace non-numeric values in numeric fields.
			}
			// If it is expected to be a string but is not, we remove it (by doing nothing here)
		}
		else // The determined data types do match and the type is numeric
		{
			strcat(formattedDataEntry, token); // Append correct token.
		}
		
		
		
		// Append delimiter except for the last token
		if (formatIndex < fieldCount - 1)
		{
			strncat(formattedDataEntry, delimiter, bufferSize - strlen(formattedDataEntry) - 1);
		}
		
		// Increment the index for the current field on this dataEntry and move to the next token(the next delimiter separated field in the dataEntry)
		formatIndex++;
		token = strtok_r(saveptr, delimiter, &saveptr);
	}
	
	
	
	free(dataCopy);
	for (int i = 0; i < fieldCount; i++)
	{
		free(dataTypes[i]);
	}
	free(dataTypes);
	
	
	
	// Reallocate formattedDataEntry to fit the actual length, if realloc fails, return the original formattedDataEntry without adjusting its size.
	char *fitOutput = realloc(formattedDataEntry, strlen(formattedDataEntry) + 1);
	if (fitOutput)
	{
		formattedDataEntry = fitOutput;
	}
	
	return formattedDataEntry;
}



/**
 * identify_plottable_fields
 *
 * This function identifies which fields from a given array of fieldNameTypePairs are plottable based on their data type.
 * It returns an array of integers where each index corresponds to a field, and the value at that index is 1 if the field is plottable (numeric),
 * and 0 otherwise. The function uses a type delimiter to parse the type from each fieldNameTypePair string.
 *
 * @param fieldNameTypePairs An array of strings where each string represents a field name paired with its data type.
 * @param fieldCount Number of fields in the data entry.
 * @param typeDelimiter The character used as a delimiter to separate the field name and type.
 * @return An array of integers where each index corresponds to a field, and the value at that index is 1 if the field is plottable (numeric), and 0 otherwise.
 */
int* identify_plottable_fields(char** fieldNameTypePairs, int fieldCount, const char *typeDelimiter)
{
	// Allocate memory for the plottability status array and handle potential memory allocation failure.
	int* plottabilityStatus = (int*)malloc(fieldCount * sizeof(int));
	if (plottabilityStatus == NULL)
	{
		perror("\n\nError: Memory allocation failed for plottabilityStatus in 'identify_plottable_fields'.\n");
		return NULL;
	}
	
	/// Iterate through each field to determine if it is plottable based on its data type.
	for (int i = 0; i < fieldCount; i++)
	{
		// Locate the delimiter in the string to separate the field name from its type.
		char* typeDelimiterPosition = strchr(fieldNameTypePairs[i], typeDelimiter[0]);
		if (typeDelimiterPosition != NULL)
		{
			// Extract the type part of the string following the delimiter.
			char* type = typeDelimiterPosition + 1;
			// Check if the type is 'numeric' or matches 'numeric' when processed.
			if (strcmp(type, "numeric") == 0 || strcmp(determine_string_representation_type(type), "numeric") == 0)
			{
				plottabilityStatus[i] = 1; // Mark the field as plottable.
			}
			else
			{
				plottabilityStatus[i] = 0; // Mark the field as not plottable.
			}
		}
		else
		{
			plottabilityStatus[i] = 0; // Default to not plottable if parsing fails.
		}
	}
	return plottabilityStatus;
}





/**
 * capture_plottable_fields
 *
 * This function processes a line from the dataset to extract and format fields that are numeric, as these are typically used for plotting.
 * It first formats the entire data entry for plotting, then iterates through each field, checking if it is numeric.
 * Numeric fields are stored in a dynamically allocated array, which is returned.
 * If no numeric fields are found, or if any memory allocation fails, the function cleans up any allocated memory and returns NULL.
 *
 * @param headerLine The header line of the dataset, used for formatting the data entry.
 * @param dataEntry A line from the dataset to be processed.
 * @param fieldCount The total number of fields in the data entry.
 * @param delimiter The delimiter used to separate fields in the data entry.
 * @return An array of strings containing only the numeric fields from the data entry, or NULL if an error occurs.
 */
char** capture_plottable_fields(char *headerLine, char* dataEntry, int fieldCount, const char *delimiter)
{
	// Format the data entry for plotting.
	char *formattedData = format_data_entry_for_plotting(headerLine, dataEntry, fieldCount, delimiter);
	//printf("\n\nformattedData: %s\n\n", formattedData);
	
	
	// Allocate memory for plottable fields.
	char** plottableFields = allocate_memory_char_ptr_ptr(0, fieldCount);
	
	// Duplicate formatted data for processing.
	char* dataCopy = strdup(formattedData);
	
	char* token = strtok(dataCopy, delimiter);
	int plottableIndex = 0;
	
	/// Extract numeric fields for plotting.
	while (token)
	{
		// Check if the field is numeric.
		if (string_is_numeric(token))
		{
			// Allocate memory for the plottable field
			plottableFields[plottableIndex] = allocate_memory_char_ptr(strlen(token) + 1);
			if (!plottableFields[plottableIndex]) // Check if memory allocation failed
			{
				// Free previously allocated strings up to the current index
				for (int j = 0; j < plottableIndex; j++)
				{
					free(plottableFields[j]);
				}
				// Free the array of pointers
				deallocate_memory_char_ptr_ptr(plottableFields, fieldCount);
				
				// Free other resources
				free(dataCopy);
				free(formattedData);
				return NULL;
			}
			
			// Copy the token to the plottable field
			strcpy(plottableFields[plottableIndex], token);
			plottableIndex++;
		}
		// Get the next token
		token = strtok(NULL, delimiter);
	}
	
	
	// Handle case where no plottable fields were found.
	if (plottableIndex == 0)
	{
		deallocate_memory_char_ptr_ptr(plottableFields, fieldCount);
		free(dataCopy);
		free(formattedData);
		return NULL;
	}
	
	
	//
	free(dataCopy);
	free(formattedData);
	return plottableFields;
}



/**
 * capture_data_set_for_plotting
 *
 * Processes a dataset for plotting by extracting and formatting numeric fields.
 * Goes through each line in the dataset, capturing only the numeric fields, and formats them for plotting.
 *
 * @param fileContents Array of strings representing the dataset.
 * @param lineCount The number of lines in the dataset.
 * @param delimiter The delimiter used in the dataset.
 * @return An array of strings, each representing a line in the dataset with only plottable fields.
 */
char** capture_data_set_for_plotting(char** fileContents, int lineCount, const char *delimiter)
{
	// Assume the first line is the header to determine the format
	char *headerLine = fileContents[0]; // The header line is the first line of the file that is assumed to contain the names of each data field seperated by the delimiter.
	char *formatLine = fileContents[1]; // The format line is the line with the first data entries and whose contents will be analyzed to serve as the basis for formatting.
	int fieldCount = count_data_fields(formatLine);
	
	// Allocate memory for the new dataset
	char** processedDataSet = allocate_memory_char_ptr_ptr(MAX_STRING_SIZE, lineCount + 1); // Assuming max length of processed line
	
	// Extracting solely the numeric values from all of the plottable fields
	for (int i = 0; i < lineCount; i++)
	{
		// Extract plottable fields from each line of the dataset.
		char** plottableFields = capture_plottable_fields(formatLine, fileContents[i], fieldCount, delimiter);
		
		// Concatenate plottable fields into a single string for each line.
		char* concatenatedFields = concatenate_string_array(plottableFields, fieldCount, delimiter);
		
		// Store the concatenated string in the processed dataset.
		processedDataSet[i] = concatenatedFields; // Shift index by 1 to leave space for the header
												  // print the i element of the processed data set
												  ////printf("\n\n processedDataSet[%d]: %s\n\n", i, processedDataSet[i]);
		
		
		// Free the array of plottable fields
		for (int j = 0; j < fieldCount && plottableFields[j] != NULL; j++)
		{
			free(plottableFields[j]);
		}
		free(plottableFields);
	}
	
	
	/// Capture the
	char **fileHeader = capture_data_set_header_for_plotting(headerLine, fileContents, delimiter);
	
	
	/// Identify plottable fields and store the plottability status in an integer array as 1 for plottable and 0 for unplottable
	int *plottabilityStatus = identify_plottable_fields(fileHeader, fieldCount, ":");
	
	
	
	
	/// Loop through all fields and count the number of plottable fields
	int plottableFieldCount = 0;
	for (int i = 0; i < fieldCount; i++)
	{
		if(plottabilityStatus[i] == 1)
		{
			plottableFieldCount++;
		}
	}
	//printf("\n\n plottableFieldCount: %d\n\n\n\n", plottableFieldCount);
	
	
	int unplottableFieldCount = 0; // Start as 0 then as unplottable fields are encountered, increment by 1. Use this variable to offset the index of the fieldNames array when assigning the field names
	
	
	
	// Capture the Data Field Names in a Single String Without the typeDelimiters
	const char *fieldNames[fieldCount]; //The string of field names as they appear in the first element of the fileContents array of strings
	for (int i = 0; i < fieldCount; i++)
	{
		// Identify the delimiting character seperating each field name from it's data type, for example,  'mass:numeric'  --->  typeDelimiter = ':'  &&  typeDelimiterToken = 'mass'
		const char *typeDelimiter = identify_delimiter(fileHeader, fieldCount);
		char* typeDelimiterToken = strtok(fileHeader[i], typeDelimiter); // Split the token by the type delimiter to get the name
																		 //printf("\n\n\n field: %d \n typeDelimiter: %s \n typeDelimiterToken: %s", i, typeDelimiter, typeDelimiterToken);
		
		
		if (typeDelimiterToken)
		{
			if(plottabilityStatus[i] == 1) /// If the token field name is not NULL and the field is plottable, duplicate the token to the field names
			{
				// Offset the index of the fieldNames array by the current number of unplottable fields which preceded this iteration of the loop when assigning the field names
				// so that the ordering of the plottable field names will correspond to the ordering of the plottable fields in the processed dataset
				fieldNames[i - unplottableFieldCount] = strdup(typeDelimiterToken);  // Duplicate the token to the field names while omitting unplottable fieldsa
			}
			else if(plottabilityStatus[i] == 0) /// the token field name is not NULL and the field is NOT plottable
			{
				unplottableFieldCount++; // Increment the count of unplottable fields
				
				
				// Remove the invalid field name from the field names array, decreasing the dimension of the array by 1 and pushing the remaining elements to the left
				//for (int j = i - unplottableFieldCount; j < fieldCount; j++)
				//{
				//    fieldNames[j] = fieldNames[j + 1];
				//}
				
			}
			
		}
	}
	
	
	
	// Concatenate field names into a single string and set as the first line of the processed dataset
	char* headerString = concatenate_string_array(fieldNames, plottableFieldCount, delimiter);
	processedDataSet[0] = headerString;
	
	
	//print_string_array(processedDataSet, lineCount, "Processed Data Set in 'capture_data_set_for_plotting'");
	
	return processedDataSet;
}




/**
 * extract_radix_sorted_data
 *
 * Extracts all plottable data fields from the dataset and radix sorts them(maintaining relative positions).
 * It processes the dataset to extract and store each plottable field in a dedicated file for easy access and plotting.
 *
 *
 * @param dataSetContents Array of strings representing the dataset.
 * @param fieldCount Number of fields in each data entry.
 * @param delimiter The delimiter used in the dataset.
 * @return A 2D array of doubles representing the extracted data fields.
 */
double **extract_radix_sorted_data(char** dataSetContents, int lineCount, int fieldCount, const char *delimiter)
{
	// Count lines and allocate memory for the dataset.
	double **radixSortedDataSet = (double**)malloc(sizeof(double*) * lineCount);
	
	
	
	// Extract each plottable field.
	for(int i = 0; i < fieldCount; i++)
	{
		radixSortedDataSet[i] = extract_plottable_data_field(dataSetContents, i, fieldCount, delimiter);
		
		double* radixSortedField = radixSortedDataSet[i];
		//printf("\n\nField %d", i);
		//print_array(fieldCount, radixSortedField, "Plottable Data Set Unsorted");
		radix_sort_doubles(radixSortedField, fieldCount);
		//print_array(fieldCount, radixSortedField, "Plottable Data Set Radix Sorted");
		
		
		radixSortedDataSet[i] = radixSortedField;
		free(radixSortedField);
	}
	
	
	return radixSortedDataSet;
}




/**
 * extract_plottable_data_field
 *
 * Extracts a specific plottable data field from the dataset for plotting purposes.
 * It allocates memory for storing the extracted numeric data and processes each line to capture the specified field.
 *
 * @param dataSetContents Array of strings representing the dataset.
 * @param fieldIndex Index of the field to be extracted.
 * @param fieldCount Number of fields in each data entry.
 * @param delimiter The delimiter used in the dataset.
 * @return An array of doubles representing the extracted data field.
 */
double *extract_plottable_data_field(char** dataSetContents, int fieldIndex, int fieldCount, const char *delimiter)
{
	// Count the number of lines in the dataset and allocate memory for plot data.
	int lineCount = count_array_strings(dataSetContents);
	double *plotData = allocate_memory_double_ptr(lineCount);
	
	// Process each line to extract the specified field.
	for (int i = 0; i < lineCount; i++)
	{
		char* dataEntry = dataSetContents[i]; // The current line of the data set.
		char* dataCopy = strdup(dataEntry);
		char* saveptr;
		char* token = strtok_r(dataCopy, delimiter, &saveptr);
		int currentIndex = 0;
		double fieldValue = 0.0;
		
		
		// Extract the field value at the specified index.
		while (token && currentIndex < fieldCount)
		{
			if (currentIndex == fieldIndex) // The field index is passed in to capture the specified data set field only.
			{
				fieldValue = atof(token); // Convert the nonnumeric token to a numeric representation as a double.
				plotData[i] = fieldValue;
				break;
			}
			
			// Not yet at the fieldIndex, iterate the currentIndex and get the next field with token.
			currentIndex++;
			token = strtok_r(NULL, delimiter, &saveptr);
		}
		
		
		// If after iterating through all field indexes, the specified fieldIndex of the data entry has not been found, set it to NULL 0.0.
		if (!token)
		{
			fieldValue = 0.0;
			plotData[i] = fieldValue;
		}
		free(dataCopy);
		//free(dataEntry);
	}
	
	
	
	
	return plotData;
}




/**
 * extract_plottable_data
 *
 * Extracts all plottable data fields from the dataset and writes them into separate files.
 * It processes the dataset to extract and store each plottable field in a dedicated file for easy access and plotting.
 *
 * @param dataSetContents Array of strings representing the dataset.
 * @param fieldCount Number of fields in each data entry.
 * @param delimiter The delimiter used in the dataset.
 * @return A 2D array of doubles representing the extracted data fields.
 */
double **extract_plottable_data(char** dataSetContents, int lineCount, int fieldCount, const char *delimiter)
{
	// Count lines and allocate memory for the dataset.
	double **plottableDataSet = (double**)malloc(sizeof(double*) * lineCount);
	
	
	
	// Extract each plottable field.
	for(int i = 0; i < fieldCount; i++)
	{
		plottableDataSet[i] = extract_plottable_data_field(dataSetContents, i, fieldCount, delimiter);
		
		double* plottableField = plottableDataSet[i];
		//printf("\n\nField %d", i);
		//print_array(fieldCount, plottableField, "Plottable Data Set Unsorted");
		radix_sort_doubles(plottableField, fieldCount);
		//print_array(fieldCount, plottableField, "Plottable Data Set Radix Sorted");
		//free(plottableField);
	}
	
	
	return plottableDataSet;
}








/**
 * write_plottable_data
 *
 * Writes the plottable data extracted from the dataset to files.
 * This function captures the dataset's plottable data and writes it into files for each field,
 * making it easier for later use in plotting or analysis.
 *
 * @param dataSetContents Array of strings representing the dataset.
 * @param headerLine The header line of the dataset.
 * @param dataDirectory Directory where the extracted data will be stored.
 * @param dataSetFileName Name of the dataset file.
 * @param delimiter The delimiter used in the dataset.
 * @return The path name of the files where the plottable data is written.
 */
const char *write_plottable_data(char** dataSetContents, char *headerLine, const char *dataDirectory, char *dataSetFileName, int lineCount, const char *delimiter)
{
	// Capture and process the dataset to extract plottable data
	int fieldCount = count_data_fields(dataSetContents[0]);
	//printf("\n\n\n\n-----fieldCount-----: %d\n\n\n", fieldCount);
	
	double **plottableDataSet = extract_plottable_data(dataSetContents, lineCount, fieldCount, delimiter);
	// perform statistical analysis on this double array array of plottable(numeric) data
	// capture the results of the analysis in a file, where the name is generated based on the test performed and the numbers are stored inside
	
	
	
	
	
	char **dataSetFieldNames = split_tokenized_string(dataSetContents[0], delimiter, fieldCount);
	char **fileHeader = capture_data_set_header_for_plotting(dataSetContents[0], dataSetContents, delimiter);
	
	//printf("\n\n\nfieldCount in 'capture_data_set_for_plotting': %d", fieldCount);
	//printf("\nheaderLine in 'capture_data_set_for_plotting': %s", headerLine);
	//printf("\ndataSetContents[0] in 'capture_data_set_for_plotting': %s\n", dataSetContents[0]);
	//print_string_array(dataSetFieldNames, fieldCount, "Data set field names in 'write_plottable_data'");
	//print_string_array(fileHeader, fieldCount, "File header in 'write_plottable_data'");
	
	
	/// Identify plottable fields and store the plottability status in an integer array as 1 for plottable and 0 for unplottable
	int *plottabilityStatus = identify_plottable_fields(dataSetFieldNames, fieldCount, ":");
	//printf("\n\n\n\n-----plottabilityStatus-----: \n");
	for (int i = 0; i < fieldCount; i++)
	{
		//printf("\nplottabilityStatus[%d]: %d", i, plottabilityStatus[i]);
	}// printf("\n\n\n-------------------------\n\n\n");
	
	
	/// Loop through all fields and count the number of plottable fields
	int plottableFieldCount = 0;
	for (int i = 0; i < fieldCount; i++)
	{
		if(plottabilityStatus[i] == 1)
		{
			plottableFieldCount++;
		}
	}
	
	
	
	
	
	// Write data fields into separate files for each plottable field.
	const char *extractedDataDirectory = dataDirectory;
	
	
	// Path to the directory in which the plottable data fields will be located(to be followed by the actual name of the file, index of the field, and the .txt extension)
	char *plottableDataFieldsDirectoryFilePath = combine_strings("/", combine_strings(dataSetFileName, "_Plottable_Field"));
	const char *plottableFieldsPathName = combine_strings(3, plottableDataFieldsDirectoryFilePath); // Full path for plottable data fields.
	
	
	for(int i = 0; i < fieldCount; i++)
	{
		// Construct file path name for each field.
		char *plottingDataFilePathName = plottableFieldsPathName;
		char *fieldExtracted = allocate_memory_char_ptr(plottableFieldCount);
		sprintf(fieldExtracted, "_%d-%s.txt", i, dataSetFieldNames[i]);
		plottingDataFilePathName = combine_strings(plottingDataFilePathName, fieldExtracted);
		
		
		// Write the data to files.
		FILE *file = fopen(plottingDataFilePathName, "w+");
		write_file_numeric_data(plottingDataFilePathName, plottableDataSet[i], lineCount, dataSetFieldNames[i]);
		fclose(file);
		free(plottingDataFilePathName);
		free(fieldExtracted);
	}
	
	
	
	
	/// Cleanup memory
	free(plottableDataFieldsDirectoryFilePath);
	//free(plottableFieldsPathName);
	deallocate_memory_char_ptr_ptr(dataSetFieldNames, fieldCount);
	deallocate_memory_char_ptr_ptr(fileHeader, fieldCount);
	deallocate_memory_double_ptr_ptr(plottableDataSet, fieldCount);
	
	
	
	
	/// OMITTED, an additional file with all the plottable fields combined into a single file, the individual fields are already stored separately
	const char *plottingDataFilePathName = combine_strings(extractedDataDirectory, ".txt");
	FILE *file = fopen(plottingDataFilePathName, "w+");
	for(int i = 0; i < plottableFieldCount; i++)
	{ write_file_numeric_data(plottingDataFilePathName, plottableDataSet[i], lineCount, dataSetFieldNames[i]); }
	fclose(file);
	
	
	
	
	return extractedDataDirectory;
}




/**
 * write_data_set
 *
 * Processes and writes a dataset to files, categorizing data into plottable and non-plottable types.
 * The function structures the dataset, creating separate files for plottable fields to aid in analysis
 * and plotting.
 *
 * @param fileContents Array of strings representing the dataset.
 * @param filePathName Path of the original dataset file.
 * @param delimiter The delimiter used in the dataset.
 * @return The directory where the processed dataset files are stored.
 */
const char *write_data_set(char** fileContents, const char *filePathName, int lineCount, int fieldCount, const char *delimiter)
{
	// Process the dataset for plotting
	/*-----------   Write the NON-Plottable Data to Files   -----------*/
	char **plottingData = capture_data_set_for_plotting(fileContents, lineCount, delimiter);  // Capture data suitable for plotting.
	
	
	
	
	/*-----------   Create the Name of a File for Plotting Data Set in Same Directory as Original File   -----------*/
	char *directoryPathName = find_file_directory_path(filePathName); // Get the directory path of the file.
	char *fileName = find_name_from_path(filePathName); // Extract the filename from the path.
	
	
	// Create a directory for plottable data fields
	char *dataDirectory = create_directory(filePathName, "_Plottable_Fields"); // Create a directory for this CSV file's plottable data fields.
	
	
	// The directory to be created at the same level and location as the passed in file pathname
	directoryPathName = combine_strings(directoryPathName, fileName);
	directoryPathName = combine_strings(directoryPathName, "_Plottable_Fields");
	
	
	
	
	
	// Write plottable fields to files, Populate the Contents of the Plotting File with the Contents of the Array of Strings(i.e., the data entries)
	// Write plottable fields to the directory at 'directoryPathName' with pathnames 'plottableFieldsPathName'(to be followed by the index of the field and the .txt extension)
	write_plottable_data(plottingData, fileContents[0], directoryPathName, fileName, lineCount, delimiter);
	
	
	
	
	/// Cleanup memory
	free(directoryPathName);
	free(fileName);
	deallocate_memory_char_ptr_ptr(plottingData, lineCount);
	return dataDirectory;
}








