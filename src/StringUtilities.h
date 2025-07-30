//  StringUtilities.h
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02
/**
 * StringUtilities code: provides a comprehensive suite of utilities for string manipulation and analysis, primarily designed for data
 * parsing and preprocessing in the context of data integration and analysis. This collection of functions is particularly useful
 * in contexts where string data must be examined, transformed, or formatted before being used in computations or visualizations.
 *
 * NOTE: Focuses on C-style strings, namely: char, char*, char*[], and char**
 **/


#ifndef StringUtilities_h
#define StringUtilities_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>




// ------------- Helper Functions for Determining Properties of Characters -------------
/// \{
bool char_is_alpha(char c); // Checks if a character is an alphabetic character.
bool char_is_digit(char c); // Checks if a character is a digit.
bool char_is_alnum(char c); // Checks if a character is alphanumeric.
bool char_is_whitespace(char c); // Checks if a character is whitespace.
bool char_is_punctuation(char c); // Checks if a character is punctuation.
bool char_is_underscore(char c); // Checks if a character is an underscore.
bool char_is_sign(char c); // Checks if a character is a sign (+, -, etc.).
bool char_is_delimiter(char c); // Checks if a character is a delimiter (non-alphanumeric and non-space)
								/// \}






// ------------- Helper Functions for Counting and Classifying Strings -------------
/// \{
size_t string_length(const char *characterString); // Returns the length of a character string.
int count_array_strings(char **stringArray); // Counts the number of strings in an array of strings.
int count_character_occurrences(const char *characterString, char c); // Counts occurrences of a specified character in a string
bool string_is_numeric(const char *characterString); // Checks if a string represents a numeric value.
bool string_is_nonnumeric(const char *characterString); // Checks if the given string contains any numeric characters.
bool string_is_numeric_with_alphanumeric(const char* characterString); // Checks if a string is alphanumeric.
bool string_is_hyphen_else_is_minus_sign(char *characterString); // Differentiates between hyphens and minus signs.
void replace_hyphen_occurences_with_underscore(char *characterString); // Detect all occurences of hyphens and replace them with underscores.
int *string_is_date_time(const char *characterString, const char *delimiter, const int fieldCount); // Analyzes a string for date/time formats.
int* string_is_unit(const char *characterString, const char *delimiter, const int fieldCount); // Analyzes a string for units/unitformats.
bool is_numeric_with_units(const char* characterStringToken, char* testUnit); // Checks a string for specified units.
const char *determine_string_representation_type(const char* token); // Determines if a string is numeric or non-numeric.
																	 /// \}






// ------------- Helper Functions for Evaluating and Searching within Strings -------------
/// \{
char find_most_common_non_alphanumeric_character(const char *characterString); // Finds the most common non-alphanumeric character in a string
char *character_last_occurence(const char *characterString, char c); // Returns a pointer to the last occurrence of a character in a string. (~= strrchr)
char *determine_most_common_string(char **stringArray, int stringCount); // Determines the most common string in an array of strings.
char* determine_most_common_delimited_string(char** stringArray, int stringCount, const char* delimiter); // Determines the most common string in an array of delimited strings.
char* current_date_time_string(void); // Returns the current local machine date/time as a string.
bool string_array_contains_date_time(char **stringArray, int stringCount, const char *delimiter); // Checks if an array of strings contains date/time strings.
char *determine_string_occurence(char *queryString, char *searchString); // Determines if a string occurs in another string.
char *find_string_in_string_array(char **stringArray, int stringCount, const char *searchString); // Finds a string in an array of strings.
																								  /// \}






// ------------- Helper Functions for Handling Delimiters in Strings -------------
/// \{
char *find_potential_delimiters(const char *characterString, int *delimiterCount); // Finds potential delimiters in a string.
char *identify_delimiter(char** stringArray, int stringCount); // Identify the most common delimiter in a string array.
char *identify_delimiter_recursive(char** stringArray, int stringCount, char **primaryDelimiters, int *primaryDelimiterCount); // Pass in primary delimiters(and count) by reference... a special case function
																															   /// \}






// ------------- Helper Functions for Comparing, Copying, and Duplicating Strings -------------
/// \{
int compare_strings(const char *characterString1, const char *characterString2); // Compares two character strings for equality. (~= strcmp)
char *duplicate_string(const char *characterString); // Duplicates a character string. (~= strdup)
char *copy_string(char *destination, const char *source); // Copies a character string. (~= strcpy)
char *copy_n_string(char *destination, const char *source, size_t n); // Copies up to n characters from the string pointed to by 'source' and to 'destination', where any n greater than the length of 'source' is padded with null characters. (~= strncpy)
																	  /// \}






// ------------- Helper Functions for Combining and Concatenating Strings -------------
/// \{
char *concatenate_string(char *destination, const char *source); // Concatenates two character strings. (~= strcat)
char *concatenate_n_string(char *destination, const char *source, size_t n); // Concatenates up to n characters from the string pointed to by 'source' to the end of the string pointed to by 'destination'. (~= strncat)
char *combine_strings(const char *characterString1, const char *characterString2);  // Returns a new string with contents of the second string appended to the first.
char **combine_string_arrays(const char **stringArray1, int stringCountArray1, const char **stringArray2, int stringCountArray2); // Combines two arrays of strings.
char *append_string_array_to_string(const char *characterString1, const char *characterStringArray[], int numArrayStrings, const char *delimiter); // Combines each string in the array with the original string to return a single string having the contents of the original string with each string from the array appended to it seperated by the specified delimiter.
char *concatenate_string_array(const char** stringArray, int stringCount, const char *delimiter); // Concatenates an array of strings into a single string with elements/fields seperated by the specified delimiter
																								  /// \}






// ------------- Helper Functions for Tokenizing and Segmenting Strings -------------
/// \{
char *tokenize_string(char *characterString, const char *delimiter); // Tokenizes a character string based on a delimiter. (~= strtok)
char* tokenize_string_r(char* str, const char* delimiters, char** saveptr);
char **split_tokenized_string(const char* characterString, const char* delimiter, int divisions); // Splits a tokenized string into an array of strings based on a given delimiter, i.e., assigns each tokenized field to an element in an array.
																								  /// \}






// ------------- Helper Functions for Pruning and Trimming Strings -------------
/// \{
char *trim_string_whitespaces(char* untrimmedString); // Trims whitespace from a string.
char *prune_string_whitespaces(char *unprunedString); // Removes all whitespaces from a string.
char *prune_repeated_delimiters_from_string(char *unprunedString, const char *delimiter);  // Handles repeated delimiters.
char *prune_and_trim_problematic_characters_from_string(char *originalString, const char *delimiter, const int fieldCount);  // Prunes and trims problematic characters.
																															 /// \}






// ------------- Helper Functions for Preprocessing Strings -------------
/// \{
char** extract_units_from_fields(const char* characterString, const char* delimiter, const int fieldCount); // Extracts units from fields.


char *replace_date_time_with_unix(char* characterString, const char *delimiter, const int fieldCount); // Replaces date/time with Unix time.
char **preprocess_string_array(char **stringArray, int stringCount, const char *delimiter); // Preprocesses an array of strings, trimming and pruning whitespaces, repeated delimiters, and standardizing some variable parameters like date/time strings which have no standard formatting and are standardized by replacing them with unix representation.
																							/// \}









#endif /* StringUtilities_h */
