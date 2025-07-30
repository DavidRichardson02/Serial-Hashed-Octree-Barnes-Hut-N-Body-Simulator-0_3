//  CommonDefinitions.h
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02
/**
 * CommonDefinitions code: The CommonDefinitions file provides a set of common constants, macros, and functions used throughout the program.
 *
 * The majority of their use is found in facilitating special case operations such as detecting occurrences of date/time strings and extracting
 * string representations of variable names from arrays.
 * The creation of this file and the inclusion of these definitions were primarily motivated by the need to distinctly separate these special case
 * definitions from the main program logic and to provide a convenient way to access these definitions across the program.
 *
 * Furthermore, while the majority of these definitions find use exclusively in special case operations, there are a number that are more general in nature,
 * and find use in a variety of contexts throughout the program. Meaning, the definitions included in this file may be broadly categorized into two groups:
 * 1. General definitions: These definitions are widely applicable and are used in a variety of contexts throughout the program.
 * 2. Special case definitions: These definitions are used in specific contexts and are not widely applicable across the program.
 *
 * The general definitions include constants such as 'MAX_STRING_SIZE', which specifies the maximum size of a string, and macros
 * such as 'ARRAY_SIZE', which computes the number of elements in a statically-allocated array. The special case definitions include
 * constants such as 'months', which is an array of the names of the twelve months of the year, and 'commonDateTimeFormats', which is an array
 * of commonly used date/time formats and is used for parsing and formatting date/time strings.
 */


#ifndef CommonDefinitions_h
#define CommonDefinitions_h

#include <stdio.h>
#include <pthread.h> // For POSIX Threads




/// ====================================   (1.) GENERAL DEFINITIONS   ====================================
/// \{
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0])) // Macro to compute the number of elements in a statically-allocated array.


static const size_t MAX_STRING_SIZE = 1000; // Maximum string size.
static const size_t MAX_NUM_FILE_LINES = 100000; // Maximum number of lines in a file.

/// \}












/// ====================================   (2.) SPECIAL CASE DEFINITIONS   ====================================
/// \{


/**@{
 ------------------ MACROS FOR EXTRACTING STRING REPRESENTATIONS FROM VARIABLE NAMES AND ARRAY ELEMENTS ------------------
 * Macros defined for extracting string representations from names of variables and array elements.
 *
 * Throughout this program, there are several instances where both the name, denoting the type of data, and the contents of the actual
 * data are stored together as string representations within the same array(a double value of 6.0 will be stored as "6.0"). In such cases,
 * it becomes necessary to have methods for separating the name of the variable stored at a particular index of the array from the
 * contents of the variable stored at that index.
 */
#define VAR_NAME_AS_STRING(var) #var // Preprocessor macro that creates a string representation of its argument, the name of the variable
#define ARRAY_ELEMENT_NAME_AS_STRING(var) #var[0] // Preprocessor macro that creates a string representation of the first element of the array in the argument
#define BUFFER_SIZE 256 // Buffer size for the string representation of an array element. Used primarily with the ARRAY_ELEMENT_TO_STRING macro to handle string conversions.


// Preprocessor macro that creates a string representation of the n-element's variable name from a string array.
#define ARRAY_ELEMENT_TO_STRING(arr, n) ({ \
static char buffer[BUFFER_SIZE]; \
n_array_element_to_string(arr, n, buffer); \
buffer; \
}) //... used to extract the name of the variable stored at the n-th index of a string array
const char* n_array_element_to_string(char** stringArray, int n, char* buffer); // Macro helper function to convert the n element of an array to a string
/** @}*/








/**@{
 ------------------ CONSTANTS FACILITATING DATE/TIME OPERATIONS INVOLVING STRINGS ------------------
 * These constants are useful for converting numeric representations of months and days into their corresponding textual names,
 * as well as for parsing and standardizing date/time strings from diverse data sources.
 **/


/**
 * This mutex should be locked before calling non-thread-safe functions that deal with localtime
 * and unlocked after their usage.
 * @example
 * pthread_mutex_lock(&localtime_mutex);
 * // Call to thread-unsafe functions
 * pthread_mutex_unlock(&localtime_mutex);
 */
static pthread_mutex_t localtime_mutex = PTHREAD_MUTEX_INITIALIZER; //Mutex for Thread Safety protecting access to localtime


/* 'months' array:
 This array contains the names of the twelve months of the year.
 It is useful for converting numeric month representations into their corresponding textual names. */
extern const char *months[12];


/* 'weekDays' array:
 This array holds the names of the seven days of the week.
 It can be used to map numeric day representations to their corresponding day names. */
extern const char *weekDays[7];


/* 'commonDateTimeFormats' array:
 An array of string literals representing commonly used date/time formats.
 These formats are crucial in parsing and standardizing date/time strings from diverse data sources,
 and follow the conventions used in 'strftime' and 'strptime' functions for date-time formatting and parsing.
 Each format specifier adheres to the C standard library conventions for date/time formatting. */
extern const char *commonDateTimeFormats[12]; // Used for converting date/time strings from datasets into standardized formats.
/** @}*/








/**@{
 ------------------ COMMON UNIT FORMATS ------------------
 * These constants represent commonly used unit formats for various measurements.
 * Cover a variety of scientific and engineering measurements, including:
 *		1. **Basic Units**: These are the standard units from the International System of Units (SI).
 *		2. **Derived Units with SI Prefixes**: These include units that are scaled by powers of ten to represent larger or smaller quantities.
 *		3. **Compound Units**: These involve combinations of basic units to represent rates, densities, and other derived quantities.
 *		4. **Specialized Units**: Specific units used in fields such as electrical engineering, physics, and chemistry.
 *		5. **Units with per time variations**: These represent rates commonly used in everyday contexts as well as scientific calculations.
 */
extern const char *commonUnitFormats[45]; // An array of string literals representing commonly used unit formats for representing physical quantities in data sets.
/** @}*/








/**@{
 ------------------ CONSTANTS FACILITATING NEW FILE NAMING ------------------
 * These constants are used for defining standard formats and identifiers for new text files.
 * They are used in the creation of new files and in the naming of files generated by the program.
 **/
/** @}*/


/// \}







#endif /* CommonDefinitions_h */
