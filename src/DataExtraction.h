//  DataExtraction.h
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02
/**
 * DataExtraction
 * DataExtraction code: The DataExtraction file provides a set of functions for processing, analyzing, and formatting data sets
 * stored in CSV format, primarily for plotting purposes. The functions are designed to facilitate data visualization and plotting by
 * preparing datasets in a format suitable for study via statistical methods and later for graphical representation. It emphasizes data
 * preprocessing to ensure compatibility with plotting tools and data analysis applications.
 *
 * It includes functionalities for parsing CSV files, counting missing values, determining data entry types, and extracting plottable fields.
 * The file also handles formatting data entries for plotting, ensuring correct data types and handling missing values. Additionally, it offers
 * functions for creating and populating formatted files from data sets, separating plottable and non-plottable data for further analysis.
 *
 * NOTE:
 * The functions in this file are primarily designed for datasets where the first line is the header, and the following lines contain data entries.
 * The file assumes CSV format, but the delimiter can be specified to accommodate different formats.
 *
 *
 * Key Features:
 *       - Parsing CSV files into structured arrays for easy manipulation.
 *       - Counting missing values to assess data quality and completeness.
 *       - Determining data types for fields to ensure correct processing and plotting.
 *       - Identifying and extracting plottable (numeric) fields for visualization.
 *       - Formatting data entries and entire datasets for consistent and accurate plotting.
 *       - Writing processed and plottable data to separate files for further analysis and visualization.
 */


#ifndef DataExtraction_h
#define DataExtraction_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>




// ------------- Helper Functions for Capturing and Formatting Data from an Existing File/Data Set -------------
/// \{
char*** parse_entire_file(char **fileContents, int lineCount, int* fieldCount, const char *delimiter);
int* count_missing_values(char** fileContents, int lineCount, int fieldCount, const char *delimiter, const char *headerLine);
char** capture_data_set_header_for_plotting(char* headerContents, char** fileContents, const char *delimiter); // Captures and formats header information for plotting purposes.
char** determine_data_entry_types(const char *dataEntry, int fieldCount, const char *delimiter); // Determines the data types of fields in a data entry.
char** determine_common_data_entry_types(const char **dataEntries, int entryCount, int fieldCount, const char *delimiter); // Determines the data types of entries in a dataset.


int* identify_plottable_fields(char** fieldNameTypePairs, int fieldCount, const char *typeDelimiter); // Identifies the plottable fields in the dataset and captures their indexes with a value of 0 if the field is not plottable, and a value of 1 if the field is plottable.
char** capture_plottable_fields(char *headerLine, char* dataEntry, int fieldCount, const char *delimiter); // Captures fields from a data entry suitable for plotting.
char *format_data_entry_for_plotting(char *headerLine, char* dataEntry, int fieldCount, const char *delimiter); // Formats a data entry for plotting, ensuring correct data type for each field.
char** capture_data_set_for_plotting(char** fileContents, int lineCount, const char *delimiter); // Processes a dataset for plotting by extracting and formatting numeric fields.
																								 /// \}






// ------------- Helper Functions for Extracting Plottable Data Fields -------------
/// \{
double **extract_radix_sorted_data(char** dataSetContents, int lineCount, int fieldCount, const char *delimiter); // Extracts all plottable data fields from the dataset and radix sorts them(maintaining relative positions).
double *extract_plottable_data_field(char** dataSetContents, int fieldIndex, int fieldCount, const char *delimiter); // Writes the plottable data extracted from the dataset to files.
double **extract_plottable_data(char** dataSetContents, int lineCount, int fieldCount, const char *delimiter); // Extracts all plottable data fields from the dataset and writes them into separate files.
																											   /// \}






// ------------- Helper Functions for Creating and Populating a Formatted File from a Data Set -------------
/// \{
const char *write_plottable_data(char** dataSetContents, char *headerLine, const char *dataDirectory, char *dataSetFileName, int lineCount, const char *delimiter); // Writes the plottable data extracted from the dataset to files, returns the pathname to the directory of contained files.
const char *write_data_set(char** fileContents, const char *filePathName, int lineCount, int fieldCount, const char *delimiter); // Processes and writes a dataset to files, separating plottable and non-plottable data.
																																 /// \}














#endif /* DataExtraction_h */
