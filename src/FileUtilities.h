//  FileUtilities.h
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02
/**
 * FileUtilities code: Provides utility functions for handling various file-related operations.
 * It includes helper functions for(but not limited to):
 *
 * - Analyzing and manipulating file paths, including extracting file extensions,
 *   identifying directory paths, and isolating filenames from paths.
 * - Counting data fields and characters within file lines, specifically targeting CSV file formats.
 * - Reading file contents into arrays of strings and writing arrays back to files,
 *   supporting both textual and numeric data.
 * - Creating and deleting directories in a standardized manner, facilitating organized data storage.
 *
 * These utilities are geared towards handling CSV and text file formats, aiding in tasks like
 * data extraction, preprocessing, and organization. They provide a suite of tools for efficiently
 * managing file and directory operations in C programming, particularly useful in data analysis and integration projects.
 */


#ifndef FileUtilities_h
#define FileUtilities_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>




// ------------- Helper Functions for Determining Various File Properties -------------
/// \{
int determine_file_depth(const char* filePathName); // Determines the depth (nested level) of a file in the file system
char* identify_file_extension(const char* filePathName); // Identifies the file extension from a file path
char* find_file_directory_path(const char* filePathName); // Finds the directory path from a full file path
char *find_directory_name_from_path(const char* directoryPathName); // Extracts the directory name from a directory path(taken from filepath using 'find_file_directory_path')
char* find_name_from_path(const char* filePathName); // Extracts the file name from a file path


int count_data_fields(char* headerLine); // Counts the number of data fields in a header line
int count_plot_data_fields(char* lineContents, const char *delimiter); // Counts the number of fields in a plot data line
int count_file_lines(const char* filePathName, int maxLines); // Counts the lines in a file up to a specified maximum
int* count_file_lines_characters(const char* filePathName, int lineCount); // Counts characters in each line of a file
int* count_characters_in_file_lines_range(const char* filePathName, int lineCount, int startLine); // Counts characters in each line in a specified range of lines of a file
int determine_maximum_line_length(const char* filePathName, int lineCount); // Determines the maximum line length in a file
																			/// \}








// ------------- Helper Functions for File I/O Operations -------------
/// \{
char** read_file_contents(const char* filePathName, int lineCount); // Blindly reads the contents of a file into a string array
char** parse_file_contents(const char* filePathName, int lineCount); // Conditionally reads the contents of a file into a string array
void write_file_contents(const char *filename, char **fileContents); // Writes content to a file from a char array
void write_file_numeric_data(const char *filename, double *data, int countDataEntries, const char *dataFieldName); // Writes data to a file from a double array
char* generate_merged_filename(const char* filePath1, const char* filePath2);
char* merge_two_files(const char* filePath1, const char* filePath2);
/// \}






// ------------- Helper Functions for File Creation and Management -------------
/// \{
char *create_directory(const char *filePathName, const char *directoryName); // Creates a new directory based on a file path
int delete_directory(const char *filePathName); // Deletes a directory based on a file path
int count_files_in_directory(const char *directoryPath); // Counts the number of files in a directory
char **get_file_pathnames_in_directory(const char *directoryPath); // Retrieves the pathnames of all files contained within a directory
char **get_file_names_in_directory(const char *directoryPath); // Retrieves the names of all files contained within a directory
															   /// \}






typedef struct
{
	const char *fileName;
	const char *fileExtension;
	const char *filePathName;
	int fileLineCount;
} FileProperties;
FileProperties extract_file_properties(const char *filePathName); // Extracts file directory information.




typedef struct
{
	const char *directoryPathName;
	const char *directoryName;
	
	int fileCount;
	FileProperties *fileProperties;
} DirectoryProperties;
DirectoryProperties extract_directory_properties(const char *directoryPathName); // Extracts directory information.
DirectoryProperties extract_directory_properties_from_file(const char *filePathName); // Extracts directory information.



void print_directory_propertiees(DirectoryProperties fileDirectory);


#endif /* FileUtilities_h */
















// ------------- File Utilities Footnotes for Streamlining Readability of Code Documentation -------------
/// \{
/** footnote 1: Counting String Length
 *     The length of a String is determined by counting all the characters in it (before the null character),
 *     looping over all non NULL characters in string and incrementing the count till the end of the string.
 *
 *     The current implementation being used, 'CountCharacters', is further configured to not include spaces between characters, ' ',
 *     in the character count, the purpose of creating this function is to allow custom configurations to any character counting
 *     needed to be done that cannot be found in the use of the built-in string.h library function strlen() or the sizeof() operator
 */
















/** footnote 2: Reading From File
 *    In the context of this program, the data extracted from files is most often read in and used to populate array(s), allowing
 *    for data to be more easily worked with, although there are still some instances where the data is processed "on the fly" , i.e.,
 *    process data to create a result directly, as opposed to storing it before reading more data.
 *
 *
 *
 *
 * ///   fgets
 *    The 'char *fgets(const char *lineString, int n, FILE *file)' function reads a line from a stream pointed to by
 *    'FILE *file' into an array pointed to by 'char *lineString' from the starting position specified by 'int n'. Additionally,
 *    the stream(file) must have been successfully opened by fopen(and not yet closed), so that reading existing data using the file object's name
 *    and/or reading user inputted data using stdin/
 *
 *    The contents of an individual file line are read into the field 'lineString' as a series of independent characters in
 *    the following format   [ '', '', ..., '\n', '\0' ], where the characters are contained in apostrophes and seperated by commas,
 *    and the reading continues until it is ended by whichever of the following criteria is first met: the number of characters read
 *    is equal to n-1, a newline character '\n' is reached, or the stream is terminated with a null character '\0'.
 *
 *    The contents of the line variable are then stored into the string/char pointer function field. The 'fgets' function will
 *    operate in the way described above, unless either there is an error parsing the file or the end-of-file condition has been met,
 *    both of which will simply return NULL.
 *
 *    Hence, to iterate through a file's lines, fgets is repeatedly called until fgets returns NULL. this is often implemented using a
 *    while loop that continues until fgets retuns NULL.
 *
 *
 *
 *
 * ///   fscanf
 *    The 'int fscanf (FILE *file, const char *formatString, argument-list)' built-in stdio.h library function fscanf reads data
 *    from a file pointed by 'file' in the form of bytes into the locations that are given by the entries in argument-list, if any.
 *    Each argument-list entry(if any) must be a pointer to a variable type that corresponds to type specifier in 'formatString'.
 *    The fscanf function will return the input string if successful, and will return zero if unsuccessful.
 *
 *
 *
 *
 *
 *
 *
 *
 * /// *----- BRIEF : Writing to a File -----*
 * ///   fprintf
 *    The 'int printf(FILE *file, const char *formatString, argument-list);' function is similar to the fgets and fscanf functions in
 *    the fields it accepts as well as how it navigates parsing data in a formatted manner, however, it differs in it's fundamental
 *    purpose, namely, while fgets and fscanf are used to read data from files, fprintf is used to write data into files.
 *    The fprintf() function writes a series of characters and values to an output stream pointed to by 'file', converting each entry in
 *    it's argument-list field, if there are any, and writes to the stream according to the corresponding format specification in the
 *    'formatString'. The fprintf function will returns the number of bytes that are printed or a negative value if an output error occurs.
 */
















/** footnote 3: Counting File Lines
 *    Count the number of lines in the file (i.e., the number of entries in the data set). The current implementation
 *    being used, 'CountFileLines', is further configured to accept a field 'int numLinesToSkip' to skip a number
 *    of specified lines at the beginning of the file, the purpose of this is to only count the lines of the
 *    file that have relevant data, for example to skip header lines or title lines. The skipping of lines at the
 *    beginning of the file is accomplished by iterating over a for loop for the number of lines to be skipped and
 *    calling fgets at each iteration of the loop.
 *
 *
 *
 *
 * /// *----- The 'CountFileLines' Function For the Purpose of Indexing the String Array of File contents -----*
 *    The 'CountFileLines' function has relevance with respect to the 'InitializeFileLinesContents' function
 *    as a consequence of the 'InitializeFileLinesContents' function being required to have indexes of the
 *    string contents of each line that reflect the actual indexes of their location in the file, so that
 *    the line count of the desired data entries can be captured appropriately so as to avoid capturing
 *    something like the header line instead of a data entry at an index of the string array of line contents.
 *
 *    For this, the 'CountFileLines' function has functionality such that the range of lines from the file
 *    that contain relevant data can be determined simply by calling 'CountFileLines' with the number of lines to be
 *    skipped specified as zero, effectively counting the total number of lines in the file, and then comparing this
 *    total number with the passed in 'lineCount'(the original count of relevant data lines) to determine if any
 *    lines need to be skipped, which will result in one of two cases. Case 1, the two values are equal meaning there
 *    were no lines skipped, or case 2, where the two values are not equal, meaning there were a non-zero number of lines
 *    skipped for the data in this file when originally counted. If the result is case 1, no actions need to be taken
 *    as the number of lines to be skipped for the purpose of indexing the string array is initialized as zero, however,
 *    if the result is case 2, the number of lines to be skipped will be set as the absolute value of the difference
 *    between the total number of lines and the number of lines that was counted originally. To assign indexes to the
 *    string array of file contents memory is allocated for a char pointer of the total number of lines, and then
 *    the data entries of the file will be captured as the elements of the string array starting from the
 *    specified number of lines to be skipped(0 if none), and all indexes preceding this starting point of data entries
 *    will be reserved for later use to capture the contents of any header lines for the purpose of data categorization,
 *    formatting, etc.
 */
















/** footnote 4: File Name/Pathname
 *    In plain c there is no way to retreive the name/pathname of a pointer FILE object(excluding OS exclusive link fetching) because
 *    it uses a file handle to read/write, which isn't associated with any named file(and hence pathname) at all.
 *
 *    Hence, in order to repeatedly perform operations on/using the same pointer FILE object without the need
 *    to create a new FILE object each time a 'const char' pointer variable is used to capture the name/pathname of a pointer FILE object,
 *    allowing easy access to the contents of the file while avoding the risk of corrupting or inadvertently modifying
 *    the FILE.
 *    For example:
 *           const char *dataFilePath = "data.txt";
 *           FILE *file = fopen(dataFilePath, "r");
 *
 *    Additionally, the use of the 'dataFilePath' variable also allows for streamlined creation of temporary/proxy pointer FILE objects
 *    based solely on a file pathname, the creation of which would otherwise introduce significant complexity to a code snippet,
 *    requiring the use of handle manipulation and file mapping to implement an operating system specific link fetch.
 */
















/** footnote 5: File Opening
 * ///   fopen
 *    The 'FILE *fopen(const char *filename, const char *mode)' function assigns a pointer FILE object to the file located at the path
 *    specified by 'filename' which is further configured with the 'mode' argument. The built-in stdio.h library function fopen opens
 *    the file whose pathname is the string pointed to by filename then associates a stream with it. When calling fopen the file will
 *    be opened in one of several modes, all of which will attempt to open an existing file, if there is no such file, some will create
 *    one at the path, and some will return NULL, the modes of the fopen function are as follows:
 *
 *           (1.)  r or rb : Read
 *                   Open existing file for reading, fails if the file does not exist or cannot be read.
 *                      r == Open for reading.
 *                      rb == Open for reading in binary mode.
 *
 *
 *            (2.)  w or wb : Truncate and Write
 *                   If the file exists, its contents are overwritten. If the file does not exist, it will be created.
 *                   Reset contents of existing file for writing, if file does not exist, it is created.
 *                      w == Open for writing.
 *                      wb == Open for writing in binary mode.
 *
 *
 *           (3.)  a or ab : Append
 *                   Open file for append, writing data at end-of-file regardless of intervening calls to fseek(),
 *                   if file does not exist, it is created.
 *                      a == Open for append.
 *                      ab == Open for append in binary mode.
 *
 *
 *           (4.)  r+ or rb+ or r+b : Update
 *                   Open file for update, able to execute both input and output, fails if the file does not exist or cannot be read.
 *                      r+ ==  Open for both reading and writing.
 *                      rb+ ==  r+b == Open for both reading and writing in binary mode.
 *
 *
 *           (5.)  w+ or wb+ or w+b : Truncate and Update
 *                   Reset contents of existing file for update, if file does not exist, it is created.
 *                      w+ ==  Open for both reading and writing.
 *                      wb+ ==  w+b == Open for both reading and writing in binary mode.
 *
 *
 *           (6.)  a+ or ab+ or a+b : Append and Update
 *                   Open file for update, able to execute both input and output. Writing at end-of-file regardless of intervening calls to fseek(),
 *                   if file does not exist, it is created.
 *                      a+ == Open for both reading and appending.
 *                      ab+ == a+b == Open for both reading and appending in binary mode.
 *
 *
 *          NOTE: For 'Update' mode, output must not be directly followed by input without an
 *          intervening call to fflush() or to a file positioning function (.Fn fseek , fsetpos() or rewind()), and input must not
 *          be directly followed by output without an intervening call to a file positioning function, unless the input operation
 *          encounters end-of-file.
 *
 *
 * /// *----- Ensuring The File is Opened Properly -----*
 *     After the initial opening of the file using fopen, a check is performed to ensure the file was successfully opened by checking
 *     to see if the file variable == NULL, and then if it is, an error has occured, and the program will be unable to continue.
 *     To address this, a debugging mode function is called to dig into the errno code to determine the reason fopen returned NULL
 *     and then present the user with relevant debug information and provide a list of possible reasons for fopen failing as well
 *     as steps for correcting potential issues.
 */
/// \}
