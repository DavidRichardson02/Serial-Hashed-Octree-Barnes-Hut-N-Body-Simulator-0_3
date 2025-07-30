//  GeneralUtilities.h
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02
/**
 * GeneralUtilities code: Provides a collection of general utility functions that support various computational
 * and data manipulation tasks across various contexts. Encompasses functions for safe memory allocation for different data types, including basic and derived types,
 * and functions for handling date/time conversions, particularly focusing on Unix time conversions.
 *
 * Additionally, it includes utilities for bitwise operations on numerical data, useful in scenarios requiring
 * direct manipulation of the binary representation of data. It also provides robust sorting algorithms, including
 * merge sort and radix sort, optimized for handling large datasets, especially useful for sorting double precision
 * floating-point numbers.
 */


#ifndef GeneralUtilities_h
#define GeneralUtilities_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>




// ------------- Helper Functions for Allocating Memory Safely For Basic Types -------------
/// \{
int *allocate_memory_int_ptr(size_t sizeI);
float *allocate_memory_float_ptr(size_t sizeF);
double *allocate_memory_double_ptr(size_t sizeD);
char *allocate_memory_char_ptr(size_t sizeC);
/// \}




// ------------- Helper Functions for Allocating and Deallocating Memory Safely For Derived Types -------------
/// \{
int **allocate_memory_int_ptr_ptr(size_t sizeI);
float **allocate_memory_float_ptr_ptr(size_t sizeF);
double **allocate_memory_double_ptr_ptr(size_t sizeD);
char **allocate_memory_char_ptr_ptr(size_t strSize, size_t numStrings);

void deallocate_memory_int_ptr_ptr(int **intPtrPtr, size_t numInts);
void deallocate_memory_float_ptr_ptr(float **floatPtrPtr, size_t numFloats);
void deallocate_memory_double_ptr_ptr(double **doublePtrPtr, size_t numDoubles);
void deallocate_memory_char_ptr_ptr(char **charPtrPtr, size_t numStrings);
/// \}






// ------------- Helper Functions for Performing Various Mathematical Operations on Containers -------------
/// \{
double minimum(double a, double b); // Returns the minimum of two values.
double maximum(double a, double b);  // Returns the maximum of two values.
double min_element(double *data, int n); // Returns the minimum element in an array of values.
double max_element(double *data, int n); // Returns the maximum element in an array of values.
int max_element_int(int *data, int n); // Returns the maximum element in an array of values.
double sum_elements(double *data, int n); // Returns the sum of all elements in an array.
/// \}






// ------------- Helper Functions for Operations with Time -------------
/// \{
time_t convert_to_unix_time(const char *dateTimeString); // Converts a date/time string into Unix time.
struct tm *thread_safe_localtime(const time_t *tim, struct tm *result); // A thread-safe wrapper around localtime.
																		/// \}






// ------------- Helper Functions for Bitwise Operations(on Binary Numerals?) -------------
/// \{
uint64_t flip_sign_bit(uint64_t value); // Helper function to flip the sign bit of the double's binary representation.
uint64_t double_to_uint64(double value); // Helper function to reinterpret a double as an uint64_t.
double uint64_to_double(uint64_t value); // Helper function to reinterpret a uint64_t as an double.
uint64_t double_to_mapped_uint64(double value); // Helper function to map a double to an uint64_t.
double mapped_uint64_to_double(uint64_t u); // Helper function to map a uint64_t to a double.
											/// \}









// ------------- Helper Functions for Sorting -------------
/// \{
void merge_data(double *unsortedData, int left, int middle, int right); // Merges two sorted subarrays into a single sorted array.
void merge_sort_data(double *unsortedData, int left, int right); // Recursively divides and sorts an array using the merge sort algorithm.
void merge_sort(double *unsortedData, const int numElements); // Sorts an array of doubles using the merge sort algorithm.


void radix_sort_doubles(double *unsortedData, const int numElements); // Sorts an array of double precision floating-point numbers by interpreting their bit patterns as unsigned 64-bit integers using Radix Sort.
																	  /// \}






// ------------- Helper Functions for Performing Memory Operations -------------
/// \{
void *set_memory_block(void *block, int c, size_t n); // Sets the first 'n' bytes of the memory 'block' to the value specified by 'c'. ( ~= memset)
void *copy_memory_block(void *destination, const void *source, size_t n); // Copies the specified number of bytes from the source memory block to the destination memory block. ( ~= memcpy)
																		  /// \}






#endif /* GeneralUtilities_h */
















// ------------- General Utilities Footnotes for Streamlining Readability of Code Documentation -------------
/** footnote 1: Memory Allocation for Basic and Derived Types
 *
 **/








/** footnote 2: Operations with Time
 *    These functions assist with converting date and time representations, particularly focusing on converting date/time strings
 *    into Unix time and providing thread-safe access to local time structures.
 *
 *    convert_to_unix_time:
 *    		This function attempts to parse a date/time string using various common date/time formats specified in an array (e.g., ISO 8601 format, RFC 2822 format).
 *    		It uses the strptime function to parse the string into a struct tm, and then converts it to Unix time (seconds since the Unix Epoch, January 1, 1970) using mktime.
 *    		If the parsing is successful, it returns the Unix time. If all parsing attempts fail, it returns -1.
 *
 *    thread_safe_localtime:
 *    		The standard localtime function is not thread-safe because it returns a pointer to a static struct tm that is shared across threads.
 *    		This function wraps localtime with a mutex to ensure that only one thread accesses the static data at a time, making it thread-safe.
 *    		It copies the result into a user-provided struct tm buffer.
 */
















/** footnote 3: Bitwise Operations on Numerical Data
 *    These functions enable manipulation of numerical data at the bit level, allowing direct interaction with the binary representation.
 *    This is particularly useful for operations like reinterpreting data types, flipping sign bits, and mapping floating-point numbers for sorting.
 *
 *    Bitwise Reinterpretation:
 *        - double_to_uint64:
 *              Reinterprets a double as a uint64_t by copying its binary representation.
 *              This does not perform any conversion; it simply treats the bits of the double as an unsigned 64-bit integer.
 *
 *        - uint64_to_double:
 *              Reinterprets a uint64_t as a double, performing the inverse operation of double_to_uint64.
 *              Allows converting a manipulated uint64_t representation back to a double.
 *
 *    Sign Bit Manipulation:
 *        - flip_sign_bit:
 *              Flips the most significant bit (sign bit) of a 64-bit unsigned integer.
 *              This is useful for changing the sign of a floating-point number represented in binary form,
 *              while preserving its magnitude and exponent.
 *
 *    Mapping for Sorting:
 *        - double_to_mapped_uint64:
 *              Converts a double to a uint64_t in such a way that when the uint64_t values are sorted,
 *              the resulting order corresponds to the numerical order of the original double values.
 *              - For negative numbers: Invert all bits using `~u`, which reverses their order in the unsigned integer space.
 *              - For positive numbers: Flip the sign bit using `u ^= (1ULL << 63)`, shifting them into the upper half of the unsigned integer space.
 *
 *        - mapped_uint64_to_double:
 *              Reverses the mapping applied by double_to_mapped_uint64 to obtain the original double value from the mapped uint64_t.
 *              - For mapped negative numbers: Invert all bits using `~u` to restore the original negative double.
 *              - For mapped positive numbers: Flip the sign bit using `u ^= (1ULL << 63)` to restore the original positive double.
 *
 *
 *
 *
 *    Usage in Sorting:
 *            Standard integer sorting algorithms (like Radix Sort) cannot directly sort double values, especially when negative numbers are involved,
 *            because the bit patterns of floating-point numbers do not correspond to their numerical order.
 *            By mapping double values to uint64_t values using double_to_mapped_uint64, we can sort the uint64_t representations using integer sorting algorithms.
 *            The mapping ensures that the numerical order of the doubles is preserved in the unsigned integer space.
 *
 *        Process:
 *            1. Map double to uint64_t: Use double_to_mapped_uint64 to map each double in the array to a uint64_t.
 *            2. Sort the uint64_t array: Use Radix Sort or another integer sorting algorithm.
 *            3. Convert back to double: Use mapped_uint64_to_double to convert the sorted uint64_t values back to their original double representations.
 *
 *    Summary:
 *        This approach ensures that the numerical order of the doubles is preserved when they are converted to and from their uint64_t representations,
 *        allowing integer sorting algorithms like Radix Sort to handle negative and positive double values correctly.
 *        After sorting, you reverse the mapping to get back the original double values.
 */

