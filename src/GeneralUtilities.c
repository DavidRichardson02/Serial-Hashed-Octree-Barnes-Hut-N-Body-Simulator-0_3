//  GeneralUtilities.c
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02


#include "GeneralUtilities.h"
#include "CommonDefinitions.h"







int *allocate_memory_int_ptr(size_t sizeI)
{
	int *intPtr = (int*)malloc(sizeI * sizeof(int));
	if (intPtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_int_ptr'.\n");
		exit(1);
	}
	
	return intPtr;
}


float *allocate_memory_float_ptr(size_t sizeF)
{
	float *floatPtr = (float*)malloc(sizeF * sizeof(float));
	if (floatPtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_float_ptr'.\n");
		exit(1);
	}
	
	return floatPtr;
}


double *allocate_memory_double_ptr(size_t sizeD)
{
	double *doublePtr = (double*)malloc(sizeD * sizeof(double));
	if (doublePtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_double_ptr'.\n");
		exit(1);
	}
	
	return doublePtr;
}


char *allocate_memory_char_ptr(size_t sizeC)
{
	char *charPtr = (char*)malloc(sizeC * sizeof(char));
	if (charPtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_char_ptr'.\n");
		exit(1);
	}
	
	return charPtr;
}




int **allocate_memory_int_ptr_ptr(size_t sizeI)
{
	int **intPtrPtr = (int**)malloc(sizeI * sizeof(int*));
	if (intPtrPtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_int_ptr_ptr'.\n");
		exit(1);
	}
	
	
	//For each pointer, allocate memory for a float pointer
	for(size_t i = 0; i < sizeI; i++)
	{
		intPtrPtr[i] = allocate_memory_int_ptr(sizeof(intPtrPtr[i]));
	}
	
	return intPtrPtr;
}


float **allocate_memory_float_ptr_ptr(size_t sizeF)
{
	float **floatPtrPtr = (float**)malloc(sizeF * sizeof(float*));
	if (floatPtrPtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_float_ptr_ptr'.\n");
		exit(1);
	}
	
	
	//For each pointer, allocate memory for a float pointer
	for(size_t i = 0; i < sizeF; i++)
	{
		floatPtrPtr[i] = allocate_memory_float_ptr(sizeof(floatPtrPtr[i]));
	}
	
	return floatPtrPtr;
}


double **allocate_memory_double_ptr_ptr(size_t sizeD)
{
	double **doublePtrPtr = (double**)malloc(sizeD * sizeof(double*));
	if (doublePtrPtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_double_ptr_ptr'.\n");
		exit(1);
	}
	
	
	//For each pointer, allocate memory for a double pointer
	for(size_t i = 0; i < sizeD; i++)
	{
		doublePtrPtr[i] = allocate_memory_double_ptr(sizeof(doublePtrPtr[i]));
	}
	
	return doublePtrPtr;
}


char **allocate_memory_char_ptr_ptr(size_t strSize, size_t numStrings) // Allocate memory for the array of char* pointers
{
	char **charPtrPtr = (char**)malloc(numStrings * sizeof(char*));
	if (charPtrPtr == NULL)
	{
		perror("\n\nError: Unable to allocate memory in 'allocate_memory_char_ptr_ptr'.\n");
		exit(1);
	}
	
	
	//For each pointer, allocate memory for a char pointer
	for(size_t i = 0; i < numStrings; i++)
	{
		charPtrPtr[i] = allocate_memory_char_ptr(strSize);
	}
	
	return charPtrPtr;
}




void deallocate_memory_int_ptr_ptr(int **intPtrPtr, size_t numInts)
{
	if (intPtrPtr != NULL)
	{
		// First, free each element
		for (size_t i = 0; i < numInts; i++)
		{
			free(intPtrPtr[i]);
		}
		// Then, free the array of pointers
		free(intPtrPtr);
	}
}


void deallocate_memory_float_ptr_ptr(float **floatPtrPtr, size_t numFloats)
{
	if (floatPtrPtr != NULL)
	{
		// First, free each element
		for (size_t i = 0; i < numFloats; i++)
		{
			free(floatPtrPtr[i]);
		}
		// Then, free the array of pointers
		free(floatPtrPtr);
	}
}


void deallocate_memory_double_ptr_ptr(double **doublePtrPtr, size_t numDoubles)
{
	if (doublePtrPtr != NULL)
	{
		// First, free each element
		for (size_t i = 0; i < numDoubles; i++)
		{
			free(doublePtrPtr[i]);
		}
		// Then, free the array of pointers
		free(doublePtrPtr);
	}
}


void deallocate_memory_char_ptr_ptr(char **charPtrPtr, size_t numStrings)
{
	if (charPtrPtr != NULL)
	{
		// First, free each string
		for (size_t i = 0; i < numStrings; i++)
		{
			free(charPtrPtr[i]);
		}
		// Then, free the array of pointers
		free(charPtrPtr);
	}
}







/**
 * minimum
 * Calculates the minimum of two double values.
 *
 * @param a The first value.
 * @param b The second value.
 * @return The minimum value.
 */
double minimum(double a, double b)
{
	return a < b ? a : b;  // If a is less than b, return a, otherwise return b.
}


/**
 * maximum
 * Calculates the maximum of two double values.
 *
 * @param a The first value.
 * @param b The second value.
 * @return The maximum value.
 */
double maximum(double a, double b)
{
	return a > b ? a : b; //  (1.) if a is greater than b, return a.  (2.) If a is NOT greater than b, then either b is greater than a or they are equal, in either case, return b
}


/**
 * min_element
 *
 * Finds the minimum element in an array of doubles.
 * Skips NaN values to find the minimum numerical value.
 *
 * @param data The array of doubles.
 * @param n The number of elements in the array.
 * @return The minimum element, or NaN if the array is NULL or only contains NaNs.
 */
double min_element(double *data, int n)
{
	if(data == NULL)
	{
		perror("\n\nError: Data to find min element from was NULL in 'min_element'.\n");
		exit(1);
	}
	
	
	double min_elem = NAN; // Initialize to NaN
	for (int i = 0; i < n; i++) // Check from the start
	{
		if (!isnan(data[i]))
		{
			min_elem = data[i]; // Found the first non-NaN element.
			break;
		}
	}
	
	
	for(int i = 1; i<n;i++) // Start from 1 since 0 is already in min_elem
	{
		if (!isnan(data[i])) // Check if the value is not NaN
		{
			min_elem = minimum(min_elem, data[i]);
		}
	}
	
	return min_elem;
}


/**
 * max_element
 *
 * Finds the maximum element in an array of doubles.
 * Skips NaN values to find the maximum numerical value.
 *
 * @param data The array of doubles.
 * @param n The number of elements in the array.
 * @return The maximum element, or NaN if the array is NULL or only contains NaNs.
 */
double max_element(double *data, int n)
{
	if(data == NULL)
	{
		perror("\n\nError: Data to find max element from was NULL in 'max_element'.\n");
		exit(1);
	}
	
	double max_elem = NAN; // Initialize to NaN
	for (int i = 0; i < n; i++) // Check from the start
	{
		if (!isnan(data[i]))
		{
			max_elem = data[i]; // Found the first non-NaN element.
			break;
		}
	}
	
	
	for(int i = 1; i < n; i++) // Continue from where we left off
	{
		if (!isnan(data[i])) // Check if the value is not NaN
		{
			max_elem = maximum(max_elem, data[i]);
		}
	}
	
	return max_elem;
}


int max_element_int(int *data, int n)
{
	if(data == NULL)
	{
		perror("\n\nError: Data to find max element from was NULL in 'max_element_int'.\n");
		exit(1);
	}
	
	int max_elem = NAN; // Initialize to the smallest possible integer value
	for (int i = 0; i < n; i++) // Check from the start
	{
		if (!isnan(data[i]))
		{
			max_elem = data[i]; // Found the first non-NaN element.
			break;
		}
	}
	
	
	for(int i = 1; i < n; i++) // Continue from where we left off
	{
		if (!isnan(data[i])) // Check if the value is not NaN
		{
			max_elem = maximum(max_elem, data[i]);
		}
	}
	
	
	return max_elem;
}

/**
 * sum_elements
 *
 * Sums the elements in an array of doubles.
 * @param data The array of doubles.
 * @param n The number of elements in the array.
 * @return The maximum element, or NaN if the array is NULL or only contains NaNs.
 */
double sum_elements(double *data, int n)
{
	if(data == NULL)
	{
		perror("\n\nError: Data to sum elements from was NULL in 'sum_elements'.\n");
		exit(1);
	}
	
	double sum = 0.0;
	for(int i = 0; i < n; i++)
	{
		sum += data[i];
	}
	
	return sum;
}






/**
 * convert_to_unix_time
 *
 * Converts a date/time string into Unix time(the number of seconds since the Unix Epoch, January 1, 1970).
 * It attempts to parse the string using various common date/time formats and returns the Unix time if successful.
 * The purpose of this function is to help in standardizing data set file contents.
 *
 * @param dateTimeString A pointer to the string containing date/time information.
 * @return Unix time as time_t. Returns -1 if conversion fails.
 */
time_t convert_to_unix_time(const char *dateTimeString)
{
	struct tm tm; // Structure to hold the broken-down time.
	char *parsed;  // A pointer to track where the parsing of the date/time string ended.
	time_t unixTime = -1;
	
	
	// Iterate through each date/time format specified in commonDateTimeFormats.
	for (int i = 0; i < 12; ++i)
	{
		memset(&tm, 0, sizeof(struct tm)); // Reset the tm structure for each iteration.
		
		// Parse the input string according to these formats.
		parsed = strptime(dateTimeString, commonDateTimeFormats[i], &tm);
		
		
		// Check if parsing was successful and the entire string was consumed.
		if (parsed != NULL && *parsed == '\0')
		{
			// Convert the parsed time (tm structure) to Unix time.
			unixTime = mktime(&tm);
			
			// If conversion is successful (unixTime is not -1), break out of the loop.
			if (unixTime != -1)
			{
				break;
			}
		}
	}
	
	// Return the resulting Unix time, or -1 if none of the formats matched.
	return unixTime;
}


/**
 * thread_safe_localtime
 *
 * Converts time_t to tm as Local Time in a thread-safe manner.
 * Uses mutexes to ensure thread safety during the conversion.
 *
 * @param tim The time_t structure to convert.
 * @param result A pointer to the struct tm where the result will be stored.
 * @return A pointer to the result if successful, NULL otherwise.
 */
struct tm *thread_safe_localtime(const time_t *tim, struct tm *result)
{
	struct tm *t = NULL;
	
	// Ensure 'result' is a non-null pointer
	if (result == NULL)
	{
		return NULL;
	}
	
	
	pthread_mutex_lock(&localtime_mutex);
	t = localtime(tim);
	if (t)
	{
		*result = *t;  // Copy the statically allocated struct
	}
	pthread_mutex_unlock(&localtime_mutex);
	
	pthread_mutex_destroy(&localtime_mutex);
	
	
	return t ? result : NULL; // Return result if localtime didn't return NULL
}








/**
 * flip_sign_bit
 *
 * Flip the Sign Bit of a 64-bit Integer.
 * Flips(toggles) the sign bit(Most Significant Bit) of a 64-bit integer, effectively changing the sign of a floating-point number
 * represented in its binary form(specifically IEEE 754 representation for doubles). This function is particularly
 * useful in scenarios involving floating-point numbers where the sign needs to be altered while preserving the
 * magnitude and exponent parts of the floating-point representation.
 *
 * @param value The 64-bit unsigned integer whose sign bit is to be flipped.
 * @return The 64-bit unsigned integer with its sign bit flipped.
 */
uint64_t flip_sign_bit(uint64_t value)
{
	// '1ull << 63' creates a 64-bit integer with only the MSB set to 1 (i.e., the 63rd bit in zero-indexed notation).
	
	// 'value ^ ...' applies the bitwise XOR operation between the input 'value' and the above-created integer.
	//   XOR with '1' toggles the corresponding bit, so if the sign bit in 'value' is 0 (positive number), it becomes 1 (negative number),
	//   and vice versa, hence this operation flips the sign bit of the input 'value'.
	return value ^ (1ull << 63);
}


/**
 * double_to_uint64
 *
 * This function reinterprets a double value as a 64-bit unsigned integer (uint64_t) by directly copying its binary representation.
 * This method of type-punning allows the bit pattern of a double to be reinterpreted as a uint64_t without any conversion,
 * facilitating direct manipulation of a double's binary representation for bit-level operations, useful in certain algorithms like
 * sorting or encoding, where floating-point operations are not required.
 *
 * @param value The double value to be reinterpreted as a uint64_t.
 * @return The binary representation of the input double as a uint64_t.
 */
uint64_t double_to_uint64(double value)
{
	uint64_t result;
	
	// Copies the binary representation of 'value' into 'result', 'sizeof(double)' ensures that exactly 8 bytes (size of a double) are copied.
	memcpy(&result, &value, sizeof(double));
	return result;
}

/**
 * uint64_to_double
 *
 * This function reinterprets a 64-bit unsigned integer (uint64_t) as a double by copying its binary representation.
 * It is the inverse operation of double_to_uint64, allowing for the conversion of a manipulated uint64_t representation
 * back to a double. Allows for conversion of a uint64_t value back to a double after performing bit-level manipulations, enabling the use
 *  or interpretation of the result as a floating-point number.
 *
 * @param value The 64-bit unsigned integer to be reinterpreted as a double.
 * @return The binary representation of the input uint64_t as a double.
 */
double uint64_to_double(uint64_t value)
{
	double result;
	
	// Copies the binary representation of 'value' into 'result', sizeof(uint64_t)' ensures that exactly 8 bytes (size of a uint64_t) are copied.
	memcpy(&result, &value, sizeof(uint64_t));
	return result;
}


/**
 Explanation:
 
 Mapping Negative Numbers: By inverting all bits of negative numbers (u = ~u), you reverse their order in the unsigned integer space, so that more negative numbers come before less negative ones when sorted.
 Mapping Positive Numbers: By flipping the sign bit of positive numbers (u ^= (1ULL << 63)), you shift them into the upper half of the unsigned integer space, preserving their order.
 Sorting: The mapped uint64_t values can be sorted using a standard radix sort for unsigned integers.
 Restoring Original Values: After sorting, you reverse the mapping to get back the original double values.
 Note: This approach ensures that the numerical order of the doubles is preserved when they are converted to and from their uint64_t representations, allowing the radix sort to handle negative values correctly.
 
 The key is to map the bit patterns of the doubles to unsigned integers in such a way that the numerical order of the doubles corresponds to the unsigned integer order after mapping.
 */

/**
 * double_to_mapped_uint64
 *
 * Converts a double value to a mapped uint64_t representation to facilitate correct sorting order of both negative and positive doubles as unsigned integers.
 * Mapping doubles to uint64_t allows the use of radix sort, which sorts numbers based on their bit patterns, ensuring that their numerical ordering is maintained correctly.
 * This function maps double values to uint64_t such that:
 *		- Negative numbers are handled by inverting all bits, placing more negative values before less negative ones in the unsigned space.
 *		- Positive numbers are handled by flipping the MSB, shifting them into the upper half of the uint64_t space to preserve their order relative to each other.
 *
 *	@param value The double value to be converted.
 *	@return The mapped uint64_t representation of the input double.
 */
uint64_t double_to_mapped_uint64(double value)
{
	uint64_t u = double_to_uint64(value);  // Convert double to uint64_t without changing bit pattern
	if (u & (1ULL << 63)) // Check if the number is negative (MSB is 1)
	{
		u = ~u; // Invert all bits for negative numbers
	}
	else // Positive number
	{
		u ^= (1ULL << 63); // Flip sign bit for positive numbers
	}
	return u;
}


/**
 * mapped_uint64_to_double
 *
 * Converts a mapped uint64_t back to its original double representation after sorting.
 * This function reverses the mapping applied in double_to_mapped_uint64:
 * 			- If the MSB is 1, it indicates a positive number in the mapped space, requiring flipping the MSB to restore the original positive double.
 * 			- If the MSB is 0, it indicates an originally negative number, requiring all bits to be inverted to restore the original negative double.
 *
 * @param u The mapped uint64_t value to be converted back to double.
 * @return The original double value from the mapped uint64_t representation.
 */
double mapped_uint64_to_double(uint64_t u)
{
	if (u & (1ULL << 63)) // MSB is 1, indicating a mapped positive number
	{
		u ^= (1ULL << 63); // Flip the sign bit to restore the original positive double
	}
	else // MSB is 0, indicating a mapped negative number
	{
		u = ~u; // Invert all bits to restore the original negative double
	}
	
	
	return uint64_to_double(u); // Convert the uint64_t back to double
}











/**
 * merge_data
 *
 * Merges two subarrays of unsortedData[].
 * This function is a part of the merge sort algorithm. It merges two sorted subarrays
 * defined by the indices [left, middle] and [middle+1, right] into a single sorted array.
 * It uses two temporary arrays to hold the subarrays and then merges them back into the
 * original array in a sorted manner.
 *
 * @param unsortedData Pointer to the array of doubles to be sorted.
 * @param left The starting index of the first subarray,  unsortedData[left..middle]
 * @param middle The ending index of the first subarray and the starting index of the second subarray minus one.
 * @param right The ending index of the second subarray,  unsortedData[middle+1..right]
 */
void merge_data(double *unsortedData, int left, int middle, int right)
{
	// Calculate the sizes of the two subarrays to be merged
	int i, j, k;
	int size1 = middle - left + 1;  // Size of the first subarray
	int size2 = right - middle;   // Size of the second subarray
	
	
	
	// Dynamically allocate memory for temporary arrays holding subarray elements
	double *leftTempData = allocate_memory_double_ptr(size1);
	double *rightTempData = allocate_memory_double_ptr(size2);
	
	
	// Copying data to the temporary arrays to isolate the subarrays
	for (i = 0; i < size1; i++)
	{
		leftTempData[i] = (unsortedData)[left + i];
	}
	for (j = 0; j < size2; j++)
	{
		rightTempData[j] = (unsortedData)[middle + 1 + j];
	}
	
	
	
	// Merging the temporary arrays back into the original array[left..right]
	i = 0; // Initial index of the first subarray
	j = 0; // Initial index of the second subarray
	k = left; // Initial index of the merged subarray
	while (i < size1 && j < size2)
	{
		if (leftTempData[i] <= rightTempData[j])
		{
			(unsortedData)[k] = leftTempData[i];
			i++;
		}
		else
		{
			(unsortedData)[k] = rightTempData[j];
			j++;
		}
		k++;
	}
	while (i < size1)  // Copy the remaining elements from leftTempData to unsortedData
	{
		(unsortedData)[k] = leftTempData[i];
		i++;
		k++;
	}
	while (j < size2)  // Copy the remaining elements from rightTempData to unsortedData
	{
		(unsortedData)[k] = rightTempData[j];
		j++;
		k++;
	}
	
	
	free(leftTempData);
	free(rightTempData);
}


/**
 * merge_sort_data
 *
 * This function implements the recursive part of the merge sort algorithm. It recursively divides
 * the array into two halves, calls itself for the two halves, and then merges the two sorted halves
 * using the merge_data function..
 *
 * @param unsortedData Pointer to the array of doubles to be sorted.
 * @param left The starting index of the array portion to be sorted.
 * @param right The ending index of the array portion to be sorted.
 */
void merge_sort_data(double *unsortedData, int left, int right)
{
	if(left < right)  // Check to ensure valid indexing bounds
	{
		// Find the middle point to divide the array into two halves
		int middle = left + (right - left) / 2; // Same as (left+right)/2, but avoids overflow for large left and right
		
		// Recursive call to sort the first half of the array until the sub-arrays are small enough to be solved directly
		merge_sort_data(unsortedData, left, middle);
		
		// Recursive call to sort the second half of the array until the sub-arrays are small enough to be solved directly
		merge_sort_data(unsortedData, middle + 1, right);
		
		// Merge the sorted halves
		merge_data(unsortedData, left, middle, right);
	}
}


/**
 * merge_sort,  a divide and conquer algorithm, O(nLogn)
 *
 * This function is the entry point for the merge sort algorithm. It checks for
 * null pointers and then calls merge_sort_data to sort the entire array.
 *
 * @param unsortedData Pointer to the array of doubles to be sorted.
 * @param numElements The number of elements in the array.
 */
void merge_sort(double *unsortedData, const int numElements)
{
	// Check for null pointers to ensure data integrity
	if(unsortedData == NULL)
	{
		perror("\n\nError: Data to be sorted was NULL in 'merge_sort'.\n");
		exit(1);
	}
	
	
	// Initial call to the recursive merge sort function on the entire array
	merge_sort_data(unsortedData, 0, numElements - 1);
}




/**
 * radix_sort_doubles
 *
 * Sorts an array of double precision floating-point numbers using Radix Sort algorithm. This implementation
 * specifically deals with the floating-point nature of the data by converting doubles to their 64-bit integer
 * representations. The sorting is performed on these integer representations, enabling the Radix Sort algorithm
 * to be applied to floating-point numbers.
 *
 *
 * Perform a radix sort on an array of double values. This function sorts an array of doubles
 * using a radix sort algorithm, which is a non-comparative integer sorting algorithm. It achieves this by
 * interpreting the bit representation of double values as unsigned 64-bit integers (uint64_t).
 *
 * The radix sort is done on the binary representation of these integers. Special handling is done for negative
 * numbers by flipping the sign bit, to ensure they are sorted correctly. Post sorting, the binary representations
 * are converted back to double values.
 *
 * This function is particularly useful in scenarios where a fast, stable sorting of a large number of floating-point numbers( > 60) is required.
 *
 * @param unsortedData A pointer to the array of double values to be sorted.
 * @param numElements The number of elements in the array.
 */
void radix_sort_doubles(double *unsortedData, const int numElements)
{
	// Allocate memory for temporary arrays used in sorting
	uint64_t* temp = (uint64_t*)malloc(numElements * sizeof(uint64_t));
	uint64_t* intValues = (uint64_t*)malloc(numElements * sizeof(uint64_t));
	
	/// Step 1: Map doubles to uint64_t for sorting, reinterpreting each double as a 64-bit integer.
	for (size_t i = 0; i < numElements; i++)
	{
		//
		/*
		 uint64_t intValue = double_to_uint64(unsortedData[i]);
		 // Check if the number is negative (indicated by the sign bit). If so, flip the sign bit to ensure correct sorting order
		 if (intValue >> 63)
		 {
		 intValue = flip_sign_bit(intValue);
		 }
		 intValues[i] = intValue;
		 //*/
		intValues[i] = double_to_mapped_uint64(unsortedData[i]);
	}
	
	
	
	
	/// Step 2: Perform Radix Sort on the mapped 64-bit integer representations.
	for (int shift = 0; shift < 64; shift += 8)
	{
		// Sorts the array byte-by-byte (8 bits at a time).
		size_t counts[256] = {0};
		size_t offsets[256] = {0};
		
		
		// Counting occurrences for each byte value
		for (size_t i = 0; i < numElements; i++)
		{
			size_t index = (intValues[i] >> shift) & 0xFF;
			counts[index]++;
		}
		
		// Calculate offsets for each byte value
		for (size_t i = 1; i < 256; i++)
		{
			offsets[i] = offsets[i - 1] + counts[i - 1];
		}
		
		// Reorder elements based on the current byte
		for (size_t i = 0; i < numElements; i++)
		{
			size_t index = (intValues[i] >> shift) & 0xFF;
			temp[offsets[index]] = intValues[i];
			offsets[index]++;
		}
		
		// Copy the sorted temporary array to the original one
		memcpy(intValues, temp, numElements * sizeof(uint64_t));
	}
	
	
	/// Step 3: Convert uint64_t back to double
	for (size_t i = 0; i < numElements; i++)
	{
		//
		/*
		 // Re-flip sign bit if necessary to restore the original sign.
		 if (intValues[i] >> 63)
		 {
		 intValues[i] = flip_sign_bit(intValues[i]);
		 }
		 unsortedData[i] = uint64_to_double(intValues[i]);
		 //*/
		unsortedData[i] = mapped_uint64_to_double(intValues[i]);
	}
	
	free(temp);
	free(intValues);
}








/**
 * set_memory_block
 *
 * This function sets the first 'n' bytes of the memory area pointed to by 'block' to the byte 'c'.
 * It optimizes the process by setting word-sized chunks of memory when possible.
 * The function first handles the alignment by setting bytes until the block pointer is aligned to the word size.
 * Once aligned, it sets memory in word-sized chunks (32-bit or 64-bit, depends on the system).
 * Lastly, it handles any remaining bytes that don't fit into a full word.
 *
 * @param block The memory block to be set.
 * @param c The byte value to be set.
 * @param n The number of bytes to set.
 * @return A pointer to the memory block.
 */
void *set_memory_block(void *block, int c, size_t n)
{
	// Cast the pointer to an unsigned char pointer for byte-wise operations
	unsigned char *p = (unsigned char *)block;
	
	
	
	// Create a word-sized value filled with the byte value 'c'
	uintptr_t wordValue = (unsigned char)c;
	
	
	
	// Fill the wordValue with the byte 'c' in all its bytes. Bitwise operations to shift the bits of wordValue to the left, effectively multiplying wordValue by 256 (2^8), 65536 (2^16) and 4294967296 (2^32), if system is 64-bit.
	wordValue |= wordValue << 8;
	wordValue |= wordValue << 16;
#if UINTPTR_MAX > 0xffffffff // Check if the system is 64-bit
	wordValue |= wordValue << 32;
#endif
	
	
	
	// Align the destination to the word size
	while (n > 0 && ((uintptr_t)p & (sizeof(uintptr_t) - 1)) != 0) // If the result is not 0, the address is not aligned to the size of a uintptr_t.
	{
		*p++ = (unsigned char)c; // Assigns the value of c cast to unsigned char to the memory location pointed to by p
		n--; // Increments p to point to the next byte
	}
	
	
	
	// Fill memory in word-sized chunks
	uintptr_t *wp = (uintptr_t *)p; // This casts p to a uintptr_t pointer and assigns it to wp
	size_t words = n / sizeof(uintptr_t); // Calculates the number of uintptr_t-sized chunks that can be filled in the remaining bytes
	while (words--)
	{
		*wp++ = wordValue; // This assigns wordValue to the uintptr_t-sized chunk of memory pointed to by wp, then increments wp by sizeof(uintptr_t), moving it to the next chunk
	}
	
	
	
	// Handle any remaining bytes
	p = (unsigned char *)wp;
	n = n % sizeof(uintptr_t);
	while (n--)
	{
		*p++ = (unsigned char)c;
	}
	
	return block;
}




/**
 * copy_memory_block
 *
 * Copies the specified number of bytes from the source memory block to the destination memory block.
 * The function first handles the alignment by copying bytes until the destination pointer is aligned to the word size.
 * Once aligned, it copies memory in word-sized chunks(32-bit or 64-bit, depends on system, then lastly,
 * it handles any remaining bytes don't fit into a full word
 * It is designed to handle unaligned memory access and optimize the copying process by copying aligned
 * memory blocks of the data.
 *
 * @param destination The destination memory block where the data will be copied.
 * @param source The source memory block from which the data will be copied.
 * @param n The number of bytes to copy from the source to the destination.
 * @return A pointer to the destination memory block.
 */
void *copy_memory_block(void *destination, const void *source, size_t n)
{
	// Cast the pointers to unsigned char pointers for byte-wise copying
	unsigned char *dst = (unsigned char *)destination;
	const unsigned char *src = (const unsigned char *)source;
	
	
	// If the source and destination are the same or n is zero, there's nothing to do
	if (dst == src || n == 0)
	{
		return destination;
	}
	
	
	
	/// Create word-sized type for copying data in word-sized chunks
	typedef uintptr_t word; // Unsigned integer type capable of storing a pointer
	const size_t wordSize = sizeof(word); // Determines the size of the word.
	const size_t wordMask = wordSize - 1; // Mask to check if the pointer is aligned to the word size
	
	
	
	// Copy bytes until the destination is aligned to the word size
	while (n > 0 && ((uintptr_t)dst & wordMask) != 0)
	{
		*dst++ = *src++;
		n--;
	}
	
	
	// Copy word-sized chunks of data
	word *wd = (word *)dst; // Cast the destination pointer to a word pointer
	const word *ws = (const word *)src; // Cast the source pointer to a word pointer
	size_t words = n / wordSize; // Calculate the number of whole words to copy
	while (words--)
	{
		*wd++ = *ws++; // Copy word-sized chunks
	}
	
	
	// Cast the pointers back to unsigned char pointers for byte-wise copying of remaining bytes
	dst = (unsigned char *)wd;
	src = (const unsigned char *)ws;
	
	
	// Calculate the number of remaining bytes to be copied from the source to the destination after copying whole words of data(wordSize bytes at a time).
	n = n % wordSize;
	
	
	// Copy the remaining bytes
	while (n--)
	{
		*dst++ = *src++;
	}
	return destination;
}







