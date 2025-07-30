//  MortonKeys.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3

/*
 * MortonKey code: for Spatial Encoding and Decoding. Used to sort elements according to the order in which they
 * would be visited by a space-filling curve.
 
 * Description:
 * This file defines the MortonKey code, which is used to encode and decode
 * 3D positions into a Morton key (a linear indexing system that preserves
 * spatial locality) and sort them into Morton Ordering. This technique is particularly useful for optimizing algorithms
 * requiring spatial coherence by reducing cache misses and improving data access patterns.
 * Performance Notes:
 *     - The "Magic Bits" method for bit shifting is used for the encoding/decoding of Morton keys.
 *     - MergeSort is faster than BubbleSort for sorting the keys.
 *     - Radix sorting is optimized
 */

/**
 The Morton key computation involves the following steps:
 
		- Coordinate Mapping and Scaling: Map the 3D coordinates from the domain (-domainSize, domainSize) to a non-negative range (0, dimKeyMax).
 
		- Conversion to Integers: Convert the scaled coordinates to integers.
 
		- Bit Interleaving (Morton Encoding): Interleave the bits of the integer coordinates to create a single Morton key.
 
		- Key Formation: Combine the interleaved bits to form the final Morton key.
 

 */
#pragma once
#include "MathUtils.hpp"
#include "SequenceContainers.hpp"
#include "ofMain.h"




// ------------- Constants and typedefs  -------------
typedef uint64_t spatialKey;
static const int MortonKeyDim = 21; //Use 21 bits for each dimension, resulting in a 63-bit Morton code (21*3 = 63).
static const spatialKey maxKeyDimension = (1ull << 21);
static const double maxKeyDimension_d = (double)(1ull << 21);
static const int NUM_BINS = 256;  // Number of bins used in binning optimization




// ------------- Masks for the morton bit swizzle. -------------
static const spatialKey sepMasks[] =
{
	0x1f00000000ffff,// shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
	0x1f0000ff0000ff,// shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
	0x100f00f00f00f00f,// shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
	0x10c30c30c30c30c3,// shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
	0x1249249249249249
};




// ------------- Helper functions for MortonKey Encoding and Decoding -------------
static inline spatialKey interleaveBits(spatialKey xInt, spatialKey yInt, spatialKey zInt); //Interleave the bits of values representing the integer coordinates
static inline spatialKey computeMortonKey(const Vec3D& _position, const double _size);   // Compute the Morton key from integer coordinates
static inline Vec3D decodeMortonKey(const spatialKey _mortonKey, const double _size); //Decoding Function




// ------------- Helper functions for MortonKey sorting -------------
static inline void Merge(spatialKey* keys, size_t left, size_t middle, size_t right); // Merges two sub-arrays of keys
static inline void MergeSort(spatialKey* keys, size_t left, size_t right); // Recursive merge sort function
static inline void MergeSortMortonKeys(const size_t numKeys, spatialKey* unsortedKeys); // Merge sorting Morton keys
static inline void RadixSortMortonKeys(spatialKey* keys, size_t numKeys); //Radix sorting of morton keys
static inline void RadixSortPass(spatialKey* keys, size_t* indexes, size_t numKeys, spatialKey* tempKeys, size_t* tempIndexes, int shift, size_t* histogram, size_t* offsets); // Radix sort pass on specified bits
static inline void ThreePassRadixSortMortonKeys(spatialKey* keys, size_t* indexes, size_t numKeys); // Three-pass radix sort with binning optimization





/**
 * Interleave the bits of three spatialKey values representing the integer coordinates in 3D space.
 *
 * This function interleaves the bits of xInt, yInt, and zInt to produce a single 63-bit Morton key, which is a linear index for 3D points that preserves spatial locality.
 * The interleaving of bits is achieved through a series of bit-wise manipulations, essentially spreading the bits of each integer across the entire 63-bit range.
 *
 *  @param xInt: X-coordinate integer representation.
 *  @param yInt: Y-coordinate integer representation.
 *  @param zInt: Z-coordinate integer representation.
 *
 *  @return Morton key interleaved from xInt, yInt, and zInt.
 */
inline spatialKey interleaveBits(spatialKey xInt, spatialKey yInt, spatialKey zInt)
{
	spatialKey _mortonKey = 0;
	
	//morton code swizzle
	xInt = (xInt | xInt << 32) & sepMasks[0];
	xInt = (xInt | xInt << 16) & sepMasks[1];
	xInt = (xInt | xInt << 8) & sepMasks[2];
	xInt = (xInt | xInt << 4) & sepMasks[3];
	xInt = (xInt | xInt << 2) & sepMasks[4];
	
	yInt = (yInt | yInt << 32) & sepMasks[0];
	yInt = (yInt | yInt << 16) & sepMasks[1];
	yInt = (yInt | yInt << 8) & sepMasks[2];
	yInt = (yInt | yInt << 4) & sepMasks[3];
	yInt = (yInt | yInt << 2) & sepMasks[4];
	
	zInt = (zInt | zInt << 32) & sepMasks[0];
	zInt = (zInt | zInt << 16) & sepMasks[1];
	zInt = (zInt | zInt << 8) & sepMasks[2];
	zInt = (zInt | zInt << 4) & sepMasks[3];
	zInt = (zInt | zInt << 2) & sepMasks[4];
	

	//cout << "xInt: " << xInt << " yInt: " << yInt << " zInt: " << zInt << endl;
	
	// Capture and configure the 64-bit key to have a left most bit(MSB) as 1, followed by the 21-bit Morton encoding for each of the x, y, and z 3D-position components.
	_mortonKey = (1ULL << 63) | (xInt << 2) | (yInt << 1) | zInt;
	return(_mortonKey);
}
/**
 * Compute the Morton key for a given 3D position and spatial size (_size)..
 *
 * This function first scales and converts the floating-point 3D position
 * values into the integer domain. After this conversion, it calls the `interleaveBits` method to create the Morton key
 * which is stored as the instance variable.
 *
 *  @param _position: 3D position to be encoded into a Morton key.
 *  @param _size: Maximum spatial extent of the positions. This value is used to scale the positions to [0, 1) range.
 *
 *  @return Encoded morton key
 */
inline spatialKey computeMortonKey(const Vec3D& _position, const double _size) //computeMortonCode based on the position of a point within a bounding box _size
{
	// Rescale the coordinates to [0, 1) relative to the relationship between the position within the
	// The effective origin is known to be located at (0,0)
	// By accounting for the size of the space enclosing these positions(by adding it to the position components) we are able to effectively treat these positions with
	// The coordinates are shifted by domainSize to map from (-domainSize, domainSize) to (0, 2*domainSize).
	double xScaled = _position.x + _size;
	double yScaled = _position.y + _size;
	double zScaled = _position.z + _size;
	
	
	// The coordinates are scaled to the range (0, dimKeyMax).
	double frac = maxKeyDimension_d / (2.0 * _size);
	xScaled *= frac;
	yScaled *= frac;
	zScaled *= frac;
	
	
	
	//std::cout << "\n\nVec3D: (" << _position.x << ", " << _position.y << ", " << _position.z << ")" << std::endl;
	//cout << "xScaled: " << xScaled << " yScaled: " << yScaled << " zScaled: " << zScaled << endl;
	
	
	
	// Convert to integer coordinates in the range [0, 2^bits)
	spatialKey xInt = (spatialKey)(xScaled); //std::min((spatialKey)(xScaled), maxKeyDimension - 1);
	spatialKey yInt = (spatialKey)(yScaled); // std::min((spatialKey)(yScaled), maxKeyDimension - 1);
	spatialKey zInt = (spatialKey)(zScaled); // std::min((spatialKey)(zScaled), maxKeyDimension - 1);
	
	
	//cout << "xInt: " << xInt << " yInt: " << yInt << " zInt: " << zInt << endl;
	
	
	spatialKey _mortonKey = interleaveBits(xInt, yInt, zInt);
	
	
	//std::cout << "Morton Key: " << _mortonKey;
	//std::bitset<64> binary(_mortonKey);
	//cout << "\nBinary MortonKey: " << binary;
	
	
	
	
	return(_mortonKey);
}



/**
 * Takes a Morton key and decodes it back into a 3D position.
 *
 * The function decodeMortonKey is designed to reverse operation performed by interleaveBits, effectively 'de-interleaving' the bits of the Morton key back into three separate x, y, and z components.
 * These integer values are then scaled back into the floating-point domain to retrieve the original 3D position.
 *
 * A loop is performed over the number of bits used for each dimension(21). For each bit index i, we extract the i-th bit from the Morton key for each dimension using bitwise AND (&) and shifting (<<, >>).
 * This is done by: bitwise ANDing (&) the Morton key with a mask that has a 1 in the position corresponding to the i-th bit of the given dimension and 0s everywhere else.
 * For x, this is (1ull << ((3*i) + 2)) - i, for y, this is (1ull << ((3*i) + 1)) - i, and for z, this is (1ull << (3*i)) - i. This operation will yield a number
 * where only the i-th bit of the dimension is potentially non-zero. Bitwise ORing (|=) the result with the accumulated result for the dimension. This 'keeps' the bits that were already set in the dimension
 * and potentially sets the i-th bit if it was set in the Morton key.
 *
 * It then scales these integer values back into the floating-point domain to retrieve the original 3D position.
 *
 *  @param _mortonKey: Morton key to be decoded.
 *  @param _size: Maximum spatial extent of the positions. This value is used to scale back the decoded position to its original range.
 *
 *  @return Decoded 3D position from the Morton key.
 */
inline Vec3D decodeMortonKey(const spatialKey _mortonKey, const double _size)
{
	spatialKey xInt = 0, yInt = 0, zInt = 0;
	for (int i = 0; i < MortonKeyDim; ++i)
	{
		//xInt |= (_mortonKey & ((uint64_t)1 << (3 * i + 2))) >> ((3 * i) + 2) - i;
		//yInt |= (_mortonKey & ((uint64_t)1 << (3 * i + 1))) >> ((3 * i) + 1) - i;
		//zInt |= (_mortonKey & ((uint64_t)1 << (3 * i))) >> (3 * i) - i;
		
		xInt |= ((_mortonKey >> (3 * i + 2)) & 1ULL) << i;
		yInt |= ((_mortonKey >> (3 * i + 1)) & 1ULL) << i;
		zInt |= ((_mortonKey >> (3 * i))     & 1ULL) << i;
	}
	double xScaled = (double)(xInt);
	double yScaled = (double)(yInt);
	double zScaled = (double)(zInt);
	xScaled = xScaled / maxKeyDimension_d * 2.0 * _size - _size;
	yScaled = yScaled / maxKeyDimension_d * 2.0 * _size - _size;
	zScaled = zScaled / maxKeyDimension_d * 2.0 * _size - _size;
	
	
	
	
	return Vec3D(xScaled, yScaled, zScaled);
}








static inline void Merge(spatialKey* keys, size_t left, size_t middle, size_t right)
{
	size_t i, j, k;
	size_t size1 = middle - left + 1;
	size_t size2 = right - middle;
	// Create temp arrays
	spatialKey* leftTemp = new spatialKey[size1];
	spatialKey* rightTemp = new spatialKey[size2];
	// Copy data to temp arrays
	for (i = 0; i < size1; i++)
	{
		leftTemp[i] = keys[left + i];
	}
	for (j = 0; j < size2; j++)
	{
		rightTemp[j] = keys[middle + 1 + j];
	}
	// Merge the temp arrays back into keys[left..right]
	i = 0;
	j = 0;
	k = left;
	while (i < size1 && j < size2)
	{
		if (leftTemp[i] <= rightTemp[j])
		{
			keys[k] = leftTemp[i];
			i++;
		}
		else
		{
			keys[k] = rightTemp[j];
			j++;
		}
		k++;
	}
	while (i < size1)  // Copy the remaining elements of leftTemp
	{
		keys[k] = leftTemp[i];
		i++;
		k++;
	}
	while (j < size2)    // Copy the remaining elements of rightTemp
	{
		keys[k] = rightTemp[j];
		j++;
		k++;
	}
	// Clean up
	delete[] leftTemp;
	delete[] rightTemp;
}
static inline void MergeSort(spatialKey* keys, size_t left, size_t right)
{
	if (left < right)
	{
		size_t middle = left + (right - left) / 2;  // Same as (left+right)/2, but avoids overflow for large left and right
													// Sort first and second halves
		MergeSort(keys, left, middle);
		MergeSort(keys, middle + 1, right);
		Merge(keys, left, middle, right);
	}
}
static inline void MergeSortMortonKeys(const size_t numKeys, spatialKey* unsortedKeys)
{
	if (unsortedKeys == nullptr)
	{
		throw std::invalid_argument("unsortedKeys is a null pointer");
	}
	MergeSort(unsortedKeys, 0, numKeys - 1);
}








/**
 * Radix sort is a non-comparative integer sorting algorithm that sorts numbers
 * digit by digit, starting from the least significant digit (LSD) to the most
 * significant digit (MSD) or vice-versa. Used here for Morton keys.
 */
static inline void RadixSortMortonKeys(spatialKey* keys, size_t numKeys)
{
	spatialKey* temp = new spatialKey[numKeys];
	for (int shift = 0; shift < 64; shift += 8)
	{
		size_t counts[256] = { 0 };
		size_t offsets[256] = { 0 };
		for (size_t i = 0; i < numKeys; i++)  // Counting occurrences
		{
			size_t index = (keys[i] >> shift) & 0xFF;
			counts[index]++;
		}
		for (size_t i = 1; i < 256; i++)   // Calculate offsets
		{
			offsets[i] = offsets[i - 1] + counts[i - 1];
		}
		for (size_t i = 0; i < numKeys; i++)    // Reorder
		{
			size_t index = (keys[i] >> shift) & 0xFF;
			temp[offsets[index]] = keys[i];
			offsets[index]++;
		}
		for (size_t i = 0; i < numKeys; i++)    // Copy the temporary array to the original one
		{
			keys[i] = temp[i];
		}
	}
	delete[] temp;
}








/**
 * Radix sort a single pass on the specified bits.
 *
 * This function handles a single pass of radix sort based on the bits indicated by the 'shift' parameter.
 *
 * It uses a histogram-based approach, calculating the occurrence of
 * each value in the histogram, then using these values to determine the
 * correct position of each key in the output.
 *
 * @param keys       The Morton keys to be sorted.
 * @param numKeys    The number of keys.
 * @param temp       Temporary storage for the sorted keys.
 * @param shift      The number of bits to shift to access the target bits.
 * @param histogram  The histogram storing the count of occurrences for each value.
 * @param offsets    The calculated positions for placing each key in the output.
 */
static inline void RadixSortPass(spatialKey* keys, size_t* indexes, size_t numKeys, spatialKey* tempKeys, size_t* tempIndexes, int shift, size_t* histogram, size_t* offsets)
{
	// Initialize histogram to zeros
	for (int i = 0; i < NUM_BINS; ++i)
	{
		histogram[i] = 0;
	}
	
	// Populate the histogram based on the values of the current bits, i.e., the bits at the 'shift' position
	// This counts how many keys belong to each bin.
	for (size_t i = 0; i < numKeys; ++i)
	{
		size_t index = (keys[i] >> shift) & 0xFF; // Extracting the relevant byte from Morton key
		histogram[index]++;
	}
	
	// Compute offsets for each bin. This determines the starting index in the sorted array for each bin.
	offsets[0] = 0;
	for (size_t i = 1; i < NUM_BINS; ++i)
	{
		offsets[i] = offsets[i - 1] + histogram[i - 1];
	}
	
	// Reorder the Morton keys and their corresponding indexes according to the current byte values based on the histogram and offsets
	for (size_t i = 0; i < numKeys; ++i)
	{
		size_t index = (keys[i] >> shift) & 0xFF;
		tempKeys[offsets[index]] = keys[i];
		tempIndexes[offsets[index]] = indexes[i]; // Sort the indexes in the same way as the keys
		offsets[index]++;
	}
	
	// Copy sorted keys from the temp array back to the original Morton keys array
	for (size_t i = 0; i < numKeys; ++i)
	{
		keys[i] = tempKeys[i];
		indexes[i] = tempIndexes[i];
	}
}




/**
 * Three-pass radix sort with binning optimization.
 *
 * This function performs a three-pass radix sort on the Morton keys.
 * Since we know that the precedence is x > y > z, we modify the
 * radix sort to prioritize bits from the x coordinate,
 * followed by bits from the y coordinate, and then bits from the z coordinate.
 * This ensures the Morton keys are sorted by x, then y, and then z.
 *
 * The function reduces cache misses by splitting the Morton code into three parts and sorting each part separately.
 * Because our Morton code is 63 bits long (manually offset by 1), each pass sorts based on 63/3 bits.
 * The first pass sorts bits 0-20 (representing x), the second sorts bits 21-41 (representing y),
 * and the third sorts bits 42-62 (representing z).
 *
 *
 * Before initiating the radix sorting passes, Morton keys are categorized into coarse bins
 * based on their high-order bits. By leveraging spatial coherence, this optimization ensures
 * nearby objects are likely to fall into the same bin, which, when sorted, will minimize cache misses.
 *
 *
 * @param keys     The Morton keys to be sorted.
 * @param numKeys  The number of keys.
 */
static inline void ThreePassRadixSortMortonKeys(spatialKey* keys, size_t* indexes, size_t numKeys) // Three-pass radix sort with binning optimization
{
	spatialKey* tempKeys = new spatialKey[numKeys];
	size_t* tempIndexes = new size_t[numKeys];  // Temporary storage for indexes
	
	// Histogram and offset tables for each bin (used in each pass)
	size_t* histograms = new size_t[NUM_BINS * NUM_BINS];
	size_t* offsets = new size_t[NUM_BINS * NUM_BINS];
	
	// Sorting based on bits representing x-coordinates (bits 0-20)
	for (int shift = 0; shift <= 20; shift += 8)
	{
		RadixSortPass(keys, indexes, numKeys, tempKeys, tempIndexes, shift, histograms, offsets);
	}
	
	// Sorting based on bits representing y-coordinates (bits 21-41)
	for (int shift = 21; shift <= 41; shift += 8)
	{
		RadixSortPass(keys, indexes, numKeys, tempKeys, tempIndexes, shift, histograms, offsets);
	}
	
	// Sorting based on bits representing z-coordinates (bits 42-62)
	for (int shift = 42; shift <= 62; shift += 8)
	{
		RadixSortPass(keys, indexes, numKeys, tempKeys, tempIndexes, shift, histograms, offsets);
	}
	
	delete[] tempKeys;
	delete[] tempIndexes;
	delete[] histograms;
	delete[] offsets;
}














/* ---------------- PERFORMANCE BENCHMARRKS ----------------     NOTE: time is recorded in seconds, benchmarks are tracked for each of the functions that sort morton keys
 
 
 ------------ BENCHMARKS FOR MORTON KEY ENCODING ------------
 --- Magic Bits method ---
 For 2097152 keys
 - Computing morton keys: 1.92861 seconds
 
 
 
 
 
 
 ------------ BENCHMARKS FOR MORTON KEY SORTING ------------
 --- MergeSortMortonKeys ---
 For 2097152 keys
 - Creating morton keys: 1.91779 seconds
 - Sorting MortonKeys Time: 3.13582 seconds
 - Creating morton keys: 1.80638 seconds
 - Sorting MortonKeys Time: 2.98662 seconds
 
 For 20971520 keys
 - Creating morton keys: 13.5163 seconds
 - Sorting MortonKeys Time: 26.4151 seconds
 
 
 
 
 --- RadixSortMortonKeys ---
 For 2097152 keys
 - Creating morton keys: 2.00174 seconds
 - Sorting MortonKeys Time: 2.17463 seconds
 - Creating morton keys: 1.81512 seconds
 - Sorting MortonKeys Time: 1.97575 seconds
 
 For 20971520 keys
 - Creating morton keys: 13.6232 seconds
 - Sorting MortonKeys Time: 15.3477 seconds
 
 
 
 
 --- ThreePassRadixSortMortonKeys ---
 For 2097152 keys
 - Creating morton keys: 1.93241 seconds
 - Sorting MortonKeys Time: 2.10793 seconds
 - Creating morton keys: 1.88512 seconds
 - Sorting MortonKeys Time: 2.04663 seconds
 
 For 20971520 keys
 - Creating morton keys: 13.653 seconds
 - Sorting MortonKeys Time: 15.35 seconds
 
 */
