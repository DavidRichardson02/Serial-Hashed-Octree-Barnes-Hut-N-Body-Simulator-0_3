//  MathUtils.cpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3


#include "MathUtils.hpp"
#include <math.h>




Vec3D::Vec3D() : x(0), y(0), z(0) {}

Vec3D::Vec3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) // using an initializer list to initialize members instead of assignment in the body of the constructor
{}

Vec3D::Vec3D(double magnitude)
{
	double theta = acos(3/sqrt(magnitude)*1);
	
	x = magnitude * cos(theta);
	y = magnitude * sin(theta);
	z = cos(theta) * sin(theta);
	
}


Vec3D& Vec3D::operator=(const Vec3D& other)
{
	if (this != &other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
	}
	return(*this);
}

bool Vec3D::operator!=(const Vec3D& other)
{
	return(x != other.x || y != other.y || z != other.z);
}


Vec3D Vec3D::operator+(const Vec3D& other)
{
	return(Vec3D(x + other.x, y + other.y, z + other.z));
}

Vec3D Vec3D::operator-(const Vec3D& other)
{
	return(Vec3D(x - other.x, y - other.y, z - other.z));
}

Vec3D Vec3D::operator*(const double scalar)
	{
	return(Vec3D(x * scalar, y * scalar, z * scalar));
}


void Vec3D::operator+=(const Vec3D& other)
	{
	x = x + other.x;
	y = y + other.y;
	z = z + other.z;
}

void Vec3D::operator-=(const Vec3D& other)
	{
	x = x - other.x;
	y = y - other.y;
	z = z - other.z;
}

void Vec3D::operator*=(const double scalar)
	{
	x = x * scalar;
	y = y * scalar;
	z = z * scalar;
}

void Vec3D::operator/=(const double scalar)
	{
	if (scalar != 0)
	{
		double inverseScalar = 1 / scalar;
		x = x * inverseScalar;
		y = y * inverseScalar;
		z = z * inverseScalar;
	}
}

Vec3D Vec3D::scaleVector(double scalar) const
	{
	Vec3D scaledVector = { x * scalar, y * scalar ,z * scalar };
	return(scaledVector);
}

double Vec3D::vectorLength() const
{
	
	return(sqrt(x * x + y * y + z * z));
}

double Vec3D::vectorSquareLength() const
{
	return(x * x + y * y + z * z);
}

Vec3D Vec3D::vectorNormalize() const
{
	double length = vectorLength();
	if (length != 0)
	{
		double inverseLength = 1 / length;
		Vec3D normalVec = { x * inverseLength , y * inverseLength , z * inverseLength };
		return(normalVec);
	}
	else
	{
		return(*this);
	}
}

double Vec3D::vectorDistance(const Vec3D& other)
{
	return(sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y) + (z - other.z) * (z - other.z)));
}


void Vec3D::reset()
{
	x = 0;
	y = 0;
	z = 0;
}







//
/*
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
 *
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
 *
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
 *
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
 *
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
 *
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


//*/
