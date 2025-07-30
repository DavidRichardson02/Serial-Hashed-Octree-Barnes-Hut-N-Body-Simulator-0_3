//  MathUtils.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3

/**
 * Vec3D Class: provides an abstraction for three-dimensional vectors, encapsulating the vector's components (x, y, z)
 *
 * This class encapsulates a 3D vector's components (x, y, z),  and provides methods
 * to perform basic vector operations and manipulate vector data.
 */


#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "ofMain.h"


// Vec3D class representing a 3D Vector
class Vec3D
{
public:
	// ------------- Constructors -------------
	Vec3D();
	Vec3D(double _x, double _y, double _z);
	
	Vec3D(double magnitude);
	
	
	
	
	// ------------- Comparison operators -------------
	bool operator!=(const Vec3D& other);
	
	
	
	
	// ------------- Accessors (Return new Vec3D based on current one) -------------
	Vec3D operator+(const Vec3D& other);
	Vec3D operator-(const Vec3D& other);
	Vec3D operator*(const double scalar);
	
	
	
	
	// ------------- Modifiers (Modify the current Vec3D and return reference) -------------
	Vec3D& operator=(const Vec3D& other);
	void operator+=(const Vec3D& other);
	void operator-=(const Vec3D& other);
	void operator*=(const double scalar);
	void operator/=(const double scalar);
	
	
	
	
	// ------------- Vector operations (do not modify the internal state) -------------
	Vec3D scaleVector(double scalar) const;
	double vectorLength() const;
	double vectorSquareLength() const;
	Vec3D vectorNormalize() const;
	double vectorDistance(const Vec3D& other);
	void reset();
	
	
	// ------------- Member variables -------------
	//private:
	double x;
	double y;
	double z;
};


static inline void ResetVec3D(Vec3D *&vec3d, int numVecs);
static inline void ParameterizeVec3D(Vec3D *&vec3d, int numVecs);
static inline void PrintVec3D(const Vec3D &vec3d);


static inline void ResetVec3D(Vec3D *&vec3d, int numVecs)
{
	if(!vec3d)
	{
		throw std::invalid_argument("vec3d is a null pointer in 'ResetVec3D'.");
	}
	
	for(int i = 0; i < numVecs; i++)
	{
		vec3d[i].reset();
		//vec3d[i].x = 0.0;
		//vec3d[i].y = 0.0;
		//vec3d[i].z = 0.0;
	}
}

static inline void ParameterizeVec3D(Vec3D *&vec3d, int numVecs)
{
	if(vec3d == nullptr)
	{
		throw std::invalid_argument("vec3d is a null pointer in 'ResetVec3D'.");
	}
	
	for(int i = 0; i < numVecs; i++)
	{
		vec3d[i].reset();
		//vec3d[i].x = 0.0;
		//vec3d[i].y = 0.0;
		//vec3d[i].z = 0.0;
	}
}


static inline void PrintVec3D(const Vec3D &vec3d)
{
	std::cout << "Vec3D: (" << vec3d.x << ", " << vec3d.y << ", " << vec3d.z << ")" << std::endl;
}


//
/*
// ------------- Helper Functions for Bitwise Operations(on Binary Numerals?) -------------
/// \{
uint64_t flip_sign_bit(uint64_t value); // Helper function to flip the sign bit of the double's binary representation.
uint64_t double_to_uint64(double value); // Helper function to reinterpret a double as an uint64_t.
double uint64_to_double(uint64_t value); // Helper function to reinterpret a uint64_t as an double.
uint64_t double_to_mapped_uint64(double value); // Helper function to map a double to an uint64_t.
double mapped_uint64_to_double(uint64_t u); // Helper function to map a uint64_t to a double.
/// \}
//*/
