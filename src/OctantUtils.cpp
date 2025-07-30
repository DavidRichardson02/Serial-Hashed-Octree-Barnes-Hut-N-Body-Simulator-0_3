//  OctantUtils.cpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3


#include "OctantUtils.hpp"




OctantBounds::OctantBounds() :center(0.0, 0.0, 0.0), size(0.0) {}
OctantBounds::OctantBounds(Vec3D _center, double _size) :center(_center), size(_size) {}
OctantBounds::OctantBounds(const OctantBounds& other) : size(other.size), center(other.center) {}


OctantBounds& OctantBounds::operator=(const OctantBounds& other)
{
	if (this != &other)
	{
		size = other.size;
		center = other.center;
	}
	return(*this);
}



OctantBounds::OctantBounds(Body* nBodies, size_t numBodies) : center (0.0,0.0,0.0 ), size(0.0)
{
	
	// Initializing the center of the global octant bounds at the average position of all bodies
	for (size_t i = 0; i < numBodies; i++)
	{
		//center += nBodies[i].position; // This was originally not commented out, unsure whether it should be or not.
		center.x += nBodies[i].position.x;
		center.y += nBodies[i].position.y;
		center.z += nBodies[i].position.z;
	}
	center /= numBodies;
	
	
	
	
	// Compute the size of the octant bounds(from the )
	// Because the center of this global boundary is located at the average position of all bodies in the simulation,
	// we assume the worst case and use the absolute position of each body relative to the center of the octant so
	// that the size of the octant bounds(from the center) will always accomodate the furthest position occupied by a body.
	Vec3D absPosition;
	double tempSize = size;
	for (size_t i = 0; i < numBodies; i++)
	{
		absPosition.x = fabs(nBodies[i].position.x);
		absPosition.y = fabs(nBodies[i].position.y);
		absPosition.z = fabs(nBodies[i].position.z);
		if ((absPosition.x - center.x) > tempSize) { tempSize = absPosition.x - center.x; }
		if ((absPosition.y - center.y) > tempSize) { tempSize = absPosition.y - center.y; }
		if ((absPosition.z - center.z) > tempSize) { tempSize = absPosition.z - center.z; }
	}
	size = tempSize;
}

OctantBounds::OctantBounds(Body* localBodies, size_t start, size_t end) : center (0.0,0.0,0.0 ), size(0.0)
{
	for (size_t i = start; i < end; i++)
	{
		center += localBodies[i].position;
	}
	center /= (end - start);
	
	Vec3D absPosition;
	double tempSize = size;
	for (size_t i = start; i < end; i++)
	{
		absPosition.x = (localBodies[i].position.x);
		absPosition.y = (localBodies[i].position.y);
		absPosition.z = (localBodies[i].position.z);
		if ((absPosition.x - center.x) > tempSize) { tempSize = absPosition.x - center.x; }
		if ((absPosition.y - center.y) > tempSize) { tempSize = absPosition.y - center.y; }
		if ((absPosition.z - center.z) > tempSize) { tempSize = absPosition.z - center.z; }
	}
	size = tempSize * 2;
	
}


