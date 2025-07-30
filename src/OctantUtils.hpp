//  OctantUtils.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3

// The 'OctantBounds' class's constructors are used primarily for constructing the octantbounds that will contain the entire tree as well as all other entities in the simulation

#pragma once
#include "MathUtils.hpp"
#include "SequenceContainers.hpp"
#include "MortonKeys.hpp"
#include "ObjectPool.hpp"
#include "Body.hpp"
#include "ofMain.h"




/* GLOBAL CONSTANTS & ENUMERATIONS */
static const int  DEFAULT_HASHED_OCTREE = 16; // Default capacity for the hashed octree, value of 16 fastest recorded(smalllll difference), followed closely by 14 and 18
static const size_t CHILD_OCTANTS = 8; // Number of child octants in an octree node
static const spatialKey ROOT_KEY = 1; // Root key for the tree


// OctantEnum: Defines the octant ordering in a 3D coordinate space based on Morton ordering (x > y > z).
enum OctantEnum
{
	Oct_0 = 0x0, //---
	Oct_1 = 0x1, //--+
	Oct_2 = 0x2, //-+-
	Oct_3 = 0x3, //-++
	Oct_4 = 0x4, //+--
	Oct_5 = 0x5, //+-+
	Oct_6 = 0x6, //++-
	Oct_7 = 0x7, //+++
};
// Lookup table for octant directions.
static const double OctantDir[] =
{
	-1.0, -1.0, -1.0,
	-1.0, -1.0,  1.0,
	-1.0,  1.0, -1.0,
	-1.0,  1.0,  1.0,
	1.0, -1.0, -1.0,
	1.0, -1.0,  1.0,
	1.0,  1.0, -1.0,
	1.0,  1.0,  1.0
};




// Helper Functions: These functions aid in manipulating and determining octant-specific properties.
static inline bool HasOctChild(uint8_t childByte, OctantEnum octant) // Checks if a specific octant is marked as having a child node. If non-zero, this octant has a child.
{
	return (childByte & (1 << octant)) != 0;
}
static inline void SetOctChild(uint8_t& childByte, OctantEnum octant) //set a specific bit in the childByte representing the existence of a child in the given octant, sets the corresponding bit to 1.
{
	childByte |= (1 << octant);
}
static inline void UnsetOctChild(uint8_t& childByte, OctantEnum octant)// unset a specific bit in the childByte, indicating the removal of a child in the given octant, sets the corresponding bit to 0
{
	childByte &= ~(1 << octant);
}
static inline OctantEnum DetermineOctant(const Vec3D& center, const Vec3D& bodyPosition) // Determines the octant for a body based on its position and the center of the current node.
{
	unsigned posX = (bodyPosition.x > center.x);
	unsigned posY = (bodyPosition.y > center.y);
	unsigned posZ = (bodyPosition.z > center.z);
	return static_cast<OctantEnum>((posX << 2) | (posY << 1) | posZ);
}


static inline Vec3D DetermineOctantCenter(const Vec3D& parentCenter, const double& parentSize, OctantEnum octant) // Retrieves the direction for a given octant based on the lookup table.
{
	double halfSize = parentSize * 0.5;
	Vec3D octantCenter;
	octantCenter.x = parentCenter.x + (halfSize * OctantDir[3 * octant]);
	octantCenter.y = parentCenter.y + (halfSize * OctantDir[3 * octant + 1]);
	octantCenter.z = parentCenter.z + (halfSize * OctantDir[3 * octant + 2]);
	return(octantCenter);
}


static inline Vec3D DetermineOctantDirection(OctantEnum octant) // Retrieves the direction for a given octant based on the lookup table.
{
	Vec3D octantDir;
	octantDir.x = OctantDir[3 * octant];
	octantDir.y = OctantDir[3 * octant + 1];
	octantDir.z = OctantDir[3 * octant + 2];
	return(octantDir);
}
static inline spatialKey GetParentKey(const spatialKey childKey) // shifts the childKey three places to the right to get the parent key
{
	const spatialKey parentKey = childKey >> 3;
	return(parentKey);
}
static inline spatialKey GetChildKey(const spatialKey parentKey, OctantEnum octant) //generates the key for a child node by shifting the parentKey three places to the left and adding the octant index
{
	const spatialKey childKey = ((parentKey) << 3 | (octant));
	return(childKey);
}
static inline OctantEnum GetOctantFromKey(spatialKey key) //Retrieves the octant enumeration from a key by taking the last three bits of the key. This is done by taking the key modulo 8, i.e., a bitwise-AND with 7, which gives the last three bits of the key.
{
	OctantEnum octant = static_cast<OctantEnum>(key & static_cast <spatialKey>(7));
	return(octant);
}





class OctantBounds
{
public:
	OctantBounds();
	OctantBounds(Body* nBodies, size_t numBodies); //returns a BoundingBox which has its center at the average position and encloses all bodies
	OctantBounds(Body* localBodies, size_t start, size_t end); //returns a BoundingBox which has its center at the average position and encloses all bodies
	OctantBounds(Vec3D _center, double _size);
	OctantBounds(const OctantBounds& other);
	OctantBounds& operator=(const OctantBounds& other);
	Vec3D center;
	double size;
};



