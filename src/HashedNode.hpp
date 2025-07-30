//  HashedNode.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3
/*
 * =====================================================================================
 * Hashed Octree Node
 * =====================================================================================
 
 * OVERVIEW:
 *
 *
 *
 *
 * CORE COMPONENTS:
 *
 * 1. STRUCTURE:
 *    - The tree follows Morton Ordering: This implementation adopts Morton ordering for its nodes, wherein octants are interpreted based on a coordinate system with ordering: x > y > z.
 *	  - The root node has index 1.
 *    - Each child node's index is the concatenation of its parent index with the octant direction, encoded over two bits.
 *
 *
 * 2. GLOBAL CONSTANTS & ENUMERATIONS:
 *    - DEFAULT_HASHED_OCTREE: Default capacity for the hashed octree.
 *    - CHILD_OCTANTS: Defines the number of child octants.
 *    - ROOT_KEY: Root key for the tree.
 *
 *
 * 3. OCTANTS:
 *    - Eight octants in 3D space, defined by OctantEnum, an enumeration representing octant ordering in 3D space based on Morton ordering.
 *    - A lookup table, OctantDir[], provides the direction of each octant.
 *
 *
 * 4. HELPER FUNCTIONS: A collection of inline functions facilitate operations on octants:
 *    - Child octant determination and manipulation.
 *    - Morton key manipulations for parent-child relationships, i.e., determining the parent or child key, setting, and unsetting children
 *
 *
 * 5. HOTNode Class: Represents an individual node in the hashed octree. Key features include:
 *    - Constructors for different use cases (e.g., node creation from parent nodes or Morton keys).
 *    - Getter and Setter methods for node attributes.
 *    - Functions for managing parent and child relationships.
 *    - Defines the core node attributes: spatial key, boundaries, barycenter, mass, number of bodies, children.
 *
 *
 */


#pragma once
#include "MathUtils.hpp"
#include "SequenceContainers.hpp"
#include "MortonKeys.hpp"
#include "ObjectPool.hpp"
#include "Body.hpp"
#include "OctantUtils.hpp"
#include "ofMain.h"




// HOTNode Class: Represents an individual node in the hashed octree structure.
class HOTNode
{
public:
	// Constructors and Destructor
	HOTNode();
	HOTNode(const HOTNode& other);
	HOTNode& operator=(const HOTNode& other);
	HOTNode(HOTNode* parentNode, OctantEnum targetOctant, const Vec3D& _baryCenter, const double& _mass);
	HOTNode(const Vec3D& _position, const double _size);
	HOTNode(const OctantBounds _nodeBounds);
	HOTNode(const OctantBounds _nodeBounds, const spatialKey rootKey);
	HOTNode(const double _size, const spatialKey rootKey);
	~HOTNode();
	
	
	bool containsBody(const Vec3D& bodyPosition);
	
	
	void reset();
	void updateCenterOfMass();
	void initializeNode(const OctantBounds _nodeBounds, const spatialKey rootKey);
	void insertBodyDirectly(spatialKey parentKey, OctantEnum targetOctant, const Vec3D bodyPosition, const double bodyMass);
	void parameterizeChildNode(HOTNode* parentNode, const OctantEnum targetOctant, const Vec3D& _baryCenter, const double& _mass);
	
	//Upper triangle of quadrupole moment tensor Q ordered as
	//Qxx, Qxy, Qxz, Qyy, Qyz, Qzz
	double quadrupoleMoment[6];
	
	Vec3D baryCenter;      double mass;// the center of mass and total mass of all bodies at or below this node
	long N; //number of bodies at or below this node.
	uint8_t childByte; //a bitfield encoding which children actually exist, each of the 8 bits can represent the existence of one of the 8 children in the octree (where a set bit indicates that the child exists and an unset bit indicates the opposite).
					   //private:
					   //OctantEnum octant; //compute the octant on the fly
	OctantBounds nodeBounds;
	spatialKey nodeKey; //the spatial key for this node/body
};










static inline void ResetHOTNodeList(HOTNode **&hotNodeList, int numBodies);





