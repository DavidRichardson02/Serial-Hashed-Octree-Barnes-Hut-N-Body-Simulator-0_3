//  HashedOctree.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3
/*
 * =====================================================================================
 * Hashed Octree Implementation Description
 * =====================================================================================
 *
 * OVERVIEW:
 * This file contains an implementation of a hashed octree data structure.
 * An octree is a tree data structure in which each internal node has eight children, typically used to partition a three-dimensional space.
 * A hashed octree stores its nodes in a hash map, indexed by their Morton codes, the main aim of which is to facilitate efficient spatial operations,
 * and this implementation particularly optimizes lookup operations by hashing nodes based on Morton ordering, where the precedence is x > y > z.
 * Morton codes, also known as Z-order codes, provide a unique identifier for each node in the octree based on its spatial location.
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
 * 6. HASHING:
 *    - A std::unordered_map is used to store nodes in the octree with Morton keys as indices, allowing for quick access and modification.
 *
 *
 * 7. HashedOctree Class: Manages nad represents the entire hashed octree structure. Key features include:
 *    - Constructors for creating trees from different initial conditions.
 *	  - Functions for traversing the tree, getting the depth, and accessing parent/child nodes.
 *	  - Methods for inserting nodes and bodies into the octree.
 *    - A hash map (std::unordered_map) storing octree nodes using Morton code as the key.
 *
 *
 */


#pragma once
#include "MathUtils.hpp"
#include "SequenceContainers.hpp"
#include "MortonKeys.hpp"
#include "ObjectPool.hpp"
#include "Body.hpp"
#include "HashedNode.hpp"
#include "OctantUtils.hpp"
#include "ofMain.h"





class LinearHashedOctree
{
public:
	LinearHashedOctree();
	LinearHashedOctree(const OctantBounds& rootBounds); //create octree based on rootBounds
	~LinearHashedOctree();
	
	
	
	// Octree Operations: Functions to interact with the tree
	size_t GetNodeTreeDepth(const spatialKey nodeKey); // Function to get the depth of the tree
	HOTNode* getParentNode(const HOTNode* childNode); // Get the parentNode of a node
	HOTNode* getChildNode(const HOTNode* parentNode); // Get the childNode of a node
	HOTNode* lookUpNode(const spatialKey code);// Lookup a node by its Morton key
	void traverseTree(HOTNode* node);
	bool descendantOf(spatialKey nodeKey, spatialKey parentKey); //checks if a key is contained within a parent key
	bool leafNodeExistsForBody(const Vec3D& bodyPosition);
	
	std::vector<OctantEnum> getPathToLeafNode(const Vec3D& bodyPosition);
	
	
	
	// Insertion Functions: Add nodes or bodies to the octree
	void insertHOTNode(HOTNode*& node);  // Insert a node into the Octree with a specified Morton code. If the node already exists, overwrite it.
	
	void insertBody(HOTNode*& node, const Vec3D bodyPosition, const double bodyMass);
	void pushHOTNodeBodyToChild(HOTNode* node);
	void processNodeInsertion(HOTNode*& node, OctantEnum& targetOctant, spatialKey& currentKey, const Vec3D bodyPosition, const double bodyMass);
	void pruneEmptyNodes(HOTNode* node);
	void deleteNode(spatialKey nodeCode);
	void clear();
	void deleteTree();
	
	void createChildNode(HOTNode*& node, OctantEnum& targetOctant, const Vec3D& bodyPosition, const double& bodyMass);
	
	void computeNodeBaryCenters(HOTNode*& node);
	void computeTreeBaryCenters(HOTNode*& node);
	
	
	bool MAC(HOTNode* node, const Vec3D& r, double theta);
	
	
	// Debug utility to print the HashedOctree's details
	void printHashedOctree();
	
	// Rendering utility to visualize the HashedOctree's node's octant bounds
	void visualizeTree();
	void visualizeNode(const HOTNode* node);
	
	
	// Traverse down the tree until we reach the lowest level or a leaf node.
	// Your hashing function should determine the path to the leaf node.
	
	
	//private:
	
	// Hash map to store the Octree nodes with Morton code as key
	std::unordered_map<spatialKey, HOTNode*> nodes;
};




static inline int BarnesHutHOTMAC(HOTNode* node, Vec3D bodyPosition, double theta);
static inline HOTNode* LookUpNode(LinearHashedOctree& HTree, spatialKey code);// Lookup a node by its Morton key

static inline void testBuildLinearHashedOctreeInPlace(LinearHashedOctree& HTree, Body*& bodies, const size_t& numBodies, OctantBounds& domainBounds);
static inline void buildLinearHashedOctreeInPlace(LinearHashedOctree& HTree, Body*& bodies, const size_t& numBodies, OctantBounds& domainBounds);
//static inline void buildHashedOctreePool(LinearHashedOctree &HTree, ObjectPool<HOTNode> nodePool, Body* bodies, const size_t numBodies, OctantBounds domainBounds);
static inline void PruneEmptyNodesFromTree(LinearHashedOctree& HTree);
static inline void ComputeHOTOctreeBaryCenters(LinearHashedOctree& HTree, HOTNode* rootNode);
static inline void ComputeHOTOctreeForce(LinearHashedOctree& HTree, Body*& bodies, Vec3D*& bodiesAccelerations, const size_t& numBodies, HOTNode**& walkList, HOTNode**& interactList, double thetaMAC);
static inline long TraverseHOTInteractionList(LinearHashedOctree& LHTree, Vec3D& bodyPosition, double theta, HOTNode**& walkList, HOTNode**& interactList);
static inline void ComputeHOTForceInteractionList(Vec3D& bodyPosition, double& bodyMass, Vec3D& acceleration, HOTNode**& interactList, long listLength);
const double SOFTENING = 0.025;


inline int BarnesHutHOTMAC(HOTNode* node, Vec3D bodyPosition, double theta)
{
	// node's diamater / distance-to-bodyPosition < theta
	Vec3D distance = node->baryCenter - bodyPosition;
	
	double distSquared = distance.x * distance.x + distance.y * distance.y + distance.z * distance.z;
	
	return(4.0 * node->nodeBounds.size * node->nodeBounds.size < distSquared * theta * theta);
}

inline HOTNode* LookUpNode(LinearHashedOctree& HTree, const spatialKey code)
{
	const auto iter = HTree.nodes.find(code);
	return (iter == HTree.nodes.end() ? nullptr : iter->second);
}



static inline void testBuildLinearHashedOctreeInPlace(LinearHashedOctree& HTree, Body*& bodies, const size_t& numBodies, OctantBounds& domainBounds)
{
	
	HTree.clear();
	//HTree = LinearHashedOctree(domainBounds);
	//HOTNode* rootNode = HTree.lookUpNode(ROOT_KEY);//LookUpNode(HTree, ROOT_KEY);//new HOTNode(domainBounds, ROOT_KEY);
	
	
	
	HOTNode* rootNode = new HOTNode(domainBounds, ROOT_KEY);
	//HTree.nodes[rootNode->nodeKey] = rootNode;
	HTree.insertHOTNode(rootNode);
	for (size_t i = 0; i < numBodies; i++)
	{
		HTree.insertBody(rootNode, bodies[i].position, bodies[i].mass);
	}
	
	cout<<"\n\nTree Constructed ofGetElapsedTimef(): " << ofGetElapsedTimef();
	
	
	PruneEmptyNodesFromTree(HTree);
	if (numBodies > 0)
	{
		ComputeHOTOctreeBaryCenters(HTree, rootNode);
	}
	
	
	cout<<"\nTree Pruned and Barycenters Computed  ofGetElapsedTimef(): " << ofGetElapsedTimef();
	
	//HTree.visualizeTree();
	
	HTree.printHashedOctree();
	
	//HTree.deleteTree();
	delete rootNode;
	
}


static inline void buildLinearHashedOctreeInPlace(LinearHashedOctree& HTree, Body*& bodies, const size_t& numBodies, OctantBounds& domainBounds)
{
	
	HTree.clear();
	//HTree = LinearHashedOctree(domainBounds);
	//HOTNode* rootNode = HTree.lookUpNode(ROOT_KEY);//LookUpNode(HTree, ROOT_KEY);//new HOTNode(domainBounds, ROOT_KEY);
	
	
	
	HOTNode* rootNode = new HOTNode(domainBounds, ROOT_KEY);
	//HTree.nodes[rootNode->nodeKey] = rootNode;
	HTree.insertHOTNode(rootNode);
	for (size_t i = 0; i < numBodies; i++)
	{
		HTree.insertBody(rootNode, bodies[i].position, bodies[i].mass);
	}
	
	
	
	PruneEmptyNodesFromTree(HTree);
	if (numBodies > 0)
	{
		ComputeHOTOctreeBaryCenters(HTree, rootNode);
	}
	
	
	//HTree.visualizeTree();
	
	//HTree.printHashedOctree();
	
	//HTree.deleteTree();
	delete rootNode;
	
}



inline void PruneEmptyNodesFromTree(LinearHashedOctree& HTree)
{
	
	HTree.pruneEmptyNodes(HTree.lookUpNode(ROOT_KEY));
}

inline void ComputeHOTOctreeBaryCenters(LinearHashedOctree& HTree, HOTNode* rootNode)
{
	HTree.computeTreeBaryCenters(rootNode);
}




inline void ComputeHOTOctreeForce(LinearHashedOctree& LHTree, Body*& bodies, Vec3D*& bodiesAccelerations, const size_t& numBodies, HOTNode**& walkList, HOTNode**& interactList, double thetaMAC)
{
	
	long interactionListLength = 0;
	for (size_t i = 0; i < numBodies; i++)
	{
		
		interactionListLength = TraverseHOTInteractionList(LHTree, bodies[i].position, thetaMAC, walkList, interactList);
		ComputeHOTForceInteractionList(bodies[i].position, bodies[i].mass, bodiesAccelerations[i], interactList, interactionListLength);
		
	}
	
	//
	/*
	 
	 omp_set_num_threads(NUM_THREADS);
	 
	 #pragma omp parallel
	 {
	 int id = omp_get_thread_num(); // current thread
	 int numThreads = omp_get_num_threads(); // Get the number of threads inside parallel section
	 
	 
	 // Calculate start and end indices based on id and numThreads
	 int bodiesPerThread = numBodies / numThreads;
	 int remainderBodies = numBodies % numThreads;//(int)numBodies - (numThreads * bodiesPerThread);
	 
	 // Thread-local copies of walkList and interactList.
	 
	 // Thread assignment: Distribute extra bodies to the first 'remainderBodies' threads.
	 size_t start = id * bodiesPerThread;
	 size_t end = start + bodiesPerThread;
	 if (id < remainderBodies)
	 {
	 start += id;  // One extra body for this thread.
	 end += (id + 1);  // End index will be one more.
	 }
	 else
	 {
	 start += remainderBodies; // All extra bodies distributed before this thread.
	 end += remainderBodies; // Same here.
	 }
	 
	 // Thread-local copies of walkList and interactList.
	 HOTNode** localWalkList = new HOTNode * [end - start];  // Set your maximum size
	 HOTNode** localInteractList = new HOTNode * [end - start];  // Set your maximum size
	 
	 #pragma omp for
	 for (int i = start; i < end; i++)
	 {
	 
	 interactionListLength = TraverseHOTInteractionList(LHTree, bodies[i].position, thetaMAC, localWalkList, localInteractList);//walkList, interactList);
	 ComputeHOTForceInteractionList(bodies[i].position, bodies[i].mass, bodiesAccelerations[i], localInteractList, interactionListLength);
	 
	 }
	 // Free thread-local storage
	 delete[] localWalkList;
	 delete[] localInteractList;
	 }
	 
	 
	 //*/
}






/*
 intended to generate the interaction list for a given bodyPosition
 The list will include either leaf nodes (individual bodies) or internal nodes (groups of bodies)
 that satisfy the MAC (multipole acceptance criterion) for approximation.
 */
static inline long TraverseHOTInteractionList(LinearHashedOctree& LHTree, Vec3D& bodyPosition, double theta, HOTNode**& walkList, HOTNode**& interactList)
{
	if (LHTree.nodes.empty())
	{
		return 0;
	}
	
	//if my parallel call of TraverseHOTInteractionList tries to call LHTree in a wrong order, the root node has not been initialized
	//if (LHTree.lookUpNode(ROOT_KEY) == nullptr)
	//{
	//    cout << "\n\nUH-OH, gay \n\n";
	//}
	
	HOTNode* node;
	HOTNode* child;
	walkList[0] = LHTree.lookUpNode(ROOT_KEY);
	long walkIdx = 1;
	long intListIdx = 0;
	while (walkIdx > 0)
	{
		node = walkList[--walkIdx];
		if (node->childByte != 0)
		{
			for (int i = 0; i < 8; i++)
			{
				if (!HasOctChild(node->childByte, static_cast<OctantEnum>(i))) //    if (!node->childByte & (1 << i))// See if the ith child exist
				{
					continue;
				}
				child = LHTree.lookUpNode(GetChildKey(node->nodeKey, static_cast<OctantEnum>(i)));
				if (LHTree.MAC(child, bodyPosition, theta))
				{
					interactList[intListIdx++] = child;
				}
				else
				{
					walkList[walkIdx++] = child;
				}
			}
		}
		else if (node->baryCenter != bodyPosition)
		{
			interactList[intListIdx++] = node;
		}
	}
	
	//delete node;
	//delete child;
	return intListIdx;
}







inline void ComputeHOTForceInteractionList(Vec3D& bodyPosition, double& bodyMass, Vec3D& acceleration, HOTNode**& interactList, long listLength)
{
	//acceleration = {0,0,0};
	
	
	HOTNode* node;
	//node = new HOTNode();
	double dx = 0, dy = 0, dz = 0, D1 = 0, D2 = 0;
	double qx = 0, qy = 0, qz = 0;
	
	
	for (int i = 0; i < listLength; i++)
	{
		//node->reset();
		node = interactList[i];
		
		
		//convert r to c.o.m. referece frame of the node.
		dx = node->baryCenter.x - bodyPosition.x;
		dy = node->baryCenter.y - bodyPosition.y;
		dz = node->baryCenter.z - bodyPosition.z;
		
		
		// m*r / |r|^3 : normal monopole part/direct interaction
		D2 = dx * dx + dy * dy + dz * dz;
		D2 += SOFTENING * SOFTENING;
		
		D1 = 1.0 / sqrt(D2); //hopefully compiler catches the one over sqrt
		D2 = 1.0 / D2;
		D1 = D1 * D2; // 1/D3
		
		//m*r / |r|^3
		acceleration.x += node->mass * dx * D1;
		acceleration.y += node->mass * dy * D1;
		acceleration.z += node->mass * dz * D1;
		
		
		if (node->N > 1)//just did monopole so now quadrupole approximate
		{
			
			//Q.r / |r|^5; recall quadMom is only the upper triangle of symmetric tensor
			qx = node->quadrupoleMoment[0] * dx + node->quadrupoleMoment[1] * dy + node->quadrupoleMoment[2] * dz;
			qy = node->quadrupoleMoment[1] * dx + node->quadrupoleMoment[3] * dy + node->quadrupoleMoment[4] * dz;
			qz = node->quadrupoleMoment[2] * dx + node->quadrupoleMoment[4] * dy + node->quadrupoleMoment[5] * dz;
			
			D1 *= D2; // 1/D5 now
			acceleration.x -= qx * D1;
			acceleration.y -= qy * D1;
			acceleration.z -= qz * D1;
			
			
			
			//5*r.Q.r*r / 2*|r|^7
			qx = dx * qx + dy * qy + dz * qz;
			qx *= 2.5;
			D1 *= D2; // 1/D7 now
			
			
			acceleration.x += qx * dx * D1;
			acceleration.y += qx * dy * D1;
			acceleration.z += qx * dz * D1;
		}
	}
	
	//delete node;
	
}



















/*
 
 bool HashedOctree::findLeafNodeAndPath(const Body& body, HOTNode*& leafNode, std::vector<OctantEnum>& path) {
 // Ensure the root node is initialized.
 assert(!nodes.empty());
 
 // Start from the root.
 spatialKey currentKey = ROOT_KEY;
 HOTNode* node = lookUpNode(currentKey);
 
 // Initialize path.
 path.clear();
 
 // Traverse the tree to find the appropriate leaf node for the body.
 while (node != nullptr) {
 // Determine the spatial subdivision for the body.
 OctantEnum targetOctant = DetermineOctant(node->nodeBounds.center, body.position);
 
 // Add to path.
 path.push_back(targetOctant);
 
 // Check if there's an existing child node in the target octant.
 if (HasOctChild(node->childByte, targetOctant)) {
 // Update the current key to the child's.
 spatialKey childKey = GetChildKey(currentKey, targetOctant);
 currentKey = childKey;
 
 // Fetch the child node from the HOT (Hashed Octree).
 node = lookUpNode(currentKey);
 }
 // If no child exists, then it's the leaf node we're looking for.
 else {
 leafNode = node;
 return true;
 }
 }
 
 // If we've reached here, no suitable leaf node was found.
 leafNode = nullptr;
 return false;
 }
 
 
 
 
 
 
 
 
 
 std::vector<OctantEnum> HashedOctree::getPathToLeafNode(const Vec3D& bodyPosition)
 {
 std::vector<OctantEnum> path;
 
 // Start with the root node
 HOTNode* currentNode = lookUpNode(ROOT_KEY);
 
 // Traverse the tree until we reach a leaf node or an empty node
 while(currentNode && currentNode->N != 0)
 {
 OctantEnum octant = DetermineOctant(currentNode->nodeBounds.center, bodyPosition);
 path.push_back(octant);
 
 spatialKey childKey = GetChildKey(currentNode->nodeKey, octant);
 
 if (nodes.find(childKey) != nodes.end())
 {
 currentNode = lookUpNode(childKey);
 }
 else
 {
 break;  // Leaf node or empty node found
 }
 }
 
 return path;
 }
 
 */



//
/*


static inline void DistanceInTermsofNodes_Test_01(HOTNode* node1, HOTNode* node2);
static inline void DistanceInTermsofNodes_Control(HOTNode* node1, HOTNode* node2);







static inline void DistanceInTermsofNodes_Control(HOTNode* node1, HOTNode* node2)
{
	double nodeDistance = 
}
//*/
