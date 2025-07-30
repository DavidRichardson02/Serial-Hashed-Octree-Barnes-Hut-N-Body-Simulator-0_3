//  HashedOctree.cpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3
// ----------------------- Overview -----------------------
/**
 * HashedOctree is a spatial partitioning data structure where space is divided into octants.
 * It uses a hashmap to store nodes, ensuring efficient memory usage and constant-time lookups.
 * Each node in the tree has a unique Morton key which preserves spatial locality.
 */






#include "HashedOctree.hpp"



/**
 * Default constructor for HashedOctree class.
 *
 * Initializes an empty HashedOctree by reserving space for the default number of nodes (2^DEFAULT_HASHED_OCTREE - 1).
 */
LinearHashedOctree::LinearHashedOctree()
{
	nodes.reserve((1 << DEFAULT_HASHED_OCTREE) - 1); //rehash
}


/**
 * Constructor with provided OctantBounds.
 *
 * This constructor initializes the tree with a HOTNode that encompasses the provided bounds and assigns it the ROOT_KEY.
 * The hash map is also prepared with a reservation for the default number of nodes.
 *
 * @param rootBounds: Bounding box that defines the space covered by the root node, i.e., the domain size of the hashed octree
 */
LinearHashedOctree::LinearHashedOctree(const OctantBounds& rootBounds)
{
	// Insert the node into the hashmap.
	//nodes.rehash((1 << DEFAULT_HASHED_OCTREE) - 1);
	//HOTNode* rootNode = new HOTNode(rootBounds, ROOT_KEY);
	//nodes[rootNode->nodeKey] = rootNode;
	//nodes.emplace(rootNode->nodeKey, rootNode);
	
}


LinearHashedOctree::~LinearHashedOctree()
{
	deleteTree();
	//clear();
}



/**
 * Compute the tree depth of a node based on its Morton key.
 * By bit-shifting the node key, each shift by 3 bits (representing one level in the octree) increments the depth counter.
 *
 * @param nodeKey: Morton key of the node.
 * @return The depth of the node in the tree.
 */
size_t LinearHashedOctree::GetNodeTreeDepth(const spatialKey nodeKey)
{
	assert(nodeKey);
	size_t depth = 0;
	for (spatialKey locationKey = nodeKey; locationKey != 1; locationKey >>= 3) //shift locationKey right by 3 bits, i.e., moving up one level in the tree until the root node is reached
	{
		depth++;  //increment depth by 1, i.e., moving down one level in the tree
	}
	return depth;
}


/**
 * Retrieve the parent node of a given node.
 * The parent's Morton key can be derived by right-shifting the child's key by 3.
 *
 * @param childNode: Pointer to the child node.
 * @return Pointer to the parent node if it exists, otherwise nullptr.
 */
HOTNode* LinearHashedOctree::getParentNode(const HOTNode* childNode)
{
	const spatialKey parentNode = childNode->nodeKey >> 3;
	return lookUpNode(parentNode);
}

/**
 * Retrieve the child node of a given node.
 * The child's Morton key can be derived by left-shifting the parent's key by 3.
 *
 * @param parentNode: Pointer to the parent node.
 * @return Pointer to the child node if it exists, otherwise nullptr.
 */
HOTNode* LinearHashedOctree::getChildNode(const HOTNode* parentNode)
{
	const spatialKey childNode = parentNode->nodeKey << 3;
	return lookUpNode(childNode);
}


/**
 * Lookup a node in the hashmap using its Morton key.
 * This method checks the existence of a node with the provided Morton key and if it exists, it returns the node, otherwise, it returns nullptr.
 *
 * @param code: Morton key of the node to be retrieved.
 * @return Pointer to the node if it exists, otherwise nullptr.
 */

HOTNode* LinearHashedOctree::lookUpNode(const spatialKey code)
{
	const auto iter = nodes.find(code);
	return (iter == nodes.end() ? nullptr : iter->second);
}

/**
 * Recursively traverse the octree, visiting all child nodes.
 * This traversal goes through all eight possible child nodes of a given node.
 *
 * @param node: Current node being traversed.
 */
void LinearHashedOctree::traverseTree(HOTNode* node)
{
	for (int i = 0; i < 8; i++) //For all eight possible children
	{
		if (node->childByte & (1 << i))// See if the ith child exist
		{
			const spatialKey locCodeChild = (node->nodeKey << 3) | i; //Compute new Morton key for the child
			auto* child = lookUpNode(locCodeChild);  //Using key, look child up in hash table and recursively visit subtree
			traverseTree(child);
		}
	}
}



/**
 * Inserts a HOTNode into the hashed octree's underlying hashmap.
 *
 * The method takes in a HOTNode and integrates it into the hashmap based on its Morton key.
 * The Morton key serves as a unique identifier, allowing for constant-time lookups in the hashmap.
 * Nodes with a key value of 1 (ROOT_KEY) are considered as root nodes.
 * Using hashmap ensures efficient space usage as only the populated parts of the octree are stored,
 * contrary to a dense tree where all possible nodes would need storage.
 *
 * Note: If the hashmap is empty, the given node is directly inserted as the root.
 * For uninitialized nodes (nullptr), an error message is displayed.
 *
 * @param node: The HOTNode to be inserted into the hashmap.
 */
void LinearHashedOctree::insertHOTNode(HOTNode*& node)
{
	///*
	if (nodes.empty()) //if this is the root node
	{
		node->nodeKey = ROOT_KEY;
		nodes[node->nodeKey] = node;
	}
	else //NOT the root node
	{
		if (nodes[node->nodeKey]) //this node already exists
		{
			//deleteNode(node->nodeKey);//delete old node
			
			nodes.erase(node->nodeKey);   // Remove the node from the map
			delete node;  // Release the node's memory
			
		}
		nodes[node->nodeKey] = node; //add node to hashmap
	}
	//*/
	
	//
	/*
	 if (nodes[node->nodeKey]) //this node already exists
	 {
	 nodes.erase(node->nodeKey);   // Remove the node from the map
	 delete node;  // Release the node's memory
	 
	 }
	 nodes[node->nodeKey] = node; //add node to hashmap
	 //*/
	
	
	
	//
	/*
	 auto result = nodes.emplace(node->nodeKey, node);
	 
	 if (!result.second) // this node already exists
	 {
	 // Overwrite the old node with the new one.
	 // Assuming ownership is transferred here, and it's safe to delete the old node.
	 delete result.first->second;  // delete old node
	 result.first->second = node;  // replace old node with new one
	 }
	 // Otherwise, the emplace has already added the new node for us
	 
	 //*/
	
	
	
	
}



/**
 * Determines whether a node, identified by its spatial key, is a descendant of another
 * node within the octree. This function is essential for verifying hierarchical
 * relationships between nodes and is used to confirm if one node exists beneath
 * another in the tree structure.
 *
 * @param nodeKey: The key of the node to be checked.
 * @param parentKey: The key of the potential parent node.
 * @return True if nodeKey is a descendant of parentKey, false otherwise.
 */
bool LinearHashedOctree::descendantOf(spatialKey nodeKey, spatialKey parentKey)
{
	if (parentKey == 1)  //only the root node has a key of 1, all nodes are descendants of root node
	{
		return(true);
	}
	if (nodeKey <= parentKey) //if this key is less than the parentKey, then it is an ancestor of the parent key, not a descendant (i.e., it resides higher up in the tree).
	{
		return(false);
	}
	// Calculate the tree depths of nodeKey and parentKey.
	size_t nodeKeyDepth = GetNodeTreeDepth(nodeKey);
	size_t parentKeyDepth = GetNodeTreeDepth(parentKey);
	
	
	// Shift the nodeKey to the same level as the parentKey.
	nodeKey >>= (nodeKeyDepth - parentKeyDepth) * 3;
	return (nodeKey == parentKey);
}




/**
 * Checks whether there exists a leaf node in the octree that directly corresponds to
 * a given body's position. This function traverses the tree from the root, following
 * the path dictated by the body's position, and determines if a leaf node (a node
 * with no children) exists at that position.
 */
bool LinearHashedOctree::leafNodeExistsForBody(const Vec3D& bodyPosition)
{
	// Initialize the traversal from the root node.
	spatialKey currentKey = ROOT_KEY;
	HOTNode* currentNode = lookUpNode(currentKey);
	
	// Traverse the tree until a leaf node or a null node is encountered.
	while (currentNode != nullptr && currentNode->N > 0)
	{
		OctantEnum targetOctant = DetermineOctant(currentNode->nodeBounds.center, bodyPosition); // Determine the octant corresponding to the body's position.
		currentKey = GetChildKey(currentKey, targetOctant); // Update the key to point to the child node in the target octant.
		currentNode = lookUpNode(currentKey); // Look up the child node using the updated key.
	}
	
	return currentNode != nullptr && currentNode->N == 0; // true if it's a leaf node
}



/**
 * This function computes and returns the path from the root of the octree to the
 * leaf node corresponding to a given body's position. It identifies the sequence
 * of octants traversed to reach a node that would either be a leaf or an empty node
 * where the body could be placed.
 */
std::vector<OctantEnum> LinearHashedOctree::getPathToLeafNode(const Vec3D& bodyPosition)
{
	// Initialize an empty vector to store the path of octants.
	std::vector<OctantEnum> path;
	
	// Start with the root node
	spatialKey currentKey = ROOT_KEY;
	HOTNode* currentNode = lookUpNode(ROOT_KEY);
	
	
	// Traverse the tree until we reach a leaf node or an empty node
	while (currentNode) //&& currentNode->N != 0
	{
		OctantEnum octant = DetermineOctant(currentNode->nodeBounds.center, bodyPosition);
		path.push_back(octant);
		
		// Calculate the spatial key for the child node in the identified octant.
		currentKey = GetChildKey(currentNode->nodeKey, octant);
		
		// Check if the child node exists in the nodes map.
		if (nodes.find(currentKey) != nodes.end())
		{
			// Update the current node to the child node.
			currentNode = lookUpNode(currentKey);
		}
		else
		{
			break;  // Exit the loop if a leaf node or an empty node is found.
		}
	}
	
	
	
	// Print the path in octant representation.
	std::cout << "\n\nPath to leaf node: \n";
	for(auto& octant : path)
	{
		std::cout << octant << " ";  // Output each octant in the path.
	}
	
	// Print the path in binary representation.
	std::cout << "\n\nBinary Path to Leaf Node: \n";
	for(auto& octant : path)
	{
		std::bitset<64> binary(octant);
		cout << binary;  // Output each octant as a binary number.
	}
	
	
	// Return the path vector containing the sequence of octants.
	return path;
}

















/**
 * Inserts a body into a given HOTNode of the Hashed Octree.
 *
 * The method starts by checking if the Hashed Octree is empty, and if so, it simply returns,
 * implying that the root node hasn't been initialized. If the provided node is nullptr,
 * this could mean one of two things:
 * 1. It's the first body being inserted, and we are at the root node.
 * 2. We are currently at the end of a branch in the octree.
 *
 * Once these preliminary checks are done, the method calculates a Morton key
 * for the given body's position, and then recursively inserts the body in
 * the appropriate node based on spatial subdivisions.
 *
 * Depending on the current node status (whether it's empty, contains a single body, or contains multiple bodies),
 * the function determines the appropriate octant for the body and inserts in
 * the appropriate node based on spatial subdivisions, potentially creating new child nodes.
 *
 * If the node is empty, then simply add the body's information to the node directly
 *
 * The Morton key is essentially a way of encoding the 3D position of the body
 * into a 1D value while preserving spatial locality. This is useful for spatial
 * searches and sorting operations.
 *
 * @param node: A reference to a pointer of the node into which the body is to be inserted.
 * @param body: The body to be inserted.
 * Note: The node parameter is passed by reference (&), allowing for direct modifications.
 *
 */

/**
 * When inserting a body into the hashed octree there are five general cases. The first two cases come from determining the starting node at which we will begin the tree traversal, 'while(node!=nullptr)'.
 * The next three general cases arise when navigating the tree to find the target node (the node that contains this body) and insert the body into the appropriate octant.
 *  For each traversal of the tree, 1 of the two first cases will be determine true, and 1 of the last three cases will be determined ture.
 *
 *
 * - **Initialization Cases**:
 FIRST TWO CASES: When traversing the hashed-octree to determine the starting point of the tree traversal, there two possible starting cases:
 
 GENERAL CASE 1. Node is nullptr
 - If the node in question is nullptr, we need to determine if we're at the root node or the end of a branch, the result of which will mean one of two things:
 1.1: We are on the first insertion and this is the root node
 1.2: We have just come from the end of a branch node (see 4.2, 5.2)
 - In both cases it means that the current key will be ROOT_KEY (1). It means the starting point of this traversal is the root node, i.e., we start the traversal from the beginning of the tree
 
 GENERAL CASE 2. Node already exists in hashed octree hash-map
 - If the node in question already exists in the hash-map it  means one of two things:
 2.1: We are in the middle of a tree traversal to find the target node (the node that contains this body) and our current node is the child of the parent node we just came from which was unable to have the body directly inserted into it (see 4.1, 5.1)
 2.2: We have encountered a collision in the hash-map
 
 *
 *
 *
 *
 LAST THREE CASES: When traversing the nodes stored in the hashmap of the hashed octree structure to find the target node, there are three general cases:
 
 GENERAL CASE 3. Node is empty, Handled in the insertBodyDirectly function
 - If the node is empty, it means we are at a leaf node or this is the first body being inserted. Just simply add the body's information to the node directly
 
 
 
 
 GENERAL CASE 4. Node has more than one body, Handled in the processMultiBodyNode function
 - If the node has more than one body, increment the body count and add the body mass to the node mass
 - Then, compute the octant of that node in which the current body resides, from which there will be two cases:
 4.1: The node's target octant exists, i.e., the desired has child node exists
 - If the node has children nodes, set the current key to the child key and set the current node to the target child node using lookup.
 - Then we will reenter the while loop from the beginning with the node in question having been set to it's pre-existing childnode in which the current body lies
 
 4.2: The node's target octant does NOT exist, i.e., the desired child node does NOT exist
 - If the target child node does not exist, a new child must be created with the current node acting as the parent node and with the bodies information (similar to the case when the node is empty, just creating the node into which it's being inserted this time)
 - Mark the existence of the child node in the target octant of the current node
 - Insert the child node into the hash-map of the hashed octree
 - Set the current node as a nullptr to indicate having reached the end of a branch (node with no children)
 - Then we will reenter the while loop from the beginning with the node in question being set the root node, meaning we are restarting from the beginning of the hashed octree and will
 
 
 
 
 GENERAL CASE 5. Node has exactly one body, Handled in the processSingleBodyNode function
 - If the Node contains only one body it's capacity has been reached, which means that the previous tree traversal for this node resulted in a body being directly inserted into the current node,
 and that we are now trying to insert another body into which also lies in the bounds of that node only it is being done on a different tree traversal(a different body being inserted),
 but because this node's capacity has been reached, we must pass the body already inside it into a child node which does not exist yet.
 - The target octant of the current node in which the current body lies is computed
 - The childNode is created with the current node acting as the parent node and with the bodies information
 - Mark the existence of the child node in the target octant of the current node
 - Insert the child node into the hash-map of the hashed octree
 - Now that space has been made for the this body in the current node(old body pushed down to child node), increment the body count and add the body mass to the node mass
 - Then, compute the octant of that node in which the current body resides, from which there will be two cases:
 5.1: The node's target octant exists, i.e., the desired has child node exists
 - If the node has children nodes, set the current key to the child key and set the current node to the target child node using lookup.
 - Then we will reenter the while loop from the beginning with the node in question having been set to it's pre-existing childnode in which the current body lies
 
 5.2: The node's target octant does NOT exist, i.e., the desired child node does NOT exist
 - If the target child node does not exist, a new child must be created with the current node acting as the parent node and with the bodies information (similar to the case when the node is empty, just creating the node into which it's being inserted this time)
 - Mark the existence of the child node in the target octant of the current node
 - Insert the child node into the hash-map of the hashed octree
 - Set the current node as a nullptr to indicate having reached the end of a branch (node with no children)
 - Then we will reenter the while loop from the beginning with the node in question being set the root node, meaning we are restarting from the beginning of the hashed octree and will
 
 
 *
 *
 *
 *
 *
 *
 *
 */
void LinearHashedOctree::insertBody(HOTNode*& node, const Vec3D bodyPosition, const double bodyMass)
{
	///*
	// Ensure the root node is initialized.
	assert(!nodes.empty());
	spatialKey currentKey;
	OctantEnum targetOctant;
	
	
	//GENERAL CASE 1
	if (node == nullptr)
	{
		currentKey = ROOT_KEY;
		node = lookUpNode(currentKey);
	}
	//GENERAL CASE 2
	else
	{
		currentKey = node->nodeKey;
	}
	
	
	
	// Traverse the tree to find the appropriate insertion point for the body.
	while (node != nullptr)
	{

		
		
		//GENERAL CASE 3: The target Node is empty
		if (node->N == 0)
		{
			node->mass = bodyMass;
			node->baryCenter = bodyPosition;
			node->N = 1;
			node->childByte = 0;
			
			node = nullptr;  // Reset the node pointer
		}
		else if (node->N > 1)         //GENERAL CASE 4: The target Node has Multiple Bodies
		{
			// Determine the spatial subdivision for the new body.
			targetOctant = DetermineOctant(node->nodeBounds.center, bodyPosition);
			node->N = node->N + 1;
			node->mass = node->mass + bodyMass;
			processNodeInsertion(node, targetOctant, currentKey, bodyPosition, bodyMass);
		}
		else if (node->N == 1)       //GENERAL CASE 5: The target Node has a Single Body
		{
			pushHOTNodeBodyToChild(node); // Push the existing body down to a child node, making space for a new body.
			targetOctant = DetermineOctant(node->nodeBounds.center, bodyPosition);
			node->N = node->N + 1;
			node->mass = node->mass + bodyMass;
			processNodeInsertion(node, targetOctant, currentKey, bodyPosition, bodyMass);
		}
	}
	
	//pruneEmptyNodes(node);
	
	//*/
}








/**
 * Transfer a body from a parent node to its appropriate child.
 *
 * This method is used when a node becomes too crowded (i.e., it contains more bodies
 * than it can accommodate) and needs to pass bodies down to its children. This maintains
 * the octree's property of having a set limit of bodies per node and ensures that bodies
 * are always found in the most specific (smallest) node that can contain them.
 *
 *
 * The method calculates which child node each body should be transferred
 * to based on its position and the spatial boundaries of the child nodes.
 * The child node is then either fetched if it already exists, or created
 * if it doesn't. Subsequently, the body is moved to the appropriate child node.
 *
 *
 * @param node: The current node that needs to push its body down to a child.
 */
void LinearHashedOctree::pushHOTNodeBodyToChild(HOTNode* node)
{
	// Check if node is valid and contains exactly one body. Exit if not.
	if (node == nullptr || node->N != 1)
	{
		return;
	}
	
	//Since N == 1, com is the position of that one body.
	//OctantEnum targetOctant = GetOctantFromKey(node->nodeKey);
	OctantEnum targetOctant = DetermineOctant(node->nodeBounds.center, node->baryCenter);
	
	// Create the child node in the determined octant with the body's mass and position.
	createChildNode(node, targetOctant, node->baryCenter, node->mass);
	
	
	//pruneEmptyNodes(node);
	//If the node is empty after pushing body to its child, delete the node
	//if (node->N == 0 && node->childByte == 0)
	//{
	//    nodes.erase(node->nodeKey);   // Remove the node from the map
	//    delete node;  // Release the node's memory
	//}
}


/**
 * Handles the insertion of a body into the correct node within the octree. This function
 * checks whether the target child node exists in the specified octant and either updates
 * the current node to that child or creates a new child node with the body's position
 * and mass. The function ensures that bodies are inserted into the most specific node
 * possible, maintaining the integrity of the octree structure.
 *
 * - node: A reference to a pointer to the current HOTNode in the insertion process.
 * - targetOctant: An OctantEnum indicating the octant into which the body should be inserted.
 * - currentKey: A reference to the spatialKey representing the current node's key.
 * - bodyPosition: A Vec3D representing the position of the body to be inserted.
 * - bodyMass: A double representing the mass of the body to be inserted.
 */
void LinearHashedOctree::processNodeInsertion(HOTNode*& node, OctantEnum& targetOctant, spatialKey& currentKey, const Vec3D bodyPosition, const double bodyMass)
{
	if (node == nullptr)
	{
		return;
	}
	
	// Check if the target child node exists in the specified octant.
	if (HasOctChild(node->childByte, targetOctant))
	{
		//updates the current key to that of the child's
		currentKey = GetChildKey(currentKey, targetOctant);
		
		//fetches the child node from the HOT
		node = lookUpNode(currentKey);
	}
	else
	{
		//create a new child node with the body's information
		createChildNode(node, targetOctant, bodyPosition, bodyMass);
		
		// Reset the node since we're now at a leaf.
		node = nullptr;
	}
	
}




/**
 * This function recursively traverses the octree, pruning nodes that are empty and
 * have no children. It helps to optimize the octree structure by removing unnecessary
 * nodes, thereby conserving memory and improving performance.
 */
void LinearHashedOctree::pruneEmptyNodes(HOTNode* node)
{
	if (node == nullptr)//|| node->nodeKey == 0)
	{
		return;
	}
	
	
	// Iterate over all possible child nodes
	for (int i = 0; i < 8; ++i)
	{
		// Check if there is a child node present at index i
		if (node->childByte & (1 << i))
		{
			spatialKey childKey = (node->nodeKey << 3) | i; // Compute the spatial key for the child node
			auto* childNode = lookUpNode(childKey); // Retrieve the child node using its spatial key
			
			// Recursively prune children of this child node
			pruneEmptyNodes(childNode);
		}
	}
	
	
	// Check if the current node is empty and has no children
	if (node->N == 0 && node->childByte == 0)
	{
		nodes.erase(node->nodeKey);   // Remove the node from the map
		delete node;  // Release the node's memory
	}
}


/**
 * This function deletes a specific node from the octree based on its spatial key.
 * It deallocates the memory associated with the node and removes the node from the
 * unordered_map that stores all nodes, ensuring no memory leaks occur.
 */
void LinearHashedOctree::deleteNode(spatialKey nodeCode)
{
	// Find the node in the nodes map using its spatial key
	auto iter = nodes.find(nodeCode);
	
	if (iter != nodes.end())
	{
		delete iter->second;  // Deallocate the memory for the node
		nodes.erase(iter);    // Remove the entry from the unordered_map
	}
}


/**
 * This function clears all nodes from the octree, effectively removing all elements
 * stored within the `nodes` map. It is used to reset the octree to an empty state
 * without deallocating memory for individual nodes, relying on the map's clear
 * function to handle this.
 */
void LinearHashedOctree::clear()
{
	nodes.clear();
}


/**
 * This function deallocates all dynamically allocated nodes in the octree and clears
 * the node map. It ensures that all memory used by the nodes is properly freed, preventing
 * memory leaks. The function iterates over all the nodes stored in the `nodes` map,
 * deletes each node, and then clears the map.
 */
void LinearHashedOctree::deleteTree()
{
	// Iterate over all node pairs in the `nodes` map
	for (auto& node_pair : nodes)
	{
		delete node_pair.second; // Deallocate the memory for the node
	}
	clear(); // Remove the entry from the unordered_map by clearing it
}




/**
 * This function creates a new child node within the octree for a specified parent node.
 * It calculates the spatial key for the child node, initializes it with given properties
 * (position and mass), updates the parent node to reflect the existence of this new child,
 * and inserts the newly created node into the octree.
 */
void LinearHashedOctree::createChildNode(HOTNode*& node, OctantEnum& targetOctant, const Vec3D& bodyPosition, const double& bodyMass)
{
	// Compute the spatial key for the new child node based on the parent's key and target octant
	spatialKey childKey = GetChildKey(node->nodeKey, targetOctant);
	
	HOTNode* childNode;  // Declare a pointer for the new child node
	childNode = new HOTNode(node, targetOctant, bodyPosition, bodyMass); // Allocate and initialize the new child node
	
	//Update the parent node to reflect the new child's existence.
	SetOctChild(node->childByte, targetOctant);
	
	
	// Insert the newly created HOTNode into the octree
	insertHOTNode(childNode);
}


/**
 * This function computes the barycenters and quadrupole moments for a node in the
 * octree. It recursively processes each child node, aggregating their masses and
 * barycenters to compute the barycenter for the parent node. The function also
 * calculates the quadrupole moments, which are essential for efficient gravitational
 * force calculations in hierarchical simulations, like those used in astrophysics.
 */
void LinearHashedOctree::computeNodeBaryCenters(HOTNode*& node)
{
	// Reset current node's mass and baryCenter.
	node->mass = 0.0;
	node->baryCenter = Vec3D(0.0, 0.0, 0.0);
	node->N = 0;
	
	//compute c.o.m. first.
	spatialKey childKey;
	HOTNode* childNode;
	double mass;
	for (size_t i = 0; i < 8; i++) //For all eight possible children
	{
		// Check if the ith child exists
		if (node->childByte & (1 << i))
		{
			childKey = (node->nodeKey << 3) | i; // Compute new Morton key for the child
			childNode = lookUpNode(childKey);            // Fetches the child node from the HOT
			computeNodeBaryCenters(childNode);		// Recursively compute barycenters for the child node
			
			
			// Accumulate mass and count for the parent node
			mass = childNode->mass;
			node->mass += mass;
			node->N += childNode->N;
			
			
			// Update barycenter weighted by child's mass
			node->baryCenter += (childNode->baryCenter * childNode->mass);
		}
		
		
		
		// Finalize barycenter calculation by dividing by the total mass
		node->baryCenter /= node->mass;
	}
	
	//compute quadrapole moment contribution
	node->quadrupoleMoment[0] = 0.0;
	node->quadrupoleMoment[1] = 0.0;
	node->quadrupoleMoment[2] = 0.0;
	node->quadrupoleMoment[3] = 0.0;
	node->quadrupoleMoment[4] = 0.0;
	node->quadrupoleMoment[5] = 0.0;
	
	Vec3D childPosition;
	double xx, yy, zz, d2;
	
	
	for (size_t i = 0; i < 8; i++) //For all eight possible children
	{
		// Check if the ith child exists
		if (node->childByte & (1 << i))
		{
			childKey = (node->nodeKey << 3) | i; //Compute new Morton key for the child
			childNode = lookUpNode(childKey);            //fetches the child node from the HOT
			
			//first, get child's pos in c.o.m. reference frame
			childPosition = childNode->baryCenter - node->baryCenter;
			
			mass = childNode->mass;
			
			
			// Calculate and add contribution to quadrupole moments
			xx = childPosition.x * childPosition.y; //xy
			node->quadrupoleMoment[1] += 3.0 * mass * xx;
			
			zz = childPosition.x * childPosition.z; //xz
			node->quadrupoleMoment[2] += 3.0 * mass * zz;
			
			
			yy = childPosition.y * childPosition.z; //yz
			node->quadrupoleMoment[4] += 3.0 * mass * yy;
			
			
			
			xx = childPosition.x * childPosition.x;
			zz = childPosition.z * childPosition.z;
			yy = childPosition.y * childPosition.y;
			
			d2 = xx + yy + zz;
			d2 *= mass;
			mass *= 3.0;
			
			
			node->quadrupoleMoment[0] += mass * xx - d2;
			node->quadrupoleMoment[3] += mass * yy - d2;
			node->quadrupoleMoment[5] += mass * zz - d2;
			
			
			//directly add child's moments also
			if (childNode->N > 1)
			{
				node->quadrupoleMoment[0] += childNode->quadrupoleMoment[0];
				node->quadrupoleMoment[1] += childNode->quadrupoleMoment[1];
				node->quadrupoleMoment[2] += childNode->quadrupoleMoment[2];
				node->quadrupoleMoment[3] += childNode->quadrupoleMoment[3];
				node->quadrupoleMoment[4] += childNode->quadrupoleMoment[4];
				node->quadrupoleMoment[5] += childNode->quadrupoleMoment[5];
			}
		}
		
		
		
	}
	
}

void LinearHashedOctree::computeTreeBaryCenters(HOTNode*& node)
{
	if (node == nullptr || nodes.empty() || node->N == 1)
	{
		return;
	}
	
	
	
	//do a depth-first search
	//recursively update children first.
	spatialKey childKey;
	for (int i = 0; i < 8; i++) //For all eight possible children
	{
		if (node->childByte & (1 << i))// See if the ith child exist
		{
			childKey = (node->nodeKey << 3) | i; //Compute new Morton key for the child
			auto* child = lookUpNode(childKey);  //Using key, look child up in hash table and recursively visit subtree
			computeTreeBaryCenters(child);
		}
	}
	
	
	// After recursing through all existing children in the tree(i.e., the last level of child nodes reached),
	computeNodeBaryCenters(node);
	
}






/*
 * This function evaluates the Multipole Acceptance Criterion (MAC) for a given node
 * in the octree with respect to a point in space. It determines whether the node can
 * be approximated as a single point mass based on its distance to the point and a
 * specified opening angle threshold.
 *
 * - node: A pointer to a HOTNode representing the node in the octree.
 * - r: A Vec3D representing the position vector to test against.
 * - theta: A double representing the opening angle threshold.
 */
bool LinearHashedOctree::MAC(HOTNode* node, const Vec3D& r, double theta)
{
	// Calculate the differences in x, y, and z coordinates between the node's barycenter and the given point
	double dx = node->baryCenter.x - r.x;
	double dy = node->baryCenter.y - r.y;
	double dz = node->baryCenter.z - r.z;
	
	double distanceSquared = dx * dx + dy * dy + dz * dz; // Compute the squared distance between the node's barycenter and the point
	double sizeSquared = node->nodeBounds.size * node->nodeBounds.size; // Compute the squared size of the node's bounding box
	
	
	// Return true if the ratio of size squared to distance squared is less than the square of theta
	return (sizeSquared / distanceSquared) < (theta * theta);
}



/**
 * This function prints detailed information about each node stored in the hashed octree.
 * It outputs the bucket count, size of the node container, and detailed data for each node,
 * including its spatial key, Morton key in binary, octant, node bounds, barycenter, mass,
 * number of bodies (N), and child byte.
 */
void LinearHashedOctree::printHashedOctree()
{
	cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nnodes.bucket_count(): " << nodes.bucket_count();
	cout << "\n\nnodes.size(): " << nodes.size();
	
	
	size_t i = 0;
	for (auto& pair : nodes)
	{
		i++;
		const spatialKey& key = pair.first;
		
		
		// Print out key
		
		std::cout << "\n\n\n\n\n\n\n\n\nKey: " << key;        cout << "                  i" << i;
		
		std::bitset<64> binary(key);
		cout << "\nBinary MortonKey: " << binary;
		
		
		
		if (pair.second != nullptr)
		{
			//OctantEnum octant1 = DetermineOctant(value->nodeBounds.center, value->baryCenter); //determine which octant the new body should be placed in
			//cout << "\n\nDetermineOctant: " << octant1;
			
			OctantEnum octant2 = GetOctantFromKey(pair.second->nodeKey);
			cout << "\n\n\n\nGetOctantFromKey: " << octant2;
			
			cout << "\nvalue->nodeBounds.size: " << pair.second->nodeBounds.size;
			cout << "\nvalue->nodeBounds.center: (" << pair.second->nodeBounds.center.x; cout << ", " << pair.second->nodeBounds.center.y; cout << ", " << pair.second->nodeBounds.center.z; cout << ")";
			
			
			cout << "\nvalue->baryCenter: (" << pair.second->baryCenter.x; cout << ", " << pair.second->baryCenter.y; cout << ", " << pair.second->baryCenter.z; cout << ")";
			cout << "\nvalue->mass: " << pair.second->mass;
			
			cout << "\nvalue->N: " << pair.second->N;
			cout << "\nvalue->childByte: " << std::bitset<8>(pair.second->childByte);
			
		}
	}
}

void LinearHashedOctree::visualizeTree()
{
	for (auto& pair : nodes)
	{
		visualizeNode(pair.second);
	}
}

void LinearHashedOctree::visualizeNode(const HOTNode* node)
{
	if (!node)
	{
		return;
	}
	
	
	ofNoFill();
	// Get the center and size of the node
	//ofVec3f center(node->nodeBounds.center.x, node->nodeBounds.center.y, node->nodeBounds.center.z);
	
	
	// Set the box color
	//ofSetColor(255, 255, 255, 63.75);
	ofSetColor(255, 255, 255, 47.8125);
	//ofSetColor(255, 255, 255, 31.875);
	
	// Draw the box
	ofDrawBox(node->nodeBounds.center.x, node->nodeBounds.center.y, node->nodeBounds.center.z, node->nodeBounds.size); //uses the passed-in coordinates as the center of the cube
	
	
	//ofSetColor(66, 66, 255);
	//ofDrawBitmapString("Mass: " + ofToString(node->mass) + "\nN: " + ofToString(node->N) + "\nQuadrupole Moment[4]: " + ofToString(node->quadrupoleMoment[4]), center);
	
	//size_t depthLevel = GetNodeTreeDepth(node->nodeKey);
	//ofDrawBitmapString("Mass: " + ofToString(node->mass) + "\nN: " + ofToString(node->N) + "\nKey: " + ofToString(node->nodeKey) + "\nDepth: " + ofToString(depthLevel), center);
}
