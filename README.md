NOTE: The documentation present in the files of this repo has been severely reduced to increase readability while I use it. Though very thorough, it is also very information dense, and makes for a cluttered work environment and will remain absent until I get around to collecting all the individual pieces of documentation spread throughout the files and implement a sort of footnote system to link code sections to the extensive documentation while minimally decreasing readability. 


This project is part of a directed study aimed at optimizing the Barnes-Hut algorithm for N-body simulations using Hashed Octrees and various performance-enhancing techniques. The N-body problem, which deals with the calculation of gravitational forces in a system of multiple bodies, faces challenges related to computational complexity, traditionally requiring O(N²) calculations. The Barnes-Hut algorithm reduces this to O(N log N) by approximating the effect of distant particles, allowing for more scalable simulations.

The primary focus of this project is the development and exploration of hashed octrees, a data structure optimized for parallel construction and traversal, improving both memory efficiency and computational speed for large-scale N-body simulations. This simulator is built using OpenFrameworks and is modularly designed to allow experimentation with different optimization strategies.

## Key Features
  - Hashed Octree Data Structure: Efficient spatial partitioning using hashed octrees instead of traditional pointer-based octrees. This improves parallelism, memory usage, and cache efficiency.
  - Barnes-Hut Algorithm: Reduces the complexity of N-body gravitational simulations from O(N²) to O(N log N).
  - Parallelization Support: Built with an emphasis on parallel computation using multi-threading techniques, making the simulation scalable for large datasets.
  - OpenFrameworks Integration: A graphical interface for visualizing the octree construction, the placement of bodies, and the results of the simulation in both 2D and 3D environments.
  - Dynamic Refinement: Adaptive refinement of the octree based on the distribution of bodies, allowing for efficient handling of simulations where bodies undergo significant movement.


## Use Cases
  - Astrophysical Simulations: Model the gravitational interactions of star clusters or galaxies over time.
  - 3D Graphics and Collision Detection: The octree data structure is applicable in fields like 3D graphics for object rendering, collision detection, and even in 3D convolutional neural networks for shape analysis.
  - Geographic Information Systems (GIS): The efficient partitioning of space using octrees can be useful in spatial data analysis for GIS applications.


## Core Logic Structure
  - HashedOctree.hpp & HashedOctree.cpp: The core of the spatial data structure, which implements the hashed octree for efficient memory management and traversal.
  - HashedNode.hpp & HashedNode.cpp: Defines the individual nodes within the octree, implementing the hashing mechanism to locate nodes efficiently.
  - Body.hpp: Represents the individual bodies in the simulation, each with properties like mass, position, velocity, and acceleration.
  - PhysicsHelpers.hpp: Contains helper functions for calculating gravitational forces, velocities, and accelerations between bodies.
  - MortonKeys.hpp: Implements the Morton encoding (Z-order curve) to linearize spatial hierarchies for efficient processing.
  - OctantUtils.hpp & OctantUtils.cpp: Utilities for managing octants (subregions) of the simulation space, helping in the recursive division of space.
