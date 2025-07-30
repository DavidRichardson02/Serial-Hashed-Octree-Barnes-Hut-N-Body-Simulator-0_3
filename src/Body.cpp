//  Body.cpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3


#include "Body.hpp"




// Default constructor
Body::Body() : position(0.0, 0.0, 0.0), velocity(0.0, 0.0, 0.0), mass(0.0), bodyKey(0), bodyColor(0,255,0), displayRadius(mass) {}

// Overloaded constructors
Body::Body(Vec3D _position, Vec3D _velocity, double _mass) : position(_position), velocity(_velocity), mass(_mass), bodyKey(0), bodyColor(0,255,0), displayRadius(_mass) {}

Body::Body(Vec3D _position, Vec3D _velocity, double _mass, ofColor _bodyColor) : position(_position), velocity(_velocity), mass(_mass), bodyKey(0), bodyColor(_bodyColor), displayRadius(_mass){}

Body::Body(Vec3D _position, Vec3D _velocity, double _mass, ofColor _bodyColor, double _displayRadius): position(_position), velocity(_velocity), mass(_mass), bodyKey(0), bodyColor(_bodyColor), displayRadius(_displayRadius){}
// Copy constructor
Body::Body(const Body& other) : position(other.position), velocity(other.velocity), mass(other.mass), bodyKey(other.bodyKey), bodyColor(other.bodyColor), displayRadius(other.displayRadius) {}

// Assignment operator
Body& Body::operator=(const Body& other)
{
	if (this != &other)
	{
		position = other.position;
		velocity = other.velocity;
		mass = other.mass;
		bodyKey = other.bodyKey;
		bodyColor = other.bodyColor;
		displayRadius = other.displayRadius;
	}
	return(*this);
}

// Setters and getters
void Body::setBodyKey(spatialKey _bodyKey) { bodyKey = _bodyKey; }
spatialKey Body::getBodyKey() const { return(bodyKey); }


// Debug utility implementation
void Body::printBodies(string label)
{
	cout << "\n" << label; cout << "_";
	cout << "    \nPosition: (" << position.x; cout << ", " << position.y; cout << ", " << position.z; cout << ")";
	cout << "    \n     Mass: " << mass;
	
	//cout << "    \n\nInteger MortonKey: " << mortonKey;
	// Convert the Morton key to binary format and print it
	//std::bitset<64> binary(mortonKey);
	//cout << "\nBinary MortonKey: " << binary;
}
