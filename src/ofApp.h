#pragma once
#include "MathUtils.hpp"
#include "SequenceContainers.hpp"
#include "MortonKeys.hpp"
#include "ObjectPool.hpp"
#include "Body.hpp"
#include "OctantUtils.hpp"
#include "HashedNode.hpp"
#include "HashedOctree.hpp"
#include "PhysicsHelpers.hpp"
#include "ofMain.h"




static inline void plummerModelInitialConditions(Body* &simBodies, size_t numBodies, double totalMass, double anchorMass, double galaxyScatter, double orbitalRadius, Vec3D galaxyAnchorPoint);
//distribute the mass of the bodies following a spherically symmetric density profile, position the bodies by sampling the cumulative distribution function derived from the Plummer density profile. Velocities of bodies are determined such that the system is in virial equilibrium, i.e., kinetic and potential energy of the system is balanced
static inline void coldStartInitialConditions(Body* &simBodies, size_t numBodies);




// Structure to hold a frustum plane
struct FrustumPlane
{
	glm::vec3 normal;
	float distance;
};

inline void calculateFrustumPlanes(double horizontalAngle, double verticalAngle, const ofVec3f& offset, double zoomScale, FrustumPlane (&planes)[6]);//const ofCamera& camera, FrustumPlane planes[6]);
bool isPointInFrustum(const glm::vec3& point, const FrustumPlane planes[6]) ;




//--------------------------------------------------------------
// Suppose you want to allow 10^+3 down to 10^-10, etc.
// We'll store them in a vector or array. You can go bigger/smaller.
static std::vector<double> ZOOM_LEVELS =
{
	1e3, 1e2, 1e1, 1.0,
	1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6,
	1e-7, 1e-8, 1e-9, 1e-10, 1e-11
	// Add more if needed
};
static int currentZoomIndex = 3; // Start at 1.0 (which is ZOOM_LEVELS[3]).
static double zoomScale = ZOOM_LEVELS[currentZoomIndex];





class ofApp : public ofBaseApp
{
	size_t numBodies;
	Body* bodies;
	OctantBounds rootNodeBounds;
	LinearHashedOctree LHTree;
	
	double theta = 1;
	long interactionCount = 0, numInteractions = 0;
	double dt = 0.001;
	
	Vec3D* bodiesAccelerations;
	
	
public:
	
	/*----------------------   2D plane coordinate system navigation  ----------------------*/
	ofRectangle canvas;
	const double zoomBuffer = 0.1; // 10% buffer
	double maxZoomScale;//0.01;  //the maximum viewable plane size will be ten times that of the initial plane size
	double minZoomScale;//0.000000000000001;  //the minimal viewable plane size will be one tenth that of the initial plane size
	//double zoomScale = 1.0;  //the variable to control the value of each scroll wheel input, i.e, how much each scroll changes the view
	double lastZoomScale = zoomScale;
	double zoomIncrement = 0.025f;//0.0000000001;  //the variable to control the value of each scroll wheel input, i.e, how much each scroll changes the view
	ofVec3f dragStartPos;
	ofVec3f offset;
	double horizontalAngle = 0.0;
	double verticalAngle = 0.0;
	bool isDragging = false;
	bool visualizeTree = false;
	
	
	float horizontalAngleAtMousePress;
	float verticalAngleAtMousePress;
	
	
	ofVec3f panningVelocity;
	double zoomVelocity = 0.0;
	double horizontalRotationVelocity = 0.0;
	double verticalRotationVelocity = 0.0;
	double interfaceInertiaDamping = 0.9;  // adjust this for stronger/weaker inertia
	
	// Add a member variable to store the offset when the mouse is pressed
	ofVec3f offsetAtMousePress;
	
	ofVec3f targetOffset;
	bool isFocusing = false;
	double focusSpeed = 0.05;  // Control how fast the focusing happens
	
	
	
	void coldStartInitialConditionzzz(Body* &simBodies, size_t numBodies);
	void setup() override;
	void update() override;
	void draw() override;
	void exit() override;
	
	void keyPressed(int key) override;
	void keyReleased(int key) override;
	void mouseMoved(int x, int y ) override;
	void mouseDragged(int x, int y, int button) override;
	void mousePressed(int x, int y, int button) override;
	void mouseReleased(int x, int y, int button) override;
	void mouseScrolled(int x, int y, float scrollX, float scrollY) override;
	void mouseEntered(int x, int y) override;
	void mouseExited(int x, int y) override;
	void windowResized(int w, int h) override;
	void dragEvent(ofDragInfo dragInfo) override;
	void gotMessage(ofMessage msg) override;
	
};









//distribute the mass of the bodies following a spherically symmetric density profile, position the bodies by sampling the cumulative distribution function derived from the Plummer density profile. Velocities of bodies are determined such that the system is in virial equilibrium, i.e., kinetic and potential energy of the system is balanced
static inline void plummerModelInitialConditions(Body* &simBodies,size_t numBodies, double totalMass, double anchorMass, double galaxyScatter, double orbitalRadius, Vec3D galaxyAnchorPoint)
{
	double mass, uniformRandom, r, _theta, phi, x, y, z, u1, u2, z1, v, vx, vy, vz;
	Vec3D position, velocity;
	mass = totalMass / numBodies + anchorMass / numBodies;
	for (size_t i = 0; i < numBodies; i++)
	{
		mass = totalMass / numBodies + anchorMass / numBodies;
		uniformRandom = ofRandom(galaxyScatter);          // Generate uniform random number between 0 and 1(between 0 and galaxyScatter)
		r = orbitalRadius * sqrt(pow(uniformRandom, -2.0 / 3) - 1);
		_theta = 2 * PI * ofRandom(1.0);      // Generate random angles uniformly distributed
		phi = acos(1 - 2 * ofRandom(1.0));
		x = r * sin(phi) * cos(_theta);
		y = r * sin(phi) * sin(_theta);
		z = r * cos(phi);
		
		
		
		position = Vec3D(x , y , z ); // Changed to ofVec3f for 3D
		
		// Generate normal distributed random numbers using Box-Muller transformation
		u1 = ofRandom(1.0);
		u2 = ofRandom(1.0);
		z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2); // Box-Muller transformation
		v = sqrt(2) * pow((1 + r * r), -0.25);
		vx = v * sqrt(2) * pow(1 + r * r, -0.25) * z1;
		vy = sqrt(v * v - vx * vx) * (ofRandom(1.0) > 0.5 ? 1 : -1);
		vz = sqrt(v * v - vx * vx - vy * vy) * (ofRandom(1.0) > 0.5 ? 1 : -1);
		
		velocity = Vec3D(vx, vy,  ofRandom(-1.0,1.0)); //vz);
		
		
		simBodies[i] = {position, velocity, 1}; // mass};
	}
}



static inline void coldStartInitialConditions(Body* &simBodies, size_t numBodies)
{
	Vec3D position, velocity;
	double mass = 1.0;
	for (size_t i = 0; i < numBodies; i++)
	{
		//position = Vec3D((0 + rand() % (1500)), (0 + rand() % (1500)), (0 + rand() % (1500)));
		//velocity = Vec3D((0 + rand() % (1)), (0 + rand() % (1)), (0 + rand() % (1)));
		
		position = Vec3D(ofRandom(-1500, 1500), ofRandom(-1500, 1500), ofRandom(-1500, 1500));
		velocity = Vec3D(ofRandom(-1, 1), ofRandom(-1, 1), ofRandom(-1, 1));
		
		
		simBodies[i] = {position, velocity, mass};
	}
	
}






/**
 If a node is entirely outside, skip its rendering.
 
 Calculate Frustum Planes:
 Each frame, calculate the six frustum planes (top, bottom, left, right, near, far) based on the camera’s position, orientation, and perspective settings.
 Frustum Check:
 For each body or node, check its position relative to the frustum planes. If it lies outside any plane, skip rendering it.

 */




/*
// Calculate frustum planes based on the camera's view-projection matrix
void calculateFrustumPlanes(const ofCamera& camera, FrustumPlane planes[6]) {
	glm::mat4 m = camera.getProjectionMatrix() * camera.getModelViewMatrix();
	
	// Left plane
	planes[0].normal = glm::vec3(m[0][3] + m[0][0], m[1][3] + m[1][0], m[2][3] + m[2][0]);
	planes[0].distance = m[3][3] + m[3][0];
	
	// Right plane
	planes[1].normal = glm::vec3(m[0][3] - m[0][0], m[1][3] - m[1][0], m[2][3] - m[2][0]);
	planes[1].distance = m[3][3] - m[3][0];
	
	// Bottom plane
	planes[2].normal = glm::vec3(m[0][3] + m[0][1], m[1][3] + m[1][1], m[2][3] + m[2][1]);
	planes[2].distance = m[3][3] + m[3][1];
	
	// Top plane
	planes[3].normal = glm::vec3(m[0][3] - m[0][1], m[1][3] - m[1][1], m[2][3] - m[2][1]);
	planes[3].distance = m[3][3] - m[3][1];
	
	// Near plane
	planes[4].normal = glm::vec3(m[0][3] + m[0][2], m[1][3] + m[1][2], m[2][3] + m[2][2]);
	planes[4].distance = m[3][3] + m[3][2];
	
	// Far plane
	planes[5].normal = glm::vec3(m[0][3] - m[0][2], m[1][3] - m[1][2], m[2][3] - m[2][2]);
	planes[5].distance = m[3][3] - m[3][2];
	
	// Normalize the planes for accurate distance checks
	for (int i = 0; i < 6; i++) {
		float length = glm::length(planes[i].normal);
		planes[i].normal /= length;
		planes[i].distance /= length;
	}
}
*/
// Use custom view-projection matrix for frustum calculation
inline void calculateFrustumPlanes(double horizontalAngle, double verticalAngle, const ofVec3f& offset, double zoomScale, FrustumPlane (&planes)[6])
{
	glm::mat4 projectionMatrix = glm::perspective(glm::radians(60.0f), (float)ofGetWidth() / (float)ofGetHeight(), 1.0f, 5000.0f); // Customize FOV, near, and far planes as needed
	
	glm::mat4 viewMatrix = glm::rotate(glm::mat4(1.0f), (float)horizontalAngle, glm::vec3(0, 1, 0))  // Horizontal rotation
	* glm::rotate(glm::mat4(1.0f), (float)verticalAngle, glm::vec3(1, 0, 0))    // Vertical rotation
	* glm::translate(glm::mat4(1.0f), -glm::vec3(offset.x, offset.y, offset.z)) // Offset/panning
	* glm::scale(glm::mat4(1.0f), glm::vec3(zoomScale));                        // Zoom scale
	
	glm::mat4 m = projectionMatrix * viewMatrix;
	
	// Calculate frustum planes using the combined view-projection matrix `m`
	// Same frustum plane extraction code as before
	planes[0].normal = glm::vec3(m[0][3] + m[0][0], m[1][3] + m[1][0], m[2][3] + m[2][0]);
	planes[0].distance = m[3][3] + m[3][0];
	planes[1].normal = glm::vec3(m[0][3] - m[0][0], m[1][3] - m[1][0], m[2][3] - m[2][0]);
	planes[1].distance = m[3][3] - m[3][0];
	planes[2].normal = glm::vec3(m[0][3] + m[0][1], m[1][3] + m[1][1], m[2][3] + m[2][1]);
	planes[2].distance = m[3][3] + m[3][1];
	planes[3].normal = glm::vec3(m[0][3] - m[0][1], m[1][3] - m[1][1], m[2][3] - m[2][1]);
	planes[3].distance = m[3][3] - m[3][1];
	planes[4].normal = glm::vec3(m[0][3] + m[0][2], m[1][3] + m[1][2], m[2][3] + m[2][2]);
	planes[4].distance = m[3][3] + m[3][2];
	planes[5].normal = glm::vec3(m[0][3] - m[0][2], m[1][3] - m[1][2], m[2][3] - m[2][2]);
	planes[5].distance = m[3][3] - m[3][2];
	
	// Normalize planes for accurate calculations
	for (int i = 0; i < 6; i++) {
		float length = glm::length(planes[i].normal);
		planes[i].normal /= length;
		planes[i].distance /= length;
	}
}




//to determine if a body is inside the frustum by testing it against each plane. If it’s outside even one plane, it’s outside the frustum.
// Check if a point (e.g., body position) is inside the frustum
bool isPointInFrustum(const glm::vec3& point, const FrustumPlane planes[6])
{
	for (int i = 0; i < 6; i++)
	{
		if (glm::dot(planes[i].normal, point) + planes[i].distance < 0)
		{
			return false; // Outside this plane
		}
	}
	return true; // Inside all planes
}





















/*
 
 to do:
 - refactor the function in my  'PhysicHelpers.f' file to scale the solar system model to be able to be visualized in a...
 while maintaining the same proportion of masses and distances, just to this much smaller scale.
 
 
 //*/


/**
 * Vec3D Class: provides an abstraction for three-dimensional vectors, encapsulating the vector's components (x, y, z)
 *
 * This class encapsulates a 3D vector's components (x, y, z),  and provides methods
 * to perform basic vector operations and manipulate vector data.
 *
 #pragma once
 
 
 
 
 // Vec3D class representing a 3D Vector
 class Vec3D
 {
 public:
 // ------------- Constructors -------------
 Vec3D();
 Vec3D(double _x, double _y, double _z);
 
 
 
 
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
 double vectorLength() const;
 double vectorSquareLength() const;
 Vec3D vectorNormalize() const;
 
 
 
 
 // ------------- Member variables -------------
 //private:
 double x;
 double y;
 double z;
 };
 
 
 
 #include "Containers.h"
 #include <math.h>
 
 Vec3D::Vec3D() : x(0), y(0), z(0) {}
 
 Vec3D::Vec3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) // using an initializer list to initialize members instead of assignment in the body of the constructor
 {}
 
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
 
 
 
 
 
 
 
 
 
 
 
 #pragma once
 #include "Vec3D.h"
 #include "ofMain.h"
 
 /**
 *  Body class representing a physical body with position, velocity, and mass properties
 *
 class Body
 {
 public:
 // ------------- Constructors and Assignment operators -------------
 Body();
 Body(Vec3D _position, Vec3D _velocity, double _mass);
 Body(Vec3D _position, Vec3D _velocity, double _mass, const double _size);  //compute bodies based on an inputted position and a defined bounding box
 Body(const Body& other);
 Body& operator=(const Body& other);
 
 
 
 
 
 
 
 // ------------- Member variables -------------
 //private:
 Vec3D position;
 Vec3D velocity;
 double mass;
 
 };
 
 
 #include "Body.h"
 
 
 // Default constructor
 Body::Body() : position(0.0, 0.0, 0.0), velocity(0.0, 0.0, 0.0), mass(0.0) {}
 
 // Overloaded constructors
 Body::Body(Vec3D _position, Vec3D _velocity, double _mass) : position(_position), velocity(_velocity), mass(_mass) {}
 Body::Body(Vec3D _position, Vec3D _velocity, double _mass, const double _size) : position(_position), velocity(_velocity), mass(_mass)
 {
 }
 
 // Copy constructor
 Body::Body(const Body& other) : position(other.position), velocity(other.velocity), mass(other.mass) {}
 
 // Assignment operator
 Body& Body::operator=(const Body& other)
 {
 if (this != &other)
 {
 position = other.position;
 velocity = other.velocity;
 mass = other.mass;
 }
 return(*this);
 }
 
 
 
 
 
 #pragma once
 #include "Vec3D.h"
 #include "Body.h"
 #include "ofMain.h"
 #include <vector>
 #include <cmath>
 
 const double PI = 3.14159265358979323846;
 
 // Function to set up the Solar System
 void setupSolarSystem(Body* &simBodies, size_t numBodies) {
 
 size_t idx = 0;
 
 // Initialize the Sun
 if (idx < numBodies) {
 simBodies[idx++] = {Vec3D(0, 0, 0), Vec3D(0, 0, 0), 1.989e30};
 }
 
 // Initialize the Planets (masses in kg, positions and velocities in m and m/s)
 struct PlanetInit {
 double mass;
 Vec3D position;
 Vec3D velocity;
 };
 
 std::vector<PlanetInit> planetData = {
 {3.3011e23, Vec3D(5.79e10, 0, 0), Vec3D(0, 47400, 0)},  // Mercury
 {4.8675e24, Vec3D(1.0821e11, 0, 0), Vec3D(0, 35000, 0)},  // Venus
 {5.97237e24, Vec3D(1.496e11, 0, 0), Vec3D(0, 29800, 0)},  // Earth
 {6.4171e23, Vec3D(2.279e11, 0, 0), Vec3D(0, 24100, 0)},  // Mars
 {1.8982e27, Vec3D(7.783e11, 0, 0), Vec3D(0, 13100, 0)},  // Jupiter
 {5.6834e26, Vec3D(1.426e12, 0, 0), Vec3D(0, 9600, 0)},  // Saturn
 {8.6810e25, Vec3D(2.871e12, 0, 0), Vec3D(0, 6800, 0)},  // Uranus
 {1.02413e26, Vec3D(4.498e12, 0, 0), Vec3D(0, 5400, 0)}  // Neptune
 };
 
 for (const auto& planet : planetData) {
 if (idx >= numBodies) break;
 simBodies[idx++] = {planet.position, planet.velocity, planet.mass};
 }
 
 // Initialize Moons (Just for Earth and Jupiter, simplified)
 // Earth's Moon
 if (idx < numBodies) {
 simBodies[idx++] = {Vec3D(1.496e11 + 3.844e8, 0, 0), Vec3D(0, 29800 + 1022, 0), 7.34767309e22};
 }
 // Jupiter's Galilean Moons (Io, Europa, Ganymede, Callisto)
 double jupiter_x = 7.783e11;
 double jupiter_v = 13100;
 std::vector<PlanetInit> moonData = {
 {8.9319e22, Vec3D(jupiter_x + 4.22e8, 0, 0), Vec3D(0, jupiter_v + 17334, 0)},  // Io
 {4.7998e22, Vec3D(jupiter_x + 6.71e8, 0, 0), Vec3D(0, jupiter_v + 13740, 0)},  // Europa
 {1.4819e23, Vec3D(jupiter_x + 1.07e9, 0, 0), Vec3D(0, jupiter_v + 10880, 0)},  // Ganymede
 {1.0759e23, Vec3D(jupiter_x + 1.88e9, 0, 0), Vec3D(0, jupiter_v + 8204, 0)}  // Callisto
 };
 
 for (const auto& moon : moonData) {
 if (idx >= numBodies) break;
 simBodies[idx++] = {moon.position, moon.velocity, moon.mass};
 }
 
 // Initialize Asteroid Belt (Approximate positions)
 for (; idx < numBodies * 0.9 && idx < numBodies; ++idx) {
 double theta = 2 * PI * ((double) rand() / RAND_MAX);
 double r = 3e11 + ((double) rand() / RAND_MAX) * 1e11;  // Between Mars and Jupiter
 simBodies[idx] = {Vec3D(r * cos(theta), r * sin(theta), 0), Vec3D(0, 0, 18000), 1e21};
 }
 
 // Initialize Oort Cloud (randomly in a sphere)
 for (; idx < numBodies; ++idx) {
 double theta = 2 * PI * ((double) rand() / RAND_MAX);
 double phi = acos(2 * ((double) rand() / RAND_MAX) - 1);
 double r = 1e13 + ((double) rand() / RAND_MAX) * 2e13;  // Up to 3e13 m
 simBodies[idx] = {Vec3D(r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)), Vec3D(0, 0, 500), 2e13};
 }
 }
 
 
 
 
 
 
 
 
 #pragma once
 #include "Vec3D.h"
 #include "Body.h"
 #include "PhysicsHelpers.h"
 #include "ofMain.h"
 
 
 class ofApp : public ofBaseApp
 {
 
 size_t numBodies;
 Body* bodies;
 OctantBounds rootNodeBounds;
 
 
 
 void setup() override;
 void update() override;
 void draw() override;
 
 };
 
 
 
 
 */
