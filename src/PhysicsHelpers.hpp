//  PhysicsHelpers.hpp
//  Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3

#pragma once
#include <vector>
#include <cmath>
#include "MathUtils.hpp"
#include "SequenceContainers.hpp"
#include "MortonKeys.hpp"
#include "Body.hpp"
//#include "DataLoader.hpp"
#include "ofMain.h"



extern "C" {
#include "CommonDefinitions.h"
#include "GeneralUtilities.h"
#include "StringUtilities.h"
#include "FileUtilities.h"
#include "DataExtraction.h"
#include "SolarSystem.h"
}





// Constants
const double G = 6.67430e-11;  // Gravitational constant in m^3 kg^-1 s^-2
const double SOLAR_MASS = 1.989e30; // Mass of the Sun in kg
const double MU = G * SOLAR_MASS; // Gravitational parameter;
const double AU = 149597870700; // Astronomical Unit

// Constants for display radius calculations
const double DISPLAY_RADIUS_SCALE_LARGE = 1000000000;
const double DISPLAY_RADIUS_SCALE_MEDIUM = 100000000;
const double DISPLAY_RADIUS_SCALE_SMALL = 10000000;
const double DISPLAY_RADIUS_BASE = 1.0e2;
//const size_t LARGE_BODY_THRESHOLD = 14;
const double LARGE_BODY_THRESHOLD = 1e25;
const double MEDIUM_BODY_THRESHOLD = 1e23;
const double SMALL_BODY_THRESHOLD = 1e15;




/* Class to hold the orbital elements for each celestial object read in to the program
 
 //*----------------------   Step 1: Source the Data Files   ---------------------- *
 * The orbital data for celestial objects is READ from the downloaded data files from an online repository
 
 * After reading from the online repository, the yet to be created specific parser function WRITES the specified data
 to a file on the computer in the data folders of this project
 
 NOTE: Use data from a reliable source. The NASA Jet Propulsion Laboratory (JPL) Horizons system is a good place to start.
 This system allows you to generate ephemerides for solar-system bodies that include all the required parameters.
 
 
 
 
 //*----------------------   Step 2: Store the Object Data   ---------------------- *
 * After the each of the object's data has been written to the on-site file, the program will then read in each object's data to be
 used as an argument list for initializing an instance of the CelestialObject class, which will be used to store the orbital state of the object for a
 specified point in time.
 
 
 
 //*----------------------   Step 3: Setup the Object Model   ---------------------- *
 * After populating the CelestialBody class instances using the data written in from the on-site file,
 compute the cartesian coordinates and velocities from each object's orbital elements
 
 //*----------------------   Step 4: Setup the Object Model   ---------------------- *
 * Use the CelestialBody object's to setup a model of the solar system that allows users to easily recgonize smaller bodies while maintaining the
 integrity of the simulation, this is accomplished by a including a seperate member variable in the Body class that represent's the radius of the body when being visualized, so that all bodies can be made an size logarithmically proportional to it's mass that makes them visible while not harming the actual physics.
 
 
 
 
 * This CelestialObject will then
 
 * Param: eccentricity = e : shape of the ellipse, describing how much it is elongated compared to a circle
 *
 * Param: semiMajorAxis = a :  the sum of the periapsis and apoapsis distances divided by two. For classic two-body orbits, the semimajor axis is the distance
 * between the centers of the bodies, not the distance of the bodies from the center of mass.
 *
 * Param: inclination = i : vertical tilt of the ellipse with respect to the reference plane, measured at the ascending node
 *
 *
 * Param: longitudeOfAscendingNode = Omega : horizontally orients the ascending node of the ellipse (where the orbit passes from south to north through the
 * reference plane, symbolized by ☊) with respect to the reference frame's vernal point (symbolized by ♈︎). This is measured in the reference plane
 *
 * Param: argumentOfPerihelion = omega : argument of periapsis, defines the orientation of the ellipse in the orbital plane, as an angle measured from the
 * ascending node to the periapsis
 *
 * Param: meanAnomaly = M : a mathematically convenient fictitious "angle" which varies linearly with time, but which does not correspond to a real geometric
 * angle. It can be converted into the true anomaly ν, which does represent the real geometric angle in the plane of the ellipse, between periapsis (closest
 * approach to the central body) and the position of the orbiting object at any given time.
 *
 */
class CelestialObject
{
public:
	double semiMajorAxis;
	double eccentricity;
	double inclination;
	double meanAnomaly;
	double argumentOfPerihelion;
	double longitudeOfAscendingNode;
	
	double mass;
	
	Vec3D position;
	Vec3D velocity;
	
	ofColor planetColor;
};


static inline void readOrbitalElementsFromFile();
static inline void SetupSolarSystem(Body* &simBodies, size_t numBodies);
static inline void initializeBody(Body &body, double semiMajorAxis, double eccentricity, double inclination, double meanAnomaly, double argumentOfPerihelion, double longitudeOfAscendingNode, double planetMass);
static inline void GetCartesianFromOrbital(Body &body, double eccentricity, double semiMajorAxis, double inclination, double longitudeOfAscendingNode, double argumentOfPerihelion, double meanAnomaly, double mu);



double toRadians(double degrees);
void calculatePositionAndVelocity(Body &body, double r, double ν, double i, double Ω, double ω, double mu, double a, double e);
static inline void GetCartesianFromOrbitalRadians(Body &body, double e, double a, double i, double Ω, double ω, double M, double mu);
static inline void GetCartesianFromOrbitalDegrees(Body &body, double e, double a, double i_deg, double Ω_deg, double ω_deg, double M, double mu);



//
/*
// Function to set up the Solar System
static inline void SetupSolarSystem(Body* &simBodies, size_t numBodies)
{
	
	
	
	
	
	size_t idx = 0;
	
	// Initialize the Sun
	if (idx < numBodies)
	{
		simBodies[idx++] = {Vec3D(0, 0, 0), Vec3D(0, 0, 0), 1.989e30, ofColor(253,184,19)};
	}
	
	// Initialize the Planets (masses in kg, positions and velocities in m and m/s)
	struct PlanetInit {
		double mass;
		Vec3D position;
		Vec3D velocity;
		ofColor planetColor;
	};
	
	std::vector<PlanetInit> planetData =
	{
		{3.3011e23 , Vec3D(5.79e10), Vec3D(47400), ofColor(205,203,202)}, //Mercury
		{4.8675e24 , Vec3D(1.0821e11), Vec3D(35000), ofColor(238,208,83)},  // Venus
		{5.97237e24 , Vec3D(1.496e11), Vec3D(29800), ofColor(52,198,144)},  // Earth
		{6.4171e23 , Vec3D(2.279e11), Vec3D(24100), ofColor(173,98,66)},  // Mars
		{1.8982e27 , Vec3D(7.783e11), Vec3D(13100), ofColor(225,225,226)},  // Jupiter
		{5.6834e26 , Vec3D(1.426e12), Vec3D(9600), ofColor(250,229,191)},  // Saturn
		{8.6810e25 , Vec3D(2.871e12), Vec3D(6800), ofColor(225,238,238)},  // Uranus
		{1.02413e26, Vec3D(4.498e12), Vec3D(5400), ofColor(91,93,223)}  // Neptune
		
	};
	
	for (const auto& planet : planetData) {
		if (idx >= numBodies) break;
		simBodies[idx++] = {planet.position, planet.velocity, planet.mass, planet.planetColor};
	}
	
	// Initialize Moons (Just for Earth and Jupiter, simplified)
	// Earth's Moon
	if (idx < numBodies) {
		simBodies[idx++] = {Vec3D(1.496e11 + 3.844e8, 0, 0), Vec3D(0, 29800 + 1022, 0), 7.34767309e22, ofColor(246,241,213)};
	}
	// Jupiter's Galilean Moons (Io, Europa, Ganymede, Callisto)
	double jupiter_x = 7.783e11;
	double jupiter_v = 13100;
	std::vector<PlanetInit> moonData = {
		{8.9319e22, Vec3D(jupiter_x + 4.22e8), Vec3D( jupiter_v + 17334), ofColor(214,232,101)},  // Io
		{4.7998e22, Vec3D(jupiter_x + 6.71e8), Vec3D( jupiter_v + 13740), ofColor(76,86,113)},  // Europa
		{1.4819e23, Vec3D(jupiter_x + 1.07e9), Vec3D( jupiter_v + 10880), ofColor(139,125,130)},  // Ganymede
		{1.0759e23, Vec3D(jupiter_x + 1.88e9), Vec3D( jupiter_v + 8204), ofColor(202,207,211)}  // Callisto
	};
	
	
	for (const auto& moon : moonData)
	{
		if (idx >= numBodies) break;
		simBodies[idx++] = {moon.position, moon.velocity, moon.mass, moon.planetColor};
	}
	
	// Initialize Asteroid Belt (Approximate positions)
	for (; idx < numBodies * 0.9 && idx < numBodies; ++idx)
	{
		double theta = 2 * PI * ((double) rand() / RAND_MAX);
		double r = 3e11 + ((double) rand() / RAND_MAX) * 1e11;  // Between Mars and Jupiter
		simBodies[idx] = {Vec3D(r * cos(theta), r * sin(theta), cos(theta)* sin(theta)), Vec3D(18000), 1e21, ofColor(246,241,213)};
	}
	
	// Initialize Oort Cloud (randomly in a sphere)
	for (; idx < numBodies; ++idx)
	{
		double theta = 2 * PI * ((double) rand() / RAND_MAX);
		double phi = acos(2 * ((double) rand() / RAND_MAX) - 1);
		double r = 1e13 + ((double) rand() / RAND_MAX) * 2e13;  // Up to 3e13 m
		simBodies[idx] = {Vec3D(r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)), Vec3D(500), 2e13};
	}
	
	
	
	// Assigning radii based on actual mass
	// This will make sure that planets and moons are visible, while not affecting the physics.
	double minMass = std::numeric_limits<double>::max();
	double maxMass = std::numeric_limits<double>::min();
	
	for (size_t i = 0; i < idx; ++i)
	{
		minMass = std::min(minMass, simBodies[i].mass);
		maxMass = std::max(maxMass, simBodies[i].mass);
	}
	
	// Set displayRadius based on a logarithmic scale
	for (size_t i = 0; i < idx; ++i)
	{
		if(i<=LARGE_BODY_THRESHOLD)
		{
			simBodies[i].displayRadius= DISPLAY_RADIUS_BASE + DISPLAY_RADIUS_SCALE_LARGE * log2(simBodies[i].mass / minMass);
		}
		else if (i>LARGE_BODY_THRESHOLD)
		{
			simBodies[i].displayRadius = DISPLAY_RADIUS_BASE + DISPLAY_RADIUS_SCALE_SMALL * log2(simBodies[i].mass / minMass); //log10
		}
	}
	
	
	cout<<"\n\nSolar System Setup ofGetElapsedTimef(): " << ofGetElapsedTimef();
}
//*/

//  We'll define some orbital elements for the main planets at a particular epoch.
//  (These are approximate or example values; use real data for best results.)
struct PlanetOrbitalElements {
	const char* name;
	double e;      // eccentricity
	double a_m;    // semi-major axis in meters (or in AU if you want to incorporate the conversion internally)
	double i_deg;  // inclination in degrees
	double Omega_deg; // longitude of ascending node
	double w_deg;     // argument of perihelion
	double M_deg;     // mean anomaly at epoch
	double mass;      // in kg
	ofColor color;
};







static inline void SetupSolarSystem(Body* &simBodies, size_t numBodies)
{
	size_t idx = 0;
	size_t lIdx = 0, lCount = 0;
	size_t mIdx = 0, mCount = 0;
	
	//
	/*
	 // Suppose you want to parse these Horizons result files:
	 const char* horizonFiles[] =
	 {
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_sun.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_mercury.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_venus.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_earth.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_mars.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_jupiter.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_saturn.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_uranus.txt",
	 "/Users/98dav/Downloads/of_v0.12.0_osx_release/apps/myApps/Serial Hashed-Octree Barnes-Hut N-Body Simulator 0_3/horizon_results/horizons_results_neptune.txt"
	 };
	 int numFiles = sizeof(horizonFiles)/sizeof(horizonFiles[0]);
	 for(int i = 0; i < numFiles; i++)
	 {
	 int recCount = 0;
	 //HorizonsRecord* records = parseHorizonsFile(horizonFiles[i], &recCount);
	 SolarBody* horizonSolarBodies = initializeSolarSystemFromCSV(horizonFiles[i], &recCount);
	 //*/
	
	
	
	// 1) The Sun at origin with zero velocity
	if(idx < numBodies)
	{
		simBodies[idx] =
		{
			Vec3D(0,0,0),
			Vec3D(0,0,0),
			SOLAR_MASS,  // ~1.989e30
			ofColor(253,184,19)
		};
	}
	idx++;
	// 2) Define orbital elements for planets (approximate placeholders!)
	//    Semi-major axis is in meters. (1 AU ~ 1.496e11 m.)
	//    The angles i, Omega, w, M are in degrees.
	//    You can fill in actual data from e.g. NASA JPL Horizons for a particular epoch.
	PlanetOrbitalElements planetData[] = {
		{
			"Mercury",
			0.205630,          // e
			0.387098 * AU, // a
			7.00487,            // i
			48.331,             // Omega
			29.124,             // w
			174.795,            // M
			3.3011e23,
			ofColor(205,203,202)
		},
		{
			"Venus",
			0.006772,
			0.723332 * AU,
			3.39471,
			76.6807,
			54.8523,
			50.416,   // mean anomaly
			4.8675e24,
			ofColor(238,208,83)
		},
		{
			"Earth",
			0.016711,
			1.000000 * AU,
			0.00005,
			348.73936,
			114.20783,
			357.51716,
			5.97237e24,
			ofColor(52,198,144)
		},
		{
			"Mars",
			0.093394,
			1.523679 * AU,
			1.85061,
			49.57854,
			286.502,
			19.3564,
			6.4171e23,
			ofColor(173,98,66)
		},
		{
			"Jupiter",
			0.0489,
			5.2026 * AU,
			1.304,
			100.464,
			273.867,
			20.02, // example M
			1.8982e27,
			ofColor(225,225,226)
		},
		{
			"Saturn",
			0.0565,
			9.5549 * AU,
			2.485,
			113.665,
			339.392,
			317.02, // example
			5.6834e26,
			ofColor(250,229,191)
		},
		{
			"Uranus",
			0.046381,
			19.2184 * AU,
			0.773,
			74.006,
			96.998857,
			142.2386,
			8.6810e25,
			ofColor(225,238,238)
		},
		{
			"Neptune",
			0.008678,
			30.1104 * AU,
			1.767,
			131.784,
			273.187,
			48.12369,
			1.02413e26,
			ofColor(91,93,223)
		}
	};
	int numPlanets = ARRAY_SIZE(planetData);//sizeof(planetData)/sizeof(planetData[0]);
	printf("\n\n\nNumber of planets: %d\n", numPlanets);
	
	// 3) For each planet, convert orbital elements to Cartesian
	for(int p=0; p<numPlanets && idx<numBodies; p++)
	{
		// Make a fresh Body
		Body planet;
		planet.mass = planetData[p].mass;
		planet.bodyColor = planetData[p].color;
		
		// Use our new approach (from the snippet) to fill in .position, .velocity
		// This is: GetCartesianFromOrbitalDegrees( Body&, e, a, i_deg, Ω_deg, ω_deg, M_deg, mu )
		GetCartesianFromOrbitalDegrees(planet,planetData[p].e, planetData[p].a_m, planetData[p].i_deg, planetData[p].Omega_deg, planetData[p].w_deg, planetData[p].M_deg, planetData[p].mass);
		
		simBodies[idx++] = planet;
	}
	
	// 4) (Optional) Moons or additional objects
	// For example, Earth's Moon:
	if(idx < numBodies)
	{
		Body moon;
		moon.mass = 7.34767309e22;
		moon.bodyColor = ofColor(246,241,213);
		
		// Example approximate orbital elements:
		double e_moon = 0.0549;
		double a_moon = 3.844e8; // ~384,400 km
		double i_deg_moon = 5.145;
		double Omega_deg_moon = 125.08;  // example
		double w_deg_moon = 318.15;      // argument of perigee
		double M_deg_moon = 135.27;      // example
										 // central mu for Earth if orbit is Earth-centered => G*(5.97237e24).
		double mu_earth = G * 5.97237e24;
		
		GetCartesianFromOrbitalDegrees(moon, e_moon, a_moon, i_deg_moon, Omega_deg_moon, w_deg_moon, M_deg_moon, mu_earth);
		
		// Shift that to Earth's heliocentric frame, or keep it Earth-centered?
		// If you want Earth + Moon in a purely heliocentric system, you must
		// add Earth's heliocentric position/velocity to the moon's Earth-centered position.
		// For simplicity, let's just do that:
		Body &earthRef = simBodies[3]; // assuming Earth is the 3rd body after Sun, Mercury, Venus
		moon.position += (earthRef.position * 1.1);  // slightly offset
		moon.velocity += earthRef.velocity;
		
		simBodies[idx++] = moon;
	}

	
	
	
	//
	/*
	// Initialize the Planets (masses in kg, positions and velocities in m and m/s)
	struct PlanetInit {
		double mass;
		Vec3D position;
		Vec3D velocity;
		ofColor planetColor;
	};
	// 5) Possibly do the same for Jupiter’s Galilean Moons or random asteroids, etc.
	// ...
	// Jupiter's Galilean Moons (Io, Europa, Ganymede, Callisto)
	double jupiter_x = 7.783e11;
	double jupiter_v = 13100;
	std::vector<PlanetInit> moonData = {
		{8.9319e22, Vec3D(jupiter_x + 4.22e8), Vec3D( jupiter_v + 17334), ofColor(214,232,101)},  // Io
		{4.7998e22, Vec3D(jupiter_x + 6.71e8), Vec3D( jupiter_v + 13740), ofColor(76,86,113)},  // Europa
		{1.4819e23, Vec3D(jupiter_x + 1.07e9), Vec3D( jupiter_v + 10880), ofColor(139,125,130)},  // Ganymede
		{1.0759e23, Vec3D(jupiter_x + 1.88e9), Vec3D( jupiter_v + 8204), ofColor(202,207,211)}  // Callisto
	};
	for (const auto& moon : moonData)
	{
		if (idx >= numBodies) break;
		simBodies[idx++] = {moon.position, moon.velocity, moon.mass, moon.planetColor};
	}
	//*/
	
	
	
	
	//
	/*
	// Initialize Asteroid Belt (Approximate positions)
	 for (; idx < numBodies * 0.9 && idx < numBodies; ++idx){
		double theta = 2 * PI * ((double) rand() / RAND_MAX);
		double r = 3e11 + ((double) rand() / RAND_MAX) * 1e11;  // Between Mars and Jupiter
		simBodies[idx] = {Vec3D(r * cos(theta), r * sin(theta), cos(theta)* sin(theta)), Vec3D(18000), 1e21, ofColor(246,241,213)};}
	// Initialize Oort Cloud (randomly in a sphere)
	 for (; idx < numBodies; ++idx){
		double theta = 2 * PI * ((double) rand() / RAND_MAX);
		double phi = acos(2 * ((double) rand() / RAND_MAX) - 1);
		double r = 1e13 + ((double) rand() / RAND_MAX) * 2e13;  // Up to 3e13 m
		simBodies[idx] = {Vec3D(r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)), Vec3D(500), 2e13};}//*/
	
	
	// Asteroid Belt initialization (improved)
	for (; idx < numBodies * 0.5 && idx < numBodies; ++idx)
	{
		// Realistic Orbital Elements for Asteroid Belt
		double e = ofRandom(0.0, 0.2);  // small eccentricities
		double a = ofRandom(2.2, 3.2) * AU; // ~2.2 to 3.2 AU
		double i_deg = ofRandom(0, 15); // low inclination
		double Omega_deg = ofRandom(0, 360);
		double omega_deg = ofRandom(0, 360);
		double M_deg = ofRandom(0, 360);
		double mass = ofRandom(1e15, 1e21); // typical asteroid mass
		
		double pos[3], vel[3];
		computeCartesianFromOrbital(a, e, toRadians(i_deg), toRadians(Omega_deg), toRadians(omega_deg), toRadians(M_deg), MU, pos, vel);
		
		
		
		simBodies[idx] = {
			Vec3D(pos[0], pos[1], pos[2]),
			Vec3D(vel[0], vel[1], vel[2]),
			mass,
			ofColor(246,241,213) // asteroid color
		};
	}
	
	// Oort Cloud initialization (improved)
	for (; idx < (numBodies-1+1); ++idx)
	{
		// Oort cloud orbits are highly elliptical and isotropic
		double e = ofRandom(0.8, 0.99);  // highly eccentric
		double a = ofRandom(1000, 100000) * AU; // 1,000 AU to 100,000 AU (~1e14m to 1.496e16 m)
		double i_deg = ofRandom(0, 180);  // fully isotropic
		double Omega_deg = ofRandom(0, 360);
		double omega_deg = ofRandom(0, 360);
		double M_deg = ofRandom(0, 360);
		double mass = ofRandom(1e12, 1e17);  // typical comet mass range
		
		double pos[3], vel[3];
		computeCartesianFromOrbital(a, e, toRadians(i_deg), toRadians(Omega_deg), toRadians(omega_deg), toRadians(M_deg), MU, pos, vel);
		
		simBodies[idx] = {
			Vec3D(pos[0], pos[1], pos[2]),
			Vec3D(vel[0], vel[1], vel[2]),
			mass,
			ofColor(180,180,240)  // Oort cloud visual distinction
		};
	}
	
	
	for(size_t i = 0; i <idx; i++)
	{
		if(simBodies[i].mass > LARGE_BODY_THRESHOLD)
		{
			lCount++;
			continue;
		}
		if(simBodies[i].mass > MEDIUM_BODY_THRESHOLD)
		{
			mCount++;
			continue;
		}
	}
	
	//instead of having sIDX, just take the difference between idx and lIdx+mIdx
	double *lIdxMasses = (double*)malloc(lCount * sizeof(double));
	double *mIdxMasses = (double*)malloc(mCount * sizeof(double));
	for(size_t i = 0; i <idx; i++)
	{
		if(simBodies[i].mass > LARGE_BODY_THRESHOLD)
		{
			//printf("\nLarge body %zu mass: %f\n", i, lIdxMasses[i]);
			lIdxMasses[lIdx++] = simBodies[i].mass;
			continue;
		}
		if(simBodies[i].mass > MEDIUM_BODY_THRESHOLD)
		{
			//printf("\nMedium body %zu mass: %f\n", i, mIdxMasses[i]);
			mIdxMasses[mIdx++] = simBodies[i].mass;
			continue;
		}
	}
	double lIdxMass = sum_elements(lIdxMasses, idx);
	double mIdxMass = sum_elements(mIdxMasses, mIdx);
	
	
	
	

	
	
	
	// 6) Now that we have assigned positions & velocities, let's assign displayRadius
	
	
	
	///*
	// Assigning radii based on actual mass DISPLAY_RADIUS_SCALE_MEDIUM SMALL_BODY_THRESHOLD
	double *bodyMasses = (double*)malloc(idx * sizeof(double));
	for(size_t i = 0; i < idx; ++i)
	{
		bodyMasses[i] = simBodies[i].mass;
		//printf("\nBody %zu mass: %f\n", i, bodyMasses[i]);
	}
	double totalMass = sum_elements(bodyMasses, idx);
	double minimumMass = min_element(bodyMasses, idx);
	double maximumMass = max_element(bodyMasses, idx);
	//printf("\n\n\nTotal mass: %f\n\nMinimum mass: %f\n\nMaximum mass: %f\n\n\n", totalMass, minimumMass, maximumMass);

	
	for(size_t iBody = 0; iBody < idx; iBody++)
	{
		if(simBodies[iBody].mass > LARGE_BODY_THRESHOLD)
		{
			//printf("\n\n\nLARGE_BODY_THRESHOLD i_%zu mass: %f\n", iBody, simBodies[iBody].mass);
			simBodies[iBody].displayRadius = DISPLAY_RADIUS_BASE +
			DISPLAY_RADIUS_SCALE_LARGE * log2(simBodies[iBody].mass / (lIdxMass-simBodies[iBody].mass));
			continue;
		}
		if(simBodies[iBody].mass > MEDIUM_BODY_THRESHOLD)
		{
			//printf("\n\n\nMEDIUM_BODY_THRESHOLD i_%zu mass: %f\n", iBody, simBodies[iBody].mass);
			simBodies[iBody].displayRadius = DISPLAY_RADIUS_BASE +
			DISPLAY_RADIUS_SCALE_MEDIUM * log2(simBodies[iBody].mass / (mIdxMass-simBodies[iBody].mass));
			continue;
		}
		else if(simBodies[iBody].mass > SMALL_BODY_THRESHOLD)
		{
			//printf("\n\n\nSMALL_BODY_THRESHOLD i_%zu mass: %f\n", iBody, simBodies[iBody].mass);
			simBodies[iBody].displayRadius = DISPLAY_RADIUS_BASE +
			DISPLAY_RADIUS_SCALE_SMALL * log2(simBodies[iBody].mass / ((totalMass-(mIdxMass+lIdxMass))-simBodies[iBody].mass));
			continue;
		}
		else
		{
			simBodies[iBody].displayRadius = DISPLAY_RADIUS_BASE +
			(DISPLAY_RADIUS_SCALE_SMALL * 0.1) * log2(simBodies[iBody].mass / totalMass);;
		}
	}
	//*/
	
	
	
	//i_0 mass: 1991668748881006311016377614336.000000
	//M mass: 1988999999999999901909255192576.000000
	//Minimum mass: 3015268630528.000000
	

	
	//
	/*
	double minMass = totalMass;//std::numeric_limits<double>::max();
	double maxMass = std::numeric_limits<double>::min();
	printf("\n\n\nMinimum mass: %f\n", minMass);
	printf("\nMaximum mass: %f\n\n\n\n\n\n\n", maxMass);
	for(size_t iBody = 0; iBody < idx; iBody++)
	{
		minMass = minimum(minMass, simBodies[iBody].mass); // std::min(
		maxMass = maximum(maxMass, simBodies[iBody].mass); // std::max(
	}
	for(size_t iBody = 0; iBody < idx; iBody++)
	{
		if(iBody <= LARGE_BODY_THRESHOLD)
		{
			simBodies[iBody].displayRadius = DISPLAY_RADIUS_BASE +
			DISPLAY_RADIUS_SCALE_LARGE * log2(simBodies[iBody].mass / minMass);
		}
		else
		{
			simBodies[iBody].displayRadius = DISPLAY_RADIUS_BASE +
			DISPLAY_RADIUS_SCALE_SMALL * log2(simBodies[iBody].mass / minMass);
		}
	}
	//*/
	
	
	

	
	// Debug message
	cout << "\n\nSolar System Setup using orbital elements. Bodies assigned: " << idx
	<< "\nCurrent ofGetElapsedTimef(): " << ofGetElapsedTimef() << endl;
}







/*
 toRadians
 Converts an angle from degrees to radians. This is a utility function that facilitates the conversion of angular measurements, which is essential for calculations that require angle inputs in radians.
 
 Parameters:
 - degrees: The angle in degrees to be converted.
 
 Returns:
 - The angle in radians.
 */
double toRadians(double degrees)
{
	/// Convert degrees to radians using the conversion factor π/180.
	/// RAD = DEG * π / 180
	/// 					1/180 = 0.0055555555556
	///						π*0.0055555555556 = 0.01745329252
	///	=> RAD = DEG * 0.01745329252
	return (degrees * 0.01745329252); // (M_PI *0.0055555555556);
}

/*
 calculatePositionAndVelocity
 Calculates the position and velocity of a celestial body in Cartesian coordinates based on its orbital parameters. The function performs coordinate transformations from the orbital plane to the inertial frame, considering angular momentum and velocity transformation.
 
 Parameters:
 - body: The celestial body whose position and velocity are being calculated. This is updated in-place.
 - r: The radial distance from the center of attraction.
 - ν: The true anomaly, representing the angle between the direction of periapsis and the current position of the body.
 - i: The inclination of the orbital plane.
 - Ω: The longitude of the ascending node.
 - ω: The argument of periapsis.
 - mu: The standard gravitational parameter.
 - a: The semi-major axis of the orbit.
 - e: The eccentricity of the orbit.
 
 Returns:
 - None. The function updates the 'body' object directly.
 */
void calculatePositionAndVelocity(Body &body, double r, double ν, double i, double Ω, double ω, double mu, double a, double e)
{
	/// Calculate position in the orbital plane using polar coordinates.
	double x_op = r * cos(ν);
	double y_op = r * sin(ν);
	
	/// Transform coordinates to the inertial frame using rotation matrices.
	body.position.x = x_op * (cos(Ω) * cos(ω) - sin(Ω) * sin(ω) * cos(i)) -
	y_op * (cos(Ω) * sin(ω) + sin(Ω) * cos(ω) * cos(i));
	body.position.y = x_op * (sin(Ω) * cos(ω) + cos(Ω) * sin(ω) * cos(i)) +
	y_op * (cos(Ω) * cos(ω) * cos(i) - sin(Ω) * sin(ω));
	body.position.z = x_op * (sin(ω) * sin(i)) + y_op * (cos(ω) * sin(i));
	
	/// Compute angular momentum and velocity components in the orbital plane.
	double h = sqrt(mu * a * (1 - e * e));
	double vx_op = -mu / h * sin(ν);
	double vy_op = mu / h * (e + cos(ν));
	
	/// Transform velocity components to the inertial frame.
	body.velocity.x = vx_op * (cos(Ω) * cos(ω) - sin(Ω) * sin(ω) * cos(i)) -
	vy_op * (cos(Ω) * sin(ω) + sin(Ω) * cos(ω) * cos(i));
	body.velocity.y = vx_op * (sin(Ω) * cos(ω) + cos(Ω) * sin(ω) * cos(i)) +
	vy_op * (cos(Ω) * cos(ω) * cos(i) - sin(Ω) * sin(ω));
	body.velocity.z = vx_op * (sin(ω) * sin(i)) + vy_op * (cos(ω) * sin(i));
}

/*
 GetCartesianFromOrbitalRadians
 Converts orbital elements expressed in radians into Cartesian coordinates for a celestial body. It uses the Newton-Raphson method to solve for the eccentric anomaly and then calculates the true anomaly, radial distance, and updates the body's position and velocity.
 
 Parameters:
 - body: The celestial body whose Cartesian coordinates are being calculated.
 - e: The eccentricity of the orbit.
 - a: The semi-major axis of the orbit.
 - i: The inclination of the orbital plane (radians).
 - Ω: The longitude of the ascending node (radians).
 - ω: The argument of periapsis (radians).
 - M: The mean anomaly.
 - mu: The standard gravitational parameter.
 
 Returns:
 - None. The function updates the 'body' object directly.
 */
static inline void GetCartesianFromOrbitalRadians(Body &body, double e, double a, double i, double Ω, double ω, double M, double mu)
{
	const double tolerance = 1e-10; // Consistent tolerance for Newton-Raphson convergence
	double E = M;
	double delta;
	int iter = 0, max_iter = 100;
	
	/// Use Newton-Raphson method to iteratively solve for the eccentric anomaly (E).
	do {
		double f = E - e * sin(E) - M;   // Kepler's equation
		double f_prime = 1 - e * cos(E); // Derivative of Kepler's equation
		delta = -f / f_prime;            // Compute the correction
		E += delta;                      // Update the estimate of E
		iter++;
	} while (fabs(delta) > tolerance && iter < max_iter);
	
	/// Calculate the true anomaly ν and radial distance r.
	double ν = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
	double r = a * (1 - e * cos(E));
	
	/// Update the body's position and velocity in Cartesian coordinates.
	calculatePositionAndVelocity(body, r, ν, i, Ω, ω, mu, a, e);
}

/*
 GetCartesianFromOrbitalDegrees
 Converts orbital elements expressed in degrees into Cartesian coordinates for a celestial body. This function first converts angle measurements from degrees to radians and then calls GetCartesianFromOrbitalRadians to perform the conversion.
 
 Parameters:
 - body: The celestial body whose Cartesian coordinates are being calculated.
 - e: The eccentricity of the orbit.
 - a: The semi-major axis of the orbit.
 - i_deg: The inclination of the orbital plane (degrees).
 - Ω_deg: The longitude of the ascending node (degrees).
 - ω_deg: The argument of periapsis (degrees).
 - M: The mean anomaly.
 - mu: The standard gravitational parameter.
 
 Returns:
 - None. The function updates the 'body' object directly.
 */
static inline void GetCartesianFromOrbitalDegrees(Body &body, double e, double a, double i_deg, double Ω_deg, double ω_deg, double M, double mu)
{
	/// Convert orbital angles from degrees to radians.
	double i_rad = toRadians(i_deg);
	double Ω_rad = toRadians(Ω_deg);
	double ω_rad = toRadians(ω_deg);
	
	/// Perform conversion using the radian-based function.
	GetCartesianFromOrbitalRadians(body, e, a, i_rad, Ω_rad, ω_rad, M, mu);
}









/*
 * // Function to initialize a Body object given its orbital elements
 * // Assumes that the central body is at the origin
 * // Based on https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
 */
static inline void initializeBody(Body &body, double eccentricity, double semiMajorAxis, double inclination, double longitudeOfAscendingNode, double argumentOfPerihelion, double meanAnomaly, double planetMass)
{
	// Error checking for orbital parameters
	if (semiMajorAxis <= 0)
	{
		throw std::invalid_argument("Semi-major axis must be positive.");
	}
	if (eccentricity < 0 || eccentricity >= 1)
	{
		throw std::invalid_argument("Eccentricity must be in the range [0, 1).");
	}
	
	
	// Calculate Mean Motion (n)
	double n = sqrt(G * (SOLAR_MASS + planetMass) / pow(semiMajorAxis, 3));
	
	// Calculate True Anomaly (nu) :  defines the position of the orbiting body along the ellipse at a specific time (the "epoch")
	double nu = meanAnomaly + 2 * eccentricity * sin(meanAnomaly) + (5.0 / 4) * pow(eccentricity, 2) * sin(2 * meanAnomaly);
	
	// Calculate r and f
	double r = semiMajorAxis * (1 - pow(eccentricity, 2)) / (1 + eccentricity * cos(nu));
	double f = nu + argumentOfPerihelion;
	
	
	// Convert to radians as C++ trigonometric functions expect arguments in radians
	double inclinationRad = ofDegToRad(inclination);
	double longitudeOfAscendingNodeRad = ofDegToRad(longitudeOfAscendingNode);
	double fRad = ofDegToRad(f);
	
	
	// Convert to Cartesian coordinates
	body.position.x = r * (cos(longitudeOfAscendingNodeRad) * cos(fRad) - sin(longitudeOfAscendingNodeRad) * sin(fRad) * cos(inclinationRad));
	body.position.y = r * (sin(longitudeOfAscendingNodeRad) * cos(fRad) + cos(longitudeOfAscendingNodeRad) * sin(fRad) * cos(inclinationRad));
	body.position.z = r * (sin(inclinationRad) * sin(fRad));
	
	// Calculate velocities
	double h = n * semiMajorAxis * sqrt(1 - pow(eccentricity, 2));
	body.velocity.x = -sin(f) * h / r;
	body.velocity.y = (cos(f) + eccentricity) * h / r;
	body.velocity.z = 0; // Assuming orbits are planar for simplicity
}






static inline void GetCartesianFromOrbital(Body &body, double eccentricity, double semiMajorAxis, double inclination, double longitudeOfAscendingNode, double argumentOfPerihelion, double meanAnomaly, double mu)
{
	/*----- 1. Calculate eccentric anomaly using Newton's method, E from the mean anomaly and eccentricity -----*/
	double E = meanAnomaly;  // Initial guess
	for (int k = 0; k < 100; ++k)
	{
		E = E - (E - eccentricity * sin(E) - meanAnomaly) / (1 - eccentricity * cos(E));
	}
	
	
	
	
	
	// Eccentric anomaly calculation with precision threshold
	const double precision = 1e-6; // Example threshold
	double E_prev;
	do {
		E_prev = E;
		E = E - (E - eccentricity * sin(E) - meanAnomaly) / (1 - eccentricity * cos(E));
	} while (std::abs(E - E_prev) > precision);
	
	
	
	
	
	
	
	/*----- 2. Compute Cartesian coordinates in a perifocal coordinate system -----*/
	// Calculate position in the orbital plane
	double xp = semiMajorAxis * (cos(E) - eccentricity);
	double yp = semiMajorAxis * sin(E) * sqrt(1 - eccentricity * eccentricity);
	
	// Calculate coordinates in 3D space
	double x = xp * (cos(longitudeOfAscendingNode) * cos(argumentOfPerihelion) - sin(longitudeOfAscendingNode) * sin(argumentOfPerihelion) * cos(inclination)) -
	yp * (cos(longitudeOfAscendingNode) * sin(argumentOfPerihelion) + sin(longitudeOfAscendingNode) * cos(argumentOfPerihelion) * cos(inclination));
	double y = xp * (sin(longitudeOfAscendingNode) * cos(argumentOfPerihelion) + cos(longitudeOfAscendingNode) * sin(argumentOfPerihelion) * cos(inclination)) +
	yp * (sin(longitudeOfAscendingNode) * sin(argumentOfPerihelion) - cos(longitudeOfAscendingNode) * cos(argumentOfPerihelion) * cos(inclination));
	double z = xp * (sin(inclination) * sin(argumentOfPerihelion)) + yp * (sin(inclination) * cos(argumentOfPerihelion));
	
	
	body.position.x = x;
	body.position.y = y;
	body.position.z = z;
	
	
	/*----- 3. Transform these coordinates to the heliocentric inertial system -----*/
	// Calculate velocity in the orbital plane
	double Edot = sqrt(mu / (semiMajorAxis * semiMajorAxis * semiMajorAxis)) / (1 - eccentricity * cos(E));
	double vxp = -semiMajorAxis * Edot * sin(E);
	double vyp = semiMajorAxis * Edot * sqrt(1 - eccentricity * eccentricity) * cos(E);
	
	// Calculate velocity in 3D space
	double vx = vxp * (cos(longitudeOfAscendingNode) * cos(argumentOfPerihelion) - sin(longitudeOfAscendingNode) * sin(argumentOfPerihelion) * cos(inclination)) -
	vyp * (cos(longitudeOfAscendingNode) * sin(argumentOfPerihelion) + sin(longitudeOfAscendingNode) * cos(argumentOfPerihelion) * cos(inclination));
	double vy = vxp * (sin(longitudeOfAscendingNode) * cos(argumentOfPerihelion) + cos(longitudeOfAscendingNode) * sin(argumentOfPerihelion) * cos(inclination)) +
	vyp * (sin(longitudeOfAscendingNode) * sin(argumentOfPerihelion) - cos(longitudeOfAscendingNode) * cos(argumentOfPerihelion) * cos(inclination));
	double vz = vxp * (sin(inclination) * sin(argumentOfPerihelion)) + vyp * (sin(inclination) * cos(argumentOfPerihelion));
	
	body.velocity.x = vx;
	body.velocity.y = vy;
	body.velocity.z = vz;
	
}













