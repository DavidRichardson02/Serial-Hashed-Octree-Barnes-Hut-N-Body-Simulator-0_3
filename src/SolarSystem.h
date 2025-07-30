#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
//#include "ofMain.h" // or your color definitions
					// Possibly your "CommonDefinitions.h" if it’s needed







//static void computeCartesianFromOrbital(double a, double e, double iRad, double OmegaRad, double omegaRad, double MRad, double centralMu, double outPos[3], double outVel[3]);
/*----------------------------------------------------------------------------
 * internal function: computeCartesianFromOrbital
 *
 * Example function that similarly does E-solve + rotation. Shown here for
 * demonstration, but if you'd prefer to unify with your PhysicsHelpers
 * approach, you can do so.
 *
 * This function is used by 'initializeSolarSystemFromCSV'.
 *
 * param a,e,iRad,OmegaRad,omegaRad,MRad,centralMu => orbital elements
 * param outPos => array of 3 doubles for position
 * param outVel => array of 3 doubles for velocity
 */
/*----------------------------------------------------------------------------
 * 3) A helper that converts orbital elements -> Cartesian (heliocentric).
 *    This function follows the standard approach:
 *      1) Solve for eccentric anomaly E from mean anomaly M
 *      2) Determine position & velocity in the orbital plane
 *      3) Rotate by (Ω, i, ω) to get heliocentric inertial coords
 *----------------------------------------------------------------------------*/
static void computeCartesianFromOrbital(double a, double e, double iRad, double OmegaRad, double omegaRad, double MRad, double centralMu, double outPos[3], double outVel[3])
{
	// 1) Solve Kepler’s Equation for E via Newton’s method
	//    E - e sin(E) = M
	double E = MRad;
	for(int iter = 0; iter < 100; iter++)
	{
		double f  = E - e*sin(E) - MRad;
		double fp = 1.0 - e*cos(E);
		E = E - f/fp;
	}
	
	// 2) Perifocal coordinates (xp,yp) in orbital plane
	double cosE = cos(E);
	double sinE = sin(E);
	double sqrt1me2 = sqrt(1.0 - e*e);
	
	// radial distance from focal point
	double r = a*(1.0 - e*cosE);
	
	double xp = a*(cosE - e);
	double yp = a* sqrt1me2 * sinE;
	
	// Perifocal velocity:
	//   Edot = sqrt( mu / a^3 ) / (1 - e cos E )
	double n = sqrt(centralMu/(a*a*a));  // mean motion
	double Edot = n/(1.0 - e*cosE);
	double vxp = -a * Edot * sinE;
	double vyp =  a * Edot * sqrt1me2 * cosE;
	
	// 3) Rotate (xp,yp) and (vxp,vyp) from orbital plane to heliocentric coords
	//    The standard rotation sequence is by Ω around Z, then i around X, then ω around Z
	//    This can be done more explicitly or with a single combined matrix.
	//
	//    We’ll do Z(Ω) X(i) Z(ω). For position:
	//
	//    x = Rz(Ω) Rx(i) Rz(ω) * [xp, yp, 0]
	//
	//    For brevity, we’ll do a direct approach below:
	//
	double cosO = cos(OmegaRad), sinO = sin(OmegaRad);
	double cosi = cos(iRad),     sini = sin(iRad);
	double cosw = cos(omegaRad), sinw = sin(omegaRad);
	
	// Position (in 3D):
	// xp*(cosw cosO - sinw sinO cosi)
	// - yp*(sinw cosO + cosw sinO cosi)
	// etc. The fully expanded form is below:
	double x = (cosO*cosw - sinO*sinw*cosi)*xp + (-cosO*sinw - sinO*cosw*cosi)*yp;
	double y = (sinO*cosw + cosO*sinw*cosi)*xp + (-sinO*sinw + cosO*cosw*cosi)*yp;
	double z = (sinw*sini)*xp + (cosw*sini)*yp;
	
	// Velocity (similarly):
	double vx = (cosO*cosw - sinO*sinw*cosi)*vxp + (-cosO*sinw - sinO*cosw*cosi)*vyp;
	double vy = (sinO*cosw + cosO*sinw*cosi)*vxp + (-sinO*sinw + cosO*cosw*cosi)*vyp;
	double vz = (sinw*sini)*vxp + (cosw*sini)*vyp;
	
	outPos[0] = x;
	outPos[1] = y;
	outPos[2] = z;
	outVel[0] = vx;
	outVel[1] = vy;
	outVel[2] = vz;
}



/**
 * \file SolarSystem.h
 * \brief Declarations for reading orbital data from CSV, reading from Horizons files, and storing results in a struct.
 */

/**
 * A local struct to hold final “body” data after conversion to Cartesian.
 * You can adapt or rename this to match your actual simulation’s needs.
 */
typedef struct SolarBody
{
	char   name[128];
	double mass;          // in kilograms
	double pos[3];        // (x,y,z) in meters
	double vel[3];        // (vx,vy,vz) in m/s
	
	double semiMajorAxis; // a
	double eccentricity;   // e
	double inclination;    // i in radians
	double lonAscending;   // Ω in radians
	double argPerihelion;  // ω in radians
	double meanAnomaly;    // M in radians
	
} SolarBody;
/**
 * This is the gravitational parameter for the Sun if you want it in this file,
 * but typically that might be in a separate definitions header.
 */
static const double muSun = 6.67430e-11 * 1.989e30;

/**
 * \brief Reads a CSV file containing orbital elements (e, a, i, Ω, ω, M, mass, etc.)
 *        and returns an array of SolarBody with positions/velocities computed.
 *
 * CSV columns might be:
 *   name, a_AU, e, i_deg, Omega_deg, argPeri, M_deg, mass_kg
 *
 * \param csvFilePath The path to the CSV file
 * \param outCount    Will be set to the number of bodies read
 * \return A dynamically allocated array of SolarBody. Caller must free().
 */
SolarBody* initializeSolarSystemFromCSV(const char* csvFilePath, int* outCount);








/**
 * A struct for a single record line from a NASA JPL Horizons text output.
 */
typedef struct
{
	double jdTDB;          /* e.g. 2460744.5 */
	char   dateString[64]; /* e.g. "2025-Mar-10 00:00:00.0000 TDB" */
	double X, Y, Z;        /* position in KM or some units (based on Horizons output) */
	double VX, VY, VZ;     /* velocity */
} HorizonsRecord;

/**
 * \brief Reads a NASA JPL Horizons text file, scanning for lines containing
 *        Julian date, X, Y, Z, VX, VY, VZ, etc. Then returns them in an array.
 *
 * \param filepath The path to the Horizons text file
 * \param outCount pointer to an int that will store how many lines were found
 * \return A dynamically allocated array of HorizonsRecord. Caller must free().
 */
HorizonsRecord* parseHorizonsFile(const char* filepath, int* outCount);















// Our ephemeris struct:
typedef struct {
	double julianDate;    // e.g., 2460744.5
	double position[3];   // X, Y, Z in kilometers
	double velocity[3];   // VX, VY, VZ in km/s
} EphemerisRecord;

/**
 * parseHorizonsFile
 *
 * Parses a NASA JPL Horizons ephemeris file. Extracts each data record’s
 * Julian date (JDTDB), plus X, Y, Z, VX, VY, VZ. Data lines are delimited by
 * the “$$SOE” and “$$EOE” markers. The function returns a dynamically allocated
 * array of EphemerisRecord, and sets *outCount to the number of records found.
 *
 * Example usage:
 *   int recordCount = 0;
 *   EphemerisRecord* ephem = parseHorizonsFile("horizons_results_earth.txt", &recordCount);
 *   if(ephem && recordCount > 0){
 *       // do something with ephem...
 *   }
 *   free(ephem);
 *
 * @param filePath:   Path to the Horizons text file
 * @param outCount:   Will be set to the number of ephemeris records parsed
 * @return            Dynamically allocated array of EphemerisRecord; caller frees
 */
EphemerisRecord* parseHorizonsFilez(const char* filePath, int* outCount);



#endif // SOLAR_SYSTEM_H














