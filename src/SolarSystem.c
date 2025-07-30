//  DataExtraction.c
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02

#include "SolarSystem.h"
#include "CommonDefinitions.h"   // if you need e.g. MAX_NUM_FILE_LINES
#include "DataExtraction.h"
#include "StringUtilities.h"
#include "FileUtilities.h"
#include <ctype.h>


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
 *----------------------------------------------------------------------------*
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
//*/


/*----------------------------------------------------------------------------
 * 4) The main function to read the CSV and fill an array of SolarBody
 *
 *   CSV columns might be (assuming 8 columns):
 *      name, a_AU, e, i_deg, Omega_deg, omega_deg, M_deg, mass_kg
 *   or in meters, or in a different order, etc.  Adapt accordingly.
 *
 *   Returns a dynamically allocated array of `SolarBody`. On success, sets
 *   *outCount = the number of bodies, or 0 on error. Caller must free result.
 *----------------------------------------------------------------------------*/
SolarBody* initializeSolarSystemFromCSV(const char* csvFilePath, int* outCount)
{
	// 1) Read lines from file using your library
	int lineCount = count_file_lines(csvFilePath, MAX_NUM_FILE_LINES);
	
	char** fileLines = read_file_contents(csvFilePath, lineCount);
	if(!fileLines || lineCount < 2)
	{
		fprintf(stderr, "Error: No data or too few lines in %s\n", csvFilePath);
		*outCount = 0;
		return NULL;
	}
	
	// 2) Parse all lines into a fieldwise container: parse_entire_file
	int    fieldCount = 0;
	char*** separatedData = parse_entire_file(fileLines, lineCount, &fieldCount, ",");
	if(!separatedData || fieldCount == 0)
	{
		fprintf(stderr, "Error: Could not parse CSV file.\n");
		*outCount = 0;
		return NULL;
	}
	
	// For example, suppose the columns are:
	//   0: name
	//   1: semiMajorAxis (AU)
	//   2: eccentricity
	//   3: inclination (deg)
	//   4: longAscNode  (deg)
	//   5: argPerihelion(deg)
	//   6: meanAnomaly  (deg)
	//   7: mass_kg
	//
	// If your CSV is in a different order, just reorder the indices below.
	
	// 3) We'll allocate an array of SolarBody to hold all rows except the header
	int bodyCount = lineCount - 1; // skip the header row
	SolarBody* bodies = (SolarBody*) malloc(bodyCount * sizeof(SolarBody));
	if(!bodies)
	{
		fprintf(stderr, "Memory allocation failed for SolarBody array.\n");
		*outCount = 0;
		return NULL;
	}
	
	// 4) Loop over each row (line) from 1..(lineCount-1)
	//    and parse the relevant columns. Then convert to Cartesian.
	for(int row = 1; row < lineCount; row++)
	{
		// index in our local bodies[] array
		int idx = row - 1;
		// For safety, initialize to zero
		memset(&bodies[idx], 0, sizeof(SolarBody));
		
		// Extract each piece from separatedData[col][row].
		// Note col is fixed, row changes.
		// You might want to check for NULL or do strtod carefully.
		
		// name (col0)
		const char* nameField = separatedData[0][row];
		strncpy(bodies[idx].name, nameField, sizeof(bodies[idx].name)-1);
		
		// read a in AU. If your CSV is already in meters, skip the conversion factor.
		const char* aField = separatedData[1][row];
		double aAU = strtod(aField, NULL);
		double a_m = aAU * 1.495978707e11; // convert AU to meters
		bodies[idx].semiMajorAxis = a_m;
		
		// e
		const char* eField = separatedData[2][row];
		double eVal = strtod(eField, NULL);
		bodies[idx].eccentricity = eVal;
		
		// i in deg
		const char* iField = separatedData[3][row];
		double iDeg = strtod(iField, NULL);
		double iRad = iDeg * (M_PI/180.0);
		bodies[idx].inclination = iRad;
		
		// Omega in deg
		const char* OField = separatedData[4][row];
		double Odeg = strtod(OField, NULL);
		double Orad = Odeg*(M_PI/180.0);
		bodies[idx].lonAscending = Orad;
		
		// omega in deg
		const char* wField = separatedData[5][row];
		double wdeg = strtod(wField, NULL);
		double wrad = wdeg*(M_PI/180.0);
		bodies[idx].argPerihelion = wrad;
		
		// M in deg
		const char* MField = separatedData[6][row];
		double Mdeg = strtod(MField, NULL);
		double Mrad = Mdeg*(M_PI/180.0);
		bodies[idx].meanAnomaly = Mrad;
		
		// mass in kg
		const char* massField = separatedData[7][row];
		double massVal = strtod(massField, NULL);
		bodies[idx].mass = massVal;
		
		// 5) Convert to Cartesian relative to the Sun
		//    (This is a purely heliocentric approach.)
		computeCartesianFromOrbital(a_m, eVal, iRad, Orad, wrad, Mrad,
									muSun,
									bodies[idx].pos,
									bodies[idx].vel);
	}
	
	// 6) Cleanup parse_entire_file data
	//    You can free the separatedData, as you no longer need it.
	//    (Implementation depends on your library’s memory model.)
	// For example:
	for(int col=0; col < fieldCount; col++)
	{
		for(int r=0; r< lineCount; r++)
		{
			free(separatedData[col][r]);
		}
		free(separatedData[col]);
	}
	free(separatedData);
	
	// free fileLines as well
	for(int i=0; i< lineCount; i++)
	{
		free(fileLines[i]);
	}
	free(fileLines);
	
	// 7) Return final array to caller
	*outCount = bodyCount;
	return bodies;
}







/*
 parseHorizonsFile:
 Reads lines from the given file, scans them for ephemeris data lines of the
 form:
 
 2460744.500000000 = A.D. 2025-Mar-10 00:00:00.0000 TDB
 X =-7.8361449E+05 Y = ...
 VX= 1.26076e-02 VY= ...
 ...
 
 Returns an array of HorizonsRecord, storing each day’s record.
 *outCount is set to the number of records. The caller must free() the result.
 Scans lines for JD plus X=...,Y=...,Z=...,VX=..., etc.
 */
HorizonsRecord* parseHorizonsFile(const char* filepath, int* outCount)
{
	*outCount = 0;
	
	// 1) Read all lines from the file
	int lineCount = count_file_lines(filepath, MAX_NUM_FILE_LINES);
	
	char** lines = read_file_contents(filepath, lineCount);
	if(!lines || lineCount < 1) {
		fprintf(stderr,"File %s empty or error reading.\n", filepath);
		return NULL;
	}
	
	// We'll store records in a dynamic array (realloc as needed)
	int capacity = 32;
	HorizonsRecord* records = (HorizonsRecord*)malloc(capacity * sizeof(HorizonsRecord));
	if(!records) { perror("malloc"); return NULL; }
	
	int recordIndex = -1;  // We'll increment when we see a new JD line
	
	// 2) Loop over each line, see if it starts with something like 2460... TDB
	for(int i=0; i<lineCount; i++)
	{
		trim_string_whitespaces(lines[i]);
		if(strlen(lines[i])<1) continue; // skip empty
		
		// A typical data line starts with "2460744.500000000 = A.D. 2025-Mar-10 00:00..."
		// We'll detect that pattern:
		double potentialJD = 0.0;
		if(isdigit(lines[i][0]))
		{
			// parse the leading double
			// e.g. "2460744.500000000 = A.D. 2025-Mar-10 00:00:00.0000 TDB"
			char temp[256];
			strcpy(temp, lines[i]);
			char *pEquals = strstr(temp, "=");
			if(pEquals) {
				*pEquals = '\0'; // cut at '='
								 // left side now e.g "2460744.500000000 "
				potentialJD = atof(temp);
				
				// Create new record
				recordIndex++;
				if(recordIndex >= capacity) {
					capacity *= 2;
					records = (HorizonsRecord*)realloc(records, capacity*sizeof(HorizonsRecord));
					if(!records) { perror("realloc"); return NULL; }
				}
				memset(&records[recordIndex], 0, sizeof(HorizonsRecord));
				records[recordIndex].jdTDB = potentialJD;
				
				// Right side is the date/time string
				char *rightSide = pEquals+1;
				while(isspace(*rightSide)) rightSide++;
				strncpy(records[recordIndex].dateString, rightSide, 63);
				records[recordIndex].dateString[63] = '\0';
			}
		}
		else if(recordIndex>=0) {
			// We are within lines for the current record’s X, Y, Z, etc.
			// e.g. "X = ... Y = ... Z = ..." or "VX= ... VY= ... VZ= ..."
			// We can parse them by scanning for lines like X =someNumber
			// Typically there's "X =", "Y =", "Z =", "VX=", "VY=", "VZ="
			// in the next 2 lines.
			
			// Just find 'X =' pattern, parse next token
			// Then find 'Y =', parse next token, etc.
			
			char *xPtr = strstr(lines[i], "X =");
			if(xPtr) {
				// parse X
				xPtr += 3; // skip "X ="
				records[recordIndex].X = atof(xPtr);
				
				char *yPtr = strstr(lines[i], "Y =");
				if(yPtr) {
					yPtr += 3;
					records[recordIndex].Y = atof(yPtr);
				}
				char *zPtr = strstr(lines[i], "Z =");
				if(zPtr) {
					zPtr += 3;
					records[recordIndex].Z = atof(zPtr);
				}
			}
			
			char *vxPtr = strstr(lines[i], "VX=");
			if(vxPtr) {
				vxPtr += 3;
				records[recordIndex].VX = atof(vxPtr);
				
				char *vyPtr = strstr(lines[i], "VY=");
				if(vyPtr) {
					vyPtr += 3;
					records[recordIndex].VY = atof(vyPtr);
				}
				char *vzPtr = strstr(lines[i], "VZ=");
				if(vzPtr) {
					vzPtr += 3;
					records[recordIndex].VZ = atof(vzPtr);
				}
			}
		}
	}
	
	// 3) Final count = recordIndex+1 if we found any lines
	int totalFound = recordIndex+1;
	*outCount = totalFound;
	
	// Cleanup
	for(int i=0; i<lineCount; i++){
		free(lines[i]);
	}
	free(lines);
	
	if(totalFound<1) {
		free(records);
		return NULL;
	}
	return records;
}




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
EphemerisRecord* parseHorizonsFilez(const char* filePath, int* outCount)
{
	*outCount = 0; // Initialize output count
	
	// 1) Read all lines from the file
	int lineCount = count_file_lines(filePath, 100000); // Just pick a large max
	if(lineCount < 1){
		fprintf(stderr, "parseHorizonsFile: file is empty or unreadable: %s\n", filePath);
		return NULL;
	}
	
	char** fileLines = read_file_contents(filePath, lineCount);
	if(!fileLines){
		fprintf(stderr, "parseHorizonsFile: unable to read lines from %s\n", filePath);
		return NULL;
	}
	
	//
	/*
	char** fileLines = (char**)malloc(lineCount * sizeof(char*));
	// Process each string in the array
	for (int i = 0; i < lineCount; i++)
	{
		fileLines[i] = prune_string_whitespaces(fileContents[i]);
	}
	//*/
	
	// 2) Identify start (“$$SOE”) and end (“$$EOE”) of ephemeris
	int startIndex = -1;
	int endIndex   = -1;
	
	
	char *found = find_string_in_string_array(fileLines, lineCount, "$$SOE");
	printf("\n\nFound: %s\n", found);
	char *indexString = "$$SOE";
	for(int i=0; i<lineCount; i++)
	{
		char *occurence = determine_string_occurence(fileLines[i], indexString);
		//printf("\n\nOccurence: %s\n", occurence);
		//printf("fileLines[i]: %s\n", fileLines[i]);
		

		
		
		if(occurence != NULL && startIndex != -1)
		{
			endIndex = i;
			break;
		}
		if(occurence != NULL && endIndex == -1)
		{
			startIndex = i+1;
			indexString = "$$EOE";
			continue;
		}
		
		
		//if(strstr(fileLines[i], "$$SOE")){startIndex = i;continue;}
		//if(strstr(fileLines[i], "$$EOE")){endIndex = i;break;}
	}
	printf("\n\nstartIndex: %d\n", startIndex);
	printf("endIndex: %d\n", endIndex);
	
	
	
	
	if(startIndex == -1 || endIndex == -1 || endIndex <= startIndex){
		fprintf(stderr, "parseHorizonsFile: Could not find $$SOE ... $$EOE in %s\n", filePath);
		// Cleanup
		for(int i=0; i<lineCount; i++) free(fileLines[i]);
		free(fileLines);
		return NULL;
	}
	
	// We’ll parse lines between startIndex+1 and endIndex-1 inclusive
	int dataLinesCount = (endIndex - (startIndex)); // should be 124 for mars
	//printf("Data Lines Count: %d\n", dataLinesCount);
	if(dataLinesCount < 1){
		fprintf(stderr, "parseHorizonsFile: No data lines found between $$SOE and $$EOE.\n");
		for(int i=0; i<lineCount; i++) free(fileLines[i]);
		free(fileLines);
		return NULL;
	}
	
	// 3) We don’t know the exact # of ephemeris records a priori, but from NASA Horizons
	//    each record is ~3 lines:
	//       1) JD + date
	//       2) X=..., Y=..., Z=...
	//       3) VX=..., VY=..., VZ=...
	//    Then a blank or next JD line, etc.
	//    So the maximum possible records is dataLinesCount / 3 (roughly).
	EphemerisRecord* records = (EphemerisRecord*)calloc(dataLinesCount, sizeof(EphemerisRecord));
	if(!records){
		perror("parseHorizonsFile: calloc() failed for records");
		// Cleanup
		for(int i=0; i<lineCount; i++) free(fileLines[i]);
		free(fileLines);
		return NULL;
	}
	
	// 4) Parse line by line inside the data region
	int recIndex = 0;
	// We'll scan in blocks of lines. The pattern is typically:
	//   <JD line>
	//   X=..., Y=..., Z=...
	//   VX=..., VY=..., VZ=...
	//   (maybe LT=..., RG=..., RR=..., then next JD line)...
	
	// We'll do a simple approach: find lines with "X =", "Y =", "Z =", etc.
	// But first, parse JD from the line that has " = A.D."
	// Then parse the next lines for X, Y, Z, etc.
	
	for(int i = startIndex; i < endIndex; /* increment inside loop */)
	{
		// 4a) Find the line that has the JD before " = A.D."
		// Example: "2460744.500000000 = A.D. 2025-Mar-10 00:00:00.0000 TDB"
		if(strstr(fileLines[i], " = A.D.")){
			// Extract the JD from start of line up to ' '
			char jdString[64] = {0};
			// The JD is typically at the line start, so we can read until first space or '='
			// e.g. "2460744.500000000"
			// We can do a simpler approach: tokenize up to '='
			char* eqPtr = strstr(fileLines[i], "=");
			if(eqPtr){
				// eqPtr points to "=", so we copy from line start up to eqPtr
				// but eqPtr might have preceding spaces
				int length = (int)(eqPtr - fileLines[i]);
				if(length > 63) length = 63;
				strncpy(jdString, fileLines[i], length);
				jdString[length] = '\0';
				// Trim trailing spaces
				// Or we can do a simpler approach, ignoring them
				
				// Validate numeric
				if(string_is_numeric_with_alphanumeric(jdString)){
					records[recIndex].julianDate = atof(jdString);
				} else {
					fprintf(stderr, "Warning: JD not numeric: %s\n", jdString);
					records[recIndex].julianDate = 0.0;
				}
			}
			i++; // Next line should have X=..., Y=..., Z=...
		}
		else {
			// Not a JD line, skip
			i++;
			continue;
		}
		
		if(i >= endIndex) break;
		
		// 4b) Parse a line that has X=..., Y=..., Z=...
		// Example: " X =-1.467683777129898E+08 Y = 2.677311511827773E+07 Z = 2.452261530162208E+04"
		if( (strstr(fileLines[i], "X =") || strstr(fileLines[i], "X = ")  || strstr(fileLines[i], "X =-") || strstr(fileLines[i], "X=")) &&
		   (strstr(fileLines[i], "Y =") || strstr(fileLines[i], "Y = ")  || strstr(fileLines[i], "Y =-")  ||  strstr(fileLines[i], "Y=")) &&
		   (strstr(fileLines[i], "Z =") || strstr(fileLines[i], "Z = ")  || strstr(fileLines[i], "Z =-")  ||  strstr(fileLines[i], "Z=")) )
		{
			double xVal=0, yVal=0, zVal=0;
			// A simple approach: break line by spaces, find tokens after 'X=', 'Y=', 'Z='
			// We'll do naive parsing. Another approach is to find the substring after 'X=' up to 'Y='.
			
			// We can do multiple calls to e.g.:
			//   char* xPtr = strstr(fileLines[i], "X");
			//   then parse the numeric substring right after '='
			// But let's do the approach of tokenizing by spaces, scanning for something that starts with 'X='
			{
				// Tokenize
				// We have a function from StringUtilities, e.g. `tokenize_string(...)`.
				char* lineCopy = strdup(fileLines[i]);
				char* token = tokenize_string(lineCopy, " \t");
				while(token){
					// e.g. "X=-1.467683777129898E+08"
					if(token[0] == 'X'){ // or starts with 'X'
						char* eq = strchr(token, '=');
						if(eq){
							eq++; // eq points to the first digit
							if(string_is_numeric_with_alphanumeric(eq)){
								xVal = atof(eq);
							}
						}
					}
					else if(token[0] == 'Y'){
						char* eq = strchr(token, '=');
						if(eq){
							eq++;
							if(string_is_numeric_with_alphanumeric(eq)){
								yVal = atof(eq);
							}
						}
					}
					else if(token[0] == 'Z'){
						char* eq = strchr(token, '=');
						if(eq){
							eq++;
							if(string_is_numeric_with_alphanumeric(eq)){
								zVal = atof(eq);
							}
						}
					}
					token = tokenize_string(NULL, " \t");
				}
				free(lineCopy);
			}
			records[recIndex].position[0] = xVal;  // kilometers
			records[recIndex].position[1] = yVal;
			records[recIndex].position[2] = zVal;
			i++;
		} else {
			// If the line doesn't have X=..., skip or i++ to find it
			i++;
			continue;
		}
		
		if(i >= endIndex) break;
		
		// 4c) Parse next line for VX=..., VY=..., VZ=...
		// Example: " VX=-5.983816950304321E+00 VY=-2.938275245040573E+01 VZ= 2.077673118884960E-03"
		if( (strstr(fileLines[i], "VX=") || strstr(fileLines[i], "VX =") || strstr(fileLines[i], "VX = ")) &&
		   (strstr(fileLines[i], "VY=") || strstr(fileLines[i], "VY =") || strstr(fileLines[i], "VY = ")) &&
		   (strstr(fileLines[i], "VZ=") || strstr(fileLines[i], "VZ =") || strstr(fileLines[i], "VZ = ")) )
		{
			double vxVal=0, vyVal=0, vzVal=0;
			// same approach as above
			char* lineCopy = strdup(fileLines[i]);
			char* token = tokenize_string(lineCopy, " \t");
			while(token){
				if(strncmp(token, "VX=", 3)==0){
					char* eq = strchr(token, '=');
					if(eq){
						eq++;
						if(string_is_numeric_with_alphanumeric(eq)){
							vxVal = atof(eq);
						}
					}
				}
				else if(strncmp(token, "VY=", 3)==0){
					char* eq = strchr(token, '=');
					if(eq){
						eq++;
						if(string_is_numeric_with_alphanumeric(eq)){
							vyVal = atof(eq);
						}
					}
				}
				else if(strncmp(token, "VZ=", 3)==0){
					char* eq = strchr(token, '=');
					if(eq){
						eq++;
						if(string_is_numeric_with_alphanumeric(eq)){
							vzVal = atof(eq);
						}
					}
				}
				token = tokenize_string(NULL, " \t");
			}
			free(lineCopy);
			
			records[recIndex].velocity[0] = vxVal; // km/s
			records[recIndex].velocity[1] = vyVal;
			records[recIndex].velocity[2] = vzVal;
			
			// Successfully parsed one record
			recIndex++;
			i++;
		}
		else {
			// Not a velocity line
			i++;
			continue;
		}
		
		// Possibly skip the next line (the line with LT=..., RG=..., RR=...) if it exists
		// or we can parse them if we want, but not needed for EphemerisRecord
	}
	
	// 5) Set outCount to the number of valid records
	*outCount = recIndex;
	
	// 6) Cleanup the file lines
	for(int i=0; i<lineCount; i++){
		free(fileLines[i]);
	}
	free(fileLines);
	
	// 7) If no records found, free records array
	if(recIndex == 0){
		free(records);
		records = NULL;
	}
	
	return records;
}

