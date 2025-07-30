#include "CommonDefinitions.h"
//  CommonDefinitions.c
//  ECE-370_Standardized_CSV_Data_Analysis
//  DavidRichardson02


#include "CommonDefinitions.h"






/**
 * Months Array:
 * Represents the names of the twelve months of the Gregorian calendar.
 * This constant array is used throughout the program wherever month names are required,
 * such as in date formatting or parsing tasks.
 */
const char *months[12] =
{
	"January", "February", "March",
	"April", "May", "June",
	"July", "August", "September",
	"October", "November", "December"
};

/**
 * Weekdays Array:
 * Contains the names of the seven days of the week.
 * This array is utilized in functions that involve day-of-week calculations, display, or parsing.
 */
const char *weekDays[7] =
{
	"Sunday", "Monday", "Tuesday", "Wednesday",
	"Thursday", "Friday", "Saturday"
};


/**
 * Common DateTime Formats Array:
 * A collection of commonly used date/time string formats.
 * These formats are used for detecting and formatting date/time fields in data set contents.
 * The array includes variations for 12-hour and 24-hour time formats, accommodating different
 * date ordering conventions (e.g., MM/DD/YYYY, YYYY/MM/DD, DD/MM/YYYY).
 *
 */
const char *commonDateTimeFormats[12] =
{
	// 12-hour time formats with different separators and precision levels
	"%m/%d/%Y %I:%M:%S %p", // MM/DD/YYYY HH:MM:SS AM/PM
	"%m/%d/%Y %I:%M %p",    // MM/DD/YYYY HH:MM AM/PM
	"%m-%d-%Y %I:%M:%S %p", // MM-DD-YYYY HH:MM:SS AM/PM
	"%m-%d-%Y %I:%M %p",    // MM-DD-YYYY HH:MM AM/PM
	
	// 24-hour time formats with Year-Month-Day ordering
	"%Y/%m/%d %H:%M:%S",    // YYYY/MM/DD HH:MM:SS
	"%Y/%m/%d %H:%M",       // YYYY/MM/DD HH:MM
	"%Y-%m-%d %H:%M:%S",    // YYYY-MM-DD HH:MM:SS
	"%Y-%m-%d %H:%M",       // YYYY-MM-DD HH:MM
	
	// 24-hour time formats with Day-Month-Year ordering
	"%d/%m/%Y %H:%M:%S",    // DD/MM/YYYY HH:MM:SS
	"%d/%m/%Y %H:%M",       // DD/MM/YYYY HH:MM
	"%d-%m-%Y %H:%M:%S",    // DD-MM-YYYY HH:MM:SS
	"%d-%m-%Y %H:%M"        // DD-MM-YYYY HH:MM
};








// ------------------ Global Constants and Macros Implementations -------------------
/// \{
const char* n_array_element_to_string(char** stringArray, int n, char* buffer)
{
	if (stringArray == NULL || buffer == NULL || BUFFER_SIZE == 0 || stringArray[n] == NULL){ return NULL; } // Validate input parameters
	snprintf(buffer, BUFFER_SIZE, "%s", stringArray[n]); // Copy the n-th element of the string array to the buffer
	return buffer; // Return the buffer containing the string
}
/// \}








/**
 * Common Unit Formats Array:
 * A collection of commonly used unit formats for representing physical quantities in data sets.
 * These unit formats cover a variety of scientific and engineering measurements.
 *
 */
const char *commonUnitFormats[] =
{
	// Basic units
	"m",    // meters
	"kg",   // kilograms
	"s",    // seconds
	"A",    // ampere
	"K",    // kelvin
	"mol",  // mole
	"cd",   // candela
	
	// Derived units with SI prefixes
	"mm",   // millimeters
	"cm",   // centimeters
	"km",   // kilometers
	"ms",   // milliseconds
	"µs",   // microseconds
	"ns",   // nanoseconds
	"MHz",  // Megahertz
	"GHz",  // Gigahertz
	"kJ",   // kilojoules
	"mW",   // milliwatts
	"kW",   // kilowatts
	
	// Compound units
	"m/s",      // meters per second
	"m/s^2",    // meters per second squared
	"kg/m^3",   // kilograms per cubic meter
	"W/m^2",    // watts per square meter
	"A/m^2",    // amperes per square meter
	"mol/m^3",  // moles per cubic meter
	"cd/m^2",   // candelas per square meter
	
	// Specialized units
	"Ohm",      // ohm
	"Pa",       // Pascal
	"N",        // Newton
	"J",        // Joule
	"Hz",       // Hertz
	"W",        // Watt
	"V",        // Volt
	"F",        // Farad
	"C",        // Coulomb
	"T",        // Tesla
	"H",        // Henry
	"lx",       // Lux
	"Bq",       // Becquerel
	"Gy",       // Gray
	"Sv",       // Sievert
	"kat",      // Katal
	
	// Units with per time variations
	"km/h",     // kilometers per hour
	"mph",      // miles per hour
	"g/cm^3",   // grams per cubic centimeter
	"l/min",    // liters per minute
};

