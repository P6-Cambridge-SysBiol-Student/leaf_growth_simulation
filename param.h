/*
 Basic diffusion and FRAP simulation
 Jonathan Ward and Francois Nedelec, Copyright EMBL 2007-2009
 */

#include <cmath>
#include <sstream>
const double SCALING_FACTOR = 100000;

// physical parameters:  ensure to add any new parameters to the readOption() function
double xBound = 50 * SCALING_FACTOR;   /// half-width of box (X) in micrometers. If this value is lower than 100 the deulaunay triangulation misses points
double yBound = xBound;   /// half-height of box (Y), is set to be equal to y for saftey
double pixel = 1;    /// size of one pixel in GL units


int nbo = 101;    /// initial number of objects (points)

// spring-physics parameters
const double fluidViscosity = 0.0016; /// Pa.s, velocity of water at 20 degrees celcius
const double mobilityCoefficient = 6 * 3.14159 * fluidViscosity;
const double breakSpringCoeff = 1.5;

// timestep parameters
double timestep = 0.00006; /// viscosity is in Pa.sec so this is seconds. 60 fps means 1sec simulated = 1.8sec realtime
int delay = 16;         /// milli-seconds between successive display
double delta = 0.00001;
unsigned long seed = 2; /// seed for random number generator

// hormone parameters

double hormone1ProdRate = 5000;
double hormone1DegRate = 50; /// not used if reaction-diffusion used
double hormone1IntroTime = 0.00006;
vector2D hormone1OriginV1 = vector2D(0.1*xBound,0.1*xBound);
double hormone1DiffCoeff = 0.005 * SCALING_FACTOR;
double hormone1DiffPro = 0.6; /// this is made up i think,
double horm1Efficacy = 0.0000;

double hormone2ProdRate = 33333;
double hormone2DegRate = 666;
double hormone2DiffCoeff = 0.0005 *SCALING_FACTOR; /// in the gray-scott model the rate of diff of horm2 is twice 1
vector2D horm2Source1 = vector2D(0.2, 0);
vector2D horm2Srouce2 = vector2D(-0.2, -0.3);

double RDfeedRate = 0.039;
double RDkillRate = 0.053;
double reactRate1to2 = 1/timestep;
double lengthOfHorm2Prod = 30*timestep;


// mitosis parameters
double baseMaxProbOfDiv = 0.0001; /// should inform this with the maximal amount of cell division that can occur, dont think this should be tunable
double baseDesiredTotalCells = 100; /// used to calculate mitosis probabilites


double realTime = 0;     /// time in the simulated world
int finalIterationNumber = 100;  /// iterations before final frame



//-----------------------------------------------------------------------------

template <typename T>
int readParameter(const char arg[], const char name[], T & ptr)
{
    if ( arg == strstr(arg, name) ) {
        std::istringstream iss(arg+strlen(name));
        iss >> ptr;
        return !iss.fail();
    }
    return 0;
}

// TODO alter readOption so it can read input variables from the genetic algorithm
int readOption(const char arg[])
{
    if ( readParameter(arg, "n=",     nbo) )    return 1;
    if ( readParameter(arg, "seed=",  seed) )   return 1;
    if ( readParameter(arg, "delay=", delay) )  return 1;
    if ( readParameter(arg, "bounds=", xBound) )  return 1;
    //if ( readParameter(arg, "timestep=", &timestep) ) return 1;

    return 0;
}

