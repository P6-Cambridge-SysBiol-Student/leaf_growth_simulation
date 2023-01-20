/*
 Basic diffusion and FRAP simulation
 Jonathan Ward and Francois Nedelec, Copyright EMBL 2007-2009
 */

#include <cmath>
#include <sstream>

// physical parameters:  ensure to add any new parameters to the readOption() function
double xBound = 200000;   /// half-width of box (X). If this value is lower than 100 the deulaunay triangulation misses points
double yBound = xBound;   /// half-height of box (Y), is set to be equal to y for saftey
double pixel = 1;    /// size of one pixel in GL units

int nbo = 200;    /// number of particles

double delta = 0.00001;   /// time-step
double repulsionRadius = xBound/15;  /// radius of repulsion around a point
double dampening = 0.95;
double velocityNoiseParam = 10;
double fluidViscosity = 1.0016; /// mPa.s, velocity of water at 20 degrees celcius

double timestep = 0.000000000000030;
int delay = 16;         /// milli-seconds between successive display
unsigned long seed = 1; /// seed for random number generator

double realTime = 0;     /// time in the simulated world



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


int readOption(const char arg[])
{
    if ( readParameter(arg, "n=",     nbo) )    return 1;
    if ( readParameter(arg, "repulsionRadius=",  repulsionRadius) )   return 1;
    if ( readParameter(arg, "delta=", delta) )  return 1;
    if ( readParameter(arg, "seed=",  seed) )   return 1;
    if ( readParameter(arg, "delay=", delay) )  return 1;
    return 0;
}

