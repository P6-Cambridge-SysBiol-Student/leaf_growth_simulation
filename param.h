/*
 Basic diffusion and FRAP simulation
 Jonathan Ward and Francois Nedelec, Copyright EMBL 2007-2009
 */

#include <cmath>
#include <sstream>

// physical parameters:  ensure to add any new parameters to the readOption() function
double xBound = 200;   /// half-width of box (X). If this value is lower than 100 the deulaunay triangulation misses points
double yBound = 200;   /// half-height of box (Y). Same as above, keep the value above 100
double pixel = 1;    /// size of one pixel in GL units

int nbo = 50;    /// number of particles

double delta = 0.001;   /// time-step
double repulsionRadius = 10;  /// radius of repulsion around a point

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

