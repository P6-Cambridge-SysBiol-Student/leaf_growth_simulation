/*
 Basic diffusion and FRAP simulation
 Jonathan Ward and Francois Nedelec, Copyright EMBL 2007-2009
 */

#include <cmath>
#include <sstream>

// physical parameters:
double xBound = 5;   // half-width of box (X)
double yBound = 5;   // half-height of box (Y)
double pixel = 1;    // size of one pixel in GL units

int nbo = 50;    // number of particles

double diff  = 1;       // diffusion constant
double delta = 0.001;   // time-step
double range = 1;       // radius of bleached zone

int delay = 16;         // milli-seconds between successive display
unsigned long seed = 1; // seed for random number generator

//derived quantities:
double alpha    = 0;     // diffusive motion within interval dt
double realTime = 0;     // time in the simulated world



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
    if ( readParameter(arg, "range=", range) )  return 1;
    if ( readParameter(arg, "n=",     nbo) )    return 1;
    if ( readParameter(arg, "diff=",  diff) )   return 1;
    if ( readParameter(arg, "delta=", delta) )  return 1;
    if ( readParameter(arg, "seed=",  seed) )   return 1;
    if ( readParameter(arg, "delay=", delay) )  return 1;
    return 0;
}

