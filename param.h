/*
 Basic diffusion and FRAP simulation
 Jonathan Ward and Francois Nedelec, Copyright EMBL 2007-2009
 */
///tesst
#include <cmath>
#include <sstream>
#include <fstream>

const double SCALING_FACTOR = 100000;

// physical parameters:  ensure to add any new parameters to the readOption() function
double xBound = 1 * SCALING_FACTOR;   /// half-width of box (X) in micrometers
double yBound = xBound;   /// half-height of box (Y), is set to be equal to y for saftey
double pixel = 1;    /// size of one pixel in GL units


int nbo = 20;    /// initial number of objects (points)

// spring-physics parameters
const double fluidViscosity = 0.0016; /// Pa.s, velocity of water at 20 degrees celcius
const double mobilityCoefficient = 6 * 3.14159 * fluidViscosity;
const double breakSpringCoeff = 1.2;

// timestep parameters
double timestep = 0.00004; /// viscosity is in Pa.sec so this is seconds. 60 fps means 1sec simulated = 1.8sec realtime
int delay = 16;         /// milli-seconds between successive display
double delta = 0.00001;
unsigned long seed = 2; /// seed for random number generator
double finalTime = 0.05;
double realTime = 0;     /// time in the simulated world
int finalIterationNumber = 100;  /// iterations before final frame
int maxFourierCoeffs = 15;

// hormone parameters

double hormone1ProdRate = 100;
double hormone1DegRate = 10; /// not used if reaction-diffusion used
double hormone1IntroTime = 0.00360;
vector2D hormone1OriginV1 = vector2D(0.1 * xBound, 0.1 * xBound);
double inputHorm1DiffCoeff = 160;
double hormone1DiffCoeff = inputHorm1DiffCoeff * SCALING_FACTOR;
double horm1Efficacy = 5;
double horm1DivOrientVertComp = 5;
double horm1DivOrientHoriComp = 0;
vector2D horm1DivOrient = vector2D(horm1DivOrientHoriComp, horm1DivOrientVertComp);


double hormone2ProdRate = 33;
double hormone2DegRate = 666;
double horm1toHorm2Ratio = 0.9375;
double hormone2DiffCoeff =
        horm1toHorm2Ratio * hormone1DiffCoeff; /// in the gray-scott model the rate of diff of horm2 is twice 1
double horm2Efficacy = 10;
double hormone2IntroTime = 0.00360;
double horm2SourceHor = 0.2;
double horm2SourceVer = 0;
vector2D horm2Source1 = vector2D(horm2SourceHor, horm2SourceVer);
double horm2DivOrientVertComp = 0;
double horm2DivOrientHoriComp = 5;
vector2D horm2DivOrient = vector2D(horm2DivOrientHoriComp, horm2DivOrientVertComp);

double RDfeedRate = 35;
double RDfeedToKillRatio = 1.14;
double RDkillRate = RDfeedRate * RDfeedToKillRatio;
double reactRate1to2 = 6400;
double lengthOfHorm2Prod = 1;


// mitosis parameters
double baseMaxProbOfDiv = 0.0005; /// should inform this with the maximal amount of cell division that can occur, dont think this should be tunable
double DesiredTotalCells = 1000; /// used to calculate mitosis probabilites




bool displayInverseFourier = true;


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
    printf("[%s]\n", arg);
    if ( readParameter(arg, "n=",     nbo) )    return 1;
    if ( readParameter(arg, "inputHorm1DiffCoeff=",  inputHorm1DiffCoeff) )   return 1;
    if ( readParameter(arg, "horm1Efficacy=", horm1Efficacy) )  return 1;
    if ( readParameter(arg, "horm1DivOrientVertComp=", horm1DivOrientVertComp) )  return 1;
    if ( readParameter(arg, "horm1DivOrientHoriComp=", horm1DivOrientHoriComp) )  return 1;

    if ( readParameter(arg, "horm1toHorm2Ratio=",     horm1toHorm2Ratio) )    return 1;
    if ( readParameter(arg, "horm2Efficacy=",  horm2Efficacy) )   return 1;
    if ( readParameter(arg, "hormone2IntroTime=", hormone2IntroTime) )  return 1;
    if ( readParameter(arg, "horm2SourceHor=", horm2SourceHor) )  return 1;
    if ( readParameter(arg, "horm2SourceVer=", horm2SourceVer) )  return 1;
    if ( readParameter(arg, "horm2DivOrientVertComp=",     horm2DivOrientVertComp) )    return 1;
    if ( readParameter(arg, "horm2DivOrientHoriComp=",  horm2DivOrientHoriComp) )   return 1;
    if ( readParameter(arg, "lengthOfHorm2Prod=",  lengthOfHorm2Prod) )   return 1;

    if ( readParameter(arg, "RDfeedRate=", RDfeedRate) )  return 1;
    if ( readParameter(arg, "RDfeedToKillRatio=", RDfeedToKillRatio) )  return 1;
    if ( readParameter(arg, "reactRate1to2=", reactRate1to2) )  return 1;
    return 0;
}

void readFile(const char path[])
{
    std::string line;
    std::ifstream is(path);
    if ( !is.good() )
        printf("File `%s' cannot be read\n", path);
    while ( is.good() )
    {
        getline(is, line);
        readOption(line.c_str());
    }
}


