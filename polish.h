//
// Created by finley on 05/01/23.
//

#ifndef FRAP_POLISH_H
#define FRAP_POLISH_H
#define WORD  unsigned long

#endif //FRAP_POLISH_H

/// hard-coded limit to the number of particles
/// size_t = size in bytes,  const means MAX is immutable

const double PI = 3.14159265358979323846;
const size_t MAX = 16384;
const bool debugStatus = 1;


/// create an array of these and allocate them the max amount of memory possibly required
/// this is global and not in the main.c file to keep it tidy
/// need to initialise the triangleIndexList pointer before delaunay triangulation

Point pointsArray[MAX];
Hormone hormoneArray[numHormones];
int numTriangleVertices = 0;
WORD* triangleIndexList;
const int NAW = 50;  /// neighbourhood array width
double currentTime = 0;   /// a tracker for how many timesteps have passed

/// window size in pixels
int winW = 1000;
int winH = 1000;


static void error(int error, const char* text)
{
    fprintf(stderr, "GLFW Error: %s\n", text);
}

void limitNbo()
{
    // limit number of particles
    if ( nbo >= MAX ) nbo = MAX-1;
#if DEBUG
    printf("The number of points used is %d \n", nbo);
#endif
    //initialize random number generator
    srandom(seed);
}

