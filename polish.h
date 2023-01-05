//
// Created by finley on 05/01/23.
//

#ifndef FRAP_POLISH_H
#define FRAP_POLISH_H
#define WORD  unsigned int

#endif //FRAP_POLISH_H

/// hard-coded limit to the number of particles
/// size_t = size in bytes,  const means MAX is immutable

const size_t MAX = 16384;
const bool debugStatus = 1;

/// create an array of these and allocate them the max amount of memory possibly required
/// this is global and not in the main.c file to keep it tidy
/// need to initialise the triangleIndexList pointer before delaunay triangulation

Point pointsArray[MAX];
int numTriangleVertices = 0;
WORD* triangleIndexList;

/// window size in pixels
int winW = 800;
int winH = 800;


