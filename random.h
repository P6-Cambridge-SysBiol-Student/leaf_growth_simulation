/*
 Basic diffusion and FRAP simulation
 Jonathan Ward and Francois Nedelec, Copyright EMBL 2007-2009
 */
#include <cstdlib>
/// signed random real in [-1, 1]
/// used to create random initial starting positions and velocities
float mySrand()
{
    const float scale = 2.0 / static_cast<float>(RAND_MAX);
    return static_cast<float>( rand() ) * scale - 1.0;
}

/// positive random real in [0, 1]
float myPrand()
{
    const float scale = 1.0 / ( 1+static_cast<float>(RAND_MAX) );
    return static_cast<float>( 1+rand() ) * scale;
}

