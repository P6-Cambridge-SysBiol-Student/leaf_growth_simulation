/*
Class which represents cells which repel and attract each other on a 2D plane
This file contains information about Object which contains
 */
#include <math.h>
/// TODO include a vector class and include vectors in object to clean code
/// for the compiler this doesn't slow down the programme

///Points class, contains info about xy displacement, velocity, and acceleration towards centre
class Point
{
public:  /// these are attributes that can be called outside of the script
    /// member variables:
    int color;
    double x, y; /// displacement
    double xvelocity, yvelocity; /// velocity
    double xSpringForce, ySpringForce;
    double extendedHooks, compressedHooks;  /// hooks constant for attracting points back to the centre
    double cellRadius;
    double cellMass;
    
    /// initialize each point in a random position with random x and y velocities
    /// currently these are set to start points randomly at the centre bottom to mimic plant leaves
    void reset()
    {
        x = (xBound * srand());
        xvelocity = 0.001 * srand();
        xSpringForce = 0;

        y = (yBound * srand());
        yvelocity = 0.001 * srand();
        ySpringForce = 0;

        extendedHooks   = 0.00003;
        compressedHooks = 0.003;
        cellRadius = 20; /// in micrometers
        cellMass = 1; /// in nanograms
        color = 1;
    }
    
    /// call initialize
    Point(){
        reset();
    }

    /// particles bounce off walls
    void bounce(){
        if ( x >  xBound )  x =  2*xBound - x, xvelocity = -xvelocity;
        if ( x < -xBound )  x = -2*xBound - x, xvelocity = -xvelocity;
        if ( y >  yBound )  y =  2*yBound - y, yvelocity = -yvelocity;
        if ( y < -yBound )  y = -2*yBound - y, yvelocity = -yvelocity;
    }


    //TODO find youngs modulus of the leaf (stiffness per unit area of the cross-section)

    /// make a step in the given direction
    void step(){
        /// change the velocity depending on the acceleration
        x += timestep * xSpringForce/(mobilityCoefficient * cellRadius/100000);
        y += timestep * ySpringForce/(mobilityCoefficient * cellRadius/100000);
    }

    /// partial display: this needs to be called between glBegin() and glEnd()
    void displayYellow(){
        /// transparency used to visualize overlapping particles
        if ( color == 1 )
            glColor4f(1.0, 1.0, 0.0, 0.5);
        else
            glColor4f(0.3, 0.3, 0.3, 0.5);
        glVertex2f(x, y);
    }

    void displayWhite()
    {
        /// transparency used to visualize overlapping particles
        if ( color == 1 )
            glColor4f(1.0, 1.0, 1.0, 0.5);
        else
            glColor4f(0.3, 0.3, 0.3, 0.5);
        glVertex2f(x, y);
    }
};
