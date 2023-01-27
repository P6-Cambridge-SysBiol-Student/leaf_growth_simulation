/*
Working simulation of points orbiting around the centre of the window as a test for later physics
This file contains information about Object which contains
 */
#include <math.h>

///Points class, contains info about xy displacement, velocity, and acceleration towards centre
class Point
{
public:  /// these are attributes that can be called outside of the script
    /// member variables:
    int color;
    double x, y; /// displacement
    double xvelocity, yvelocity; /// velocity
    double xacceleration, yacceleration; /// acceleration
    double xSpringForce, ySpringForce;
    double xStokesDrag, yStokesDrag;
    double extendedHooks, compressedHooks;  /// hooks constant for attracting points back to the centre
    double cellRadius;
    double cellMass;
    
    /// initialize each point in a random position with random x and y velocities
    /// currently these are set to start points randomly at the centre bottom to mimic plant leaves
    void reset()
    {
        x = (xBound * srand());
        xvelocity = 0.001 * srand();
        xacceleration = 0;
        xSpringForce = 0;
        xStokesDrag = 0;

        y = (yBound * srand());
        xvelocity = 0.001 * srand();
        yacceleration = 0;
        ySpringForce = 0;
        yStokesDrag = 0;

        extendedHooks   = 0.000003;
        compressedHooks = 0.00001;
        cellRadius = repulsionRadius; /// in micrometers
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

    void calcStokesDrag(){ /// returns the drag force of X component
        xStokesDrag = 6 * 3.14159 * cellRadius/1000000 * fluidViscosity * xvelocity;
        yStokesDrag = 6 * 3.14159 * cellRadius/1000000 * fluidViscosity * yvelocity;
        /// velocity and cellRadius in micrometres
    }

    void calculateVelocity(){
        xvelocity += timestep * (xSpringForce - xStokesDrag)/(cellMass/1000000000);
        yvelocity += timestep * (ySpringForce - yStokesDrag)/(cellMass/1000000000);
    }

    /// make a step in the given direction
    void step(){
#if ( 1 )
        /// change the velocity depending on the acceleration
        x += timestep * xvelocity;
        y += timestep * yvelocity;
#else
        //use 2 Gaussian random numbers in polar coordinates:
        float angle = 2 * PI * prand();
        float norm = alpha * sqrt( -log( prand() ));
        x += norm * cos(angle);
        y += norm * sin(angle);
#endif
        bounce();
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
