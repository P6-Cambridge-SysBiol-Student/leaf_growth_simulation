/*
Class which represents cells which repel and attract each other on a 2D plane, based on a delaunay triangulation.
Points contain hooks constants which are used to simulate spring-like physics
This file contains information about Object which contains
 */
#include <math.h>
/// for the compiler this doesn't slow down the programme

///Points class, contains info about xy displacement, velocity, and acceleration towards centre
class Point
{
public:  /// these are attributes that can be called outside of the script
    /// member variables:
    int color;
    vector2D disVec = vector2D(double (xBound*srand()), double (yBound*srand())); /// sets x and y values randomly
    vector2D velVec = vector2D(0.0001, 0.0001); /// initial velocities set to very small, prevents bugs
    vector2D springVec = vector2D(0, 0);  /// would be set (0, 0) by default but just in case
    double extendedHooks, compressedHooks;  /// hooks constant for attracting points back to the centre
    double cellRadius;
    double cellMass;
    
    /// initialize each point in a random position with random x and y velocities
    /// currently these are set to start points randomly at the centre bottom to mimic plant leaves
    void reset()
    {
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
        if ( disVec.x >  xBound )  disVec.x =  2*xBound - disVec.x, velVec.x = -velVec.x;
        if ( disVec.x >  xBound )  disVec.x =  2*xBound - disVec.x, velVec.x = -velVec.x;
        if ( disVec.y >  xBound )  disVec.y =  2*xBound - disVec.y, velVec.y = -velVec.x;
        if ( disVec.y >  xBound )  disVec.y =  2*xBound - disVec.y, velVec.y = -velVec.x;
    }


    //TODO find youngs modulus of the leaf (stiffness per unit area of the cross-section)

    /// make a step in the given direction
    void step(){
        /// change the velocity depending on the acceleration
        /// TODO the value given to
        disVec += (timestep/(mobilityCoefficient * cellRadius/100000)) * springVec;
    }

    /// partial display: this needs to be called between glBegin() and glEnd()
    void displayYellow(){
        /// transparency used to visualize overlapping particles
        if ( color == 1 )
            glColor4f(1.0, 1.0, 0.0, 0.5);
        else
            glColor4f(0.3, 0.3, 0.3, 0.5);
        glVertex2f(disVec.x, disVec.y);
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
