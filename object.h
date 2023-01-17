/*
Working simulation of points orbiting around the centre of the window as a test for later physics
This file contains information about Object which contains
 */

///Points class, contains info about xy displacement, velocity, and acceleration towards centre
class Point
{
public:  /// these are attributes that can be called outside of the script
    /// member variables:
    int color;
    float x, y; /// displacement
    float xvelocity, yvelocity; /// velocity
    float xacceleration, yacceleration; /// acceleration
    float extendedHooks, compressedHooks;  /// hooks constant for attracting points back to the centre
    float cellRadius;
    
    /// initialize each point in a random position with random x and y velocities
    /// currently these are set to start points randomly at the centre bottom to mimic plant leaves
    void reset()
    {
        x = 0.1*(xBound * srand());
        // xvelocity = xBound/1500 * srand();
        xvelocity = 0;
        xacceleration = 0;

        y = 0.1*(yBound * srand());
        // yvelocity = xBound/1500 * srand();
        xvelocity = 0;
        yacceleration = 0;

        extendedHooks = 0.005;
        compressedHooks = 1;

        cellRadius = xBound / 30;
        color = 1;
    }
    
    /// call initialize
    Point()
    {
        reset();
    }

    /// particles bounce off walls
    void bounce()
    {
        if ( x >  xBound )  x =  2*xBound - x, xvelocity = -xvelocity;
        if ( x < -xBound )  x = -2*xBound - x, xvelocity = -xvelocity;
        if ( y >  yBound )  y =  2*yBound - y, yvelocity = -yvelocity;
        if ( y < -yBound )  y = -2*yBound - y, yvelocity = -yvelocity;
    }
    
    /// make a step in the given direction
    void step()
    {
#if ( 1 )
        /// change the velocity depending on the acceleration
        /*xvelocity += xacceleration;
        yvelocity += yacceleration;*/ /// dont think i need this due to high viscotiy
        /// change displacement by a unit of velocity
        x += xvelocity*timestep;
        y += yvelocity*timestep;
#else
        //use 2 Gaussian random numbers in polar coordinates:
        float angle = 2 * PI * prand();
        float norm = alpha * sqrt( -log( prand() ));
        x += norm * cos(angle);
        y += norm * sin(angle);
#endif
        bounce();
    }

    double stokesDrag(float inputVelocity) /// returns the drag force
    {
        double dragForce = 6 * 3.14159 * cellRadius * fluidViscosity * inputVelocity;
        return dragForce;
    }

    void accelerateToCentre() /// change the velocity of each point proportional to the distance from centre and hooks constant
    {
        xvelocity = xvelocity -x*extendedHooks;
        yvelocity = yvelocity -y*extendedHooks;
    }
    
    /// partial display: this needs to be called between glBegin() and glEnd()
    void displayYellow()
    {
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
            glColor4f(1.0, 1.0, 1.0, 1);
        else
            glColor4f(0.3, 0.3, 0.3, 1);
        glVertex2f(x, y);
    }
};
