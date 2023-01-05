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
    float hooks;  /// hooks constant for attracting points back to the centre
    
    /// initialize each point in a random position with random x and y velocities
    /// currently these are set to start points randomly at the centre bottom to mimic plant leaves
    void reset()
    {
        x = 0.25*(xBound * srand());
        xvelocity = xBound/1500 * srand();
        xacceleration = 0;

        y = -abs(yBound * srand());
        yvelocity = xBound/1500 * srand();
        yacceleration = 0;

        hooks = 0.05;
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

    void accelerateToCentre() /// change the velocity of each point proportional to the distance from centre and hooks constant
    {
        xvelocity = xvelocity -x*hooks;
        yvelocity = yvelocity -y*hooks;
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
            glColor4f(1.0, 1.0, 1.0, 0.5);
        else
            glColor4f(0.3, 0.3, 0.3, 0.5);
        glVertex2f(x, y);
    }
};
