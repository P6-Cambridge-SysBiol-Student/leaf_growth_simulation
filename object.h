/*
 Basic diffusion and FRAP simulation
 Francois Nedelec, Copyright Cambridge University 2021
 */
#include <math.h>

float timegap = 0.16;

///Object class
class pointClass
{
public:  // these are attributes i think can be called outside of the script
    /// member variables:
    int color;
    float x, y; ///< position
    float xvelocity, yvelocity;
    
    /// initialize in a random position
    void reset()
    {
        color = 1;
        x = xBound * srand();
        y = yBound * srand();
        xvelocity = 0.2 * srand();
        yvelocity = 0.2 * srand();

    }
    
    /// call initialize
    pointClass()
    {
        reset();
    }

/*    /// particles bounce off walls
    void bounce()
    {
        if ( x >  xBound )  x =  2*xBound - x, xvelocity = -xvelocity;
        if ( x < -xBound )  x = -2*xBound - x, xvelocity = -xvelocity;
        if ( y >  yBound )  y =  2*yBound - y, yvelocity = -yvelocity;
        if ( y < -yBound )  y = -2*yBound - y, yvelocity = -yvelocity;
    }
*/
    /// make a step in a specified
    void step()
    {
        int num = 2;
        float hooks = 2;
        float xAcceleration = 0.5 * hooks * x;  /// modelling attraction to centre as a spring
        float yAcceleration = 0.5 * hooks * y;

#if ( 1 )
        //use 2 uniform random numbers
        x += xvelocity*timegap - xAcceleration*pow(timegap, 2);
        y += yvelocity*timegap - yAcceleration*pow(timegap, 2);
#endif
        num = 3;
    }


    /// particles are constantly attracted to the centre
    void attract()
    {
        float XdispFromCentre = x - 0;  /// these are redundant but I'm keeping for clarity / for working with pow()
        float YdispFromCentre = y - 0;
        float hooksConstant = 0.01;

        xvelocity = xvelocity - 1/hooksConstant*XdispFromCentre;  /// accelerate the velocities towards the centre
        yvelocity = yvelocity - 1/YdispFromCentre*YdispFromCentre;
    }


    /// partial display: this needs to be called between glBegin() and glEnd()
    void display()
    {
        //we use transparency to visualize overlapping particles
        if ( color == 1 )
            glColor4f(1.0, 1.0, 0.0, 0.5);
        else
            glColor4f(0.3, 0.3, 0.3, 0.5);
        glVertex2f(x, y);
    }
};
