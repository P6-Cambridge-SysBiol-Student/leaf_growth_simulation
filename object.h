/*
Working simulation of points orbiting around the centre
 */

///Object class
class Object
{
public:  // these are attributes i think can be called outside of the script
    /// member variables:
    int color;
    float x, y; /// displacement
    float xvelocity, yvelocity; /// velocity
    float hooks;
    
    /// initialize in a random position
    void reset()
    {
        color = 1;
        x = xBound * srand();
        xvelocity = 0.2 * srand();
        y = yBound * srand();
        yvelocity = 0.2 * srand();
        hooks = 0.005;
    }
    
    /// call initialize
    Object()
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
        //use 2 uniform random numbers
        x += xvelocity;
        y += yvelocity;
#else
        //use 2 Gaussian random numbers in polar coordinates:
        float angle = 2 * PI * prand();
        float norm = alpha * sqrt( -log( prand() ));
        x += norm * cos(angle);
        y += norm * sin(angle);
#endif
        bounce();
    }

    void accelerate()
    {
        xvelocity = xvelocity -x*hooks;
        yvelocity = yvelocity -y*hooks;
    }
    
    /// partial display: this needs to be called between glBegin() and glEnd()
    void display()
    {
        //we use transparency to visualize overlapping particules
        if ( color == 1 )
            glColor4f(1.0, 1.0, 0.0, 0.5);
        else
            glColor4f(0.3, 0.3, 0.3, 0.5);
        glVertex2f(x, y);
    }
};
