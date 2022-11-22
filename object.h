/*
 Basic diffusion and FRAP simulation
 Francois Nedelec, Copyright Cambridge University 2021
 */

///Object class
class Object
{
public:  // these are attributes i think can be called outside of the script
    /// member variables:
    int color;
    float x, y; ///< position
    
    /// initialize in a random position
    void reset()
    {
        color = 1;
        x = xBound * srand();
        y = yBound * srand();
    }
    
    /// call initialize
    Object()
    {
        reset();
    }
    
    /// true if particle is within distance R of (xc, yc)
    bool within(float R, float xc=0, float yc=0)
    {
        return ( (x-xc)*(x-xc) + (y-yc)*(y-yc) < R*R );
    }
    
    /// particles bounce off walls
    void bounce()
    {
        if ( x >  xBound )  x =  2*xBound - x;
        if ( x < -xBound )  x = -2*xBound - x;
        if ( y >  yBound )  y =  2*yBound - y;
        if ( y < -yBound )  y = -2*yBound - y;
    }
    
    /// make a step in a random direction
    void step()
    {
#if ( 1 )
        //use 2 uniform random numbers
        x += alpha * srand();
        y += alpha * srand();
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
