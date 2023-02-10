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
    vector2D disVec = vector2D(double (0.8*xBound*mySrand()), double (0.8*yBound*mySrand())); /// sets x and y values randomly
    vector2D velVec = vector2D(0.0001, 0.0001); /// initial velocities set to very small, prevents bugs
    vector2D springVec = vector2D(0, 0);  /// would be set (0, 0) by default but just in case
    vector2D mitosisOrient = vector2D(1, 1);

    double extendedHooks, compressedHooks;  /// hooks constant for attracting points back to the centre
    double cellRadius;
    double cellMass;
    int color;

    /// members related to hormone function

    bool isHormoneProducer = false;
    double myTotalHormone = 0;
    double myMaxTotalHorm  = 2;  /// TODO i need a max value to scale
    double myDiffCoeff = 0;
    double myRateOfProd = 0;
    double myRateOfDeg = 0;
    double myExpandEffect = 0;
    double myHormoneSensitivity = 0;

    /// initialize each point in a random position with random x and y velocities
    void reset()
    {
      extendedHooks   = 0.3;
      compressedHooks = 3;
      cellRadius = 10 * 100000; /// in micrometers
      cellMass = 1; /// in nanograms
      color = 1;
    }
    
    /// call initialize
    Point(){
        reset();
    }

    /// particles bounce off walls
    void bounce(){
        if ( disVec.xx >  xBound )  disVec.xx =  2*xBound - disVec.xx, velVec.xx = -velVec.xx;
        if ( disVec.xx >  xBound )  disVec.xx =  2*xBound - disVec.xx, velVec.xx = -velVec.xx;
        if ( disVec.yy >  xBound )  disVec.yy =  2*xBound - disVec.yy, velVec.yy = -velVec.yy;
        if ( disVec.yy >  xBound )  disVec.yy =  2*xBound - disVec.yy, velVec.yy = -velVec.yy;
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
        glVertex2f(disVec.xx, disVec.yy);
    }

    void displayWhite()
    {
        /// transparency used to visualize overlapping particles
        if ( color == 1 )
            glColor4f(1.0, 1.0, 1.0, 0.5);
        else
            glColor4f(0.3, 0.3, 0.3, 0.5);
        glVertex2f(disVec.xx, disVec.yy);
    }

    void produceHormone(int inputProdRate){
        myRateOfProd = inputProdRate;
        myTotalHormone += myRateOfProd;

    }

    void degradeHormone(){
        myTotalHormone += - myTotalHormone*myRateOfDeg;
    }

    void displayHormone() {
        if (isHormoneProducer == true) {
            glColor4f(0, 1, 0, 1);
        } // TODO fix the lack of graphics for showing hormones
        else if(isHormoneProducer == false){
            double normHormConc = (1 + myTotalHormone)/myMaxTotalHorm;
            double sigmoidDenominator = pow(exp, -(myTotalHormone-5));
            double squishedHormoneValue = 1 / (1 + pow(exp, -(myTotalHormone-5)));
            glColor4f((normHormConc), 0, (1 - normHormConc), 1);
        }
        glVertex2f(disVec.xx, disVec.yy);
    }


    /* ignoring mitosis for now, focusing on morphogens
    void wasMotherMitosisDisplace(){
        mitosisOrient.normalise();
        disVec += 0.25*cellRadius*mitosisOrient; /// the cell that inherits mothers indexing
    }

    void isNewCellMitosisDisplace(){ /// is called when definiting position of newly created cell
        mitosisOrient.normalise();
    } */
};
