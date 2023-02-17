/*
Class which represents cells which repel and attract each other on a 2D plane, based on a delaunay triangulation.
Points contain hooks constants which are used to simulate spring-like physics
This file contains information about Object which contains
 */
#include <math.h>
#include "sigmoid.h"
/// for the compiler this doesn't slow down the programme

///Points class, contains info about xy displacement, velocity, and acceleration towards centre
class Point
{
public:  /// these are attributes that can be called outside of the script
    /// member variables:
    vector2D disVec = vector2D(double (0.1*xBound*mySrand()), double (0.1*yBound*mySrand())); /// sets x and y values randomly
    vector2D velVec = vector2D(0.0001, 0.0001); /// initial velocities set to very small, prevents bugs
    vector2D springVec = vector2D(0, 0);  /// would be set (0, 0) by default but just in case
    vector2D mitosisOrient = vector2D(1, 1);

    double extendedHooks, compressedHooks, innerMultiplier, innerCompressedHooks;  /// hooks constant for attracting points back to the centre
    double cellRadiusBase, cellRadius;
    double cellMass;
    double probOfDividing;
    int color;

    /// members related to hormone function
    bool isHormoneProducer = false;
    double myTotalHormone = 0;
    double myHormConc = 0; // TODO add function to find concentration, check units in diff equation
    double myDeltaHormone = 0; /// keeps track of amount of hormone gained/lost
    double myDiffCoeff = 0;
    double myRateOfProd = 0;
    double myRateOfDeg = 0.1;
    double myExpandEffect = 0;
    double myHormoneSensitivity = 0;

    /// members related to cell division
    bool wasMotherCell = true;
    bool newDaughterCell = false;

    /// initialize each point in a random position with random x and y velocities
    void reset()
    {
        extendedHooks = 10;
        compressedHooks = 30;
        innerMultiplier = 2;
        innerCompressedHooks = innerMultiplier * compressedHooks;
        cellRadiusBase = 4 * SCALING_FACTOR; /// in micrometers
        cellRadius = cellRadiusBase;
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
        disVec += (timestep/(mobilityCoefficient * cellRadius/SCALING_FACTOR)) * springVec;
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

    void produceHormone(double inputProdRate){
        myRateOfProd = inputProdRate;
        myTotalHormone += myRateOfProd;
    }

    void degradeHormone(double inputDegRate){
        myTotalHormone -= myTotalHormone*inputDegRate;
    }


    void sigmoidDisplayHormone() {
        double sigmoidHormConc = sigmoid((12*myTotalHormone)-5); /// this shifts curve so that color scales from 0 to 1
        glColor4f((sigmoidHormConc), 0, (1 - sigmoidHormConc), 1);
        glVertex2f(disVec.xx, disVec.yy);
    }

    void linearDisplayHormone() {
        double linearHormConc = myTotalHormone / 1.5;
        glColor4f((linearHormConc), (0.5 - 0.5*linearHormConc), (1 - linearHormConc), 1);
        glVertex2f(disVec.xx, disVec.yy);
    }

    double divisionProb(double maxProbOfDiv, int numCurrentCells, int finalTotCells) { /// each cell has a p(mitosis) varied by number of existing points, cell size etc.
        double divisionProb = -(maxProbOfDiv / finalTotCells) * numCurrentCells + maxProbOfDiv; /// just a linear equation
        return divisionProb;
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
