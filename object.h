/*
Class which represents cells which repel and attract each other on a 2D plane, based on a delaunay triangulation.
Points contain hooks constants which are used to simulate spring-like physics
This file contains information about Object which contains
 */
#include <math.h>
#include "sigmoid.h"
/// for the compiler this doesn't slow down the programme

class Point
{
public:  /// these are attributes that can be called outside of the script
    /// member variables:
    vector2D disVec = vector2D(double (0.7*xBound*mySrand()), double (0.7*yBound*mySrand())); /// sets x and y values randomly
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
    double myTotalHormone1 = 0;
    double myDeltaHormone1 = 0; /// keeps track of amount of hormone gained/lost
    double myRateOfProd1 = 0;
    double myExpandEffect = 0;
    double myHormoneSensitivity = 0;

    double myTotalHormone2 = 0;
    double myDeltaHormone2 = 0;
    double myRateOfProd2 = 0;

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
        cellRadiusBase = 3 * SCALING_FACTOR; /// in micrometers
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

    void displayGreen(){
        if (color == 1)
            glColor4f(0, 1, 0, 1);
        else
            glColor4f(0, 1, 0, 1);
        glVertex2f(disVec.xx, disVec.yy);
    }

    void produceHormone1(double inputProdRate){
        myRateOfProd1 = inputProdRate;
        myDeltaHormone1 = timestep*(myRateOfProd1);
        printf("Locally myDeltaHormone1 is %f\n", myDeltaHormone1);
    }

    void produceHormone2(double inputProdRate){
        myRateOfProd2 = inputProdRate;
        myDeltaHormone2 = timestep*(myRateOfProd2);
    }

    void degradeHormone1(double inputDegRate){
        myDeltaHormone1 = -timestep*(myTotalHormone1*inputDegRate);
    }

    void degradeHormone2(double inputDegRate){
        myDeltaHormone2 = -timestep*(myDeltaHormone2*inputDegRate);
    }

    void react1With2(double input1and2ReactRate){
        double deltaHormoneReaction = timestep*(input1and2ReactRate*myTotalHormone1*myTotalHormone2*myTotalHormone2);
        myDeltaHormone1 -= deltaHormoneReaction;
        myDeltaHormone2 += deltaHormoneReaction;
    }

    void updateTotalHormone(){
        myTotalHormone1 += myDeltaHormone1;
        myTotalHormone2 += myDeltaHormone2;
    }

    void sigmoidDisplayHormone() {
        double sigmoidHormConc = sigmoid((12*myTotalHormone1)-5); /// this shifts curve so that color scales from 0 to 1
        glColor4f((sigmoidHormConc), 0, (1 - sigmoidHormConc), 1);
        glVertex2f(disVec.xx, disVec.yy);
    }

    void linearDisplayHormone() {
        double linearHormConc = myTotalHormone1 / 1.5;
        glColor4f((linearHormConc), (0.5 - 0.5*linearHormConc), (1 - linearHormConc), 1);
        glVertex2f(disVec.xx, disVec.yy);
    }

    double divisionProb(double maxProbOfDiv, int numCurrentCells, int finalTotCells) { /// each cell has a p(mitosis) varied by number of existing points, cell size etc.
        double divisionProb = -(maxProbOfDiv / finalTotCells) * numCurrentCells + maxProbOfDiv + (myTotalHormone1 * horm1Efficacy); /// just a linear equation
        return divisionProb;
    }
};
