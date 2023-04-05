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
    vector2D disVec = vector2D(double (0.01*xBound*mySrand()), double (0.01*yBound*mySrand())); /// sets x and y values randomly
    vector2D velVec = vector2D(0.0001, 0.0001); /// initial velocities set to very small, prevents bugs
    vector2D springVec = vector2D(0, 0);  /// would be set (0, 0) by default but just in case
    vector2D mitosisOrient = vector2D(1, 1);

    double extendedHooks, compressedHooks, innerMultiplier, innerCompressedHooks;  /// hooks constant for attracting points back to the centre
    double cellRadiusBase, cellRadius;
    double cellMass;
    double probOfDividing;
    int color;

    /// members related to hormone function
    bool isHormone1Producer = false;
    double myTotalHormone1 = 0;
    double myDeltaHormone1 = 0; /// keeps track of amount of hormone gained/lost
    double myRateOfProd1 = 0;
    double myRateOfDeg1 = 0;
    double myExpandEffect = 0;
    double myHormoneSensitivity = 0;

    bool isHormone2Producer = false;
    double myTotalHormone2 = 0;
    double myDeltaHormone2 = 0;
    double myRateOfProd2 = 0;

    /// members related to cell division

    /// initialize each point in a random position with random x and y velocities
    void reset()
    {
        extendedHooks = 10;
        compressedHooks = 30;
        innerMultiplier = 2;
        innerCompressedHooks = innerMultiplier * compressedHooks;
        cellRadiusBase = 1 * SCALING_FACTOR; /// in micrometers
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
/// BD here represents Birth-death process, need new functions for reaction-diffusion
    void produceHormone1BD(double inputProdRate){
        myRateOfProd1 = inputProdRate;
        myDeltaHormone1 += myRateOfProd1;
    }

    void degradeHormone1BD(double inputDegRate){
        myRateOfDeg1 = inputDegRate;
        myDeltaHormone1 += myRateOfDeg1*myTotalHormone1;
    }

    void produceHormone1ReactD(double inputFeedRate){
        myDeltaHormone1 += inputFeedRate*(1-myTotalHormone1);
    }

    void produceHormone2ForInit(double inputFeedRate){
        myDeltaHormone2 += inputFeedRate;
    }

    void productHormone2ReactD(double inputFeedRate){
        myDeltaHormone2 += inputFeedRate;
    }

    void degradeHormone2ReactD(double inputKillRate, double inputFeedRate) {
        myDeltaHormone2 += -(inputFeedRate + inputKillRate) * myTotalHormone2;
        /// feedrate added to killrate so killrate is never < feedrate
    }

    void react1With2(double input1and2ReactRate){
        double deltaHormoneReaction = (input1and2ReactRate*myTotalHormone1*myTotalHormone2*myTotalHormone2);
        myDeltaHormone1 -= deltaHormoneReaction;
        myDeltaHormone2 += deltaHormoneReaction;
    }

    void updateTotalHormone(){
        myTotalHormone1 += timestep*(myDeltaHormone1);
        myTotalHormone2 += timestep*(myDeltaHormone2);
        myDeltaHormone1 = 0;
        myDeltaHormone2 = 0;
        if (myTotalHormone1 < 0){
            myTotalHormone1 = 0;
        }
        if (myTotalHormone2 < 0){
            myTotalHormone2 = 0;
        }
    }

    void sigmoidDisplayHormone() {
        double sigmoidHormConc = sigmoid((12*myTotalHormone1)-5); /// this shifts curve so that color scales from 0 to 1
        glColor4f((sigmoidHormConc), 0, (1 - sigmoidHormConc), 1);
        glVertex2f(disVec.xx, disVec.yy);
    }

    void linearDisplayHormone2(double inputMaxHormone) {
        double linearHormConc = myTotalHormone2 / inputMaxHormone;
        glColor4f((linearHormConc), (0.5 - 0.5*linearHormConc), (1 - linearHormConc), 1);
        glVertex2f(disVec.xx, disVec.yy);
    }

    void linearDisplayHormone1(double inputMaxHormone) {
        double linearHormConc = myTotalHormone1 / 0.5*inputMaxHormone;
        glColor4f((linearHormConc), (0.5 - 0.5*linearHormConc), (1 - linearHormConc), 1);
        glVertex2f(disVec.xx, disVec.yy);
    }


    double divisionProb(double maxProbOfDiv, int numCurrentCells, int finalTotCells) { /// each cell has a p(mitosis) varied by number of existing points, cell size etc.
        double divisionProb = -(maxProbOfDiv / finalTotCells) * numCurrentCells + maxProbOfDiv + (myTotalHormone1 * horm1Efficacy); /// just a linear equation
        return divisionProb;
    }
};
