//
// Created by finley on 09/02/23.
//

#ifndef FRAP_HORMONE_H
#define FRAP_HORMONE_H

#endif //FRAP_HORMONE_H

/// includes parameters relevant to diffusion and behaviour of hormones
class Hormone{
public:
    /// member variables
    double myDiffCoeff, myRateOfProd, myRateOfDeg, myExpandEffect, myEfficacy, myMitosisRateEffect, myMitosisOrientEffect;

    /// default constructor
    Hormone() :
            myDiffCoeff(0.0),
            myRateOfProd(0.0),
            myRateOfDeg(0.0),
            myExpandEffect(0.0),
            myMitosisRateEffect(0.0),
            myMitosisOrientEffect(0.0),
            myEfficacy(0.0)
    {}

    /// constructor with arguments
    Hormone(double inputDiffCoeff, double inputProdRate, double inputExpansion, double inputEfficacy,
            double inputDegRate, double inputMitosisRateEffect, double inputMitosisOrientEffect) :
            myDiffCoeff(inputDiffCoeff),
            myRateOfProd(inputProdRate),
            myRateOfDeg(inputDegRate),
            myExpandEffect(inputExpansion),
            myMitosisRateEffect(inputMitosisRateEffect),
            myMitosisOrientEffect(inputMitosisOrientEffect),
            myEfficacy(inputEfficacy)
    {}
};
