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
    double DiffCoeff, rateOfProd, rateOfDeg, expandEffect, efficacy, currentConcn;

    /// constructor
    Hormone(double inputDiffCoeff, double inputProdRate, double inputExpansion, double inputEfficacy, double inputDegRate) :
            DiffCoeff(inputDiffCoeff),
            rateOfProd(inputProdRate),
            rateOfDeg(inputDegRate),
            expandEffect(inputExpansion),
            efficacy(inputEfficacy)
            {}
};

