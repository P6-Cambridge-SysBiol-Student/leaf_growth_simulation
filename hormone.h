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
    double DiffCoeff, StartTime, initConcn, hormoneExpansionEffect, efficacy;

    /// constructor
    Hormone(double inputDiffCoeff, double inputStartTime,
            double inputInitConcn, double inputExpansion, double inputEfficacy) :

            DiffCoeff(inputDiffCoeff),
            StartTime(inputStartTime),
            initConcn(inputInitConcn),
            hormoneExpansionEffect(inputExpansion),
            efficacy(inputEfficacy) {}
};

