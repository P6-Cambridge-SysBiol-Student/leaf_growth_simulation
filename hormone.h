//
// Created by finley on 09/02/23.
//

#ifndef FRAP_HORMONE_H
#define FRAP_HORMONE_H
#include <float.h>

#endif //FRAP_HORMONE_H

void startHormoneBD(double inputStartTime){
    static bool flag = false;
    if ((currentTime > inputStartTime) and (flag == false)){
        flag = true;
        /// find the point closest to the hormone Origin
        int closest_point_index = -1;
        double squareMinDist = 1000*1000*xBound;
        for (int i = 0; i < nbo; i++) {
            double squareDisFromOrigin = (pointsArray[i].disVec - horm2Source1).magnitude_squared();
            if (squareDisFromOrigin < squareMinDist) {
                squareMinDist = squareDisFromOrigin;
                closest_point_index = i;
            }
        }
        /// set this point as the hormone producer
        pointsArray[closest_point_index].isHormone1Producer = true;
        printf("Closest point is point %d\n", closest_point_index);
    }
    else{
    }
}

void calcHormBirthDeath(){
    for (int i = 0; i < nbo; i++){
        Point& cell = pointsArray[i]; /// alias for pointsArray[i]
        /// calculate amount of hormone made by producers
        if ((cell.isHormone1Producer == true)){
            cell.produceHormone1BD(hormone1ProdRate);
            cell.degradeHormone1BD(hormone1DegRate);
        }
        else{
            cell.degradeHormone1BD(hormone1DegRate);
        }
    }
}

void hormReactDiffuse(double inputStartTime) {
    static bool flag = false;
    if ((currentTime > inputStartTime) and (flag == false)) {
        flag = true;
        /// find the point closest to the hormone Origin
        int closest_point_source1_index = -1;
        int closest_point_source2_index = -1;
        double squareMinDist1 = 1000 * 1000 * xBound;
        for (int i = 0; i < nbo; i++) {
            double squareDisFromOrigin = (pointsArray[i].disVec - horm2Source1).magnitude_squared();
            if (squareDisFromOrigin < squareMinDist1) {
                squareMinDist1 = squareDisFromOrigin;
                closest_point_source1_index = i;
            }
        }
        double squareMinDist2 = 1000 * 1000 * xBound;
        for (int k = 0; k < nbo; k++) {
            double squareDisFromOrigin = (pointsArray[k].disVec - horm2Source1).magnitude_squared();
            if (squareDisFromOrigin < squareMinDist1) {
                squareMinDist2 = squareDisFromOrigin;
                closest_point_source2_index = k;
            }
            /// set this point as the hormone producer
            pointsArray[closest_point_source1_index].isHormone2Producer = true;
            pointsArray[closest_point_source2_index].isHormone2Producer = true;
        }
    }
    for (int i = 0; i < nbo; i++) {
        Point &cell = pointsArray[i]; /// alias for pointsArray[i]
        /// in reaction diffusion all cells produce horm1
        if (cell.isHormone2Producer == true) {
            cell.produceHormone1ReactD( RDfeedRate);
            cell.productHormone2ReactD( 1*RDfeedRate);
            cell.react1With2( reactRate1to2);
            cell.degradeHormone2ReactD( RDkillRate,  RDfeedRate);
            //printf("Point %d is a horm2 producer\n", i);
        } else {
            cell.produceHormone1ReactD( RDfeedRate);
            cell.react1With2( reactRate1to2);
            cell.degradeHormone2ReactD( RDkillRate,  RDfeedRate);
        }
    }
}

void v1DiffuseHorm(int** neighbourhoods) {

    for (int i = 0; i < nbo; i++) { ///for each primary point in pointsArray (iterates through each point using i)
        Point &centre = pointsArray[i]; /// alias for pointsArray[i]
        for (int l = 0; l < NAW; l++) {
            Point &neighbour = pointsArray[neighbourhoods[i][l]];
            if (neighbourhoods[i][l] != -1) {
                /// using squared magnitudes here is computationally faster
                if ((neighbour.disVec - centre.disVec).magnitude_squared() <
                    (0.2 * centre.cellRadius * 0.2 * centre.cellRadius)) {
                } /// stops diffusion if points overlap
                else {
                    /// find the magnitude of distance between the neighbouring point and the central point
                    double magnitudeOfDistance = (centre.disVec - neighbour.disVec).magnitude(); // m

                    /// find difference in hormone amount between cells
                    double hormone1ConcnDiff = centre.myTotalHormone1 - neighbour.myTotalHormone1;  //n / m
                    double hormone2ConcnDiff = centre.myTotalHormone2 - neighbour.myTotalHormone2;

                    double hormone1ConcnGrad = hormone1ConcnDiff / (magnitudeOfDistance * magnitudeOfDistance); //n / m^2
                    double hormone2ConcnGrad = hormone2ConcnDiff / (magnitudeOfDistance * magnitudeOfDistance);
                    /// diffuse the hormone from the centre to neighbour
                    neighbour.myDeltaHormone1 += timestep*(hormone1DiffCoeff * hormone1ConcnGrad * centre.cellRadius); //  n = t * (m^2/t * n/m * m)
                    centre.myDeltaHormone1 -= timestep*(hormone1DiffCoeff * hormone1ConcnGrad * centre.cellRadius);

                    neighbour.myDeltaHormone2 += timestep*(hormone2DiffCoeff * hormone2ConcnGrad * centre.cellRadius); //  n = t * (m^2/t * n/m * m)
                    centre.myDeltaHormone2 -= timestep*(hormone2DiffCoeff * hormone2ConcnGrad * centre.cellRadius);
                }
            }
        }
    }
    double sumHorm1 = 0;
    double sumHorm2 = 0;

    for (int j = 0; j < nbo; j++) {
        Point &cell = pointsArray[j];

        sumHorm1 += cell.myTotalHormone1;
        sumHorm2 += cell.myTotalHormone2;
    }
#if DEBUG
printf("The sum of hormone1 is %f\nThe sum of hormone 2 is %f \n", sumHorm1, sumHorm2); /// test conservation of hormone
#endif
}

double sumHormone2(){
    double sumHorm1 = 0;
    double sumHorm2 = 0;

    for (int j = 0; j < nbo; j++) {
        Point &cell = pointsArray[j];

        sumHorm1 += cell.myTotalHormone1;
        sumHorm2 += cell.myTotalHormone2;
    }
    return sumHorm2;
}

double findMaxHormone1(){
    double maxHormone = 0;
    for (int i = 0; i<nbo; i++){
        if (pointsArray[i].myTotalHormone1 > maxHormone){
            maxHormone = pointsArray[i].myTotalHormone1;
        }
    }
    return maxHormone;
}

double findMinHormone1(){
    double minHormone = DBL_MAX;
    for (int i = 0; i<nbo; i++){
        if (pointsArray[i].myTotalHormone1 < minHormone){
            minHormone = pointsArray[i].myTotalHormone1;
        }
    }
    return minHormone;
}

double findMaxHormone2(){
    double maxHormone = 0;
    for (int i = 0; i<nbo; i++){
        if (pointsArray[i].myTotalHormone2 > maxHormone){
            maxHormone = pointsArray[i].myTotalHormone2;
        }
    }
    return maxHormone;
}

double findMinHormone2(){
    double minHormone = DBL_MAX;
    for (int i = 0; i<nbo; i++){
        if (pointsArray[i].myTotalHormone2 < minHormone){
            minHormone = pointsArray[i].myTotalHormone2;
        }
    }
    return minHormone;
}

void globalUpdateHormone(){
    for (int i = 0; i<nbo; i++){
        pointsArray[i].updateTotalHormone();
    }
}

void hormoneExpandEffect(){
    for (int i = 0; i < nbo; i++){
        Point& centre = pointsArray[i];
        centre.cellRadius = centre.cellRadiusBase + (horm1Efficacy * centre.myTotalHormone1 * SCALING_FACTOR);
    }
}

