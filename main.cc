/*
 * Elementary simulation using GLFW + OpenGL for display modified for use in PlantSim by Finley Webb
 * Original code by Francois J Nedelec, Cambridge University, 13 Nov 2021, 11 Oct 2022
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define DEBUG false
#define DISPLAY true /// set to true to display
#define BENCHMARK false /// set to true to benchmark (not bottlenecked by printing or displaying)
#define GLAD_GL_IMPLEMENTATION
#include <glad/gl.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <vector>

#include "random.h"
#include "arrays.h"
#include "vector.h"
#include "param.h"
#include "hormone.h"
#include "object.h"
#include "polish.h"
#include "Clarkson-Delaunay.cpp"  /// this is slightly odd, would be better to compile them seperately and link together
#include "graphics.h"

/// TODO Hormones simulation (simple & reaction-diffusion)


///-----------------------------------------------------------------------------


/// evolves system, stepping points forward and accelerating velocity
static void animate(){
    realTime += delta;
    for ( int i = 0; i < nbo; ++i ) {
        pointsArray[i].step();
    }
}

/// creates an array of xy co-ords for the delaunay triangulation function, then execute it
void create_triangles_list(){
    float xyValuesArray[nbo][2];
    for ( int i = 0; i < nbo; i++){
        xyValuesArray[i][0] = pointsArray[i].disVec.xx;
        xyValuesArray[i][1] = pointsArray[i].disVec.yy;
    }

    numTriangleVertices = 0;
    triangleIndexList = BuildTriangleIndexList((void*)xyValuesArray, (float)1.0, nbo, (int)2, (int)1, &numTriangleVertices);

#if DEBUG
    printf("\nThere are %d points moving around \n", nbo);
    printf("\nThe number of vertices defined by numTriangleVertices is %d\n", numTriangleVertices);
    /*printf("int has size %ld \n", sizeof(int ));
    printf("\ntriangleIndexList contains the values: ");
    for (int i = 0; i < numTriangleVertices; i++)
        printf("%u, ", triangleIndexList[i]);
    printf("\n");*/
#endif ///DEBUG


}

/// checks through the pointsConnected array to see if secondary point is present
bool noDuplicateCheck(int indexValueToCheck, int arrayToCheck[], int max){
    static bool unique = true;

    for (int i = 0; i <= max; i++) {  /// going up to total instead of over all the array is faster
        if (arrayToCheck[i] == indexValueToCheck) {
            unique = false;
        }
    }

    if ((unique = true)){
        return true;
    }
    else{
        return false;
    }
}

/// create a neighbourhood array (seperated out from CalculateSpringForces function)
void fill2DArrayNeighbourhoods(int** neighbourhoods, int* total, int rows){
    for (int v = 0; v < numTriangleVertices; v+=3) {
        /// add neighbours of the first value of the triangle to its row in the neighbourhood array
        total[triangleIndexList[v]]++;
        neighbourhoods[triangleIndexList[v]][total[triangleIndexList[v]]] = triangleIndexList[v + 1];
        total[triangleIndexList[v]]++;
        neighbourhoods[triangleIndexList[v]][total[triangleIndexList[v]]] = triangleIndexList[v + 2];
        /// add neighbours of the second value in triangleIndex List to its row in the neighbour array
        total[triangleIndexList[v + 1]]++;
        neighbourhoods[triangleIndexList[v + 1]][total[triangleIndexList[v + 1]]] = triangleIndexList[v];
        total[triangleIndexList[v + 1]]++;
        neighbourhoods[triangleIndexList[v + 1]][total[triangleIndexList[v + 1]]] = triangleIndexList[v + 2];
        /// add neighbours of the third value
        total[triangleIndexList[v + 2]]++;
        neighbourhoods[triangleIndexList[v + 2]][total[triangleIndexList[v + 2]]] = triangleIndexList[v];
        total[triangleIndexList[v + 2]]++;
        neighbourhoods[triangleIndexList[v + 2]][total[triangleIndexList[v + 2]]] = triangleIndexList[v + 1];
    }
#if DEBUG
        printf("Neighbourhood array BEFORE cleaning: \n");
    for (int n = 0; n < nbo; n++){
        printf("nbo %d:  ", n);
        for (int i = 0; i < NAW; i++){
            printf(" %d", neighbourhoods[n][i]);
        }
        printf("\n");
    }
    printf("\n\n");
#endif
        /// remove duplicates
        for (int i = 0; i < nbo; i++) {    /// for each row
            for (int j = 0; j < (total[i]+1); j++) {    /// for each value in the row
                for (int k = j+1; k < (total[i]+1); k++) {    /// for each subsequent value until the end of filled elements
                    if (neighbourhoods[i][j] == neighbourhoods[i][k]) {    /// check if there are pairwise duplicates
                        for (int l = k; l < (total[i]+1); l++) {
                            neighbourhoods[i][l] = neighbourhoods[i][l+1];  /// override duplicate and shift all values left
                            neighbourhoods[i][total[i]] = -1;  /// left most value is converted to empty element
                        }
                        total[i]--; /// decrement pointer for the array to account for this
                        j--;    /// pointer to possible duplicates shifts left one
                    }
                }
            }
        }
#if DEBUG
    printf("Neighbourhood array AFTER cleaning: \n");
    for (int n = 0; n < nbo; n++){
        printf("nbo %d:  ", n);
        for (int i = 0; i < NAW; i++){
            printf(" %d", neighbourhoods[n][i]);
        }
        printf("\n");
    }
    printf("\n\n");
#endif
}

/// repels/attracts points to each other dependent on relative displacement
/// currently only v3 has aliases
void v3CalcSprings(int** neighbourhoods){
    for(int i = 0; i < nbo; i++) { ///for each primary point in pointsArray (iterates through each point using i)
        Point& centre = pointsArray[i]; /// alias for pointsArray[i]
        pointsArray[i].springVec.setZeros(); /// set spring forces to 0

        for (int l = 0; l < NAW; l++) {
            if (neighbourhoods[i][l] != -1){
                Point& neighbour = pointsArray[neighbourhoods[i][l]]; /// alias for pointsArray[neighbourhoods[i][l]]
                /// find the magnitude of distance between the neighbouring point and the central point
                double magnitudeOfDistance = (neighbour.disVec - centre.disVec).magnitude();
                double deltaMagnitude = magnitudeOfDistance - centre.cellRadius;
#if DEBUG
                printf("deltaMag for %d to %d is %f \n", i, (neighbourhoods[i][l]), magnitudeOfDistance);
#endif
                if ((deltaMagnitude > breakSpringCoeff*centre.cellRadius)) {
                    /// do nothing, the connection is ignored (need to show this in graphics somehow)
                }
                else if ((deltaMagnitude > 0)){
                    /// aka point exists outside of the repulsion radius of neighbour it is attracted
                    centre.springVec += (neighbour.disVec - (centre.disVec))
                                        * (deltaMagnitude/magnitudeOfDistance) * centre.extendedHooks;  /// deltaMag/Mag is needed to scale the x component to only that outside the radius of equilibrium
                    neighbour.springVec -= (neighbour.disVec - (centre.disVec))
                                           * (deltaMagnitude/magnitudeOfDistance) * centre.extendedHooks;
                }
                else if ((deltaMagnitude < 0) and (deltaMagnitude > -0.95*centre.cellRadius)){
                    /// aka point exists just within the radius of the neighbouring point and is repelled
                    centre.springVec -= ((neighbour.disVec) - (centre.disVec)) * centre.compressedHooks;
                    neighbour.springVec += ((neighbour.disVec) - (centre.disVec)) * centre.compressedHooks;
                }
                else if ((deltaMagnitude < 0) and (deltaMagnitude < -0.95*centre.cellRadius)) {
                    /// aka point exists just very far within the radius of the neighbouring point and is repelled strongly
                    centre.springVec -= ((neighbour.disVec) - (centre.disVec)) * centre.innerCompressedHooks;
                    neighbour.springVec += ((neighbour.disVec) - (centre.disVec)) * centre.innerCompressedHooks;
                }
            }
        }
    }
}

void v2CalcSprings(){

    /// number of triangle vertices seems to average at 6 per point, setting to 15 for saftey
    /// NAW = neighbourhood array width
    int neighbourhoods[nbo][NAW];
    memset(neighbourhoods, (int)-1, nbo * NAW * sizeof(int));

    int total[nbo];
    memset(total, (int)-1, nbo*sizeof(int));

    /// now to fill the neighbourhood array. nb only goes over triangleIndexList once, rather than once per point
    for (int v = 0; v < numTriangleVertices; v+=3){

        /// add neighbours of the first value of the triangle to its row in the neighbourhood array
        total[triangleIndexList[v]]++;
        neighbourhoods[triangleIndexList[v]][total[triangleIndexList[v]]] = triangleIndexList[v+1];
        total[triangleIndexList[v]]++;
        neighbourhoods[triangleIndexList[v]][total[triangleIndexList[v]]] = triangleIndexList[v+2];

        /// add neighbours of the second value in triangleIndex List to its row in the neighbour array
        total[triangleIndexList[v+1]]++;
        neighbourhoods[triangleIndexList[v+1]][total[triangleIndexList[v+1]]] = triangleIndexList[v];
        total[triangleIndexList[v+1]]++;
        neighbourhoods[triangleIndexList[v+1]][total[triangleIndexList[v+1]]] = triangleIndexList[v+2];

        /// add neighbours of the third value
        total[triangleIndexList[v+2]]++;
        neighbourhoods[triangleIndexList[v+2]][total[triangleIndexList[v+2]]] = triangleIndexList[v];
        total[triangleIndexList[v+2]]++;
        neighbourhoods[triangleIndexList[v+2]][total[triangleIndexList[v+2]]] = triangleIndexList[v+1];
    }


#if DEBUG
    printf("Neighbourhood array BEFORE cleaning: \n");
    for (int n = 0; n < nbo; n++){
        printf("nbo %d:  ", n);
        for (int i = 0; i < NAW; i++){
            printf(" %d", neighbourhoods[n][i]);
        }
        printf("\n");
    }
    printf("\n\n");
#endif

    /// need to remove duplicates from each row, so that each interaction is only present once
    for (int i = 0; i < nbo; i++) {    /// for each row
        for (int j = 0; j < (total[i]+1); j++) {    /// for each value in the row
            for (int k = j+1; k < (total[i]+1); k++) {    /// for each subsequent value until the end of filled elements
                if (neighbourhoods[i][j] == neighbourhoods[i][k]) {    /// check if there are pairwise duplicates
                    for (int l = k; l < (total[i]+1); l++) {
                        neighbourhoods[i][l] = neighbourhoods[i][l+1];  /// override duplicate and shift all values left
                        neighbourhoods[i][total[i]] = -1;  /// left most value is converted to empty element
                    }
                    total[i]--; /// decrement pointer for the array to account for this
                    j--;    /// pointer to possible duplicates shifts left one
                }
            }
        }
    }

#if DEBUG
    printf("Neighbourhood array AFTER cleaning: \n");
    for (int n = 0; n < nbo; n++){
        printf("nbo %d:  ", n);
        for (int i = 0; i < NAW; i++){
            printf(" %d", neighbourhoods[n][i]);
        }
        printf("\n");
    }
    printf("\n\n");
#endif


    for(int i = 0; i < nbo; i++) { ///for each primary point in pointsArray (iterates through each point using i)

        pointsArray[i].springVec.setZeros(); /// set spring forces to 0

        for (int l = 0; l < NAW; l++) {
            if (neighbourhoods[i][l] != -1){
                /// find the magnitude of distance between the neighbouring point and the central point
                double magnitudeOfDistance = (pointsArray[neighbourhoods[i][l]].disVec - pointsArray[i].disVec).magnitude();
                double deltaMagnitude = magnitudeOfDistance - pointsArray[i].cellRadius;
#if DEBUG
                printf("deltaMag for %d to %d is %f \n", i, (neighbourhoods[i][l]), magnitudeOfDistance);
#endif
                if ((deltaMagnitude > 3*pointsArray[i].cellRadius)) {
                /// do nothing, the connection is ignored (need to show this in graphics somehow)
                }
                else if ((deltaMagnitude > 0)){
                /// aka point exists outside of the repulsion radius of neighbour it is attracted
                pointsArray[i].springVec += (pointsArray[neighbourhoods[i][l]].disVec - (pointsArray[i].disVec))
                                          * (deltaMagnitude/magnitudeOfDistance) * pointsArray[i].extendedHooks;  /// deltaMag/Mag is needed to scale the x component to only that outside the radius of equilibrium
                }
                else if ((deltaMagnitude < 0) and (deltaMagnitude > -0.5*pointsArray[i].cellRadius)){
                    /// aka point exists just within the radius of the neighbouring point and is repelled
                    pointsArray[i].springVec -= ((pointsArray[neighbourhoods[i][l]].disVec) - (pointsArray[i].disVec)) * pointsArray[i].compressedHooks;
                }
                else if ((deltaMagnitude < 0) and (deltaMagnitude < -0.1*pointsArray[i].cellRadius)) {
                    /// aka point exists just very far within the radius of the neighbouring point and is repelled strongly
                    pointsArray[i].springVec -=
                            ((pointsArray[neighbourhoods[i][l]].disVec) - (pointsArray[i].disVec)) * 20 *
                            pointsArray[i].compressedHooks;
                }
            }
        }
    }

} /// currently deprecated (needs pairwise interactions & aliases)

void v1CalcSprings(){  /// currently deprecated (needs pairwise interactions & aliases)
    for(int i = 0; i < nbo; i++) { ///for each primary point in pointsArray (iterates through each point using i)

        int pointsConnected[MAX]; /// create an array for the neighbours of a primary point
        int xtotal = 0; /// is a pointer for the pointsConnected array
        pointsArray[i].springVec.setZeros();

        for(int j = 0; j < numTriangleVertices; j +=3){ /// we step through triangleIndexList in 3's
            if((triangleIndexList[j] == i) ||
               (triangleIndexList[j+1] == i) ||
               (triangleIndexList[j+2] == i)){  /// iterates through each triangle to check if any of the vertices are the primary point
                for(int k = 0; k < 3; k++){  /// triangle iterated through using k
                    /// below checks the secondary point isn't the same as the primary point and has not been referenced before
                    if((triangleIndexList[k+j] != i) and (noDuplicateCheck(triangleIndexList[j+k], pointsConnected, xtotal) == true)){
                        pointsConnected[xtotal] = triangleIndexList[j+k];  /// adds the connected points to the array
                        xtotal++;  /// increments the pointer of the pointsConnected Array
                    }
                }
            }
        }

        /// now we've gotten all the connected points we need to change the velocity of each central point in turn
        for (int l = 0; l < xtotal; l++) {
            /// find the magnitude of distance between the neighbouring point and the central point
            /// TODO this stuff is very wrong, ignoring for now
            double magnitudeOfDistance = pointsArray[i].disVec.magnitude();
            double deltaMagnitude = magnitudeOfDistance - pointsArray[i].cellRadius;

            if ((deltaMagnitude > 0)){
                /// aka point exists outside of the repulsion radius of neighbour it is attracted
                pointsArray[i].springVec += (pointsArray[pointsConnected[l]].disVec - (pointsArray[i].disVec))
                                            * (deltaMagnitude/magnitudeOfDistance) * pointsArray[i].extendedHooks;  /// deltaMag/Mag is needed to scale the x component to only that outside the radius of equilibrium
            }
            else if ((deltaMagnitude < 0)){
                /// aka point exists within the radius of the neighbouring point and is repelled
                pointsArray[i].springVec -= ((pointsArray[pointsConnected[l]].disVec) - (pointsArray[i].disVec)) * pointsArray[i].compressedHooks;
            }
        }
    }
}

void iterateDisplace(){
    for(int i = 0; i<nbo; i++){
        pointsArray[i].step();
    }
}

double trackTime(){
    return currentTime += timestep;
}

void startHormone(double inputStartTime){
    static bool flag = false;
    if ((currentTime > inputStartTime) and (flag == false)){
        flag = true;
        /// find the point closest to the hormone Origin
        int closest_point_index = -1;
        double squareMinDist = 1000*1000*xBound;
        for (int i = 0; i < nbo; i++) {
            double squareDisFromOrigin = (pointsArray[i].disVec - hormone1Origin).magnitude_squared();
            if (squareDisFromOrigin < squareMinDist) {
                squareMinDist = squareDisFromOrigin;
                closest_point_index = i;
            }
        }
    /// set this point as the hormone producer
    pointsArray[closest_point_index].isHormoneProducer = true;
    }
    else{
    }
}

/// every cell degrades hormone, only produces produce it
/// this calculates the amount of production / degredation within cells
void hormBirthDeath(double inputStartTime){
    for (int i = 0; i < nbo; i++){
        Point& cell = pointsArray[i]; /// alias for pointsArray[i]
        /// calculate amount of hormone made by producers
        if ((cell.isHormoneProducer == true) and (currentTime < 1000*inputStartTime)){
            printf("Hormone is being produced!\n");
            cell.produceHormone1(hormone1ProdRate);
            cell.degradeHormone1(hormone1DegRate);
        }
        else{
            cell.degradeHormone1(hormone1DegRate);
        }

    }
}

void hormReactDiffuse(double inputStartTime){
    for (int i = 0; i < nbo; i++){
        Point& cell = pointsArray[i]; /// alias for pointsArray[i]
        /// calculate amount of hormone made by producers
        if ((cell.isHormoneProducer == true) and (currentTime < 1000*inputStartTime)){
            printf("Hormone is being produced!\n");
            cell.produceHormone1(hormone1ProdRate);
            cell.produceHormone2(hormone2ProdRate);
            cell.degradeHormone2(hormone2DegRate);
            //cell.react1With2(); /// TODO write member function for this
        }
        else{
            cell.degradeHormone2(hormone2DegRate);
        }

    }
}

void v1DiffuseHorm(int** neighbourhoods) {

    for (int n = 0; n < nbo; n++) {
        pointsArray[n].myDeltaHormone1 = 0;
    }

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
                    double hormoneConcnDiff = centre.myTotalHormone1 - neighbour.myTotalHormone1;  //n / m

                    double hormoneConcnGrad = hormoneConcnDiff / magnitudeOfDistance; //n / m^2
                    /// diffuse the hormone from the centre to neighbour
                    neighbour.myDeltaHormone1 += timestep*(hormone1DiffCoeff * hormoneConcnGrad * centre.cellRadius); //  n = t * (m^2/t * n/m * m)
                    centre.myDeltaHormone1 -= timestep*(hormone1DiffCoeff * hormoneConcnGrad * centre.cellRadius);
                }
            }
        }
    }
    double sum = 0;

    for (int j = 0; j < nbo; j++) {
        Point &cell = pointsArray[j];
        cell.myTotalHormone1 += cell.myDeltaHormone1;
        sum += cell.myTotalHormone1;
        if (cell.myTotalHormone1 < 0) {
            cell.myTotalHormone1 = 0;
        }
    }
printf("The sum of myTotalHormone is %f\n", sum); /// test conservation of hormone
}

int findMaxHormone(){
    int maxPointer = 0;
    for (int i = 0; i<nbo; i++){
        if (pointsArray[i].myTotalHormone1 > pointsArray[maxPointer].myTotalHormone1){
            maxPointer = i;
        }
    }
    return maxPointer;
}

void hormoneExpandEffect(){
    for (int i = 0; i < nbo; i++){
    Point& centre = pointsArray[i];
    centre.cellRadius = centre.cellRadiusBase + (horm1Efficacy * centre.myTotalHormone1 * SCALING_FACTOR);
    }
}

// TODO add a check so that cells cannot divide immediately after dividing again
void calcMitosis(){
    for (int i = 0; i < nbo; i++){
    Point &motherCell = pointsArray[i];
        if (myPrand() < motherCell.divisionProb(baseMaxProbOfDiv, nbo, baseDesiredTotalCells)){

            nbo++; /// MAX points already exist, need to increase pointer by one to access new cell

            Point& daughterCell = pointsArray[nbo-1];
            vector2D normOrient = vector2D(3*mySrand(), mySrand()).normalise();

            vector2D displaceVec = 0.15 * motherCell.cellRadius * normOrient;
            daughterCell.disVec = motherCell.disVec + displaceVec; /// change daughter cell to inherit mother cell position + random orientation
            motherCell.disVec -= displaceVec;  /// mother cell displaced in opposite direction
            // TODO add something that can alter orientation of division
            // TODO maybe add something at causes orientation of division to align with tension
        }
    }
}

void computerDiscreteFourierCoeffs(int iteration, int finalIterationInput){

    if(iteration == finalIterationInput){
        float xx[nbo], yy[nbo];

        for ( int i = 0; i < nbo; i++){
            xx[i] = pointsArray[i].disVec.xx;
            yy[i] = pointsArray[i].disVec.yy;
        }

        double t = 1.0; /// period length
        double sampleFreq = nbo / t;
        double nyquistLim = sampleFreq / 2; /// coefficients representing sample distances greater than 2 per period result in aliasing

        double fourierCoeffs[nbo][2]; /// array of real-valued Fourier coefficients, a read and imaginary component per Coeff

        /// Compute DFT
        for (int k = 0; k < nyquistLim; k++) {
            double sum_re = 0, sum_im = 0;
            for (int n = 0; n < nbo; n++) {
                double angle = 2 * M_PI * k * n / nbo;  /// M_PI is math.h definition of PI to high precision
                sum_re += 2*(xx[n] * cos(angle) + yy[n] * sin(angle)); /// multiplied by 2 to account for nyquist lim
                sum_im += 2*(yy[n] * cos(angle) - xx[n] * sin(angle));
            }
            fourierCoeffs[k][0] = sum_re / nbo;
            fourierCoeffs[k][1] = sum_im / nbo;
        }

        // Print coefficients
        for (int k = 0; k < nyquistLim; k++) {
            double re = fourierCoeffs[k][0];
            double im = fourierCoeffs[k][1];
            printf("c[%d] = %f + %fi\n", k, re, im);
        }
    }
};

int* findAlphaShapePoints(int** neighbourhoods){

    int* concaveHullPointsIndices = create1Darray(nbo);
    init1DArray(concaveHullPointsIndices, nbo, -1);
    int currentHullArrayPointer = 0;


    /// find the point with most negative y co-ord
    int firstPointIndex = -1;
    double currentMinY = yBound;
    for(int i = 0; i < nbo; i++){
        Point& centre = pointsArray[i];
        if(centre.disVec.yy < currentMinY){
            currentMinY = centre.disVec.yy;
            firstPointIndex = i;
            printf("Current Minimum Point is %d with a Y co-ord of %f\n", firstPointIndex, centre.disVec.yy);
        }
    }

    /// fill the array of concaveHullPoints with the initial point's index
    concaveHullPointsIndices[0] = firstPointIndex;
    currentHullArrayPointer++;
    printf("concaveHullPointsIndices[0] =: %d\n", firstPointIndex);

    /// init alias and the vector used for the first comparison
    int currentCentreIndex = firstPointIndex;
    Point &central = pointsArray[currentCentreIndex];
    vector2D initialComparisonVector = vector2D(0, -1);
    int nextCentreIndex = -1;
    double currentMinAngle = 10; /// above max possible value of 2pi

    /// finds next central point in the concave hull
    for(int k =0; k<NAW; k++){
        if (neighbourhoods[firstPointIndex][k] != -1){
            Point &neighbour = pointsArray[neighbourhoods[firstPointIndex][k]]; /// alias for neighbour
            double squareMagDistance = (central.disVec - neighbour.disVec).magnitude_squared();
            if (squareMagDistance < (central.cellRadius*breakSpringCoeff*central.cellRadius*breakSpringCoeff)){
                double angBetweenInitVecAndNeighbour = angleBetweenVecs(initialComparisonVector,
                                                                        (central.disVec - neighbour.disVec));
                if (angBetweenInitVecAndNeighbour < currentMinAngle){
                    nextCentreIndex = neighbourhoods[firstPointIndex][k];
                    currentMinAngle = angBetweenInitVecAndNeighbour;
                }
            }
        }
    }
    concaveHullPointsIndices[1] = nextCentreIndex;
    currentHullArrayPointer++;

    /// change reference points for the new current centre and track the old centre
    int previousCentreIndex = firstPointIndex;
    currentCentreIndex = nextCentreIndex;

    /// carry on finding next point in concave hull
    while(true){

        /// check if the last central point was the same as the first, if so the concave hull has been found
        if (currentCentreIndex == firstPointIndex) {
            printf("Hull finding algorithm has gotten back to the first point\n");
            break;
        }
        else{
            int iterationOfHullFinding = 1; /// for debugging
            currentMinAngle = 10;
            for (int m = 0; m<NAW; m++) {
                Point &neighbour = pointsArray[neighbourhoods[currentCentreIndex][m]];
                Point &previousCentre = pointsArray[previousCentreIndex];
                Point &currentCentre = pointsArray[currentCentreIndex];

                double squareMagDistance = (currentCentre.disVec - neighbour.disVec).magnitude_squared();
                vector2D vectorFromPrevPoint = (currentCentre.disVec - previousCentre.disVec);
                if (neighbourhoods[currentCentreIndex][m] != -1){
                    if (squareMagDistance < (currentCentre.cellRadius * breakSpringCoeff * currentCentre.cellRadius * breakSpringCoeff)) {
                        double angBetweenPrevPointAndNeighbour = angleBetweenVecs(vectorFromPrevPoint,
                                                                                  (neighbour.disVec - currentCentre.disVec));
                        if (angBetweenPrevPointAndNeighbour < currentMinAngle) {
                            previousCentreIndex = currentCentreIndex;
                            currentCentreIndex = neighbourhoods[currentCentreIndex][m];
                            currentMinAngle = angBetweenPrevPointAndNeighbour; // update the value of currentMinAngle
                        }
                    }
                }
            }
            concaveHullPointsIndices[currentHullArrayPointer] = currentCentreIndex;
            currentHullArrayPointer++;
            iterationOfHullFinding++;
            printf("%d points have been included in the concave hull\n", iterationOfHullFinding);
        }
    }

    /*
    printf("Concave Hull is made of points:");
    for(int ii = 0; ii<nbo; ii++){
        while(ii != -1){
            printf("%d", concaveHullPoints[ii]);
        }
    }
    printf("\n");
    */

    return concaveHullPointsIndices;
}

void speedTest(int iterationNumber, int versionOfAlgoUsed, int nboDesired){
    double now =glfwGetTime();
    for (int i = 0; i < iterationNumber; i++)
        {
            create_triangles_list();
            if (versionOfAlgoUsed == 1){
                v1CalcSprings();
                iterateDisplace();
            }
            else if (versionOfAlgoUsed == 2){
                v2CalcSprings();
                iterateDisplace();
            }
            else if (versionOfAlgoUsed == 3){
                int** neighbourhoods = create2Darray(nbo, NAW);
                init2DArray(neighbourhoods, nbo, NAW, -1);
                int* totalArray = create1Darray(nbo);
                init1DArray(totalArray, nbo, -1);
                fill2DArrayNeighbourhoods(neighbourhoods, totalArray, NAW);
                v3CalcSprings(neighbourhoods);
                iterateDisplace();
                free(neighbourhoods);
                free(totalArray);
            }
        }
    double cpu = glfwGetTime() - now;
    printf("Iterations = %d\n Time taken = %f \n", iterationNumber, cpu);

}


/* program entry */
/// argc is the number of arguements, argv    y = yBound * srand(); is pointer to array of strings
int main(int argc, char *argv[]){
    for ( int i=1; i<argc; ++i ) {  /// iterates through arguements
       if ( 0 == readOption(argv[i]) )
           printf("Argument '%s' was ignored\n", argv[i]);
    }
    limitNbo();

    if ( !glfwInit() )
    {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return EXIT_FAILURE;
    }
    glfwSetErrorCallback(error);

    glfwWindowHint(GLFW_DEPTH_BITS, 0);
    //glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GLFW_TRUE);
    //glfwWindowHint(GLFW_CONTEXT_CREATION_API, GLFW_NATIVE_CONTEXT_API);

    GLFWwindow* win = glfwCreateWindow(winW, winH, "FRAP", NULL, NULL);
    if (!win)
    {
        fprintf(stderr, "Failed to open GLFW window\n");
        glfwTerminate();
        return EXIT_FAILURE;
    }
    init(win);

#if BENCHMARK
    for (int i = 1; i < 11; i++) {
        nbo = 100 * i;
        printf("Points to be simulated: %d\n", nbo);
        speedTest(1000, false, 10);
        printf("\n");
    }
#endif

#if DISPLAY
    double next = 0;
    while( !glfwWindowShouldClose(win) )
    {
        static int iterationNumber = 1;
        double now = glfwGetTime();
        if ( now > next)
        {
            while (iterationNumber <= finalIterationNumber){
                iterationNumber++;
                next += delay/100000;
                trackTime();

                printf("nbo is %d\n", nbo);
                create_triangles_list();
                int** neighbourhoods = create2Darray(nbo, NAW);
                init2DArray(neighbourhoods, nbo, NAW, -1);
                int* totalArray = create1Darray(nbo);
                init1DArray(totalArray, nbo, -1);
                fill2DArrayNeighbourhoods(neighbourhoods, totalArray, NAW);
                v3CalcSprings(neighbourhoods);

                iterateDisplace();
                startHormone(hormone1IntroTime);
                v1DiffuseHorm(neighbourhoods);
                hormBirthDeath(hormone1IntroTime);
                calcMitosis();
                //hormoneExpandEffect();

                drawTrianglesAndPoints();
                printf("This is iteration: %d \n\n\n", iterationNumber);
                computerDiscreteFourierCoeffs(iterationNumber, finalIterationNumber);
                int* alphaShapePoints = findAlphaShapePoints(neighbourhoods);
                drawConcaveHull(alphaShapePoints);


                free(triangleIndexList);
                free(neighbourhoods);
                free(totalArray);
                free(alphaShapePoints);
                glfwSwapBuffers(win);
            }
        }
        glfwPollEvents();
    }

    glfwDestroyWindow(win);
    glfwTerminate();
#endif
}
