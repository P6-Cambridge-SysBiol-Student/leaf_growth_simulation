//
// Created by finley on 03/03/23.
//

#ifndef FRAP_SPRINGS_H
#define FRAP_SPRINGS_H

#endif //FRAP_SPRINGS_H


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
