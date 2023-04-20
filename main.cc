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
#define REGULAR_LATTICE false
#define MOVING_POINTS true
#define GLAD_GL_IMPLEMENTATION
#include <glad/gl.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <vector>

#include "random.h"
#include "vector.h"
#include "param.h"
#include "object.h"
#include "polish.h"
#include "arrays.h"
#include "hormone.h"
#include "Clarkson-Delaunay.cpp"
#include "springs.h"
#include "graphics.h"
#include "fitness.h"


///-----------------------------------------------------------------------------

/// evolves system, stepping points forward and accelerating velocity
static void animate(){
    realTime += delta;
    for ( int i = 0; i < nbo; ++i ) {
        pointsArray[i].step();
    }
}

void initRegularTriangularLattice() {
    int index = 0;
    bool isSqrtNBOWhole = !fmod(sqrt(nbo), 1);

    int numPointsX = sqrt(nbo);
    int numPointsY = sqrt(nbo);
    double spacing = pointsArray[0].cellRadius*2;
    double xSum = 0.0, ySum = 0.0;
    int numPoints = numPointsX * numPointsY;
    if (isSqrtNBOWhole){
        for (int i = 0; i < numPointsX; i++) {
            for (int j = 0; j < numPointsY; j++) {
                double x = i * spacing + ((j % 2 == 0) ? 0 : spacing / 2.0);
                double y = j * spacing * sin(M_PI / 3.0);
                Point& p = pointsArray[index];
                p.disVec = vector2D(x, y);
                xSum += x;
                ySum += y;
                index++;
            }
        }
        double xCenter = xSum / numPoints;
        double yCenter = ySum / numPoints;
        for (int i = 0; i < numPoints; i++) {
            pointsArray[i].disVec -= vector2D(xCenter, yCenter);
        }
    }
    else{
        printf("NBO DOES NOT EQUAL numPointsX * numPointsY\n");
    }
}

void initPerfectCircle(double circleRadius) {
    int index = 0;
    double angleSpacing = 2 * M_PI / nbo;
    double xSum = 0.0, ySum = 0.0;

    for (int i = 0; i < nbo; i++) {
        double angle = 2 * i * angleSpacing;
        double x = circleRadius * cos(angle);
        double y = circleRadius * sin(angle);
        Point& p = pointsArray[index];
        p.disVec = vector2D(x, y);
        xSum += x;
        ySum += y;
        index++;
    }

    double xCenter = xSum / nbo;
    double yCenter = ySum / nbo;
    for (int i = 0; i < nbo; i++) {
        pointsArray[i].disVec -= vector2D(xCenter, yCenter);
    }
}

void initHollowSquare(double sideLength, int nbo) {
    int pointsPerSide = nbo / 4;
    double spacing = sideLength / (pointsPerSide);
    double xSum = 0.0, ySum = 0.0;

    for (int i = 0; i < nbo; i++) {
        int side = i / pointsPerSide;
        int j = i % pointsPerSide;
        double x, y;

        if (side == 0) {
            x = j * spacing - sideLength / 2;
            y = -sideLength / 2;
        } else if (side == 1) {
            x = sideLength / 2;
            y = j * spacing - sideLength / 2;
        } else if (side == 2) {
            x = sideLength / 2 - j * spacing;
            y = sideLength / 2;
        } else { // side == 3
            x = -sideLength / 2;
            y = sideLength / 2 - j * spacing;
        }

        Point& p = pointsArray[i];
        p.disVec = vector2D(x, y);
        xSum += x;
        ySum += y;
    }

    double xCenter = xSum / nbo;
    double yCenter = ySum / nbo;
    for (int i = 0; i < nbo; i++) {
        pointsArray[i].disVec -= vector2D(xCenter, yCenter);
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
                    }
                    total[i]--; /// decrement pointer for the array to account for this
                    j--;    /// pointer to possible duplicates shifts left one
                    neighbourhoods[i][total[i] + 1] = -1;  /// set the correct value to -1 (empty element)
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

void iterateDisplace(){
    for(int i = 0; i<nbo; i++){
        pointsArray[i].step();
    }
}

double trackTime(){
    return currentTime += timestep;
}

// TODO add a check so that cells cannot divide immediately after dividing again
void calcMitosis(){
    for (int i = 0; i < nbo; i++){
        Point &motherCell = pointsArray[i];
        if (myPrand() < motherCell.divisionProb(baseMaxProbOfDiv, nbo, DesiredTotalCells)){

            nbo++; /// MAX points already exist, need to increase pointer by one to access new cell

            Point& daughterCell = pointsArray[nbo-1]; /// nbo-1 as C is 0 indexed
            vector2D OrientVec = vector2D(mySrand(), mySrand())
                               + (motherCell.myTotalHormone1*horm1Efficacy*horm1DivOrient)
                               + (motherCell.myTotalHormone2*horm2Efficacy*horm2DivOrient);
            vector2D normOrient = OrientVec.normalise();

            vector2D displaceVec = 0.5 * motherCell.cellRadius * normOrient;
            daughterCell.disVec = motherCell.disVec + displaceVec; /// change daughter cell to inherit mother cell position + random orientation
            motherCell.disVec -= displaceVec;  /// mother cell displaced in opposite direction
        }
    }
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
int main(int argc, char *argv[]) {
    bool cym_file_found = false;

    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        size_t n = strlen(arg);
        if (n > 4 && strcmp(arg + n - 4, ".cym") == 0) {
            cym_file_found = true;
            readFile(arg);
            break;
        }
    }

    if (!cym_file_found) {
        printf(".cym file not found\n Using Defaults\nF");
    }

    if (!glfwInit()) { // Call glfwInit() before using any other GLFW functions
        fprintf(stderr, "Failed to initialize GLFW\n");
        return EXIT_FAILURE;
    }
    glfwSetErrorCallback(error);

    glfwWindowHint(GLFW_DEPTH_BITS, 0);
    //glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GLFW_TRUE);
    //glfwWindowHint(GLFW_CONTEXT_CREATION_API, GLFW_NATIVE_CONTEXT_API);

#if DISPLAY
    GLFWwindow *win = glfwCreateWindow(winW, winH, "LifeSim", NULL, NULL);
    if (!win) {
        fprintf(stderr, "Failed to open GLFW window\n");
        glfwTerminate();
        return EXIT_FAILURE;
    }
    init(win);
#endif
#if BENCHMARK
    for (int i = 1; i < 11; i++) {
        nbo = 100 * i;
        printf("Points to be simulated: %d\n", nbo);
        speedTest(1000, false, 10);
        printf("\n");
    }
#endif
    double next = 0;
    bool shouldTerminate = false;
#if DISPLAY
    while (!glfwWindowShouldClose(win) && !shouldTerminate) {
#else
        while (!shouldTerminate) {
#endif
        static int iterationNumber = 1;
        double now = glfwGetTime();
        if (now > next) {
            static double currentTime = 0;
            while (currentTime <= finalTime + timestep) {
                currentTime += timestep;
#if REGULAR_LATTICE
                if (iterationNumber == 1) {
                    //initPerfectCircle(20*SCALING_FACTOR);
                    //initHollowSquare(20 * SCALING_FACTOR, nbo);
                    initRegularTriangularLattice();
                }
#endif
#if DISPLAY
                glfwPollEvents();
                glClear(GL_COLOR_BUFFER_BIT);
#endif
                iterationNumber++;
                next += delay / 100000;
                trackTime();printf("%d cells exist\n", nbo);
                printf("Current time = %f\n", currentTime);
                calcMitosis();

                create_triangles_list();
                int **neighbourhoods = create2Darray(nbo, NAW); /// malloc empty nbo * NAW array
                init2DArray(neighbourhoods, nbo, NAW, -1); /// fill it with -1s
                int *totalArray = create1Darray(nbo); /// create empty nbo array
                init1DArray(totalArray, nbo, -1);  /// fill it with -1s
                fill2DArrayNeighbourhoods(neighbourhoods, totalArray, NAW); /// fill neighbourhood aray

#if MOVING_POINTS
                v3CalcSprings(neighbourhoods);
#endif
                iterateDisplace();
                calcHormBirthDeath();
                v1DiffuseHorm(neighbourhoods);
                hormReactDiffuse(hormone2IntroTime);
                globalUpdateHormone();
                double globalHorm2 = sumHormone2();
                double globalHorm1 = sumHormone1();
                //printf("Global Horm 2 is %f: Average Horm2 is %f\n", globalHorm2, globalHorm2/nbo);

                if ((currentTime > hormone2IntroTime) and (isnan(globalHorm2) or isnan(globalHorm1))){
                    shouldTerminate = true;
                    printf("A Hormone Has Exploded!\n");
                    break;
                }

#if DISPLAY
                double maxHormone2 = findMaxHormone2();
                double minHormone2 = findMinHormone2();
                drawPointsHorm2(maxHormone2); // calls
#endif
                free(triangleIndexList);
                free(totalArray);
                for (int i = 0; i < nbo; i++) {
                    free(neighbourhoods[i]);
                }
                free(neighbourhoods);

                if (currentTime >= finalTime) {
                    int fourierCoeffsNum = 0.5*nbo;
                    if (nbo > 2*maxFourierCoeffs){
                        fourierCoeffsNum = maxFourierCoeffs;
                    }
                    double **fourierCoeffs = computeDeltaFourierCoeffs(fourierCoeffsNum);
                    outputFourierToFile(fourierCoeffs, fourierCoeffsNum, "outputFourierCoeffs.csv");
#if DISPLAY
                    if (displayInverseFourier) {
                        glfwPollEvents();
                        reconstructShape(fourierCoeffs, fourierCoeffsNum);
                    }
#endif
                    for (int i = 0; i < fourierCoeffsNum; i++){
                        free(fourierCoeffs[i]);
                    }
                    free(fourierCoeffs);

                    printf("Fourier Coefficients Saved!\n");
                    break;

                }

                free(out_of_flat_p_neigh.basis); /// not freed in delaunay.cpp by default, causes memory leak if not
                shouldTerminate = true;
#if DISPLAY

                glFlush();
                glfwSwapBuffers(win);
#endif

            }
        }
#if DISPLAY
        if (glfwGetKey(win, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(win, GLFW_TRUE);
        }
#endif
    }
#if DISPLAY
    glfwDestroyWindow(win);
    glfwTerminate();
#endif
}
