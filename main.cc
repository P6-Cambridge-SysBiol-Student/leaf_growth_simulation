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
#include "param.h"
#include "vector.h"
#include "object.h"
#include "polish.h"
#include "Clarkson-Delaunay.cpp"  /// this is slightly odd, would be better to compile them seperately and link together (don't know how to do that lol)


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

/// repels/attracts points to each other dependent on relative displacement
void newCalculateSpringForces(){

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
                double deltaMagnitude = magnitudeOfDistance - repulsionRadius;
#if DEBUG
                printf("deltaMag for %d to %d is %f \n", i, (neighbourhoods[i][l]), magnitudeOfDistance);
#endif
                if ((deltaMagnitude > 0)){
                    /// aka point exists outside of the repulsion radius of neighbour it is attracted
                    pointsArray[i].springVec += (pointsArray[neighbourhoods[i][l]].disVec - (pointsArray[i].disVec))
                                                   * (deltaMagnitude/magnitudeOfDistance) * pointsArray[i].extendedHooks;  /// deltaMag/Mag is needed to scale the x component to only that outside the radius of equilibrium
                }
                else if ((deltaMagnitude < 0)){
                    /// aka point exists within the radius of the neighbouring point and is repelled
                    pointsArray[i].springVec -= ((pointsArray[neighbourhoods[i][l]].disVec) - (pointsArray[i].disVec)) * pointsArray[i].compressedHooks;
                }
            }
        }
    }
}

void oldCalculateSpringForces(){
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
            double deltaMagnitude = magnitudeOfDistance - repulsionRadius;

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

/// draws the square in the window that contains the poinst
void drawSquare(float w, float h){
    glColor3f(0.5, 0.5, 0.5);
    glLineWidth(3);
    glBegin(GL_LINE_LOOP);
    glVertex2f(-w, -h);
    glVertex2f( w, -h);
    glVertex2f( w,  h);
    glVertex2f(-w,  h);
    glEnd();
}

/// draws the points as single points in the window
static void drawPoints(){
    glClear(GL_COLOR_BUFFER_BIT);
    
    // draw system's edges
    drawSquare(xBound, yBound);
    
    // draw particles as points:
    glPointSize(6);
    glBegin(GL_POINTS);
    for ( size_t i = 0; i < nbo; ++i )
        pointsArray[i].displayYellow();
    glEnd();
    
    //printf("draw @ %f\n", realTime);
    glFlush();
}

/// connects 3 consecutive points into triangles, if using the indexed list of points created
/// by the delaunay triangulation it should produce the triangulation
static void drawTrianglesAndPoints(){
    glClear(GL_COLOR_BUFFER_BIT);

    /// draw system's edges
    drawSquare(xBound, yBound);

    /// draw particles as Triangles:
    for(int i = 0; i < numTriangleVertices; i += 3) { /// iterates through 3 points at a time, needed as it loops back to the start of each polygon
        glBegin(GL_LINE_LOOP);
        for (int j = 0; j <= 2; ++j) {
            pointsArray[triangleIndexList[i + j]].displayYellow();

        }
glEnd();
        /// overlay points on top
        for(int k = 0; k < nbo; k += 1) {
            glPointSize(10);
            glBegin(GL_POINTS);
            pointsArray[k].displayWhite();
            glEnd();
        }
    }

    printf("draw @ %f\n", realTime);
    glFlush();
}

static void drawTriangles(){
    glClear(GL_COLOR_BUFFER_BIT);

    // draw system's edges
    drawSquare(xBound, yBound);

    // draw particles as Triangles:
    for(int i = 0; i < numTriangleVertices; i += 3) { /// iterates through 3 points at a time, needed as it loops back to the start of each polygon
        glBegin(GL_LINE_LOOP);
        for (int j = 0; j <= 2; ++j) {
            pointsArray[triangleIndexList[i + j]].displayYellow();

        }
        glEnd();
    }

    printf("draw @ %f\n", realTime);
    glFlush();
}

/// some more graphics stuff
/* change view angle, exit upon ESC */
void key(GLFWwindow* win, int k, int s, int action, int mods){
    if ( action != GLFW_PRESS )
        return;
    
    switch (k)
    {
        case GLFW_KEY_ESCAPE:
            glfwSetWindowShouldClose(win, GLFW_TRUE);
            break;
        case GLFW_KEY_UP:
            break;
        case GLFW_KEY_DOWN:
            break;
        case GLFW_KEY_LEFT:
            break;
        case GLFW_KEY_RIGHT:
            break;
        default:
            return;
    }
}

/* change window size, adjust display to maintain isometric axes */
void reshape(GLFWwindow* win, int W, int H){
    glfwGetWindowSize(win, &winW, &winH);
    //printf("window size %i %i buffer : %i %i\n", winW, winH, W, H);

    pixel = 2 * std::min(xBound/winW, yBound/winH);

    glViewport(0, 0, W, H);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // buffer size:
    double mag = std::min(xBound/W, yBound/H);
    glOrtho(-mag * W, mag * W, -mag * H, mag * H, -1, 1);
    
    //printf("window size %f %f\n", midW, midH);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/* program & OpenGL initialization */
static void init(GLFWwindow* win){
    // Set GLFW callback functions
    glfwSetFramebufferSizeCallback(win, reshape);
    glfwSetKeyCallback(win, key);
    
    glfwMakeContextCurrent(win);
    gladLoadGL(glfwGetProcAddress);
    glfwSwapInterval(1);
    
    int W = winW, H = winH;
    glfwGetFramebufferSize(win, &W, &H);
    reshape(win, W, H);

    // Init OpenGL rendering
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glDisable(GL_DEPTH_TEST);
}

void speedTest(int iterationNumber, bool useOldAlgo, int nboDesired){
    double now =glfwGetTime();
    for (int i = 0; i < iterationNumber; i++)
        {
            create_triangles_list();
            if (useOldAlgo == true){
                oldCalculateSpringForces();
            }
            else{
                newCalculateSpringForces();
            }
            iterateDisplace();
            free(triangleIndexList);
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
    for (int i = 11; i < 16; i++) {
        nbo = 200;
        printf("Points to be simulated: %d\n", nbo);
        speedTest(1000, false, 10);
        printf("\n");
        /// TODO plot these values, take derivative to find exponent (derivatve may also be exponenetial
    }
#endif

#if DISPLAY
    double next = 0;
    while( !glfwWindowShouldClose(win) )
    {
        static int interationNumber = 1;
        double now = glfwGetTime();
        if ( now > next)
        {
            interationNumber++;
            next += delay/100000;
            create_triangles_list();
            newCalculateSpringForces();
            iterateDisplace();
            drawTrianglesAndPoints();
            printf("This is iteration: %d \n", interationNumber);
            free(triangleIndexList);
            glfwSwapBuffers(win);
        }
        glfwPollEvents();
    }

    glfwDestroyWindow(win);
    glfwTerminate();
#endif
}
