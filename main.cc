/*
 * Elementary simulation using GLFW + OpenGL for display modified for use in PlantSim by Finley Webb
 * Original code by Francois J Nedelec, Cambridge University, 13 Nov 2021, 11 Oct 2022
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define DEBUG TRUE
#define GLAD_GL_IMPLEMENTATION
#include <glad/gl.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <vector>

#include "param.h"
#include "random.h"
#include "object.h"
#include "polish.h"
#include "Clarkson-Delaunay.cpp"  /// this is slightly odd, would be better to compile them seperately and link together (don't know how to do that lol)


///-----------------------------------------------------------------------------

static void error(int error, const char* text)
{
    fprintf(stderr, "GLFW Error: %s\n", text);
}

void polish()
{
    // limit number of particles
    if ( nbo >= MAX ) nbo = MAX-1;
    printf("The number of points used is %d", nbo);

    //initialize random number generator
    srandom(seed);
}


/// evolves system, stepping points forward and accelerating velocity
static void animate()
{
    realTime += delta;
    for ( int i = 0; i < nbo; ++i ) {
        pointsArray[i].step();
        /// pointsArray[i].accelerateToCentre();
    }
}

/// creates an array of xy co-ords for the delaunay triangulation function, then execute it
void create_triangles_list()
{
    float xyValuesArray[nbo][2];
    for ( int i = 0; i < nbo; i++){
        xyValuesArray[i][0] = pointsArray[i].x;
        xyValuesArray[i][1] = pointsArray[i].y;
    }

    numTriangleVertices = 0;
    triangleIndexList = BuildTriangleIndexList((void*)xyValuesArray, (float)1.0, nbo, (int)2, (int)1, &numTriangleVertices); /// having issues calling the delaunay tessleation function here

#if DEBUG
    printf("\nThere are %d points moving around", nbo);
    printf("\nThe number of vertices defined by numTriangleVertices is %d", numTriangleVertices);
    printf("\ntriangleIndexList contains the values: ");
#endif ///DEBUG

    for (int i = 0; i < numTriangleVertices; i++)
        printf("%u, ", triangleIndexList[i]);
}


bool noDuplicateCheck(int indexValueToCheck, int arrayToCheck[], int max) /// checks through the pointsConnected array to see if secondary point is present
{
    static bool found = false;

    for (int i = 0; i <= max; i++) {  /// going up to total instead of over all the array is faster
        if (arrayToCheck[i] == indexValueToCheck) {
            found = true;
            break;
        }
    }

    if (found = false){
        return true;
    }
    else{
        return false;
    }
}

/// finds all the points connected to the primary point through delaunay triangulations (IS REDUNDANT)
int* connectedSecondaryPoints(int indexOfInterest, int numberOfTriangleVertices)
{
    static int pointsConnected[MAX];  /// important to be static, so it can be called in pointAttractOrRepel function. Cannot be nbo sized as this is dynamic
    static int total = 0; /// is a pointer for the pointsConnected array

    for(int i = 0; i < numTriangleVertices; i +=3){
        if((triangleIndexList[i] == indexOfInterest) || \
           (triangleIndexList[i+1] == indexOfInterest) || \
           (triangleIndexList[i+2] == indexOfInterest)){  /// iterates through each triangle in turn
            for(int j = 0; j < j+3; j++){  /// iterates through each point the current triangle
                /// below checks the secondary point isn't the same as the primary point and has not been referenced before
                if((triangleIndexList[j] != indexOfInterest) and (noDuplicateCheck(triangleIndexList[j], pointsConnected, total) == true)){
                    pointsConnected[total] = triangleIndexList[j];  /// adds the connected points to the array
                    total++;
                }
            }
        }
    }

    return pointsConnected;
}

/// repels/attracts points to each other dependent on relative displacement (CURRENTLY PSEUDOCODE)
void pointsAttractOrRepel()
{
    for(int i = 0; i < nbo; i++) { ///for each primary point in pointsArray (a primary point is identified as i)

        /// create an array of all neighbouring points (a secondary point is identified as j)
        int pointsConnected[MAX];  /// important to be static, so it can be called in pointAttractOrRepel function. Cannot be nbo sized as this is dynamic
        int total = 0; /// is a pointer for the pointsConnected array

        for(int j = 0; j < numTriangleVertices; j +=3){ /// we step through triangleIndexList in 3's
            if((triangleIndexList[j] == i) || \
           (triangleIndexList[j+1] == i) || \
           (triangleIndexList[j+2] == i)){  /// iterates through each triangle to check if any of the vertices are out primary point
                for(int k = 0; k < (j+3) ; k++){  /// iterates through each point the current triangle if our primary point is found
                    /// below checks the secondary point isn't the same as the primary point and has not been referenced before
                    if((triangleIndexList[k] != i) and (noDuplicateCheck(triangleIndexList[j+k], pointsConnected, total) == true)){
                        pointsConnected[total] = triangleIndexList[j+k];  /// adds the connected points to the array
                        total++;  /// increments the pointer of the poitnsConnected Array
                    }
                }
            }
        }

        /// now we've gotten all the connected points we need to change the velocities of the central point
        for (int j = 0; j <= total; j++) {
            double magnitudeOfDistance = sqrt((pow((pointsArray[pointsConnected[j]].x) - (pointsArray[i].x), 2)\
                                               + pow((pointsArray[pointsConnected[j]].y) - (pointsArray[i].y), 2) ));
            double deltaMagnitude = magnitudeOfDistance - repulsionRadius;
            if (deltaMagnitude > 0){ /// aka point exists outside of the repulsion radius and are attracted
                pointsArray[i].xvelocity += (pointsArray[pointsConnected[j]].x) * (deltaMagnitude/magnitudeOfDistance) * pointsArray[i].hooks;  /// deltaMag/Mag is needed to scale the x component to only that outside the radius of equilibrium
                pointsArray[i].yvelocity += (pointsArray[pointsConnected[j]].y) * (deltaMagnitude/magnitudeOfDistance) * pointsArray[i].hooks;
            }
            else{ /// aka point exists within the radius and is repelled
                pointsArray[i].xvelocity -= (pointsArray[pointsConnected[j]].x) - (pointsArray[i].x) * pointsArray[i].hooks;
                pointsArray[i].yvelocity -= (pointsArray[pointsConnected[j]].y) - (pointsArray[i].y) * pointsArray[i].hooks;
            }
        }
    }
}


void drawSquare(float w, float h)  /// draws the square in the window that contains the poinst
{
    glColor3f(0.5, 0.5, 0.5);
    glLineWidth(3);
    glBegin(GL_LINE_LOOP);
    glVertex2f(-w, -h);
    glVertex2f( w, -h);
    glVertex2f( w,  h);
    glVertex2f(-w,  h);
    glEnd();
}

/* draw System */
static void draw()  /// draws the points as single points in the window
{
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


/// connects 3 consective points into triangles, if using the indexed list of points created
/// by the delaunay triangulation it should produce the triangulation
static void drawTrianglesAndPoints()
{
    glClear(GL_COLOR_BUFFER_BIT);

    /// draw system's edges
    drawSquare(xBound, yBound);

    /// draw particles as Triangles:
    for(int i = 0; i < numTriangleVertices; i += 3) { /// iterates through 3 points at a time, needed as it loops back to the start of each polygon
        glBegin(GL_LINE_LOOP);
        for (int j = 0; j <= 2; ++j) {
            pointsArray[triangleIndexList[i + j]].displayYellow();
#if DEBUG
            printf("I is %d, J is %d \n", i, j);
            printf("The IndexListofTriangles[i+j] is %d\n", IndexListofTriangles[i + j]);
            printf("The x value of the point is %f\n", pointsArray[i + j].x);
            printf("The y value of the point is %f\n", pointsArray[i + j].y);
#endif ///DEBUG
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

static void drawTriangles()
{
    glClear(GL_COLOR_BUFFER_BIT);

    // draw system's edges
    drawSquare(xBound, yBound);

    // draw particles as Triangles:
    for(int i = 0; i < numTriangleVertices; i += 3) { /// iterates through 3 points at a time, needed as it loops back to the start of each polygon
        glBegin(GL_LINE_LOOP);
        for (int j = 0; j <= 2; ++j) {
            pointsArray[triangleIndexList[i + j]].displayYellow();

#if DEBUG
            printf("I is %d, J is %d \n", i, j);
            printf("The IndexListofTriangles[i+j] is %d\n", IndexListofTriangles[i + j]);
            printf("The x value of the point is %f\n", pointsArray[i + j].x);
            printf("The y value of the point is %f\n", pointsArray[i + j].y);
#endif ///DEBUG

        }
        glEnd();
    }

    printf("draw @ %f\n", realTime);
    glFlush();
}

/// some more graphics stuff
/* change view angle, exit upon ESC */
void key(GLFWwindow* win, int k, int s, int action, int mods)
{
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
void reshape(GLFWwindow* win, int W, int H)
{
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
static void init(GLFWwindow* win)
{
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


/* program entry */
/// argc is the number of arguements, argv    y = yBound * srand(); is pointer to array of strings
int main(int argc, char *argv[])
{
    for ( int i=1; i<argc; ++i ) {  /// iterates through arguements
       if ( 0 == readOption(argv[i]) )
           printf("Argument '%s' was ignored\n", argv[i]);
    }
    printf("I am here \n");
    polish();
    
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


    double next = 0;
    while( !glfwWindowShouldClose(win) )
    {
        static int interationNumber = 1;
        double now = glfwGetTime();
        if ( now > next)
        {
            interationNumber++;
            next += delay/100000;
            pointsAttractOrRepel();
            animate();
            create_triangles_list();
            drawTrianglesAndPoints();
            printf("This is iteration: %d \n", interationNumber);
            free(triangleIndexList);
            glfwSwapBuffers(win);
        }
        glfwPollEvents();
    }
    
    glfwDestroyWindow(win);
    glfwTerminate();
}

