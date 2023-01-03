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
#include "Clarkson-Delaunay.cpp"  /// this is slightly ood, would be better to compile them seperately and link together (don't know how to do that lol)

/// hard-coded limit to the number of particles
/// size_t = size in bytes,  const means MAX is immutable

const size_t MAX = 16384;
const bool debugStatus = 1;

/// create an array of these and allocate them the max amount of memory possibly required
/// this is global, isn't on the local stack which is good
Point pointsArray[MAX];
int numTriangleVertices = 0;


/// window size in pixels
int winW = 800;
int winH = 800;

///-----------------------------------------------------------------------------

static void error(int error, const char* text)
{
    fprintf(stderr, "GLFW Error: %s\n", text);
}

/// ensure nbo is a multiple of 3 for triangle drawing
void polish()
{
    // limit number of particles
    if ( nbo >= MAX ) nbo = MAX-1;
    printf("The number of points used is %d", nbo);

    //initialize random number generator
    srandom(seed);
}


/// evolves system, steping points forward and accelerating velocity
static void animate()
{
    realTime += delta;
    for ( int i = 0; i < nbo; ++i ) {
        pointsArray[i].step();
        pointsArray[i].accelerate();
    }
}

/// creates an array of xy co-ords for the delaunay triangulation function, then execute it
static WORD* create_triangles_list()
{
    float xyValuesArray[nbo][2];
    for ( int i = 0; i < nbo; i++){
        xyValuesArray[i][0] = pointsArray[i].x;
        xyValuesArray[i][1] = pointsArray[i].y;
    }

    numTriangleVertices = 0;
    WORD* triangleIndexList = NULL;
    triangleIndexList = BuildTriangleIndexList((void*)xyValuesArray, (float)1.0, nbo, (int)2, (int)1, &numTriangleVertices); /// having issues calling the delaunay tessleation function here

#if DEBUG
    printf("\nThere are %d points moving around", nbo);
    printf("\nThe number of vertices defined by numTriangleVertices is %d", numTriangleVertices);
    printf("\ntriangleIndexList contains the values: ");
#endif ///DEBUG

    for (int i = 0; i < numTriangleVertices; i++)
        printf("%u, ", triangleIndexList[i]);
    return triangleIndexList;  /// the triangle index list created is to be used by the draw_triangles function
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
static void draw_triangles(WORD* IndexListofTriangles)
{
    glClear(GL_COLOR_BUFFER_BIT);

    // draw system's edges
    drawSquare(xBound, yBound);

    // draw particles as Triangles:
    for(int i = 0; i < numTriangleVertices; i += 3) { /// iterates through 3 points at a time, needed as it loops back to the start of each polygon
        glBegin(GL_LINE_LOOP);
        for (int j = 0; j <= 2; ++j) {
            pointsArray[IndexListofTriangles[i + j]].displayYellow();
#if DEBUG
            printf("I is %d, J is %d \n", i, j);
            printf("The IndexListofTriangles[i+j] is %d\n", IndexListofTriangles[i + j]);
            printf("The x value of the point is %f\n", pointsArray[i + j].x);
            printf("The y value of the point is %f\n", pointsArray[i + j].y);
#endif ///DEBUG
        }
glEnd();
        for (int k = 0; k < nbo; k += 1) {
            glPointSize(10);
            glBegin(GL_POINTS);
            pointsArray[k].displayWhite();
            glEnd();
        }
    }
#if DEBUG
    glPointSize(6);
    glBegin(GL_POINTS);
    for ( size_t i = 0; i < nbo; ++i )
        pointsArray[i].display();
    glEnd();
    }
#endif ///DEBUG

    printf("draw @ %f\n", realTime);
    glFlush();  // i need to set this to not work if both draw_points and draw_triangles is true

    // EDITS BY FIN
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
    WORD* triangleIndexList;
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
            next += delay/1000;
            animate();
            triangleIndexList = create_triangles_list();
            draw_triangles(triangleIndexList);
            printf("This is iteration: %d \n", interationNumber);
            free(triangleIndexList);
            glfwSwapBuffers(win);
        }
        glfwPollEvents();
    }
    
    glfwDestroyWindow(win);
    glfwTerminate();
}

