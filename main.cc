/*
 * Elementary simulation using GLFW + OpenGL for display
 * Francois J Nedelec, Cambridge University, 13 Nov 2021, 11 Oct 2022
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define GLAD_GL_IMPLEMENTATION
#include <glad/gl.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <vector>

#include <GL/freeglut_std.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "param.h"
#include "random.h"
#include "object.h"
#include "delaunay.hpp"


// hard-coded limit to the number of particles
const size_t MAX = 16384;

Object obj[MAX];

// window size in pixels
int winW = 800;
int winH = 800;

//-----------------------------------------------------------------------------

static void error(int error, const char* text)
{
    fprintf(stderr, "GLFW Error: %s\n", text);
}

// calculate derived parameters
void polish()
{
    // limit number of particles:
    if ( nbo >= MAX ) nbo = MAX-1;
    // calibrate diffusion:
    alpha = sqrt( 2 * diff * delta );
    //initialize random number generator
    srandom(seed);
}

/* evolve System */
static void animate()
{
    realTime += delta;
    for ( int i = 0; i < nbo; ++i )
        obj[i].step();
}

void drawSquare(float w, float h)
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
static void draw()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
    // draw system's edges
    drawSquare(xBound, yBound);
    
    // draw particles as points:
    glPointSize(20);
    glBegin(GL_POINTS);
    for ( size_t i = 0; i < nbo; ++i )
        obj[i].display();
    glEnd();
    
    //printf("draw @ %f\n", realTime);
    glFlush();
}

static void draw_triangles()
{
    glClear(GL_COLOR_BUFFER_BIT);

    // draw system's edges
    drawSquare(xBound, yBound);

    // draw particles as Triangles:
    glPointSize(8);
    int i;  // not quite sure why it wanted me to initialise i here before using it in the loop
    for(i = 0; i < nbo; i += 3) {
        glBegin(GL_TRIANGLES);
        for (size_t j = 0; j <= 2; ++j)  // especially when it didn't need me to initialise j
            obj[(i+j)].display();
        glEnd();
    }
    printf("draw @ %f\n", realTime);
    glFlush();  // i need to set this to not work if both draw_points and draw_triangles is true

    // EDITS BY FIN
}

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

using namespace delaunay;

namespace context {
    std::vector<Point<float>> points;
} /* namespace context */



/* program entry */
int main(int argc, char *argv[])
{
    for ( int i=1; i<argc; ++i ) {
       if ( 0 == readOption(argv[i]) )
           printf("Argument '%s' was ignored\n", argv[i]);
    }
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
        double now = glfwGetTime();
        if ( now > next )
        {
            next += 0.05; // will give 20 frames/second
            animate();
            draw_triangles();
            glfwSwapBuffers(win);
        }
        glfwPollEvents();
    }
    
    glfwDestroyWindow(win);
    glfwTerminate();
}

