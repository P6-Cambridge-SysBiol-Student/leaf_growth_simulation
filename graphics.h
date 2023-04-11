//
// Created by finley on 26/02/23.
//

#ifndef FRAP_GRAPHICS_H
#define FRAP_GRAPHICS_H
#include <GLFW/glfw3.h>

#endif //FRAP_GRAPHICS_H

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
static void drawPointsHorm1(int inputMaxHormone){
    glClear(GL_COLOR_BUFFER_BIT);

    // draw system's edges
    drawSquare(xBound, yBound);

    // draw particles as points:
    glPointSize(3);
    glBegin(GL_POINTS);
    for(int k = 0; k < nbo; k += 1) {
        pointsArray[k].linearDisplayHormone1(inputMaxHormone);
    }
    glEnd();

    //printf("draw @ %f\n", realTime);
}

/// draws the points as single points in the window
static void drawPointsHorm2(double inputMaxHormone){
    glClear(GL_COLOR_BUFFER_BIT);

    // draw system's edges
    drawSquare(xBound, yBound);

    // draw particles as points:
    glPointSize(8);
    glBegin(GL_POINTS);
    for(int k = 0; k < nbo; k += 1) {
        pointsArray[k].linearDisplayHormone2(inputMaxHormone);
    }
    glEnd();

    //printf("draw @ %f\n", realTime);
}



/// connects 3 consecutive points into triangles, if using the indexed list of points created
/// by the delaunay triangulation it should produce the triangulation
static void drawTrianglesAndPoints(double inputMaxHormone){
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
        glPointSize(12);
        glBegin(GL_POINTS);
        for(int k = 0; k < nbo; k += 1) {
            pointsArray[k].linearDisplayHormone1(inputMaxHormone);
        }
        glEnd();
    }

    printf("\ndraw @ %f\n", realTime);
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

static void drawConcaveHull(int* inputConcaveHullArray){
    for(int i = 0; i < nbo; i++){
        glBegin(GL_LINE_LOOP);
        Point& p = pointsArray[inputConcaveHullArray[i]];
        p.displayGreen(); // call displayGreen() using an object of the Point class
    }
    glEnd();
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
    glClearColor(1, 1, 1, 1);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glDisable(GL_DEPTH_TEST);
}