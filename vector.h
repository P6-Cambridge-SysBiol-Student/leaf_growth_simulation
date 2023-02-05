//
// Created by finley on 01/02/23.
//
#include <math.h>
#ifndef FRAP_VECTOR_H
#define FRAP_VECTOR_H

#endif //FRAP_VECTOR_H

class vector2D
{
public:
    double x, y; /// x and y components of the vector

    vector2D() : x(xBound*srand()), y(xBound*srand()) {} /// default constructor, sets x and y to random positions

    vector2D(double x, double y) : x(x), y(y) {} /// constructor, sets x and y to the specified values

    /// overload the + operator to add two vectors
    vector2D operator+(const vector2D &v) const
    {
        return vector2D(x + v.x, y + v.y);
    }

    /// overload the - operator to subtract two vectors
    vector2D operator-(const vector2D &v) const
    {
        return vector2D(x - v.x, y - v.y);
    }

    /// overload the * operator to multiply a vector by a scalar
    vector2D operator*(double k) const
    {
        return vector2D(k * x, k * y);
    }

    /// returns the magnitude (length) of the vector
    double magnitude() const
    {
        return sqrt(x * x + y * y);
    }

    /// normalizes the vector, i.e. scales it to have a magnitude of 1
    void normalize()
    {
        double mag = magnitude();
        x /= mag;
        y /= mag;
    }
};