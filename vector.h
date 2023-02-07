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
    vector2D(double x, double y) : x(x), y(y) {} /// constructor, sets x and y to the specified values if desired

    /// overload the + operator to add two vectors
    vector2D operator+(const vector2D &v) /// all these overloadings used to have const, removed
    {
        return vector2D(x + v.x,
                        y + v.y);
        /// v1 = (1, 2)  v2 = (3, 4)  v3 = v1 + v2; is now value (v3 = 3, 6)
    }

    /// overload the += operator to increment on a vector
    vector2D& operator+=(vector2D &v)
    {
        x += v.x;
        y += v.y;
        return *this; /// *this is a pointer to the current object being operated on
        /// v4 = (5, 6)  v5 = (7, 8)   v4 += v5;  v4 now evaluates to (12, 14)
    }

    /// overload the - operator to subtract two vectors
    vector2D operator-(const vector2D &v)
    {
        return vector2D(x - v.x,
                        y - v.y);
        /// v4 = (5, 6)  v5 = (7, 8)   v6 = v4 - v5;  v6 now evaluates to (-2, -2)
    }

    vector2D& operator-=(const vector2D &v)
    {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    /// overload the * operator to multiply a vector by a scalar
    vector2D operator*(double k)
    {
        return vector2D(k * x,
                        k * y);
    }

    vector2D operator*=(double k)
    {
        return vector2D(x *= k, y *= k);
    }

    /// overload the / operator to divide a vector by a scalar
    vector2D operator/(double k)
    {
        return vector2D(x/k, y/k);
    }

    /// overlaod the /= operator
    vector2D operator/=(double k)
    {
        x/= k;
        y/= k;
        return *this;
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