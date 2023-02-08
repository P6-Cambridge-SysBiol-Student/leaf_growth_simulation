//
// Created by finley on 01/02/23.
//
#include <math.h>
#include "random.h"
#ifndef FRAP_VECTOR_H
#define FRAP_VECTOR_H

#endif //FRAP_VECTOR_H

/// TODO change variables to conventions
/// TODO maybe use a typedef to define size of elements
/// want it to be efficient
typedef double elem_type;

class vector2D
{
public:
    elem_type xx, yy; /// x and y components of the vector

    vector2D() : xx(0), yy(0) {}  /// default is to set x and y to 0
    vector2D(elem_type inputX, elem_type inputY) : xx(inputX), yy(inputY) {} /// constructor, sets x and y to the specified values if desired

    /// overload the + operator to add two vectors
    vector2D operator+(const vector2D& vec) /// pointer used to pass vector by reference, function automatiically derefernce the address?
    {
        return vector2D(xx + vec.xx, yy + vec.yy);
        /// v1 = (1, 2)  v2 = (3, 4)  v3 = v1 + v2; is now value (v3 = 3, 6)
    }

    /// overload the += operator to increment on a vector
    void operator+=(const vector2D& vec)
    {
        xx += vec.xx;
        yy += vec.yy;
        /// *this is a pointer to the current object being operated on
        /// v4 = (5, 6)  v5 = (7, 8)   v4 += v5;  v4 now evaluates to (12, 14)
    }

    /// overload the - operator to subtract two vectors
    vector2D operator-(const vector2D& vec)
    {
        return vector2D(xx - vec.xx, yy - vec.yy);
        /// v4 = (5, 6)  v5 = (7, 8)   v6 = v4 - v5;  v6 now evaluates to (-2, -2)
    }

    vector2D& operator-=(const vector2D& vec)
    {
        xx -= vec.xx;
        yy -= vec.yy;
        return *this;
    }

    /// overload the * operator to multiply a vector by a scalar
    /// here vec is the vector in the operator, no matter what side it is on
    friend vector2D operator*(vector2D& vec, const double scalar)
    {
        return vector2D(vec.xx * scalar, vec.yy * scalar);
    }

    void operator*=(double const scalar)
    {
        xx *= scalar;
        yy *= scalar;
    }

    /// overload the / operator to divide a vector by a scalar
    vector2D operator/(double const scalar)
    {
        return vector2D(xx / scalar, yy / scalar);
    }

    /// overlaod the /= operator
    void operator/=(double const scalar)
    {
        xx /= scalar;
        yy /= scalar;
    }

    /// returns the magnitude (length) of the vector
    double magnitude_squared() const
    {
        return (xx * xx + yy * yy);
    }

    double magnitude() const
    {
        return (xx * xx + yy * yy);
    }

    /// normalises the vector, i.e. scales it to have a magnitude of 1
    void normalise()
    {
        double mag = magnitude();
        xx /= mag;
        yy /= mag;
    }
};