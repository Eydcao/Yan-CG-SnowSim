// source: MIT OCW Computer Graphics 6.837
// Creative Common's License

// Modified by Yadi on 4/28/20 for snow simulation

#ifndef EXTRA_H
#define EXTRA_H

// #ifdef WIN32
// #include <windows.h>
// #endif
#include <eigen3/Eigen/Eigen>
using namespace Eigen;
#include <GL/gl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #ifndef M_PI
// #define M_PI 3.14159265358979
// #endif

// Inline functions to help with drawing
inline void glVertex(const Vector3f &a)
{
    float target[3] = {a.x(), a.y(), a.z()};
    glVertex3fv(target);
}

inline void glNormal(const Vector3f &a)
{
    float target[3] = {a.x(), a.y(), a.z()};
    glNormal3fv(target);
}

inline void glLoadMatrix(const Matrix4f &m)
{
    float target[16] = {m(0, 0), m(1, 0), m(2, 0), m(3, 0), m(0, 1), m(1, 1),
                        m(2, 1), m(3, 1), m(0, 2), m(1, 2), m(2, 2), m(3, 2),
                        m(0, 3), m(1, 3), m(2, 3), m(3, 3)};
    glLoadMatrixf(target);
}

inline void glMultMatrix(const Matrix4f &m)
{
    float target[16] = {m(0, 0), m(1, 0), m(2, 0), m(3, 0), m(0, 1), m(1, 1),
                        m(2, 1), m(3, 1), m(0, 2), m(1, 2), m(2, 2), m(3, 2),
                        m(0, 3), m(1, 3), m(2, 3), m(3, 3)};
    glMultMatrixf(target);
}

#endif
