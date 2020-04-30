//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_INTERSECTION_H
#define RAYTRACING_INTERSECTION_H
#include "Shape.hpp"
#include <eigen3/Eigen/Eigen>
#include <limits>

using namespace Eigen;

class Shape;

struct Intersection
{
    Intersection()
    {
        happened = false;
        coords = Vector3f();
        normal = Vector3f();
        distance = std::numeric_limits<double>::max();
        shape = nullptr;
        // m = nullptr;
    }
    bool happened;
    Vector3f coords;
    // Vector3f tcoords;
    Vector3f normal;
    // Vector3f emit;
    double distance;
    Shape* shape;
    // Material* m;
};
#endif  // RAYTRACING_INTERSECTION_H
