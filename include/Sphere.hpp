//
// Created by LEI XU on 5/13/19.
//

// Modified by Yadi on 4/28/20 for snow simulation

#ifndef SNOWSIM_SPHERE
#define SNOWSIM_SPHERE

#include "Bounds3.hpp"
#include "Shape.hpp"

class Sphere : public Shape
{
   public:
    Vector3f center;
    float radius, radius2;
    Sphere(const Vector3f &c, const float &r)
        : center(c), radius(r), radius2(r * r)
    {
    }
    Bounds3 getBounds()
    {
        return Bounds3(Vector3f(center.x() - radius, center.y() - radius,
                                center.z() - radius),
                       Vector3f(center.x() + radius, center.y() + radius,
                                center.z() + radius));
    }
    int generateParticlesInside(int lNumberDensity,
                                std::vector<Vector3f> &contain)
    {
        // TODO generate particles pos inside
        contain.clear();
        return 0;
    }
};

#endif  // SNOWSIM_SPHERE
