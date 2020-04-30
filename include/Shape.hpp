#ifndef SNOWSIM_SHAPE
#define SNOWSIM_SHAPE

#include "Bounds3.hpp"
#include "Intersection.hpp"
#include "Ray.hpp"

class Intersection;

class Shape
{
   public:
    Shape()
    {
    }
    virtual ~Shape()
    {
    }
    virtual float getArea() = 0;
    virtual float getVolume() = 0;
    virtual Bounds3 getBounds() = 0;
    virtual Intersection getIntersection(Ray ray) = 0;
    virtual bool isInside(const Vector3f& pos) = 0;
    virtual int generateParticlesInside(int lNumberDensity,
                                        std::vector<Vector3f>& contain) = 0;
};

#endif  // SNOWSIM_SHAPE
