#ifndef SNOWSIM_SHAPE
#define SNOWSIM_SHAPE

#include "Bounds3.hpp"

class Shape
{
   public:
    Shape()
    {
    }
    virtual ~Shape()
    {
    }
    virtual Bounds3 getBounds() = 0;
    virtual bool isInside(const Vector3f& pos) = 0;
    virtual int generateParticlesInside(int lNumberDensity,
                                        std::vector<Vector3f>& contain) = 0;
};

#endif  // SNOWSIM_SHAPE
