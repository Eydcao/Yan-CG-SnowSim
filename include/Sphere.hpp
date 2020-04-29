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
    bool isInside(const Vector3f &pos)
    {
        return (pos - center).squaredNorm() <= radius2;
    }
    int generateParticlesInside(int lNumberDensity,
                                std::vector<Vector3f> &contain)
    {
        // generate particles pos inside
        contain.clear();
        // std::cout << " size is " << contain.size() << std::endl;
        // fill the bbox first, can be represented by maxi,maxj,maxk,dijk
        int numi, numj, numk;
        float di, dj, dk;
        Bounds3 bbox = getBounds();
        Vector3f diag = bbox.Diagonal();
        numi = fmax(lNumberDensity * diag.x(), 1);
        numj = fmax(lNumberDensity * diag.y(), 1);
        numk = fmax(lNumberDensity * diag.z(), 1);
        di = (numi > 2) ? diag.x() / (numi - 1.) : diag.x();
        dj = (numj > 2) ? diag.y() / (numj - 1.) : diag.y();
        dk = (numk > 2) ? diag.z() / (numk - 1.) : diag.z();
        // std::cout << "num i is" << numi << std::endl;
        // std::cout << "num j is" << numj << std::endl;
        // std::cout << "num k is" << numk << std::endl;
        for (int i = 0; i < numi; i++)
        {
            for (int j = 0; j < numj; j++)
            {
                for (int k = 0; k < numk; k++)
                {
                    Vector3f tempPos =
                        bbox.pMin + Vector3f(di * i + 0.5 * di * (numi == 1),
                                             dj * j + 0.5 * dj * (numj == 1),
                                             dk * k + 0.5 * dk * (numk == 1));
                    if (isInside(tempPos)) contain.push_back(tempPos);
                    // std::cout << " size is " << contain.size() << std::endl;
                }
            }
        }
        // filter only keep those inside ones
        // std::cout << " size is " << contain.size() << std::endl;
        return contain.size();
    }
};

#endif  // SNOWSIM_SPHERE
