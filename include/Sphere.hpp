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
    float getArea()
    {
        return 4. * M_PI * radius2;
    }
    float getVolume()
    {
        return 4. / 3. * M_PI * radius2 * radius;
    }
    Intersection getIntersection(Ray ray)
    {
        Intersection result;
        result.happened = false;
        Vector3f L = ray.origin - center;
        float a = ray.direction.squaredNorm();
        float b = 2. * ray.direction.dot(L);
        float c = L.squaredNorm() - radius2;
        float t0, t1;
        if (!solveQuadratic(a, b, c, t0, t1)) return result;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return result;
        result.happened = true;

        result.coords = ray.origin + ray.direction * t0;
        result.normal = (result.coords - center).normalized();
        // result.m = this->m;
        result.shape = this;
        result.distance = t0;
        return result;
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
