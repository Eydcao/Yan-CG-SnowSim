// Created by Yadi on 4/28/20 for snow simulation

#ifndef SNOWSIM_RECTANGULAR
#define SNOWSIM_RECTANGULAR

#include "Bounds3.hpp"
#include "Shape.hpp"

class Rectangular : public Shape
{
   public:
    Vector3f pMin, pMax;
    Rectangular(const Vector3f &p1, const Vector3f &p2)
    {
        pMin = Vector3f(fmin(p1.x(), p2.x()), fmin(p1.y(), p2.y()),
                        fmin(p1.z(), p2.z()));
        pMax = Vector3f(fmax(p1.x(), p2.x()), fmax(p1.y(), p2.y()),
                        fmax(p1.z(), p2.z()));
    }
    Bounds3 getBounds()
    {
        return Bounds3(pMin, pMax);
    }
    bool isInside(const Vector3f &pos)
    {
        return (pos.x() >= pMin.x() && pos.x() <= pMax.x() &&
                pos.y() >= pMin.y() && pos.y() <= pMax.y() &&
                pos.z() >= pMin.z() && pos.z() <= pMax.z());
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
        // Bounds3 bbox = getBounds();
        Vector3f diag = pMax - pMin;
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
                        pMin + Vector3f(di * i + 0.5 * di * (numi == 1),
                                        dj * j + 0.5 * dj * (numj == 1),
                                        dk * k + 0.5 * dk * (numk == 1));
                    contain.push_back(tempPos);
                    // std::cout << " size is " << contain.size() << std::endl;
                }
            }
        }
        // filter only keep those inside ones
        // std::cout << " size is " << contain.size() << std::endl;
        return contain.size();
    }
};

#endif  // SNOWSIM_RECTANGULAR
