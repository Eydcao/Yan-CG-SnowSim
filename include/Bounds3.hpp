//
// Created by LEI XU on 5/16/19.
//

// Modified by Yadi on 4/28/20 for snow simulation

#ifndef SNOWSIM_BOUNDS3
#define SNOWSIM_BOUNDS3
#include <array>
#include <eigen3/Eigen/Eigen>
#include <limits>

using namespace Eigen;

class Bounds3
{
   public:
    Vector3f pMin, pMax;  // two points to specify the bounding box
    Bounds3()
    {
        double minNum = std::numeric_limits<double>::lowest();
        double maxNum = std::numeric_limits<double>::max();
        pMax = Vector3f(minNum, minNum, minNum);
        pMin = Vector3f(maxNum, maxNum, maxNum);
    }
    Bounds3(const Vector3f p) : pMin(p), pMax(p)
    {
    }
    Bounds3(const Vector3f p1, const Vector3f p2)
    {
        pMin = Vector3f(fmin(p1.x(), p2.x()), fmin(p1.y(), p2.y()),
                        fmin(p1.z(), p2.z()));
        pMax = Vector3f(fmax(p1.x(), p2.x()), fmax(p1.y(), p2.y()),
                        fmax(p1.z(), p2.z()));
    }

    Vector3f Diagonal() const
    {
        return pMax - pMin;
    }
    int maxExtent() const
    {
        Vector3f d = Diagonal();
        if (d.x() > d.y() && d.x() > d.z())
            return 0;
        else if (d.y() > d.z())
            return 1;
        else
            return 2;
    }

    double SurfaceArea() const
    {
        Vector3f d = Diagonal();
        return 2 * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
    }

    Vector3f Centroid()
    {
        return 0.5 * pMin + 0.5 * pMax;
    }
    Bounds3 Intersect(const Bounds3& b)
    {
        return Bounds3(
            Vector3f(fmax(pMin.x(), b.pMin.x()), fmax(pMin.y(), b.pMin.y()),
                     fmax(pMin.z(), b.pMin.z())),
            Vector3f(fmin(pMax.x(), b.pMax.x()), fmin(pMax.y(), b.pMax.y()),
                     fmin(pMax.z(), b.pMax.z())));
    }

    Vector3f Offset(const Vector3f& p) const
    {
        Vector3f o = p - pMin;
        if (pMax.x() > pMin.x()) o.x() /= pMax.x() - pMin.x();
        if (pMax.y() > pMin.y()) o.y() /= pMax.y() - pMin.y();
        if (pMax.z() > pMin.z()) o.z() /= pMax.z() - pMin.z();
        return o;
    }

    bool Overlaps(const Bounds3& b1, const Bounds3& b2)
    {
        bool x = (b1.pMax.x() >= b2.pMin.x()) && (b1.pMin.x() <= b2.pMax.x());
        bool y = (b1.pMax.y() >= b2.pMin.y()) && (b1.pMin.y() <= b2.pMax.y());
        bool z = (b1.pMax.z() >= b2.pMin.z()) && (b1.pMin.z() <= b2.pMax.z());
        return (x && y && z);
    }

    bool Inside(const Vector3f& p, const Bounds3& b)
    {
        return (p.x() >= b.pMin.x() && p.x() <= b.pMax.x() &&
                p.y() >= b.pMin.y() && p.y() <= b.pMax.y() &&
                p.z() >= b.pMin.z() && p.z() <= b.pMax.z());
    }
    inline const Vector3f& operator[](int i) const
    {
        return (i == 0) ? pMin : pMax;
    }
};

inline Bounds3 Union(const Bounds3& b1, const Bounds3& b2)
{
    Bounds3 ret;
    ret.pMin =
        Vector3f(fmin(b1.pMin.x(), b2.pMin.x()), fmin(b1.pMin.y(), b2.pMin.y()),
                 fmin(b1.pMin.z(), b2.pMin.z()));
    ret.pMax =
        Vector3f(fmax(b1.pMax.x(), b2.pMax.x()), fmax(b1.pMax.y(), b2.pMax.y()),
                 fmax(b1.pMax.z(), b2.pMax.z()));
    return ret;
}

inline Bounds3 Union(const Bounds3& b, const Vector3f& p)
{
    Bounds3 ret;
    ret.pMin = Vector3f(fmin(b.pMin.x(), p.x()), fmin(b.pMin.y(), p.y()),
                        fmin(b.pMin.z(), p.z()));
    ret.pMax = Vector3f(fmax(b.pMax.x(), p.x()), fmax(b.pMax.y(), p.y()),
                        fmax(b.pMax.z(), p.z()));
    return ret;
}

#endif  // SNOWSIM_BOUNDS3
