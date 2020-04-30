// Modified by Yadi on 4/28/20 for snow simulation
// Original creator not pointed out by source: Yan CG assignment7's triangle.hpp

#ifndef SNOWSIM_TRIANGLE
#define SNOWSIM_TRIANGLE

#include "BVH.hpp"
#include "Intersection.hpp"
#include "OBJ_Loader.hpp"
#include "Shape.hpp"
#include <array>
#include <cassert>

using namespace Eigen;

class Triangle : public Shape
{
   public:
    Vector3f v0, v1, v2;  // vertices A, B ,C , counter-clockwise order
    Vector3f e1, e2;      // 2 edges v1-v0, v2-v0;
    Vector3f normal;
    float area;

    Triangle(Vector3f _v0, Vector3f _v1, Vector3f _v2)
        : v0(_v0), v1(_v1), v2(_v2)
    {
        e1 = v1 - v0;
        e2 = v2 - v0;
        normal = e1.cross(e2).normalized();
        area = e1.cross(e2).norm() * 0.5f;
    }
    Bounds3 getBounds()
    {
        return Union(Bounds3(v0, v1), v2);
    }
    bool isInside(const Vector3f& pos)
    {
        return false;
    }
    int generateParticlesInside(int lNumberDensity,
                                std::vector<Vector3f>& contain)
    {
        contain.clear();
        return 0;
    }
    float getArea()
    {
        return area;
    }
    float getVolume()
    {
        return 0.;
    }
    Intersection getIntersection(Ray ray)
    {
        Intersection inter;

        if (ray.direction.dot(normal) > 0) return inter;
        double u, v, t_tmp = 0;
        Vector3f pvec = ray.direction.cross(e2);
        double det = e1.dot(pvec);
        if (fabs(det) < EPSILON) return inter;

        double det_inv = 1. / det;
        Vector3f tvec = ray.origin - v0;
        u = tvec.dot(pvec) * det_inv;
        if (u < 0. || u > 1.) return inter;
        Vector3f qvec = tvec.cross(e1);
        v = ray.direction.dot(qvec) * det_inv;
        if (v < 0. || u + v > 1.) return inter;
        t_tmp = e2.dot(qvec) * det_inv;

        // find ray triangle intersection
        if (t_tmp > ray.t_min - EPSILON && t_tmp < ray.t_max + EPSILON)
        {
            inter.coords =
                u * this->v1 + v * this->v2 + (1.f - u - v) * this->v0;
            inter.normal = this->normal;
            inter.distance = t_tmp;
            inter.shape = this;
            // inter.m = this->m;
            inter.happened = true;
        }
        return inter;
    }

    Vector3f centroid()
    {
        return (v0 + v1 + v2) / 3.;
    }
};

class MeshTriangle : public Shape
{
   public:
    MeshTriangle(const std::string& filename)
    {
        objl::Loader loader;
        loader.LoadFile(filename);
        area = 0;
        assert(loader.LoadedMeshes.size() == 1);
        auto mesh = loader.LoadedMeshes[0];

        Vector3f min_vert = Vector3f{std::numeric_limits<float>::infinity(),
                                     std::numeric_limits<float>::infinity(),
                                     std::numeric_limits<float>::infinity()};
        Vector3f max_vert = Vector3f{-std::numeric_limits<float>::infinity(),
                                     -std::numeric_limits<float>::infinity(),
                                     -std::numeric_limits<float>::infinity()};
        for (int i = 0; i < mesh.Vertices.size(); i += 3)
        {
            std::array<Vector3f, 3> face_vertices;

            for (int j = 0; j < 3; j++)
            {
                auto vert = Vector3f(mesh.Vertices[i + j].Position.X,
                                     mesh.Vertices[i + j].Position.Y,
                                     mesh.Vertices[i + j].Position.Z);
                face_vertices[j] = vert;

                min_vert = Vector3f(std::min(min_vert.x(), vert.x()),
                                    std::min(min_vert.y(), vert.y()),
                                    std::min(min_vert.z(), vert.z()));
                max_vert = Vector3f(std::max(max_vert.x(), vert.x()),
                                    std::max(max_vert.y(), vert.y()),
                                    std::max(max_vert.z(), vert.z()));
            }

            triangles.emplace_back(face_vertices[0], face_vertices[1],
                                   face_vertices[2]);
        }

        bounding_box = Bounds3(min_vert, max_vert);

        std::vector<Shape*> ptrs;

        std::cout << "how many triangles " << triangles.size() << std::endl;

        volume = 0.;

        for (auto& tri : triangles)
        {
            ptrs.push_back(&tri);
            area += tri.area;
            float volumeGuassian =
                tri.centroid().dot(tri.normal) * tri.area / 3.;
            volume += volumeGuassian;
        }
        bvh = new BVHAccel(ptrs);
    }

    // bool intersect(const Ray& ray, float& tnear, uint32_t& index) const
    // {
    //     bool intersect = false;
    //     for (uint32_t k = 0; k < numTriangles; ++k)
    //     {
    //         const Vector3f& v0 = vertices[vertexIndex[k * 3]];
    //         const Vector3f& v1 = vertices[vertexIndex[k * 3 + 1]];
    //         const Vector3f& v2 = vertices[vertexIndex[k * 3 + 2]];
    //         float t, u, v;
    //         if (rayTriangleIntersect(v0, v1, v2, ray.origin, ray.direction,
    //         t,
    //                                  u, v) &&
    //             t < tnear)
    //         {
    //             tnear = t;
    //             index = k;
    //             intersect |= true;
    //         }
    //     }

    //     return intersect;
    // }

    Bounds3 getBounds()
    {
        return bounding_box;
    }
    bool isInside(const Vector3f& pos)
    {
        Vector3f outerPos = bounding_box.pMax + bounding_box.Diagonal();
        Ray ray(pos, (outerPos - pos).normalized(), 0.0);
        // Vector3f dir = Vector3f(1., 1., 1.).normalized();
        // Ray ray(pos, dir, 0.0);
        std::vector<Intersection> tempInters = {};
        this->bvh->allIntersect(ray, tempInters);
        // connect a ray between pos and outer pos
        // use bvh to traverse all possible intersections (ts)
        std::vector<float> ts;
        bool isNearParticle = false;
        for (auto& oneInter : tempInters)
        {
            if (std::abs(oneInter.distance) < 0.00866) isNearParticle = true;
            ts.push_back(oneInter.distance);
        }
        if (isNearParticle) return true;
        int size = ts.size();
        // std::cout << "intersection size is " << size << std::endl;
        // TODO consider very close ts, if they are the same intersect or
        // tangential to the surface
        // int size = tempInters.size();
        // if (size == 2)
        // {
        //     std::cout << "intersection size is " << size << std::endl;
        //     std::cout << "t1 is " << ts[0] << std::endl;
        //     if (size == 2) std::cout << "t2 is " << ts[1] << std::endl;
        // }
        // std::cout << " intersection for this particle is " << size <<
        // std::endl;

        if (size % 2 == 0)
            return false;
        else
            return true;
    }
    int generateParticlesInside(int lNumberDensity,
                                std::vector<Vector3f>& contain)
    {
        // generate particles pos inside
        contain.clear();
        // std::cout << " size is " << contain.size() << std::endl;
        // fill the bbox first, can be represented by maxi,maxj,maxk,dijk
        int numi, numj, numk;
        float di, dj, dk;
        Bounds3 bbox = getBounds();
        Vector3f diag = bbox.Diagonal();
        // std::cout << "bbox diag is" << diag << std::endl;
        numi = fmax(lNumberDensity * diag.x(), 1);
        numj = fmax(lNumberDensity * diag.y(), 1);
        numk = fmax(lNumberDensity * diag.z(), 1);
        // std::cout << "num i is" << numi << std::endl;
        // std::cout << "num j is" << numj << std::endl;
        // std::cout << "num k is" << numk << std::endl;
        di = (numi > 2) ? diag.x() / (numi - 1.) : diag.x();
        dj = (numj > 2) ? diag.y() / (numj - 1.) : diag.y();
        dk = (numk > 2) ? diag.z() / (numk - 1.) : diag.z();
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
        // std::cout << " size is " << contain.size() << std::endl;
        return contain.size();
    }

    Intersection getIntersection(Ray ray)
    {
        Intersection intersec;

        if (bvh)
        {
            intersec = bvh->Intersect(ray);
        }

        return intersec;
    }

    float getArea()
    {
        return area;
    }

    float getVolume()
    {
        return volume;
    }

    Bounds3 bounding_box;
    uint32_t numTriangles;

    std::vector<Triangle> triangles;

    BVHAccel* bvh;
    float area;
    float volume;
};

// inline Intersection Triangle::getIntersection(Ray ray)
// {
//     Intersection inter;

//     if (dotProduct(ray.direction, normal) > 0) return inter;
//     double u, v, t_tmp = 0;
//     Vector3f pvec = crossProduct(ray.direction, e2);
//     double det = dotProduct(e1, pvec);
//     if (fabs(det) < EPSILON) return inter;

//     double det_inv = 1. / det;
//     Vector3f tvec = ray.origin - v0;
//     u = dotProduct(tvec, pvec) * det_inv;
//     if (u < 0 || u > 1) return inter;
//     Vector3f qvec = crossProduct(tvec, e1);
//     v = dotProduct(ray.direction, qvec) * det_inv;
//     if (v < 0 || u + v > 1) return inter;
//     t_tmp = dotProduct(e2, qvec) * det_inv;

//     // find ray triangle intersection
//     if (t_tmp > ray.t_min && t_tmp < ray.t_max)
//     {
//         inter.coords = u * this->v1 + v * this->v2 + (1.f - u - v) *
//         this->v0; inter.normal = this->normal; inter.distance = t_tmp;
//         inter.obj = this;
//         inter.m = this->m;
//         inter.happened = true;
//     }
//     return inter;
// }

#endif  // SNOWSIM_TRIANGLE
