//
// Created by LEI XU on 5/16/19.
//

// Modified by Yadi on 4/28/20 for snow simulation

#ifndef SNOWSIM_BVH
#define SNOWSIM_BVH

#include "Bounds3.hpp"
#include "Intersection.hpp"
#include "Ray.hpp"
#include "Shape.hpp"
#include <atomic>
#include <ctime>
#include <memory>
#include <vector>

struct BVHBuildNode;
// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;

// BVHAccel Declarations
inline int leafNodes, totalLeafNodes, totalPrimitives, interiorNodes;
class BVHAccel
{
   public:
    // BVHAccel Public Types
    enum class SplitMethod
    {
        NAIVE,
        SAH
    };

    // BVHAccel Public Methods
    BVHAccel(std::vector<Shape *> p, int maxPrimsInNode = 1,
             SplitMethod splitMethod = SplitMethod::NAIVE);
    Bounds3 WorldBound() const;
    ~BVHAccel();

    Intersection Intersect(const Ray &ray) const;
    int allIntersect(const Ray &ray, std::vector<Intersection> &contain) const;
    Intersection getIntersection(BVHBuildNode *node, const Ray &ray) const;
    int getAllIntersection(BVHBuildNode *node, const Ray &ray,
                           std::vector<Intersection> &contain) const;
    // bool IntersectP(const Ray &ray) const;
    BVHBuildNode *root;

    // BVHAccel Private Methods
    BVHBuildNode *recursiveBuild(std::vector<Shape *> objects);

    // BVHAccel Private Data
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    std::vector<Shape *> primitives;

    // void getSample(BVHBuildNode *node, float p, Intersection &pos, float
    // &pdf); void Sample(Intersection &pos, float &pdf);
};

struct BVHBuildNode
{
    Bounds3 bounds;
    BVHBuildNode *left;
    BVHBuildNode *right;
    Shape *shape;
    float area;

   public:
    int splitAxis = 0, firstPrimOffset = 0, nPrimitives = 0;
    // BVHBuildNode Public Methods
    BVHBuildNode()
    {
        bounds = Bounds3();
        left = nullptr;
        right = nullptr;
        shape = nullptr;
    }
};

#endif  // SNOWSIM_BVH
