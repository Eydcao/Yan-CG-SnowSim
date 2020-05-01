#include "BVH.hpp"
#include <algorithm>
#include <cassert>

BVHAccel::BVHAccel(std::vector<Shape*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty()) return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    // printf(
    //     "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i
    //     secs\n\n", hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Shape*> shapes)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < shapes.size(); ++i)
        bounds = Union(bounds, shapes[i]->getBounds());
    if (shapes.size() == 1)
    {
        // Create leaf _BVHBuildNode_
        node->bounds = shapes[0]->getBounds();
        node->shape = shapes[0];
        node->left = nullptr;
        node->right = nullptr;
        node->area = shapes[0]->getArea();
        return node;
    }
    else if (shapes.size() == 2)
    {
        node->left = recursiveBuild(std::vector{shapes[0]});
        node->right = recursiveBuild(std::vector{shapes[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        node->area = node->left->area + node->right->area;
        return node;
    }
    else
    {
        Bounds3 centroidBounds;
        for (int i = 0; i < shapes.size(); ++i)
            centroidBounds =
                Union(centroidBounds, shapes[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim)
        {
            case 0:
                std::sort(shapes.begin(), shapes.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().x() <
                           f2->getBounds().Centroid().x();
                });
                break;
            case 1:
                std::sort(shapes.begin(), shapes.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().y() <
                           f2->getBounds().Centroid().y();
                });
                break;
            case 2:
                std::sort(shapes.begin(), shapes.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().z() <
                           f2->getBounds().Centroid().z();
                });
                break;
        }

        auto beginning = shapes.begin();
        auto middling = shapes.begin() + (shapes.size() / 2);
        auto ending = shapes.end();

        auto leftshapes = std::vector<Shape*>(beginning, middling);
        auto rightshapes = std::vector<Shape*>(middling, ending);

        assert(shapes.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
        node->area = node->left->area + node->right->area;
    }

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root) return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

int BVHAccel::allIntersect(const Ray& ray,
                           std::vector<Intersection>& contain) const
{
    contain.clear();
    if (!root) return 0;
    getAllIntersection(root, ray, contain);
    return contain.size();
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // Traverse the BVH to find intersection
    // child case
    if (!node->left && !node->right)
    {
        assert(node->shape);
        Intersection inter = node->shape->getIntersection(ray);
        return inter;
    }
    // ray inverse dir and isDirNeg
    Vector3f dir = ray.direction;
    Vector3f invDir = ray.direction_inv;
    std::array<int, 3> dirIsNeg = {int(dir.x() > 0), int(dir.y() > 0),
                                   int(dir.z() > 0)};
    // bb judge for left
    Intersection interL;
    if (node->left->bounds.IntersectP(ray, invDir, dirIsNeg))
    {
        interL = getIntersection(node->left, ray);
    }
    // then for right
    Intersection interR;
    if (node->right->bounds.IntersectP(ray, invDir, dirIsNeg))
    {
        interR = getIntersection(node->right, ray);
    }
    Intersection inter = (interL.distance < interR.distance) ? interL : interR;
    return inter;
}

int BVHAccel::getAllIntersection(BVHBuildNode* node, const Ray& ray,
                                 std::vector<Intersection>& contain) const
{
    // Traverse the BVH to find ALL intersection
    // child case
    if (!node->left && !node->right)
    {
        assert(node->shape);
        Intersection inter = node->shape->getIntersection(ray);
        if (inter.happened) contain.push_back(inter);
        return contain.size();
    }
    // ray inverse dir and isDirNeg
    Vector3f dir = ray.direction;
    Vector3f invDir = ray.direction_inv;
    std::array<int, 3> dirIsNeg = {int(dir.x() > 0), int(dir.y() > 0),
                                   int(dir.z() > 0)};
    // bb judge for left
    if (node->left->bounds.IntersectP(ray, invDir, dirIsNeg))
    {
        getAllIntersection(node->left, ray, contain);
    }
    // then for right
    if (node->right->bounds.IntersectP(ray, invDir, dirIsNeg))
    {
        getAllIntersection(node->right, ray, contain);
    }
    return contain.size();
}

// void BVHAccel::getSample(BVHBuildNode* node, float p, Intersection& pos,
//                          float& pdf)
// {
//     if (node->left == nullptr || node->right == nullptr)
//     {
//         node->shape->Sample(pos, pdf);
//         pdf *= node->area;
//         return;
//     }
//     if (p < node->left->area)
//         getSample(node->left, p, pos, pdf);
//     else
//         getSample(node->right, p - node->left->area, pos, pdf);
// }

// void BVHAccel::Sample(Intersection& pos, float& pdf)
// {
//     float p = std::sqrt(get_random_float()) * root->area;
//     getSample(root, p, pos, pdf);
//     pdf /= root->area;
// }