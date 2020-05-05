#ifndef SNOWSIM_GRID
#define SNOWSIM_GRID

#include "SnowParticle.hpp"
#include "global.hpp"
#include <cstring>
#include <math.h>
#include <stdio.h>
using namespace Eigen;

const float BSPLINE_EPSILON = 1e-4;
const int BSPLINE_RADIUS = 2;

class SnowParticleMaterial;
class SnowParticle;
class SnowParticleSet;

struct GridCell
{
    // if cell has snow particle, then active
    bool active;
    float mass;
    Vector3f velocity;
    // star is super script, meaning temp or intermedia velocity, under explicit
    // mode velocity star will be new velocity
    Vector3f velocityStar;

    // the following pars are only used in implicit mode
    bool implicitSolve;
    Vector3f force;
    Vector3f error;
    Vector3f residual;
    // TODO fig out what is p, gradient of resiudal or error?
    Vector3f p;
    // TODO fig out what EP ER is
    Vector3f Ep;
    Vector3f Er;
    // dotProduct(r, Er)
    float rEr;
};

class GridMesh
{
   private:
   public:
    // Vector3f origin;
    // Vector3f size;
    Bounds3 bbox;
    // Vector3f cellsize;
    // how many cells in x, y, and z direction
    Vector3i cellNum;
    Vector3f cellSize;
    SnowParticleSet* SPS;
    float eachCellVolume;
    // Nodes: use (y*size[0] + x) to index, where zero is the bottom-left corner
    // (e.g. like a cartesian grid)
    int totalCellNum;
    std::vector<GridCell*> cells;

    // Grid be at least one cell; there must be one layer of cells surrounding
    // all particles
    GridMesh(const Bounds3& bbox, const Vector3f& cellSize,
             SnowParticleSet* SPS);
    // GridMesh(const GridMesh& anotherGridMesh);
    ~GridMesh();

    // Map particles to grid
    void initializeMass();
    void initializeVelocities();
    // Map grid volumes back to particles (first timestep only)
    void calculateVolumes() const;
    // Compute grid velocities
    void explicitVelocities(const Vector3f& gravity);
    // #if ENABLE_IMPLICIT
    void implicitVelocities();
    void recomputeImplicitForces();
    // #endif
    // Map grid velocities back to particles
    void updateVelocities() const;

    // Collision detection
    void collisionGrid();
    void collisionParticles() const;

    // TODO change the corresponding spline curve in 3D
    // Cubic B-spline shape/basis/interpolation function
    // A smooth curve from (0,1) to (1,0)
    static float bspline(float x)
    {
        x = fabs(x);
        float w;
        if (x < 1)
            w = x * x * (x / 2 - 1) + 2 / 3.0;
        else if (x < 2)
            w = x * (x * (-x / 6 + 1) - 2) + 4 / 3.0;
        else
            return 0;
        // Clamp between 0 and 1... if needed
        if (w < BSPLINE_EPSILON) return 0;
        return w;
    }
    // Slope of interpolation function
    static float bsplineSlope(float x)
    {
        float abs_x = fabs(x), w;
        if (abs_x < 1)
            return 1.5 * x * abs_x - 2 * x;
        else if (x < 2)
            return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
        else
            return 0;
        // Clamp between -2/3 and 2/3... if needed
    }
};

#endif  // SNOWSIM_GRID
