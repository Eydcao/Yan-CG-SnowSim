#ifndef SNOWSIM_SIMDOMAIN
#define SNOWSIM_SIMDOMAIN

#include "Grid.hpp"
#include "SnowParticle.hpp"
#include "global.hpp"
using namespace Eigen;

class SimDomain
{
   private:
   public:
    float currentTime;
    SnowParticleSet* SPS;
    GridMesh* gridMesh;
    void initializeSimulator();
    // void restartSimulator();
    void oneTimeSimulate();
    // void simulate();
    float correctedDeltaT();
    // void saveBuffer(int time);
    SimDomain(SnowParticleSet* SPS, GridMesh* gridMesh);
    ~SimDomain();
};

#endif  // SNOWSIM_SIMDOMAIN
