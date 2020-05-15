#include "simDomain.hpp"
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

SimDomain::SimDomain(SnowParticleSet* SPS, GridMesh* gridMesh)
{
    currentTime = 0;
    this->SPS = SPS;
    this->gridMesh = gridMesh;
}

SimDomain::~SimDomain()
{
}

void SimDomain::initializeSimulator()
{
    gridMesh->initializeGridMeshActiveMassAndMomentum();
    gridMesh->calculateParticleVolume();
}

void SimDomain::oneTimeSimulate()
{
    gridMesh->initializeGridMeshActiveMassAndMomentum();
    // correct delta T after initialize grid mass
    // this is because the first one has to be serial
    // so why not just get the maxV here
    deltaT = correctedDeltaT();
    // gravity
    Vector3f g(0, GRAVITY, 0);
    gridMesh->updateVelocityInGrids(g);
    // Map back to particles
    gridMesh->mapVelocityToSPS();
    currentTime += deltaT;
}

float SimDomain::correctedDeltaT()
{
    float prevMaxVelocity = SPS->maxVelocity;
    float dt;
    if (prevMaxVelocity > 1.e-8)
    {
        float minCellSize =
            std::min(std::min(gridMesh->cellSize[0], gridMesh->cellSize[1]),
                     gridMesh->cellSize[2]);
        float dt = CFL * minCellSize / prevMaxVelocity;
        dt = dt > 1. / FRAMERATE ? 1. / FRAMERATE : dt;
    }
    else
    {
        dt = 1. / FRAMERATE;
    }
    return dt > MAX_TIMESTEP ? MAX_TIMESTEP : dt;
}