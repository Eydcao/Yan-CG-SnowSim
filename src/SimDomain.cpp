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
    gridMesh->initializeMass();
    gridMesh->calculateVolumes();
}

void SimDomain::oneTimeSimulate()
{
    deltaT = correctedDeltaT();
    gridMesh->initializeMass();
    gridMesh->initializeVelocities();
    // gravity
    Vector3f g(0, GRAVITY, 0);
    gridMesh->explicitVelocities(g);
    // TODO after finish explicit then do implicit
    // #if ENABLE_IMPLICIT
    //     if (IMPLICIT_RATIO > 0) grid->implicitVelocities();
    // #endif
    // Map back to particles
    gridMesh->updateVelocities();
    // Update particle data
    SPS->update();
    currentTime += deltaT;
}

float SimDomain::correctedDeltaT()
{
    float prevMaxVelocity = SPS->maxVelocity;
    float dt;
    if (prevMaxVelocity > 1.e-8)
    {
        // We should really take the min(cellsize) I think, if the grid is not
        // square
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