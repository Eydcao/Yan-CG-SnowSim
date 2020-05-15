#include "Grid.hpp"

GridCell::GridCell()
{
    active = false;
    // the mass will be mapped from adj pars each time of simulation
    mass = 0;
    velocity.setZero();
    velocityStar.setZero();
    force.setZero();
    affectedParticles = {};
    affectedParticlesWeightID = {};

    error.setZero();
    residual.setZero();
    p.setZero();
    Ep.setZero();
    Er.setZero();
    rEr = 0;
}

GridCell::~GridCell()
{
}

void GridCell::resetCell()
{
    active = false;
    // the mass will be mapped from adj pars each time of simulation
    mass = 0;
    velocity.setZero();
    velocityStar.setZero();
    force.setZero();
    affectedParticles = {};
    affectedParticlesWeightID = {};

    error.setZero();
    residual.setZero();
    p.setZero();
    Ep.setZero();
    Er.setZero();
    rEr = 0;
}

GridMesh::GridMesh(const Bounds3& bbox, const Vector3f& cellSize,
                   SnowParticleSet* SPS)
{
    this->SPS = SPS;
    this->cellSize = cellSize;
    eachCellVolume = cellSize.x() * cellSize.y() * cellSize.z();
    // the real bbox may be bigger than the input one
    int numx, numy, numz;
    Vector3f oldDiag = bbox.Diagonal();
    Vector3f tempCellNum = oldDiag.cwiseProduct(cellSize.cwiseInverse());
    numx = std::max((int)tempCellNum.x() + 1, 3);
    numy = std::max((int)tempCellNum.y() + 1, 3);
    numz = std::max((int)tempCellNum.z() + 1, 3);
    cellNum = Vector3i(numx, numy, numz);
    totalCellNum = numx * numy * numz;
    // calibrate bbox
    Vector3f newDiag =
        Vector3f((float)numx * cellSize.x(), (float)numy * cellSize.y(),
                 (float)numz * cellSize.z());
    Vector3f halfDiffDiag = .5 * (newDiag - oldDiag);
    Vector3f newPmin = bbox.pMin - halfDiffDiag;
    Vector3f newPmax = bbox.pMax + halfDiffDiag;
    this->bbox = Bounds3(newPmin, newPmax);
    // new cells
    cells.resize(totalCellNum);
    for (int i = 0; i < numx; i++)
    {
        for (int j = 0; j < numy; j++)
        {
            for (int k = 0; k < numz; k++)
            {
                // why write 3 loops instead of one
                // just in case sth related to xyz dimension will be calculated
                GridCell* tempCell = new GridCell();
                tempCell->i = i;
                tempCell->j = j;
                tempCell->k = k;
                cells[i * numy * numz + j * numz + k] = tempCell;
            }
        }
    }
    // now we dont know how many cells effective yet
    totalEffectiveCellNum = 0;
    effectiveCells = {};
}

GridMesh::~GridMesh()
{
    // TODO check if this destroy way is okay
    for (auto& oneCell : cells)
    {
        delete (oneCell);
    }
}

// Maps mass and velocity to the grid
void GridMesh::initializeGridMeshActiveMassAndMomentum()
{
    // this loop will reset the old EFFECTIVE cell; TODO trivial parfor
    for (int iCell = 0; iCell < totalEffectiveCellNum; iCell++)
    {
        GridCell* cell = effectiveCells[iCell];
        cell->resetCell();
    }
    effectiveCells.clear();

    // reset SPSbbox
    SPSbbox = Bounds3(SPS->particles[0]->position);

    // This initial mapping has to serial loop with particles
    // after this is done create a 'effective cells'
    // this is for the sake of parfor efficiency later (avoid uneffective cells)

    float& maxVelocity = SPS->maxVelocity;
    maxVelocity = 0;
    for (int i = 0; i < SPS->particles.size(); i++)
    {
        SnowParticle* p = SPS->particles[i];
        // Solve for grid internal forces
        // Matrix3f energyD = p->energyDerivative();
        Vector3f pos = p->position;
        // greedy to get the max vel
        Vector3f vel = p->velocity;
        float v = vel.norm();
        if (v > maxVelocity) maxVelocity = v;
        SPSbbox = Union(SPSbbox, pos);
        Vector3f posRel =
            (pos - bbox.pMin).cwiseProduct(cellSize.cwiseInverse());
        for (int countForAdjGrid = 0, ipx = (int)posRel.x() - 1;
             ipx < (int)posRel.x() + 3; ipx++)
        {
            for (int ipy = (int)posRel.y() - 1; ipy < (int)posRel.y() + 3;
                 ipy++)
            {
                for (int ipz = (int)posRel.z() - 1; ipz < (int)posRel.z() + 3;
                     ipz++)
                {
                    // see if it is inside the bbox
                    // no: this weight and weight grad is 0, nothing else
                    // yes: can calculate the weights and weight grads
                    //      then push back this p and the count to corresponding
                    //      cell for record; mark that cell to active; addup
                    //      cell's mass
                    if (ipx >= 0 && ipy >= 0 && ipz >= 0 && ipx < cellNum.x() &&
                        ipy < cellNum.y() && ipz < cellNum.z())
                    {
                        // get the cell
                        int iCell = ipx * cellNum.y() * cellNum.z() +
                                    ipy * cellNum.z() + ipz;
                        GridCell* cell = cells[iCell];
                        // hat means this is non dimensional
                        float xHat = (float)ipx - posRel.x();
                        float yHat = (float)ipy - posRel.y();
                        float zHat = (float)ipz - posRel.z();
                        float wx = bspline(xHat);
                        float wy = bspline(yHat);
                        float wz = bspline(zHat);
                        float w = wx * wy * wz;
                        if (w > BSPLINE_EPSILON)
                        {
                            // get grad by chain's rule (NOTICE that xyz hat is
                            // dimensional)
                            Vector3f wGrad(
                                wy * wz * bsplineSlope(xHat) / cellSize.x(),
                                wz * wx * bsplineSlope(yHat) / cellSize.y(),
                                wx * wy * bsplineSlope(zHat) / cellSize.z());
                            p->weights[countForAdjGrid] = w;
                            p->weight_gradient[countForAdjGrid] = wGrad;
                            // push to cell
                            cell->affectedParticles.push_back(p);
                            cell->affectedParticlesWeightID.push_back(
                                countForAdjGrid);
                            // this cell is active now
                            cell->active = true;
                            // add mass contribution
                            cell->mass += w * p->mass;
                            // add momentum contribution into velocity
                            // later devide by mass then get velocity
                            cell->velocity += w * p->mass * p->velocity;
                            // if (ipx == 25 && ipy == 1 && ipz == 30)
                            // {
                            //     std::cout << " w is " << w << std::endl;
                            //     std::cout << " wx is " << wx << std::endl;
                            //     std::cout << " wy is " << wy << std::endl;
                            //     std::cout << " wz is " << wz << std::endl;
                            //     std::cout << " p mass " << p->mass <<
                            //     std::endl; std::cout << " p vel " <<
                            //     p->velocity
                            //               << std::endl;
                            //     std::cout << " cell id is " << cell->i << ",
                            //     "
                            //               << cell->j << ", " << cell->k
                            //               << std::endl;
                            // }
                        }
                        else
                        {
                            p->weights[countForAdjGrid] = 0;
                            p->weight_gradient[countForAdjGrid].setZero();
                        }
                    }
                    else
                    {
                        p->weights[countForAdjGrid] = 0;
                        p->weight_gradient[countForAdjGrid].setZero();
                    }
                    // do it or not, it has goto next cell
                    countForAdjGrid++;
                }
            }
        }
    }

    // std::cout << " the SPS bbox Pmin is " << SPSbbox.pMin.x() << ", "
    //           << SPSbbox.pMin.y() << ", " << SPSbbox.pMin.z() << std::endl;
    // std::cout << " the SPS bbox Pmax is " << SPSbbox.pMax.x() << ", "
    //           << SPSbbox.pMax.y() << ", " << SPSbbox.pMax.z() << std::endl;

    // then fill the new EFFECTIVE cells
    for (int iCell = 0; iCell < totalCellNum; iCell++)
    {
        GridCell* cell = cells[iCell];
        if (cell->active) effectiveCells.push_back(cell);
    }
    totalEffectiveCellNum = effectiveCells.size();
}

void GridMesh::updateVelocityInGrids(const Vector3f& gravity)
{
#pragma omp parallel for
    for (int iCell = 0; iCell < totalEffectiveCellNum; iCell++)
    {
        // get vel
        {
            GridCell* cell = effectiveCells[iCell];
            cell->velocity /= cell->mass;
            // test
            // if (isnan(cell->velocity.norm()))
            // {
            //     std::cout << " cell velocity old is " << cell->velocity
            //               << std::endl;
            //     std::cout << " cell mass " << cell->mass << std::endl;
            //     std::cout << " cell id is " << cell->i << ", " << cell->j <<
            //     ", "
            //               << cell->k << std::endl;
            //     exit(0);
            // }
        }
        // get vel star
        {
            // First, compute the forces
            GridCell* cell = effectiveCells[iCell];
            Vector3f force = Vector3f::Zero();
            for (int i = 0; i < cell->affectedParticles.size(); i++)
            {
                SnowParticle* p = cell->affectedParticles[i];
                int countForAdjGrid = cell->affectedParticlesWeightID[i];
                float w = p->weights[countForAdjGrid];
                Vector3f wGrad = p->weight_gradient[countForAdjGrid];
                Matrix3f energyD = p->energyDerivative();
                if (w > BSPLINE_EPSILON)
                {
                    force += energyD * wGrad;
                    // std::cout << " force + energyD is " << energyD <<
                    // std::endl; std::cout << " force + wGrad is " << wGrad <<
                    // std::endl;
                }
            }
            // std::cout << " force is " << force << std::endl;
            // std::cout << " g is " << gravity << std::endl;
            // std::cout << " mass is " << cell->mass << std::endl;
            // force+gravity->a->dv->vStar
            cell->velocityStar =
                cell->velocity + deltaT * (gravity - force / cell->mass);

            // if (isnan(cell->velocityStar.norm()))
            // {
            //     std::cout << " cell velocity star is " << cell->velocityStar
            //               << std::endl;
            //     std::cout << " cell velocity old is " << cell->velocity
            //               << std::endl;
            //     std::cout << " cell velocity mass " << cell->mass <<
            //     std::endl; std::cout << " cell id is " << cell->i << ", " <<
            //     cell->j << ", "
            //               << cell->k << std::endl;
            //     exit(0);
            // }
            // if (cell->i == 0 && cell->j == 0 && cell->k == 16)
            // {
            //     std::cout << " cell velocity star is " << cell->velocityStar
            //               << std::endl;
            //     std::cout << " cell velocity old is " << cell->velocity
            //               << std::endl;
            //     std::cout << " cell velocity mass " << cell->mass <<
            //     std::endl; std::cout << " cell id is " << cell->i << ", " <<
            //     cell->j << ", "
            //               << cell->k << std::endl;
            //     // exit(0);
            // }
        }
        // todo consider implicit

        // collision test
        {
            // int i = iCell / (cellNum.y() * cellNum.z());
            // int j = (iCell % (cellNum.y() * cellNum.z())) / cellNum.z();
            // int k = iCell % cellNum.z();

            // int iCell = i * cellNum.y() * cellNum.z() + j * cellNum.z() + k;
            GridCell* cell = cells[iCell];
            int i = cell->i;
            int j = cell->j;
            int k = cell->k;
            // collision test
            Vector3f oldCellPos((i + 0.5) * cellSize.x(),
                                (j + 0.5) * cellSize.y(),
                                (k + 0.5) * cellSize.z());
            Vector3f tempPos = oldCellPos + deltaT * cell->velocityStar;
            // TODO consider wrap this fixing vel as a function
            // x dir
            if (tempPos.x() > bbox.pMax.x() || tempPos.x() < bbox.pMin.x())
            {
                Vector3f VSolid(0, 0, 0);
                Vector3f VRel = cell->velocityStar - VSolid;
                Vector3f Vn = Vector3f(VRel.x(), 0, 0);
                Vector3f Vt = Vector3f(0, VRel.y(), VRel.z());
                float vn = Vn.norm();
                float vt = Vt.norm();
                // TODO this mu calculation should be cleverer
                float mu = SPS->particles[0]->m->sticky;
                if (vt < mu * vn || vt < 1.e-12)
                {
                    VRel.setZero();
                }
                else
                {
                    VRel = (1. - mu * vn / vt) * Vt;
                }
                cell->velocityStar = VRel + VSolid;
                // std::cout << " cell velocity star is " << cell->velocityStar
                //           << std::endl;
            }
            // y dir
            if (tempPos.y() > bbox.pMax.y() || tempPos.y() < bbox.pMin.y())
            {
                Vector3f VSolid(0, 0, 0);
                Vector3f VRel = cell->velocityStar - VSolid;
                Vector3f Vn = Vector3f(0, VRel.y(), 0);
                Vector3f Vt = Vector3f(VRel.x(), 0, VRel.z());
                float vn = Vn.norm();
                float vt = Vt.norm();
                // TODO this mu calculation should be cleverer
                float mu = SPS->particles[0]->m->sticky;
                if (vt < mu * vn || vt < 1.e-12)
                {
                    VRel.setZero();
                }
                else
                {
                    VRel = (1. - mu * vn / vt) * Vt;
                }
                cell->velocityStar = VRel + VSolid;
                // std::cout << " cell velocity star is " << cell->velocityStar
                //           << std::endl;
            }
            // z dir
            if (tempPos.z() > bbox.pMax.z() || tempPos.z() < bbox.pMin.z())
            {
                Vector3f VSolid(0, 0, 0);
                Vector3f VRel = cell->velocityStar - VSolid;
                Vector3f Vn = Vector3f(0, 0, VRel.z());
                Vector3f Vt = Vector3f(VRel.x(), VRel.y(), 0);
                float vn = Vn.norm();
                float vt = Vt.norm();
                // TODO this mu calculation should be cleverer
                float mu = SPS->particles[0]->m->sticky;
                if (vt < mu * vn || vt < 1.e-12)
                {
                    VRel.setZero();
                }
                else
                {
                    VRel = (1. - mu * vn / vt) * Vt;
                }
                cell->velocityStar = VRel + VSolid;
                // std::cout << " cell velocity star is " << cell->velocityStar
                //           << std::endl;
            }
            // TODO there are more collisions
            // Finish collision test

            // if (isnan(cell->velocityStar.norm()))
            // {
            //     std::cout << " cell velocity star is " << cell->velocityStar
            //               << std::endl;
            //     std::cout << " cell velocity old is " << cell->velocity
            //               << std::endl;
            //     std::cout << " cell velocity mass " << cell->mass <<
            //     std::endl; std::cout << " cell id is " << cell->i << ", " <<
            //     cell->j << ", "
            //               << cell->k << std::endl;
            //     exit(0);
            // }
        }
    }
}

// void GridMesh::initializeGridVelocity()
// {
//     // TODO trivial parfor
//     for (int iCell = 0; iCell < totalEffectiveCellNum; iCell++)
//     {
//         GridCell* cell = effectiveCells[iCell];
//         cell->velocity /= cell->mass;
//         // test
//         // if (isnan(cell->velocity.norm()))
//         // {
//         //     std::cout << " cell velocity old is " << cell->velocity
//         //               << std::endl;
//         //     std::cout << " cell mass " << cell->mass << std::endl;
//         //     std::cout << " cell id is " << cell->i << ", " << cell->j <<
//         ", "
//         //               << cell->k << std::endl;
//         //     exit(0);
//         // }
//     }
// }

// Maps volume from the grid to particles
// This should only be called once, at the beginning of the simulation
void GridMesh::calculateParticleVolume() const
{
    // Volume is calculated once and stored
    // density, however, will update each time for viewing purposes
    // NB YADI : seems that for viewing purpose it only needs density
    //           so vol still only updates once
    for (int i = 0; i < SPS->particles.size(); i++)
    {
        SnowParticle* p = SPS->particles[i];
        // Solve for grid internal forces
        // Matrix3f energyD = p->energyDerivative();
        Vector3f pos = p->position - bbox.pMin;
        p->density = 0;
        // temp sum the mass in p density
        for (int countForAdjGrid = 0, px = (int)(pos.x() / cellSize.x()),
                 ipx = px - 1, pxEnd = px + 3;
             ipx < pxEnd; ipx++)
        {
            if (ipx >= 0 && ipx < cellNum.x())
            {
                for (int py = (int)(pos.y() / cellSize.y()), ipy = py - 1,
                         pyEnd = py + 3;
                     ipy < pyEnd; ipy++)
                {
                    if (ipy >= 0 && ipy < cellNum.y())
                    {
                        for (int pz = (int)(pos.z() / cellSize.z()),
                                 ipz = pz - 1, pzEnd = pz + 3;
                             ipz < pzEnd; ipz++)
                        {
                            if (ipz >= 0 && ipz < cellNum.z())
                            {
                                float w = p->weights[countForAdjGrid];
                                if (w > BSPLINE_EPSILON)
                                {
                                    int iCell =
                                        ipx * cellNum.y() * cellNum.z() +
                                        ipy * cellNum.z() + ipz;
                                    p->density += w * cells[iCell]->mass;
                                    // std::cout << " here " << std::endl;
                                }
                            }
                            // do it or not, it has goto next cell
                            countForAdjGrid++;
                        }
                    }
                    else
                    {
                        // do nothing and skip this z edge
                        countForAdjGrid += 4;
                    }
                }
            }
            else
            {
                // do nothing and skip this yz face
                countForAdjGrid += 16;
            }
        }
        p->density /= eachCellVolume;
        p->volume = p->mass / p->density;
        // std::cout << " rho is" << p->density << std::endl;
        // std::cout << " vol is" << p->volume << std::endl;
    }
}

// // Calculate next timestep velocities using pure explicit method
// void GridMesh::updateGridVelocityStar(const Vector3f& gravity)
// {
//     // TODO trivial parfor
//     for (int iCell = 0; iCell < totalEffectiveCellNum; iCell++)
//     {
//         // First, compute the forces
//         GridCell* cell = effectiveCells[iCell];
//         Vector3f force = Vector3f::Zero();
//         for (int i = 0; i < cell->affectedParticles.size(); i++)
//         {
//             SnowParticle* p = cell->affectedParticles[i];
//             int countForAdjGrid = cell->affectedParticlesWeightID[i];
//             float w = p->weights[countForAdjGrid];
//             Vector3f wGrad = p->weight_gradient[countForAdjGrid];
//             Matrix3f energyD = p->energyDerivative();
//             if (w > BSPLINE_EPSILON)
//             {
//                 force += energyD * wGrad;
//                 // std::cout << " force + energyD is " << energyD <<
//                 std::endl;
//                 // std::cout << " force + wGrad is " << wGrad << std::endl;
//             }
//         }
//         // std::cout << " force is " << force << std::endl;
//         // std::cout << " g is " << gravity << std::endl;
//         // std::cout << " mass is " << cell->mass << std::endl;
//         // force+gravity->a->dv->vStar
//         cell->velocityStar =
//             cell->velocity + deltaT * (gravity - force / cell->mass);

//         // if (isnan(cell->velocityStar.norm()))
//         // {
//         //     std::cout << " cell velocity star is " << cell->velocityStar
//         //               << std::endl;
//         //     std::cout << " cell velocity old is " << cell->velocity
//         //               << std::endl;
//         //     std::cout << " cell velocity mass " << cell->mass <<
//         std::endl;
//         //     std::cout << " cell id is " << cell->i << ", " << cell->j <<
//         ", "
//         //               << cell->k << std::endl;
//         //     exit(0);
//         // }
//         // if (cell->i == 0 && cell->j == 0 && cell->k == 16)
//         // {
//         //     std::cout << " cell velocity star is " << cell->velocityStar
//         //               << std::endl;
//         //     std::cout << " cell velocity old is " << cell->velocity
//         //               << std::endl;
//         //     std::cout << " cell velocity mass " << cell->mass <<
//         std::endl;
//         //     std::cout << " cell id is " << cell->i << ", " << cell->j <<
//         ", "
//         //               << cell->k << std::endl;
//         //     // exit(0);
//         // }
//     }

//     // TODO the above loop shall be combined with collision Grid when Parfor
//     gridCollisionTest();
// }

// #if ENABLE_IMPLICIT
// Solve linear system for implicit velocities
void GridMesh::implicitUpdateGridVelocity()
{
    // TODO
}

void GridMesh::recomputeImplicitForces()
{
    // TODO
}
// #endif

// Map grid velocities back to particles
void GridMesh::mapVelocityToSPS() const
{
#pragma omp parallel for
    for (int i = 0; i < SPS->particles.size(); i++)
    {
        // map back to particles
        {
            SnowParticle* p = SPS->particles[i];
            // position
            Vector3f pos = p->position - bbox.pMin;
            // We calculate Vpic and VFLIP velocities separately
            Vector3f Vpic = Vector3f::Zero();
            Vector3f Vflip = p->velocity;
            // Also keep track of velocity gradient
            Matrix3f& grad = p->velocityGradient;
            grad.setZero();
            // VISUALIZATION PURPOSES ONLY:
            // Recompute density
            p->density = 0;

            for (int countForAdjGrid = 0, px = (int)(pos.x() / cellSize.x()),
                     ipx = px - 1, pxEnd = px + 3;
                 ipx < pxEnd; ipx++)
            {
                if (ipx >= 0 && ipx < cellNum.x())
                {
                    for (int py = (int)(pos.y() / cellSize.y()), ipy = py - 1,
                             pyEnd = py + 3;
                         ipy < pyEnd; ipy++)
                    {
                        if (ipy >= 0 && ipy < cellNum.y())
                        {
                            for (int pz = (int)(pos.z() / cellSize.z()),
                                     ipz = pz - 1, pzEnd = pz + 3;
                                 ipz < pzEnd; ipz++)
                            {
                                if (ipz >= 0 && ipz < cellNum.z())
                                {
                                    float w = p->weights[countForAdjGrid];
                                    Vector3f wGrad =
                                        p->weight_gradient[countForAdjGrid];
                                    if (w > BSPLINE_EPSILON)
                                    {
                                        GridCell* cell = cells[(
                                            int)(ipx * cellNum.y() *
                                                     cellNum.z() +
                                                 ipy * cellNum.z() + ipz)];
                                        // Particle in cell
                                        Vpic += cell->velocityStar * w;
                                        // Fluid implicit particle
                                        Vflip += (cell->velocityStar -
                                                  cell->velocity) *
                                                 w;
                                        // Velocity gradient
                                        grad += cell->velocityStar *
                                                wGrad.transpose();
                                        // VISUALIZATION ONLY: Update density
                                        p->density += w * cell->mass;
                                        // test
                                        // if (isnan(Vflip.norm()))
                                        // {
                                        //     std::cout << " vel flip is " <<
                                        //     Vflip
                                        //               << std::endl;
                                        //     std::cout << " cell active "
                                        //               << cell->active <<
                                        //               std::endl;
                                        //     std::cout << " cell id is " <<
                                        //     ipx
                                        //               << ", " << ipy << ", "
                                        //               << ipz
                                        //               << std::endl;
                                        //     std::cout << " cell id is " <<
                                        //     cell->i
                                        //               << ", " << cell->j <<
                                        //               ", "
                                        //               << cell->k <<
                                        //               std::endl;
                                        //     std::cout << " cell velocity star
                                        //     is
                                        //     "
                                        //               << cell->velocityStar
                                        //               << std::endl;
                                        //     std::cout << " cell velocity old
                                        //     is "
                                        //               << cell->velocity
                                        //               << std::endl;
                                        //     exit(0);
                                        // }
                                        // if (isnan(Vpic.norm()))
                                        // {
                                        //     std::cout << " vel pic is " <<
                                        //     Vpic
                                        //               << std::endl;
                                        // }
                                    }
                                }
                                // do it or not, it has goto next cell
                                countForAdjGrid++;
                            }
                        }
                        else
                        {
                            // do nothing and skip this z edge
                            countForAdjGrid += 4;
                        }
                    }
                }
                else
                {
                    // do nothing and skip this yz face
                    countForAdjGrid += 16;
                }
            }
            // Final velocity is a linear combination of Vpic and VFLIP
            // components
            p->velocity = Vflip * FLIP_PERCENT + Vpic * (1. - FLIP_PERCENT);
            // if (isnan(Vflip.norm()))
            // {
            //     std::cout << " vel flip is " << Vflip << std::endl;
            // }
            // if (isnan(Vpic.norm()))
            // {
            //     std::cout << " vel pic is " << Vpic << std::endl;
            // }
            // VISUALIZATION: Update density
            p->density /= eachCellVolume;
        }
        // collision test
        {
            // collision test
            SnowParticle* p = SPS->particles[i];
            Vector3f tempPos = p->position + deltaT * p->velocity;
            // x dir
            if (tempPos.x() > bbox.pMax.x() || tempPos.x() < bbox.pMin.x())
            {
                Vector3f VSolid(0, 0, 0);
                Vector3f VRel = p->velocity - VSolid;
                Vector3f Vn = Vector3f(VRel.x(), 0, 0);
                Vector3f Vt = Vector3f(0, VRel.y(), VRel.z());
                float vn = Vn.norm();
                float vt = Vt.norm();
                float mu = p->m->sticky;
                if (vt < mu * vn || vt < 1.e-12)
                {
                    VRel.setZero();
                }
                else
                {
                    VRel = (1. - mu * vn / vt) * Vt;
                }
                p->velocity = VRel + VSolid;
            }
            // y dir
            if (tempPos.y() > bbox.pMax.y() || tempPos.y() < bbox.pMin.y())
            {
                Vector3f VSolid(0, 0, 0);
                Vector3f VRel = p->velocity - VSolid;
                Vector3f Vn = Vector3f(0, VRel.y(), 0);
                Vector3f Vt = Vector3f(VRel.x(), 0, VRel.z());
                float vn = Vn.norm();
                float vt = Vt.norm();
                float mu = p->m->sticky;
                if (vt < mu * vn || vt < 1.e-12)
                {
                    VRel.setZero();
                }
                else
                {
                    VRel = (1. - mu * vn / vt) * Vt;
                }
                p->velocity = VRel + VSolid;
            }
            // z dir
            if (tempPos.z() > bbox.pMax.z() || tempPos.z() < bbox.pMin.z())
            {
                Vector3f VSolid(0, 0, 0);
                Vector3f VRel = p->velocity - VSolid;
                Vector3f Vn = Vector3f(0, 0, VRel.z());
                Vector3f Vt = Vector3f(VRel.x(), VRel.y(), 0);
                float vn = Vn.norm();
                float vt = Vt.norm();
                float mu = p->m->sticky;
                if (vt < mu * vn || vt < 1.e-12)
                {
                    VRel.setZero();
                }
                else
                {
                    VRel = (1. - mu * vn / vt) * Vt;
                }
                p->velocity = VRel + VSolid;
            }
            // TODO there are more collisions
            // Finish collision test
            // Now push Particles Postion
            p->updatePos();
            // update pure elastic deformation
            p->updatePureElasticGradient();
            // push excessive part to plastic deformation
            p->updateCombinedPElasticGradient();
        }
    }
}

// void GridMesh::gridCollisionTest()
// {
//     for (int iCell = 0; iCell < totalEffectiveCellNum; iCell++)
//     {
//         // int i = iCell / (cellNum.y() * cellNum.z());
//         // int j = (iCell % (cellNum.y() * cellNum.z())) / cellNum.z();
//         // int k = iCell % cellNum.z();

//         // int iCell = i * cellNum.y() * cellNum.z() + j * cellNum.z() + k;
//         GridCell* cell = cells[iCell];
//         int i = cell->i;
//         int j = cell->j;
//         int k = cell->k;
//         // collision test
//         Vector3f oldCellPos((i + 0.5) * cellSize.x(), (j + 0.5) *
//         cellSize.y(),
//                             (k + 0.5) * cellSize.z());
//         Vector3f tempPos = oldCellPos + deltaT * cell->velocityStar;
//         // TODO consider wrap this fixing vel as a function
//         // x dir
//         if (tempPos.x() > bbox.pMax.x() || tempPos.x() < bbox.pMin.x())
//         {
//             Vector3f VSolid(0, 0, 0);
//             Vector3f VRel = cell->velocityStar - VSolid;
//             Vector3f Vn = Vector3f(VRel.x(), 0, 0);
//             Vector3f Vt = Vector3f(0, VRel.y(), VRel.z());
//             float vn = Vn.norm();
//             float vt = Vt.norm();
//             // TODO this mu calculation should be cleverer
//             float mu = SPS->particles[0]->m->sticky;
//             if (vt < mu * vn || vt < 1.e-12)
//             {
//                 VRel.setZero();
//             }
//             else
//             {
//                 VRel = (1. - mu * vn / vt) * Vt;
//             }
//             cell->velocityStar = VRel + VSolid;
//             // std::cout << " cell velocity star is " << cell->velocityStar
//             //           << std::endl;
//         }
//         // y dir
//         if (tempPos.y() > bbox.pMax.y() || tempPos.y() < bbox.pMin.y())
//         {
//             Vector3f VSolid(0, 0, 0);
//             Vector3f VRel = cell->velocityStar - VSolid;
//             Vector3f Vn = Vector3f(0, VRel.y(), 0);
//             Vector3f Vt = Vector3f(VRel.x(), 0, VRel.z());
//             float vn = Vn.norm();
//             float vt = Vt.norm();
//             // TODO this mu calculation should be cleverer
//             float mu = SPS->particles[0]->m->sticky;
//             if (vt < mu * vn || vt < 1.e-12)
//             {
//                 VRel.setZero();
//             }
//             else
//             {
//                 VRel = (1. - mu * vn / vt) * Vt;
//             }
//             cell->velocityStar = VRel + VSolid;
//             // std::cout << " cell velocity star is " << cell->velocityStar
//             //           << std::endl;
//         }
//         // z dir
//         if (tempPos.z() > bbox.pMax.z() || tempPos.z() < bbox.pMin.z())
//         {
//             Vector3f VSolid(0, 0, 0);
//             Vector3f VRel = cell->velocityStar - VSolid;
//             Vector3f Vn = Vector3f(0, 0, VRel.z());
//             Vector3f Vt = Vector3f(VRel.x(), VRel.y(), 0);
//             float vn = Vn.norm();
//             float vt = Vt.norm();
//             // TODO this mu calculation should be cleverer
//             float mu = SPS->particles[0]->m->sticky;
//             if (vt < mu * vn || vt < 1.e-12)
//             {
//                 VRel.setZero();
//             }
//             else
//             {
//                 VRel = (1. - mu * vn / vt) * Vt;
//             }
//             cell->velocityStar = VRel + VSolid;
//             // std::cout << " cell velocity star is " << cell->velocityStar
//             //           << std::endl;
//         }
//         // TODO there are more collisions
//         // Finish collision test

//         // if (isnan(cell->velocityStar.norm()))
//         // {
//         //     std::cout << " cell velocity star is " << cell->velocityStar
//         //               << std::endl;
//         //     std::cout << " cell velocity old is " << cell->velocity
//         //               << std::endl;
//         //     std::cout << " cell velocity mass " << cell->mass <<
//         std::endl;
//         //     std::cout << " cell id is " << cell->i << ", " << cell->j <<
//         ", "
//         //               << cell->k << std::endl;
//         //     exit(0);
//         // }
//     }
//     // if (cell->i == 0 && cell->j == 0 && cell->k == 16)
//     // {
//     //     std::cout << " cell velocity star is " << cell->velocityStar
//     //               << std::endl;
//     //     std::cout << " cell velocity old is " << cell->velocity <<
//     std::endl;
//     //     std::cout << " cell velocity mass " << cell->mass << std::endl;
//     //     std::cout << " cell id is " << cell->i << ", " << cell->j << ", "
//     //               << cell->k << std::endl;
//     //     // exit(0);
//     // }
// }

// void GridMesh::particleCollisionTest() const
// {
//     // Try to write all particles' update here
//     // instead of creating new for loops later, better for parfor
//     // Now collision only happens in the bbox of the grid

//     for (int i = 0; i < SPS->particles.size(); i++)
//     {
//         // collision test
//         SnowParticle* p = SPS->particles[i];
//         Vector3f tempPos = p->position + deltaT * p->velocity;
//         // x dir
//         if (tempPos.x() > bbox.pMax.x() || tempPos.x() < bbox.pMin.x())
//         {
//             Vector3f VSolid(0, 0, 0);
//             Vector3f VRel = p->velocity - VSolid;
//             Vector3f Vn = Vector3f(VRel.x(), 0, 0);
//             Vector3f Vt = Vector3f(0, VRel.y(), VRel.z());
//             float vn = Vn.norm();
//             float vt = Vt.norm();
//             float mu = p->m->sticky;
//             if (vt < mu * vn || vt < 1.e-12)
//             {
//                 VRel.setZero();
//             }
//             else
//             {
//                 VRel = (1. - mu * vn / vt) * Vt;
//             }
//             p->velocity = VRel + VSolid;
//         }
//         // y dir
//         if (tempPos.y() > bbox.pMax.y() || tempPos.y() < bbox.pMin.y())
//         {
//             Vector3f VSolid(0, 0, 0);
//             Vector3f VRel = p->velocity - VSolid;
//             Vector3f Vn = Vector3f(0, VRel.y(), 0);
//             Vector3f Vt = Vector3f(VRel.x(), 0, VRel.z());
//             float vn = Vn.norm();
//             float vt = Vt.norm();
//             float mu = p->m->sticky;
//             if (vt < mu * vn || vt < 1.e-12)
//             {
//                 VRel.setZero();
//             }
//             else
//             {
//                 VRel = (1. - mu * vn / vt) * Vt;
//             }
//             p->velocity = VRel + VSolid;
//         }
//         // z dir
//         if (tempPos.z() > bbox.pMax.z() || tempPos.z() < bbox.pMin.z())
//         {
//             Vector3f VSolid(0, 0, 0);
//             Vector3f VRel = p->velocity - VSolid;
//             Vector3f Vn = Vector3f(0, 0, VRel.z());
//             Vector3f Vt = Vector3f(VRel.x(), VRel.y(), 0);
//             float vn = Vn.norm();
//             float vt = Vt.norm();
//             float mu = p->m->sticky;
//             if (vt < mu * vn || vt < 1.e-12)
//             {
//                 VRel.setZero();
//             }
//             else
//             {
//                 VRel = (1. - mu * vn / vt) * Vt;
//             }
//             p->velocity = VRel + VSolid;
//         }
//         // TODO there are more collisions
//         // Finish collision test
//         // Now push Particles Postion
//         p->updatePos();
//         // update pure elastic deformation
//         p->updatePureElasticGradient();
//         // push excessive part to plastic deformation
//         p->updateCombinedPElasticGradient();
//     }
// }