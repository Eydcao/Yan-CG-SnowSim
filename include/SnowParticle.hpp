#ifndef SNOWSIM_SNOWPARTICLES
#define SNOWSIM_SNOWPARTICLES

#include "Shape.hpp"
#include <iostream>
using namespace Eigen;

class SnowParticleMaterial
{
   private:
   public:
    // TODO add more related data inside later, give them defualt values
    // TODO fig out if this will be updated or only volume needs to
    float initialDensity = 400.;
    // line number density when generating particles
    // initial value 138 = 0.0072 diameter
    int lNumDensity = 138;
    float criticalStress = 1. - 2.5e-2;
    float criticalStretch = 1. + 7.5e-3;
    float hardening = 10.;
    float youngsModule = 1.4e5;
    float PoissonsRatio = .2;
    // Collision stickiness (lower = stickier)
    float sticky = .9;
    // Lame parameters
    float mu;
    float lambda;

    // etc

    SnowParticleMaterial();
    ~SnowParticleMaterial();
};

class SnowParticle
{
   private:
   public:
    float volume;
    float mass;
    // TODO fig out if this will be updated or only volume needs to
    float density;
    Vector3f position;
    Vector3f velocity;
    Matrix3f velocityGradient;
    Matrix3f deformationGradientElastic;
    Matrix3f deformationGradientPlastic;
    // Cached SVD's for elastic deformation gradient
    Matrix3f SVDW;
    Matrix3f SVDV;
    Vector3f SVDE;
    // Cached polar decomposition
    // TODO fig out in 3d if the polar phi is also needed
    Matrix3f polarR;
    Matrix3f polarTheta;
    Matrix3f polarPhi;
    // Grid interpolation weights
    // TODO fig out if grid position can be replaced by id or sth
    Vector3f grid_position;
    // TODO fig out what is the size of adj grids in 3d
    Vector3f weight_gradient[16];
    float weights[16];

    SnowParticleMaterial* m;

    SnowParticle();
    SnowParticle(const Vector3f& pos, SnowParticleMaterial* material);
    ~SnowParticle();

    // Update position, based on velocity
    void updatePos();
    // Update deformation gradient
    void updatePureElasticGradient();
    void updateCombinedPElasticGradient();
    // Compute stress tensor
    const Matrix3f energyDerivative();

    // Computes stress force delta, for implicit velocity update
    const Vector3f deltaForce(const Vector2f& u, const Vector2f& weight_grad);
};

class SnowParticleSet
{
   private:
   public:
    std::vector<SnowParticle*> particles;
    float maxVelocity;

    SnowParticleSet();
    ~SnowParticleSet();
    void addParticle(SnowParticle* sp);
    void addParticle(const Vector3f& pos, SnowParticleMaterial* m);
    void addParticlesInAShape(Shape* s, SnowParticleMaterial* m);
    // void appendSet(const SnowParticleSet& anotherSet);
    // inline SnowParticleSet unionSet(const SnowParticleSet& set1,
    //                                 const SnowParticleSet& set2);
    // inline SnowParticleSet unionSet(const std::vector<SnowParticleSet>&
    // sets);
    void update();
};

#endif  // SNOWSIM_SNOWPARTICLES
