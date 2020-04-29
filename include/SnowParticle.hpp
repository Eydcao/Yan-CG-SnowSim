#ifndef SNOWSIM_SNOWPARTICLES
#define SNOWSIM_SNOWPARTICLES

#include "Shape.hpp"
#include <iostream>
using namespace Eigen;

struct SnowParticleMaterial
{
    // TODO add more related data inside later, give them defualt values
    float density = 1.;
    float youngsModule = 1.;
    int lNumDensity = 1;
    // etc
};

class SnowParticle
{
   private:
   public:
    Vector3f position;
    Vector3f velocity;
    Vector3f deformationGradient;

    SnowParticleMaterial* m;

    SnowParticle();
    SnowParticle(const Vector3f& pos, SnowParticleMaterial* material);
    ~SnowParticle();
};

SnowParticle::SnowParticle()
    : position(), velocity(), deformationGradient(), m(nullptr)
{
}

SnowParticle::SnowParticle(const Vector3f& pos, SnowParticleMaterial* material)
    : position(pos), m(material)
{
}

SnowParticle::~SnowParticle()
{
}

class SnowParticleSet
{
   private:
   public:
    std::vector<SnowParticle*> particles;

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
};

SnowParticleSet::SnowParticleSet() : particles()
{
}

SnowParticleSet::~SnowParticleSet()
{
    for (auto& oneParticle : particles)
    {
        delete (oneParticle);
    }
}

void SnowParticleSet::addParticle(SnowParticle* sp)
{
    particles.push_back(sp);
}

void SnowParticleSet::addParticle(const Vector3f& pos, SnowParticleMaterial* m)
{
    SnowParticle* sp = new SnowParticle(pos, m);
    particles.push_back(sp);
}

void SnowParticleSet::addParticlesInAShape(Shape* s, SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    // std::cout << " size is " << temp << std::endl;
    if (temp > 0)
    {
        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, m);
        }
    }
}

#endif  // SNOWSIM_SNOWPARTICLES
