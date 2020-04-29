#ifndef SNOWSIM_SNOWPARTICLES
#define SNOWSIM_SNOWPARTICLES

#include "Shape.hpp"
using namespace Eigen;

struct SnowParticleMaterial
{
    // TODO add more related data inside later, give them defualt values
    float density = 1.;
    float youngsModule = 1.;
    int lNumDensity = 10;
    // etc
};

class SnowParticle
{
   private:
    SnowParticleMaterial* m;
    Vector3f position;
    Vector3f velocity;
    Vector3f deformationGradient;

   public:
    SnowParticle();
    SnowParticle(const Vector3f& pos, SnowParticleMaterial* material);
    ~SnowParticle();
};

SnowParticle::SnowParticle() : position(), m(nullptr)
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
    std::vector<SnowParticle*> particles;

   public:
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

SnowParticleSet::SnowParticleSet()
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
    if (s->generateParticlesInside(m->lNumDensity, tempPos) > 0)
    {
        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, m);
        }
    }
}

#endif  // SNOWSIM_SNOWPARTICLES
