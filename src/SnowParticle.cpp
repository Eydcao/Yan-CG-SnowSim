#include "SnowParticle.hpp"

SnowParticleMaterial::SnowParticleMaterial()
{
    mu = youngsModule / (2. + 2. * PoissonsRatio);
    lambda = youngsModule * PoissonsRatio /
             ((1. + PoissonsRatio) * (1. - 2. * PoissonsRatio));
}

SnowParticleMaterial::~SnowParticleMaterial()
{
}

SnowParticle::SnowParticle() : m(nullptr)
{
}

SnowParticle::SnowParticle(const Vector3f& pos, SnowParticleMaterial* material)
    : position(pos), m(material)
{
}

SnowParticle::~SnowParticle()
{
}

SnowParticleSet::SnowParticleSet() : particles()
{
}

SnowParticleSet::~SnowParticleSet()
{
    // TODO check if this destroy way is okay
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