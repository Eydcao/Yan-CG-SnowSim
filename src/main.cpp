
// #include "Renderer.hpp"
// #include "Scene.hpp"
// #include "Triangle.hpp"
// #include "Sphere.hpp"
// #include "Vector.hpp"
// #include "global.hpp"
// #include <chrono>
#include "Rectangular.hpp"
#include "SnowParticle.hpp"
#include "Sphere.hpp"
#include "Triangle.hpp"
#include <GL/glut.h>
#include <iostream>

// Main function for snow simulation
int main(int argc, char** argv)
{
    // tests
    {
        std::cout << "hello snow simulation tests" << std::endl;
    }
    {
        std::cout << "test p generation in a sphere" << std::endl;
        SnowParticleMaterial m;
        m.lNumDensity = 2;
        Sphere Omega(Vector3f(0.0, 0.0, 0.0), 10.0);
        SnowParticleSet spSet;
        spSet.addParticlesInAShape(&Omega, &m);
        assert(spSet.particles.size() == 30976);
        for (auto& anyParticle : spSet.particles)
        {
            assert(anyParticle->m == &m);
        }
    }
    {
        std::cout << "test p generation in a box" << std::endl;
        SnowParticleMaterial m;
        m.lNumDensity = 2;
        Rectangular Box(Vector3f(0.0, 0.0, 0.0), Vector3f(10.0, 10.0, 10.0));
        SnowParticleSet spSet;
        spSet.addParticlesInAShape(&Box, &m);
        assert(spSet.particles.size() == 8000);
        for (auto& anyParticle : spSet.particles)
        {
            assert(anyParticle->m == &m);
        }
    }
    {
        std::cout << "test p generation in a closed tri mesh" << std::endl;
        SnowParticleMaterial m;
        m.lNumDensity = 10;
        MeshTriangle cow("../media/spot_triangulated_good.obj");
        SnowParticleSet spSet;
        spSet.addParticlesInAShape(&cow, &m);
        std::cout << " the vol ratio is "
                  << cow.getVolume() / (cow.getBounds().volume()) << std::endl;
        std::cout << " the size ratio is "
                  << spSet.particles.size() / 9. / 16. / 17. << std::endl;
        // assert(spSet.particles.size() == 8000);
        for (auto& anyParticle : spSet.particles)
        {
            assert(anyParticle->m == &m);
        }
    }
    {
        std::cout << "all snow simulation tests passed" << std::endl;
    }
    // end tests

    return 0;
}