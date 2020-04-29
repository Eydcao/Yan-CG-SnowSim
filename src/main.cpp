
// #include "Renderer.hpp"
// #include "Scene.hpp"
// #include "Triangle.hpp"
// #include "Sphere.hpp"
// #include "Vector.hpp"
// #include "global.hpp"
// #include <chrono>
#include "SnowParticle.hpp"
#include "Sphere.hpp"
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
        std::cout << "all snow simulation tests passed" << std::endl;
    }
    // end tests
    return 0;
}