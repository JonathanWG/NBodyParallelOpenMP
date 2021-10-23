#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

namespace NBODY {

    typedef struct {
        double x, y, z;
        double mass;
        } Particle;
    typedef struct {
        double xold, yold, zold;
        double fx, fy, fz;
        } ParticleV;

    static long seed = DEFAULT;

    class NBody {

    public:

        static double Random(void);
        void InitParticles(Particle[], ParticleV [], int );
        double ComputeForces( Particle [], Particle [], ParticleV [], int );
        double ComputeNewPos( Particle [], ParticleV [], int, double);

        void ParallelInitParticles(Particle[], ParticleV [], int );
        double ParallelComputeForces( Particle [], Particle [], ParticleV [], int );
        double ParallelComputeNewPos( Particle [], ParticleV [], int, double);

    };
}
