#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <chrono>
#include <string>
#include <stdio.h>
#include <omp.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789
#define PAD 8

namespace NBODY {
static long seed = DEFAULT;

using namespace std;

    typedef struct {
        double x, y, z;
        double mass;
        } Particle;
    typedef struct {
        double xold, yold, zold;
        double fx, fy, fz;
        } ParticleV;


    class NBody {

    public:
        NBody();
        static double Random(void);
        void InitParticles(Particle[], ParticleV [], int );
        double ComputeForces( Particle [], Particle [], ParticleV [], int );
        double ComputeNewPos( Particle [], ParticleV [], int, double);
        //parallel
        void ParallelInitParticles(Particle[], ParticleV [], int );
        double ParallelComputeForces( Particle [], Particle [], ParticleV [], int );
        double ParallelComputeNewPos( Particle [], ParticleV [], int, double);
        //utils
        static Particle** PadParticle(Particle parts[],int pad_size);
        static ParticleV** PadParticleV(ParticleV p_v[],int pad_size);
        void setNumThreads(int n);

    private:
        int num_threads;

    };
}
