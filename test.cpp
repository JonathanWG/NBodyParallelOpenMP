#include <iostream>
#include <omp.h>
#include <stdio.h>
#include "nbody_parallel/nbody_parallel.hpp"

using namespace std;
using namespace NBODY;



int serial_test(NBody &nbdObj) {
    double time;
    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int         npart, i, j;
    int         cnt;         /* number of times in loop */
    double      sim_t;       /* Simulation time */
    int tmp;
    tmp = fscanf(stdin,"%d\n",&npart);
    tmp = fscanf(stdin,"%d\n",&cnt);
/* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle)*npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);
/* Generate the initial values */
    nbdObj.InitParticles( particles, pv, npart);
    sim_t = 0.0;

    while (cnt--) {
      double max_f;
      /* Compute forces (2D only) */
      max_f = nbdObj.ComputeForces( particles, particles, pv, npart );
      /* Once we have the forces, we compute the changes in position */
      sim_t += nbdObj.ComputeNewPos( particles, pv, npart, max_f);
    }
    for (i=0; i<npart; i++)
      fprintf(stdout,"%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z);
    return 0;
}



int parallel_test(NBody &nbdObj) {
    double time;
    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int         npart, i, j;
    int         cnt;         /* number of times in loop */
    double      sim_t;       /* Simulation time */
    int tmp;
    tmp = fscanf(stdin,"%d\n",&npart);
    tmp = fscanf(stdin,"%d\n",&cnt);
/* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle)*npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);
/* Generate the initial values */
    nbdObj.ParallelInitParticles( particles, pv, npart);
    sim_t = 0.0;

    while (cnt--) {
      double max_f;
      /* Compute forces (2D only) */
      max_f = nbdObj.ParallelComputeForces( particles, particles, pv, npart );
      /* Once we have the forces, we compute the changes in position */
      sim_t += nbdObj.ParallelComputeNewPos( particles, pv, npart, max_f);
    }
    for (i=0; i<npart; i++)
      fprintf(stdout,"%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z);
    return 0;
}


int main()
{
    NBody nbdObj;

    serial_test(nbdObj);

    parallel_test(nbdObj);

}
