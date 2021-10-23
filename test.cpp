#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "nbody_parallel/nbody_parallel.hpp"

using namespace std;
using namespace NBODY;


std::ifstream in_ifs;
typedef struct unitTest{
    int num_part;
    int num_inter;
    unitTest(int n1, int n2) {
        num_part = n1;
        num_inter = n2;
    }
} unitTest;


int serial_test(NBody &nbdObj,int npart,int cnt) { /* number of particles, number of times in loop */
    double time;
    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int          i, j;
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



int parallel_test(NBody &nbdObj,int npart,int cnt) { /* number of particles, number of times in loop */
    double time;
    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int i, j;
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

    in_ifs.open("nbody.in",std::ifstream::in);
    if(!in_ifs){
        cout << "Error while trying read nbody.in file." << endl;
        return 0;
    }

    std::string line;
    int test_counter = 0;
    std::vector<unitTest> tests_vect;

    while (std::getline(in_ifs, line))
    {
        cout << "Reading .in line " << test_counter << endl;
        std::istringstream iss(line);
        int n_part, n_inter;
        if (!(iss >> n_part >> n_inter)) { break; } // error
        tests_vect.push_back(unitTest(n_part,n_inter));
        cout << "n_part << " << n_part <<  "  |  n_part <<  " << n_part  << endl;
        test_counter++;
    }


    serial_test(nbdObj,tests_vect.at(0).num_part,tests_vect.at(0).num_inter);

    parallel_test(nbdObj,tests_vect.at(0).num_part,tests_vect.at(0).num_inter);

}
