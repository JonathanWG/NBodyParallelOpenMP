#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include "nbody_parallel/nbody_parallel.hpp"
#define NUMBER_OF_TESTS 20

using namespace std;
using namespace NBODY;
bool log_details = false;


std::ifstream in_ifs;
typedef struct unitTest{
    int num_part;
    int num_inter;
    int n_threads;

    unitTest(int n1, int n2,int n3) {
        num_part = n1;
        num_inter = n2;
        n_threads = n3;
    }
} unitTest;

const auto p1 = std::chrono::system_clock::now();
std::string filename;
std::ofstream log_file;
clock_t clk_start, clk_end;

double serial_test(NBody &nbdObj,int npart,int cnt) {
    /* number of particles, number of times in loop */
     clk_start = clock();

    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int          i, j;
    double      sim_t;       /* Simulation time */


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
    clk_end = clock();
    double time_taken = double(clk_end - clk_start) / double(CLOCKS_PER_SEC);

    if(log_details)
        log_file << "---- Serial Execution  ----- " << endl;

    if(log_details)
        log_file << "\n time spended : " << time_taken << endl;
    if(log_details)
        for (i=0; i<npart; i++)
          log_file << ("%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z) << endl;

   if(log_details) log_file << "------------------------   ";
    return time_taken;
}



double parallel_test(NBody &nbdObj,int npart,int cnt) { /* number of particles, number of times in loop */
    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int i, j;
    double      sim_t;       /* Simulation time */

    clk_start = clock();

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


    clk_end = clock();
    double time_taken = double(clk_end - clk_start) / double(CLOCKS_PER_SEC);

    if(log_details) log_file << "---- Parallel Execution  ----- " << endl;
   if(log_details) log_file << "\n time spended : " << time_taken << endl;

    if(log_details)
        for (i=0; i<npart; i++)
          log_file << ("%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z) << endl;

    if(log_details)log_file << "------------------------\n\n";
    return time_taken;
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
    int count = 0;
    filename = "_testOutput.txt";


    while (std::getline(in_ifs, line))
    {
        cout << "Reading .in line " << test_counter << endl;
        std::istringstream iss(line);

        int n_part, n_inter,n_threads;

        if (!(iss >> n_part >> n_inter >> n_threads)) { break; } // error
        tests_vect.push_back(unitTest(n_part,n_inter,n_threads));
        test_counter++;
    }

    log_file.open(filename,std::ios_base::out);
    log_file.precision(5);

    log_file << "-- TEST BATERY. NUMBER OF TESTS : " << test_counter << endl << endl;

    do {
       //for each test do 20 runs
       double global_serial_time_sum  = 0.0, global_parallel_time_sum = 0.0,global_speedup = 0.0;
       for(int i =0; i< NUMBER_OF_TESTS;i++) {
            if(log_details)
                log_file << "--- UNIT TEST   "<< count <<   "-----" <<  endl;

            nbdObj.setNumThreads(tests_vect.at(count).n_threads);
            if(log_details)
                log_file << "--------------------------------------\n" <<  endl;
            double serial_time = serial_test(nbdObj,tests_vect.at(count).num_part,tests_vect.at(count).num_inter);

            double parallel_time = parallel_test(nbdObj,tests_vect.at(count).num_part,tests_vect.at(count).num_inter);

            global_serial_time_sum += serial_time;
            global_parallel_time_sum += parallel_time;
        }
       log_file << "--------------------------------------\n" <<  endl;
       log_file  << "n_part << " << tests_vect.at(count).num_part<<  "  |  n_part <<  " << tests_vect.at(0).num_inter
                 << " |  nthreads  " << tests_vect.at(count).n_threads<< endl;
       log_file << "SERIAL average of 20 executions this test : " << (global_serial_time_sum/NUMBER_OF_TESTS)<< endl;
       log_file << "PARALLEL average of 20 executions this test : " << (global_parallel_time_sum/NUMBER_OF_TESTS)<< endl;
       log_file << "SPEEDUP average of 20 executions this test : " << (global_serial_time_sum/global_parallel_time_sum)<< endl;
        count++;
        log_file << "--------------------------------------\n" <<  endl;
     }while(count < test_counter);
}
