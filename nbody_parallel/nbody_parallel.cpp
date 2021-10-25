#include "nbody_parallel.hpp"
using namespace NBODY;
using namespace std;

 NBody::NBody() {

}


double  NBody::Random(void)
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed % Q) - R * (seed / Q);
  if (t > 0)
    seed = t;
  else
    seed = t + MODULUS;
  return ((double) seed / MODULUS);
}

void NBody::InitParticles( Particle particles[], ParticleV pv[], int npart )
{
    int i;
//    cout << " " << endl;
    for (i=0; i<npart; i++) {
    particles[i].x	  = this->Random();
    particles[i].y	  = this->Random();
    particles[i].z	  = this->Random();
    particles[i].mass = 1.0;
    pv[i].xold	  = particles[i].x;
    pv[i].yold	  = particles[i].y;
    pv[i].zold	  = particles[i].z;
    pv[i].fx	  = 0;
    pv[i].fy	  = 0;
    pv[i].fz	  = 0;
    }
}

double NBody::ComputeForces( Particle myparticles[], Particle others[], ParticleV pv[], int npart )
{
  double max_f;
  int i;
  max_f = 0.0;
  for (i=0; i<npart; i++) {
    int j;
    double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
    rmin = 100.0;
    xi   = myparticles[i].x;
    yi   = myparticles[i].y;
    fx   = 0.0;
    fy   = 0.0;
    for (j=0; j<npart; j++) {
      rx = xi - others[j].x;
      ry = yi - others[j].y;
      mj = others[j].mass;
      r  = rx * rx + ry * ry;
      /* ignore overlap and same particle */
      if (r == 0.0) continue;
      if (r < rmin) rmin = r;
      r  = r * sqrt(r);
      fx -= mj * rx / r;
      fy -= mj * ry / r;
    }
    pv[i].fx += fx;
    pv[i].fy += fy;
    fx = sqrt(fx*fx + fy*fy)/rmin;
    if (fx > max_f) max_f = fx;
  }
  return max_f;
}


double NBody::ComputeNewPos( Particle particles[], ParticleV pv[], int npart, double max_f)
{
  int i;
  double a0, a1, a2;
  static double dt_old = 0.001, dt = 0.001;
  double dt_new;
  a0	 = 2.0 / (dt * (dt + dt_old));
  a2	 = 2.0 / (dt_old * (dt + dt_old));
  a1	 = -(a0 + a2);
  for (i=0; i<npart; i++) {
    double xi, yi;
    xi	           = particles[i].x;
    yi	           = particles[i].y;
    particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
    particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
    pv[i].xold     = xi;
    pv[i].yold     = yi;
    pv[i].fx       = 0;
    pv[i].fy       = 0;
  }
  dt_new = 1.0/sqrt(max_f);
  /* Set a minimum: */
  if (dt_new < 1.0e-6) dt_new = 1.0e-6;
  /* Modify time step */
  if (dt_new < dt) {
    dt_old = dt;
    dt     = dt_new;
  }
  else if (dt_new > 4.0 * dt) {
    dt_old = dt;
    dt    *= 2.0;
  }
  return dt_old;
}


//parallel


void NBody::ParallelInitParticles( Particle particles[], ParticleV pv[], int npart )
{
    int i;
    #pragma omp parallel for num_threads(num_threads)
        for (i=0; i<npart; i++) {
            particles[i].x	  = this->Random();
            particles[i].y	  = this->Random();
            particles[i].z	  = this->Random();
            particles[i].mass = 1.0;
            pv[i].xold	  = particles[i].x;
            pv[i].yold	  = particles[i].y;
            pv[i].zold	  = particles[i].z;
            pv[i].fx	  = 0;
            pv[i].fy	  = 0;
            pv[i].fz	  = 0;
        }
}

double NBody::ParallelComputeForces( Particle myparticles[], Particle others[], ParticleV pv[], int npart )
{
  double max_f;
  int i;
  max_f = 0.0;
#pragma omp parallel for reduction(+ : max_f) collapse(2) num_threads(num_threads)

  for (i=0; i<npart; i++) {
        int j;
        double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
        rmin = 100.0;
        xi   = myparticles[i].x;
        yi   = myparticles[i].y;
        fx   = 0.0;
        fy   = 0.0;
        for (j=0; j<npart; j++) {
          rx = xi - others[j].x;
          ry = yi - others[j].y;
          mj = others[j].mass;
          r  = rx * rx + ry * ry;
          /* ignore overlap and same particle */
          if (r == 0.0) continue;
          if (r < rmin) rmin = r;
          r  = r * sqrt(r);
          fx -= mj * rx / r;
          fy -= mj * ry / r;
        }
        pv[i].fx += fx;
        pv[i].fy += fy;
        fx = sqrt(fx*fx + fy*fy)/rmin;
        if (fx > max_f) max_f = fx;
      }
  return max_f;
}


double NBody::ParallelComputeNewPos(Particle particles[], ParticleV pv[], int npart, double max_f)
{
  int i;
  double a0, a1, a2;
  static double dt_old = 0.001, dt = 0.001;
  double dt_new;
  a0	 = 2.0 / (dt * (dt + dt_old));
  a2	 = 2.0 / (dt_old * (dt + dt_old));
  a1	 = -(a0 + a2);

  #pragma omp parallel for num_threads(num_threads)
  for (i=0; i<npart; i++) {
    double xi, yi;
    xi	           = particles[i].x;
    yi	           = particles[i].y;
    particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
    particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
    pv[i].xold     = xi;
    pv[i].yold     = yi;
    pv[i].fx       = 0;
    pv[i].fy       = 0;
  }

  dt_new = 1.0/sqrt(max_f);
  /* Set a minimum: */
  if (dt_new < 1.0e-6) dt_new = 1.0e-6;
  /* Modify time step */
  if (dt_new < dt) {
    dt_old = dt;
    dt     = dt_new;
  }
  else if (dt_new > 4.0 * dt) {
    dt_old = dt;
    dt    *= 2.0;
  }
  return dt_old;
}


//utils

void NBody::setNumThreads(int n) {
    num_threads = n;
}


