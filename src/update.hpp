#ifndef UPDATENEWH
#define UPDATENEWH
#include "mesh2d.hpp"
#include "sweep.hpp"

// void MUSCLHancock1D(
//     double* &rho, double* &E, double* &momN, double* &momT,
//     int ni, int iUpper, double gamma, double dt, double dx);

// double getTimestep(Mesh2D);
void sweepX(Mesh2D &mesh, double dt);
void sweepY(Mesh2D &mesh, double dt);

#endif
