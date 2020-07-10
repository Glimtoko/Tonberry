#ifndef UPDATENEWH
#define UPDATENEWH
#include "mesh2d.hpp"
#include "sweep.hpp"

void sweep1D(
    QUANT_1D &rho, QUANT_1D &E, QUANT_1D &momN, QUANT_1D &momT,
    int ni, int iUpper, double gamma, double dt, double dx);

double getTimestep(Mesh2D);
void sweepX(Mesh2D &mesh, double dt);
void sweepY(Mesh2D &mesh, double dt);

#endif
