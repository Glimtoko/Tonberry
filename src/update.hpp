#ifndef UPDATENEWH
#define UPDATENEWH
#include "mesh/mesh2d.hpp"
#include "sweep.hpp"

// double getTimestep(Mesh2D);
void sweepX(Mesh2D &mesh, double dt);
void sweepY(Mesh2D &mesh, double dt);

#endif
