#ifndef UPDATEH
#define UPDATEH
#include "mesh2d.hpp"
#include "sweep.hpp"

// enum Sweep { x, y };

double getTimestep(Mesh2D);
void sweep(
//     QUANT_2D &,  QUANT_2D &,  QUANT_2D &,  QUANT_2D &,
    Mesh2D &,
//     int, int, int, int,
    Sweep,
    double);
#endif
