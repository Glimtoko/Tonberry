#ifndef UPDATEH
#define UPDATEH
#include "mesh2d.hpp"

enum Sweep { x, y };

double getTimestep(Mesh2D);
void sweep(
    QUANT_2D &,  QUANT_2D &,  QUANT_2D &,  QUANT_2D &,
    int, int, int, int,
    Sweep,
    double, double, double, double);
#endif
