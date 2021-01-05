#ifndef MESH2D_H
#define MESH2D_H
// #include "sweep.hpp"/**/

#include "iris.hpp"
#include <vector>

#define GET(M, D, J, I) M.D[(J)*M.niGhosts + I]
#define _LGET(D, J, I) D[(J)*niGhosts + I]
#define _PGET(M, D, J, I) M->D[(J)*M->niGhosts + I]

class Mesh2D {
public:
    const int nghosts = 2;
    int niGhosts;
    int njGhosts;
    int iUpper;
    int jUpper;
    int nreg;

    int dumpStateNo = 0;

    double gamma;
    double dtmax;
    double dx, dy;
    double cfl;

    double* meshBoundaryUD;
    double* meshBoundaryLR;

    double* x;
    double* y;

    int * materials;
    int * regions;
    int * cell_map;

    double* rho;
    double* momU;
    double* momV;
    double* E;

    // Constructor to just produce a Sod mesh
    Mesh2D(int, int, int);

    void Kill();
    void setBoundaries();
    void dumpToSILO(double, int, iris::Logger*);

    void sweepX(double dt,
                void(*flux)(double*, double*, double*, double*, int, int, double, double, double)
    );

    void sweepY(double dt,
                void(*flux)(double*, double*, double*, double*, int, int, double, double, double)
    );

//     double* &getMomentum(Sweep, Direction);
};
#endif
