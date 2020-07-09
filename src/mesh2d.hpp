#ifndef MESH2D_H
#define MESH2D_H
#include <vector>

typedef std::vector<double> QUANT_1D;
typedef std::vector<std::vector<double> > QUANT_2D;

class Mesh2D {
public:
    const int nghosts = 2;
    int niGhosts;
    int njGhosts;
    int iUpper;
    int jUpper;

    double gamma;
    double dtmax;
    double dx, dy;
    double cfl;

    QUANT_1D x;
    QUANT_1D y;

    QUANT_2D rho;
    QUANT_2D momU;
    QUANT_2D momV;
    QUANT_2D E;

    // Constructor to just produce a Sod mesh
    Mesh2D(int, int, int);
    Mesh2D(int, int, double);
    void setBoundaries();
    void dumpToNetCDF();
    void dumpToNetCDF_NG();
};
#endif
