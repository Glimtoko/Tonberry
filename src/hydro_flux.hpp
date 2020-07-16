#ifndef FLUXH
#define FLUXH

#include "mesh2d.hpp"

namespace Hydro {
    struct Flux {
        double rho;
        double momU;
        double momV;
        double E;
    };

    Flux getFluxHLLC(
        double uL, double vL, double rhoL, double pL,
        double uR, double vR, double rhoR, double pR,
        double gamma
    );

    void MUSCLHancock1D(
        double* &rho, double* &E, double* &momN, double* &momT,
        int ni, int iUpper, double gamma, double dt, double dx
    );

    double getTimestep(Mesh2D);
}

#endif
