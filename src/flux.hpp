#ifndef FLUXH
#define FLUXH

#include "mesh2d.hpp"

namespace Flux {
    struct Flux {
        double rho;
        double momU;
        double momV;
        double E;
    };

    Flux getFluxHLLC(
    double uL, double vL, double rhoL, double pL,
    double uR, double vR, double rhoR, double pR,
    double gamma);
}

#endif
