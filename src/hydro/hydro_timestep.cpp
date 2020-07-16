#include "hydro.hpp"
#include <math.h>

// double Hydro::getTimestep(Mesh2D mesh) {
//     double Sx = 0.0, Sy = 0.0;
//     for (int j=0; j<mesh.njGhosts-1; j++) {
//         for (int i=0; i<mesh.niGhosts-1; i++) {
//             double u = GET(mesh, momU, j, i)/GET(mesh, rho, j, i);
//             double v = GET(mesh, momV, j, i)/GET(mesh, rho, j, i);
//             double p = (gamma - 1.0)*(
//                 GET(mesh, E, j, i) - 0.5*GET(mesh, rho, j, i)*u*u
//                                    - 0.5*GET(mesh, rho, j, i)*v*v
//             );
//             double a = sqrt((gamma*p)/GET(mesh, rho, j, i));
//             Sx = std::max(Sx, a + abs(u));
//             Sy = std::max(Sy, a + abs(v));
//         }
//     }
//     return std::min(mesh.dtmax, std::min(mesh.cfl*mesh.dx/Sx, mesh.cfl*mesh.dy/Sy));
// }

double Hydro::getTimestep(
    double *rho, double *momU, double *momV, double *E,
    double dx, double dy,
    int nCells,
    double gamma, double cfl, double dtmax
) {
    double Sx = 0.0, Sy = 0.0;
    for (int i=0; i<nCells; i++) {
        double u = momU[i]/rho[i];
        double v = momV[i]/rho[i];
        double p = (gamma - 1.0)*(E[i] - 0.5*rho[i]*u*u - 0.5*rho[i]*u*u);
        double a = sqrt((gamma*p)/rho[i]);

        Sx = std::max(Sx, a + abs(u));
        Sy = std::max(Sy, a + abs(v));
    }

    return std::min(dtmax, std::min(cfl*dx/Sx, cfl*dy/Sy));
}

