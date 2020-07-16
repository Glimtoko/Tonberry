#include "hydro_flux.hpp"
#include <math.h>

double Hydro::getTimestep(Mesh2D mesh) {
    double Sx = 0.0, Sy = 0.0;
    for (int j=0; j<mesh.njGhosts-1; j++) {
        for (int i=0; i<mesh.niGhosts-1; i++) {
            double u = GET(mesh, momU, j, i)/GET(mesh, rho, j, i);
            double v = GET(mesh, momV, j, i)/GET(mesh, rho, j, i);
            double p = (mesh.gamma - 1.0)*(
                GET(mesh, E, j, i) - 0.5*GET(mesh, rho, j, i)*u*u
                                   - 0.5*GET(mesh, rho, j, i)*v*v
            );
            double a = sqrt((mesh.gamma*p)/GET(mesh, rho, j, i));
            Sx = std::max(Sx, a + abs(u));
            Sy = std::max(Sy, a + abs(v));
        }
    }
    return std::min(mesh.dtmax, std::min(mesh.cfl*mesh.dx/Sx, mesh.cfl*mesh.dy/Sy));
}

