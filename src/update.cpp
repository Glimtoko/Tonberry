#include "update.hpp"
#include "flux.hpp"

#include <math.h>
#include <iostream>

double getLimiter(double di1, double di2, double omega) {
    double xi;
    // Slope limiter - Van Leer
    if (di2 == 0) {
        xi = 0.0;
    } else {
        double r = di1/di2;
        if (r <= 0.0) {
            xi = 0.0;
        } else {
            double xiR = 2.0/(1.0 - omega + (1 + omega)*r);
            xi = std::min(2*r/(1+r), xiR);
        }
    }
    return xi;
}

double getSlope(QUANT_2D U, int i, int j, int ud, int lr, double omega) {
    double di1 = U[j][i] - U[j-ud][i-lr];
    double di2 = U[j+ud][i+lr] - U[j][i];

    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    diU *= getLimiter(di1, di2, omega);

    return diU;
}


inline double calc_flux_rho(double rho, double u) {
    return rho*u;
}


inline double calc_flux_mom(double rho, double u, double v, double p) {
    return rho*u*v + p;
}


inline double calc_flux_E(double u, double E, double p) {
    return u*(E + p);
}


double getTimestep(Mesh2D mesh) {
    double Sx = 0.0, Sy = 0.0;
    for (int j=0; j<mesh.njGhosts-1; j++) {
        for (int i=0; i<mesh.niGhosts-1; i++) {
            double u = mesh.momU[j][i]/mesh.rho[j][i];
            double v = mesh.momV[j][i]/mesh.rho[j][i];
            double p = (mesh.gamma - 1.0)*(
                mesh.E[j][i] - 0.5*mesh.rho[j][i]*u*u
                             - 0.5*mesh.rho[j][i]*v*v
            );
            double a = sqrt((mesh.gamma*p)/mesh.rho[j][i]);
            Sx = std::max(Sx, a + u);
            Sy = std::max(Sy, a + v);
        }
    }
    return std::min(mesh.dtmax, std::min(mesh.cfl*mesh.dx/Sx, mesh.cfl*mesh.dy/Sy));
}


// momU is always momentum in sweep direction, momV is tangential
void sweep(
    QUANT_2D &rho, QUANT_2D &momU, QUANT_2D &momV, QUANT_2D &E,
    int ni, int nj, int iUpper, int jUpper,
    Sweep sweep,
    double dx, double omega, double dt, double gamma)
{

    int di, dj;
    switch (sweep) {
        case x:
            di = 1;
            dj = 0;
            break;
        case y:
            di = 0;
            dj = 1;
            break;
    }

    Mesh2D meshLD = Mesh2D(ni, nj, 0.0);
    Mesh2D meshRU = Mesh2D(ni, nj, 0.0);

    std::vector<std::vector<Flux::Flux>> flux;
    flux.resize(nj, std::vector<Flux::Flux>(ni));

    // Data reconstruction
    for (int j=1; j<jUpper+1; j++) {
        for (int i=1; i<iUpper+1; i++) {
            double dRho = 0.5*getSlope(rho, i, j, dj, di, omega);
            meshLD.rho[j][i] = rho[j][i] - dRho;
            meshRU.rho[j][i] = rho[j][i] + dRho;

            double dMomU = 0.5*getSlope(momU, i, j, dj, di, omega);
            meshLD.momU[j][i] = momU[j][i] - dMomU;
            meshRU.momU[j][i] = momU[j][i] + dMomU;

            double dMomV = 0.5*getSlope(momV, i, j, dj, di, omega);
            meshLD.momV[j][i] = momV[j][i] - dMomV;
            meshRU.momV[j][i] = momV[j][i] + dMomV;

            double dE = 0.5*getSlope(E, i, j, dj, di, omega);
            meshLD.E[j][i] = E[j][i] - dE;
            meshRU.E[j][i] = E[j][i] + dE;
        }
    }


    // Data evolution to half timestep
    double f = 0.5*dt/dx;
    for (int j=1; j<jUpper+1; j++) {
        for (int i=1; i<iUpper+1; i++) {
            double uLD = meshLD.momU[j][i]/meshLD.rho[j][i];
            double uRU = meshRU.momU[j][i]/meshRU.rho[j][i];

            double vLD = meshLD.momV[j][i]/meshLD.rho[j][i];
            double vRU = meshRU.momV[j][i]/meshRU.rho[j][i];

            double pLD = (gamma - 1.0)*(
                meshLD.E[j][i] - 0.5*meshLD.rho[j][i]*uLD*uLD
                            - 0.5*meshLD.rho[j][i]*vLD*vLD
            );

            double pRU = (gamma - 1.0)*(
                meshRU.E[j][i] - 0.5*meshRU.rho[j][i]*uRU*uRU
                            - 0.5*meshRU.rho[j][i]*vRU*vRU
            );

            double dF_rho = f*(calc_flux_rho(meshLD.rho[j][i], uLD) -
                               calc_flux_rho(meshRU.rho[j][i], uRU));

            double dF_momU = f*(calc_flux_mom(meshLD.rho[j][i], uLD, uLD, pLD) -
                               calc_flux_mom(meshRU.rho[j][i], uRU, uRU, pRU));

            double dF_momV = f*(calc_flux_mom(meshLD.rho[j][i], uLD, vLD, 0.0) -
                               calc_flux_mom(meshRU.rho[j][i], uRU, vRU, 0.0));

            double dF_E = f*(calc_flux_E(uLD, meshLD.E[j][i], pLD) -
                             calc_flux_E(uRU, meshRU.E[j][i], pRU));

            meshLD.rho[j][i] += dF_rho;
            meshRU.rho[j][i] += dF_rho;
            meshLD.momU[j][i] += dF_momU;
            meshRU.momU[j][i] += dF_momU;
            meshLD.momV[j][i] += dF_momV;
            meshRU.momV[j][i] += dF_momV;
            meshLD.E[j][i] += dF_E;
            meshRU.E[j][i] += dF_E;
        }
    }

    // Get intercell fluxes from Riemann solver
    for (int j=1; j<jUpper+1; j++) {
        for (int i=1; i<iUpper+1; i++) {
            double rhoLD_cell = meshRU.rho[j][i];
            double uLD_cell = meshRU.momU[j][i]/meshRU.rho[j][i];
            double vLD_cell = meshRU.momV[j][i]/meshRU.rho[j][i];
            double pLD_cell = (gamma - 1.0)*(
                meshRU.E[j][i] - 0.5*meshRU.rho[j][i]*uLD_cell*uLD_cell
                               - 0.5*meshRU.rho[j][i]*vLD_cell*vLD_cell
            );

            double rhoRU_cell = meshRU.rho[j+dj][i+di];
            double uRU_cell = meshRU.momU[j+dj][i+di]/meshRU.rho[j+dj][i+di];
            double vRU_cell = meshRU.momV[j+dj][i+di]/meshRU.rho[j+dj][i+di];
            double pRU_cell = (gamma - 1.0)*(
                meshLD.E[j+dj][i+di] - 0.5*meshLD.rho[j+dj][i+di]*uRU_cell*uRU_cell
                                     - 0.5*meshLD.rho[j+dj][i+di]*vRU_cell*vRU_cell
            );

            flux[j][i] = Flux::getFluxHLLC(
                uLD_cell, vLD_cell, rhoLD_cell, pLD_cell,
                uRU_cell, vRU_cell, rhoRU_cell, pRU_cell,
                gamma);

        }
    }

    // Update
    f = dt/dx;
    for (int j=2; j<jUpper; j++) {
        for (int i=2; i<iUpper; i++) {
            rho[j][i] += f*(flux[j-dj][i-di].rho - flux[j][i].rho);
            momU[j][i] += f*(flux[j-dj][i-di].momU - flux[j][i].momU);
            momV[j][i] += f*(flux[j-dj][i-di].momV - flux[j][i].momV);
            E[j][i] += f*(flux[j-dj][i-di].E - flux[j][i].E);
        }
    }
}
