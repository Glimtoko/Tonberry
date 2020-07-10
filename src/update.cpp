#include "update.hpp"
#include "flux.hpp"

#include <math.h>
#include <iostream>

double getLimiter1D(double di1, double di2, double omega) {
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

double getSlope1D(QUANT_1D U, int i, double omega) {
    double di1 = U[i] - U[i-1];
    double di2 = U[i+1] - U[i];

    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    diU *= getLimiter1D(di1, di2, omega);

    return diU;
}


inline double calcFluxRho(double rho, double u) {
    return rho*u;
}


inline double calcFluxMom(double rho, double u, double v, double p) {
    return rho*u*v + p;
}


inline double calcFluxE(double u, double E, double p) {
    return u*(E + p);
}


void sweep1D(
    QUANT_1D &rho, QUANT_1D &E, QUANT_1D &momN, QUANT_1D &momT,
    int ni, int iUpper, double gamma, double dt, double dx)
{
    double omega = 0.0;

    QUANT_1D rhoL;
    QUANT_1D rhoR;
    QUANT_1D momTL;
    QUANT_1D momTR;
    QUANT_1D momNL;
    QUANT_1D momNR;
    QUANT_1D EL;
    QUANT_1D ER;

    rhoL.assign(ni, 0.0);
    rhoR.assign(ni, 0.0);
    momTL.assign(ni, 0.0);
    momTR.assign(ni, 0.0);
    momNL.assign(ni, 0.0);
    momNR.assign(ni, 0.0);
    EL.assign(ni, 0.0);
    ER.assign(ni, 0.0);

    std::vector<Flux::Flux> flux;
    flux.resize(ni);

    // Data reconstruction
    for (int i=1; i<iUpper+1; i++) {
        double dRho = 0.5*getSlope1D(rho, i, omega);
        rhoL[i] = rho[i] - dRho;
        rhoR[i] = rho[i] + dRho;

        double dMomU = 0.5*getSlope1D(momN, i, omega);
        momNL[i] = momN[i] - dMomU;
        momNR[i] = momN[i] + dMomU;

        double dMomV = 0.5*getSlope1D(momT, i, omega);
        momTL[i] = momT[i] - dMomV;
        momTR[i] = momT[i] + dMomV;

        double dE = 0.5*getSlope1D(E, i, omega);
        EL[i] = E[i] - dE;
        ER[i] = E[i] + dE;
    }



    // Data evolution to half timestep
    double f = 0.5*dt/dx;
    for (int i=1; i<iUpper+1; i++) {
        double uLD = momNL[i]/rhoL[i];
        double uRU = momNR[i]/rhoR[i];

        double vLD = momTL[i]/rhoL[i];
        double vRU = momTR[i]/rhoR[i];

        double pLD = (gamma - 1.0)*(
            EL[i] - 0.5*rhoL[i]*uLD*uLD
                        - 0.5*rhoL[i]*vLD*vLD
        );

        double pRU = (gamma - 1.0)*(
            ER[i] - 0.5*rhoR[i]*uRU*uRU
                        - 0.5*rhoR[i]*vRU*vRU
        );

        double dF_rho = f*(calcFluxRho(rhoL[i], uLD) -
                            calcFluxRho(rhoR[i], uRU));

        double dF_momN = f*(calcFluxMom(rhoL[i], uLD, uLD, pLD) -
                            calcFluxMom(rhoR[i], uRU, uRU, pRU));

        double dF_momT = f*(calcFluxMom(rhoL[i], uLD, vLD, 0.0) -
                            calcFluxMom(rhoR[i], uRU, vRU, 0.0));

        double dF_E = f*(calcFluxE(uLD, EL[i], pLD) -
                            calcFluxE(uRU, ER[i], pRU));

        rhoL[i] += dF_rho;
        rhoR[i] += dF_rho;
        momNL[i] += dF_momN;
        momNR[i] += dF_momN;
        momTL[i] += dF_momT;
        momTR[i] += dF_momT;
        EL[i] += dF_E;
        ER[i] += dF_E;
    }

    // Get intercell fluxes from Riemann solver
    for (int i=1; i<iUpper+1; i++) {
        double rhoLD_cell = rhoR[i];
        double uLD_cell = momNR[i]/rhoR[i];
        double vLD_cell = momTR[i]/rhoR[i];
        double pLD_cell = (gamma - 1.0)*(
            ER[i] - 0.5*rhoR[i]*uLD_cell*uLD_cell
                            - 0.5*rhoR[i]*vLD_cell*vLD_cell
        );

        double rhoRU_cell = rhoR[i+1];
        double uRU_cell = momNR[i+1]/rhoR[i+1];
        double vRU_cell = momTR[i+1]/rhoR[i+1];
        double pRU_cell = (gamma - 1.0)*(
            EL[i+1] - 0.5*rhoL[i+1]*uRU_cell*uRU_cell
                                    - 0.5*rhoL[i+1]*vRU_cell*vRU_cell
        );

        flux[i] = Flux::getFluxHLLC(
            uLD_cell, vLD_cell, rhoLD_cell, pLD_cell,
            uRU_cell, vRU_cell, rhoRU_cell, pRU_cell,
            gamma);

    }

    // Update
    f = dt/dx;
    for (int i=2; i<iUpper; i++) {
        rho[i] += f*(flux[i-1].rho - flux[i].rho);
        momN[i] += f*(flux[i-1].momU - flux[i].momU);
        momT[i] += f*(flux[i-1].momV - flux[i].momV);
        E[i] += f*(flux[i-1].E - flux[i].E);
    }
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


void sweepX(Mesh2D &mesh, double dt) {
    int ni = mesh.niGhosts;
    int iUpper = mesh.iUpper;
    int jUpper = mesh.jUpper;
    double gamma = mesh.gamma;
    double dx = mesh.dx;

    QUANT_1D rho;
    QUANT_1D momN;
    QUANT_1D momT;
    QUANT_1D E;

    rho.assign(ni, 0.0);
    momN.assign(ni, 0.0);
    momT.assign(ni, 0.0);
    E.assign(ni, 0.0);

    for (int j=2; j<jUpper; j++) {
        // Pack 1D data
        for (int i=0; i<ni; i++) {
            rho[i] = mesh.rho[j][i];
            momN[i] = mesh.momU[j][i];
            momT[i] = mesh.momV[j][i];
            E[i] = mesh.E[j][i];
        }
        // Sweep
        sweep1D(rho, E, momN, momT, ni, iUpper, gamma, dt, dx);

        // Unpack 1D data
        for (int i=0; i<ni; i++) {
            mesh.rho[j][i] = rho[i];
            mesh.momU[j][i] = momN[i];
            mesh.momV[j][i] = momT[i];
            mesh.E[j][i] = E[i];
        }
    }

    // Boundary update
    mesh.setBoundaries();
}


void sweepY(Mesh2D &mesh, double dt) {
    int nj = mesh.njGhosts;
    int iUpper = mesh.iUpper;
    int jUpper = mesh.jUpper;
    double gamma = mesh.gamma;
    double dy = mesh.dy;

    QUANT_1D rho;
    QUANT_1D momN;
    QUANT_1D momT;
    QUANT_1D E;

    rho.assign(nj, 0.0);
    momN.assign(nj, 0.0);
    momT.assign(nj, 0.0);
    E.assign(nj, 0.0);

    for (int i=2; i<iUpper; i++) {
        // Pack 1D data
        for (int j=0; j<nj; j++) {
            rho[j] = mesh.rho[j][i];
            momN[j] = mesh.momV[j][i];
            momT[j] = mesh.momU[j][i];
            E[j] = mesh.E[j][i];
        }
        // Sweep
        sweep1D(rho, E, momN, momT, nj, jUpper, gamma, dt, dy);

        // Unpack 1D data
        for (int j=0; j<nj; j++) {
            mesh.rho[j][i] = rho[j];
            mesh.momV[j][i] = momN[j];
            mesh.momU[j][i] = momT[j];
            mesh.E[j][i] = E[j];
        }
    }

    // Boundary update
    mesh.setBoundaries();
}
