#include "update.hpp"
#include "flux.hpp"

#include <math.h>
#include <iostream>
#include <cstdlib>

#include "mpi.h"

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

    QUANT_1D rhoL = new double[ni];
    QUANT_1D rhoR = new double[ni];
    QUANT_1D momTL = new double[ni];
    QUANT_1D momTR = new double[ni];
    QUANT_1D momNL = new double[ni];
    QUANT_1D momNR = new double[ni];
    QUANT_1D EL = new double[ni];
    QUANT_1D ER = new double[ni];

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

        double rhoRU_cell = rhoL[i+1];
        double uRU_cell = momNL[i+1]/rhoL[i+1];
        double vRU_cell = momTL[i+1]/rhoL[i+1];
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

    delete[] rhoL;
    delete[] rhoR;
    delete[] momTL;
    delete[] momTR;
    delete[] momNL;
    delete[] momNR;
    delete[] EL;
    delete[] ER;
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
            Sx = std::max(Sx, a + abs(u));
            Sy = std::max(Sy, a + abs(v));
        }
    }
    return std::min(mesh.dtmax, std::min(mesh.cfl*mesh.dx/Sx, mesh.cfl*mesh.dy/Sy));
}


void sweepX(Mesh2D &mesh, double dt) {
    int ni = mesh.niGhosts;
    int nj = mesh.jUpper - 2;
    int iUpper = mesh.iUpper;
    int jUpper = mesh.jUpper;
    double gamma = mesh.gamma;
    double dx = mesh.dx;

    // MPI environment
    int nprocs, myrank, error;
    nprocs = MPI::COMM_WORLD.Get_size();
    myrank = MPI::COMM_WORLD.Get_rank();

    // Get decomposition size
    int njDecomp = nj/nprocs;

    // Create a full 1D representation of the data
    QUANT_1D rhoAll = new double[ni*nj];
    QUANT_1D momNAll = new double[ni*nj];
    QUANT_1D momTAll = new double[ni*nj];
    QUANT_1D EAll = new double[ni*nj];

    int index = 0;
    for (int j=2; j<jUpper; j++) {
        for (int i=0; i<ni; i++) {
            rhoAll[index] = mesh.rho[j][i];
            momNAll[index] = mesh.momU[j][i];
            momTAll[index] = mesh.momV[j][i];
            EAll[index] = mesh.E[j][i];

            index += 1;
        }
    }

    // Local 1D representation of data
    QUANT_1D rhoLocal = new double[ni*njDecomp];
    QUANT_1D momNLocal = new double[ni*njDecomp];
    QUANT_1D momTLocal = new double[ni*njDecomp];
    QUANT_1D ELocal = new double[ni*njDecomp];

    MPI::COMM_WORLD.Scatter(rhoAll, ni*njDecomp, MPI::DOUBLE, rhoLocal, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momNAll, ni*njDecomp, MPI::DOUBLE, momNLocal, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momTAll, ni*njDecomp, MPI::DOUBLE, momTLocal, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(EAll, ni*njDecomp, MPI::DOUBLE, ELocal, ni*njDecomp, MPI::DOUBLE, 0);

    QUANT_1D rho = new double[ni];
    QUANT_1D momN = new double[ni];
    QUANT_1D momT = new double[ni];
    QUANT_1D E = new double[ni];

    index = 0;
    int index2 = 0;
    for (int j=0; j<njDecomp; j++) {
        // Pack 1D data
        for (int i=0; i<ni; i++) {
            rho[i] = rhoLocal[index];
            momN[i] = momNLocal[index];
            momT[i] = momTLocal[index];
            E[i] = ELocal[index];
            index++;
        }
        // Sweep
        sweep1D(rho, E, momN, momT, ni, iUpper, gamma, dt, dx);

        // Unpack 1D data
        for (int i=0; i<ni; i++) {
            rhoLocal[index2] = rho[i];
            momNLocal[index2] = momN[i];
            momTLocal[index2] = momT[i];
            ELocal[index2] = E[i];
            index2++;
        }
    }

    MPI::COMM_WORLD.Gather(rhoLocal, ni*njDecomp, MPI::DOUBLE, rhoAll, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Gather(momNLocal, ni*njDecomp, MPI::DOUBLE, momNAll, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Gather(momTLocal, ni*njDecomp, MPI::DOUBLE, momTAll, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Gather(ELocal, ni*njDecomp, MPI::DOUBLE, EAll, ni*njDecomp, MPI::DOUBLE, 0);

    index = 0;
    for (int j=2; j<jUpper; j++) {
        for (int i=0; i<ni; i++) {
            mesh.rho[j][i] = rhoAll[index];
            mesh.momU[j][i] = momNAll[index];
            mesh.momV[j][i] = momTAll[index];
            mesh.E[j][i] = EAll[index];
            index += 1;
        }
    }

    delete[] rho;
    delete[] momN;
    delete[] momT;
    delete[] E;

    delete[] rhoAll;
    delete[] momNAll;
    delete[] momTAll;
    delete[] EAll;

    delete[] rhoLocal;
    delete[] momNLocal;
    delete[] momTLocal;
    delete[] ELocal;

    // Boundary update
    mesh.setBoundaries();
}


void sweepY(Mesh2D &mesh, double dt) {
    int nj = mesh.njGhosts;
    int ni = mesh.iUpper - 2;
    int iUpper = mesh.iUpper;
    int jUpper = mesh.jUpper;
    double gamma = mesh.gamma;
    double dy = mesh.dy;

    // MPI environment
    int nprocs, myrank, error;
    nprocs = MPI::COMM_WORLD.Get_size();
    myrank = MPI::COMM_WORLD.Get_rank();

    // Get decomposition size
    int niDecomp = ni/nprocs;

    // Create a full 1D representation of the data
    QUANT_1D rhoAll = new double[ni*nj];
    QUANT_1D momNAll = new double[ni*nj];
    QUANT_1D momTAll = new double[ni*nj];
    QUANT_1D EAll = new double[ni*nj];

    int index = 0;
    for (int i=2; i<iUpper; i++) {
        for (int j=0; j<nj; j++) {
            rhoAll[index] = mesh.rho[j][i];
            momNAll[index] = mesh.momV[j][i];
            momTAll[index] = mesh.momU[j][i];
            EAll[index] = mesh.E[j][i];

            index += 1;
        }
    }

    // Local 1D representation of data
    QUANT_1D rhoLocal = new double[nj*niDecomp];
    QUANT_1D momNLocal = new double[nj*niDecomp];
    QUANT_1D momTLocal = new double[nj*niDecomp];
    QUANT_1D ELocal = new double[nj*niDecomp];

    MPI::COMM_WORLD.Scatter(rhoAll, nj*niDecomp, MPI::DOUBLE, rhoLocal, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momNAll, nj*niDecomp, MPI::DOUBLE, momNLocal, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momTAll, nj*niDecomp, MPI::DOUBLE, momTLocal, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(EAll, nj*niDecomp, MPI::DOUBLE, ELocal, nj*niDecomp, MPI::DOUBLE, 0);

    QUANT_1D rho = new double[nj];
    QUANT_1D momN = new double[nj];
    QUANT_1D momT = new double[nj];
    QUANT_1D E = new double[nj];

    int index2 = 0;
    index = 0;
    for (int i=0; i<niDecomp; i++) {
        // Pack 1D data
        for (int j=0; j<nj; j++) {
            rho[j] = rhoLocal[index];
            momN[j] = momNLocal[index];
            momT[j] = momTLocal[index];
            E[j] = ELocal[index];
            index++;
        }
        // Sweep
        sweep1D(rho, E, momN, momT, nj, jUpper, gamma, dt, dy);

        // Unpack 1D data
        for (int j=0; j<nj; j++) {
            rhoLocal[index2] = rho[j];
            momNLocal[index2] = momN[j];
            momTLocal[index2] = momT[j];
            ELocal[index2] = E[j];
            index2++;
        }
    }

    MPI::COMM_WORLD.Gather(rhoLocal, nj*niDecomp, MPI::DOUBLE, rhoAll, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Gather(momNLocal, nj*niDecomp, MPI::DOUBLE, momNAll, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Gather(momTLocal, nj*niDecomp, MPI::DOUBLE, momTAll, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Gather(ELocal, nj*niDecomp, MPI::DOUBLE, EAll, nj*niDecomp, MPI::DOUBLE, 0);

    index = 0;
    for (int i=2; i<iUpper; i++) {
        for (int j=0; j<nj; j++) {
            mesh.rho[j][i] = rhoAll[index];
            mesh.momU[j][i] = momTAll[index];
            mesh.momV[j][i] = momNAll[index];
            mesh.E[j][i] = EAll[index];
            index += 1;
        }
    }

    delete[] rho;
    delete[] momN;
    delete[] momT;
    delete[] E;

    delete[] rhoAll;
    delete[] momNAll;
    delete[] momTAll;
    delete[] EAll;

    delete[] rhoLocal;
    delete[] momNLocal;
    delete[] momTLocal;
    delete[] ELocal;

    // Boundary update
    mesh.setBoundaries();
}
