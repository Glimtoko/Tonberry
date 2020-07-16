#include "hydro.hpp"

#include <math.h>
#include <iostream>
#include <cstdlib>

double getLimiter1D(double di1, double di2, double omega);
double getSlope1D(double* U, int i, double omega);
double calcFluxRho(double rho, double u);
double calcFluxMom(double rho, double u, double v, double p);
double calcFluxE(double u, double E, double p);

void Hydro::MUSCLHancock1D(
    double* &rho, double* &E, double* &momN, double* &momT,
    int ni, int iUpper, double gamma, double dt, double dx)
{
    double omega = 0.0;

    double* rhoL = new double[ni];
    double* rhoR = new double[ni];
    double* momTL = new double[ni];
    double* momTR = new double[ni];
    double* momNL = new double[ni];
    double* momNR = new double[ni];
    double* EL = new double[ni];
    double* ER = new double[ni];

    std::vector<Hydro::Flux> flux;
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

        flux[i] = Hydro::getFluxHLLC(
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

double getSlope1D(double* U, int i, double omega) {
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
