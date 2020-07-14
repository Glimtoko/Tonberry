#include "flux.hpp"
#include <math.h>
#include <iostream>

Flux::Flux Flux::getFluxHLLC(
    double uL, double vL, double rhoL, double pL,
    double uR, double vR, double rhoR, double pR,
    double gamma) {

    // Pressure estimate from PVRS solver
    double aL = sqrt((gamma*pL)/rhoL);
    double aR = sqrt((gamma*pR)/rhoR);

    double rho_bar = 0.5*(rhoL + rhoR);
    double a_bar = 0.5*(aL + aR);

    double p_guess = 0.5*(pL + pR) - 0.5*(uR - uL)*(rho_bar*a_bar);


    double qL;
    if (p_guess <= pL) {
        qL = 1.0;
    } else {
        qL = sqrt(1.0 + (gamma + 1.0)/(2.0*gamma)*(p_guess/pL - 1.0));
    }

    double qR;
    if (p_guess <= pR) {
        qR = 1.0;
    } else {
        qR = sqrt(1.0 + (gamma + 1.0)/(2.0*gamma)*(p_guess/pR - 1.0));
    }

    // Estimate wave speeds
    double SL = uL - aL*qL;
    double SR = uR + aR*qR;
    double Sstar = pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR);
    Sstar /= (rhoL*(SL - uL) - rhoR*(SR - uR));

    // Get energy and momenta on boundaries
    double eL = pL/((gamma - 1.0)*rhoL);
    double EL = rhoL*(0.5*uL*uL + 0.5*vL*vL + eL);
    double momUL = uL*rhoL;
    double momVL = vL*rhoL;

    double eR = pR/((gamma - 1.0)*rhoR);
    double ER = rhoR*(0.5*uR*uR + 0.5*vR*vR + eR);
    double momUR = uR*rhoR;
    double momVR = vR*rhoR;

    // Funcion return
    Flux flux;

    // Get flux
    if (0 <= SL) {
        // Flux = F(UL)
        flux.rho = rhoL*uL;
        flux.momU = rhoL*uL*uL + pL;
        flux.momV = rhoL*vL*uL;
        flux.E = uL*(EL + pL);
    } else if (SL < 0 && 0 <= Sstar) {
        // Flux = F(*L)
        double pLR = 0.5*(pL + pR + rhoL*(SL - uL)*(Sstar - uL) + rhoR*(SR - uR)*(Sstar - uR));
        double d = SL - Sstar;

        // Set F(UL)
        double fL_rho = rhoL*uL;
        double fL_momU = rhoL*uL*uL + pL;
        double fL_momV = rhoL*uL*vL;
        double fL_E = uL*(EL + pL);

        // Find F(*L)
        flux.rho = (Sstar*(SL*rhoL - fL_rho))/d;
        flux.momU = (Sstar*(SL*momUL - fL_momU) + SL*pLR)/d;
        flux.momV = Sstar*(SL*momVL - fL_momV)/d;
        flux.E = (Sstar*(SL*EL - fL_E) + SL*pLR*Sstar)/d;
    } else if (Sstar < 0 && 0 <= SR) {
        // Flux = F(*R)
        double pLR = 0.5*(pL + pR + rhoL*(SL - uL)*(Sstar - uL) + rhoR*(SR - uR)*(Sstar - uR));
        double d = SR - Sstar;

        // Set F(UR)
        double fR_rho = rhoR*uR;
        double fR_momU = rhoR*uR*uR + pR;
        double fR_momV = rhoR*uR*vR;
        double fR_E = uR*(ER + pR);

        // Find F(*R)
        flux.rho = (Sstar*(SR*rhoR - fR_rho))/d;
        flux.momU = (Sstar*(SR*momUR - fR_momU) + SR*pLR)/d;
        flux.momV = Sstar*(SR*momVR - fR_momV)/d;
        flux.E = (Sstar*(SR*ER - fR_E) + SR*pLR*Sstar)/d;
    } else {
        // Flux = F(UR)
        flux.rho = rhoR*uR;
        flux.momU = rhoR*uR*uR + pR;
        flux.momV = rhoR*vR*uR;
        flux.E = uR*(ER + pR);
    }

    return flux;
}
