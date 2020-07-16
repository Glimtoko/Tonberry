#include "update.hpp"
#include "hydro/hydro.hpp"

#include <math.h>
#include <iostream>
#include <cstdlib>

#include "mpi.h"



// double getTimestep(Mesh2D mesh) {
//     double Sx = 0.0, Sy = 0.0;
//     for (int j=0; j<mesh.njGhosts-1; j++) {
//         for (int i=0; i<mesh.niGhosts-1; i++) {
//             double u = GET(mesh, momU, j, i)/GET(mesh, rho, j, i);
//             double v = GET(mesh, momV, j, i)/GET(mesh, rho, j, i);
//             double p = (mesh.gamma - 1.0)*(
//                 GET(mesh, E, j, i) - 0.5*GET(mesh, rho, j, i)*u*u
//                                    - 0.5*GET(mesh, rho, j, i)*v*v
//             );
//             double a = sqrt((mesh.gamma*p)/GET(mesh, rho, j, i));
//             Sx = std::max(Sx, a + abs(u));
//             Sy = std::max(Sy, a + abs(v));
//         }
//     }
//     return std::min(mesh.dtmax, std::min(mesh.cfl*mesh.dx/Sx, mesh.cfl*mesh.dy/Sy));
// }


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
    double* rhoAll = new double[ni*nj];
    double* momNAll = new double[ni*nj];
    double* momTAll = new double[ni*nj];
    double* EAll = new double[ni*nj];

    int index = 0;
    for (int j=2; j<jUpper; j++) {
        for (int i=0; i<ni; i++) {
            rhoAll[index] = GET(mesh, rho, j, i);
            momNAll[index] = GET(mesh, momU, j, i);
            momTAll[index] = GET(mesh, momV, j, i);
            EAll[index] = GET(mesh, E, j, i);

            index += 1;
        }
    }

    // Local 1D representation of data
    double* rhoLocal = new double[ni*njDecomp];
    double* momNLocal = new double[ni*njDecomp];
    double* momTLocal = new double[ni*njDecomp];
    double* ELocal = new double[ni*njDecomp];

    MPI::COMM_WORLD.Scatter(rhoAll, ni*njDecomp, MPI::DOUBLE, rhoLocal, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momNAll, ni*njDecomp, MPI::DOUBLE, momNLocal, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momTAll, ni*njDecomp, MPI::DOUBLE, momTLocal, ni*njDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(EAll, ni*njDecomp, MPI::DOUBLE, ELocal, ni*njDecomp, MPI::DOUBLE, 0);

    double* rho = new double[ni];
    double* momN = new double[ni];
    double* momT = new double[ni];
    double* E = new double[ni];

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
        Hydro::MUSCLHancock1D(rho, E, momN, momT, ni, iUpper, gamma, dt, dx);

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
            GET(mesh, rho, j, i) = rhoAll[index];
            GET(mesh, momU, j, i) = momNAll[index];
            GET(mesh, momV, j, i) = momTAll[index];
            GET(mesh, E, j, i) = EAll[index];
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
    double* rhoAll = new double[ni*nj];
    double* momNAll = new double[ni*nj];
    double* momTAll = new double[ni*nj];
    double* EAll = new double[ni*nj];

    int index = 0;
    for (int i=2; i<iUpper; i++) {
        for (int j=0; j<nj; j++) {
            rhoAll[index] = GET(mesh, rho, j, i);
            momNAll[index] = GET(mesh, momV, j, i);
            momTAll[index] = GET(mesh, momU, j, i);
            EAll[index] = GET(mesh, E, j, i);

            index += 1;
        }
    }

    // Local 1D representation of data
    double* rhoLocal = new double[nj*niDecomp];
    double* momNLocal = new double[nj*niDecomp];
    double* momTLocal = new double[nj*niDecomp];
    double* ELocal = new double[nj*niDecomp];

    MPI::COMM_WORLD.Scatter(rhoAll, nj*niDecomp, MPI::DOUBLE, rhoLocal, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momNAll, nj*niDecomp, MPI::DOUBLE, momNLocal, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(momTAll, nj*niDecomp, MPI::DOUBLE, momTLocal, nj*niDecomp, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Scatter(EAll, nj*niDecomp, MPI::DOUBLE, ELocal, nj*niDecomp, MPI::DOUBLE, 0);

    double* rho = new double[nj];
    double* momN = new double[nj];
    double* momT = new double[nj];
    double* E = new double[nj];

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
        Hydro::MUSCLHancock1D(rho, E, momN, momT, nj, jUpper, gamma, dt, dy);

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
            GET(mesh, rho, j, i) = rhoAll[index];
            GET(mesh, momU, j, i) = momTAll[index];
            GET(mesh, momV, j, i) = momNAll[index];
            GET(mesh, E, j, i) = EAll[index];
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
