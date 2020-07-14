#include "mesh2d.hpp"
#include "update.hpp"

#include <iostream>
#include "mpi.h"

#include <fenv.h>

int main(int argc, char* argv[]) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);

    // Mesh sizes - hardcoded for now
    const int ni = 250;
    const int nj = 250;
    const int problem = 2;

    const double dtOut = 0.05;
    const double tEnd = 2.0;

    // Initialise MPI
    MPI::Init(argc, argv);

    // MPI environment
    int nprocs, myrank, error;
    nprocs = MPI::COMM_WORLD.Get_size();
    myrank = MPI::COMM_WORLD.Get_rank();

    // 100x100 mesh using spherical Sod set-up
    Mesh2D mesh(ni, nj, problem);

    if (myrank == 0) {
        mesh.dumpToSILO(0.0);
    }

    double t = 0.0;
    double outNext = t + dtOut;
    int step = 0;

    for( ; ; ) {
        step++;

        double dt = getTimestep(mesh);
        dt = std::min(dt, outNext - t);

        MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);

        if (myrank == 0) {
            std::cout << "Step: " << step;
            std::cout << ", time = " << t;
            std::cout << ", dt = " << dt << std::endl;
        }

        // Sweep is half X, then Y, then half X
        sweepX(mesh, dt/2.0);
        sweepY(mesh, dt);
        sweepX(mesh, dt/2.0);
//         mesh.setBoundaries();

        t += dt;
        if (t >= outNext) {
            outNext += dtOut;
            mesh.dumpToSILO(t);
        }
        if (t > tEnd || step > 10000) break;
    }

    std::cout << step << std::endl;
    if (myrank == 0) {
//         mesh.dumpToNetCDF();
//         mesh.dumpToNetCDF_NG();
//         mesh.dumpToSILO(5.0);
    }

//     mesh.Kill();

//     for (int i=2; i<mesh.jUpper; i++) {
//         std::cout << i << " " << mesh.y[i] << " " << mesh.rho[i][2] << " " << mesh.momV[i][2] << std::endl;
//     }

    MPI::Finalize();
    return 0;
}
