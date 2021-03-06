#include "mesh/mesh2d.hpp"
#include "hydro/hydro.hpp"

#include <iostream>
#include "mpi.h"

#include <fenv.h>

int main(int argc, char* argv[]) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);

    // Mesh sizes - hardcoded for now
    const int ni = 600;
    const int nj = 200;
    const int problem = 4;

    const double dtOut = 0.2;
    const double tEnd = 80.0;

    // Initialise MPI
    MPI::Init(argc, argv);

    // MPI environment
    int nprocs, myrank, error;
    nprocs = MPI::COMM_WORLD.Get_size();
    myrank = MPI::COMM_WORLD.Get_rank();

    // 100x100 mesh using spherical Sod set-up
    Mesh2D mesh(ni, nj, problem);

    if (myrank == 0) {
        mesh.dumpToSILO(0.0, 0);
    }

    double t = 0.0;
    double outNext = t + dtOut;
    int step = 0;

    for( ; ; ) {
        step++;

        double dt = Hydro::getTimestep(
            mesh.rho,  mesh.momU,  mesh.momV,  mesh.E,
            mesh.dx, mesh.dy,
            mesh.niGhosts*mesh.njGhosts,
            mesh.gamma, mesh.cfl, mesh.dtmax
        );

        dt = std::min(dt, outNext - t);

        MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);

        if (myrank == 0) {
            std::cout << "Step: " << step;
            std::cout << ", time = " << t;
            std::cout << ", dt = " << dt << std::endl;
        }

        // Sweep is half X, then Y, then half X
        mesh.sweepX(dt/2.0, Hydro::MUSCLHancock1D);
        mesh.sweepY(dt, Hydro::MUSCLHancock1D);
        mesh.sweepX(dt/2.0, Hydro::MUSCLHancock1D);

        t += dt;
        if (t >= outNext) {
            outNext += dtOut;
            if (myrank == 0) mesh.dumpToSILO(t, step);
        }
        if (t > tEnd || step > 10000) break;
    }

    mesh.Kill();

    MPI::Finalize();
    return 0;
}
