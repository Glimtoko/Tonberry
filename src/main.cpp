#include "mesh/mesh2d.hpp"
#include "hydro/hydro.hpp"

#include <iostream>
#include <memory>

#include "mpi.h"

#include "iris.hpp"
#include "arcus.hpp"

#include <fenv.h>

int main(int argc, char* argv[]) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);

    // Mesh sizes - hardcoded for now
    const int ni = 60;
    const int nj = 20;
    const int problem = 4;

    const double dtOut = 1.0;
    const double tEnd = 10.0;  // N.b. this *should* be 80+

    // Initialise MPI
    MPI::Init(argc, argv);

    // MPI environment
    int nprocs, myrank, error;
    nprocs = MPI::COMM_WORLD.Get_size();
    myrank = MPI::COMM_WORLD.Get_rank();

    // 100x100 mesh using spherical Sod set-up
    Mesh2D mesh(ni, nj, problem);

    // Set logger registry
    const int niG = ni + 4;
    const int njG = nj + 4;
    std::shared_ptr<iris::arcus::LogRegistry> log_reg(
        iris::arcus::LogRegistry::Create(
        mesh.regions, mesh.materials, ni*nj, ni*nj
        )
    );

    double t = 0.0;
    double outNext = t + dtOut;
    int step = 0;
    double dt;

    log_reg->AddQuant("density", mesh.rho, iris::arcus::Centring::cell);
    log_reg->AddQuant("momv", mesh.momV, iris::arcus::Centring::cell);
    log_reg->AddQuant("momu", mesh.momU, iris::arcus::Centring::cell);
    log_reg->AddQuant("e", mesh.E, iris::arcus::Centring::cell);

    log_reg->AddScalarDouble("t", &t);
    log_reg->AddScalarDouble("dt", &dt);
    log_reg->AddScalarInt("step", &step);

    std::shared_ptr<iris::Logger> coreLogger(
        new iris::Logger("BaseLogger", iris::Levels::Info)
    );

    coreLogger->addConsoleAppender("::CORE:: %m%n");
    coreLogger->setRegistry(log_reg);

    std::shared_ptr<iris::Logger> spLogger(
        new iris::Logger("SPrintLogger", iris::Levels::Info)
    );

    spLogger->addConsoleAppender("::SPRINT:: %m%n");
    spLogger->addFileAppender("::SPRINT:: %m%n", "sprint.dat");
    spLogger->setRegistry(log_reg);

    if (myrank == 0) {
        mesh.dumpToSILO(0.0, 0, coreLogger);
    }

    for( ; ; ) {
        step++;

        dt = Hydro::getTimestep(
            mesh.rho,  mesh.momU,  mesh.momV,  mesh.E,
            mesh.dx, mesh.dy,
            mesh.niGhosts*mesh.njGhosts,
            mesh.gamma, mesh.cfl, mesh.dtmax
        );

        dt = std::min(dt, outNext - t);

        MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);

        if (myrank == 0) 
            coreLogger->logvar(iris::Levels::Info, "Step: @{step} :: Time = @{t}, dt = @{dt}");

        // Sweep is half X, then Y, then half X
        mesh.sweepX(dt/2.0, Hydro::MUSCLHancock1D);
        mesh.sweepY(dt, Hydro::MUSCLHancock1D);
        mesh.sweepX(dt/2.0, Hydro::MUSCLHancock1D);

        t += dt;
        if (t >= outNext) {
            outNext += dtOut;
            if (myrank == 0) {
                mesh.dumpToSILO(t, step, coreLogger);

                spLogger->logvar(iris::Levels::Info, "");
                spLogger->logvar(iris::Levels::Info, "Short Print at t = @{t}");
                spLogger->logvar(iris::Levels::Info, "Region  Rho          momU         momV         E");
                for (int r=0; r<mesh.nreg; r++) {
                    spLogger->logvar(
                        iris::Levels::Info,
                        "@{regno}       @{density.region.n.sum}  @{momu.region.n.sum}  @{momv.region.n.sum}  @{e.region.n.sum}",
                        r+1
                    );
                }
                spLogger->logvar(iris::Levels::Info, "");
                spLogger->logvar(
                    iris::Levels::Info,
                    "Total   @{density.sum}  @{momu.sum}  @{momv.sum}  @{e.sum}"
                );
                spLogger->logvar(iris::Levels::Info, "");
            }
        }
        if (t > tEnd || step > 10000) break;
    }

    mesh.Kill();

    MPI::Finalize();
    return 0;
}
