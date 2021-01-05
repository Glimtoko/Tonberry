#include "mesh/mesh2d.hpp"
#include "hydro/hydro.hpp"

#include <iostream>
#include <memory>

#include "mpi.h"

#include "iris.hpp"
#include "arcus.hpp"
#include "arcus_c.h"
#include "iris_c.h"

#include <fenv.h>

int main(int argc, char* argv[]) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);

    // Mesh sizes - hardcoded for now
    const int ni = 60;
    const int nj = 20;
    const int problem = 4;

    const double dtOut = 1.0;
    const double tEnd = 1.0;  // N.b. this *should* be 80+

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


    log_reg->SetCellMap(mesh.cell_map);

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

    // ------------------------------------------------------------------------
    // Test of C interface
    arcus_reg_t c_reg = arcus_create_registry(
         mesh.regions, mesh.materials, ni*nj, ni*nj
    );

    int status;
    status = arcus_set_cell_map(c_reg, mesh.cell_map);
    
    char qname[] = "density";
    status = arcus_add_quant(c_reg, qname, mesh.rho, ARCUS_CENTRE_CELL);

    char iname[] = "t";
    status = arcus_add_scalar_double(c_reg, iname, &t);

    char lname[] = "clogger";
    iris_log_t c_log = iris_create_logger(lname, IRIS_LEVEL_INFO);

    status = iris_add_registry(c_log, c_reg);

    char fmt[] = "%y %m%n";
    status = iris_add_console_appender(c_log, fmt);


    // ------------------------------------------------------------------------


    iris::Logger * coreLogger = new iris::Logger("BaseLogger", iris::Levels::Info);

    coreLogger->addConsoleAppender("C %m%n");
    coreLogger->setRegistry(log_reg);
    // coreLogger->setMPIAll(false);

    iris::Logger * spLogger = new iris::Logger("SPrintLogger", iris::Levels::Info);

    spLogger->addConsoleAppender(":S: %m%n");
    spLogger->addFileAppender(":S: %m%n", "sprint.dat");
    spLogger->setRegistry(log_reg);
    spLogger->setDefaultLevel(iris::Levels::Info);

    mesh.dumpToSILO(0.0, 0, coreLogger);

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

        coreLogger->logvar(iris::Levels::Info, "Step: @{step:%-5d} :: Time = @{t}, dt = @{dt:%e}");

        char msg[] = "t = @{t:%e}";
        status = iris_log(c_log, IRIS_LEVEL_INFO, msg);

        // Sweep is half X, then Y, then half X
        mesh.sweepX(dt/2.0, Hydro::MUSCLHancock1D);
        mesh.sweepY(dt, Hydro::MUSCLHancock1D);
        mesh.sweepX(dt/2.0, Hydro::MUSCLHancock1D);

        t += dt;
        if (t >= outNext) {
            outNext += dtOut;
            if (myrank == 0) {
                mesh.dumpToSILO(t, step, coreLogger);
            }
            spLogger->GenerateRegionalPrint(
                "regno|Region:%6d density|Rho momu|momU momv|momV e|Energy:%12.6e"
            );

            spLogger->GenerateMaterialPrint(
                "matno|Material:%8d density|Rho momu|momU momv|momV e|Energy"
            );
        }
        if (t > tEnd || step > 10000) break;
    }

    mesh.Kill();

    MPI::Finalize();
    return 0;
}
