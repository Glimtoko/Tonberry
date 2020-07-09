#include "mesh2d.hpp"
#include "update.hpp"

#include <iostream>

int main() {


    // 100x100 mesh using spherical Sod set-up
    Mesh2D mesh(100, 100, 2);

    mesh.dumpToNetCDF();

    double t = 0.0;
    int step = 0;

    for( ; ; ) {
        step++;

        double dt = getTimestep(mesh);
        std::cout << "Step: " << step;
        std::cout << ", time = " << t;
        std::cout << ", dt = " << dt << std::endl;

        // Sweep is half X, then Y, then half X
        sweep(mesh, Sweep::x, dt/2.0);
        sweep(mesh, Sweep::y, dt);
        sweep(mesh, Sweep::x, dt/2.0);

        t += dt;
        if (t > 0.25 || step > 2000) break;
    }

    std::cout << step << std::endl;
    mesh.dumpToNetCDF_NG();

//     for (int i=2; i<mesh.jUpper; i++) {
//         std::cout << i << " " << mesh.y[i] << " " << mesh.rho[i][2] << " " << mesh.momV[i][2] << std::endl;
//     }
    return 0;
}