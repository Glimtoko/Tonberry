#include "mesh2d.hpp"
#include "update.hpp"

#include <iostream>

int main() {


//     Mesh2D mesh(3, 100, 0);
    Mesh2D mesh(100, 100, 2);

    mesh.dumpToNetCDF();
//     mesh.dumpToNetCDF_NG();


    double t = 0.0;
    int step = 0;

    for( ; ; ) {
        step++;

        double dt = getTimestep(mesh);
        std::cout << "Step: " << step;
        std::cout << ", time = " << t;
        std::cout << ", dt = " << dt << std::endl;

        sweep(
            mesh.rho, mesh.momU, mesh.momV, mesh.E,
            mesh.niGhosts, mesh.njGhosts, mesh.iUpper, mesh.jUpper,
            Sweep::x,
            mesh.dx, 0.0, dt, mesh.gamma);
        mesh.setBoundaries();

        sweep(
            mesh.rho, mesh.momV, mesh.momU, mesh.E,
            mesh.niGhosts, mesh.njGhosts, mesh.iUpper, mesh.jUpper,
            Sweep::y,
            mesh.dy, 0.0, dt, mesh.gamma);
        mesh.setBoundaries();

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
