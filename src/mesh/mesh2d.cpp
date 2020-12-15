#include <valarray>

#include "ncFile.h"
#include "ncDim.h"
#include "ncVar.h"

#include "silo.h"
#include "mpi.h"

#include <math.h>
#include <string.h>

#include "mesh2d.hpp"

Mesh2D::Mesh2D(int ni, int nj, int xy) {
    // Set X and Y lengths for Sod
    double xLength = 1.0;
    double yLength = 1.0;

    if (xy == 3) {
        xLength = 20.0;
        yLength = 10.0;
    }

    if (xy == 4) {
        xLength = 40.0;
        yLength = 10.0;
    }

    // Number of ghosts on each mesh side
    niGhosts = ni + 2*nghosts;
    njGhosts = nj + 2*nghosts;
    iUpper = ni + nghosts;
    jUpper = nj + nghosts;

    // Set coordinate arrays. Since this is a cartesian mesh, these need only
    // be 1D arrays, as, for example, all cells with the same j-index will have
    // the same y position.
    dx = xLength/ni;
    x = new double[niGhosts];

    // Ghost/boundary cells mean that the first "data" cell is at index 2.
    // Coordinates are for cell centres, so start at dx/2, not 0
    for (int i=0; i<niGhosts; i++) {
        x[i] = (0.5 + i - 2)*dx;
    }

    // Same for y
    dy = yLength/nj;
    y = new double[njGhosts];
    for (int j=0; j<njGhosts; j++) {
        y[j] = (0.5 + j - 2)*dy;
    }

    // Set physical arrays
    rho = new double[njGhosts*niGhosts];
    momU = new double[njGhosts*niGhosts];
    momV = new double[njGhosts*niGhosts];
    E = new double[njGhosts*niGhosts];

    // Set material and region arrays
    materials = new int[njGhosts*niGhosts];
    regions = new int[njGhosts*niGhosts];
    
    // Set initial cell map
    cell_map = new int[njGhosts*niGhosts];
    for (int i=0; i<niGhosts; i++) {
        for (int j=0; j<njGhosts; j++) {
            _LGET(cell_map, j, i) = 0;
        }
    }


    // Boundary values
    double bL = 1.0;
    double bR = 1.0;
    double bU = 1.0;
    double bD = 1.0;

    // Physical parameters for Sod
    const double x0 = 0.5;
    const double uL = 0.0;
    const double uR = 0.0;
    const double rhoL = 1.0;
    const double rhoR = 0.125;
    const double pL = 1.0;
    const double pR = 0.1;
    gamma = 1.4;
    dtmax = 0.1;
    cfl = 0.6;

    if (xy == 0) {
        nreg = 2;
        for (int i=2; i<iUpper; i++) {
            double cellRho = x[i] < x0 ? rhoL : rhoR;
            double cellMom = x[i] < x0 ? rhoL * uL : rhoR * uR;
            double e = x[i] < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
            // E has a w**2 term, but it's always zero here... (actually, so is u...)
            double cellE = x[i] < x0 ? rhoL*(0.5*uL*uL + e) : rhoR*(0.5*uR*uR + e);

            int regmat = x[i] < x0 ? 1 : 2;

            for (int j=2; j<jUpper; j++) {
                _LGET(cell_map, j, i) = 1;
                _LGET(rho, j, i) = cellRho;
                _LGET(momU, j, i) = cellMom;
                _LGET(momV, j, i) = 0.0;
                _LGET(E, j, i) = cellE;
                _LGET(materials, j, i) = regmat;
                _LGET(regions, j, i) = regmat;
            }
        }
    } else if (xy==1) {
        nreg = 2;
        for (int j=2; j<jUpper; j++) {
            double cellRho = y[j] < x0 ? rhoL : rhoR;
            double cellMom = y[j] < x0 ? rhoL * uL : rhoR * uR;
            double e = y[j] < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
            // E has a w**2 term, but it's always zero here... (actually, so is u...)
            double cellE = y[j] < x0 ? rhoL*(0.5*uL*uL + e) : rhoR*(0.5*uR*uR + e);

            int regmat = x[j] < x0 ? 1 : 2;

            for (int i=2; i<iUpper; i++) {
                _LGET(cell_map, j, i) = 1;
                _LGET(rho, j, i) = cellRho;
                _LGET(momV, j, i) = cellMom;
                _LGET(momU, j, i) = 0.0;
                _LGET(E, j, i) = cellE;
                _LGET(materials, j, i) = regmat;
                _LGET(regions, j, i) = regmat;
            }
        }
    } else if (xy == 2) {
        nreg = 2;
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
                _LGET(cell_map, j, i) = 1;
                double r = sqrt(y[j]*y[j] + x[i]*x[i]);
                double cellRho = r < x0 ? rhoL : rhoR;

                double e = r < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
                double cellE = r < x0 ? rhoL*e : rhoR*e;

                int regmat = r < x0 ? 1 : 2;

                _LGET(rho, j, i) = cellRho;
                _LGET(momV, j, i) = 0.0;
                _LGET(momU, j, i) = 0.0;
                _LGET(E, j, i) = cellE;
                _LGET(materials, j, i) = regmat;
                _LGET(regions, j, i) = regmat;
            }
        }
    } else if (xy == 3) {
        double pi = 3.14159265359;
        double B = 5.0;
        double x0 = 5, y0 = 5;
        double u0 = 1.0;
        double v0 = 0.0;
        nreg = 1;
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
                _LGET(cell_map, j, i) = 1;
                double r = sqrt(pow(y[j] - y0, 2) + pow(x[i] - x0, 2));
                double f = B/(2*pi) * exp((1-r*r)/2);
                double u = u0 + (x[i] - x0)*f;
                double v = v0 + (y[j] - y0)*f;
                double T = 1 - ((gamma-1)*B*B)/(8*gamma*pi*pi)*exp(1-r*r);

                _LGET(rho, j, i) = pow(T, (1/(gamma-1)));
                _LGET(momU, j, i) = _LGET(rho, j, i)*u;
                _LGET(momV, j, i) = _LGET(rho, j, i)*v;

                double p = _LGET(rho, j, i)*T;
                double e = p/((gamma - 1)*_LGET(rho, j, i));
                _LGET(E, j, i) = _LGET(rho, j, i)*(e + 0.5*u*u + 0.5*v*v);
                _LGET(materials, j, i) = 1;
                _LGET(regions, j, i) = 1;
            }
        }
    } else if (xy == 4) {
        double x0 = 0, y0 = 5, x1 = 10.0, x2 = 20.0, y1 = 3.0, y2 = 7.0;
        double p;
        nreg = 3;
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
                _LGET(cell_map, j, i) = 1;
                double r = sqrt(pow(y[j] - y0, 2) + pow(x[i] - x0, 2));
                if (r <= 2.0) {
                    _LGET(rho, j, i) = 1.0;
                    p = 1.0;
                } else {
                    _LGET(rho, j, i) = 0.125;
                    p = 0.1;
                }

                r = sqrt(pow(y[j] - y0, 2) + pow(x[i] - x1, 2));
                if (r <= 2.0) _LGET(rho, j, i) = 0.5;

                r = sqrt(pow(y[j] - y1, 2) + pow(x[i] - x2, 2));
                if (r <= 1.5) _LGET(rho, j, i) = 0.5;

                r = sqrt(pow(y[j] - y2, 2) + pow(x[i] - x2, 2));
                if (r <= 1.5) _LGET(rho, j, i) = 0.5;

                double e = p/((gamma - 1.0)*_LGET(rho, j, i));
                _LGET(E, j, i) = e*_LGET(rho, j, i);
                _LGET(momU, j, i) = 0.0;
                _LGET(momV, j, i) = 0.0;
                _LGET(materials, j, i) = 1;

                if (x[i] < x1) {
                    _LGET(regions, j, i) = 1;
                } else if (x[i] < x2) {
                    _LGET(regions, j, i) = 2;
                } else {
                    _LGET(regions, j, i) = 3;
                }
            }
        }
        // Set reflective boundaries top, bottom and right
        bU = -1.0;
        bD = -1.0;
        bR = -1.0;
    }

    // Set boundary factor arrays
    meshBoundaryLR = new double[2];
    meshBoundaryUD = new double[2];

    // Initialise all boundaries to be transmissive
    meshBoundaryLR[0] = bL;
    meshBoundaryLR[1] = bR;

    meshBoundaryUD[0] = bD;
    meshBoundaryUD[1] = bU;

    // Set boundary conditions
    setBoundaries();

}

void Mesh2D::setBoundaries() {
    for (int i=2; i<iUpper; i++) {
        _LGET(rho, 1, i) = _LGET(rho, 2, i);
        _LGET(rho, 0, i) = _LGET(rho, 3, i);
        _LGET(rho, njGhosts-2, i) = _LGET(rho, njGhosts-3, i);
        _LGET(rho, njGhosts-1, i) = _LGET(rho, njGhosts-4, i);

        _LGET(momU, 1, i) = _LGET(momU, 2, i);
        _LGET(momU, 0, i) = _LGET(momU, 3, i);
        _LGET(momU, njGhosts-2, i) = _LGET(momU, njGhosts-3, i);
        _LGET(momU, njGhosts-1, i) = _LGET(momU, njGhosts-4, i);

        _LGET(momV, 1, i) = meshBoundaryUD[0]*_LGET(momV, 2, i);
        _LGET(momV, 0, i) = meshBoundaryUD[0]*_LGET(momV, 3, i);
        _LGET(momV, njGhosts-2, i) = meshBoundaryUD[1]*_LGET(momV, njGhosts-3, i);
        _LGET(momV, njGhosts-1, i) = meshBoundaryUD[1]*_LGET(momV, njGhosts-4, i);

        _LGET(E, 1, i) = _LGET(E, 2, i);
        _LGET(E, 0, i) = _LGET(E, 3, i);
        _LGET(E, njGhosts-2, i) = _LGET(E, njGhosts-3, i);
        _LGET(E, njGhosts-1, i) = _LGET(E, njGhosts-4, i);
    }

    for (int j=2; j<jUpper; j++) {
        _LGET(rho, j, 1) = _LGET(rho, j, 2);
        _LGET(rho, j, 0) = _LGET(rho, j, 3);
        _LGET(rho, j, niGhosts-2) = _LGET(rho, j, niGhosts-3);
        _LGET(rho, j, niGhosts-1) = _LGET(rho, j, niGhosts-4);

        _LGET(momU, j, 1) = meshBoundaryLR[0]*_LGET(momU, j, 2);
        _LGET(momU, j, 0) = meshBoundaryLR[0]*_LGET(momU, j, 3);
        _LGET(momU, j, niGhosts-2) = meshBoundaryLR[1]*_LGET(momU, j, niGhosts-3);
        _LGET(momU, j, niGhosts-1) = meshBoundaryLR[1]*_LGET(momU, j, niGhosts-4);

        _LGET(momV, j, 1) = _LGET(momV, j, 2);
        _LGET(momV, j, 0) = _LGET(momV, j, 3);
        _LGET(momV, j, niGhosts-2) = _LGET(momV, j, niGhosts-3);
        _LGET(momV, j, niGhosts-1) = _LGET(momV, j, niGhosts-4);

        _LGET(E, j, 1) = _LGET(E, j, 2);
        _LGET(E, j, 0) = _LGET(E, j, 3);
        _LGET(E, j, niGhosts-2) = _LGET(E, j, niGhosts-3);
        _LGET(E, j, niGhosts-1) = _LGET(E, j, niGhosts-4);
    }

}

void Mesh2D::dumpToSILO(double time, int step, std::shared_ptr<iris::Logger> logger) {
    char stateNo[4];
    this->dumpStateNo++;
    sprintf(stateNo, "%03d", this->dumpStateNo);

    char fileName[20];
    strcpy(fileName, "tonberry");
    strcat(fileName, stateNo);
    strcat(fileName, ".silo");

    logger->log(iris::Levels::Info, "Outputting to " + std::string(fileName));

    DBfile *file = NULL;
    file = DBCreate(fileName, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);


    // Size of mesh, ignoring boundary cells
    int jSize = this->jUpper - this->nghosts;
    int iSize = this->iUpper - this->nghosts;

    // Create node centred Coordinates
    int nodeDims[2] = {iSize + 1, jSize + 1};
    double *x = new double[nodeDims[0]];
    double *y = new double[nodeDims[1]];

    double *coordinates[2];
    char *coordnames[2];

    for (int i=0; i<nodeDims[0]; i++) {
        x[i] = (i)*this->dx;
    }

    for (int j=0; j<nodeDims[1]; j++) {
        y[j] = (j)*this->dy;
    }

    coordinates[0] = x;
    coordinates[1] = y;
    coordnames[0] = strdup("X");
    coordnames[1] = strdup("Y");

    // Create option list
    DBoptlist *optList = DBMakeOptlist(2);
    DBAddOption(optList, DBOPT_DTIME, &time);
    DBAddOption(optList, DBOPT_CYCLE, &step);

    // Write mesh to file
    DBPutQuadmesh(file, "mesh1", coordnames, coordinates, nodeDims, 2, DB_DOUBLE, DB_COLLINEAR, optList);

    // Set quant storage
    int cellDims[2] = {iSize, jSize};
    double *data = new double[cellDims[0]*cellDims[1]];


    // Write density
    int index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            data[index++] = _PGET(this, rho, j, i);
        }
    }
    DBPutQuadvar1(file, "Density", "mesh1", data, cellDims, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

    // Write pressure
    index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            double rho = _PGET(this, rho, j, i);
            double u = _PGET(this, momU, j, i)/rho;
            double v = _PGET(this, momV, j, i)/rho;
            double p = (this->gamma - 1.0)*(_PGET(this, E, j, i) - 0.5*v*v - 0.5*u*u);

            data[index++] = p;
        }
    }
    DBPutQuadvar1(file, "Pressure", "mesh1", data, cellDims, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

    // Close file
    DBClose(file);

    // Free storage
    delete[] data;
}


void Mesh2D::Kill() {
    delete[] x;
    delete[] y;

    delete[] rho;
    delete[] momU;
    delete[] momV;
    delete[] E;
}


void Mesh2D::sweepX(double dt,
            void(*flux)(double*, double*, double*, double*, int, int, double, double, double)
){
    int ni = this->niGhosts;
    int nj = this->jUpper - 2;
    int iUpper = this->iUpper;
    int jUpper = this->jUpper;
    double gamma = this->gamma;
    double dx = this->dx;

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
            rhoAll[index] = _PGET(this, rho, j, i);
            momNAll[index] = _PGET(this, momU, j, i);
            momTAll[index] = _PGET(this, momV, j, i);
            EAll[index] = _PGET(this, E, j, i);

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
        flux(rho, E, momN, momT, ni, iUpper, gamma, dt, dx);

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
            _PGET(this, rho, j, i) = rhoAll[index];
            _PGET(this, momU, j, i) = momNAll[index];
            _PGET(this, momV, j, i) = momTAll[index];
            _PGET(this, E, j, i) = EAll[index];
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
    this->setBoundaries();
}


void Mesh2D::sweepY(double dt,
            void(*flux)(double*, double*, double*, double*, int, int, double, double, double)
) {
    int nj = this->njGhosts;
    int ni = this->iUpper - 2;
    int iUpper = this->iUpper;
    int jUpper = this->jUpper;
    double gamma = this->gamma;
    double dy = this->dy;

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
            rhoAll[index] = _PGET(this, rho, j, i);
            momNAll[index] = _PGET(this, momV, j, i);
            momTAll[index] = _PGET(this, momU, j, i);
            EAll[index] = _PGET(this, E, j, i);

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
        flux(rho, E, momN, momT, nj, jUpper, gamma, dt, dy);

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
            _PGET(this, rho, j, i) = rhoAll[index];
            _PGET(this, momU, j, i) = momTAll[index];
            _PGET(this, momV, j, i) = momNAll[index];
            _PGET(this, E, j, i) = EAll[index];
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
    this->setBoundaries();
}
