#include <valarray>

#include "ncFile.h"
#include "ncDim.h"
#include "ncVar.h"

#include <silo.h>

#include <math.h>
#include <string.h>

#include "mesh2d.hpp"

#define WRITEFULL(A,B) for (int j=0; j<this->njGhosts; j++) {for (int i=0; i<this->niGhosts; i++) {test[j][i] = B[j][i];}} A.putVar(test);

#define WRITENG(A,B) for (int j=2; j<this->jUpper; j++) {for (int i=2; i<this->jUpper; i++) {test[j-2][i-2] = B[j][i];}} A.putVar(test);



Mesh2D::Mesh2D(int ni, int nj, double v) {
    // Set X and Y lengths for Sod
    double xLength = 1.0;
    double yLength = 1.0;

    // Number of ghosts on each mesh side
    niGhosts = ni + 2*nghosts;
    njGhosts = nj + 2*nghosts;
    iUpper = ni + nghosts;
    jUpper = nj + nghosts;

    // Set coordinate arrays. Since this is a cartesian mesh, these need only
    // be 1D arrays, as, for example, all cells with the same j-index will have
    // the same y position.
    double dx = xLength/ni;
    x = new double(niGhosts);

    // Ghost/boundary cells mean that the first "data" cell is at index 2.
    // Coordinates are for cell centres, so start at dx/2, not 0
    for (int i=0; i<niGhosts; i++) {
        x[i] = (0.5 + i - 2)*dx;
    }

    // Same for y
    double dy = yLength/ni;
    y = new double(njGhosts);
    for (int j=0; j<njGhosts; j++) {
        y[j] = (0.5 + j - 2)*dy;
    }

    // Set physical arrays
    ALLOCATE(rho, niGhosts, njGhosts);
    ALLOCATE(momU, niGhosts, njGhosts);
    ALLOCATE(momV, niGhosts, njGhosts);
    ALLOCATE(E, niGhosts, njGhosts);
}


Mesh2D::Mesh2D(int ni, int nj, int xy) {
    // Set X and Y lengths for Sod
    double xLength = 1.0;
    double yLength = 1.0;

    if (xy == 3) {
        xLength = 20.0;
        yLength = 10.0;
    }

    if (xy == 4) {
        xLength = 30.0;
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
    ALLOCATE(rho, niGhosts, njGhosts);
    ALLOCATE(momU, niGhosts, njGhosts);
    ALLOCATE(momV, niGhosts, njGhosts);
    ALLOCATE(E, niGhosts, njGhosts);

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
        for (int i=2; i<iUpper; i++) {
            double cellRho = x[i] < x0 ? rhoL : rhoR;
            double cellMom = x[i] < x0 ? rhoL * uL : rhoR * uR;
            double e = x[i] < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
            // E has a w**2 term, but it's always zero here... (actually, so is u...)
            double cellE = x[i] < x0 ? rhoL*(0.5*uL*uL + e) : rhoR*(0.5*uR*uR + e);

            for (int j=2; j<jUpper; j++) {
                rho[j][i] = cellRho;
                momU[j][i] = cellMom;
                momV[j][i] = 0.0;
                E[j][i] = cellE;
            }
        }
    } else if (xy==1) {
        for (int j=2; j<jUpper; j++) {
            double cellRho = y[j] < x0 ? rhoL : rhoR;
            double cellMom = y[j] < x0 ? rhoL * uL : rhoR * uR;
            double e = y[j] < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
            // E has a w**2 term, but it's always zero here... (actually, so is u...)
            double cellE = y[j] < x0 ? rhoL*(0.5*uL*uL + e) : rhoR*(0.5*uR*uR + e);
            for (int i=2; i<iUpper; i++) {
                rho[j][i] = cellRho;
                momV[j][i] = cellMom;
                momU[j][i] = 0.0;
                E[j][i] = cellE;
            }
        }
    } else if (xy == 2) {
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
                double r = sqrt(y[j]*y[j] + x[i]*x[i]);
                double cellRho = r < x0 ? rhoL : rhoR;

                double e = r < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
                double cellE = r < x0 ? rhoL*e : rhoR*e;

                rho[j][i] = cellRho;
                momV[j][i] = 0.0;
                momU[j][i] = 0.0;
                E[j][i] = cellE;
            }
        }
    } else if (xy == 3) {
        double pi = 3.14159265359;
        double B = 5.0;
        double x0 = 5, y0 = 5;
        double u0 = 1.0;
        double v0 = 0.0;
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
                double r = sqrt(pow(y[j] - y0, 2) + pow(x[i] - x0, 2));
                double f = B/(2*pi) * exp((1-r*r)/2);
                double u = u0 + (x[i] - x0)*f;
                double v = v0 + (y[j] - y0)*f;
                double T = 1 - ((gamma-1)*B*B)/(8*gamma*pi*pi)*exp(1-r*r);

                rho[j][i] = pow(T, (1/(gamma-1)));
                momU[j][i] = rho[j][i]*u;
                momV[j][i] = rho[j][i]*v;

                double p = rho[j][i]*T;
                double e = p/((gamma - 1)*rho[j][i]);
                E[j][i] = rho[j][i]*(e + 0.5*u*u + 0.5*v*v);
            }
        }
    } else if (xy == 4) {
        double x0 = 0, y0 = 5, x1 = 10.0;
        double p;
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
                double r = sqrt(pow(y[j] - y0, 2) + pow(x[i] - x0, 2));
                if (r <= 2.0) {
                    rho[j][i] = 1.0;
                    p = 1.0;
                } else {
                    rho[j][i] = 0.125;
                    p = 0.1;
                }

                r = sqrt(pow(y[j] - y0, 2) + pow(x[i] - x1, 2));
                if (r <= 2.0) rho[j][i] = 0.5;

                double e = p/((gamma - 1.0)*rho[j][i]);
                E[j][i] = e*rho[j][i];
            }
        }
        // Set reflective boundaries top and bottom
        bU = -1.0;
        bD = -1.0;
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
        double f = -1.0;
        rho[1][i] = rho[2][i];
        rho[0][i] = rho[3][i];
        rho[njGhosts-2][i] = rho[njGhosts-3][i];
        rho[njGhosts-1][i] = rho[njGhosts-4][i];

        momU[1][i] = momU[2][i];
        momU[0][i] = momU[3][i];
        momU[njGhosts-2][i] = momU[njGhosts-3][i];
        momU[njGhosts-1][i] = momU[njGhosts-4][i];

        momV[1][i] = meshBoundaryUD[0]*momV[2][i];
        momV[0][i] = meshBoundaryUD[0]*momV[3][i];
        momV[njGhosts-2][i] = meshBoundaryUD[1]*momV[njGhosts-3][i];
        momV[njGhosts-1][i] = meshBoundaryUD[1]*momV[njGhosts-4][i];

        E[1][i] = E[2][i];
        E[0][i] = E[3][i];
        E[njGhosts-2][i] = E[njGhosts-3][i];
        E[njGhosts-1][i] = E[njGhosts-4][i];
    }

    for (int j=2; j<jUpper; j++) {
        double f = 1.0;
        rho[j][1] = rho[j][2];
        rho[j][0] = rho[j][3];
        rho[j][niGhosts-2] = rho[j][niGhosts-3];
        rho[j][niGhosts-1] = rho[j][niGhosts-4];

        momU[j][1] = meshBoundaryLR[0]*momU[j][2];
        momU[j][0] = meshBoundaryLR[0]*momU[j][3];
        momU[j][niGhosts-2] = meshBoundaryLR[1]*momU[j][niGhosts-3];
        momU[j][niGhosts-1] = meshBoundaryLR[1]*momU[j][niGhosts-4];

        momV[j][1] = momV[j][2];
        momV[j][0] = momV[j][3];
        momV[j][niGhosts-2] = momV[j][niGhosts-3];
        momV[j][niGhosts-1] = momV[j][niGhosts-4];

        E[j][1] = E[j][2];
        E[j][0] = E[j][3];
        E[j][niGhosts-2] = E[j][niGhosts-3];
        E[j][niGhosts-1] = E[j][niGhosts-4];
    }

}

template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;

	std::vector<T> vec(first, last);
	return vec;
}

void Mesh2D::dumpToSILO(double time) {
    char stateNo[4];
    this->dumpStateNo++;
    sprintf(stateNo, "%03d", this->dumpStateNo);

    char fileName[20];
    strcpy(fileName, "tonberry");
    strcat(fileName, stateNo);
    strcat(fileName, ".silo");

    std::cout << "Outputting to " << fileName << std::endl;

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

    // Write mesh to file
    DBPutQuadmesh(file, "mesh1", coordnames, coordinates, nodeDims, 2, DB_DOUBLE, DB_COLLINEAR, NULL);

    // Test writing a quant
    int cellDims[2] = {iSize, jSize};
    double *data = new double[cellDims[0]*cellDims[1]];
    int index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            data[index++] = this->rho[j][i];
        }
    }

    DBPutQuadvar1(file, "Rho", "mesh1", data, cellDims, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);


    // Close file
    DBClose(file);
}

QUANT_2D &Mesh2D::getMomentum(Sweep sweep, Direction direction) {
    switch (sweep) {
        case Sweep::x:
            switch (direction) {
                case normal:
                    return momU;
                default:
                    return momV;
            }
        default:
            switch (direction) {
                case normal:
                    return momV;
                default:
                    return momU;
            }
    }
}

void Mesh2D::Kill() {
    delete[] x;
    delete[] y;

    delete[] rho;
    delete[] momU;
    delete[] momV;
    delete[] E;
}
