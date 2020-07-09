#include <valarray>

#include "ncFile.h"
#include "ncDim.h"
#include "ncVar.h"

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
    x.resize(niGhosts);

    // Ghost/boundary cells mean that the first "data" cell is at index 2.
    // Coordinates are for cell centres, so start at dx/2, not 0
    for (int i=0; i<niGhosts; i++) {
        x[i] = (0.5 + i - 2)*dx;
    }

    // Same for y
    double dy = yLength/ni;
    y.resize(njGhosts);
    for (int j=0; j<njGhosts; j++) {
        y[j] = (0.5 + j - 2)*dy;
    }

    // Set physical arrays
    rho.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, v));
    momU.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, v));
    momV.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, v));
    E.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, v));
}


Mesh2D::Mesh2D(int ni, int nj, int xy) {
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
    dx = xLength/ni;
    x.resize(niGhosts);

    // Ghost/boundary cells mean that the first "data" cell is at index 2.
    // Coordinates are for cell centres, so start at dx/2, not 0
    for (int i=0; i<niGhosts; i++) {
        x[i] = (0.5 + i - 2)*dx;
    }

    // Same for y
    dy = yLength/nj;
    y.resize(njGhosts);
    for (int j=0; j<njGhosts; j++) {
        y[j] = (0.5 + j - 2)*dy;
    }

    // Set physical arrays
    rho.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, 0.0));
    momU.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, 0.0));
    momV.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, 0.0));
    E.resize(nj+2*nghosts, std::vector<double>(ni+2*nghosts, 0.0));

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
    }

    // Set boundary conditions
    setBoundaries();

}

void Mesh2D::setBoundaries() {
    for (int j=2; j<jUpper; j++) {
        rho[j][1] = rho[j][2];
        rho[j][0] = rho[j][3];
        rho[j][niGhosts-2] = rho[j][niGhosts-3];
        rho[j][niGhosts-1] = rho[j][niGhosts-4];

        momU[j][1] = momU[j][2];
        momU[j][0] = momU[j][3];
        momU[j][niGhosts-2] = momU[j][niGhosts-3];
        momU[j][niGhosts-1] = momU[j][niGhosts-4];

        momV[j][1] = momV[j][2];
        momV[j][0] = momV[j][3];
        momV[j][niGhosts-2] = momV[j][niGhosts-3];
        momV[j][niGhosts-1] = momV[j][niGhosts-4];

        E[j][1] = E[j][2];
        E[j][0] = E[j][3];
        E[j][niGhosts-2] = E[j][niGhosts-3];
        E[j][niGhosts-1] = E[j][niGhosts-4];
    }

    for (int i=2; i<iUpper; i++) {
        rho[1][i] = rho[2][i];
        rho[0][i] = rho[3][i];
        rho[njGhosts-2][i] = rho[njGhosts-3][i];
        rho[njGhosts-1][i] = rho[njGhosts-4][i];

        momU[1][i] = momU[2][i];
        momU[0][i] = momU[3][i];
        momU[njGhosts-2][i] = momU[njGhosts-3][i];
        momU[njGhosts-1][i] = momU[njGhosts-4][i];

        momV[1][i] = momV[2][i];
        momV[0][i] = momV[3][i];
        momV[njGhosts-2][i] = momV[njGhosts-3][i];
        momV[njGhosts-1][i] = momV[njGhosts-4][i];

        E[1][i] = E[2][i];
        E[0][i] = E[3][i];
        E[njGhosts-2][i] = E[njGhosts-3][i];
        E[njGhosts-1][i] = E[njGhosts-4][i];
    }
}

void Mesh2D::dumpToNetCDF() {
    // Create a 2D array for holding data
    double test[this->njGhosts][this->niGhosts];

    netCDF::NcFile meshFile("initial.nc", netCDF::NcFile::FileMode::replace);

    // Define axis sizes
    netCDF::NcDim yDim = meshFile.addDim("Y", this->njGhosts);
    netCDF::NcDim xDim = meshFile.addDim("X", this->niGhosts);

    // Define variables to hold axis values
    netCDF::NcVar yVar = meshFile.addVar("Y", netCDF::ncDouble, yDim);
    netCDF::NcVar xVar = meshFile.addVar("X", netCDF::ncDouble, xDim);

    xVar.putVar(this->x.data());
    yVar.putVar(this->y.data());

    xVar.putAtt("units", "cm");
    yVar.putAtt("units", "cm");

    std::vector<netCDF::NcDim> dims;
    dims.push_back(yDim);
    dims.push_back(xDim);

    // Write out a variable
    netCDF::NcVar rhoVar = meshFile.addVar("Rho", netCDF::ncDouble, dims);
    WRITEFULL(rhoVar, this->rho);
    rhoVar.putAtt("units", "g/cc");

    netCDF::NcVar momUVar = meshFile.addVar("Momentum", netCDF::ncDouble, dims);
    WRITEFULL(momUVar, this->momU);

    netCDF::NcVar EVar = meshFile.addVar("Energy", netCDF::ncDouble, dims);
    WRITEFULL(EVar, this->E);
}

template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;

	std::vector<T> vec(first, last);
	return vec;
}

void Mesh2D::dumpToNetCDF_NG() {
    // Create a 2D array for holding data
    int jSize = this->jUpper - this->nghosts;
    int iSize = this->iUpper - this->nghosts;
    double test[jSize][iSize];

    netCDF::NcFile meshFile("initial_ng.nc", netCDF::NcFile::FileMode::replace);

    // Define axis sizes
    netCDF::NcDim yDim = meshFile.addDim("Y", jSize);
    netCDF::NcDim xDim = meshFile.addDim("X", iSize);

    // Define variables to hold axis values
    netCDF::NcVar yVar = meshFile.addVar("Y", netCDF::ncDouble, yDim);
    netCDF::NcVar xVar = meshFile.addVar("X", netCDF::ncDouble, xDim);

    xVar.putVar(slice(this->x, 2, iUpper).data());
    yVar.putVar(slice(this->y, 2, jUpper).data());

    xVar.putAtt("units", "cm");
    yVar.putAtt("units", "cm");

    std::vector<netCDF::NcDim> dims;
    dims.push_back(yDim);
    dims.push_back(xDim);

    // Write out a variable
    netCDF::NcVar rhoVar = meshFile.addVar("Rho", netCDF::ncDouble, dims);
    WRITENG(rhoVar, this->rho);
    rhoVar.putAtt("units", "g/cc");

    netCDF::NcVar momUVar = meshFile.addVar("Momentum", netCDF::ncDouble, dims);
    WRITENG(momUVar, this->momU);

    netCDF::NcVar EVar = meshFile.addVar("Energy", netCDF::ncDouble, dims);
    WRITENG(EVar, this->E);
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
