#include <valarray>

#ifdef HasSILO
#include "silo.h"
#endif

#ifdef HasTIO
#include "typhonio.h"
#endif

#include <math.h>
#include <string.h>
#include <iostream>

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
    x.reserve(niGhosts);

    // Ghost/boundary cells mean that the first "data" cell is at index 2.
    // Coordinates are for cell centres, so start at dx/2, not 0
    for (int i=0; i<niGhosts; i++) {
        x[i] = (0.5 + i - 2)*dx;
    }

    // Same for y
    dy = yLength/nj;
    y.reserve(njGhosts);
    for (int j=0; j<njGhosts; j++) {
        y[j] = (0.5 + j - 2)*dy;
    }

    // Set physical arrays
    rho.reserve(njGhosts*niGhosts);
    momU.reserve(njGhosts*niGhosts);
    momV.reserve(njGhosts*niGhosts);
    E.reserve(njGhosts*niGhosts);

    for (int i=0; i<njGhosts*niGhosts; i++) {
        rho[i] = 0.0001;
        momU[i] = 0.0001;
        momV[i] = 0.0001;
        E[i] = 0.0001;
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
        for (int i=2; i<iUpper; i++) {
            double cellRho = x[i] < x0 ? rhoL : rhoR;
            double cellMom = x[i] < x0 ? rhoL * uL : rhoR * uR;
            double e = x[i] < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
            // E has a w**2 term, but it's always zero here... (actually, so is u...)
            double cellE = x[i] < x0 ? rhoL*(0.5*uL*uL + e) : rhoR*(0.5*uR*uR + e);

            for (int j=2; j<jUpper; j++) {
                _LGET(rho, j, i) = cellRho;
                _LGET(momU, j, i) = cellMom;
                _LGET(momV, j, i) = 0.0;
                _LGET(E, j, i) = cellE;
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
                _LGET(rho, j, i) = cellRho;
                _LGET(momV, j, i) = cellMom;
                _LGET(momU, j, i) = 0.0;
                _LGET(E, j, i) = cellE;
            }
        }
    } else if (xy == 2) {
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
                double r = sqrt(y[j]*y[j] + x[i]*x[i]);
                double cellRho = r < x0 ? rhoL : rhoR;

                double e = r < x0 ? pL/((gamma - 1.0)*rhoL) : pR/((gamma - 1.0)*rhoR);
                double cellE = r < x0 ? rhoL*e : rhoR*e;

                _LGET(rho, j, i) = cellRho;
                _LGET(momV, j, i) = 0.0;
                _LGET(momU, j, i) = 0.0;
                _LGET(E, j, i) = cellE;
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

                _LGET(rho, j, i) = pow(T, (1/(gamma-1)));
                _LGET(momU, j, i) = _LGET(rho, j, i)*u;
                _LGET(momV, j, i) = _LGET(rho, j, i)*v;

                double p = _LGET(rho, j, i)*T;
                double e = p/((gamma - 1)*_LGET(rho, j, i));
                _LGET(E, j, i) = _LGET(rho, j, i)*(e + 0.5*u*u + 0.5*v*v);
            }
        }
    } else if (xy == 4) {
        double x0 = 0, y0 = 5, x1 = 10.0, x2 = 20.0, y1 = 3.0, y2 = 7.0;
        double p;
        for (int j=2; j<jUpper; j++) {
            for (int i=2; i<iUpper; i++) {
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
            }
        }
        // Set reflective boundaries top, bottom and right
        bU = -1.0;
        bD = -1.0;
        bR = -1.0;
    }

    // Set boundary factor arrays
    meshBoundaryLR.reserve(2);
    meshBoundaryUD.reserve(2);

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

#ifdef HasTIO
void Mesh2D::dumpToTIO(double time, int step) {
    char stateNo[4];
    this->dumpStateNoTIO++;
    sprintf(stateNo, "%03d", this->dumpStateNoTIO);

    char fileName[20];
    strcpy(fileName, "pom");
    strcat(fileName, stateNo);
    strcat(fileName, ".h5");

    std::cout << "Outputting to " << fileName << std::endl;

    TIO_t stat;

    // Create a new file. We're doing one file per state, as it's just easier
    TIO_File_t fileID;
    stat = TIO_Create(
        fileName, 
        &fileID, 
        TIO_ACC_REPLACE, 
        "pom", 
        "0.1", 
        "N/A", 
        "N/A", 
        TIO_NULL, 
        TIO_NULL, 
        TIO_NULL
    );
    std::cout << "Stat: " << stat << std::endl;

    // Create a state
    char stateName[20];
    strcpy(stateName, "state_");
    strcat(stateName, stateNo);

    TIO_Object_t stateID;
    stat = TIO_Create_State(fileID, stateName, &stateID, this->dumpStateNoTIO, time, "s");
    std::cout << "Stat: " << stat << std::endl;

    // Size of mesh, ignoring boundary cells
    int jSize = this->jUpper - this->nghosts;
    int iSize = this->iUpper - this->nghosts;

    // Create the mesh
    TIO_Object_t meshID;
    stat = TIO_Create_Mesh( 
        fileID,
        stateID,
        "mesh",
        &meshID,
        TIO_MESH_QUAD_COLINEAR,
        TIO_COORD_CARTESIAN,
        TIO_FALSE,
        "",
        0,
        TIO_INT,
        TIO_DOUBLE,
        TIO_2D,
        iSize+1,
        jSize+1,
        0,
        0,
        1,
        "cm",
        "cm",
        "",
        "X",
        "Y",
        "" 
    );
    std::cout << "Stat: " << stat << std::endl;

    // Create node centred Coordinates
    int nodeDims[2] = {iSize + 1, jSize + 1};
    std::vector<double> x = new double[nodeDims[0]];
    std::vector<double> y = new double[nodeDims[1]];

    for (int i=0; i<nodeDims[0]; i++) {
        x[i] = (i)*this->dx;
    }

    for (int j=0; j<nodeDims[1]; j++) {
        y[j] = (j)*this->dy;
    }

    // Set the chunk...
    TIO_Set_Quad_Chunk(
        fileID,
        meshID,
        0,
        TIO_2D,
        0,
        iSize,
        0,
        jSize,
        0,
        0,
        0,
        0 
    );

    // Write the mesh
    stat = TIO_Write_QuadMesh_All( 
        fileID,
        meshID,
        TIO_DOUBLE,
        x,
        y,
        NULL 
    );
    std::cout << "Stat: " << stat << std::endl;

    // We need a material, so just create an array of ones...
    int *mat = new int[iSize*jSize];
    for (int i = 0; i < iSize*jSize; i++) {
        mat[i] = 1;
    }

    TIO_Object_t materialID;
    TIO_Create_Material(
        fileID,
        meshID,
        "material",
        &materialID,
        TIO_INT,
        1,
        0,
        TIO_FALSE,
        TIO_INT,
        TIO_INT,
        TIO_DOUBLE
    );

    TIO_Write_QuadMaterial_Chunk( 
        fileID,
        materialID,
        0,
        TIO_XFER_NULL,
        TIO_INT,
        mat,
        TIO_INT,
        TIO_INT,
        TIO_DOUBLE,
        NULL,
        NULL,
        NULL,
        NULL
    );

    // Quants
    TIO_Object_t quantID;

    // Set quant storage
    std::vector<double> data = new double[iSize*jSize];

    // Create a quant for Density
    TIO_Create_Quant( 
        fileID,
        meshID,
        "density",
        &quantID,
        TIO_DOUBLE,
        TIO_CENTRE_CELL,
        0,
        TIO_FALSE,
        "g/cc"
    );

    // Write density into storage
    int index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            data[index++] = _PGET(this, rho, j, i);
        }
    }

    // Write quant
    TIO_Write_QuadQuant_Chunk(
        fileID,
        quantID,
        0,
        TIO_XFER_NULL,
        TIO_DOUBLE,
        data,
        NULL
    );

    TIO_Close_Quant(
        fileID,
        quantID 
    );


    // Create a quant for momU
    TIO_Create_Quant( 
        fileID,
        meshID,
        "momentum - U",
        &quantID,
        TIO_DOUBLE,
        TIO_CENTRE_CELL,
        0,
        TIO_FALSE,
        "gcm/s"
    );

    // Write density into storage
    index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            data[index++] = _PGET(this, momU, j, i);
        }
    }

    // Write quant
    TIO_Write_QuadQuant_Chunk(
        fileID,
        quantID,
        0,
        TIO_XFER_NULL,
        TIO_DOUBLE,
        data,
        NULL
    );

    TIO_Close_Quant(
        fileID,
        quantID 
    );


    // Create a quant for momV
    TIO_Create_Quant( 
        fileID,
        meshID,
        "momentum - V",
        &quantID,
        TIO_DOUBLE,
        TIO_CENTRE_CELL,
        0,
        TIO_FALSE,
        "gcm/s"
    );

    // Write density into storage
    index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            data[index++] = _PGET(this, momV, j, i);
        }
    }

    // Write quant
    TIO_Write_QuadQuant_Chunk(
        fileID,
        quantID,
        0,
        TIO_XFER_NULL,
        TIO_DOUBLE,
        data,
        NULL
    );

    TIO_Close_Quant(
        fileID,
        quantID 
    );

    // Create a quant for pressure
    TIO_Create_Quant( 
        fileID,
        meshID,
        "pressure",
        &quantID,
        TIO_DOUBLE,
        TIO_CENTRE_CELL,
        0,
        TIO_FALSE,
        ""
    );

    // Write density into storage
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

    // Write quant
    TIO_Write_QuadQuant_Chunk(
        fileID,
        quantID,
        0,
        TIO_XFER_NULL,
        TIO_DOUBLE,
        data,
        NULL
    );

    TIO_Close_Quant(
        fileID,
        quantID 
    );

    TIO_Close_Mesh(
        fileID,
        meshID 
    );

    TIO_Close_State(
        fileID,
        stateID 
    );

    TIO_Close(fileID);
}
#endif

#ifdef HasSILO
void Mesh2D::dumpToSILO(double time, int step) {
    char stateNo[4];
    this->dumpStateNoSILO++;
    sprintf(stateNo, "%03d", this->dumpStateNoSILO);

    char fileName[20];
    strcpy(fileName, "pom");
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
    std::vector<double> x = new double[nodeDims[0]];
    std::vector<double> y = new double[nodeDims[1]];

    std::vector<double> coordinates[2];
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
    std::vector<double> data = new double[cellDims[0]*cellDims[1]];
    std::vector<double> data2 = new double[cellDims[0]*cellDims[1]];


    // Write density
    int index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            data[index++] = _PGET(this, rho, j, i);
        }
    }
    DBPutQuadvar1(file, "Density", "mesh1", data, cellDims, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

    // Write momenta
    index = 0;
    for (int j=2; j<this->jUpper; j++) {
        for (int i=2; i<this->iUpper; i++) {
            data[index] = _PGET(this, momU, j, i);
            data2[index++] = _PGET(this, momV, j, i);
        }
    }
    DBPutQuadvar1(file, "MomU", "mesh1", data, cellDims, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
    DBPutQuadvar1(file, "MomV", "mesh1", data2, cellDims, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

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
    delete[] data2;
}
#endif

void Mesh2D::Kill() {
    // delete[] rho;
    // delete[] momU;
    // delete[] momV;
    // delete[] E;
}
