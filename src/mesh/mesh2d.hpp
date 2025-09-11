#ifndef MESH2D_H
#define MESH2D_H
// #include "sweep.hpp"/**/

#include <vector>

#define GET(M, D, J, I) M.D[(J)*M.niGhosts + I]
#define _LGET(D, J, I) D[(J)*niGhosts + I]
#define _PGET(M, D, J, I) M->D[(J)*M->niGhosts + I]

class Mesh2D {
public:
    const int nghosts = 2;
    int niGhosts;
    int njGhosts;
    int iUpper;
    int jUpper;

    int dumpStateNoSILO = 0;
    int dumpStateNoTIO = 0;

    double gamma;
    double dtmax;
    double dx, dy;
    double cfl;

    std::vector<double> meshBoundaryUD;
    std::vector<double> meshBoundaryLR;

    std::vector<double> x, y;
    // std::vector<double> x;
    // std::vector<double> y;

    std::vector<double> rho;
    std::vector<double> momU;
    std::vector<double> momV;
    std::vector<double> E;

    // Constructor to just produce a Sod mesh
    Mesh2D(int, int, int);

    void Kill();
    void setBoundaries();

#ifdef HasSILO
    void dumpToSILO(double, int);
#endif
#ifdef HasTIO
    void dumpToTIO(double, int);
#endif
};
#endif
