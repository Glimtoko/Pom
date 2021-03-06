#ifndef FLUXH
#define FLUXH

// #include "mesh2d.hpp"

namespace Hydro {
    struct Flux {
        double rho;
        double momU;
        double momV;
        double E;
    };

    Flux getFluxHLLC(
        double uL, double vL, double rhoL, double pL,
        double uR, double vR, double rhoR, double pR,
        double gamma
    );

    void MUSCLHancock1D(
        double* rho, double* E, double* momN, double* momT,
        int ni, int iUpper, double gamma, double dt, double dx
    );

    void MUSCLHancock2D(
        double* rhoOld, double* EOld, double* momUOld, double* momVOld,
        int i, int j, int k, int niGhosts,
        double gamma, double dt, double dx, double dy,
        double* rhoNew, double* ENew, double* momUNew, double* momVNew
    );

    double getTimestep(
        double *rho, double *momU, double *momV, double *E,
        double dx, double dy,
        int nCells,
        double gamma, double cfl, double dtmax
    );
}

#endif
