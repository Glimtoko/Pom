#ifndef FLUXH
#define FLUXH

// #include "mesh2d.hpp"
#include <vector>

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
        std::vector<double> rho, std::vector<double> E, std::vector<double> momN, std::vector<double> momT,
        int ni, int iUpper, double gamma, double dt, double dx
    );

    void MUSCLHancock2D(
        std::vector<double> rhoOld, std::vector<double> EOld, std::vector<double> momUOld, std::vector<double> momVOld,
        int i, int j, int k, int niGhosts,
        double gamma, double dt, double dx, double dy,
        std::vector<double> rhoNew, std::vector<double> ENew, std::vector<double> momUNew, std::vector<double> momVNew
    );

    double getTimestep(
        std::vector<double> rho, std::vector<double> momU, std::vector<double> momV, std::vector<double> E,
        double dx, double dy,
        int nCells,
        double gamma, double cfl, double dtmax
    );
}

#endif
