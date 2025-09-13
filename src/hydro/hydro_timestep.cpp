#include "hydro.hpp"
#include <math.h>
#include <iostream>

double Hydro::getTimestep(
    std::vector<double>& rho, std::vector<double>& momU, std::vector<double>& momV, std::vector<double>& E,
    double dx, double dy,
    int nCells,
    double gamma, double cfl, double dtmax
) {
    double Sx = 0.0, Sy = 0.0;
    for (int i=0; i<nCells; i++) {
        double u = momU[i]/rho[i];
        double v = momV[i]/rho[i];
        double p = (gamma - 1.0)*(E[i] - 0.5*rho[i]*u*u - 0.5*rho[i]*u*u);
        double a = sqrt((gamma*p)/rho[i]);

        Sx = std::max(Sx, a + abs(u));
        Sy = std::max(Sy, a + abs(v));
    }

    return std::min(dtmax, std::min(cfl*dx/Sx, cfl*dy/Sy));
}

