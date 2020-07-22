#include "hydro.hpp"

#include <math.h>
#include <iostream>

double getLimiter(double di1, double di2, double omega);
double getSlopeX(double* U, int i, int j, int niGhosts, double omega);
double getSlopeY(double* U, int i, int j, int niGhosts, double omega);

double calcFluxRho(double rho, double u);
double calcFluxMom(double rho, double u, double v, double p);
double calcFluxE(double u, double E, double p);

#define GET(D, I, J) D[(J)*niGhosts + I]

void Hydro::MUSCLHancock2D(
    double* rhoOld, double* EOld, double* momUOld, double* momVOld,
    int iIndex, int jIndex, int kIndex, int niGhosts,
    double gamma, double dt, double dx, double dy,
    double& rhoNew, double& ENew, double& momUNew, double& momVNew
) {

    double omega = 0.0;

    // Data reconstruction in X
    double rhoL[3];
    double rhoR[3];
    double momUL[3];
    double momUR[3];
    double momVL[3];
    double momVR[3];
    double EL[3];
    double ER[3];

    int n = 0;
    for (int i=iIndex-1; i<iIndex+2; i++) {
        double di = 0.5*getSlopeX(rhoOld, i, jIndex, niGhosts, omega);
        rhoL[n] = GET(rhoOld, i, jIndex) - di;
        rhoR[n] = GET(rhoOld, i, jIndex) + di;

        di = 0.5*getSlopeX(momUOld, i, jIndex, niGhosts, omega);
        momUL[n] = GET(momUOld, i, jIndex) - di;
        momUR[n] = GET(momUOld, i, jIndex) + di;

        di = 0.5*getSlopeX(momVOld, i, jIndex, niGhosts, omega);
        momVL[n] = GET(momVOld, i, jIndex) - di;
        momVR[n] = GET(momVOld, i, jIndex) + di;

        di = 0.5*getSlopeX(EOld, i, jIndex, niGhosts, omega);
        EL[n] = GET(EOld, i, jIndex) - di;
        ER[n] = GET(EOld, i, jIndex) + di;

        n++;
    }

    // Data reconstruction in Y
    double rhoD[3];
    double rhoU[3];
    double momUD[3];
    double momUU[3];
    double momVD[3];
    double momVU[3];
    double ED[3];
    double EU[3];

    n = 0;
    for (int j=jIndex-1; j<jIndex+2; j++) {
        double di = 0.5*getSlopeY(rhoOld, iIndex, j, niGhosts, omega);
        rhoD[n] = GET(rhoOld, iIndex, j) - di;
        rhoU[n] = GET(rhoOld, iIndex, j) + di;

        di = 0.5*getSlopeY(momUOld, iIndex, j, niGhosts, omega);
        momUD[n] = GET(momUOld, iIndex, j) - di;
        momUU[n] = GET(momUOld, iIndex, j) + di;

        di = 0.5*getSlopeY(momVOld, iIndex, j, niGhosts, omega);
        momVD[n] = GET(momVOld, iIndex, j) - di;
        momVU[n] = GET(momVOld, iIndex, j) + di;

        di = 0.5*getSlopeY(EOld, iIndex, j, niGhosts, omega);
        ED[n] = GET(EOld, iIndex, j) - di;
        EU[n] = GET(EOld, iIndex, j) + di;

        n++;
    }

    double fx = 0.5*dt/dx;
    double fy = 0.5*dt/dy;
    for (int i=0; i<3; i++) {
        double uL = momUL[i]/rhoL[i];
        double uR = momUR[i]/rhoR[i];

        double vL = momVL[i]/rhoL[i];
        double vR = momVR[i]/rhoR[i];

        double pL = (gamma - 1.0)*(
            EL[i] - 0.5*rhoL[i]*uL*uL
                  - 0.5*rhoL[i]*vL*vL
        );

        double pR = (gamma - 1.0)*(
            ER[i] - 0.5*rhoR[i]*uR*uR
                  - 0.5*rhoR[i]*vR*vR
        );


        double dFx_rho = fx*(calcFluxRho(rhoL[i], uL) -
                             calcFluxRho(rhoR[i], uR));

        double dFx_momN = fx*(calcFluxMom(rhoL[i], uL, uL, pL) -
                              calcFluxMom(rhoR[i], uR, uR, pR));

        double dFx_momT = fx*(calcFluxMom(rhoL[i], uL, vL, 0.0) -
                              calcFluxMom(rhoR[i], uR, vR, 0.0));

        double dFx_E = fx*(calcFluxE(uL, EL[i], pL) -
                           calcFluxE(uR, ER[i], pR));


        double uD = momUD[i]/rhoD[i];
        double uU = momUU[i]/rhoU[i];

        double vD = momVD[i]/rhoD[i];
        double vU = momVU[i]/rhoU[i];

        double pD = (gamma - 1.0)*(
            ED[i] - 0.5*rhoD[i]*uD*uD
                  - 0.5*rhoD[i]*vD*vD
        );

        double pU = (gamma - 1.0)*(
            EU[i] - 0.5*rhoU[i]*uU*uU
                  - 0.5*rhoU[i]*vU*vU
        );


        double dFy_rho = fy*(calcFluxRho(rhoD[i], vD) -
                             calcFluxRho(rhoU[i], vU));

        double dFy_momN = fy*(calcFluxMom(rhoD[i], vD, vD, pD) -
                              calcFluxMom(rhoU[i], vU, vU, pU));

        double dFy_momT = fy*(calcFluxMom(rhoD[i], uD, vD, 0.0) -
                              calcFluxMom(rhoU[i], uU, vU, 0.0));

        double dFy_E = fy*(calcFluxE(vD, ED[i], pD) -
                           calcFluxE(vU, EU[i], pU));


        rhoL[i] += dFx_rho + dFy_rho;
        rhoR[i] += dFx_rho + dFy_rho;
        momUL[i] += dFx_momN + dFy_momT;
        momUR[i] += dFx_momN + dFy_momT;
        momVL[i] += dFx_momT + dFy_momN;
        momVR[i] += dFx_momT + dFy_momN;
        EL[i] += dFx_E + dFy_E;
        ER[i] += dFx_E + dFy_E;

        rhoD[i] += dFx_rho + dFy_rho;
        rhoU[i] += dFx_rho + dFy_rho;
        momUD[i] += dFx_momN + dFy_momT;
        momUU[i] += dFx_momN + dFy_momT;
        momVD[i] += dFx_momT + dFy_momN;
        momVU[i] += dFx_momT + dFy_momN;
        ED[i] += dFx_E + dFy_E;
        EU[i] += dFx_E + dFy_E;

    }


    Flux fluxL, fluxR, fluxU, fluxD;
    int i = 1;
    double rhoX1 = rhoR[i-1];
    double uX1 = momUR[i-1]/rhoR[i-1];
    double vX1 = momVR[i-1]/rhoR[i-1];
    double pX1 = (gamma - 1.0)*(
        ER[i-1] - 0.5*rhoR[i-1]*uX1*uX1
                - 0.5*rhoR[i-1]*vX1*vX1
    );

    double rhoX2 = rhoL[i];
    double uX2 = momUL[i]/rhoL[i];
    double vX2 = momVL[i]/rhoL[i];
    double pX2 = (gamma - 1.0)*(
        EL[i] - 0.5*rhoL[i]*uX2*uX2
              - 0.5*rhoL[i]*vX2*vX2
    );

    double rhoX3 = rhoR[i];
    double uX3 = momUR[i]/rhoR[i];
    double vX3 = momVR[i]/rhoR[i];
    double pX3 = (gamma - 1.0)*(
        ER[i] - 0.5*rhoR[i]*uX3*uX3
              - 0.5*rhoR[i]*vX3*vX3
    );

    double rhoX4 = rhoL[i+1];
    double uX4 = momUL[i+1]/rhoL[i+1];
    double vX4 = momVL[i+1]/rhoL[i+1];
    double pX4 = (gamma - 1.0)*(
        EL[i+1] - 0.5*rhoL[i+1]*uX4*uX4
                - 0.5*rhoL[i+1]*vX4*vX4
    );

    fluxL = Hydro::getFluxHLLC(
        uX1, vX1, rhoX1, pX1,
        uX2, vX2, rhoX2, pX2,
        gamma);

    fluxR = Hydro::getFluxHLLC(
        uX3, vX3, rhoX3, pX3,
        uX4, vX4, rhoX4, pX4,
        gamma);


    // Y
    double rhoY1 = rhoU[i-1];
    double uY1 = momUU[i-1]/rhoU[i-1];
    double vY1 = momVU[i-1]/rhoU[i-1];
    double pY1 = (gamma - 1.0)*(
        EU[i-1] - 0.5*rhoU[i-1]*uY1*uY1
                - 0.5*rhoU[i-1]*vY1*vY1
    );

    double rhoY2 = rhoD[i];
    double uY2 = momUD[i]/rhoD[i];
    double vY2 = momVD[i]/rhoD[i];
    double pY2 = (gamma - 1.0)*(
        ED[i] - 0.5*rhoD[i]*uY2*uY2
              - 0.5*rhoD[i]*vY2*vY2
    );

    double rhoY3 = rhoU[i];
    double uY3 = momUU[i]/rhoU[i];
    double vY3 = momVU[i]/rhoU[i];
    double pY3 = (gamma - 1.0)*(
        EU[i] - 0.5*rhoU[i]*uY3*uY3
              - 0.5*rhoU[i]*vY3*vY3
    );

    double rhoY4 = rhoD[i+1];
    double uY4 = momUD[i+1]/rhoD[i+1];
    double vY4 = momVD[i+1]/rhoD[i+1];
    double pY4 = (gamma - 1.0)*(
        ED[i+1] - 0.5*rhoD[i+1]*uY4*uY4
                - 0.5*rhoD[i+1]*vY4*vY4
    );

    fluxD = Hydro::getFluxHLLC(
        vY1, uY1, rhoY1, pY1,
        vY2, uY2, rhoY2, pY2,
        gamma);

    fluxU = Hydro::getFluxHLLC(
        vY3, uY3, rhoY3, pY3,
        vY4, uY4, rhoY4, pY4,
        gamma);

    fx = dt/dx;
    fy = dt/dy;

    rhoNew = GET(rhoOld, iIndex, jIndex) + fx*(fluxL.rho - fluxR.rho);// + fy*(fluxD.rho - fluxU.rho);
    momUNew = GET(momUOld, iIndex, jIndex) + fx*(fluxL.momU - fluxR.momU);// + fy*(fluxD.momV - fluxU.momV);
    momVNew = GET(momVOld, iIndex, jIndex) + fx*(fluxL.momV - fluxR.momV);// + fy*(fluxD.momU - fluxU.momU);
    ENew = GET(EOld, iIndex, jIndex) + fx*(fluxL.E - fluxR.E);// + fy*(fluxD.E - fluxU.E);

}

double getSlopeX(double* U, int i, int j, int niGhosts, double omega) {
    double di1 = GET(U, i, j) - GET(U, i-1, j);
    double di2 = GET(U, i+1, j) - GET(U, i, j);

    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    diU *= getLimiter(di1, di2, omega);

    return diU;
}

double getSlopeY(double* U, int i, int j, int niGhosts, double omega) {
    double di1 = GET(U, i, j) - GET(U, i, j-1);
    double di2 = GET(U, i, j+1) - GET(U, i, j);

    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    diU *= getLimiter(di1, di2, omega);

    return diU;
}

double getLimiter(double di1, double di2, double omega) {
    double xi;
    // Slope limiter - Van Leer
    if (di2 == 0) {
        xi = 0.0;
    } else {
        double r = di1/di2;
        if (r <= 0.0) {
            xi = 0.0;
        } else {
            double xiR = 2.0/(1.0 - omega + (1 + omega)*r);
            xi = std::min(2*r/(1+r), xiR);
        }
    }
    return xi;
}

inline double calcFluxRho(double rho, double u) {
    return rho*u;
}


inline double calcFluxMom(double rho, double u, double v, double p) {
    return rho*u*v + p;
}


inline double calcFluxE(double u, double E, double p) {
    return u*(E + p);
}
