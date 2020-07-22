#include "mesh/mesh2d.hpp"
#include "hydro/hydro.hpp"

#include <iostream>
#include <omp.h>

#include <fenv.h>

int main(int argc, char* argv[]) {
    //feenableexcept(FE_INVALID | FE_OVERFLOW);

    // Mesh sizes - hardcoded for now
    const int ni = 300;
    const int nj = 100;
    const int problem = 3;

    const double dtOut = 1.0;
    const double tEnd = 20.0;

    // MPI environment
    int nprocs, myrank, error;
    nprocs = 1;
    myrank = 0;

    // 100x100 mesh using spherical Sod set-up
    Mesh2D mesh(ni, nj, problem);

    if (myrank == 0) {
        mesh.dumpToSILO(0.0, 0);
    }

    double t = 0.0;
    double outNext = t + dtOut;
    int step = 0;

    double *rhoNew = new double[mesh.niGhosts*mesh.njGhosts];
    double *momUNew = new double[mesh.niGhosts*mesh.njGhosts];
    double *momVNew = new double[mesh.niGhosts*mesh.njGhosts];
    double *ENew = new double[mesh.niGhosts*mesh.njGhosts];

    for (int i=0; i<mesh.njGhosts*mesh.niGhosts; i++) {
        rhoNew[i] = 0.0001;
        momUNew[i] = 0.0001;
        momVNew[i] = 0.0001;
        ENew[i] = 0.0001;
    }

    mesh.setBoundaries();

    for( ; ; ) {
        step++;

        double dt = Hydro::getTimestep(
            mesh.rho,  mesh.momU,  mesh.momV,  mesh.E,
            mesh.dx, mesh.dy,
            mesh.niGhosts*mesh.njGhosts,
            mesh.gamma, mesh.cfl, mesh.dtmax
        );

        dt = std::min(dt, outNext - t);

        if (myrank == 0) {
            std::cout << "Step: " << step;
            std::cout << ", time = " << t;
            std::cout << ", dt = " << dt << std::endl;
        }

#pragma omp parallel for collapse(2)
        for (int j = 2; j<mesh.jUpper; j++) {
            for (int i = 2; i<mesh.iUpper; i++) {
                Hydro::MUSCLHancock2D(
                    mesh.rho, mesh.E, mesh.momU, mesh.momV,
                    i, j, 0, mesh.niGhosts,
                    mesh.gamma, dt, mesh.dx, mesh.dy,
                    rhoNew[j*mesh.niGhosts + i],
                    ENew[j*mesh.niGhosts + i],
                    momUNew[j*mesh.niGhosts + i],
                    momVNew[j*mesh.niGhosts + i]
                );
            }
        }

//         for (int i=0; i<mesh.niGhosts*mesh.njGhosts; i++) {
//             mesh.rho[i] = rhoNew[i];
//             mesh.momU[i] = momUNew[i];
//             mesh.momV[i] = momVNew[i];
//             mesh.E[i] = ENew[i];
//         }

        for (int j = 2; j<mesh.jUpper; j++) {
            for (int i = 2; i<mesh.iUpper; i++) {
                int index = j*mesh.niGhosts + i;
                mesh.rho[index] = rhoNew[index];
                mesh.E[index] = ENew[index];
                mesh.momU[index] = momUNew[index];
                mesh.momV[index] = momVNew[index];
            }
        }

//         mesh.setBoundaries();

        t += dt;
        if (t >= outNext) {
            outNext += dtOut;
            if (myrank == 0) mesh.dumpToSILO(t, step);
        }
        if (t > tEnd || step > 100000) break;
    }

    mesh.Kill();

    return 0;
}
