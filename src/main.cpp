#include "mesh/mesh2d.hpp"
#include "hydro/hydro.hpp"

#include <iostream>
#include <omp.h>

#include <fenv.h>

void copy(double *rhoOut, double *momUOut, double *momVOut, double *EOut,
          double *rhoIn, double *momUIn, double *momVIn, double *EIn,
          int nCells)
{
    #pragma omp parallel for
    for (int n=0; n<nCells; n++) {
        rhoOut[n] = rhoIn[n];
        EOut[n] = EIn[n];
        momUOut[n] = momUIn[n];
        momVOut[n] = momVIn[n];
    }
}


void evolve(double *rhoIn, double *momUIn, double *momVIn, double *EIn,
          double *rhoOut, double *momUOut, double *momVOut, double *EOut,
          double dt, double dx, double dy, double gamma,
          int nCells, int niGhosts, int ni)
{
        #pragma omp parallel for
        for (int n =0; n<nCells; n++) {
            int j = n/ni;
            int i = n - (ni*j);

            Hydro::MUSCLHancock2D(
                rhoIn, EIn, momUIn, momVIn,
                i+2, j+2, 0, niGhosts,
                gamma, dt, dx, dy,
                rhoOut, EOut, momUOut, momVOut
            );
        }
}


int main(int argc, char* argv[]) {
#ifdef DEBUG
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif

    // Mesh sizes - hardcoded for now
    const int ni = 600;
    const int nj = 240;
    const int problem = 4;

    const double dtOut = 2.0;
    const double tEnd = 80.0;


    const int nCells = ni*nj;

    // MPI environment
    int nprocs, myrank, error;
    nprocs = 1;
    myrank = 0;

    // 100x100 mesh using bubbles set-up
    Mesh2D mesh(ni, nj, problem);

    if (myrank == 0) {
        mesh.dumpToSILO(0.0, 0);
        mesh.dumpToTIO(0.0, 0);
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

        // Evolve the solution to next time step
        evolve(
            mesh.rho, mesh.momU, mesh.momV, mesh.E,
            rhoNew, momUNew, momVNew, ENew,
            dt, mesh.dx, mesh.dy, mesh.gamma,
            nCells, mesh.niGhosts, ni
        );

        // Copy into main data arrays
        copy(
            mesh.rho, mesh.momU, mesh.momV, mesh.E,
            rhoNew, momUNew, momVNew, ENew,
            mesh.niGhosts*mesh.njGhosts
        );

        mesh.setBoundaries();

        t += dt;
        if (t >= outNext) {
            outNext += dtOut;
            if (myrank == 0) {
                mesh.dumpToSILO(t, step);
                mesh.dumpToTIO(t, step);
            }
        }
        if (t > tEnd || step > 100000) break;
    }

    mesh.Kill();

    return 0;
}
