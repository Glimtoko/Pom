#include "mesh/mesh2d.hpp"
#include "hydro/hydro.hpp"

#include <iostream>
#include <vector>

#include <omp.h>

#include <fenv.h>

void copy(std::vector<double> rhoOut, std::vector<double> momUOut, std::vector<double> momVOut, std::vector<double> EOut,
          std::vector<double> rhoIn, std::vector<double> momUIn, std::vector<double> momVIn, std::vector<double> EIn,
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


void evolve(std::vector<double> rhoIn, std::vector<double> momUIn, std::vector<double> momVIn, std::vector<double> EIn,
          std::vector<double> rhoOut, std::vector<double> momUOut, std::vector<double> momVOut, std::vector<double> EOut,
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


int main() {
#ifdef DEBUG
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif

    // Mesh sizes - hardcoded for now
    const int ni = 600;
    const int nj = 240;
    const int problem = 4;

    const double dtOut = 2.0;
    const double tEnd = 20.0;
    int maxSteps = 100000;  // Set to a big number


    const int nCells = ni*nj;

    // MPI environment
    int nprocs, myrank, error;
    nprocs = 1;
    myrank = 0;

    // OpenMP
    #pragma omp parallel
    if (omp_get_thread_num() == 0) printf("This is POM running on %d threads\n", omp_get_num_threads());

    printf("Mesh size: %dx%d = %d\n", ni, nj, ni*nj);

    double tStart = omp_get_wtime();
    double tNow;
    double tLast = tStart;

    // 100x100 mesh using bubbles set-up
    Mesh2D mesh(ni, nj, problem);
    printf("Setup Complete");

    if (myrank == 0) {
#ifdef HasTIO
        mesh.dumpToSILO(0.0, 0);
#endif
#ifdef HasSILO
        mesh.dumpToTIO(0.0, 0);
#endif
    }

    double t = 0.0;
    double outNext = t + dtOut;
    int step = 0;

    std::vector<double> rhoNew;
    std::vector<double> momUNew;
    std::vector<double> momVNew;
    std::vector<double> ENew;

    rhoNew.reserve(mesh.niGhosts*mesh.njGhosts);
    momUNew.reserve(mesh.niGhosts*mesh.njGhosts);
    momVNew.reserve(mesh.niGhosts*mesh.njGhosts);
    ENew.reserve(mesh.niGhosts*mesh.njGhosts);

    // std::vector<double> rhoNew = new double[mesh.niGhosts*mesh.njGhosts];
    // std::vector<double> momUNew = new double[mesh.niGhosts*mesh.njGhosts];
    // std::vector<double> momVNew = new double[mesh.niGhosts*mesh.njGhosts];
    // std::vector<double> ENew = new double[mesh.niGhosts*mesh.njGhosts];

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

        tNow = omp_get_wtime();
        if (myrank == 0) {
            printf(
                "Step %6d, time = %f, dt = %f, grind = %f (per-step = %e)\n",
                step, t, dt, tNow - tLast, (tNow-tLast)/(ni*nj)
            );
        }
        tLast = tNow;

        t += dt;
        if (t >= outNext) {
            outNext += dtOut;
            if (myrank == 0) {
#ifdef HasTIO
                mesh.dumpToTIO(t, step);
#endif
#ifdef HasSILO
                mesh.dumpToSILO(t, step);
#endif
            }
        }
        if (t > tEnd || step >= maxSteps) break;
    }

    mesh.Kill();
    tNow = omp_get_wtime();
    std::cout << "Calculation completed in " << tNow - tStart << "s" << std::endl;

    return 0;
}
