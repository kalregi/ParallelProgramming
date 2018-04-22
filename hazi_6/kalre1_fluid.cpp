#include <chrono>
#include <iostream>
#include <math.h>
#include <numeric>
#include <sstream>
#include "mpi_navier.h"

// Global constants in the equations are
double rkold[3];
int nx1;
double rinv8;
double rinv9;
double Minf;
double rinv1;
double rinv4;
double rinv5;
double Pr;
double rinv12;
double rinv13;
double rinv10;
double rinv11;
double rknew[3];
double rc6;
double rc7;
double rc0;
double rc2;
double rc3;
int nx0;
double deltai1;
double deltai0;
double Re;
double deltat;
double gama;
int itercount;
double bandwidth;
double executionCount;

std::stringstream result;

// main program start
int main(int argc, char **argv) {
    
    // Initialising global constants
    nx0 = 257;
    if (argc > 1)
        nx0 = atoi(argv[1]);
    nx1 = 257;
    if (argc > 2)
        nx1 = atoi(argv[2]);
    
    mpi_setup(argc, argv, &nx0, &nx1);
    
    gama = 1.40000000000000;
    Pr = 0.710000000000000;
    Re = 1600;
    deltat = 0.000846250000000000;
    Minf = 0.100000000000000;
    rc6 = 5.0 / 2.0;
    rc7 = 4.0 / 3.0;
    rc0 = 1.0 / 2.0;
    rc2 = 1.0 / 12.0;
    rc3 = 2.0 / 3.0;
    rkold[0] = 1.0 / 4.0;
    rkold[1] = 3.0 / 20.0;
    rkold[2] = 3.0 / 5.0;
    rknew[0] = 2.0 / 3.0;
    rknew[1] = 5.0 / 12.0;
    rknew[2] = 3.0 / 5.0;
    rinv12 = pow(Minf, -2);
    rinv13 = 1.0 / (gama * pow(Minf, 2));
    rinv10 = 1.0 / Pr;
    rinv11 = 1.0 / (gama - 1);
    rinv9 = 1.0 / Re;
    deltai1 = (1.0 / (nx1 - 1.0)) * M_PI;
    deltai0 = (1.0 / (nx0 - 1.0)) * M_PI;
    rinv8 = pow(deltai1, -2);
    rinv1 = 1.0 / deltai1;
    rinv4 = 1.0 / deltai0;
    rinv5 = pow(deltai0, -2);
    itercount = 20;
    
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    
    std::cout << "Running on a " << nx0 + 4 << "x" << nx1 + 4 << " mesh for "
    << itercount << " iterations\n";
    
    // Allocating mesh
    double *rho = new double[(nx0 + 4) * (nx1 + 4)];
    double *rhou0 = new double[(nx0 + 4) * (nx1 + 4)];
    double *rhou1 = new double[(nx0 + 4) * (nx1 + 4)];
    double *rhoE = new double[(nx0 + 4) * (nx1 + 4)];
    double *rho_old = new double[(nx0 + 4) * (nx1 + 4)];
    double *rhou0_old = new double[(nx0 + 4) * (nx1 + 4)];
    double *rhou1_old = new double[(nx0 + 4) * (nx1 + 4)];
    double *rhoE_old = new double[(nx0 + 4) * (nx1 + 4)];
    double *T = new double[(nx0 + 4) * (nx1 + 4)];
    double *u0 = new double[(nx0 + 4) * (nx1 + 4)];
    double *u1 = new double[(nx0 + 4) * (nx1 + 4)];
    double *p = new double[(nx0 + 4) * (nx1 + 4)];
    double *wk0 = new double[(nx0 + 4) * (nx1 + 4)];
    double *wk1 = new double[(nx0 + 4) * (nx1 + 4)];
    double *wk2 = new double[(nx0 + 4) * (nx1 + 4)];
    double *wk3 = new double[(nx0 + 4) * (nx1 + 4)];
    
    result << "Name" << "\t" << "Elapsed time (s)" << "\t" << "Execution count" << "\t" << "Achieved bandwidth (GB/s)" << "\n";
    
    //Fill up dat_ptrs
    dat_ptrs[0] = rho;
    dat_ptrs[1] = rhou0;
    dat_ptrs[2] = rhou1;
    dat_ptrs[3] = rhoE;
    dat_ptrs[4] = T;
    dat_ptrs[5] = u0;
    dat_ptrs[6] = u1;
    dat_ptrs[7] = p;
    
    // Initialisation
    // writing dataset rho with (i,j) access
    // writing dataset rhou0 with (i,j) access
    // writing dataset rhou1 with (i,j) access
    // wirting dataset rhoE with (i,j) access
    start = std::chrono::system_clock::now();
    for (int j = 0; j < nx1 + 4; j++) {
        for (int i = 0; i < nx0 + 4; i++) {
            double x = deltai0 * (i - 2);
            double y = deltai1 * (j - 2);
            double u = sin(x) * cos(y);
            double v = -cos(x) * sin(y);
            double p = 1.0 * rinv13 + 0.25 * (sin(2.0 * x) + sin(2.0 * y));
            double r = gama * pow(Minf, 2) * p;
            rho[(j + 0) * (nx0 + 4) + (i + 0)] = r;
            rhou0[(j + 0) * (nx0 + 4) + (i + 0)] = r * u;
            rhou1[(j + 0) * (nx0 + 4) + (i + 0)] = r * v;
            rhoE[(j + 0) * (nx0 + 4) + (i + 0)] =
            rinv11 * p + 0.5 * r * (pow(u, 2) + pow(v, 2));
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    executionCount = (nx1 + 4)*(nx0 + 4);
    bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
    
    result << "Initialisation" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
    
    
    // Apply boundary conditions
    // Left
    // writing dataset rho with (i-1,j), (i-2,j) access
    // writing dataset rhou0 with (j-1,j), (i-2,j) access
    // writing dataset rhou1 with (i-1,j), (i-2,j) access
    // writing dataset rhoE with (i-1,j), (i-2,j) access
    // reading dataset rho with (i+1,j), (i+2,j) access
    // reading dataset rhou0 with (i+1,j), (i+2,j) access
    // reading dataset rhou1 with (i+1,j), (i+2,j) access
    // reading dataset rhoE with (i+1,j), (i+2,j) access
    start = std::chrono::system_clock::now();
    if(prev_x == MPI_PROC_NULL) {
        for (int j = 0; j < nx1 + 4; j++) {
            for (int i = 2; i < 3; i++) {
                rho[(j + 0) * (nx0 + 4) + (i - 1)] = rho[(j + 0) * (nx0 + 4) + (i + 1)];
                rho[(j + 0) * (nx0 + 4) + (i - 2)] = rho[(j + 0) * (nx0 + 4) + (i + 2)];
                rhou0[(j + 0) * (nx0 + 4) + (i - 1)] =
                rhou0[(j + 0) * (nx0 + 4) + (i + 1)];
                rhou0[(j + 0) * (nx0 + 4) + (i - 2)] =
                rhou0[(j + 0) * (nx0 + 4) + (i + 2)];
                rhou1[(j + 0) * (nx0 + 4) + (i - 1)] =
                rhou1[(j + 0) * (nx0 + 4) + (i + 1)];
                rhou1[(j + 0) * (nx0 + 4) + (i - 2)] =
                rhou1[(j + 0) * (nx0 + 4) + (i + 2)];
                rhoE[(j + 0) * (nx0 + 4) + (i - 1)] = rhoE[(j + 0) * (nx0 + 4) + (i + 1)];
                rhoE[(j + 0) * (nx0 + 4) + (i - 2)] = rhoE[(j + 0) * (nx0 + 4) + (i + 2)];
            }
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    executionCount = (nx1 + 4)*(3-2);
    bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
    
    result << "ABC Left" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
    
    // Right
    // writing dataset rho with (i+1,j), (i+2,j) access
    // writing dataset rhou0 with (j+1,j), (i+2,j) access
    // writing dataset rhou1 with (i+1,j), (i+2,j) access
    // writing dataset rhoE with (i+1,j), (i+2,j) access
    // reading dataset rho with (i-1,j), (i-2,j) access
    // reading dataset rhou0 with (i-1,j), (i-2,j) access
    // reading dataset rhou1 with (i-1,j), (i-2,j) access
    // reading dataset rhoE with (i-1,j), (i-2,j) access
    start = std::chrono::system_clock::now();
    if (next_x == MPI_PROC_NULL) {
        for (int j = 0; j < nx1 + 4; j++) {
            for (int i = nx0 + 1; i < nx0 + 2; i++) {
                rho[(j + 0) * (nx0 + 4) + (i + 1)] = rho[(j + 0) * (nx0 + 4) + (i - 1)];
                rho[(j + 0) * (nx0 + 4) + (i + 2)] = rho[(j + 0) * (nx0 + 4) + (i - 2)];
                rhou0[(j + 0) * (nx0 + 4) + (i + 1)] =
                rhou0[(j + 0) * (nx0 + 4) + (i - 1)];
                rhou0[(j + 0) * (nx0 + 4) + (i + 2)] =
                rhou0[(j + 0) * (nx0 + 4) + (i - 2)];
                rhou1[(j + 0) * (nx0 + 4) + (i + 1)] =
                rhou1[(j + 0) * (nx0 + 4) + (i - 1)];
                rhou1[(j + 0) * (nx0 + 4) + (i + 2)] =
                rhou1[(j + 0) * (nx0 + 4) + (i - 2)];
                rhoE[(j + 0) * (nx0 + 4) + (i + 1)] = rhoE[(j + 0) * (nx0 + 4) + (i - 1)];
                rhoE[(j + 0) * (nx0 + 4) + (i + 2)] = rhoE[(j + 0) * (nx0 + 4) + (i - 2)];
            }
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    executionCount = (nx1 + 4)*(nx0 + 2-(nx0 + 1));
    bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
    
    result << "ABC Right" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
    
    
    // Top
    // writing dataset rho with (i,j-1), (i,j-2) access
    // writing dataset rhou0 with (i,j-1), (i,j-2) access
    // writing dataset rhou1 with (i,j-1), (i,j-2) access
    // writing dataset rhoE with (i,j-1), (i,j-2) access
    // reading dataset rho with (i,j+1), (i,j+2) access
    // reading dataset rhou0 with (i,j+1), (i,j+2) access
    // reading dataset rhou1 with (i,j+1), (i,j+2) access
    // reading dataset rhoE with (i,j+1), (i,j+2) access
    start = std::chrono::system_clock::now();
    if (prev_y == MPI_PROC_NULL) {
        for (int j = 2; j < 3; j++) {
            for (int i = 0; i < nx0 + 4; i++) {
                rho[(j - 1) * (nx0 + 4) + (i + 0)] = rho[(j + 1) * (nx0 + 4) + (i + 0)];
                rho[(j - 2) * (nx0 + 4) + (i + 0)] = rho[(j + 2) * (nx0 + 4) + (i + 0)];
                rhou0[(j - 1) * (nx0 + 4) + (i + 0)] =
                rhou0[(j + 1) * (nx0 + 4) + (i + 0)];
                rhou0[(j - 2) * (nx0 + 4) + (i + 0)] =
                rhou0[(j + 2) * (nx0 + 4) + (i + 0)];
                rhou1[(j - 1) * (nx0 + 4) + (i + 0)] =
                rhou1[(j + 1) * (nx0 + 4) + (i + 0)];
                rhou1[(j - 2) * (nx0 + 4) + (i + 0)] =
                rhou1[(j + 2) * (nx0 + 4) + (i + 0)];
                rhoE[(j - 1) * (nx0 + 4) + (i + 0)] = rhoE[(j + 1) * (nx0 + 4) + (i + 0)];
                rhoE[(j - 2) * (nx0 + 4) + (i + 0)] = rhoE[(j + 2) * (nx0 + 4) + (i + 0)];
            }
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    executionCount = (3-2)*(nx0 + 4);
    bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
    
    result << "ABC Top" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
    
    
    // Bottom
    // writing dataset rho with (i,j+1), (i,j+2) access
    // writing dataset rhou0 with (i,j+1), (i,j+2) access
    // writing dataset rhou1 with (i,j+1), (i,j+2) access
    // writing dataset rhoE with (i,j+1), (i,j+2) access
    // reading dataset rho with (i,j-1), (i,j-2) access
    // reading dataset rhou0 with (i,j-1), (i,j-2) access
    // reading dataset rhou1 with (i,j-1), (i,j-2) access
    // reading dataset rhoE with (i,j-1), (i,j-2) access
    start = std::chrono::system_clock::now();
    if (next_y == MPI_PROC_NULL) {
        for (int j = nx1 + 1; j < nx1 + 2; j++) {
            for (int i = 0; i < nx0 + 4; i++) {
                rho[(j + 1) * (nx0 + 4) + (i + 0)] = rho[(j - 1) * (nx0 + 4) + (i + 0)];
                rho[(j + 2) * (nx0 + 4) + (i + 0)] = rho[(j - 2) * (nx0 + 4) + (i + 0)];
                rhou0[(j + 1) * (nx0 + 4) + (i + 0)] =
                rhou0[(j - 1) * (nx0 + 4) + (i + 0)];
                rhou0[(j + 2) * (nx0 + 4) + (i + 0)] =
                rhou0[(j - 2) * (nx0 + 4) + (i + 0)];
                rhou1[(j + 1) * (nx0 + 4) + (i + 0)] =
                rhou1[(j - 1) * (nx0 + 4) + (i + 0)];
                rhou1[(j + 2) * (nx0 + 4) + (i + 0)] =
                rhou1[(j - 2) * (nx0 + 4) + (i + 0)];
                rhoE[(j + 1) * (nx0 + 4) + (i + 0)] = rhoE[(j - 1) * (nx0 + 4) + (i + 0)];
                rhoE[(j + 2) * (nx0 + 4) + (i + 0)] = rhoE[(j - 2) * (nx0 + 4) + (i + 0)];
            }
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    executionCount = (nx1 + 2 -(nx1 + 1))*(nx0 + 4);
    bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
    
    result << "ABC Bottom" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
    
    //setting dirty the written arrays
    set_dirty(rho);
    set_dirty(rhou0);
    set_dirty(rhou1);
    set_dirty(rhoE);
    
    // Record start time
    auto startMain = std::chrono::high_resolution_clock::now();
    
    // Main time iteration loop
    for (int iteration = 0; iteration < itercount; iteration++) {
        
        //exchange halo
        exchange_halo(nx0, nx1, rho);
        exchange_halo(nx0, nx1, rhou0);
        exchange_halo(nx0, nx1, rhou1);
        exchange_halo(nx0, nx1, rhoE);
        
        // Save equations
        // writing dataset rho_old with (i,j) access
        // writing dataset rhou0_old with (i,j) access
        // writing dataset rhou1_old with (i,j) access
        // writing dataset rhoE_old with (i,j) access
        // reading dataset rho with (i,j) access
        // reading dataset rhou0 with (i,j) access
        // reading dataset rhou1 with (i,j) access
        // reading dataset rhoE with (i,j) access
        if(iteration==0) start = std::chrono::system_clock::now();
        for (int j = 0; j < nx1 + 4; j++) {
            for (int i = 0; i < nx0 + 4; i++) {
                rho_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                rho[(j + 0) * (nx0 + 4) + (i + 0)];
                rhou0_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                rhou0[(j + 0) * (nx0 + 4) + (i + 0)];
                rhou1_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                rhou1[(j + 0) * (nx0 + 4) + (i + 0)];
                rhoE_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                rhoE[(j + 0) * (nx0 + 4) + (i + 0)];
            }
        }
        if(iteration==0) {
            end = std::chrono::system_clock::now();
            elapsed_seconds = end-start;
            executionCount = (nx0 + 4)*(nx0 + 4);
            bandwidth = executionCount * 8 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
            
            result << "Save equations" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
        }
        // Runge-Kutta time-stepper
        for (int stage = 0; stage < 3; stage++) {
            
            // Grouped Formula Evaluation
            // writing dataset T with (i,j) access
            // wrrting dataset p with (i,j) access
            // writing dataset u1 with (i,j) access
            // writing dataset u0 with (i,j) access
            // reading dataset rhou0 with (i,j) access
            // reading dataset rhou1 with (i,j) access
            // reading dataset rho with (i,j) access
            // reading dataset rhoE with (i,j) access
            start = std::chrono::system_clock::now();
            for (int j = 0; j < nx1 + 4; j++) {
                for (int i = 0; i < nx0 + 4; i++) {
                    T[(j + 0) * (nx0 + 4) + (i + 0)] =
                    gama * (gama - 1) *
                    ((-rc0 * pow(rhou0[(j + 0) * (nx0 + 4) + (i + 0)], 2) -
                      rc0 * pow(rhou1[(j + 0) * (nx0 + 4) + (i + 0)], 2)) /
                     rho[(j + 0) * (nx0 + 4) + (i + 0)] +
                     rhoE[(j + 0) * (nx0 + 4) + (i + 0)]) *
                    pow(Minf, 2) / rho[(j + 0) * (nx0 + 4) + (i + 0)];
                    p[(j + 0) * (nx0 + 4) + (i + 0)] =
                    (gama - 1) *
                    ((-rc0 * pow(rhou0[(j + 0) * (nx0 + 4) + (i + 0)], 2) -
                      rc0 * pow(rhou1[(j + 0) * (nx0 + 4) + (i + 0)], 2)) /
                     rho[(j + 0) * (nx0 + 4) + (i + 0)] +
                     rhoE[(j + 0) * (nx0 + 4) + (i + 0)]);
                    u1[(j + 0) * (nx0 + 4) + (i + 0)] =
                    rhou1[(j + 0) * (nx0 + 4) + (i + 0)] /
                    rho[(j + 0) * (nx0 + 4) + (i + 0)];
                    u0[(j + 0) * (nx0 + 4) + (i + 0)] =
                    rhou0[(j + 0) * (nx0 + 4) + (i + 0)] /
                    rho[(j + 0) * (nx0 + 4) + (i + 0)];
                }
            }
            
            //setting dirty the written arrays
            set_dirty(T);
            set_dirty(p);
            set_dirty(u0);
            set_dirty(u1);
            
            //exchange halo
            exchange_halo(nx0, nx1, T);
            exchange_halo(nx0, nx1, p);
            exchange_halo(nx0, nx1, u0);
            exchange_halo(nx0, nx1, u1);

            
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (nx1 + 4)*(nx0 + 4);
                bandwidth = executionCount * 8 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "Grouped Formula Evaluation" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
                
            }
            
            // Residual of equation
            // writing dataset wk0 with (i,j) access
            // writing dataset wk1 with (i,j) access
            // writing dataset wk2 with (i,j) access
            // writing dataset wk3 with (i,j) access
            // reading dataset rhoE with (i,j),(i,j-2),(i,j-1),(i,j+1),(i,j+2),(i-2,j),(i-1,j),(i+1,j),(i+2,j) access
            // reading dataset T with (i,j),(i-2,j),(i-1,j),(i+1,j),(i+2,j),(i,j-2),(i,j-1),(i,j+1),(i,j+2) access
            // reading dataset rhou1 with (i,j),(i,j-2),(i,j-1),(i,j+1),(i,j+2),(i-2,j),(i-1,j),(i+1,j),(i+2,j) access
            // reading dataset rho with (i,j),(i,j-2),(i,j-1),(i,j+1),(i,j+2),(i-2,j),(i-1,j),(i+1,j),(i+2,j) access
            // reading dataset rhou0 with (i,j),(i-2,j),(i-1,j),(i+1,j),(i+2,j),(i,j-2),(i,j-1),(i,j+1),(i,j+2) access
            // reading dataset p with (i,j-2),(i,j-1),(i,j+1),(i,j+2),(i-2,j),(i-1,j),(i+1,j),(i+2,j) access
            // reading dataset u1 with (i,j),(i,j-2),(i,j-1),(i,j+1),(i,j+2),(i-2,j),(i-1,j),(i+1,j),(i+2,j),
            // (i-2,j-2),(i-1,j-2),(i+1,j-2),(i+2,j-2),(i-2,j-1),(i-1,j-1),(i+1,j-1),(i+2,j-1),(i-2,j+1),(i-1,j+1),
            // (i+1,j+1),(i+2,j+1),(i-2,j+2),(i-1,j+2),(i+1,j+2),(i+2,j+2) access
            // reading dataset u0 with (i,j),(i-2,j),(i-1,j),(i+1,j),(i+2,j),(i,j-2),(i,j-1),(i,j+1),(i,j+2),
            // (i-2,j-2),(i-1,j-2),(i+1,j-2),(i+2,j-2),(i-2,j-1),(i-1,j-1),(i+1,j-1),(i+2,j-1),(i-2,j+1),(i-1,j+1),
            // (i+1,j+1),(i+2,j+1),(i-2,j+2),(i-1,j+2),(i+1,j+2),(i+2,j+2) access
            start = std::chrono::system_clock::now();
            for (int j = 2; j < nx1 + 2; j++) {
                for (int i = 2; i < nx0 + 2; i++) {
                    double temp_eval0 =
                    rinv1 * ((rc2)*rhoE[(j - 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rhoE[(j - 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rhoE[(j + 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rhoE[(j + 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval1 =
                    rinv4 * ((rc2)*rhoE[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rhoE[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rhoE[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rhoE[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval2 =
                    rinv4 * ((rc2)*rhoE[(j + 0) * (nx0 + 4) + (i - 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rhoE[(j + 0) * (nx0 + 4) + (i - 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rhoE[(j + 0) * (nx0 + 4) + (i + 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rhoE[(j + 0) * (nx0 + 4) + (i + 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval3 = rinv1 * ((rc2)*u1[(j - 2) * (nx0 + 4) + (i + 0)] -
                                                 rc3 * u1[(j - 1) * (nx0 + 4) + (i + 0)] +
                                                 (rc3)*u1[(j + 1) * (nx0 + 4) + (i + 0)] -
                                                 rc2 * u1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval4 = rinv5 * (-rc6 * T[(j + 0) * (nx0 + 4) + (i + 0)] -
                                                 rc2 * T[(j + 0) * (nx0 + 4) + (i - 2)] +
                                                 (rc7)*T[(j + 0) * (nx0 + 4) + (i - 1)] +
                                                 (rc7)*T[(j + 0) * (nx0 + 4) + (i + 1)] -
                                                 rc2 * T[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval5 =
                    rinv1 * ((rc2)*rhou1[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rhou1[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rhou1[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rhou1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval6 =
                    rinv5 * (-rc6 * u0[(j + 0) * (nx0 + 4) + (i + 0)] -
                             rc2 * u0[(j + 0) * (nx0 + 4) + (i - 2)] +
                             (rc7)*u0[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc7)*u0[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * u0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval7 =
                    rinv1 * ((rc2)*rhou1[(j - 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rhou1[(j - 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rhou1[(j + 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rhou1[(j + 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval8 =
                    rinv1 * ((rc2)*rho[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rho[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rho[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rho[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval9 =
                    rinv4 * ((rc2)*rhou0[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rhou0[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rhou0[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rhou0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval10 =
                    rinv1 * ((rc2)*u0[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * u0[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*u0[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * u0[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval11 = rinv1 * ((rc2)*p[(j - 2) * (nx0 + 4) + (i + 0)] -
                                                  rc3 * p[(j - 1) * (nx0 + 4) + (i + 0)] +
                                                  (rc3)*p[(j + 1) * (nx0 + 4) + (i + 0)] -
                                                  rc2 * p[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval12 =
                    rinv4 * ((rc2)*u1[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * u1[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*u1[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * u1[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval13 =
                    rinv4 * ((rc2)*rhou1[(j + 0) * (nx0 + 4) + (i - 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rhou1[(j + 0) * (nx0 + 4) + (i - 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rhou1[(j + 0) * (nx0 + 4) + (i + 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rhou1[(j + 0) * (nx0 + 4) + (i + 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval14 =
                    rinv4 * ((rc2)*rho[(j + 0) * (nx0 + 4) + (i - 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rho[(j + 0) * (nx0 + 4) + (i - 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rho[(j + 0) * (nx0 + 4) + (i + 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rho[(j + 0) * (nx0 + 4) + (i + 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval15 =
                    rinv4 * ((rc2)*rho[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rho[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rho[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rho[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval16 =
                    rinv1 * ((rc2)*rhou0[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rhou0[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rhou0[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rhou0[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval17 =
                    rinv5 * (-rc6 * u1[(j + 0) * (nx0 + 4) + (i + 0)] -
                             rc2 * u1[(j + 0) * (nx0 + 4) + (i - 2)] +
                             (rc7)*u1[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc7)*u1[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * u1[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval18 =
                    rinv4 * ((rc2)*rhou0[(j + 0) * (nx0 + 4) + (i - 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rhou0[(j + 0) * (nx0 + 4) + (i - 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rhou0[(j + 0) * (nx0 + 4) + (i + 1)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rhou0[(j + 0) * (nx0 + 4) + (i + 2)] *
                             u0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval19 =
                    rinv1 * ((rc2)*rhou0[(j - 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rhou0[(j - 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rhou0[(j + 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rhou0[(j + 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval20 =
                    rinv8 * (-rc6 * u1[(j + 0) * (nx0 + 4) + (i + 0)] -
                             rc2 * u1[(j - 2) * (nx0 + 4) + (i + 0)] +
                             (rc7)*u1[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc7)*u1[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * u1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval21 = rinv4 * ((rc2)*p[(j + 0) * (nx0 + 4) + (i - 2)] -
                                                  rc3 * p[(j + 0) * (nx0 + 4) + (i - 1)] +
                                                  (rc3)*p[(j + 0) * (nx0 + 4) + (i + 1)] -
                                                  rc2 * p[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval22 =
                    rinv8 * (-rc6 * u0[(j + 0) * (nx0 + 4) + (i + 0)] -
                             rc2 * u0[(j - 2) * (nx0 + 4) + (i + 0)] +
                             (rc7)*u0[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc7)*u0[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * u0[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval23 =
                    rinv1 * ((rc2)*rhoE[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rhoE[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rhoE[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rhoE[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval24 =
                    rinv8 * (-rc6 * T[(j + 0) * (nx0 + 4) + (i + 0)] -
                             rc2 * T[(j - 2) * (nx0 + 4) + (i + 0)] +
                             (rc7)*T[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc7)*T[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * T[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval25 =
                    rinv4 * ((rc2)*rhou1[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * rhou1[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*rhou1[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * rhou1[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval26 = rinv4 * ((rc2)*p[(j + 0) * (nx0 + 4) + (i - 2)] *
                                                  u0[(j + 0) * (nx0 + 4) + (i - 2)] -
                                                  rc3 * p[(j + 0) * (nx0 + 4) + (i - 1)] *
                                                  u0[(j + 0) * (nx0 + 4) + (i - 1)] +
                                                  (rc3)*p[(j + 0) * (nx0 + 4) + (i + 1)] *
                                                  u0[(j + 0) * (nx0 + 4) + (i + 1)] -
                                                  rc2 * p[(j + 0) * (nx0 + 4) + (i + 2)] *
                                                  u0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval27 =
                    rinv4 * ((rc2)*u0[(j + 0) * (nx0 + 4) + (i - 2)] -
                             rc3 * u0[(j + 0) * (nx0 + 4) + (i - 1)] +
                             (rc3)*u0[(j + 0) * (nx0 + 4) + (i + 1)] -
                             rc2 * u0[(j + 0) * (nx0 + 4) + (i + 2)]);
                    double temp_eval28 = rinv1 * ((rc2)*p[(j - 2) * (nx0 + 4) + (i + 0)] *
                                                  u1[(j - 2) * (nx0 + 4) + (i + 0)] -
                                                  rc3 * p[(j - 1) * (nx0 + 4) + (i + 0)] *
                                                  u1[(j - 1) * (nx0 + 4) + (i + 0)] +
                                                  (rc3)*p[(j + 1) * (nx0 + 4) + (i + 0)] *
                                                  u1[(j + 1) * (nx0 + 4) + (i + 0)] -
                                                  rc2 * p[(j + 2) * (nx0 + 4) + (i + 0)] *
                                                  u1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval29 =
                    rinv1 * ((rc2)*rho[(j - 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 2) * (nx0 + 4) + (i + 0)] -
                             rc3 * rho[(j - 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j - 1) * (nx0 + 4) + (i + 0)] +
                             (rc3)*rho[(j + 1) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 1) * (nx0 + 4) + (i + 0)] -
                             rc2 * rho[(j + 2) * (nx0 + 4) + (i + 0)] *
                             u1[(j + 2) * (nx0 + 4) + (i + 0)]);
                    double temp_eval30 =
                    rinv1 * ((rc2)*rinv4 * ((rc2)*u0[(j - 2) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u0[(j - 2) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u0[(j - 2) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u0[(j - 2) * (nx0 + 4) + (i + 2)]) -
                             rc3 * rinv4 * ((rc2)*u0[(j - 1) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u0[(j - 1) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u0[(j - 1) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u0[(j - 1) * (nx0 + 4) + (i + 2)]) +
                             (rc3)*rinv4 * ((rc2)*u0[(j + 1) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u0[(j + 1) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u0[(j + 1) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u0[(j + 1) * (nx0 + 4) + (i + 2)]) -
                             rc2 * rinv4 * ((rc2)*u0[(j + 2) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u0[(j + 2) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u0[(j + 2) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u0[(j + 2) * (nx0 + 4) + (i + 2)]));
                    double temp_eval31 =
                    rinv1 * ((rc2)*rinv4 * ((rc2)*u1[(j - 2) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u1[(j - 2) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u1[(j - 2) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u1[(j - 2) * (nx0 + 4) + (i + 2)]) -
                             rc3 * rinv4 * ((rc2)*u1[(j - 1) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u1[(j - 1) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u1[(j - 1) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u1[(j - 1) * (nx0 + 4) + (i + 2)]) +
                             (rc3)*rinv4 * ((rc2)*u1[(j + 1) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u1[(j + 1) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u1[(j + 1) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u1[(j + 1) * (nx0 + 4) + (i + 2)]) -
                             rc2 * rinv4 * ((rc2)*u1[(j + 2) * (nx0 + 4) + (i - 2)] -
                                            rc3 * u1[(j + 2) * (nx0 + 4) + (i - 1)] +
                                            (rc3)*u1[(j + 2) * (nx0 + 4) + (i + 1)] -
                                            rc2 * u1[(j + 2) * (nx0 + 4) + (i + 2)]));
                    wk0[(j + 0) * (nx0 + 4) + (i + 0)] =
                    -0.5 * temp_eval14 -
                    0.5 * temp_eval15 * u0[(j + 0) * (nx0 + 4) + (i + 0)] -
                    0.5 * temp_eval29 -
                    0.5 * temp_eval8 * u1[(j + 0) * (nx0 + 4) + (i + 0)] -
                    0.5 * (temp_eval27 + temp_eval3) *
                    rho[(j + 0) * (nx0 + 4) + (i + 0)];
                    wk1[(j + 0) * (nx0 + 4) + (i + 0)] =
                    -0.5 * temp_eval16 * u1[(j + 0) * (nx0 + 4) + (i + 0)] -
                    0.5 * temp_eval18 - 0.5 * temp_eval19 - temp_eval21 -
                    0.5 * temp_eval9 * u0[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rinv9 * (temp_eval22 + temp_eval31) +
                    rinv9 * (-rc3 * temp_eval31 + (rc7)*temp_eval6) -
                    0.5 * (temp_eval27 + temp_eval3) *
                    rhou0[(j + 0) * (nx0 + 4) + (i + 0)];
                    wk2[(j + 0) * (nx0 + 4) + (i + 0)] =
                    -temp_eval11 - 0.5 * temp_eval13 -
                    0.5 * temp_eval25 * u0[(j + 0) * (nx0 + 4) + (i + 0)] -
                    0.5 * temp_eval5 * u1[(j + 0) * (nx0 + 4) + (i + 0)] -
                    0.5 * temp_eval7 + rinv9 * (temp_eval17 + temp_eval30) +
                    rinv9 * ((rc7)*temp_eval20 - rc3 * temp_eval30) -
                    0.5 * (temp_eval27 + temp_eval3) *
                    rhou1[(j + 0) * (nx0 + 4) + (i + 0)];
                    wk3[(j + 0) * (nx0 + 4) + (i + 0)] =
                    -0.5 * temp_eval0 -
                    0.5 * temp_eval1 * u0[(j + 0) * (nx0 + 4) + (i + 0)] +
                    temp_eval10 * rinv9 * (temp_eval10 + temp_eval12) +
                    temp_eval12 * rinv9 * (temp_eval10 + temp_eval12) -
                    0.5 * temp_eval2 -
                    0.5 * temp_eval23 * u1[(j + 0) * (nx0 + 4) + (i + 0)] +
                    temp_eval24 * rinv10 * rinv11 * rinv12 * rinv9 - temp_eval26 +
                    temp_eval27 * rinv9 * ((rc7)*temp_eval27 - rc3 * temp_eval3) -
                    temp_eval28 +
                    temp_eval3 * rinv9 * (-rc3 * temp_eval27 + (rc7)*temp_eval3) +
                    temp_eval4 * rinv10 * rinv11 * rinv12 * rinv9 +
                    rinv9 * (temp_eval17 + temp_eval30) *
                    u1[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rinv9 * ((rc7)*temp_eval20 - rc3 * temp_eval30) *
                    u1[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rinv9 * (temp_eval22 + temp_eval31) *
                    u0[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rinv9 * (-rc3 * temp_eval31 + (rc7)*temp_eval6) *
                    u0[(j + 0) * (nx0 + 4) + (i + 0)] -
                    0.5 * (temp_eval27 + temp_eval3) *
                    rhoE[(j + 0) * (nx0 + 4) + (i + 0)];
                }
            }
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (nx1)*(nx0);
                bandwidth = executionCount * 12 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "Residual of equation" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
                
            }
            
            // RK new (subloop) update
            // writing dataset rho with(i,j) access
            // writing dataset rhou0 with(i,j) access
            // writing dataset rhou1 with(i,j) access
            // writing dataset rhoE with(i,j) access
            // reading dataset wk0 with (i,j) access
            // reading dataset rho_old with (i,j) access
            // reading dataset wk1 with (i,j) access
            // reading dataset rhou0_old with (i,j) access
            // reading dataset wk2 with (i,j) access
            // reading dataset rhou1_old with (i,j) access
            // reading dataset wk3 with (i,j) access
            // reading dataset rhoE_old with (i,j) access
            start = std::chrono::system_clock::now();
            
            for (int j = 0; j < nx1 + 4; j++) {
                for (int i = 0; i < nx0 + 4; i++) {
                    rho[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rknew[0] * wk0[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rho_old[(j + 0) * (nx0 + 4) + (i + 0)];
                    rhou0[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rknew[0] * wk1[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rhou0_old[(j + 0) * (nx0 + 4) + (i + 0)];
                    rhou1[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rknew[0] * wk2[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rhou1_old[(j + 0) * (nx0 + 4) + (i + 0)];
                    rhoE[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rknew[0] * wk3[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rhoE_old[(j + 0) * (nx0 + 4) + (i + 0)];
                }
            }
            
            //setting dirty the written arrays
            set_dirty(rho);
            set_dirty(rhou0);
            set_dirty(rhou1);
            set_dirty(rhoE);
            
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (nx1)*(nx0);
                bandwidth = executionCount * 8 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "RK new (subloop) update" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
                
            }
            
            // RK old update
            // writing dataset rho_old with(i,j) access
            // writing dataset rhou0_old with(i,j) access
            // writing dataset rhou1_old with(i,j) access
            // writing dataset rhoE_old with(i,j) access
            // reading dataset wk0 with (i,j) access
            // reading dataset rho_old with (i,j) access
            // reading dataset wk1 with (i,j) access
            // reading dataset rhou0_old with (i,j) access
            // reading dataset wk2 with (i,j) access
            // reading dataset rhou1_old with (i,j) access
            // reading dataset wk3 with (i,j) access
            // reading dataset rhoE_old with (i,j) access
            start = std::chrono::system_clock::now();
            for (int j = 0; j < nx1 + 4; j++) {
                for (int i = 0; i < nx0 + 4; i++) {
                    rho_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rkold[0] * wk0[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rho_old[(j + 0) * (nx0 + 4) + (i + 0)];
                    rhou0_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rkold[0] * wk1[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rhou0_old[(j + 0) * (nx0 + 4) + (i + 0)];
                    rhou1_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rkold[0] * wk2[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rhou1_old[(j + 0) * (nx0 + 4) + (i + 0)];
                    rhoE_old[(j + 0) * (nx0 + 4) + (i + 0)] =
                    deltat * rkold[0] * wk3[(j + 0) * (nx0 + 4) + (i + 0)] +
                    rhoE_old[(j + 0) * (nx0 + 4) + (i + 0)];
                }
            }
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (nx1+4)*(nx0+4);
                bandwidth = executionCount * 8 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "RK old update" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
                
            }
            
            // Apply boundary conditions
            
            // Left
            // writing dataset rho with (i-1,j), (i-2,j) access
            // writing dataset rhou0 with (i-1,j), (i-2,j) access
            // writing dataset rhou1 with (i-1,j), (i-2,j) access
            // writing dataset rhoE with (i-1,j), (i-2,j) access
            // reading dataset rho with (i+1,j), (i+2,j) access
            // reading dataset rhou0 with (i+1,j), (i+2,j) access
            // reading dataset rhou1 with (i+1,j), (i+2,j) access
            // reading dataset rhoE with (i+1,j), (i+2,j) access
            start = std::chrono::system_clock::now();
            if (prev_x == MPI_PROC_NULL) {
                for (int j = 0; j < nx1 + 4; j++) {
                    for (int i = 2; i < 3; i++) {
                        rho[(j + 0) * (nx0 + 4) + (i - 1)] =
                        rho[(j + 0) * (nx0 + 4) + (i + 1)];
                        rho[(j + 0) * (nx0 + 4) + (i - 2)] =
                        rho[(j + 0) * (nx0 + 4) + (i + 2)];
                        rhou0[(j + 0) * (nx0 + 4) + (i - 1)] =
                        rhou0[(j + 0) * (nx0 + 4) + (i + 1)];
                        rhou0[(j + 0) * (nx0 + 4) + (i - 2)] =
                        rhou0[(j + 0) * (nx0 + 4) + (i + 2)];
                        rhou1[(j + 0) * (nx0 + 4) + (i - 1)] =
                        rhou1[(j + 0) * (nx0 + 4) + (i + 1)];
                        rhou1[(j + 0) * (nx0 + 4) + (i - 2)] =
                        rhou1[(j + 0) * (nx0 + 4) + (i + 2)];
                        rhoE[(j + 0) * (nx0 + 4) + (i - 1)] =
                        rhoE[(j + 0) * (nx0 + 4) + (i + 1)];
                        rhoE[(j + 0) * (nx0 + 4) + (i - 2)] =
                        rhoE[(j + 0) * (nx0 + 4) + (i + 2)];
                    }
                }
            }
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (nx1+4)*(3-2);
                bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "ABC2 Left" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
            }
            
            // Right
            // writing dataset rho with (i+1,j), (i+2,j) access
            // writing dataset rhou0 with (i+1,j), (i+2,j) access
            // writing dataset rhou1 with (i+1,j), (i+2,j) access
            // writing dataset rhoE with (i+1,j), (i+2,j) access
            // reading dataset rho with (i-1,j), (i-2,j) access
            // reading dataset rhou0 with (i-1,j), (i-2,j) access
            // reading dataset rhou1 with (i-1,j), (i-2,j) access
            // reading dataset rhoE with (i-1,j), (i-2,j) access
            start = std::chrono::system_clock::now();
            if (next_x == MPI_PROC_NULL) {
                for (int j = 0; j < nx1 + 4; j++) {
                    for (int i = nx0 + 1; i < nx0 + 2; i++) {
                        rho[(j + 0) * (nx0 + 4) + (i + 1)] =
                        rho[(j + 0) * (nx0 + 4) + (i - 1)];
                        rho[(j + 0) * (nx0 + 4) + (i + 2)] =
                        rho[(j + 0) * (nx0 + 4) + (i - 2)];
                        rhou0[(j + 0) * (nx0 + 4) + (i + 1)] =
                        rhou0[(j + 0) * (nx0 + 4) + (i - 1)];
                        rhou0[(j + 0) * (nx0 + 4) + (i + 2)] =
                        rhou0[(j + 0) * (nx0 + 4) + (i - 2)];
                        rhou1[(j + 0) * (nx0 + 4) + (i + 1)] =
                        rhou1[(j + 0) * (nx0 + 4) + (i - 1)];
                        rhou1[(j + 0) * (nx0 + 4) + (i + 2)] =
                        rhou1[(j + 0) * (nx0 + 4) + (i - 2)];
                        rhoE[(j + 0) * (nx0 + 4) + (i + 1)] =
                        rhoE[(j + 0) * (nx0 + 4) + (i - 1)];
                        rhoE[(j + 0) * (nx0 + 4) + (i + 2)] =
                        rhoE[(j + 0) * (nx0 + 4) + (i - 2)];
                    }
                }
            }
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (nx1+4)*(nx0+2-(nx0+1));
                bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "ABC2 Right" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
            }
            
            // Top
            // writing dataset rho with (i,j-1), (i,j-2) access
            // writing dataset rhou0 with (i,j-1), (i,j-2) access
            // writing dataset rhou1 with (i,j-1), (i,j-2) access
            // writing dataset rhoE with (i,j-1), (i,j-2) access
            // reading dataset rho with (i,j+1), (i,j+2) access
            // reading dataset rhou0 with (i,j+1), (i,j+2) access
            // reading dataset rhou1 with (i,j+1), (i,j+2) access
            // reading dataset rhoE with (i,j+1), (i,j+2) access
            start = std::chrono::system_clock::now();
            if (prev_y == MPI_PROC_NULL) {
                for (int j = 2; j < 3; j++) {
                    for (int i = 0; i < nx0 + 4; i++) {
                        rho[(j - 1) * (nx0 + 4) + (i + 0)] =
                        rho[(j + 1) * (nx0 + 4) + (i + 0)];
                        rho[(j - 2) * (nx0 + 4) + (i + 0)] =
                        rho[(j + 2) * (nx0 + 4) + (i + 0)];
                        rhou0[(j - 1) * (nx0 + 4) + (i + 0)] =
                        rhou0[(j + 1) * (nx0 + 4) + (i + 0)];
                        rhou0[(j - 2) * (nx0 + 4) + (i + 0)] =
                        rhou0[(j + 2) * (nx0 + 4) + (i + 0)];
                        rhou1[(j - 1) * (nx0 + 4) + (i + 0)] =
                        rhou1[(j + 1) * (nx0 + 4) + (i + 0)];
                        rhou1[(j - 2) * (nx0 + 4) + (i + 0)] =
                        rhou1[(j + 2) * (nx0 + 4) + (i + 0)];
                        rhoE[(j - 1) * (nx0 + 4) + (i + 0)] =
                        rhoE[(j + 1) * (nx0 + 4) + (i + 0)];
                        rhoE[(j - 2) * (nx0 + 4) + (i + 0)] =
                        rhoE[(j + 2) * (nx0 + 4) + (i + 0)];
                    }
                }
            }
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (3-2)*(nx0+4);
                bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "ABC2 Top" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
            }
            
            // Bottom
            // writing dataset rho with (i,j+1), (i,j+2) access
            // writing dataset rhou0 with (i,j+1), (i,j+2) access
            // writing dataset rhou1 with (i,j+1), (i,j+2) access
            // writing dataset rhoE with (i,j+1), (i,j+2) access
            // reading dataset rho with (i,j-1), (i,j-2) access
            // reading dataset rhou0 with (i,j-1), (i,j-2) access
            // reading dataset rhou1 with (i,j-1), (i,j-2) access
            // reading dataset rhoE with (i,j-1), (i,j-2) access
            start = std::chrono::system_clock::now();
            if (next_y == MPI_PROC_NULL) {
                for (int j = nx1 + 1; j < nx1 + 2; j++) {
                    for (int i = 0; i < nx0 + 4; i++) {
                        rho[(j + 1) * (nx0 + 4) + (i + 0)] =
                        rho[(j - 1) * (nx0 + 4) + (i + 0)];
                        rho[(j + 2) * (nx0 + 4) + (i + 0)] =
                        rho[(j - 2) * (nx0 + 4) + (i + 0)];
                        rhou0[(j + 1) * (nx0 + 4) + (i + 0)] =
                        rhou0[(j - 1) * (nx0 + 4) + (i + 0)];
                        rhou0[(j + 2) * (nx0 + 4) + (i + 0)] =
                        rhou0[(j - 2) * (nx0 + 4) + (i + 0)];
                        rhou1[(j + 1) * (nx0 + 4) + (i + 0)] =
                        rhou1[(j - 1) * (nx0 + 4) + (i + 0)];
                        rhou1[(j + 2) * (nx0 + 4) + (i + 0)] =
                        rhou1[(j - 2) * (nx0 + 4) + (i + 0)];
                        rhoE[(j + 1) * (nx0 + 4) + (i + 0)] =
                        rhoE[(j - 1) * (nx0 + 4) + (i + 0)];
                        rhoE[(j + 2) * (nx0 + 4) + (i + 0)] =
                        rhoE[(j - 2) * (nx0 + 4) + (i + 0)];
                    }
                }
            }
            if(iteration==0 && stage == 0) {
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                executionCount = (nx1+2-(nx1+1))*(nx0+4);
                bandwidth = executionCount * 4 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
                
                result << "ABC2 Bottom" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
            }
        } // End of stage loop
        
        double sum = 0.0;
        double sum2 = 0.0;
        
        // reading dataset rho with (i,j) access
        // reading dataset p with (i,j) access
        start = std::chrono::system_clock::now();
        double glb_sum;
        double glb_sum2;
        for (int j = 0; j < nx1 + 4; j++) {
            for (int i = 0; i < nx0 + 4; i++) {
                sum += rho[j * (nx0 + 4) + i] * rho[j * (nx0 + 4) + i];
                sum2 += p[j * (nx0 + 4) + i] * p[j * (nx0 + 4) + i];
            }
        }
        MPI_Allreduce(&sum, &glb_sum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        sum = glb_sum;
        MPI_Allreduce(&sum2, &glb_sum2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        sum2 = glb_sum2;
        
        if(iteration==0) {
            end = std::chrono::system_clock::now();
            elapsed_seconds = end-start;
            executionCount = (nx1+4)*(nx0+4);
            bandwidth = executionCount * 2 * sizeof(double) * (1/elapsed_seconds.count())/1073741824;
            
            result << "Sum" << "\t" << elapsed_seconds.count() << "\t" << executionCount << " \t" << bandwidth << " \n";
        }
        std::cout << "Checksums: " << sqrt(sum) << " " << sqrt(sum2) << "\n";
        
    } // End of time loop
    
    // Record end time
    auto endMain = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = endMain - startMain;
    
    std::cout << "\nTimings are:\n";
    std::cout << "-----------------------------------------\n";
    // TODO: per-loop statistics come here
    std::cout << result.str() << "\n";
    std::cout << "Total Wall time " << diff.count() << " seconds\n";
    

    
    delete[] rho;
    delete[] rhou0;
    delete[] rhou1;
    delete[] rhoE;
    delete[] rho_old;
    delete[] rhou0_old;
    delete[] rhou1_old;
    delete[] rhoE_old;
    delete[] T;
    delete[] u0;
    delete[] u1;
    delete[] p;
    delete[] wk0;
    delete[] wk1;
    delete[] wk2;
    delete[] wk3;
    
    MPI_Finalize();
    
    return 0;
}
