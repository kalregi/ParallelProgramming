/* Copyright (c) 2012, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include "mpi_navier.h"

double simplest_checksum(int imax, int jmax, double *in) {
  double checksum = 0.0;
  int i;
  int j;
  for (j=1; j<(jmax+1); j++){
    for (i=1; i<(imax+1); i++){
      checksum+=in[(j)*(imax+2)+i]*((double)((j)*jmax)+(i));
    }
  }
  
  return checksum;
}

void print_file(int imax, int jmax, double *in) {
  char fname[10];
  sprintf(fname,"out_%d.txt",0);
  FILE *f=fopen(fname,"w");
  int i;
  int j;
  for (j=0; j<(jmax+2); j++){
    for (i=0; i<(imax+2); i++){
      fprintf(f,"%g ",in[(j)*(imax+2)+i]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
}


int main(int argc, char** argv)
{
    //Size along y
    int jmax = 4094;
    //Size along x
    int imax = 4094;
    int iter_max = 1000;
   
    mpi_setup(argc, argv, &imax, &jmax);
 
    const double pi  = 2.0 * asin(1.0);
    const double tol = 1.0e-5;
    double error     = 1.0;

    double *A;
    double *Anew;
    double *y0;

    A    = (double *)malloc((imax+4) * (jmax+4) * sizeof(double));
    Anew = (double *)malloc((imax+4) * (jmax+4) * sizeof(double));
    y0   = (double *)malloc((imax+4) * sizeof(double));

    memset(A, 0, (imax+4) * (jmax+4) * sizeof(double));
    
    // set boundary conditions
    //TODO: this loop should only done by processes who have no previous neighbour in the y direction
    if (prev_y == MPI_PROC_NULL)
    for (int i = 0; i < imax+4; i++)
      A[(0)*(imax+4)+i]   = 0.0;

    //TODO: this loop should only done by processes who have no next neighbour in the y direction
    if (next_y == MPI_PROC_NULL)
    for (int i = 0; i < imax+4; i++)
      A[(jmax+3)*(imax+4)+i] = 0.0;
    
    //TODO: this loop should only done by processes who have no previous neighbour in the x direction
    if (prev_x == MPI_PROC_NULL)
    for (int j = 0; j < jmax+4; j++)
    {
        y0[j] = sin(pi * (gbl_j_begin+j) / (jmax+3)); //TODO: within sin(), j is a global index and jmax is the global size
                                          //need to offset j by the beginning j index of this partition, and use the full jmax size
        A[(j)*(imax+4)+0] = y0[j];
    }

    //TODO: this loop should only done by processes who have no next neighbour in the x direction
    if (next_x == MPI_PROC_NULL)
    for (int j = 0; j < imax+4; j++)
    {
        y0[j] = sin(pi * (gbl_j_begin+j)/ (jmax+3)); //TODO: within sin(), j is a global index and jmax is the global size
                                          //need to offset j by the beginning j index of this partition, and use the full jmax size
        A[(j)*(imax+4)+imax+3] = y0[j]*exp(-pi);
    }
    
    //TODO: Need to set A to dirty
    set_dirty(A);

    //TODO: Only process 0 should print this, and with the full imax and jmax sizes
    printf("Jacobi relaxation Calculation: %d x %d mesh\n", imax+2, jmax+2);
    
    double t1 = omp_get_wtime();
    int iter = 0;
    
    //TODO: this loop should only done by processes who have no previous neighbour in the y direction
    if (prev_y == MPI_PROC_NULL)
    for (int i = 1; i < imax+4; i++)
       Anew[(0)*(imax+4)+i]   = 0.0;

    //TODO: this loop should only done by processes who have no next neighbour in the y direction
    if (next_y == MPI_PROC_NULL)
    for (int i = 1; i < imax+4; i++)
       Anew[(jmax+3)*(imax+4)+i] = 0.0;

    //TODO: this loop should only done by processes who have no previous neighbour in the x direction
    if (prev_x == MPI_PROC_NULL)
    for (int j = 1; j < jmax+4; j++)
        Anew[(j)*(imax+4)+0]   = y0[j];

    //TODO: this loop should only done by processes who have no next neighbour in the x direction
    if (next_x == MPI_PROC_NULL)
    for (int j = 1; j < jmax+4; j++)
        Anew[(j)*(imax+4)+jmax+3] = y0[j]*expf(-pi);
    
    //TODO: Need to set Anew to dirty
    set_dirty(Anew);

    while ( error > tol && iter < iter_max )
    {
        error = 0.0;
        //TODO: need to make sure stencil accesses to A read correct data
        exchange_halo(imax, jmax, A);

#pragma omp parallel for reduction(max:error)
        for( int j = 1; j < jmax+3; j++ )
        {
            for( int i = 1; i < imax+3; i++)
            {
                Anew[(j)*(imax+4)+i] = 0.25f * ( A[(j)*(imax+4)+i+1] + A[(j)*(imax+4)+i-1]
                                     + A[(j-1)*(imax+4)+i] + A[(j+1)*(imax+4)+i]);
                error = fmax( error, fabs(Anew[(j)*(imax+4)+i]-A[(j)*(imax+4)+i]));
            }
        }
        //TODO: need global reduction for error
        double gbl_err;
        MPI_Allreduce(&error, &gbl_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        error = gbl_err;

        //TODO: Need to set Anew to dirty
        set_dirty(Anew);
        
        //No stencil accesses to Anew, no halo exchange necessary
#pragma omp parallel for 
        for( int j = 1; j < jmax+3; j++ )
        {
            for( int i = 1; i < imax+3; i++)
            {
                A[(j)*(imax+4)+i] = Anew[(j)*(imax+4)+i];
            }
        }
        //TODO: Need to set A to dirty
        set_dirty(A);
        //TODO: Only process 0 should print this
        if(my_rank == 0 && iter % 100 == 0) printf("%5d, %0.6f\n", iter, error);
        
        iter++;
    }

    double runtime = omp_get_wtime()-t1;
    MPI_Finalize();
    printf("checksum: %d \n", simplest_checksum(imax,jmax,Anew));
    printf(" total: %f s\n", runtime);
}
