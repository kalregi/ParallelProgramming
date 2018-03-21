#include <random>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <sstream>

int main(int argc, char *argv[])
{
    int N = 4;
    
    float A[16] = {16.0,2.0,3.0,13.0,5.0,11.0,10.0,8.0,9.0,7.0,6.0,12.0,4.0,14.0,15.0,1.0};
    float B[16] = {9.0,7.0,6.0,12.0,5.0,11.0,10.0,8.0,16.0,2.0,3.0,13.0,4.0,14.0,15.0,1.0};
    float C[16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


    int rc;
    rc = MPI_Init(&argc,&argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
        
    int p; //number of tasks
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    printf ("Number of tasks = %d. My rank = %d\n", p, rank);

    
    MPI_Status status;
    float a[(N/2)*(N/2)];
    float b[(N/2)*(N/2)];
    if (rank == 0){
        //Printing A
        std::ostringstream ssA;
        for(int i = 0; i<N*N;i++){
            ssA << " " << A[i];
            
        }
        std::cout <<"A: " << ssA.str() << "\n";
        //Printing B
        std::ostringstream ssB;
        for(int i = 0; i<N*N;i++){
            ssB << " " << B[i];
            
        }
        std::cout <<"B: " << ssB.str() << "\n";
        
        for(int i = 0; i<N/2; i++){
            for(int j = 0; j<N/2; j++){
                //A[1-2,1-2]
                a[i*(N/2)+j] = A[i*N+j];
                //B[1-2,1-2]
                b[i*(N/2)+j] = B[i*N+j];
            }
        }
        
        //Counting C[1-2,1-2], part 1
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[i*N+j]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        printf ("Thread %d is sending a to thread 1\n",rank);
        //Sending a matrix to thread 1
        MPI_Send(&a,N/2*N/2,MPI_FLOAT,1,22,MPI_COMM_WORLD);
        //Receive a matrix from thread 1
        MPI_Recv(&a,N/2*N/2,MPI_FLOAT,1,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received a from thread 1\n",rank);

        printf ("Thread %d is sending b to thread 2\n",rank);
        //Sending b matrix to thread 2
        MPI_Send(&b,N/2*N/2,MPI_FLOAT,2,22,MPI_COMM_WORLD);
        //Receive b matrix from thread 2
        MPI_Recv(&b,N/2*N/2,MPI_FLOAT,2,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received b from thread 1\n",rank);
        
        //Counting C[1-2,1-2], part 2
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[i*N+j]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        //Assembling C from parts
        
        //Getting C from thread 1
        float C1[16];
        MPI_Recv(&C1,N*N,MPI_FLOAT,1,22,MPI_COMM_WORLD,&status);
        
        //Getting C from thread 2
        float C2[16];
        MPI_Recv(&C2,N*N,MPI_FLOAT,2,22,MPI_COMM_WORLD,&status);
        
        //Getting C from thread 3
        float C3[16];
        MPI_Recv(&C3,N*N,MPI_FLOAT,3,22,MPI_COMM_WORLD,&status);
        
        std::ostringstream ss;
        for(int i = 0; i<N*N;i++){
            C[i]+=C1[i]+C2[i]+C3[i];
            ss << " " << C[i];

        }
        std::cout <<"C: " << ss.str() << "\n";
    }
    if (rank == 1){
        for(int i = 0; i<N/2; i++){
            for(int j = 0; j<N/2; j++){
                //A[1-2,3-4]
                a[i*(N/2)+j] = A[i*N+j+(N/2)];
                //B[3-4,3-4]
                b[i*(N/2)+j] = B[(i+N/2)*N+j+(N/2)];
            }
        }
        //Counting C[1-2,3-4], part 1
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[i*N+j+N/2]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        printf ("Thread %d is sending a to thread 0\n",rank);
        //Sending a matrix to thread 0
        MPI_Send(&a,N/2*N/2,MPI_FLOAT,0,22,MPI_COMM_WORLD);
        //Receive a matrix from thread 0
        MPI_Recv(&a,N/2*N/2,MPI_FLOAT,0,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received a from thread 0\n",rank);

        printf ("Thread %d is sending b to thread 3\n",rank);
        //Sending b matrix to thread 3
        MPI_Send(&b,N/2*N/2,MPI_FLOAT,3,22,MPI_COMM_WORLD);
        //Receive b matrix from thread 3
        MPI_Recv(&b,N/2*N/2,MPI_FLOAT,3,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received b from thread 3\n",rank);


        //Counting C[1-2,3-4], part 2
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[i*N+j+N/2]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        //Sending result to thread 0
        MPI_Send(&C,N*N,MPI_FLOAT,0,22,MPI_COMM_WORLD);

    }
    if (rank == 2){
        for(int i = 0; i<N/2; i++){
            for(int j = 0; j<N/2; j++){
                //A[3-4,3-4]
                a[i*(N/2)+j] = A[(i+N/2)*N+j+(N/2)];
                //B[3-4,1-2]
                b[i*(N/2)+j] = B[(i+N/2)*N+j];
            }
        }
        
        //Counting C[3-4,1-2], part 1
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[(i+(N/2))*N+j]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        printf ("Thread %d is sending a to thread 3\n",rank);
        //Sending a matrix to thread 3
        MPI_Send(&a,N/2*N/2,MPI_FLOAT,3,22,MPI_COMM_WORLD);
        //Receive a matrix from thread 3
        MPI_Recv(&a,N/2*N/2,MPI_FLOAT,3,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received a from thread 3\n",rank);

        printf ("Thread %d is sending b to thread 0\n",rank);
        //Sending b matrix to thread 0
        MPI_Send(&b,N/2*N/2,MPI_FLOAT,0,22,MPI_COMM_WORLD);
        //Receive b matrix from thread 0
        MPI_Recv(&b,N/2*N/2,MPI_FLOAT,0,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received b from thread 3\n",rank);

        
        //Counting C[3-4,1-2], part 2
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[(i+(N/2))*N+j]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        //Sending result to thread 0
        MPI_Send(&C,N*N,MPI_FLOAT,0,22,MPI_COMM_WORLD);

    }
    if (rank == 3){
        for(int i = 0; i<N/2; i++){
            for(int j = 0; j<N/2; j++){
                //A[3-4,1-2]
                a[i*(N/2)+j] = A[(i+N/2)*N+j];
                //B[1-2,3-4]
                b[i*(N/2)+j] = B[i*N+j+(N/2)];
            }
        }
        //Counting C[3-4,3-4], part 1
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[(i+(N/2))*N+j+(N/2)]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        printf ("Thread %d is sending a to thread 2\n",rank);
        //Sending a matrix to thread 2
        MPI_Send(&a,N/2*N/2,MPI_FLOAT,2,22,MPI_COMM_WORLD);
        //Receive a matrix from thread 2
        MPI_Recv(&a,N/2*N/2,MPI_FLOAT,2,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received a from thread 2\n",rank);

        printf ("Thread %d is sending b to thread 1\n",rank);
        //Sending b matrix to thread 1
        MPI_Send(&b,N/2*N/2,MPI_FLOAT,1,22,MPI_COMM_WORLD);
        //Receive b matrix from thread 1
        MPI_Recv(&b,N/2*N/2,MPI_FLOAT,1,22,MPI_COMM_WORLD,&status);
        printf ("Thread %d has received b from thread 1\n",rank);

        //Counting C[3-4,3-4], part 2
        for (int i = 0; i<N/2; i++){
            for (int j = 0; j<N/2; j++){
                for (int k = 0; k<N/2; k++){
                    C[(i+(N/2))*N+j+(N/2)]+=a[i*(N/2)+k]*b[k*(N/2)+j];
                }
            }
        }
        
        //Sending result to thread 0
        MPI_Send(&C,N*N,MPI_FLOAT,0,22,MPI_COMM_WORLD);
        
    }
    MPI_Finalize();
    
    return 0;
}



