/*
 * ------------------------------------------------------------------------------------
 * Program for Prim's Minimum Spanning Tree (MST) algorithm. 
 * The program is for adjacency matrix representation of the graph.
 * ------------------------------------------------------------------------------------

 *  Configure:          source /opt/nfs/config/source_mpich32.sh
 *  Compile:            mpicc main.c -o main 
 *  Claster setup:      /opt/nfs/config/station206_name_list.sh 1 16 > nodes
 *  Run:                mpiexec -f nodes -n 16 ./main filename
 *                          -proc - number of processes
 *                          -filename - name of the file which contains input adjency matrix
 *                      for example date: mpiexec -f nodes -n 16 ./main data01.txt
 * ------------------------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "mpi.h"

#define MASTER 0            // Master node number
int N;                      // number of vertices



int main(int argc, char *argv[])
{
    char* fileName = argv[1];
    FILE* file = fopen(fileName, "r");
    int** graph;            // adjency matrix of graph
    int i, j, v, k;

    int rank;               // rank of current processor
    int procNum;            // number of processors
    int* chunkSize;         // specifying the number of elements to send to each processor
                            // liczba kolumn
    int* displs;            // specifies the displacement (starting index) needed in MPI_Scatterv
    int* chunk;             // chunk of matrix (columns) that belongs to each processor
                            // wektor elementów należących do danego procesu (rozplaszczona macierz kolumn które ma przetworzyc dany proces)
    int* MST;               // minimum spanning tree with selected edges: MST[v1] = v2

    // struct contains tuple of minimal weight and rank for row
    struct { 
        int minWeight; 
        int rank; 
    } localRow, globalRow;

    // struct contains tuple of vertices
    struct { 
        int v1; 
        int v2; 
    } edge;                   

    int globalMinWeight, minWeight, min, v1, v2; 


    MPI_Init(&argc, &argv);
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &procNum);


    // master processor read adjency matrix data from file to 2D matrix (graph)
    if (rank == MASTER) {
        printf("[INFO] Master rank is %d and number of processors is %d\n", rank, procNum);

        // Preparing command for counting number of lines
        char* grep = "grep \"\" -c ";
        char* bashCommand = malloc(strlen(fileName) + strlen(grep));
        strcpy(bashCommand, grep);
        strcat(bashCommand, fileName);

        // Checking if file exists
        if (file == NULL)
            printf("Could not open this file. Choose again proper filename..\n");
        else {
            // Reading number of lines from Linux command
            FILE* res = popen(bashCommand, "r");
            if (res)
                fscanf(res, "%d", (int*)(&N));

            fclose(res);
            free(bashCommand);

            graph = calloc(N, sizeof(int*));
            for(i=0; i < N; i++)    
                graph[i] = calloc(N, sizeof(int));


            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    fscanf(file, "%d", &v);
                    graph[i][j] = v;
                }
            }
            fclose(file);

            if (procNum > N) {
                procNum = N;
                printf("[INFO] Proc number > graph size.. reducing to %d\n", procNum);
            }

        } // end of else
    } // end of master node tasks

    // printf("rank = %d\t procNum = %d \t N = %d\n", rank, procNum, N);

    // if(rank >= procNum - N) {
    //     printf("rank = %d\t procNum = %d\n", rank, procNum);
    //     // MPI_Finalize();
    // }

    MPI_Bcast(&N, 1, MPI_INT, MASTER, world);
    MPI_Bcast(&procNum, 1, MPI_INT, MASTER, world);
    
    if(rank != MASTER)
    {
        graph = calloc(N, sizeof(int*));
        for(i=0; i < N; i++)    
            graph[i] = calloc(N, sizeof(int));
    }

    for (i = 0; i < N; i++) 
        MPI_Bcast((int **)&(graph[i][0]), N, MPI_INT, MASTER, MPI_COMM_WORLD);
    

    
    chunkSize = calloc(procNum, sizeof(int));
    displs = calloc(procNum, sizeof(int));

    int remains = (N % procNum);  // if number of vertices is not a multiply of number of processors
    displs[0] = 0;
    
    for (i = 0; i < procNum; i++) {
        chunkSize[i] = N / procNum;
        if (i < remains) 
            ++chunkSize[i];
    }


    for (i = 1; i < procNum; i++) {
        displs[i] = displs[i-1] + chunkSize[i-1];
    }
    
    chunk = calloc(chunkSize[rank]*N, sizeof(int));     // each processor needs its own chunk of data
    

    for (i = 0; i < procNum; ++i) {
        chunkSize[i] *= N;
        displs[i] *= N;
    }


    int* flatten_graph = calloc(N*N, sizeof(int));
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            flatten_graph[N*i+j] = graph[i][j];
        }
    }

    // here the chunk each processor needs will be scatter to it
    MPI_Scatterv(&(flatten_graph[0]), chunkSize, displs, MPI_INT, chunk, chunkSize[rank], MPI_INT, MASTER, world);

 
    // if (rank == 4) {
    //     printf("rank %d \t chunk = [", rank);
    //     for (i = 0; i < chunkSize[rank]; ++i)
    //         printf("%d, ", chunk[i]);
    //     printf("]\n");
    //     printf("-------------------\n");
    // }

    for (i = 0; i < procNum; ++i)
        chunkSize[i] /= N;

    MST = calloc(N, sizeof(int));   // max size is number of vertices

    MST[0] = 0;
    for (i = 1; i < N; ++i)
    {
        MST[i] = -1;
    }

    globalMinWeight = 0;
    v1 = v2 = 0;

    
    // iterates over vertices
    for (k = 0; k < N-1; ++k) 
    {
        minWeight = INT_MAX;
        
        // iterates over vertices in chunk
        for (i = 0; i < chunkSize[rank]; ++i)
        {
            // check if beginning of edge has not been visited
            if (MST[i + displs[rank]] != -1) 
            {
                // iterates over vertices (end of our edge which we are looking for by minimal weight)
                for (j = 0; j < N; ++j)
                {
                    // check if end of edge has not been visited
                    if (MST[j] == -1) 
                    {
                        // if the chunk[N*i+j] is less than minWeight value
                        if (chunk[N*i+j] < minWeight && chunk[N*i+j] != 0)
                        {
                            minWeight = chunk[N*i+j];   // current minimal weight
                            v1 = i;                     // beginning of current edge with minimal weight
                            v2 = j;                     // end of current edge with minimal weight
                        }
                    }
                } // end of loop over j
            }
        } // end of loop over vertices in chunk

        localRow.minWeight = minWeight;
        localRow.rank = rank;

        // printf("[DEBUG] rank=%d\t minWeight=%d\t v1=%d\t v2=%d\n", 
            // rank, minWeight, v1, v2);
        
        // each process have to send its min weight row to others processes
        MPI_Allreduce(&localRow, &globalRow, 1, MPI_2INT, MPI_MINLOC, world); 
        edge.v1 = v1 + displs[rank];
        edge.v2 = v2;

        // broadcasts informations from the rank process to all other processes of the communicator
        MPI_Bcast(&edge, 1, MPI_2INT, globalRow.rank, world);

        MST[edge.v2] = edge.v1;
        globalMinWeight += globalRow.minWeight;
    } // end of loop over vertices

    MPI_Barrier(world);

    if (rank == MASTER) {
        printf("---------------------------------\n");;
        printf("MST weight: %d\n", globalMinWeight);
        // dodać czas działania
        // zapisanie macierzy sąsiedztwa MST do pliku
        // printowanie krawędzi ?
    }


    // dealocating memory
    for(i = 0; i < N; ++i)
        free(graph[i]);
    free(graph);
    free(flatten_graph);
    free(chunkSize);
    free(displs);
    free(chunk);
    free(MST);

    MPI_Finalize();

    return 0;
}