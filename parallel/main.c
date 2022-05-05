/*
 * ------------------------------------------------------------------------------------
 * Program for Prim's Minimum Spanning Tree (MST) algorithm. 
 * The program is for adjacency matrix representation of the graph.
 * ------------------------------------------------------------------------------------

 *  Configure:          source /opt/nfs/config/source_mpich32.sh
 *  Compile:            mpicc main.c -o main 
 *  Claster setup:      /opt/nfs/config/station206_name_list.sh 1 16 > nodes
 *  Run:                mpiexec -f nodes -n 16 ./main input_filename [results_filename]
 *                          -proc - number of processes
 *                          -input_filename - name of the file which contains input adjency matrix
 *                          -[results_filename] - name of the file for MST adjency matrix results
 *                      for example date: mpiexec -f nodes -n 5 ./main data01.txt mst.txt
 * ------------------------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "mpi.h"

#define MASTER 0    // master node number
int N;              // number of vertices



int main(int argc, char *argv[])
{
    char* fileName = argv[1];
    FILE* file = fopen(fileName, "r");
    char* fileNameForResults = argv[2];
    FILE* resFile;
    
    int i, j, v, k;
    int** graph;                // adjency matrix of graph
    int* flatten_graph;         // graph converted to vector
    int rank;                   // rank of current processor
    int procNum;                // number of processors
    int* chunkSize;             // specifying the number of elements to send to each processor
    int* displs;                // specifies the displacement (starting index) needed in MPI_Scatterv
    int* chunk;                 // chunk of matrix (columns) that belongs to each processor (flatten vector of columns)
    int* mstEdges;              // minimum spanning tree by edges: mstEdges[v1] = v2
    int** MST;                  // adjency matrix of minimum spanning tree (MST) found in input graph
    double startTime, endTime;  // to measure Prim algo compute time for each process

    // struct contains tuple of minimal weight and rank for row
    struct { 
        int minWeight; 
        int rank; 
    } localRow, globalRow;
                  
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
                printf("[INFO] Proc number > graph size.. reducing... will be computed only by %d procceses.\n", N);
                printf("--------------------------------------------------------------------------------\n\n");
            }
        } // end of else
    } // end of master node tasks

    MPI_Bcast(&N, 1, MPI_INT, MASTER, world);
    MPI_Bcast(&procNum, 1, MPI_INT, MASTER, world);


    if (rank < procNum) {
        MST = calloc(N, sizeof(int*));
        for(i=0; i < N; i++)    
            MST[i] = calloc(N, sizeof(int));
    
        if(rank != MASTER)
        {
            graph = calloc(N, sizeof(int*));
            for(i=0; i < N; i++)    
                graph[i] = calloc(N, sizeof(int));
        }

        for (i = 0; i < N; i++) 
            MPI_Bcast((int **)&(graph[i][0]), N, MPI_INT, MASTER, world);
        

        
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


        flatten_graph = calloc(N*N, sizeof(int));
        
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                flatten_graph[N*i+j] = graph[i][j];
            }
        }

        // here the chunk each processor needs will be scatter to it
        MPI_Scatterv(&(flatten_graph[0]), chunkSize, displs, MPI_INT, chunk, chunkSize[rank], MPI_INT, MASTER, world);


        for (i = 0; i < procNum; ++i){
            chunkSize[i] /= N;
            displs[i] /= N;
        }

        mstEdges = calloc(N, sizeof(int));   // max size is number of vertices

        
        mstEdges[0] = 0;
        for (i = 1; i < N; ++i)
        {
            mstEdges[i] = -1;
        }

        globalMinWeight = 0;
        v1 = v2 = 0;

        startTime = MPI_Wtime();    // start measuring compute time

        // iterates over vertices
        for (k = 0; k < N-1; ++k) 
        {
            minWeight = INT_MAX;
            
            // iterates over vertices in chunk
            for (i = 0; i < chunkSize[rank]; ++i)
            {
                // check if beginning of edge has not been visited
                if (mstEdges[i + displs[rank]] != -1) 
                {
                    // iterates over vertices (end of our edge which we are looking for by minimal weight)
                    for (j = 0; j < N; ++j)
                    {
                        // check if end of edge has not been visited
                        if (mstEdges[j] == -1) 
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
            
            // each process have to send its min weight row to others processes
            MPI_Allreduce(&localRow, &globalRow, 1, MPI_2INT, MPI_MINLOC, world); 
            v1 = v1 + displs[rank];
            v2 = v2;

            // broadcasts informations from the rank process to all other processes of the communicator
            MPI_Bcast(&v1, 1, MPI_INT, globalRow.rank, world);
            MPI_Bcast(&v2, 1, MPI_INT, globalRow.rank, world);
            
            mstEdges[v2] = v1;
            MST[v1][v2] = globalRow.minWeight;
            MST[v2][v1] = globalRow.minWeight;

            globalMinWeight += globalRow.minWeight;
            
        } // end of loop over vertices

        endTime = MPI_Wtime();
        if (rank < N)
            printf("Compute time from %d process: %.2f ms\n", rank, (endTime-startTime)*1000);

        MPI_Barrier(world);

        if (rank == MASTER) {
            printf("\nMST edges with weights:\n");
            printf("------------------------\n");

            for (i = 1; i < N; i++) 
                printf(" %d --- %d  \t%d\n", i, mstEdges[i], MST[i][mstEdges[i]]);

            printf("\n\nMST weight: %d\n", globalMinWeight);
            printf("--------------------------------------------------------------------------------\n");

            if (fileNameForResults != NULL) {
                printf("[INFO] Saving results to file...\n");
                resFile=fopen(fileNameForResults,"w+");
                if (resFile==NULL)
                {
                    printf("Couldn't save results..%s file could not be opened", fileNameForResults);
                } else {
                    for (i = 0; i < N; i++)
                    {   
                        for (j = 0; j < N; j++)
                            fprintf(resFile,"%d ", MST[i][j]);
                        fprintf(resFile,"\n");
                    }
                    fclose(resFile);
                    printf("[INFO] Adjency matrix of found MST saved to file: %s.\n", fileNameForResults);
                }
            } else {
                printf("[INFO] If you want to save MST adjency matrix to file add filename to args.\n");
            }
        }
  
    }   

    // dealocating memory
    for(i = 0; i < N; ++i) {
        free(graph[i]);
        free(MST[i]);
    }
    free(graph);
    free(MST);
    free(flatten_graph);
    free(chunkSize);
    free(displs);
    free(chunk);
    free(mstEdges);

    MPI_Finalize();

    return 0;
}