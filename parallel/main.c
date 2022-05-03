/*
 * ------------------------------------------------------------------------------------
 * Program for Prim's Minimum Spanning Tree (MST) algorithm. 
 * The program is for adjacency matrix representation of the graph.
 * ------------------------------------------------------------------------------------
 * Compile: gcc main.c -o main 
 * Run:     mpicc -n proc ./main filename
 *          where: 
 *                  -proc - number of processes
 *                  -filename - name of the file which contains input adjency matrix
 *          for example date: ./main data01.txt
 * ------------------------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <mpi.h>

#define MASTER 0            // Master node number
int N;                      // number of vertices



int main(int argc, char *argv[])
{
    char* fileName = argv[1];
    FILE* file = fopen(fileName, "r");
    int** graph;            // adjency matrix of graph
    int i, j, v;

    int rank;               // rank of current processor
    int procNum;            // number of processors
    int* chunkSize;         // specifying the number of elements to send to each processor     // chunkSize
                            // liczba kolumn
    int* displ;             // specifies the displacement (starting index) needed in MPI_Scatterv
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

    int minWeight, v1, v2;    // for fiding localRow.minWeight
    


    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &procNum);
    
    chunkSize = (int*)malloc(sizeof(int) * procNum);
    displs = (int*)malloc(sizeof(int) * procNum);

    int remains = procNum - (N % procNum);  // if number of vertices is not a multiply of number of processors
    chunkSize[0] = N / procNum;
    displs[0] = 0;

    for (i = 1; i < procNum; i++) {
        chunkSize[i] = N / procNum;
        if (i < remains) 
            ++chunkSize[i];
        displs[i] = displs[i-1] + chunkSize[i-1];
    }


    // master processor read adjency matrix data from file to 2D matrix (graph)
    if (rank == MASTER) {
        printf("rank is %d and number of processors is %d\n", rank, procNum);

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

        free(bashCommand);

        graph = calloc(N, sizeof(int *));
        for(i=0; i < N; i++)    
            graph[i] = calloc(N, sizeof(int));


        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                fscanf(file, "%d", &v);
                graph[i][j] = v;
            }
        }
        fclose(file);
    } // end of master node tasks


    chunk = (int*)malloc(chunkSize[rank]*N*sizeof(int)); // each processor needs its own chunk of data

    MPI_Datatype vectorType; 
    MPI_Type_contiguous(N, MPI_INT, &vectorType);
    MPI_Type_commit(&vectorType);

    // here the chunk each processor needs will be scatter to it
    MPI_Scatterv(matrix, chunkSize, displs, vectorType, chunk, chunkSize[rank], vectorType, 0, world);
 

    MST = (int*)malloc(N*sizeof(int)); // max size is number of vertices

    MST[0] = 0;
    minWeight = 0;
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
        
        // each process have to send its min weight row to others processes
        MPI_Allreduce(&localRow, &globalRow, 1, MPI_2INT, MPI_MINLOC, world); 
        edge.v1 = v1 + displs[rank];
        edge.v2 = v2;

        // broadcasts informations from the rank process to all other processes of the communicator
        MPI_Bcast(&edge, 1, MPI_2INT, globalRow.rank, world);

        MST[edge.v2] = edge.v1;
        minWeight += globalRow.minWeight;
    } // end of loop over vertices


    if (rank == MASTER) {
        printf("MST weight: %d\n", minWeight);
        // dodać czas działania
        // zapisanie macierzy sąsiedztwa MST do pliku
        // printowanie krawędzi ?
    }


    // dealocating memory
    for(i = 0; i < N; ++i)
        free(graph[i]);
    free(graph);
    free(chunkSize);
    free(displ);
    free(chunk);
    free(MST);

    MPI_Finalize();

    return 0;
}