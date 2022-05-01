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
#include <stdbool.h>
#include <string.h>
#include <mpi.h>
#define MASTER 0        // Master node number
 
int N;

// Find the vertex with minimum key value, 
// from the set of vertices not yet included in MST
int minKey(int* key, bool* mstSet)
{
    // Initialize min value
    int min = INT_MAX, min_index;
    for (int v = 0; v < N; v++) {
        if (mstSet[v] == false && key[v] < min)
            min = key[v], min_index = v;
    }

    return min_index;
}
 

// Print constructed MST stored in parent[]
int printMST(int* parent, int** graph)
{
    printf("Edge \tWeight\n");
    for (int i = 1; i < N; i++)
        printf("%d - %d \t%d \n", i+1, parent[i]+1, graph[i][parent[i]]);
}
 

// Construct and print MST for a graph represented 
// using adjacency matrix representation
void primMST(int** graph)
{
    int* parent = calloc(N, sizeof(int));        // Array to store constructed MST
    int* key = calloc(N-1, sizeof(int));         // Key values used to pick minimum weight edge in cut
    bool* mstSet = calloc(N, sizeof(bool));      // To represent set of vertices included in MST
 
    // Initialize all keys as INFINITE
    for (int i = 0; i < N; i++)
        key[i] = INT_MAX, mstSet[i] = false;
 
    // Include first 1st vertex in MST.
    // Make key 0 so that this vertex is picked as first vertex.
    key[0] = 0;
    parent[0] = -1; // First node is always root of MST
 

    for (int count = 0; count < N - 1; count++) {

        // Pick the minimum key vertex from the set of vertices not yet included in MST
        int u = minKey(key, mstSet);
 
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
 
        // Update key value and parent index of the adjacent vertices of the picked vertex.
        // Consider only those vertices which are not yet included in MST
        for (int v = 0; v < N; v++)
 
            // graph[u][v] is non zero only for adjacent vertices of m
            // mstSet[v] is false for vertices not yet included in MST
            // Update the key only if graph[u][v] is smaller than key[v]
            if (graph[u][v] && mstSet[v] == false && graph[u][v] < key[v])
                parent[v] = u, key[v] = graph[u][v];
    }

    printMST(parent, graph);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
 
int main(int argc, char *argv[])
{
    FILE* file = fopen(argv[1], "r");
    char* fileName = argv[1];
    int** graph;        // adjency matrix of graph
    int i, j, v;

    int rank;           // rank of current processor
    int procNum;        // number of processors
    int* chunkSizes;    // specifying the number of elements to send to each processor 
    int* displ;         // specifies the displacement (starting index) needed in MPI_Scatterv
    int* chunk;         // chunk of matrix (columns) that belongs to each processor

    MPI_Comm world world = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &procNum);
    
    chunkSizes = (int*)malloc(sizeof(int) * procNum);
    displs = (int*)malloc(sizeof(int) * procNum);

    chunkSize = N / procNum;                // chunk size for each processor
    int remains = procNum - (N % procNum);  // if number of vertices is not a multiply of number of processors
    chunkSizes[0] = chunkSize;
    displs[0] = 0;

    for (i = 1; i < procNum; i++) {
        chunkSizes[i] = chunkSize;
        if (i < remains) 
            ++chunkSizes[i];
        displs[i] = displs[i-1] + chunkSizes[i-1];
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
        for(i=0; i< N; i++)    
            graph[i] = calloc(N, sizeof(int));


        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                fscanf(file, "%d", &v);
                graph[i][j] = v;
            }
        }
        fclose(file);
    } // end of master node tasks


    // constructed MST and print results
    primMST(graph);
    }

    // dealocating memory
    for(i = 0; i < N; ++i)
        free(graph[i]);
    free(graph);
    free(chunkSizes);
    free(displ);

    return 0;
}