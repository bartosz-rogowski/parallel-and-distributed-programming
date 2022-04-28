/*
 * ------------------------------------------------------------------------------------
 * Program for Prim's Minimum Spanning Tree (MST) algorithm. 
 * The program is for adjacency matrix representation of the graph.
 * ------------------------------------------------------------------------------------
 * Compile: gcc main.c -o main 
 * Run:     ./main filename N
 *          where: -filename - name of the file which contains input adjency matrix
 *                 -N        - size of adjency matrix (2D NxN) 
 *          for example date: ./main data01.txt 5
 * ------------------------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
 
int N;


// Find the vertex with minimum key value, 
// from the set of vertices not yet included in MST
int minKey(int key[], bool mstSet[])
{
    // Initialize min value
    int min = INT_MAX, min_index;
 
    for (int v = 0; v < N; v++)
        if (mstSet[v] == false && key[v] < min)
            min = key[v], min_index = v;
 
    return min_index;
}
 

// Print constructed MST stored in parent[]
int printMST(int parent[], int graph[N][N])
{
    printf("Edge \tWeight\n");
    for (int i = 1; i < N; i++)
        printf("%d - %d \t%d \n", i+1, parent[i]+1, graph[i][parent[i]]);
}
 

// Construct and print MST for a graph represented 
// using adjacency matrix representation
void primMST(int graph[N][N])
{
    int parent[N];  // Array to store constructed MST
    int key[N];     // Key values used to pick minimum weight edge in cut
    bool mstSet[N]; // To represent set of vertices included in MST
 
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

 
int main(int argc, char *argv[])
{
  char* fn = argv[1];
  N = atoi(argv[2]);
  FILE* file = fopen(fn, "r");
  int graph[N][N];
  int MST[N];
  int i, j, v;


  if (file == NULL)
    printf("Could not open this file. Choose again proper filename..\n");
  else {
    
    // read adjency matrix data from file to 2D matrix (graph)
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
          fscanf(file, "%d", &v);
          graph[i][j] = v;
      }
    }

    // print input read data
    printf("Input data (adjency matrix): \n");
    printf("------------\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%d ", graph[i][j]);
        }
        printf("\n");
    }
    printf("------------\n\n");

    // constructed MST and print results
    primMST(graph);
  }

  fclose(file);
 
  return 0;
}