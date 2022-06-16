/*
 * ------------------------------------------------------------------------------------
 * Program for Prim's Minimum Spanning Tree (MST) algorithm. 
 * The program is for adjacency matrix representation of the graph.
 * ------------------------------------------------------------------------------------
 *  Configure:          source /opt/nfs/config/source_upcxx_2022.3.sh
 *  Compile:            UPCXX_GASNET_CONDUIT=udp upcxx -O2 main.cpp -o main 
 *  Claster setup:      /opt/nfs/config/station206_name_list.sh 1 16 > nodes
 *  Run:                upcxx-run -shared-heap 256M -n 4 $(upcxx-nodes nodes) .main input_filename [results_filename] 
 *                          -proc - number of processes
 *                          -input_filename - name of the file which contains input adjency matrix
 *                          -[results_filename] - name of the file for MST adjency matrix results
 *                      for example date: upcxx-run -shared-heap 256M -n 4 $(upcxx-nodes nodes) main ./example/data01.txt mst.txt
 * ------------------------------------------------------------------------------------
 */

#include <iostream>
#include <stdlib.h>
#include <upcxx/upcxx.hpp>
#include <limits.h>
#include <string>
#include <fstream>

#define MASTER 0    // master node number
int N;              // number of vertices


using namespace std;

int main(int argc, char *argv[])
{
    char* fileName = argv[1];
    FILE* file = fopen(fileName, "r");
    char* fileNameForResults = argv[2];
    FILE* resFile;
    
    int i, j, v, k;
    int** graph;                // adjency matrix of graph
    int rank;                   // rank of current processor
    int procNum;                // number of processors
    int* chunkSize;             // specifying the number of elements to send to each processor
    int* displs;                // specifies the displacement (starting index) needed in MPI_Scatterv
    int* chunk;                 // chunk of matrix (columns) that belongs to each processor (flatten vector of columns)
    int* mstEdges;              // minimum spanning tree by edges: mstEdges[v1] = v2
    int** MST;                  // adjency matrix of minimum spanning tree (MST) found in input graph
    double startTime, endTime;  // to measure Prim algo compute time for each process

	int* localRow = new int[2];  	// [0] - minimal weight, [1] - rank
    int* globalRow = new int[2];	// [0] - minimal weight, [1] - rank
    int minWeight, min, v1, v2;
    upcxx::global_ptr<int> globalMinWeight;

    upcxx::init();

    rank = upcxx::rank_me();
    procNum = upcxx::rank_n();

    // // master processor read adjency matrix data from file to 2D matrix (graph)
    if (rank == MASTER) {
        printf("[INFO] Master rank is %d and number of processors is %d\n", rank, procNum);    // for example date: mpiexec -f nodes -n 5 ./main data01.txt mst.txt

        // Preparing command for counting number of lines
        string grep = "grep \"\" -c ";
        string bashCommand = grep + fileName;

        // Checking if file exists
        if (file == NULL) 
            printf("Could not open this file. Choose again proper filename...\n");
        else {
            // Reading number of lines from Linux command
            FILE* res = popen(&bashCommand[0], "r");
            if (res)
                if(fscanf(res, "%d", (int*)(&N))){};

            fclose(res);

            graph = (int**)calloc(N, sizeof(int*));
            for(i=0; i < N; i++)    
                graph[i] = (int*)calloc(N, sizeof(int));


            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    if(fscanf(file, "%d", &v)){};;
                    graph[i][j] = v;
                }
            }
            fclose(file);

        } // end of else
    } // end of master node tasks

    N = upcxx::broadcast(N, MASTER).wait();

	MST = (int**)calloc(N, sizeof(int*));
	for(i=0; i < N; i++)    
		MST[i] = (int*)calloc(N, sizeof(int));

	if(rank != MASTER)
	{
		graph = (int**)calloc(N, sizeof(int*));
		for(i=0; i < N; i++)    
			graph[i] = (int*)calloc(N, sizeof(int));
	}

	for (i = 0; i < N; i++) {
		upcxx::broadcast((int **)&(graph[i][0]), N, MASTER).wait();
	}
	

	chunkSize = (int*)calloc(procNum, sizeof(int));
	displs = (int*)calloc(procNum, sizeof(int));

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
	
	chunk = (int*)calloc(chunkSize[rank]*N, sizeof(int));     // each processor needs its own chunk of data

	int t = 0;
	for(i = displs[rank]; i < displs[rank] + chunkSize[rank]; i++){
		for(j = 0; j < N; j++) {
			chunk[t++] = graph[i][j];
		}
	}

	for(i = 0; i < chunkSize[rank]*N; i++) {
        cout << "chunk[i=" << i << "]=" << chunk[i] << endl;;
    }
	cout << endl;


	mstEdges = (int*)calloc(N, sizeof(int));   // max size is number of vertices
	
	mstEdges[0] = 0;
	for (i = 1; i < N; ++i)
	{
		mstEdges[i] = -1;
	}

	globalMinWeight = upcxx::new_<int>(0);
	v1 = v2 = 0;

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

		localRow[0] = minWeight;
		localRow[1] = rank;
            
        // each process have to send its min weight row to others processes
		upcxx::reduce_all(localRow, globalRow, 2, upcxx::op_fast_min).wait();
        v1 = v1 + displs[rank];
        v2 = v2;

        // broadcasts informations from the rank process to all other processes of the communicator
        v1 = upcxx::broadcast(v1, globalRow[1]).wait();
		v2 = upcxx::broadcast(v2, globalRow[1]).wait();
		
        mstEdges[v2] = v1;
        MST[v1][v2] = globalRow[0];
        MST[v2][v1] = globalRow[0];

		upcxx::future<int> fut = upcxx::rget(globalMinWeight);
		int weight = fut.wait() + globalRow[0];
		rput(weight, globalMinWeight).wait();
            
    } // end of loop over vertices

    upcxx::barrier();

    if (rank == MASTER) 
    {
        cout << endl << "MST edges with weights:" << endl;
        cout << "------------------------" << endl;

        for (i = 1; i < N; i++) 
            cout << " " << i << " --- " << mstEdges[i] << "  \t" << MST[i][mstEdges[i]] << endl;

		upcxx::future<int> future = upcxx::rget(globalMinWeight);
		int resultMinWeight = future.wait();

        cout << endl << endl << "MST weight: " << resultMinWeight << endl;
        cout << "--------------------------------------------------------------------------------" << endl;

        if (fileNameForResults != NULL) {
            cout << "[INFO] Saving results to file..." << endl;
            resFile=fopen(fileNameForResults,"w+");
            if (resFile==NULL)
            {
                cout << "Couldn't save results.. " << fileNameForResults << " file could not be opened" << endl;
            } else {
                for (i = 0; i < N; i++)
                {   
                    for (j = 0; j < N; j++)
                        fprintf(resFile,"%d ", MST[i][j]);
                    fprintf(resFile,"\n");
                }
                fclose(resFile);
                cout << "[INFO] Adjency matrix of found MST saved to file: " << fileNameForResults << "." << endl;
            }
        } else {
            cout << "[INFO] If you want to save MST adjency matrix to file add filename to args." << endl;
        }
    }

    // // dealocating memory
    for(i = 0; i < N; ++i) {
        free(graph[i]);
        free(MST[i]);
    }
    free(graph);
    free(MST);
    free(chunkSize);
    free(displs);
    free(chunk);
    free(mstEdges);

    if (rank == MASTER) {  
        cout << "[DEBUG] END OF PROGRAM" << endl;
    }
    
	upcxx::delete_(globalMinWeight);
    upcxx::finalize();
    
    return 0;
}