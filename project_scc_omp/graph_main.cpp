#include "graph_scc_fwbw.h"
//#include "graph_scc_parallel.h"
//#include "graph_scc_serial.h"
//#include "graph_for_loop.h"
#include <omp.h>
#include <iostream>
#include <time.h>
using namespace std;
// Driver program to test above functions
int main()
{
    // create the graph given in above fugure
    //int V = 5;
    time_t start, finish;
    struct Graph* graph = init_graph_from_file();   
    cout << "Created graph"  << endl;
    // print the adjacency list representation of the above graph
    //printGraph(graph);
    //start timer for serial_coloring_scc


    time(&start);
    scc_fwbw(graph);
    //serial_color_scc(graph);
    time(&finish);
    cout<< "Serial coloring took " << difftime(finish, start) << " seconds." << endl;
    system ("pause");
    return 0;
}