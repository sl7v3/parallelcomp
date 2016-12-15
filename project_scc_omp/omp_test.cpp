#include <omp.h>
#include <iostream>
#include <time.h>
using namespace std;
// Driver program to test above functions
int main()
{
    cout << "Maximum number of threads available in OpenMP is " << omp_get_max_threads() << endl;
    int count = 20;
    int i = 0;
    for (i; i < count; ++i)
    {
      if ( i == 9){
        break;
      }
      /* code */
    }
    cout << "I is " << i << endl;
    // print the adjacency list representation of the above graph
    //printGraph(graph);
    //start timer for serial_coloring_scc
    system ("pause");
    return 0;
}