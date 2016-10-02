//
// Copyright (c) 2004-2006 The Trustees of Indiana University and Indiana
//                         University Research and Technology
//                         Corporation.  All rights reserved.
// Copyright (c) 2006      Cisco Systems, Inc.  All rights reserved.
//
// Sample MPI "hello world" application in C++
//
// NOTE: The MPI C++ bindings were deprecated in MPI-2.2 and removed
// from the standard in MPI-3.  Open MPI still provides C++ MPI
// bindings, but they are no longer built by default (and may be
// removed in a future version of Open MPI).  You must
// --enable-mpi-cxx when configuring Open MPI to enable the MPI C++
// bindings.
//
#include <cstdlib>
#include <mpi.h>
#include <iomanip>
#include <iostream>
#define proc_A 0
#define proc_B 1
#define ping 101
#define pong 101
int buffer [101];
void processor_A (void);
void processor_B (void);
using namespace std;
int main(int argc, char **argv)
{
    int rank, size, len;
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == proc_A) processor_A();
    else if (rank == proc_B) processor_B();
    MPI_Finalize();
    return 0;
}
void processor_A ( void){
MPI_Status status;
double start, finish, time;
extern int buffer [101];
int i, length;
  for (length =1; length<= 101; length +=10){
    start = MPI_Wtime();
    for (i = 1; i<=10; i++){
      MPI_Send(buffer, length, MPI_INT, proc_B, ping, MPI_COMM_WORLD);
      MPI_Recv(buffer, length, MPI_INT, proc_B, pong, MPI_COMM_WORLD, &status);
    }
    finish = MPI_Wtime();
    time = finish - start;

    cout << "Elapsed time for ping pong with length " << length << " was " << time << " seconds." << endl;

  }  
}
void processor_B ( void){
MPI_Status status;
double start, finish, time;
extern int buffer [101];
int i, length;
  for (length =1; length<= 101; length +=10){
    for (i = 1; i<=10; i++){
      MPI_Recv(buffer, length, MPI_INT, proc_A, ping, MPI_COMM_WORLD, &status);
      MPI_Send(buffer, length, MPI_INT, proc_A, pong, MPI_COMM_WORLD);
    }
  }  
}
