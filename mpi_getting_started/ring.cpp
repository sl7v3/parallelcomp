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
#define ping 101
#define pong 101
int buffer [101];
void processor_A (void);
void processor_B (void);
void processor_C (void);
using namespace std;
int main(int argc, char **argv)
{
    int rank, size, len;
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int proc_A = 0;
    int proc_C = size-1;
    if (rank == proc_A) processor_A();
    else if (rank == proc_C) processor_C();
    else processor_B();
    MPI_Finalize();
    return 0;
}
void processor_A ( void){
MPI_Status status;
//double start, finish, time;
int rank, size, sum = 0;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
extern int buffer [101];
buffer[0] = rank;
int i, length = 1;
  for (int sends =0; sends< size; sends++){
    MPI_Send(buffer, length, MPI_INT, rank+1, ping, MPI_COMM_WORLD);
    MPI_Recv(buffer, length, MPI_INT, size-1, pong, MPI_COMM_WORLD, &status);
    sum+= buffer[0];
  }
  cout << "Sum of all ranks is received by processor at rank " <<rank<< " is "  << sum << endl;
}
void processor_B ( void){
MPI_Status status;
//double start, finish, time;
int rank, size, sum = 0;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
extern int buffer [101];
buffer[0] = rank;
int i, length = 1;
  for (int sends =0; sends< size; sends++){
    MPI_Send(buffer, length, MPI_INT, rank+1, ping, MPI_COMM_WORLD);
    MPI_Recv(buffer, length, MPI_INT, rank-1, pong, MPI_COMM_WORLD, &status);
    sum+= buffer[0];
  }
  cout << "Sum of all ranks is received by processor at rank " <<rank<< " is "  << sum << endl;
}

void processor_C ( void){
MPI_Status status;
//double start, finish, time;
int rank, size, sum = 0;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
extern int buffer [101];
buffer[0] = rank;
int i, length = 1;
  for (int sends =0; sends< size; sends++){
    MPI_Send(buffer, length, MPI_INT, 0, ping, MPI_COMM_WORLD);
    MPI_Recv(buffer, length, MPI_INT, rank-1, pong, MPI_COMM_WORLD, &status);
    sum+= buffer[0];
  }
  cout << "Sum of all ranks is received by processor at rank " <<rank<< " is "  << sum << endl;
}
