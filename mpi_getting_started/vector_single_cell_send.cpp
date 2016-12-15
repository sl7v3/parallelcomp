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
#define NBLOCKS 2
#define ROWS 10
#define COLS 15
int my_ary [ROWS][COLS];
void print_ary(int ary[][COLS]);
void init_ary(int ary[][COLS]);
void processor_A (int ary[][COLS], int rank);
void processor_B (int ary[][COLS], int rank);
using namespace std;
int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  init_ary(my_ary);
  int rank, size, len;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int proc_A = 0;
  if (rank == proc_A) {
    cout << my_ary[5][5] << " is initial value for my_ary[5][5]" << endl;
    cout << my_ary[9][5] << " is initial value for my_ary[9][5]" << endl;
    print_ary(my_ary);
    processor_A(my_ary, rank);
  }
  else processor_B(my_ary, rank);
  MPI_Finalize();
  return 0;
}

void print_ary(int ary[][COLS]){

 for (int i = 0;i<ROWS; i++){
    for (int j = 0; j< COLS; j++){
      cout << ary[i][j] << " ";
    }
    cout << endl;
  }      
}

void init_ary(int ary [][COLS]){
 for (int i  = 0;i<ROWS; i++){
    for (int j = 0; j<COLS; j++){
      ary[i][j] = i*j;
    }
  }      
}
void processor_A (int ary[][COLS], int rank){
MPI_Status status;
int length =1;
MPI_Datatype newType;
MPI_Type_vector(ROWS, 1, COLS, MPI_INT, &newType);MPI_Type_commit (&newType);
cout << my_ary[5][5]<< " is the initial value for my_ary[5][5] from " << rank <<"th processor."<<endl;
MPI_Send(&(my_ary[5][5]), length, newType, rank+1, ping, MPI_COMM_WORLD);
MPI_Recv(&(my_ary[5][5]), length, newType, rank+1, pong, MPI_COMM_WORLD, &status);
cout << my_ary[5][5]<< " is the final value for my_ary[5][5] from " << rank <<"th processor."<<endl;
}
void processor_B (int ary[][COLS], int rank){
MPI_Datatype newType;
MPI_Type_vector(ROWS, 1, COLS, MPI_INT, &newType);MPI_Type_commit (&newType);
MPI_Status status;
int length =1;
cout << my_ary[9][5]<< " is the initial value for my_ary[9][5] from " << rank <<"th processor."<<endl;
MPI_Send(&(my_ary[9][5]), length, newType, rank-1, ping, MPI_COMM_WORLD);
MPI_Recv(&(my_ary[9][5]), length, newType, rank-1, pong, MPI_COMM_WORLD, &status);
cout << my_ary[9][5]<< " is the final value for my_ary[9][5] from " << rank <<"th processor."<<endl;
}
