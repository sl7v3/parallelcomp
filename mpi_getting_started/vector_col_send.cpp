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
#define ROWS 10
#define COLS 18
void print_ary(int ary[][COLS]);
void init_ary(int (&ary)[ROWS][COLS], int rank);
void processor_A (int rank);
void processor_B (int rank);
using namespace std;
int main(int argc, char **argv)
{
  int rank, size, len;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int proc_A = 0;
  if (rank == proc_A) {
    processor_A(rank);
  }
  else {
    processor_B(rank);
  }
  MPI_Finalize();
  return 0;
}

void print_ary(int ary[ROWS][COLS]){

 for (int i = 0;i<ROWS; i++){
    for (int j = 0; j< COLS; j++){
      cout << ary[i][j] << " ";
    }
    cout << endl;
  }      
}

void init_ary(int (&ary)[ROWS][COLS], int rank){
 for (int i  = 0;i<ROWS; i++){
    for (int j = 0; j<COLS; j++){
      ary[i][j] = i*j*(rank+1);
    }
  }      
}
void processor_A ( int rank){
int my_ary [ROWS][COLS];
init_ary(my_ary, rank);
cout<< "Initial array for " << rank << "th processor is below "<<endl;
print_ary(my_ary);
MPI_Status status;
int ierr, length =1;
MPI_Datatype newType;
MPI_Type_vector(ROWS+2, length, COLS, MPI_INT, &newType);MPI_Type_commit (&newType);
MPI_Send(&(my_ary[0][0]), ROWS, newType, rank+1, ping, MPI_COMM_WORLD);
MPI_Recv(&(my_ary[0][COLS-1]), ROWS, newType, rank+1, pong, MPI_COMM_WORLD, &status);
ierr = MPI_Type_free(&newType);
cout<< "Final array for " << rank << "th processor is below "<<endl;
print_ary(my_ary);
}
void processor_B (int rank){
int ierr, length =1, my_ary [ROWS][COLS];
MPI_Status status;
init_ary(my_ary, rank);
cout<< "Initial array for " << rank << "th processor is below "<<endl;
print_ary(my_ary);
MPI_Datatype newType;
/*Tried different combinations of block count for MPI_Type_vector for column sends  before settling on ROWS+2. It works for all columns. Just trying ROWS shifts values in to the right of the last column while ROWS+1 works in all cases except the last column where it//creates a segmentation fault.*/
MPI_Type_vector(ROWS+2, length, COLS, MPI_INT, &newType);MPI_Type_commit (&newType);
MPI_Send(&(my_ary[0][2]), ROWS, newType, rank-1, ping, MPI_COMM_WORLD);
MPI_Recv(&(my_ary[0][COLS-1]), ROWS, newType, rank-1, pong, MPI_COMM_WORLD, &status);
ierr = MPI_Type_free(&newType);
cout<< "Final array for " << rank << "th processor is below "<<endl;
print_ary(my_ary);
}
