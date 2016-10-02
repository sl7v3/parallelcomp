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
struct packet {
   int intRank; 
   float floatRank;
};
void processor_A (void);
void processor_B (void);
void processor_C (void);
using namespace std;
MPI_Datatype rankPacketType;
int main(int argc, char **argv)
{
    int rank, size, len;
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
int rank, size, sum = 0;
float floatSum = 0;
MPI_Datatype rankPacketType;
MPI_Datatype array_of_types[NBLOCKS] = {MPI_INT, MPI_FLOAT};
int array_of_blocklengths[NBLOCKS] = {1,1};
MPI_Aint array_of_displacements[NBLOCKS], address1, address2;
packet rankPacket;
MPI_Get_address(&rankPacket, &address1);
MPI_Get_address(&rankPacket.floatRank, &address2);
array_of_displacements[1] = address2 - address1;
array_of_displacements[0] = 0;
MPI_Type_create_struct (2, array_of_blocklengths, array_of_displacements, array_of_types, &rankPacketType);
MPI_Type_commit (&rankPacketType);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
rankPacket.intRank = rank;
rankPacket.floatRank = (float)rank*2.150;
 
int i, length = 1;
  for (int sends =0; sends< size; sends++){
    MPI_Send(&rankPacket, length, rankPacketType, rank+1, ping, MPI_COMM_WORLD);
    MPI_Recv(&rankPacket, length, rankPacketType, size-1, pong, MPI_COMM_WORLD, &status);
    sum+= rankPacket.intRank;
    floatSum+= rankPacket.floatRank;
  }
  cout << "Float sum of all ranks received by processor at rank " <<rank<< " is "  << floatSum << endl;
  cout << "Integer sum of all ranks received by processor at rank " <<rank<< " is "  << sum << endl;
}
void processor_B ( void){
MPI_Status status;
int rank, size, sum = 0;
float floatSum = 0;
MPI_Datatype rankPacketType;
MPI_Datatype array_of_types[NBLOCKS] = {MPI_INT, MPI_FLOAT};
int array_of_blocklengths[NBLOCKS] = {1,1};
MPI_Aint array_of_displacements[NBLOCKS], address1, address2;
packet rankPacket;
MPI_Get_address(&rankPacket, &address1);
MPI_Get_address(&rankPacket.floatRank, &address2);
array_of_displacements[1] = address2 - address1;
array_of_displacements[0] = 0;
MPI_Type_create_struct (2, array_of_blocklengths, array_of_displacements, array_of_types, &rankPacketType);
MPI_Type_commit (&rankPacketType);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
rankPacket.intRank = rank;
rankPacket.floatRank = (float)rank*2.150;
int i, length = 1;
for (int sends =0; sends< size; sends++){
    MPI_Send(&rankPacket, length, rankPacketType, rank+1, ping, MPI_COMM_WORLD);
    MPI_Recv(&rankPacket, length, rankPacketType, rank-1, pong, MPI_COMM_WORLD, &status);
    sum+= rankPacket.intRank;
    floatSum+= rankPacket.floatRank;
  }
  cout << "Float sum of all ranks received by processor at rank " <<rank<< " is "  << floatSum << endl;
  cout << "Integer sum of all ranks received by processor at rank " <<rank<< " is "  << sum << endl;
}

void processor_C ( void){
MPI_Status status;
int rank, size, sum = 0;
float floatSum = 0;
MPI_Datatype rankPacketType;
MPI_Datatype array_of_types[NBLOCKS] = {MPI_INT, MPI_FLOAT};
int array_of_blocklengths[NBLOCKS] = {1,1};
MPI_Aint array_of_displacements[NBLOCKS], address1, address2;
packet rankPacket;
MPI_Get_address(&rankPacket, &address1);
MPI_Get_address(&rankPacket.floatRank, &address2);
array_of_displacements[1] = address2 - address1;
array_of_displacements[0] = 0;
MPI_Type_create_struct (2, array_of_blocklengths, array_of_displacements, array_of_types, &rankPacketType);
MPI_Type_commit (&rankPacketType);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
rankPacket.intRank = rank;
rankPacket.floatRank = (float)rank*2.150;
int i, length = 1;
for (int sends =0; sends< size; sends++){
    MPI_Send(&rankPacket, length, rankPacketType, 0, ping, MPI_COMM_WORLD);
    MPI_Recv(&rankPacket, length, rankPacketType, rank-1, pong, MPI_COMM_WORLD, &status);    
    sum+= rankPacket.intRank;
    floatSum+= rankPacket.floatRank;
  }
  cout << "Float sum of all ranks received by processor at rank " <<rank<< " is "  << floatSum << endl;
  cout << "Integer sum of all ranks received by processor at rank " <<rank<< " is "  << sum << endl;
}

