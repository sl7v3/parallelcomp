#include <cstdlib>
#include <mpi.h>
#include <stack>
#include <queue>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stddef.h>
#include "cq_dynamic.h"
#define work 101
#define max 001
#define worker_finish 000
#define MULTIPLIER 50
#define NBLOCKS 4
#define MAX_ARRAY 810
double f (double x);
const double eps = pow(10, -6);
const double s = 12;
double a = 1, b = 100;
void set_initial_queue_size(int size, int& initial_queue_size);
void get_current_max(double& curr_max, CQ::BoundaryVals current_interval);
void proc_dfs(CQ::BoundaryVals(&local_queue)[MAX_ARRAY], double& local_max, int arr_len, int rank);
void bfs_queue(CQ::Queue& interval_queue, double& local_max, int q_size, double& start);
void get_new_intervals(CQ::BoundaryVals working_range, CQ::BoundaryVals& first_new_interval, CQ::BoundaryVals& second_new_interval);
MPI_Request max_request;

using namespace std;
int main(int argc, char **argv)
{
  std::cout.precision(8);
  int initial_queue_size, rank, size, root = 0, total_sent = 0;
  const int root_array_definition = 16; //Setting this to the maximum array size this HW will deal with
  double local_max = 0, start = 0, time_elapsed = 0;
  CQ::BoundaryVals initial_range{ a, b, f(a), f(b) };
  CQ::Queue interval_queue;
  CQ::BoundaryVals  local_buffer[MAX_ARRAY], local_queue[MAX_ARRAY];
  double rcv_array[root_array_definition];
  MPI_Request send_request_array[root_array_definition], rcv_request_array[root_array_definition];
  MPI_Status send_status_array[root_array_definition], rcv_status_array[root_array_definition], max_status;
  bool time_to_count = true;
  int flag_array[root_array_definition] = { false };
  std::vector<int> working_processors;
  interval_queue.push(initial_range);   
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  set_initial_queue_size(size, initial_queue_size);
  int send_measure = (2) * size, length = initial_queue_size / send_measure, term_code = 0, global_term_code = 1;
  if (size == 1) { length = initial_queue_size; }
  MPI_Datatype boundaryValPacketType;
  MPI_Datatype array_of_types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
  int array_of_blocklengths[4] = { 1, 1, 1, 1 };
  MPI_Aint array_of_displacements[4];
  array_of_displacements[0] = offsetof(CQ::BoundaryVals, lower_bound);
  array_of_displacements[1] = offsetof(CQ::BoundaryVals, upper_bound);
  array_of_displacements[2] = offsetof(CQ::BoundaryVals, f_lower_bound);
  array_of_displacements[3] = offsetof(CQ::BoundaryVals, f_upper_bound);
  MPI_Type_create_struct(4, array_of_blocklengths, array_of_displacements, array_of_types, &boundaryValPacketType);
  MPI_Type_commit(&boundaryValPacketType);
  if (rank == root){
	  bfs_queue(interval_queue, local_max, initial_queue_size, start);	  
	  for (int i = 0; i < initial_queue_size; i++){
		  CQ::BoundaryVals rvPacket = interval_queue.pop();
		  local_buffer[i] = rvPacket;
	  }	
	  if (size == 1){
		 proc_dfs(local_buffer, local_max, length, rank);
		 time_elapsed = MPI_Wtime() - start;
		 cout << "The number of processors is " << size << " ,the global max is " << local_max << " and it took " << time_elapsed << " seconds." << endl;
		 MPI_Finalize();
		 return 0;
	  }
  }
  MPI_Bcast(&local_max, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Scatter(&local_buffer, length, boundaryValPacketType, &local_queue, length, boundaryValPacketType, root, MPI_COMM_WORLD);  
  if (rank == root){
	  total_sent += length * (size-1);
  }
  MPI_Barrier(MPI_COMM_WORLD); 
  if (rank == root){
	  for (int i = 1; i < size; i++){
		  MPI_Irecv(&rcv_array[i], 1, MPI_DOUBLE, i, worker_finish, MPI_COMM_WORLD, &rcv_request_array[i]);
	  }
  }
  proc_dfs(local_queue, local_max, length, rank);
  do {
	  if (rank == root){
		  std::vector<int> recently_sent_proc;
		  for (int i = 1; i < size; i++){
			  MPI_Test(&rcv_request_array[i], &flag_array[i], &rcv_status_array[i]);
			  if (flag_array[i] == true){
				  local_max = fmax(local_max, rcv_array[i]);
				  if (total_sent < initial_queue_size){
					  MPI_Send(&local_max, 1, MPI_DOUBLE, i, max, MPI_COMM_WORLD);
					  if (total_sent < (size*length)){
						  MPI_Isend(&local_buffer, length, boundaryValPacketType, i, work, MPI_COMM_WORLD, &send_request_array[i]);
					  }
					  else {
						  MPI_Isend(&local_buffer[total_sent], length, boundaryValPacketType, i, work, MPI_COMM_WORLD, &send_request_array[i]);
					  }					  
					  recently_sent_proc.push_back(i);
					  MPI_Irecv(&rcv_array[i], 1, MPI_DOUBLE, i, worker_finish, MPI_COMM_WORLD, &rcv_request_array[i]);
					  total_sent += length;
					  if ((total_sent >= initial_queue_size)){
						  i = 0;
					  }					  
					 }
				  else{
					    MPI_Send(&global_term_code, 1, MPI_INT, i, max, MPI_COMM_WORLD);
					}
			  }
			  else {
				  if (time_to_count && (total_sent >= initial_queue_size)){
					  working_processors.push_back(i);
				}
			  }
			  
		  }
		  //Waiting for all sent messages
		  int recently_sent_num = recently_sent_proc.size(), temp_rank = 0;
		  for (int i = 0; i < recently_sent_num; i++){
			  temp_rank = recently_sent_proc.at(i);
			  MPI_Wait(&send_request_array[temp_rank], &send_status_array[temp_rank]);			  
		  }
		  if (total_sent < initial_queue_size){
			  //if (size <= 2){
				  if (total_sent < (size*length)){
					  for (int i = 0; i < length; i++){
						  local_queue[i] = local_buffer[i];
					  }
				  }
				  else {
					  for (int i = total_sent; i < (total_sent + length); i++){
						  local_queue[i - total_sent] = local_buffer[total_sent];
					  }
				  }
				  total_sent += length;
			//}			  
		  }
		  else {
			  term_code = global_term_code;
		  }
		  if (time_to_count && (term_code == global_term_code)){
			  int workers_num = working_processors.size(), temp_rank = 0;
			  for (int i = 0; i < workers_num; i++){
				  temp_rank = working_processors.at(i);
				  MPI_Wait(&rcv_request_array[temp_rank], &rcv_status_array[temp_rank]);
				  local_max = fmax(local_max, rcv_array[temp_rank]);
				  MPI_Send(&global_term_code, 1, MPI_INT, temp_rank, max, MPI_COMM_WORLD);
			  }
			  time_to_count = false;
		  }
	  }	   
	  if (rank != root){
		  //Expect a msg, probe it and save it where it belongs.
		  int msg_size_count;
		  MPI_Probe(root, max, MPI_COMM_WORLD, &max_status);
		  MPI_Get_count(&max_status, MPI_INT, &msg_size_count);
		  if (msg_size_count != 1){
			 MPI_Recv(&local_max, 1, MPI_DOUBLE, root, max, MPI_COMM_WORLD, &max_status);
		  }
		  else {
			  MPI_Recv(&term_code, 1, MPI_INT, root, max, MPI_COMM_WORLD, &max_status);
			  time_to_count = false;		
		  }
		  if (term_code != global_term_code){
			  MPI_Recv(&local_queue, length, boundaryValPacketType, root, work, MPI_COMM_WORLD, &rcv_status_array[rank]);
		  }		  
	  }
	  if (term_code != global_term_code){			  
		 /* if (size > 2){
			  if (rank != root){
				  proc_dfs(local_queue, local_max, length, rank);
			  }
		  }
		  else {*/
			  proc_dfs(local_queue, local_max, length, rank);
		  //}
		  
	  }
  } while (time_to_count);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == root){
	  time_elapsed = MPI_Wtime() - start;
	  cout << "The number of processors is " << size << " ,the global max is " << local_max << " and it took " << time_elapsed << " seconds." << endl;
  }
  MPI_Finalize();
  return 0;
}

//Function that is all about nothing
double f(double x){
	double  sum_all = 0;
	for (int i = 100; i >= 1; i--){
		double jsum_all = 0;
		for (int j = i; j >= 1; j--){
			jsum_all += pow((x + (0.5*j)), (-3.3));
		}
		sum_all += (sin(x + jsum_all) / pow(1.3, i));
	}
	return sum_all;
}



double intervalMax(CQ::BoundaryVals range_object, double s){
	return (((range_object.f_lower_bound + range_object.f_upper_bound) + (s*(range_object.upper_bound - range_object.lower_bound))) / 2);
}

void proc_dfs(CQ::BoundaryVals(&local_queue)[MAX_ARRAY], double& local_max, int arr_len, int rank){
	local_max = fmax(local_max, fmax(local_queue[0].f_lower_bound, local_queue[0].f_upper_bound));
	for (int arr_front = 0; arr_front < arr_len; arr_front++){
		std::stack<CQ::BoundaryVals> proc_stack;
		CQ::BoundaryVals stack_interval;
		bool stack_interval_available;
		stack_interval = local_queue[arr_front];
		do {
			stack_interval_available = false;
			get_current_max(local_max, stack_interval);
			if (intervalMax(stack_interval, s) >= (local_max + eps)){
				CQ::BoundaryVals first_new_interval, second_new_interval;
				get_new_intervals(stack_interval, first_new_interval, second_new_interval);
				proc_stack.push(second_new_interval);
				stack_interval = first_new_interval;
				stack_interval_available = true;
			}
			else {
				if (!proc_stack.empty()){
					stack_interval = proc_stack.top();
					proc_stack.pop();
					stack_interval_available = true;
				}
			}
			
		} while (stack_interval_available);
	}
	if (rank != 0){
		MPI_Isend(&local_max, 1, MPI_DOUBLE, 0, worker_finish, MPI_COMM_WORLD, &max_request);
	}
}

void get_current_max(double& curr_max, CQ::BoundaryVals current_interval){
	//double new_max = f(current_interval.lower_bound + current_interval.upper_bound) / 2;
	//double new_max = fmax(current_interval.f_lower_bound, current_interval.f_upper_bound);
	double new_max = fmax(current_interval.f_lower_bound, current_interval.f_upper_bound);
	curr_max = fmax(curr_max, new_max);
}

void get_new_intervals(CQ::BoundaryVals working_range, CQ::BoundaryVals& first_new_interval, CQ::BoundaryVals& second_new_interval){
	double new_bound = (working_range.lower_bound + working_range.upper_bound) / 2;
	double f_new_bound = f(new_bound);
	CQ::BoundaryVals temp_first_new_interval{ working_range.lower_bound, new_bound, working_range.f_lower_bound, f_new_bound };
	CQ::BoundaryVals temp_second_new_interval{ new_bound, working_range.upper_bound, f_new_bound, working_range.f_upper_bound };
	{
		first_new_interval = temp_first_new_interval;
		second_new_interval = temp_second_new_interval;
	}
}

void set_initial_queue_size(int size, int& initial_queue_size){
	initial_queue_size = size * MULTIPLIER;
}

void bfs_queue(CQ::Queue& interval_queue, double& local_max, int q_size, double& start){
	CQ::BoundaryVals first_new_interval, second_new_interval;
	bool working_range_available;
	CQ::BoundaryVals working_range = interval_queue.pop();
	do {
		get_current_max(local_max, working_range); //assign max to maximum of current_max and the boundaries of the interval
		if (intervalMax(working_range, s) >= (local_max + eps)){
			if ((working_range.lower_bound == a) && (working_range.upper_bound == b)) {
				start = MPI_Wtime(); //start the clock if this is the first time intervals are being divided
			}
			get_new_intervals(working_range, first_new_interval, second_new_interval); //gets new intervals (high and low) from current interval
			interval_queue.push(first_new_interval);
			interval_queue.push(second_new_interval);
		}

		if (interval_queue.get_num_items() < q_size){
			if (!interval_queue.empty()){
				working_range = interval_queue.pop(); //pop the front item in queue and repeat
				working_range_available = true;
			}
			else {
				working_range_available = false;
			}
		}
		else {
			working_range_available = false; //this is set to false to signal that we are done creating the queue
		}
	} while (working_range_available);
}