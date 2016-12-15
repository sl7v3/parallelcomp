#include <iostream>
#include <math.h>
#include <time.h>
# include <omp.h>
#include "cq_dynamic.h"
#include <stack>
using namespace std;
#define NUM_THREADS 1
const double eps = pow(10, -6);
const double s = 12;
double f (double x);
double a = 1, b = 100;
omp_lock_t my_lock;
double intervalMax(CQ::BoundaryVals range_object, double s);
void set_initial_queue_size(int& initial_queue_size);
int num_threads = int(NUM_THREADS);
void get_current_global_max(double& max, CQ::BoundaryVals current_interval);
void get_current_global_max(double& max, CQ::BoundaryVals current_interval, omp_lock_t &my_lock);
void thread_dfs(CQ::Queue& interval_queue, double& max, omp_lock_t &my_lock, int stack_size);
void single_thread_dfs(CQ::BoundaryVals initial_range, double& max, time_t& timer);
void bfs_queue(CQ::Queue& interval_queue, double& max, int q_size, time_t& timer);
void spawn_all_threads(CQ::Queue& interval_queue, double& max, omp_lock_t &my_lock, int initial_queue_size, int stack_size);
void get_new_intervals(CQ::BoundaryVals working_range, CQ::BoundaryVals& first_new_interval, CQ::BoundaryVals& second_new_interval);
int main()
{
  std::cout.precision(8);
  CQ::BoundaryVals initial_range { a, b, f(a), f(b)};
  CQ::BoundaryVals working_range = initial_range;
  omp_init_lock(&my_lock);
  int initial_queue_size, stack_size = 4;
  if (num_threads == 1){ stack_size += 2; };
  set_initial_queue_size(initial_queue_size);
  double max = fmax (initial_range.f_lower_bound, initial_range.f_upper_bound);
  std::cout << "Number of threads is " << num_threads << endl; 
  std::cout << "The initial max value is " << (max) << " and the boundary intervals are " << working_range.lower_bound << " and "  << working_range.upper_bound << endl;
  CQ::Queue interval_queue;
  interval_queue.push(initial_range);
  time_t timer;
  	  bfs_queue(interval_queue, max, initial_queue_size, timer); //queue of desired size created
	  cout << "Initial queue size is " << initial_queue_size << endl;
	  do {
		  spawn_all_threads(interval_queue, max, my_lock, initial_queue_size, stack_size); //now spawn parallel 
		  if (interval_queue.get_num_items() > 0){
			  if (num_threads == 1){
				 bfs_queue(interval_queue, max, initial_queue_size, timer);
			  }
			  else {
				  if (interval_queue.get_num_items() < num_threads){
					  bfs_queue(interval_queue, max, num_threads, timer);
					  initial_queue_size = num_threads;
				  }
			  }
		  }
	} while (interval_queue.get_num_items() > 0);  
    double seconds = omp_get_wtime() - timer;
  omp_destroy_lock(&my_lock);  
  std::cout << "The final max value is " << (max) << " and it took " << fixed << seconds << " seconds." << endl;
  system("pause");
  return 0;
}
double intervalMax(CQ::BoundaryVals range_object, double s){
  return (((range_object.f_lower_bound +range_object.f_upper_bound) + (s*(range_object.upper_bound - range_object.lower_bound)))/2);
}
double f (double x){
  double  sum_all = 0;
  for (int i = 100; i>=1; i--){
    double jsum_all = 0;
    for (int j = i; j >= 1; j--){
        jsum_all += pow((x + (0.5*j)),(-3.3));
      }
	  sum_all += (sin(x + jsum_all) / pow(1.3, i));
    }
  return sum_all;
}

void thread_dfs(CQ::Queue& interval_queue, double& max, omp_lock_t &my_lock, int stack_limit){
  std::stack<CQ::BoundaryVals> thread_stack;
  CQ::BoundaryVals stack_interval;
  bool stack_interval_available;
  #pragma omp critical(queueupdate)
  {
	  stack_interval = interval_queue.pop();
  }   
  do {
	      stack_interval_available = false;
	      get_current_global_max(max, stack_interval, my_lock);
		  if (intervalMax(stack_interval, s) >= (max + eps)){
			  CQ::BoundaryVals first_new_interval, second_new_interval;
			  get_new_intervals(stack_interval, first_new_interval, second_new_interval);
			  thread_stack.push(second_new_interval);
			  
				  if (thread_stack.size() >= stack_limit){
                  #pragma omp critical(queueupdate)
				{
					interval_queue.push(first_new_interval);
					int stack_size = thread_stack.size();
					for (int i = 0; i <stack_size; i++){
						stack_interval = thread_stack.top();
						thread_stack.pop();
						interval_queue.push(stack_interval); 
					}
				}
				  }
				  else {
					  stack_interval = first_new_interval;
					  stack_interval_available = true;
				  }
			  
		 } 
		  else {
			  if (!thread_stack.empty()){
				  stack_interval = thread_stack.top();
				  thread_stack.pop();
				  stack_interval_available = true;
			  }
		  }
     } while (stack_interval_available);
}
void spawn_all_threads(CQ::Queue& interval_queue, double& max, omp_lock_t &my_lock, int initial_queue_size, int stack_size){
	//stop enqueueing, get each thread to pop from queue dynamically
	//compute ranges using a stack for each thread
	#pragma omp parallel for schedule(dynamic)   
	for (int n = 0; n < int(initial_queue_size); ++n){
		thread_dfs(interval_queue, max, my_lock, stack_size);
	}
  
}

void get_current_global_max(double& max, CQ::BoundaryVals current_interval){
	double new_max = fmax(current_interval.f_lower_bound, current_interval.f_upper_bound);
	max = fmax(max, new_max);
}
void get_current_global_max(double& max, CQ::BoundaryVals current_interval, omp_lock_t &my_lock){
	double new_max = fmax(current_interval.f_lower_bound, current_interval.f_upper_bound);
	if (new_max > max){
		omp_set_lock(&my_lock);
		{
			if (new_max > max){
				max = new_max;
			}
		}
		omp_unset_lock(&my_lock);
	}
}
	

void get_new_intervals(CQ::BoundaryVals working_range, CQ::BoundaryVals& first_new_interval, CQ::BoundaryVals& second_new_interval){
  double new_bound = (working_range.lower_bound + working_range.upper_bound)/2;
  double f_new_bound = f(new_bound);
  CQ::BoundaryVals temp_first_new_interval { working_range.lower_bound, new_bound, working_range.f_lower_bound, f_new_bound}; 
  CQ::BoundaryVals temp_second_new_interval { new_bound, working_range.upper_bound,  f_new_bound, working_range.f_upper_bound};
  first_new_interval = temp_first_new_interval;
  second_new_interval = temp_second_new_interval;
}

void set_initial_queue_size(int& initial_queue_size){
	int multiplier = 1;
	if (NUM_THREADS == 1){
		multiplier = 56;
	}
	if (NUM_THREADS == 2){
		multiplier = 28;
	}
	else if (NUM_THREADS == 4){
		multiplier = 14;
	}
	else if (NUM_THREADS == 16){
		multiplier = 4;
	}
	else if (NUM_THREADS == 28){
		multiplier = 2;
	}
	initial_queue_size = int(NUM_THREADS) * multiplier;
}

void bfs_queue(CQ::Queue& interval_queue, double& max, int q_size, time_t& timer){
	CQ::BoundaryVals first_new_interval, second_new_interval;
	bool working_range_available;
	CQ::BoundaryVals working_range = interval_queue.pop();
	do {
		get_current_global_max(max, working_range); //assign max to maximum of current_max and the boundaries of the interval
		if (intervalMax(working_range, s) >= (max + eps)){
			if ((working_range.lower_bound == a) && (working_range.upper_bound == b)) {
				timer = omp_get_wtime(); //start the clock if this is the first time intervals are being divided
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


