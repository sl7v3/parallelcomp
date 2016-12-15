#include <cstdlib>
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#define ping 101
#define pong 101
#define ROWS 2000
#define COLS 500
void setup_dimensions (int rank, int size, int sqrt_size, int& proc_rows, int& proc_cols, int row_max_per_proc, int col_max_per_proc, int& offset_row, int& offset_col, int& boundary_row, int& boundary_col, int& start_col, int& start_row, int& row_span, int& col_span );
double f (double x);
using namespace std;
int main(int argc, char **argv)
{
  int rank, size, proc_rows, proc_cols, row_max_per_proc, col_max_per_proc, len, row_span, col_span, root =0, offset_col, offset_row, start_col, start_row, boundary_row, boundary_col, sqrt_size;
  double local_sum = 0, local_sum_of_squares = 0, global_sum_of_squares = 0,  global_sum =0, start = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  sqrt_size = (int)(sqrt(size));
  MPI_Status status;
  col_max_per_proc = ceil((double)COLS/sqrt_size);
  row_max_per_proc = ceil((double)ROWS/sqrt_size);
  setup_dimensions (rank, size, sqrt_size, proc_rows, proc_cols, row_max_per_proc, col_max_per_proc, offset_row, offset_col, boundary_row, boundary_col, start_col, start_row, row_span, col_span);
  double arr[proc_rows][proc_cols], arrB[proc_rows][proc_cols];
  MPI_Datatype rowType, colType;
  MPI_Type_vector(proc_rows, 1, proc_cols, MPI_DOUBLE, &colType);MPI_Type_commit(&colType);
  MPI_Type_vector(1, proc_cols, 0, MPI_DOUBLE, &rowType);MPI_Type_commit(&rowType);
  
//Initialize arrays based on processor rank
  for (int row_count = offset_row; row_count < row_span; row_count++){
    start_col = (rank%sqrt_size) * col_max_per_proc;
    for (int col_count = offset_col; col_count < col_span; col_count++){
      arr[row_count][col_count] = start_row*sin(start_row) + start_col*cos(start_col) + sqrt(start_row+start_col);
      start_col++;
    }
    start_row++;
  }
 //Make sure all processors have gotten to this point
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == root){
    start = MPI_Wtime();
  }
//Iterative step
  for (int iter = 0; iter < 10; iter++){
    start_col = (rank%sqrt_size) * col_max_per_proc;
    start_row = (rank/sqrt_size) * row_max_per_proc;
    if(size>1){
  	  if (!(start_row== 0)){
  	    MPI_Send((&arr[1][0]), 1, rowType, rank-sqrt_size, ping, MPI_COMM_WORLD);
  	  }

      if (!(boundary_row== ROWS)){
        MPI_Recv((&arr[proc_rows-1][0]), 1, rowType,  rank+sqrt_size, pong, MPI_COMM_WORLD, &status);
      }
      if (!(start_row== ((sqrt_size-1)*row_max_per_proc))){
  	    MPI_Send((&arr[proc_rows -2][0]), 1, rowType, rank+sqrt_size, ping, MPI_COMM_WORLD);
  	  }
      if (!(start_row== 0)){
        MPI_Recv((&arr[0][0]), 1, rowType, rank-sqrt_size, pong, MPI_COMM_WORLD, &status);
      }

      if (!(start_col== 0)){
  	    MPI_Send((&arr[0][1]), 1, colType, rank-1, ping, MPI_COMM_WORLD);
  	  }
      
      if (!(boundary_col== COLS)){
        MPI_Recv((&arr[0][proc_cols-1]), 1, colType, rank+1, pong, MPI_COMM_WORLD, &status);
      }

      if (!(start_col== ((sqrt_size-1)*col_max_per_proc))){
  	    MPI_Send((&arr[0][proc_cols-2]), 1, colType, rank+1, ping, MPI_COMM_WORLD);
  	  }
  	 if (!(start_col== 0)){
  	    MPI_Recv((&arr[0][0]), 1, colType, rank-1, pong, MPI_COMM_WORLD, &status);
  	  }
  }
    start_row = (rank/sqrt_size) * row_max_per_proc;
    std::copy(&arr[0][0], &arr[0][0]+proc_rows*proc_cols,&arrB[0][0]); 
    for (int row_count = offset_row; row_count < row_span; row_count++){
      start_col = (rank%sqrt_size) * col_max_per_proc;
	    for (int col_count = offset_col; col_count < col_span; col_count++){
		     if (!(start_col == 0 || start_col == COLS-1 || start_row == 0 || start_row == ROWS-1 )){ 
            arrB[row_count][col_count] = f(arr[row_count][col_count]);
			    }
	      start_col++;
	    }
      start_row++;
    }
    start_row = (rank/sqrt_size) * row_max_per_proc;
    for (int row_count = offset_row; row_count < row_span; row_count++){
      start_col = (rank%sqrt_size) * col_max_per_proc;
      for (int col_count = offset_col; col_count < col_span; col_count++){
         if (!(start_col == 0 || start_col == COLS-1 || start_row == 0 || start_row == ROWS-1 )){ 
            double z = arrB[row_count-1][col_count] + arrB[row_count+1][col_count]+ arrB[row_count][col_count-1] + arrB[row_count][col_count+1] + arrB[row_count][col_count];
            z/=5;
            arr[row_count][col_count] = fmax (double(-100), fmin(double(100),z));                
          }
        start_col++;
      }
      start_row++;
    }
  }

  //Calculate local sum for each processor.
   for (int row_count = offset_row; row_count < row_span; row_count++){
	    for (int col_count = offset_col; col_count < col_span; col_count++){
	      local_sum += arr[row_count][col_count];
	      local_sum_of_squares += pow(arr[row_count][col_count],2);
	    }
   }
  //Reduce the sum into the root processor and print time taken for entire process
  if (size > 1){
     MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&local_sum_of_squares, &global_sum_of_squares, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
     global_sum = local_sum;
     global_sum_of_squares = local_sum_of_squares;
  }
 if (rank ==root) {
    cout << "Size " << size << ", root rank " << rank << " .Time elapsed:  " << MPI_Wtime()-start << " seconds, global sum: " << global_sum << ", global sum of squares: " << global_sum_of_squares << endl;
  }
  MPI_Finalize();
  return 0;
}
//Function that is all about nothing
double f (double x){
  int i;
  double y = x;
  for (i = 1; i<11; i++){
    y = y + sin(x*i)/pow(2.0,i);
  }
  return y;
}
//Function that computes dimensions of arrays in each processor
void setup_dimensions (int rank, int size, int sqrt_size, int& proc_rows, int& proc_cols, int row_max_per_proc, int col_max_per_proc, int& offset_row, int& offset_col, int& boundary_row, int& boundary_col, int& start_col, int& start_row, int& row_span, int& col_span ){
  offset_col = 0;
  offset_row = 0;
  start_col = (rank%sqrt_size) * col_max_per_proc;
  start_row = (rank/sqrt_size) * row_max_per_proc;
  if (!(start_row== ((sqrt_size-1)*row_max_per_proc))){
    boundary_row = ((rank/sqrt_size)+1)*row_max_per_proc;
  } else boundary_row = ROWS;
  
  if (!(start_col== ((sqrt_size-1)*col_max_per_proc))){
    boundary_col = ((rank%sqrt_size)+1)*col_max_per_proc;
  } else boundary_col  = COLS;
  proc_rows = row_span = boundary_row - start_row;
  proc_cols = col_span = boundary_col - start_col;
  if (size >1) {
    if (!(rank%(sqrt_size)==0)){ proc_cols+=1; offset_col+=1; col_span+=1;};
    if (!(rank%(sqrt_size)==(sqrt_size-1))){ proc_cols+=1;};
    if (!(rank/(sqrt_size)==0)){ proc_rows+=1; offset_row+=1; row_span+=1;};
    if (!(rank/(sqrt_size)==(sqrt_size-1))){ proc_rows+=1;};
  }
}

