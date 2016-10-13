#include <cstdlib>
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#define ping 101
#define pong 101
#define ROWS 100
#define COLS 40
void setup_dimensions (int rank, int size, int& proc_rows, int& proc_cols, int& row_col_max_per_proc );
double f (double x);
using namespace std;
int main(int argc, char **argv)
{
  int rank, size, proc_rows, proc_cols, row_col_max_per_proc, len, row_span, col_span, root =0, offset_col, offset_row, start_col, start_row, boundary_row, boundary_col;
  double local_sum = 0, local_sum_of_squares = 0, global_sum_of_squares = 0,  global_sum =0, start, time_elapsed;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;
  if (size ==1){
    proc_cols = COLS;
    proc_rows = ROWS;
    double arr[proc_rows][proc_cols], arrB[proc_rows][proc_cols];
    for (int row = 0; row < proc_rows; row++){
       for (int col = 0; col < proc_cols; col++){
	      arr[row][col] = row*sin(row) + col*cos(col) + sqrt(row+col);
       }
    }
    start = MPI_Wtime();
    for (int i = 0; i<10; i++){
    std::copy(&arr[0][0], &arr[0][0]+proc_rows*proc_cols,&arrB[0][0]); 
      for (int row = 0; row < proc_rows; row++){
	    for (int col =0; col < proc_cols; col++){
	      if (!(col == 0 || col == proc_cols-1 || row == 0 || row == proc_rows-1 )){ 
                 double z = f(arrB[row-1][col]) + f(arrB[row+1][col])+ f(arrB[row][col-1]) + f(arrB[row][col+1]) + f(arrB[row][col]);
                 z/=5;
                 arr[row][col] =  fmax (double(-100), fmin(double(100),z));                
	      }
	    }
      }
    }
    for (int r = 0; r < proc_rows; r++){
        for (int c = 0; c<proc_cols; c++){
          local_sum += arr[r][c];
          local_sum_of_squares += pow(arr[r][c],2);
        }
     }
    time_elapsed = MPI_Wtime() - start;
    cout << "The parallel time elapsed for printing and calculating sums and sums of squares is " << time_elapsed << " seconds." << endl;
    cout << "The global sum printed by rank " << rank << " is " << local_sum << endl;
    cout << "The global sum of squares printed by rank " << rank << " is " << local_sum_of_squares << endl;
  
  } else {
          setup_dimensions(rank, size, proc_rows, proc_cols, row_col_max_per_proc);
          double arr[proc_rows][proc_cols], arrB[proc_rows][proc_cols];
          MPI_Datatype rowType, colType;
	  //Initialize arrays based on processor rank and row vs col comparison
	  if (ROWS>COLS){
	     MPI_Type_vector(proc_rows, 1, proc_cols, MPI_DOUBLE, &colType);MPI_Type_commit(&colType);
		  if (rank < size -1){
		    boundary_col = ((rank+1)*row_col_max_per_proc);
		  } else boundary_col = COLS;
		  if (rank > 0){
		    offset_col = 1;
		  } else offset_col = 0;
		  start_col  = rank*row_col_max_per_proc; 
		  col_span = (boundary_col - start_col) + offset_col;
		  for (int row = 0; row < proc_rows; row++){
		    start_col  = rank*row_col_max_per_proc; 
		    for (int col_count = offset_col; col_count < col_span; col_count++){
		      arr[row][col_count] = row*sin(row) + start_col*cos(start_col) + sqrt(row+start_col);
		      start_col++;
		    }
		  }
          }
	  else {
	        MPI_Type_vector(1, COLS, 0, MPI_DOUBLE, &rowType);MPI_Type_commit(&rowType);
		if (rank < size -1){
		  boundary_row = ((rank+1)*row_col_max_per_proc);
		} else boundary_row = ROWS;
		if (rank > 0){
		  offset_row = 1;
		} else offset_row = 0;
		start_row  = rank*row_col_max_per_proc;
		row_span = (boundary_row - start_row) + offset_row;
		for (offset_row; offset_row < row_span; offset_row++){
		  for (int col = 0; col < proc_cols; col++){
		    arr[offset_row][col] = start_row*sin(start_row) + col*cos(col) + sqrt(start_row+col);
		  }
		  start_row++;
		}
	  }
	//Make sure all processors have gotten to this point
cout << "The value in the 5th row and 21st column is "<< arr[4][20]<<endl;
	  MPI_Barrier(MPI_COMM_WORLD);
	  if (rank == root){
	    start = MPI_Wtime();
	  }
	//Iterative step - prepare a MPI_Vector to send or receive from neighboring processors
	  for (int iter = 0; iter < 10; iter++){
	   if (ROWS>COLS){
	     if (rank > 0 && rank < size -1){
	       MPI_Send((&arr[0][1]),1, colType, rank-1, ping, MPI_COMM_WORLD);
	       MPI_Send((&arr[0][proc_cols-2]),1, colType, rank+1, ping, MPI_COMM_WORLD);
	       MPI_Recv((&arr[0][0]),1, colType, rank-1, pong, MPI_COMM_WORLD, &status);
	       MPI_Recv((&arr[0][proc_cols-1]),1, colType, rank+1, pong, MPI_COMM_WORLD, &status);
	     }
	     else if (rank == 0){
	       MPI_Send((&arr[0][proc_cols-2]),1, colType, rank+1, ping, MPI_COMM_WORLD);
	       MPI_Recv((&arr[0][proc_cols-1]),1, colType, rank+1, pong, MPI_COMM_WORLD, &status);
	     }
	     else {
	       MPI_Send((&arr[0][1]),1, colType, rank-1, ping, MPI_COMM_WORLD);
	       MPI_Recv((&arr[0][0]),1, colType, rank-1, pong, MPI_COMM_WORLD, &status);
	     }
	     std::copy(&arr[0][0], &arr[0][0]+proc_rows*proc_cols,&arrB[0][0]); 
		  start_col  = rank*row_col_max_per_proc; 
		  for (int row = 0; row < proc_rows; row++){
		    start_col  = rank*row_col_max_per_proc; 
		    for (int col_count = offset_col; col_count < col_span; col_count++){
		      if (!(start_col == 0 || start_col == COLS-1 || row == 0 || row == ROWS-1 )){ 
			 double z = f(arrB[row-1][col_count]) + f(arrB[row+1][col_count])+ f(arrB[row][col_count-1]) + f(arrB[row][col_count+1]) + f(arrB[row][col_count]);
			 z/=5;
			 arr[row][col_count] = fmax (double(-100), fmin(double(100),z));                
		      }
		      start_col++;
		    }
	          } 
         } 
	   else {
	      if (rank > 0 && rank < size -1){
		MPI_Send((&arr[1][0]),len, rowType, rank-1, ping, MPI_COMM_WORLD);
		MPI_Send((&arr[proc_rows -2][0]),len, rowType, rank+1, ping, MPI_COMM_WORLD);
		MPI_Recv((&arr[0][0]),len, rowType, rank-1, pong, MPI_COMM_WORLD, &status);
		MPI_Recv((&arr[proc_rows -1][0]),len, rowType, rank+1, pong, MPI_COMM_WORLD, &status);
	      }
	      else if (rank == 0){
		MPI_Send((&arr[proc_rows -2][0]),len, rowType, rank+1, ping, MPI_COMM_WORLD);
		MPI_Recv((&arr[proc_rows -1][0]),len, rowType, rank+1, pong, MPI_COMM_WORLD, &status);
	      }
	      else {
		MPI_Send((&arr[1][0]),len, rowType, rank-1, ping, MPI_COMM_WORLD);
		MPI_Recv((&arr[0][0]),len, rowType, rank-1, pong, MPI_COMM_WORLD, &status);
	      }  
	      std::copy(&arr[0][0], &arr[0][0]+proc_rows*proc_cols,&arrB[0][0]); 
              if (rank > 0){
	        offset_row = 1;
	      } else offset_row = 0;
	      start_row  = rank*row_col_max_per_proc;
		for (offset_row; offset_row < row_span; offset_row++){
		  for (int col = 0; col < proc_cols; col++){
		    if (!(start_row == 0 || start_row == ROWS-1 || col == 0 || col == COLS-1 )){ 
			 double z = f(arrB[offset_row-1][col]) + f(arrB[offset_row+1][col])+ f(arrB[offset_row][col-1]) + f(arrB[offset_row][col+1]) + f(arrB[offset_row][col]);
			 z/=5;
			 arr[offset_row][col] = fmax (double(-100), fmin(double(100),z));                
		    }
		  }
		 start_row++;
	      }
	  }
      }
	  //Calculate local sum for each processor depending on the strip that each processor holds.
	  if (ROWS > COLS){
	    for (int row = 0; row < proc_rows; row++){
	      int col_count = offset_col;
	      for (col_count; col_count < col_span; col_count++){
		local_sum += arr[row][col_count];
		local_sum_of_squares += pow(arr[row][col_count],2);
	      }
	    }
	  }
	  else {
	    if (rank > 0){
	      offset_row = 1;
	    } else offset_row = 0;
	   for (offset_row; offset_row < row_span; offset_row++){
	      for (int col = 0; col < proc_cols; col++){
		local_sum += arr[offset_row][col];
		local_sum_of_squares += pow(arr[offset_row][col],2);
	      }
	    }
	  }
	  //Reduce the sum into the root processor and print time taken for entire process
	  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&local_sum_of_squares, &global_sum_of_squares, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  if (rank ==root) {
	    time_elapsed = MPI_Wtime() - start;
	    cout << "The parallel time elapsed for printing and calculating local sums and local sums of squares is " << time_elapsed << " seconds." << endl;
	    cout << "The global sum printed by rank " << rank << " is " << global_sum << endl;
	    cout << "The global sum of squares printed by rank " << rank << " is " << global_sum_of_squares << endl;
	  }
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
void setup_dimensions (int rank, int size, int& proc_rows, int& proc_cols, int& row_col_max_per_proc ){
  if (ROWS>COLS){
     row_col_max_per_proc = ceil((double)COLS/size);
     proc_rows = ROWS;
     if (rank < size -1)
     {
        proc_cols = row_col_max_per_proc;
     }
     else proc_cols = COLS - ((size-1)*row_col_max_per_proc);
     if (rank >0 && rank < size-1){
       proc_cols+=2;
     } else proc_cols+=1;
  }
  else {
     row_col_max_per_proc = ceil((double)ROWS/size);
     proc_cols = COLS;
     if (rank < size -1)
     {  
       proc_rows = row_col_max_per_proc;
     }
     else proc_rows = ROWS - ((size-1)*row_col_max_per_proc);
     if (rank >0 && rank < size-1){
       proc_rows+=2;
     } else proc_rows+=1;
 
  }
}

