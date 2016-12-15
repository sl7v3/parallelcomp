#ifdef _OPENMP
# include <omp.h>
#endif
namespace CQ{
  struct BoundaryVals
  { 
    double lower_bound, upper_bound, f_lower_bound, f_upper_bound;
    /*void operator=(const BoundaryVals &B )
      { 
        lower_bound = B.lower_bound;
        upper_bound = B.upper_bound;
        f_lower_bound = B.f_lower_bound;
        f_upper_bound = B.f_upper_bound;
      }*/
  };


	class Queue {
		public:
		Queue();
		void push(const BoundaryVals& item);
		BoundaryVals pop();
		int get_num_items();
		int get_capacity();
		bool empty();
    void reallocate();
		private:
		static const size_t DEFAULT_CAPACITY = 800000;
    size_t capacity;
    size_t num_items;
    size_t front_index;
    size_t rear_index;
    BoundaryVals* the_data;
	};
	Queue::Queue(){
	num_items = 0;
	front_index = 0;
	rear_index = -1;
	capacity= DEFAULT_CAPACITY;
  the_data = new BoundaryVals[DEFAULT_CAPACITY];
    
	}
	void CQ::Queue::push(const BoundaryVals& item) {
    if (num_items == capacity) {
      reallocate();
    }
		  num_items++;
		  rear_index = (rear_index + 1) % capacity;
		  the_data[rear_index] = item;
	}
	CQ::BoundaryVals CQ::Queue::pop() {
		BoundaryVals front_interval = the_data[front_index];
		front_index = (front_index + 1) % capacity;
		num_items--;
		return front_interval;
	}
	int CQ::Queue::get_num_items() {
		return num_items;
	}
	int CQ::Queue::get_capacity() {
		return capacity;
	}
	bool Queue::empty() {
	  return (num_items == 0);
	}

  void Queue::reallocate() {
	  std::cout << "Started reallocation " << std::endl;
    size_t new_capacity = 2 * capacity;
    BoundaryVals* new_data = new BoundaryVals[new_capacity];
    BoundaryVals* temp_data;
    size_t j = front_index;
    for (size_t i = 0; i < num_items; i++) {
      new_data[i] = the_data[j];
      j = (j + 1) % capacity;
    }
    front_index = 0;
    rear_index = num_items - 1;
    capacity = new_capacity;
    temp_data = the_data;
    the_data = new_data;
    delete[] temp_data;
	std::cout << "Queue size is " << new_capacity << std::endl;
  }

 
}
