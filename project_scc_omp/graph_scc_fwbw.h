#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <omp.h>
#include <time.h>
#include <utility>
#include <algorithm>
using namespace std;


 
// A structure to represent an adjacency list node
struct AdjListNode
{
    int dest;
    struct AdjListNode* next;
    bool valid_coloring_root;
    bool valid_for_coloring;
};
 
// A structure to represent an adjacency list
struct AdjList
{
    int node_id;
    bool list_used;
    struct AdjListNode* head;  // pointer to head node of list
};
 
// A structure to represent a graph. A graph is an array of adjacency lists.
// Size of array will be V (number of vertices in graph)
struct Graph
{
    int vertices;
    int edges;
    std::vector<AdjList*> array;
    std::vector <int> colors;
    std::vector <bool> unique_color;
    std::vector <bool> valid_vertex;
    ~Graph (){;}
};
 
// A utility function to create a new adjacency list node
struct AdjListNode* newAdjListNode(int dest)
{
    struct AdjListNode* newNode =
            (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}
 
// A utility function that creates a graph of V vertices
/*struct Graph* createGraph(int V)
{
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->vertices = V;
 
    // Create an array of adjacency lists.  Size of array will be V
    graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));
 
     // Initialize each adjacency list as empty by making head as NULL
    int i;
    for (i = 0; i < V; ++i){
        
        graph->array[i].head = NULL;
    }
        
    return graph;
}
*/

// Adds an edge to an undirected graph
void addEdge(struct Graph* graph, int src, int dest)
{
    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    struct AdjListNode* newNode = newAdjListNode(dest);
    newNode->next = graph->array.at(src)->head;
    graph->array.at(src)->head = newNode;
    if (graph->array.at(src)->list_used ==  false){
      graph->array.at(src)->list_used = true;
    }
    if (graph->array.at(dest)->list_used ==  false){
      graph->array.at(dest)->list_used = true;
    }
 
    /*
    // Since graph is undirected, add an edge from dest to src also
    newNode = newAdjListNode(src);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
    */
}
/*
void addEdgetoList(struct AdjList* list, int dest)
{
    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    struct AdjListNode* newNode = newAdjListNode(dest);
    newNode->next = list.head;
    list.head = newNode;
 
    
    // Since graph is undirected, add an edge from dest to src also
    newNode = newAdjListNode(src);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
    
}
*/

struct Graph* init_graph_from_file(){
  int vertices, edges;
  int vertex, adj_vertex;
  string junk;
  size_t num_start, num_end;
  char colon_delimeter = ':', space_delimeter = ' ';
  //string graph_file = "test_graph_mini.txt";
  //string graph_file = "wiki_vote_mini.txt";
  string graph_file = "input_graphs//soc-Epinions1.txt";
  //string graph_file = "Wiki-Vote.txt";
  ifstream fin;
  fin.open(graph_file);
  if (fin.fail()){
    cout << "Could not open file" << endl;    
  }
  std::getline(fin, junk);
  std::getline(fin, junk); //Read first two lines and save in junk
  std::getline(fin, junk); //Read third line and extract number of vertices and edges
  num_start = junk.find (colon_delimeter, 0);
  num_end = junk.find (space_delimeter, num_start+2);
  vertices = std::stoi(junk.substr(num_start+2, (num_end-(num_start+2))));
  num_start = junk.find (colon_delimeter, num_end);
  edges = std::stoi(junk.substr(num_start+2));
  struct Graph* graph = new Graph();
  // Initialize each adjacency list as empty by making head as NULL
  //Initialize colors vector to unique identifier of its vertex id
  int i;
  for (i = 0; i < vertices; ++i){
    AdjList* new_list = new AdjList();
    new_list->node_id = i;
    new_list->list_used = false;
    new_list->head = NULL;
    graph->array.push_back(new_list);
    graph->colors.push_back(-1);
    graph->unique_color.push_back(true);
    graph->valid_vertex.push_back(true);
  }    
  graph->vertices = vertices;
  graph->edges = edges;
  std::getline(fin, junk); //Read another line in junk
  while( fin >> vertex >> adj_vertex){
    addEdge(graph, vertex, adj_vertex);    
  }
  return graph;
}

// A utility function to print the adjacenncy list representation of graph
void printGraph(struct Graph* graph)
{
    int v;
    for (v = 0; v < graph->vertices; ++v)
    {
        struct AdjListNode* pCrawl = graph->array[v]->head;
        if (pCrawl !=NULL){
            printf("\n Adjacency list of vertex %d\n head ", v);
            while (pCrawl)
            {
                printf("-> %d", pCrawl->dest);
                pCrawl = pCrawl->next;
            }
            printf("\n");
        }
        
    }
}

bool isconnected (int vertex_num, int neighbor_num, std::vector<AdjList*>& current_vec){
  struct AdjListNode* vertex_crawler = current_vec.at(vertex_num)->head;
  if (vertex_crawler !=NULL){
            while (vertex_crawler)
            {
              if (vertex_crawler->dest == neighbor_num){
                return true;
              } else vertex_crawler = vertex_crawler->next;
            }
  }
  return false;
}
void print_all_scc_vectors(std::vector<pair<std::vector<int>, int> >& all_scc_vector, struct Graph*& graph){
  cout << "The size of the scc_vector is " << all_scc_vector.size() << endl;
  int num_of_scc = 0;
  int max_size_scc = -1*(numeric_limits<int>::max()), max_size_scc_root = -1*(numeric_limits<int>::max());
  for (int i = 0; i < all_scc_vector.size(); i++){
    if (graph->array.at(all_scc_vector.at(i).first.front())->list_used == true){
      num_of_scc++;
      int curr_vec_size = all_scc_vector.at(i).second;
      if (curr_vec_size > max_size_scc){
        max_size_scc = curr_vec_size;
        max_size_scc_root = all_scc_vector.at(i).first.front();
      }
      //cout << "SCC for vertex " << all_scc_vector.at(i).front() << endl;
      /*for (int j = 0; j< curr_vec_size; j++){
        //cout << all_scc_vector.at(i).at(j) << " ";
        if (j == curr_vec_size-1){
        //  cout<<endl;
        }
      }*/
    }
  }
  cout << "The number of SCC in graph, maximum SCC size and maximum SCC root are " << num_of_scc << ", " << max_size_scc << " and " 
  << max_size_scc_root << " respectively." << endl;
}

int get_pivot(int last_pivot, struct Graph*& graph, bool& v_size_changed, size_t& v_size){
  for (int vertex_counter = last_pivot; vertex_counter >=0 ; vertex_counter--)
  {
    if (graph->array.at(vertex_counter)->list_used == true && graph->valid_vertex.at(vertex_counter)==true){
      if (v_size > 0) { 
        //cout << "Old V_Size is: " << v_size << endl;
        v_size--;
        //cout << "New V_Size is: " << v_size << endl;
      }
      if (v_size_changed == false ) 
      { 
        v_size_changed = true;
      }
      graph->valid_vertex.at(vertex_counter)=false;
      return vertex_counter;
    }    
  }
  return -1;
}
 
void scc_fwbw (struct Graph*& graph){
  std::vector<pair<std::vector<int>, int> > all_scc_vector;
  bool v_size_changed;
  size_t v_size = graph->valid_vertex.size();
  int vertex_to_start_choosing_pivot_from = v_size-1;
  int pivot;
  pivot = get_pivot(vertex_to_start_choosing_pivot_from, graph, v_size_changed, v_size);
  vertex_to_start_choosing_pivot_from = pivot;
  omp_lock_t writelock;
  omp_init_lock(&writelock);
  do
  {
    v_size_changed = false;
    std::vector<int> descendant_set;
    descendant_set.push_back(pivot);
    int new_desc_set_start = 0, new_desc_set_end = 1;
    bool desc_set_member_found = false;
    do
    {  
      int num_desc_set_members_found = 0;
      desc_set_member_found = false;
      #pragma omp parallel for schedule(static)
      for (int i = new_desc_set_start; i< new_desc_set_end; i++){
        struct AdjListNode* pCrawl = graph->array[descendant_set.at(i)]->head;
        if (pCrawl != NULL){
          while (pCrawl)
          {
            omp_set_lock(&writelock); 
              std::vector<int>::iterator it = find (descendant_set.begin(), descendant_set.end(), pCrawl->dest);
            omp_unset_lock(&writelock); 
            if (graph->valid_vertex.at(pCrawl->dest)==true && it == descendant_set.end()){
                omp_set_lock(&writelock);
                  descendant_set.push_back(pCrawl->dest);
                  num_desc_set_members_found++;
                  if (desc_set_member_found == false)
                  {
                    desc_set_member_found = true;
                  } 
                omp_unset_lock(&writelock);          
            }
            omp_set_lock(&writelock);
              pCrawl = pCrawl->next; 
            omp_unset_lock(&writelock);     
          }
        }
      }
      new_desc_set_start = new_desc_set_end;
      new_desc_set_end+=num_desc_set_members_found; 
    } while (desc_set_member_found == true);
    
    std::vector<int> predecessor_set;
    std::vector<int> new_scc_set;
    predecessor_set.push_back(pivot);
    new_scc_set.push_back(pivot);
    int new_pred_set_start = 0, new_pred_set_end = 1;
    bool pred_set_member_found = false;
    do
    {
      pred_set_member_found = false;
      int num_pred_set_members_found = 0;
      for (int i = new_pred_set_start; i < new_pred_set_end; ++i)
      {
        #pragma omp parallel for schedule(static)
        for (int vertex_counter = 0; vertex_counter < graph->array.size(); vertex_counter++)
        {
          if (graph->array.at(vertex_counter)->list_used == true && graph->valid_vertex.at(vertex_counter)==true){
            std::vector<int>::iterator it = find (predecessor_set.begin(), predecessor_set.end(), vertex_counter);
            if (it == predecessor_set.end()){      
              if (isconnected (vertex_counter, predecessor_set.at(i), graph->array)){
                std::vector<int>::iterator it = find (descendant_set.begin(), descendant_set.end(), vertex_counter);
                if (it != descendant_set.end()){
                  omp_set_lock(&writelock); 
                    new_scc_set.push_back(vertex_counter);
                    graph->valid_vertex.at(vertex_counter)=false;
                    v_size--;
                    predecessor_set.push_back(vertex_counter);
                    num_pred_set_members_found++;
                    if (pred_set_member_found == false)
                    {
                      pred_set_member_found = true;
                    }  
                  omp_unset_lock(&writelock); 
                }              
              }
            }
          } 
        } 
      }
      new_pred_set_start = new_pred_set_end;
      new_pred_set_end+=num_pred_set_members_found;
    } while (pred_set_member_found);
    all_scc_vector.push_back(std::make_pair(new_scc_set, new_scc_set.size()));
    pivot = get_pivot(vertex_to_start_choosing_pivot_from, graph, v_size_changed, v_size);
    vertex_to_start_choosing_pivot_from = pivot;
    //cout << "Pivot is " << pivot << ", has v_size changed? " << v_size_changed << endl;
    //cout << "V_size is " << v_size << endl;
  } while (v_size > 0 && pivot >= 0 && v_size_changed == true );
  print_all_scc_vectors(all_scc_vector, graph);
}
/*











  //std::vector<pair<std::vector<int>, int> > myVec
  //auto p1 = std::make_pair(n, a[1]);
  std::vector<AdjList*> current_vec = graph->array;
  size_t v_size = graph->valid_vertex.size();
  std::vector<pair<std::vector<int>, int> > all_scc_vector;
  bool v_size_changed;
  //cout << "Started parallel color" << endl;
  //while there are vertices left in the graph to compute scc
  time_t start, finish;
  time(&start);
  do{
    //cout << "V_size is " << v_size << endl;
    v_size_changed = false;    
    omp_lock_t writelock;
    omp_init_lock(&writelock);  
    #pragma omp parallel for schedule(static)
      for (int vertex_counter = 0; vertex_counter < graph->array.size(); vertex_counter++){        
        while (vertex_counter < graph->array.size() && graph->valid_vertex.at(vertex_counter)==false){
          //omp_set_lock(&writelock);
            //cout << "OMP Thread " << omp_get_thread_num() << " is working." << endl;
            vertex_counter++;
          //omp_unset_lock(&writelock);
        }
        if (vertex_counter < graph->array.size() && graph->valid_vertex.at(vertex_counter)==true){
          int u_value = current_vec.at(vertex_counter)->node_id;
          graph->colors[u_value] = u_value;
          graph->unique_color[u_value] = true;
        }      
        //vertex_counter++;
        //cout << "Ended parallel color by thread "<< omp_get_thread_num() <<". Vertex counter is " << vertex_counter << endl;
      }
    //time(&finish);
    //cout<< "Assigning color took " << difftime(finish, start) << " seconds." << endl;
    bool vertex_changed_color;
    //cout << "All coloring done" << endl;
    do{
      vertex_changed_color = false;
      //time_t start, finish;
      //time(&start);
      #pragma omp parallel for schedule(static)
        for (int vert_iter = 0; vert_iter < graph->array.size(); vert_iter++){
          while (vert_iter < graph->array.size() && graph->valid_vertex.at(vert_iter)==false){
            //omp_set_lock(&writelock);
              vert_iter++;
            //omp_unset_lock(&writelock);
          }
          if (vert_iter < graph->array.size() && graph->valid_vertex.at(vert_iter)==true){
            struct AdjListNode* vertex_crawler = current_vec.at(vert_iter)->head;
            if (vertex_crawler !=NULL){
                int u_value = current_vec.at(vert_iter)->node_id;
                while (vertex_crawler)
                {
                  int v_value = vertex_crawler->dest;
                  if (graph->colors[u_value]>graph->colors[v_value]){
                    graph->colors[v_value]=graph->colors[u_value];
                    graph->unique_color[v_value] = false;
                    if (vertex_changed_color == false){
                      vertex_changed_color = true;
                    }
                  }
                  vertex_crawler = vertex_crawler->next;
                }
            }
          }
          //cout << "The number of threads is " << omp_get_num_threads() << endl;
          //cout << "Color propagation by thread "<< omp_get_thread_num() <<". Vertex iterator is " << vertex_counter << endl;        
        }
      //time(&finish);
      //cout<< "Propagating color took " << difftime(finish, start) << " seconds." << endl;
    } while (vertex_changed_color == true);
    
    //time(&start);
    #pragma omp parallel for schedule(guided)
      for (int color_iter = 0; color_iter < graph->colors.size(); color_iter++){
        if (graph->valid_vertex.at(color_iter)==true && graph->unique_color.at(color_iter)==true && graph->array.at(color_iter)->list_used == true){
          std::vector<int> new_vc;
          std::vector<bool> new_vc_added_to_scv;
          new_vc.push_back(color_iter);
          new_vc_added_to_scv.push_back(true);
          for (int new_vc_counter = 0; new_vc_counter < graph->array.size(); new_vc_counter++){
            if ((graph->valid_vertex.at(new_vc_counter)==true)&&(graph->colors[new_vc_counter]==color_iter) && (graph->unique_color[new_vc_counter] == false)){
              new_vc.push_back(new_vc_counter);
              new_vc_added_to_scv.push_back(false);
            }
          }
          graph->valid_vertex.at(color_iter) = false;
          v_size--;
          if (v_size_changed == false ) { v_size_changed = true;}
          int scv_new_set_start = 0, scv_new_set_end = 1;
          int first_scv_entry = 0, last_scv_entry = 0;
          bool new_scv_member_found;        
          do{
              new_scv_member_found = false;
              int num_members_found = 0;
              int new_vc_iter = last_scv_entry+1;
              for (new_vc_iter; new_vc_iter<new_vc.size(); new_vc_iter++){
                int scv_counter = scv_new_set_start;
                while (scv_counter < scv_new_set_end){
                  if (isconnected (new_vc.at(new_vc_iter), new_vc.at(scv_counter), current_vec)){
                    last_scv_entry++;
                    std::swap(new_vc.at(last_scv_entry),new_vc.at(new_vc_iter));
                    graph->valid_vertex.at(new_vc.at(last_scv_entry)) = false;
                    new_vc_added_to_scv.at(last_scv_entry)=true;
                    v_size--;
                    num_members_found++;
                    if (new_scv_member_found == false){
                      new_scv_member_found = true;
                    }
                    break;
                  }
                  scv_counter++;
                }
              }
              scv_new_set_start = scv_new_set_end;
              scv_new_set_end += num_members_found;
          } while (new_scv_member_found);
          omp_set_lock(&writelock);
          all_scc_vector.push_back(std::make_pair(new_vc, scv_new_set_end));
          omp_unset_lock(&writelock);
        }
      }
    //time(&finish);
    //cout<< "Populating SCC took " << difftime(finish, start) << " seconds." << endl;
    //cout << " .The number of vertices is " << v_size << endl;
  } while (v_size>0 && (v_size_changed == true));
  time(&finish);
  cout<< "Parallel coloring took " << difftime(finish, start) << " seconds." << endl;  
  print_all_scc_vectors(all_scc_vector, graph);
}

*/