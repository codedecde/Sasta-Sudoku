#include <stdio.h>  // for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include "sudoku.h"
// #define SIZE 9
#define BOARD_SIZE (SIZE*SIZE + 1)

struct cell_t; // r, c, allowed, n_allowed
struct pool_t; // Thread's access this... array of board's, FIFO/LIFO/PQ behaviour to be defined.
  struct pool_node_t; // private struct for pool_t really.
struct heap_t;

typedef struct cell_t cell_t;
typedef struct pool_t pool_t;
typedef struct pool_node_t pool_node_t;
typedef struct heap_t heap_t;
// Functions of : cell_t
  #define NOT_ALLOWED 37
  int is_allowed(cell_t*, int);
  void set_allowed(cell_t*, int);
  void unset_allowed(cell_t*, int);

// Functions of : pool_t
  int is_empty(pool_t* this);
  void free_pool(pool_t* this);
  void pushf(pool_t*, cell_t*);
  void pushb(pool_t*, cell_t*);
  cell_t* popf(pool_t*);
  cell_t* popb(pool_t*);

// Function of : heap_t
  void heap_init(heap_t* heap,cell_t* init_array);
  void build_heap(heap_t* heap);
  void heapify(heap_t* heap,int i);
  cell_t* heap_top(heap_t* heap);
  void heap_destructor(heap_t*);
// Definitions Now.
struct cell_t {
  int row; //Row
  int col; //Col
  int allowed; //Bitmask of values allowed
  int n_allowed; //Number of values allowed
  int value; // value of cell, if filled
};  // cell_t ke functions.
  int is_allowed(cell_t* this, int val){
    //0 if not allowed
    return this->allowed & 1<<val;
  }

  void set_allowed(cell_t* this, int val){
    /*
    Check if already set. 
    If not update allowed and n_allowed
    */
    if( (this->allowed & 1<<val) ){
        return;
    }
    else{
        this->allowed |= 1<<val;
        this->n_allowed ++;
    }
  }

  void unset_allowed(cell_t* this, int val){
    /*
    Check if already set.
    If set, then unset and update n_allowed . 
    Handle special case of n_allowed = 0
    Else return silently 
    */
    if ((this->allowed & 1<<val)){
        this->allowed &= ~(1<<val);
        this->n_allowed--;
        if(!(this->n_allowed)){
            (this->n_allowed) = NOT_ALLOWED;
        }
    }
  }
  // END OF WHAT WOULD HAVE BEEN A BEAUTIFUL cell_t definition

struct pool_node_t {
  cell_t* board_config; // malloc-ing
  pool_node_t* next;
  pool_node_t* prev;
};
struct pool_t { // A doubly linked list. Allows for FIFO/LIFO.
  pool_node_t* hd;
  pool_node_t* tl;
  int sz;
};
  int is_empty(pool_t* this) {
    return (this->sz > 0 ? 1 : 0);
  }

  void free_pool(pool_t* this) {
    while (this->sz) {
      cell_t* gonna_die_soon = popf(this); // searched
      free(gonna_die_soon); // dieded.
    }
  }
  /*
    Invocation time state : [p,q,r,s..]
    pushf(x)
    Return time state: [x,p,q,r,s...]
   */
  void pushf(pool_t* this, cell_t* new_board) { // ASSERT: cell_t* is a fixed size board etc.
    struct pool_node_t* new_node = (pool_node_t*) malloc(sizeof(pool_node_t));
    new_node->board_config = (cell_t*)malloc((BOARD_SIZE)*sizeof(cell_t));
    memcpy((void*)new_node->board_config, (void*)new_board, (BOARD_SIZE)*sizeof(cell_t));
    
    new_node->next = this->hd;
    if(this->hd) { // this->hd is not NULL
      this->hd->prev = new_node;
    }
    this->hd = new_node;
    this->sz++;
    if(this->sz == 1) {
      this->tl = this->hd;
    }
  }
  /*
    Invocation time state : [...,p,q,r,s]
    pushb(x)
    Return time state: [...,p,q,r,s,x]
  */
  void pushb(pool_t* this, cell_t* new_board) {
    struct pool_node_t* new_node = (pool_node_t*) malloc(sizeof(pool_node_t));
    new_node->board_config = (cell_t*)malloc((BOARD_SIZE)*sizeof(cell_t));
    memcpy((void*)new_node->board_config, (void*)new_board, (BOARD_SIZE)*sizeof(cell_t));

    new_node->next = this->tl;
    if(this->tl) { // this->tl is not NULL
      this->tl->next = new_node;
    }
    this->tl = new_node;
    this->sz++;
    if(this->sz == 1) {
      this->hd = this->tl;
    }
  }
  /*
    Invocation time state : [p,q,r,s,...]
    return [p]
    
  */
  cell_t* popf(pool_t* this) {
    if (!this->hd) { // Empty pool.
      return NULL;
    } else if (this->hd == this->tl) { // only one element in pool
      this->tl = NULL;
    }
    pool_node_t* ret_node = this->hd;
    this->hd = ret_node->next;
    if(this->hd) {
      this->hd->prev = NULL;
    }
    this->sz--;
    cell_t* retval = ret_node->board_config;
    free(ret_node); // Doesn't free retval...or what its pointing at.
    return retval;
  }
  /*
    Invocation time state : [...,p,q,r,s]
    return s
  */  
  cell_t* popb(pool_t* this) {
    if(!this->tl) {
      return NULL;
    } else if ( this->hd == this->tl ) {
      this->hd = NULL;
    }
    pool_node_t* ret_node = this->tl;
    this->tl = ret_node->prev;
    if(this->tl) {
      this->tl->next = NULL;
    }
    this->sz--;
    cell_t* retval = ret_node->board_config;
    free(ret_node);
    return retval;
  }
  // END OF WHAT WOULD HAVE BEEN A BEAUTIFUL pool_t* 


struct heap_t{
  cell_t** base_array;
}; // START OF heapt_t functions.
  void heap_init(heap_t* heap,cell_t* init_array){
    heap->base_array = malloc(BOARD_SIZE * sizeof(cell_t*));
    heap->base_array[0] = NULL;
    int i = 1;
    for(;i<BOARD_SIZE;i++){
      heap->base_array[i] = &init_array[i];
    }
  }
  void build_heap(heap_t* heap){
    int n = BOARD_SIZE / 2;
    while(n > 0){
      heapify(heap,n); //Heapify from bottom up. Happens in O(n) That is the heap magic !!
      n--;
    }
  }
  void heapify(heap_t* heap,int i){
    /*
      min stores the index in the base_array, which points to a cell_t of lowest n_allowed
      check parent has idx == min 
      else swap and recurse 
    */
    int min = i;    
    int lchild = 2*i; 
    int rchild = 2*i+1; 
    if(lchild < BOARD_SIZE){ 
      if( ( (*( heap->base_array[lchild]) ).n_allowed) < ( (*( heap->base_array[min]) ).n_allowed ) ){
        min = lchild; 
      }
    }
    if(rchild < BOARD_SIZE){
      if( ( (*( heap->base_array[rchild]) ).n_allowed)  < ( (*( heap->base_array[min]) ).n_allowed) ){
        min = rchild;
      }
    }
    if(min != i){
      cell_t* temp = heap->base_array[i];
      heap->base_array[i] = heap->base_array[min];
      heap->base_array[min] = temp;
      heapify(heap,min);
    }
  }
  cell_t* heap_top(heap_t* heap){
    if(heap->base_array)
      return heap->base_array[1];
    else
      return NULL;
  }
  void heap_destructor(heap_t* heap){
    //Because ~ is too mainstream for creators of c
    free(heap->base_array);
  }
// END OF WHAT WOULD'VE BEEN A GORGEOUS heap_t