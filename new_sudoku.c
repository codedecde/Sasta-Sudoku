#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include "sudoku.h"
// #define SIZE 9
#define BOARD_SIZE (SIZE*SIZE + 1)
#define NO_SOLN 2
#define PARTIAL_SOLN 1
#define SOLVED 0
#define MAX_INT 2*SIZE

struct cell_t; // r, c, allowed, n_allowed
typedef struct cell_t cell_t;
	// Functions of cell_t
	void cell_init(cell_t* this);
	int is_allowed(cell_t* this, int val);
	void add_allowed(cell_t*, int);
	void remove_allowed(cell_t* this, int val);
	int set_value(cell_t*, int);
	int cell_destructor(cell_t*);

struct pool_t; // Thread's access this... array of board's, FIFO/LIFO/PQ behaviour to be defined.
	struct pool_node_t; // private struct for pool_t really.
typedef struct pool_t pool_t;
	typedef struct pool_node_t pool_node_t;
	
	// Functions of : pool_t
	int is_empty(pool_t* this);
	void free_pool(pool_t* this);
	void pushf(pool_t*, cell_t*);
	void pushb(pool_t*, cell_t*);
	cell_t* popf(pool_t*);
	cell_t* popb(pool_t*);

struct cell_t {
	int row;
	int col;
	int* list;
	int* base_array;
	int n_allowed;
	int value;
};
	// RELATED FUNCTIONS;
	void cell_init(cell_t* this) {
		this->list = (int*) malloc( (SIZE+1)*sizeof(int) ); // idx0 corresponds to zilch.
		this->base_array = (int*) malloc( (SIZE+1)*sizeof(int) ); // idx0 is useless... value 0 means not allowed.
		this->n_allowed = 0;
		this->value = 0;
	}

	int is_allowed(cell_t* this, int val) {
		return (this->base_array[val]?1:0);
	}
	void add_allowed(cell_t* this, int val) {
		if(this->base_array[val]) {
			return; // Do nothing;	
		} else {
			this->n_allowed++; // This must happen first because our list and base_array are 1-indexed
			this->list[this->n_allowed] = val;
			this->base_array[val] = this->n_allowed;
		}
	}

	void remove_allowed(cell_t* this, int val) {
		int idx_list = this->base_array[val];
		if (idx_list) { // it is present in the array.
			// Case 1: It is the tail of the list.
			if (this->n_allowed == idx_list) { // 1-indexed means that their equality is what we care about.
				this->n_allowed--;
				this->base_array[val] = 0; // The value is not allowed.
			} else {
				// Case 2: It is not the tail of the list
					// Subcase 1: It is the only element...in which case it'll also be the tail of the list.
					// Subcase 2: There is some other element that is the tail of the array.	
				// We only need to worry about subcase 2.
				// Swap the items at idx_list and this->n_allowed;
				int temp_val;
				this->base_array[val] = 0;
				this->list[idx_list] = this->list[this->n_allowed]; // New value.
				this->base_array[ this->list[this->n_allowed] ] = idx_list;
				this->n_allowed--;
			}
		}
	}

	int set_value(cell_t* this, int val) {
		if(this->value == 0) {
			this->value = val;
			this->n_allowed = 0;
			int idx_base_array = 1;
			for(; idx_base_array <= SIZE; ++idx_base_array) {
				this->base_array[ idx_base_array ] = 0; // Say that the value is not allowed.
					// This is the equivalent of free-ing the base_array in old cell_t
			}
			return 1; // On success?
		}
		else{
			this->value = val;
		}
		return -1;
	}
	int cell_destructor(cell_t* this) {
		free(this->list);
		free(this->base_array);
	}
	// WHAT WOULD HAVE BEEN A GORGEOUS cell_t struct :'(


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
	// pool_t functions begin here.
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
	// END OF WHAT WOULD'VE BEEN AN AMAZING (dead)pool.
/*
	Convert a cell_t* board to an int** board
	return: int** 
*/
	int** cell2board(cell_t* board){
		int** new_board = malloc(SIZE * sizeof(int*));
		int i;
		for(i = 0;i < SIZE;i++){
			new_board[i] = malloc(SIZE*sizeof(int));
			int j;
			for(j = 0;j < SIZE;j++){
				int idx = SIZE*i + j + 1;
				new_board[i][j] = board[idx].value;
			}
		}
		return new_board;
	}
/*
  Function call made from main reaches here
*/
int** solveSudoku(int** originalGrid){
  cell_t board[SIZE*SIZE + 1];
  int r,c;
  for(r = 0;r < SIZE;r++){
	for(c = 0; c < SIZE; c++){
	  int idx = r*SIZE + c + 1;
	  cell_init(&board[idx]);
	  if(originalGrid[r][c] != 0){
	  	set_value( &board[idx], originalGrid[r][c] );
	  }
	  else{
		int bitmask = (1<<(SIZE+1) ) - 2;
		int jdx = 0;
		int r_prime = r - (r % MINIGRIDSIZE);
		int c_prime = c - (c % MINIGRIDSIZE);
		for(;jdx<SIZE;jdx++){
			bitmask &= ~(1<<originalGrid[r][jdx]);
			bitmask &= ~(1<<originalGrid[jdx][c]);
			bitmask &= ~(1<<originalGrid[r_prime + (jdx/MINIGRIDSIZE)][c_prime + (jdx % MINIGRIDSIZE)]);
		}
		jdx = 1;
		for(;jdx<=SIZE;jdx++){
			if(bitmask & 1<<jdx){
				add_allowed(&board[idx],jdx);
			}
		}
	  }
	}
  }

  int error = dfs(board);
  return cell2board(board);

}

/*
  Heuristic solve: Applies elimination (TODO: OTHER STUFF)
  returns 0 for unsolvable, 
  return 1 for partially_solved,
  return 2 for completely_solved
*/
int heuristic_solve(cell_t* board){
	int flag = 1;
  	int err_code = -1;
  	while(flag){
	  	flag = 0;
	  	err_code = SOLVED;
	  	int r,c;
	  	for(r=0;r<SIZE;r++){
		  	for(c=0;c<SIZE;c++){ // c++ ! Hah! Get it? :D
			  	int idx = r*SIZE + c + 1;
			  	if(board[idx].n_allowed == 0){
				  	if(board[idx].value == 0){ // value is still unfilled.
					  	err_code = NO_SOLN;
					  	return err_code;
				  	}
			  	}
			  	else if(board[idx].n_allowed == 1){
				  	flag = 1;
				  	int val = board[idx].list[1]; // 1 indexed array
				  	set_value(&board[idx],val);
				  	int r_prime = r - (r % MINIGRIDSIZE);
				  	int c_prime = c - (c % MINIGRIDSIZE);
				  	int jdx = 0;
				  	for(;jdx<SIZE;++jdx){
					  	int jdx_r = r*SIZE + jdx + 1;
					  	int jdx_c = jdx*SIZE + c + 1;
					  	int jdx_b = ((r_prime + (jdx / MINIGRIDSIZE))*SIZE) + (c_prime + (jdx%MINIGRIDSIZE) ) + 1;
					  	remove_allowed(&board[jdx_r],val);
					  	remove_allowed(&board[jdx_c],val); 
					  	remove_allowed(&board[jdx_b],val);
				  	}
			  	}
			  	else if( board[idx].n_allowed > 1){ // Made it a check for sanity purposes.
				  	err_code = PARTIAL_SOLN;
			  	}
		  	}
	  	}
  	}
	return err_code;
}

int dfs(cell_t* board){
	/*
	Currently iterate through the entire board
	Can we do better ?? 
	*/
	int min_val = MAX_INT;
	int min_idx = -1;
	int i = 1;
	for(;i<BOARD_SIZE;i++){
		if(board[i].n_allowed != 0){
			if(min_val > board[i].n_allowed){
				min_val = board[i].n_allowed;
				min_idx = i;
			}
		}
		else if(board[i].n_allowed == 0){
			if(board[i].value == 0){
				return NO_SOLN;
			}
		}
	}
	if(min_idx == -1)
		return SOLVED;
	int num_allowed = board[min_idx].n_allowed;
	int r = (min_idx - 1) / SIZE;
	int c = (min_idx - 1) % SIZE;
	int r_prime = r - (r % MINIGRIDSIZE);
	int c_prime = c - (c % MINIGRIDSIZE);
	int base_array[SIZE+1];
	base_array[0] = 0;
	for(i = 1; i <= num_allowed;i++){
		base_array[i] = board[min_idx].base_array[i];
	}
	for(i = 1; i <= num_allowed;i++){
		int val = board[min_idx].list[i];
		set_value(&board[min_idx],val);
		int jdx = 0;
		for(;jdx < SIZE;jdx++){
			int jdx_r = r*SIZE + jdx + 1;
		  	int jdx_c = jdx*SIZE + c + 1;
		  	int jdx_b = ((r_prime + (jdx / MINIGRIDSIZE))*SIZE) + (c_prime + (jdx%MINIGRIDSIZE) ) + 1;
		  	remove_allowed(&board[jdx_r],val);
		  	remove_allowed(&board[jdx_c],val); 
		  	remove_allowed(&board[jdx_b],val);
		}
		int error_code = dfs(board);
		if(error_code == SOLVED){
			return SOLVED;
		}
		for(;jdx < SIZE;jdx++){
			int jdx_r = r*SIZE + jdx + 1;
		  	int jdx_c = jdx*SIZE + c + 1;
		  	int jdx_b = ((r_prime + (jdx / MINIGRIDSIZE))*SIZE) + (c_prime + (jdx%MINIGRIDSIZE) ) + 1;
		  	add_allowed(&board[jdx_r],val);
		  	add_allowed(&board[jdx_c],val); 
		  	add_allowed(&board[jdx_b],val);
		}
	}
	board[min_idx].value = 0;
	board[min_idx].n_allowed = num_allowed;
	for(i = 1; i <= num_allowed;i++){
		board[min_idx].base_array[i] = base_array[i];
	}
	return NO_SOLN;
}

