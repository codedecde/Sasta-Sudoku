#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include "sudoku.h"
// #define SIZE 9
#define BOARD_SIZE (SIZE*SIZE + 1)
#define NO_SOLN 2
#define PARTIAL_SOLN 1
#define SOLVED 0

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
			if (this->n_allowed = idx_list) { // 1-indexed means that their equality is what we care about.
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
			int idx_list = 1;
			for(; this->list[idx_list]; ++idx_list) {
				this->base_array[ this->list[idx_list] ] = 0; // Say that the value is not allowed.
					// This is the equivalent of free-ing the base_array in old cell_t
			}
			return 1; // On success?
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