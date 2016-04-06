#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include <assert.h>
#include <limits.h>
#include <omp.h>	// OpenMP
#include "sudoku.h"

#define BOARD_SIZE (SIZE*SIZE + 1)

#define NO_SOLN 2
#define PARTIAL_SOLN 1
#define SOLVED 0

#define INITIAL_WORKPOOL_SZ 30
#define MAX_WORKPOOL_SZ 40
#define INITIAL_STORE_SZ 100

struct cell_t {
	int row;
	int col;
	int *list;
	int *base_array;
	int n_allowed;
	int value;
};
typedef struct cell_t cell_t;
	int is_allowed(cell_t* this, int val);
	int add_allowed(cell_t*, int); // returns -1 on behaviour dependent on this having a non-zero value
	int remove_allowed(cell_t* this, int val); // returns -1 on behaviour dependent on this having a non-zero value
	int set_value(cell_t*, int);
// cell_t functions.
	void cell_init(cell_t* t) { // TODO: Remove this once safe.
		printf("cell_init is no longer a valid function\n"); 
		assert(0);
	}
	void cell_destructor(cell_t* this) {
		free(this->list);
		free(this->base_array);
		assert("cell_destructor is ambiguously supported\n");
	}
	void print_cell(cell_t* this) {
		printf("cell(%d,%d)=%d\n", this->row, this->col, this->value);
		printf("list: ");
		int i;
		for(i=0; i<(SIZE+1); i++) {
			printf("%d ", this->list[i]);
		}
		printf("\nbase_array: ");
		for(i=0; i<(SIZE+1); i++) {
			printf("%d ", this->base_array[i]);
		}
	}
	int is_allowed(cell_t* this, int val) {
		return (((val>0) && (this->base_array[val] > 0)) ? 1 : 0 );
	}
	int add_allowed(cell_t* this, int val) {
		if (this->value == 0) {
			if (this->base_array[val] > 0) { // val is already allowed
				return 0;
			} else {
				this->n_allowed++;
				this->list[this->n_allowed] = val;
				this->base_array[val] = this->n_allowed;
				return 1;
			}		
		} else { // value is nonzero, nothing is changed. base_array and nallowed must remain as earlier.
			return -1;
		}
	}
	int remove_allowed(cell_t* this, int val) {
		int idx_list = this->base_array[val];
		if ( this->value != 0 ) {// This is the case where nothing is allowed anyway. It should be covered by other case.
			return -1; // Special behaviour, no?
		} 
		if ( idx_list==0 ) { // This is the case where the given val is not actually allowed.)
			return 0;
		}
		if (this->n_allowed == idx_list) { // 1-indexed means that their equality is what we care about.
			this->n_allowed--;
			assert(this->n_allowed >= 0);
			this->base_array[val] = 0; // The value is not allowed.
			return 1;
			// We make no guarantees about the this->list values at indices > n_allowed.
		} else {
			// Case 2: It is not the tail of the list
				// Subcase 1: It is the only element...in which case it'll also be the tail of the list.
				// Subcase 2: There is some other element that is the tail of the array.	
			// We only need to worry about subcase 2.
			// Swap the items at idx_list and this->n_allowed;
			// Note: THIS REORDERS THE FRIGGIN LIST. however, if my value has been set, it shouldn't matter.
			this->base_array[val] = 0;
			this->list[idx_list] = this->list[this->n_allowed]; // New value.
			this->base_array[ this->list[this->n_allowed] ] = idx_list;
			this->n_allowed--;
			assert(this->n_allowed >= 0 );
			return 1;
		}
	}
	int set_value(cell_t* this, int val) {
		if(this->value == 0) {
			this->value = val;
			this->n_allowed = 0;
			
			memset(this->base_array,0,(SIZE + 1)*sizeof(int));

			return 1; // On successfully unallowing things.
		} else {
			this->value = val;
			assert(this->n_allowed == 0);
			return 0; // 
		}
		return -1;
	}

	int set_all_allowed(cell_t* this) {
		// Set all possible values to "allowed"...
		if(this->value ==0) {
			int val;
			for(val=1; val<=SIZE; ++val) {
				this->base_array[val] = val;
				this->list[val] = val;
			}
			this->n_allowed = SIZE;
		} else {
			return -1; // Doesn't make sense to set everything to allowed.
		}
	}
	// END OF cell_t functions... except cell_init

	// Functions relating to a board.
	/*
		print_board function launches a nuke. Obviously.
	 */
	void print_board(cell_t* board) {
		int r,c;
		printf("\n\tprinting board\n");
		for(c=-1; c<SIZE; c++) {
			if(c==-1) printf("   ");
			else printf("%02d ", c);
		} 
		printf("\n");
		for(r=0; r<SIZE; r++) {
			printf("%02d ",r);
			for(c=0; c<SIZE; c++) {
				printf("%02d ", board[r*SIZE + c + 1].value);
			}
			printf("\n");
		}
		printf("\n\n");
	}
	int set_board_value_rc(cell_t* board, int r, int c, int val, int* altered_indices, int* ptr_naltered) {
		// ASSERT: altered_indices must be set correctly. This function is rather unsafe :)
		int idx, rp, cp, i, ir, ic, ib;
		rp = r - (r%MINIGRIDSIZE);
		cp = c - (c%MINIGRIDSIZE);
		idx = r*SIZE + c + 1;
		memset(altered_indices, 0, BOARD_SIZE*sizeof(int));
		if( altered_indices!=NULL && ptr_naltered!=NULL ) {
			*ptr_naltered = -1;
			set_value(&board[idx], val);
			for(i=0; i<SIZE; ++i) {
				ir = rp*SIZE + i + 1;
				ic = i*SIZE  + cp + 1;
				ib = (rp - (i/MINIGRIDSIZE))*SIZE + (cp - (i%MINIGRIDSIZE)) + 1;
				if( ir!=idx && remove_allowed(&board[ir], val)==1 ) {
					(*ptr_naltered)++;
					altered_indices[*ptr_naltered] = ir;
				}
				if( ic!=idx && remove_allowed(&board[ic], val)==1 ) {
					(*ptr_naltered)++;
					altered_indices[*ptr_naltered] = ic;
				}
				if( ib!=idx && remove_allowed(&board[ib], val)==1 ) {
					(*ptr_naltered)++;
					altered_indices[*ptr_naltered] = ib;
				}
			}
		} else { // ASSERT: altered_indices and ptr_naltered are NULL.
			set_value(&board[idx], val);
			for(i=0; i<SIZE; ++i) {
				ir = rp*SIZE + i + 1;
				ic = i*SIZE  + cp + 1;
				ib = (rp - (i/MINIGRIDSIZE))*SIZE + (cp - (i%MINIGRIDSIZE)) + 1;
				remove_allowed(&board[ir], val);
				remove_allowed(&board[ic], val);
				remove_allowed(&board[ib], val);
			}
		}
		return 1; // This function won't really fail, so its ok.
	}
	int set_board_value_idx(cell_t* board, int idx, int val, int* altered_indices, int* ptr_naltered) {
		return set_board_value_rc(board, ((idx-1)/SIZE), ((idx-1)%SIZE), val, altered_indices, ptr_naltered);
	}
	cell_t* init_board(cell_t* board) {
		// initialize the board... allocate it, and its list and base_array
		// Note: boards are 1 indexed. Hence, board[0] = NULL
		board = malloc(BOARD_SIZE*sizeof(cell_t));
		// Initialize the cells now.
		int i;
		// board[0] = NULL; ... How do I do this?
		for(i = 1; i<BOARD_SIZE; ++i) {
			board[i].row = (i-1)/SIZE;
			board[i].col = (i-1)%SIZE;
			board[i].value = 0;
			board[i].n_allowed = 0;
			board[i].list = malloc((SIZE+1)*sizeof(int));
			memset(board[i].list, 0, (SIZE+1)*sizeof(int));
			board[i].base_array = malloc((SIZE+1)*sizeof(int));	
			memset(board[i].base_array, 0, (SIZE+1)*sizeof(int));
		}
		return board;
	}

	int destroy_board(cell_t* board) { // Goes in deep and destroys the board entirely.
		int i;
		for(i=1; i<BOARD_SIZE; ++i) {
			free(board[i].list);
			free(board[i].base_array);
			free(&board[i]);
		}
		free(board);
		return 1;
	}

	int copy_board(cell_t* dst, cell_t* src) {
		// ASSERT : dst and src have both been allocated already.
		int i;
		for(i = 1; i<BOARD_SIZE; ++i) {
			// Deep copy src into dst.
			dst[i].row = src[i].row;
			dst[i].col = src[i].col;
			dst[i].value = src[i].value;
			dst[i].n_allowed = src[i].n_allowed;
			// ASSERT: dst[i].list and dst[i].base_array have been allocated.
			memcpy(dst[i].list, src[i].list, (SIZE+1)*sizeof(int));
			memcpy(dst[i].base_array, src[i].base_array, (SIZE+1)*sizeof(int));
		}
	}
	/*
		For outputting purposes.
	 */
	int** board2mat(cell_t* board) {
		int** mat = malloc(SIZE*sizeof(int*));
		int i,j;
		for(i = 0; i<SIZE; ++i) {
			mat[i] = malloc(SIZE*sizeof(int));
			for(j=0; j<SIZE; j++) {
				mat[i][j] = board[i*SIZE + j + 1].value;
			}
		}
		return mat;
	}

// pool struct
struct pool_t;
typedef struct pool_t pool_t;
struct pool_t {
	// IMPORTANT: arr_boards is 0 indexed
	cell_t** arr_boards; // It is an array of pointers, hence moving them around should be easier.
	int sz_alloc; 		 // Total size allocated
	int sz;				 // Present size (occupied)
	// INVAR: arr_boards[0...(sz-1)] are all initialized, no guarantees on the others.
};
// Creating a distinction for safer coding.
typedef pool_t workpool_t;
typedef pool_t board_store_t;
workpool_t* workpool; // Global
board_store_t* store; // Global

void destory_pool(workpool_t* pool) {
	int i;
	for(i=0; i<pool->sz; ++i) {
		destroy_board(pool->arr_boards[i]);
	}
	free(pool);
}

workpool_t* init_workpool(workpool_t* pool, int maxsz) {
	pool = malloc(sizeof(workpool_t));
	pool->arr_boards = malloc(maxsz*sizeof(cell_t*));
	pool->sz = 0;
	pool->sz_alloc = maxsz;
	return pool;
}
void destory_workpool(workpool_t* pool) {
	destory_pool((pool_t*)pool);
}
board_store_t* init_store(board_store_t* pool, int initsz) {
	pool = malloc(sizeof(board_store_t));
	pool->sz_alloc = initsz;
	pool->arr_boards = malloc(initsz*sizeof(cell_t*));
	int i;
	for(i=0; i<initsz; ++i) {
		pool->arr_boards[i] = init_board(pool->arr_boards[i]); // NOTE: I'm not very sure about C syntax here.
	}
	return pool;
}

void destory_store(board_store_t* pool) {
	destory_pool((pool_t*)pool);
}

int not_full(pool_t* pool) { // Applies to both, board_store_t and workpool_t
	return ((pool->sz < pool->sz_alloc)?1:0);
}

cell_t* pop_board(workpool_t* pool) {
	if(pool->sz>0) {
		pool->sz--;
		return pool->arr_boards[pool->sz];
	} else {
		return NULL;
	}
}

cell_t* pop_mem(board_store_t* pool) {
	if( not_full(pool) ) {
		printf("Store is empty\n");
		return pop_board(pool);
	} else {
		// Make a new board...allocate it, and return it.
		cell_t* newboard;
		newboard = init_board(newboard);
		return newboard;
	}
}
cell_t* push_board(workpool_t* pool, cell_t* board) {
	// IMPORTANT: The difference between push_board and push_mem is the copy that push_board performs, but push_mem doesnt
	// Case 1: Need to grow... 
	// Case 2: Don't need to grow... 
	// Return ptr to the board that was pushed.
	cell_t * retval;
	if(pool->sz >= pool->sz_alloc) { // CASE 1 ... time to grow.
		cell_t** oldarr = pool->arr_boards;
		pool->arr_boards = malloc(2*pool->sz_alloc*sizeof(cell_t*));
		pool->sz_alloc = 2*pool->sz_alloc;
		// IMPORTANT: I don't do a memset here because I make no guarantees on pool->arr_boards[pool->sz...] anyway
		memcpy(pool->arr_boards, oldarr, pool->sz_alloc*sizeof(cell_t*));
		free(oldarr); // Cleaning up after yourself.
	}
	retval = pop_mem(store);
	pool->arr_boards[pool->sz] = retval;
	copy_board(retval, board); // Go in deep and copy stuff. O(n) operation man.
	pool->sz++;
	return retval;
}
int push_mem(board_store_t* pool, cell_t* board) {
	// IMPORTANT: The difference between push_board and push_mem is the copy that push_board performs, but push_mem doesnt
	// Case1: Had to grow pool... return 1
	// Case2: Didn't have to grow... return 0
	int retval = 0;
	if( !not_full(pool) ) { // is full :P
		cell_t** oldarr = pool->arr_boards;
		pool->arr_boards = malloc(2*pool->sz_alloc*sizeof(cell_t*));
		// IMPORTANT: I don't allocate the new boards... I expect them to be malloced elsewhere.
		// 	Rationale:
		// 		the board must've come from somewhere.
		// 		Initially, all boards are born from the store.
		// 		That would mean that somebody is mallocing boards, Me allocing boards would hence be pointless.
		memcpy(pool->arr_boards, oldarr, pool->sz_alloc*sizeof(cell_t*));
		free(oldarr);
		retval = 1;
	}
	pool->arr_boards[pool->sz] = board;
	pool->sz++;
	return retval;
}
int dispose_board(cell_t* board, board_store_t* pool) {
	if(pool==NULL) {
		return destroy_board(board);
	} else {
		return push_mem(pool, board);
	}
}

// END OF BACK-END FUNCTIONS.


// What follows is the badass code.
// TODO: Everything from here on out.
// Functions used for initialization.
int copy_mat2board(int** mat, cell_t* board); // ASSERT: board and mat have been initialized.
int restore_invariants(cell_t* board);	// Invokes the wrath of god by calling meteor showers on all involved parties. Obviously.

int dfs(cell_t* board); // Take board ka dfs solution to completion
int heuristic_solve(cell_t* board); // Take board and solve it heuristically as much as possible.
int branch(cell_t* board, workpool_t* pool, int max_branches); // Find the children of board. Push STRICTLY LESS(EQ) than max_branches

int copy_mat2board(int** mat, cell_t* board) {
	// ASSERT: board has already been allocated.
	int r,c,idx;
	for(r=0; r<SIZE; ++r) {
		for(c=0; c<SIZE; ++c) {
			idx = r*SIZE + c + 1;
			board[idx].value = mat[r][c];
			board[idx].row = r;
			board[idx].col = c;
		}
	}
	return 1; // always succeed, really.
}
int restore_invariants(cell_t* board) {
	int i, r,c, rp,cp;
	int j, jr, jc, jb;
	for(i=1; i<BOARD_SIZE; ++i) {
		set_all_allowed(&board[i]);
	}

	for(i=1; i<BOARD_SIZE; ++i) {
		if(board[i].value == 0) continue; // This imposes no constraints on anybody.
		r = i/SIZE;
		c = i%SIZE;
		rp = r - (r%MINIGRIDSIZE);
		cp = c - (c%MINIGRIDSIZE);
		for(j=0; j<SIZE;++j) {
			jr = r*SIZE + j + 1;
			jc = j*SIZE + r + 1;
			jb = (rp + (j/MINIGRIDSIZE))*SIZE + (cp +(j%MINIGRIDSIZE) ) + 1;
			remove_allowed(&board[jr], board[i].value);
			remove_allowed(&board[jc], board[i].value);
			remove_allowed(&board[jb], board[i].value);
		}
	}
	return 1;
}
/*
	Behaves as my main. Pretty much.
 */
int** solveSudoku(int** originalGrid) {
	// Step1: Convert int** to orig_board.
	// Step2: initialize workpool, board_store_t.
	// Step 2.5: Attempt to solve the orig_board using heuristics. If SOLVED or NO_SOLN, return
	// ELSE :
	// <Parallelization etc follows>
	// branch() and put children of orig_board onto workpool. 
	
	// Workpool execution:
	// pop from workpool. 
	//  put it through heuristic solve. 
	//  If
	//  	SOLVED: stop execution and return
	//  	NO_SOLN: free(onto store) the board and repeat process.
	//  	PARTIAL_SOLN: try to branch and push onto workpool.
	//  					If fail, then try dfs...
	//  						If SOLVED: stop execution and return
	//  						If NO_SOLN: free(onto store) and repeat process
	//  					If succeed, then free(onto store) and repeat process.
	
	cell_t* solved_board;		// global/shared solution.
	int solved;				// global/shared boolean

	cell_t* orig_board;
	orig_board = init_board(orig_board); // NOTE: Is this correct C Syntax?
	copy_mat2board(originalGrid, orig_board); // Simply copies values etc.
	restore_invariants(orig_board);

	// TESTING CODE IN HERE
	
	// TESTING CODE ENDS HERE
	
	solved = heuristic_solve(orig_board);
	if( solved != PARTIAL_SOLN ) {
		return board2mat(orig_board);
	}

	workpool = init_workpool(workpool, INITIAL_WORKPOOL_SZ); // NOTE: These two #define macros need to be set.
	store 	 = init_store(store, INITIAL_STORE_SZ);

	int thr_solved; // Private variable
	cell_t* thr_board;
	#pragma omp parallel shared(solved, solved_board) private(thr_solved, thr_board)
	{
		// TODO: the entire thread algorithm
	}
	return board2mat(solved_board);
}

int dfs(cell_t* board) {
	// TODO: Write dfs.
}

int heuristic_solve(cell_t* board) {
	// TODO
}

int branch(cell_t* board, workpool_t* pool, int maxpushes) {
	// Initialize the list of children you can branch on.
	// 		The array of pointers to stacks is arr_bhelp
	// 		The array of sizes of aforementioned stacks in arr_bhsz
	int** arr_bhelp;	// points to a stack of max-size BOARDSIZE... ( don't really need that much... BUT consider empty board.)
	int*  arr_bhsz;		// points at the topmost NON-garbage value in arr_bhelp;
	int idx_h;
	arr_bhelp = malloc((SIZE+1)*sizeof(int*));
	arr_bhsz  = malloc((SIZE+1)*sizeof(int));
	for(idx_h=0; idx_h < (SIZE+1); ++idx_h) {
		arr_bhsz[idx_h] = -1; // Must do a ++ AFTER pushing.
		arr_bhelp[idx_h] = malloc(BOARD_SIZE*sizeof(int)); // Allocating a lot of memory here. O(size^3)
		// memset(arr_bhelp[idx_h], 0, BOARD_SIZE*sizeof(int)); // Don't need this. Not in my invariant.
	}
	// Now to scour the board and update values.
	
	for(idx_h=1; idx_h<BOARD_SIZE; ++idx_h) {
		// ... scour the board to initialize arr_bhelp and arr_bhsz
		int nal = board[idx_h].n_allowed;
		if(nal==0) continue; // Don't push these.
		else {
			arr_bhsz[nal]++;
			arr_bhelp[nal][arr_bhsz[nal]] = idx;
		}
	} // ASSERT: Done initializing.


	// Once initialized, greedy branch on indices with least amount of possible values which appear first.
	// However, I must also keep in mind the number of pushes I am making.
	int idx, branch_on, mypushes;
	mypushes = 0;
	// To branch on a given index on the board, we need to iterate through its allowed values
	// 		For each allowed value, we push a newboard onto the workpool.
	// 		This new board would first be a copy of the board.
	// 		We would then edit the board to set the value of the board correctly.
	// 			DO A SET VALUE ON THE BOARD.
	int iter_val;
	cell_t* new_board;
	for(idx = 0; idx<(SIZE+1); ++idx) {
		while(arr_bhsz[idx]>=0) {
			branch_on = arr_bhelp[idx][ arr_bhsz[idx] ];
			for(iter_val = 1; iter_val< board[branch_on].n_allowed; ++iter_val) { // Silly 1 indexed array.
				new_board = push_board(pool, board);
				mypushes++;
				set_board_value_idx(new_board, branch_on, new_board[branch_on].list[iter_val], NULL, NULL);
			
			}
			arr_bhsz[idx]--; // Reduce stack pointer.
			// ASSERT: pushed all children at branch_on	.
			if( mypushes>=maxpushes || pool->sz >= MAX_WORKPOOL_SZ ) { break; }
		}
		if( mypushes>=maxpushes || pool->sz >= MAX_WORKPOOL_SZ ) { break; } // Inner break only breaks out of while loop.
	}
}