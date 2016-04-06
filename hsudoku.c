#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include <assert.h>
#include <limits.h>
#include <omp.h>	// OpenMP
#include "sudoku.h"

int mypow(int base, int power) { // Because C confused me .
	int i;
	int acc;
	acc = 1;
	for(i=0; i<power; ++i) {
		acc = acc*base;
	}
	return acc;
} 


#define BOARD_SIZE (SIZE*SIZE + 1)

#define NO_SOLN 2
#define PARTIAL_SOLN 1
#define SOLVED 0
#define TIME_TO_LEAVE 42

#define INITIAL_WORKPOOL_SZ 6
#define MAX_WORKPOOL_SZ 40
#define INITIAL_STORE_SZ 100*BOARD_SIZE
#define NCELL_POPULATE_WORKPOOL 2

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
		if(this->value == 0) {
			int val;
			for(val=1; val<=SIZE; ++val) {
				this->base_array[val] = val;
				this->list[val] = val;
			}
			this->n_allowed = SIZE;
			return 0;
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
		
		if( altered_indices!=NULL && ptr_naltered!=NULL ) {
			memset(altered_indices, 0, BOARD_SIZE*sizeof(int));
			*ptr_naltered = -1;
			set_value(&board[idx], val);
			for(i=0; i<SIZE; ++i) {
				ir = r*SIZE + i + 1;
				ic = i*SIZE  + c + 1;
				ib = (rp + (i/MINIGRIDSIZE))*SIZE + (cp + (i%MINIGRIDSIZE)) + 1;
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
				ir = r*SIZE + i + 1;
				ic = i*SIZE  + c + 1;
				ib = (rp + (i/MINIGRIDSIZE))*SIZE + (cp + (i%MINIGRIDSIZE)) + 1;
				remove_allowed(&board[ir], val);
				remove_allowed(&board[ic], val);
				remove_allowed(&board[ib], val);
			}
		}
		return 1; // This function won't really fail, so its ok.
	}
	int set_board_value_idx(cell_t* board, int idx, int val, int* altered_indices, int* ptr_naltered) {
		assert(idx>=1);
		return set_board_value_rc(board, ((idx-1)/SIZE), ((idx-1)%SIZE), val, altered_indices, ptr_naltered);
	}
	int undo_alterations(cell_t* board, int val, int* altered_indices, int* ptr_naltered) {
		// ASSERT: altered_indices is a stack of things
		// *ptr_naltered is the list of values that has been changed. Their val was removed from allowed-ness.
		
		while(*ptr_naltered>=0) {
			add_allowed(&board[ altered_indices[*ptr_naltered] ], val);
			(*ptr_naltered)--;
		}
		// free(altered_indices);
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
			// free(&board[i]);
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
		return 0;
	}
	/*
		For outputting purposes.
	 */
	int** board2mat(cell_t* board) {
		int** mat = malloc(SIZE*sizeof(int*));
		int i,j;
		for(i = 0; i<SIZE; ++i) {
			mat[i] = malloc(SIZE*sizeof(int));
			for(j = 0; j<SIZE; j++) {
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
cell_t* solved_board;		// global/shared solution.
int solved;				// global/shared boolean


void destory_pool(workpool_t* pool) {
	int i;
	for(i=0; i<pool->sz; ++i) {
		destroy_board(pool->arr_boards[i]);
	}
	free(pool->arr_boards);
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
	pool->sz = initsz;
	return pool;
}

void destory_store(board_store_t* pool) {
	destory_pool((pool_t*)pool);
}

int not_full(pool_t* pool) { // Applies to both, board_store_t and workpool_t
	return ((pool->sz < pool->sz_alloc)?1:0);
}

int not_empty(pool_t* pool) {
	return (pool->sz > 0? 1:0);
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
	if( not_empty(pool) ) {
		return pop_board(pool);
	} else {
		// Make a new board...allocate it, and return it.
		printf("Store is empty\n");		
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
		// IMPORTANT: I don't do a memset here because I make no guarantees on pool->arr_boards[pool->sz...] anyway
		memcpy(pool->arr_boards, oldarr, pool->sz_alloc*sizeof(cell_t*));
		pool->sz_alloc = 2*pool->sz_alloc;
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
		pool->sz_alloc = 2*pool->sz_alloc;
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

int get_branchon(cell_t* board); // Behaviour not specificed... returns index of cell to branch on.
int populate_workpool(workpool_t* pool, cell_t* board, int nsize);
int hdfs(cell_t** ptr_board); // As discussed on 6th april.
int dfs_recursive(cell_t* board); // Take board ka dfs solution to completion
int dfs_iterative(cell_t* board);
int heuristic_solve(cell_t* board); // Take board and solve it heuristically as much as possible.
int branch(cell_t* board, workpool_t* pool, int max_branches); // Find the children of board. Push STRICTLY LESS(EQ) than max_branches

int get_branchon(cell_t* board) {
	int i;
	for(i=1; i<BOARD_SIZE; ++i) {
		// if board[i] is fillable.
		if( board[i].n_allowed > 0 ) {
			if( board[i].value == 0) {
				return i;
			} else {
				printf("Invalid board... in get_branchon\n");
				assert(0);
			}
		}
	}
	return -1;
}

int populate_workpool(workpool_t* pool, cell_t* board, int ncells) {
	assert(pool->sz == 0); // ASSUMES that board is empty.
	push_board(pool, board);
	int iter; // Run 0...(ncells-1) iterations.

	int stacktop =-1;
	cell_t** boardstack;
	boardstack = malloc( mypow(SIZE, ncells)*sizeof(cell_t*) );
	cell_t* topush, *tobranch, *newbranch;
	int bon;
	int val_iter,nal;
	for(iter=0; iter<ncells; ++iter) {
		while( pool->sz > 0) {
			topush = pop_board(pool); // No copying happens here, everything is cool.
			
			stacktop++;	boardstack[stacktop] = topush; // This is a push.
		}

		// ASSERT: pool is empty right now.
		while( stacktop > -1) { // taking stuff from stack, branching on each and pushing onto workpool.
			tobranch = boardstack[stacktop]; stacktop--; // This is a pop.
			bon = get_branchon(tobranch);
			if(bon == -1) {
				// We cn't branch no more...
				// Just put the stack back into the pool.
				push_board(pool, tobranch); // This is slow, TODO: think about replacing it by push_mem.
				push_mem(store, tobranch);	// IMPORTANT: This must accompany the push_board;

				while(stacktop > -1) { // Push all the others too.
					tobranch = boardstack[stacktop]; stacktop--;
					push_board(pool, tobranch);	
					push_mem(store, tobranch); // IMPORTANT: This must accompany push_board.
				}
				free(boardstack);
				return -1; // A problem arose.
			}
			// Do da branch.
			nal = tobranch[bon].n_allowed;
			for(val_iter =1 ; val_iter<=nal; ++val_iter) {
				// Alter tobranch[bon]
				newbranch = push_board(pool, tobranch);		//Does a copy and whatnot.
				set_board_value_idx(newbranch, bon, newbranch[bon].list[val_iter], NULL, NULL);
			}
			push_mem(store, tobranch); // Free da branch.
		}
	}
	free(boardstack);
	return pool->sz;
}

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
		r = (i-1)/SIZE;
		c = (i-1)%SIZE;
		rp = r - (r%MINIGRIDSIZE);
		cp = c - (c%MINIGRIDSIZE);
		for(j=0; j<SIZE;++j) {
			jr = r*SIZE + j + 1;
			jc = j*SIZE + c + 1;
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
	
	workpool = init_workpool(workpool, INITIAL_WORKPOOL_SZ); // NOTE: These two #define macros need to be set.
	store 	 = init_store(store, INITIAL_STORE_SZ);
	
	// Done only on a board that has had heuristic_solve run on it.
	
	cell_t* orig_board;
	orig_board = init_board(orig_board); // NOTE: Is this correct C Syntax?
	copy_mat2board(originalGrid, orig_board); // Simply copies values etc.
	restore_invariants(orig_board);
	// TESTING CODE IN HERE
	
	// solved = hdfs(&orig_board); // Solves it.
	// if( solved != PARTIAL_SOLN ) {
	// 	return board2mat(orig_board);
	// }
	// printf("Holy fuck no\n");
	// assert(0);
	
	// TESTING CODE ENDS HERE
	
	// TODO: Do we want to run heuristic before branches?
	// solved = heuristic_solve(orig_board);
	// if( solved != PARTIAL_SOLN ) {
	// 	return board2mat(orig_board);
	// }
	
	//TODO: Test populate_workpool
	populate_workpool(workpool, orig_board, NCELL_POPULATE_WORKPOOL); // Branch on first ncells.
	
	push_mem(store, orig_board);

	int thr_solved; // Private variable
	cell_t* thr_board;
	if( omp_get_num_threads() > workpool->sz ) { // We don't need as many threads.
		omp_set_num_threads(workpool->sz); // NOTE: High hopes, brun ... maybe even low hopes.
	}
	
	solved_board = NULL;
	#pragma omp parallel shared(solved, solved_board) private(thr_solved, thr_board)
	{
		// TODO: the entire thread algorithm
		while(not_empty(workpool)) { // The number of threads is dynamic, I'd rather avoid that for now.

			#pragma omp critical
			{
				thr_board = pop_board(workpool);
			}
			if (thr_board == NULL) {
				break; // Out of loop.
			}
			if(solved == SOLVED)
				break;
			int error_code = hdfs(&thr_board);
			if(error_code == SOLVED){
				#pragma omp critical
				{
					if(solved_board == NULL){
						solved_board = thr_board;
					}
				}
				solved = SOLVED;
				break;
			}
		}
	}
	
	destory_workpool(workpool);
	destory_store(store); // TODO: Check if positioning is appropriate.
	if (solved_board==NULL) {
		return originalGrid;
	} else {
		return board2mat(solved_board);
	}
}

// TODO: Test this.
// save oldstate (use store)
// run heuristic on board
// find index/cell to branch on - called bon
// 		for each value that bon can take:
// 			setval(bon,value)
// 			hdfs(board)
// 			if solved - woohoo.
// 				free oldstate ( use store )
// 				return board
// 			restore board to oldstate ( use deepcopy without malloc... copy_board )
// 		Get here means that no value could solve the board.
// 		restore oldstate...
// 			do free_mem( board) ... use store
// 			board = oldstate.
int hdfs(cell_t** ptr_board) {
	cell_t* oldboard;
	oldboard = pop_mem(store); // Global variable?
	
	copy_board(oldboard, *ptr_board); // returns 0... pretty much always
	int solv_state = -1;
	solv_state = heuristic_solve(*ptr_board);
	
	if(solv_state == SOLVED) {
		// Get rid of oldboard
		// printf("This happened\n");
		push_mem(store, oldboard);
		return SOLVED;
	} else if(solv_state == NO_SOLN) {

		// Restore oldboard and return.
		// printf("That happened\n");
		push_mem(store, *ptr_board);
		*ptr_board = oldboard;
		return NO_SOLN;
	}
	// PARTIAL_SOLN
	int bon, nal, lst_iter;
	bon = get_branchon(*ptr_board);
	nal = ((*ptr_board)[bon]).n_allowed;
	// ASSERT: bon can be branched on.
	int *alters, nalters; //nalters is an int.
	alters = malloc(BOARD_SIZE*sizeof(int));
	for(lst_iter=1; lst_iter<=nal; ++lst_iter) {
		nalters = -1; // These are used recursively. Hence, they can't be declared globally.
		
		int val = (*ptr_board)[bon].list[lst_iter];
		set_board_value_idx(*ptr_board, bon, val , alters, &nalters);
		
		// Before you recurse, check if somebody else has solved it.
		if (solved == SOLVED) { // 'Tis a global value.
			push_mem(store, *ptr_board);
			*ptr_board = oldboard;
			free(alters);
			return TIME_TO_LEAVE;
		}
		solv_state = hdfs(ptr_board);

		if(solv_state == SOLVED) {
			push_mem(store, oldboard);
			free(alters);
			return SOLVED;
		} else if(solv_state == TIME_TO_LEAVE) {
			// Free memory and leave.
			push_mem(store, *ptr_board);
			*ptr_board = oldboard;
			free(alters);
			return TIME_TO_LEAVE;
		}
		// copy_board(*ptr_board, oldboard); ITS WRONG ... hdfs restores to the board it found.
		undo_alterations(*ptr_board, val, alters, &nalters);
	}
	free(alters);
	push_mem(store, *ptr_board);
	*ptr_board = oldboard;
	return NO_SOLN;
}

int heuristic_solve(cell_t* board) {
	int flag = 1;
  	int err_code = -1;
  	while(flag){ 
	  	// While any changes are possible.
		// Eliminate
		// And lone range
	  	flag = 0;
	  	err_code = SOLVED;
	  	// ELIMINATION
	  	int r,c,idx;
	  	for(r=0;r<SIZE;r++){
		  	for(c=0;c<SIZE;c++){ // c++ ! Hah! Get it? :D
			  	int idx = r*SIZE + c + 1;
			  	if(board[idx].n_allowed == 0){
				  	if(board[idx].value == 0){ // value is still unfilled.
					  	err_code = NO_SOLN;
					  	return err_code;
				  	}
			  	} else if(board[idx].n_allowed == 1){
				  	flag = 1;
				  	int val = board[idx].list[1]; // 1 indexed array
				  	set_board_value_idx(board, idx, val, NULL, NULL);
			  	} else if( board[idx].n_allowed > 1){ // Made it a check for sanity purposes.
				  	err_code = PARTIAL_SOLN;
			  	}
		  	}
	  	}
  		
  		int row_count,col_count,box_count;
  		int row_idx,col_idx,box_idx;
  		int val;
  		for(val = 1;val <= SIZE;val++){
  			
	  		//Row iterator
  			// int ii = 0, jj = 0;
  			for(r = 0 ;r < SIZE; r++){
  				row_count = 0;
  				row_idx = 0;
	  			for(c = 0; c < SIZE;c++){
	  				int r_idx = r*SIZE  + c + 1;
	  				if(board[r_idx].value == val) {
	  					// row_idx = 0;
	  					break;
	  				} else if( is_allowed(&board[r_idx], val) ){
	  					row_count++;
	  					if(row_count == 2){
	  						// row_idx = 0;
	  						break;
	  					} else{ // row_count == 1 essentially.
	  						row_idx = r_idx;
	  					}
	  				}
	  			}
	  			if( row_count == 1 ){
	  				flag = 1;
	  				set_board_value_idx(board, row_idx, val, NULL, NULL);
	  			}
	  		}
	  		//Column Iterator
	  		for(c = 0; c < SIZE ; ++c){
	  			col_count = 0;
	  			col_idx = 0;
	  			for(r = 0; r < SIZE;r++){
	  				int c_idx = r*SIZE  + c + 1;
	  				if(board[c_idx].value == val) {
	  					break;
	  				} else if( is_allowed(&board[c_idx], val) ){
	  					col_count ++;
	  					if(col_count == 2){
	  						col_idx = 0;
	  						break;
	  					}else{ // col_count == 1
	  						col_idx = c_idx;
	  					}
	  				}
	  			}
	  			if(col_count == 1){
	  				flag = 1;
	  				set_board_value_idx(board, col_idx, val, NULL, NULL);
	  			}
	  		}
	  		//Box iterator
	  		int ii, jj;
	  		int r_prime, c_prime, box_count, box_idx, b_idx;
	  		for(ii = 0; ii < SIZE;ii ++ ){
	  			r_prime = ii - (ii % MINIGRIDSIZE);
	  			c_prime = MINIGRIDSIZE * (ii % MINIGRIDSIZE);
	  			box_count = 0;
	  			box_idx = 0;
	  			for(jj = 0; jj< SIZE;jj++){
	  				r = r_prime + (jj / MINIGRIDSIZE);
	  				c = c_prime + (jj % MINIGRIDSIZE);
	  				b_idx = (SIZE*r) + c + 1;
	  				if(board[b_idx].value == val) {
	  					break;
	  				}
	  				if( is_allowed( &board[b_idx],val) ){
		  				box_count ++;
		  				if(box_count == 2){
		  					// box_idx = 0;
		  					break;
		  				} else {
		  					box_idx = b_idx;
		  				}
				  	}
	  			}
	  			if(box_count ==1 ){
	  				flag = 1;
	  				// set_value(&board[box_idx],val);
	  				set_board_value_idx(board, box_idx, val, NULL, NULL);
	  			}
	  		}
  		}
  	}

	return err_code;
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
			arr_bhelp[nal][arr_bhsz[nal]] = idx_h;
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
			if( mypushes>=maxpushes ) { break; }
			if( pool->sz >= MAX_WORKPOOL_SZ ) { return -1; }
		}
		if( mypushes>=maxpushes || pool->sz >= MAX_WORKPOOL_SZ ) { break; } // Inner break only breaks out of while loop.
	}
	free(arr_bhsz);
	free(arr_bhelp);
	return 0;
}

// UNUSED CODE.
int dfs_recursive(cell_t* board){
	/*
	Currently iterate through the entire board
	Can we do better ?? ... TODO : yes... look at brute_force.c
	*/
	// ASSERT : input, i.e., board is valid.
	int min_val = SIZE+1; // Good enough.
	int min_idx = 0;
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
	// printf("#dfs: %d\n", min_idx);
	if(min_idx == 0){
		return SOLVED;
	}
	
	int num_allowed = board[min_idx].n_allowed;
	int r = (min_idx - 1) / SIZE;
	int c = (min_idx - 1) % SIZE;
	int r_prime = r - (r % MINIGRIDSIZE);
	int c_prime = c - (c % MINIGRIDSIZE);
	int base_array[SIZE+1];
	
	// print_cell(&board[min_idx]);
	memcpy((void*)base_array,(void*)board[min_idx].base_array,(SIZE+1)*sizeof(int));
	
	// printf("n_allowed=%d\n", num_allowed);
	// for(i=0; i<(SIZE+1); i++) {
	// 	printf("%d ", base_array[i]);
	// }
	// printf("\n");
	// base_array[0] = 0;

	int lst_indices[BOARD_SIZE];
	memset(lst_indices,0, BOARD_SIZE*sizeof(int));
	int lst_idx = -1; // It is set withing the loop too...
	// lst_idx = 0; ... initialized inside forloop.
	for(i = 1; i <= num_allowed;i++){
		int val = board[min_idx].list[i];
		int jdx;
		
		// ASSERTION CODE FOLLOWS
				int fl_valid_set_value = 1;
				for(jdx = 0; jdx<SIZE;++jdx) {
					int jdx_r = r*SIZE + jdx + 1;
				  	int jdx_c = jdx*SIZE + c + 1;
				  	int jdx_b = ((r_prime + (jdx / MINIGRIDSIZE))*SIZE) + (c_prime + (jdx%MINIGRIDSIZE) ) + 1;
				  	
				  	if( (board[jdx_r].value == val && (jdx_r != min_idx)/*!= 0 && is_allowed(&board[min_idx], board[jdx_r].value)*/) ||
				  		(board[jdx_c].value == val && (jdx_c != min_idx)/*!= 0 && is_allowed(&board[min_idx], board[jdx_c].value)*/) ||
				  		(board[jdx_b].value == val && (jdx_b != min_idx)/*!= 0 && is_allowed(&board[min_idx], board[jdx_b].value)*/) ) {
				  		
				  		fl_valid_set_value = 0; break;
				  	}
				}
				if( fl_valid_set_value!=1 ) {
					printf("\tASSERTION FAIL :<\n");
					print_board(board);
					printf("min_idx=%d ; val=%d\n", min_idx, val);
					print_cell(&board[min_idx]);
					int extremely_temp_iter;
					printf("base_array original:\n");
					for(extremely_temp_iter = 0; extremely_temp_iter<(SIZE+1); extremely_temp_iter++) 
						printf("%d ", base_array[i]);
					printf("\n");
					assert(fl_valid_set_value == 1);
				}
		
		set_value(&board[min_idx],val);
		jdx = 0;
		lst_idx = -1;
		for(;jdx < SIZE;jdx++){
			int jdx_r = r*SIZE + jdx + 1;
		  	int jdx_c = jdx*SIZE + c + 1;
		  	int jdx_b = ((r_prime + (jdx / MINIGRIDSIZE))*SIZE) + (c_prime + (jdx%MINIGRIDSIZE) ) + 1;
		  	if( is_allowed(&board[jdx_r],val) ) {
		  		assert(jdx_r != min_idx); // won't get here becase if condition stops it. set_value works.
		  		lst_idx++;
		  		lst_indices[lst_idx] = jdx_r;
		  		remove_allowed(&board[jdx_r],val);
		  	}
		  	if( is_allowed(&board[jdx_c],val) ) {
		  		assert(jdx_c != min_idx); 
		  		lst_idx++;
		  		lst_indices[lst_idx] = jdx_c;
		  		remove_allowed(&board[jdx_c],val);
		  	}
		  	if( is_allowed(&board[jdx_b],val) ) {
		  		assert(jdx_b != min_idx); 
		  		lst_idx++;
		  		lst_indices[lst_idx] = jdx_b;
		  		remove_allowed(&board[jdx_b],val);
		  	}
		}
		int error_code = dfs_recursive(board);
		if(error_code == SOLVED){
			return SOLVED;
		}
		for(; lst_idx>-1; lst_idx--) {
			add_allowed(&board[ lst_indices[lst_idx] ], val);
		}
	}
	board[min_idx].value = 0;
	board[min_idx].n_allowed = num_allowed;
	
	memcpy(board[min_idx].base_array,base_array,(SIZE+1)*sizeof(int));
	// free(base_array); ... doesn't need to be freed. Its a stack object.
	return NO_SOLN;
}


/*
int dfs_iterative(cell_t* board) {
	// TODO: Write iterative dfs.
	// INIT:
	// 		list of indices that need to be filled.
	//		depth which points to the index in <unfilled> that was filled.
	// 		list of value_iter
	// 		List of indices affected by 
	printf("def_iterative hasn't been implemented\n");
	assert(0);
	int* unfilled, value_iter;
	int max_depth;
	int** alterations; int* naltered;
	int depth, iter_depth; // Points directly at the top of the stack.
	int idx; 	// loop iterator.
	unfilled 	= malloc(BOARD_SIZE*sizeof(int)); // parameterized by iter_depth
	value_iter 	= malloc(BOARD_SIZE*sizeof(int)); // parameterized by iter_depth... remember that list is 1 indexed
	alterations = malloc(BOARD_SIZE*sizeof(int*));
	naltered 	= malloc(BOARD_SIZE*sizeof(int));
	depth = -1;
	max_depth = 0;
	// IMPORTANT: depth < max_depth 
	// STEP1: initialize unfilled.
	for(idx=1; idx<BOARD_SIZE; ++idx) {
		if( board[idx].n_allowed>0 ) {
			max_depth++;
			depth++;
			unfilled[depth] = idx;
		}
	}
	for(idx=0; idx<max_depth; ++idx) {
		naltered[idx] = -1;
		alterations[idx] = malloc(BOARD_SIZE*sizeof(int));
	}

	// ASSERT: INIT complete.
	// iterate through unfilled.
	iter_depth = 0;
	// 
	int branch_on =-1;
	for( ;iter_depth<=max_depth; ) {
		branch_on = unfilled[iter_depth];
		if( value_iter[iter_depth] == 0 ) {
			// Time to branch on this.
			
		} else if( value_iter[iter_depth] <= board[branch_on].n_allowed ) {
			// Time to move to next branch.
			// The present value didn't work, time to try the next value.
		} else {
			// Time to backtrack...
			// No value worked here.
		}
		iter_depth++;
	}
}
*/