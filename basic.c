/*2013CS10225_2013CS10773*/
#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include <assert.h>
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
#define BASE_SIZE SIZE*(SIZE+1)

/*
	TIME_TO_LEAVE is returned by HDFS whenever somebody has already found a solution. Hence, recursing is pointless.
 */
#define NO_SOLN -1
#define PARTIAL_SOLN 0
#define SOLVED 1
#define TIME_TO_LEAVE 42

#define INITIAL_STORE_SZ 100
#define NCELL_POPULATE_WORKPOOL 2
#define INITIAL_WORKPOOL_SZ 40

typedef struct {
	int val;
	int nal;
	int *lst;
	int *arr;
} cell_t;
void cell_init(cell_t* this) {
	this->lst = malloc((SIZE+1)*sizeof(int));
	this->arr = malloc((SIZE+1)*sizeof(int));
	this->nal = 0;
	this->val = 0;
}
int is_allowed(cell_t* this, int val) {
	return (((this->val == val) || (this->arr[val] > 0)) ? 1 : 0 );
}
int add_allowed(cell_t* this, int val) {
	if (this->val == 0) {
		if (this->arr[val] > 0) { // val is already allowed
			return 0;
		} else {
			this->nal++;
			this->lst[this->nal] = val;
			this->arr[val] = this->nal;
			return 1;
		}		
	} else { // value is nonzero, nothing is changed. base_array and nallowed must remain as earlier.
		return -1;
	}
}
int remove_allowed(cell_t* this, int val) {
	int idx_list = this->arr[val];
	if ( this->val != 0 ) {// This is the case where nothing is allowed anyway. It should be covered by other case.
		return -1; // Special behaviour, no?
	} 
	if ( idx_list==0 ) { // This is the case where the given val is not actually allowed.)
		return 0;
	}
	if (this->nal == idx_list) { // 1-indexed means that their equality is what we care about.
		this->nal--;
		this->arr[val] = 0; // The value is not allowed.
		return 1;
		// We make no guarantees about the this->list values at indices > n_allowed.
	} else {
		// Case 2: It is not the tail of the list
			// Subcase 1: It is the only element...in which case it'll also be the tail of the list.
			// Subcase 2: There is some other element that is the tail of the array.	
		// We only need to worry about subcase 2.
		// Swap the items at idx_list and this->n_allowed;
		// Note: THIS REORDERS THE FRIGGIN LIST. however, if my value has been set, it shouldn't matter.
		this->arr[val] = 0;
		this->lst[idx_list] = this->lst[this->nal]; // New value.
		this->arr[ this->lst[this->nal] ] = idx_list;
		this->nal--;
		if(this->nal <0 ) {
			printf("assertion problem tid=%d\n", omp_get_thread_num());
			assert(this->nal >= 0 );
		}
		return 1;
	}
}
int set_value(cell_t* this, int val) {
	if(this->val == 0) {
		this->val = val;
		this->nal = 0;
		memset(this->arr,0,(SIZE + 1)*sizeof(int));
		return 1; // On successfully unallowing things.
	} else {
		this->val = val;
		return 0; // 
	}
	return -1; // Should never get here.
}

int set_all_allowed(cell_t* this) { // Set all possible values to "allowed"...
	if(this->val == 0) {
		int val;
		for(val=1; val<=SIZE; ++val) {
			this->arr[val] = val;
			this->lst[val] = val;
		}
		this->nal = SIZE;
		return 0;
	} else {
		return -1; // Doesn't make sense to set everything to allowed.
	}
}
void print_board(cell_t* board) { // This function summons a level SSS vampire. Obviously.
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
			printf("%02d ", board[r*SIZE + c + 1].val);
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
cell_t* init_board(cell_t** board) {
	// initialize the board... allocate it, and its list and base_array
	// Note: boards are 1 indexed. Hence, board[0] = NULL
	*board = malloc(BOARD_SIZE*sizeof(cell_t));
	// Initialize the cells now.
	int i;
	// board[0] = NULL; ... How do I do this?
	for(i = 1; i<BOARD_SIZE; ++i) {
		(*board)[i].val = 0;
		(*board)[i].nal = 0;
		(*board)[i].lst = malloc((SIZE+1)*sizeof(int));
		memset((*board)[i].lst, 0, (SIZE+1)*sizeof(int));
		(*board)[i].arr = malloc((SIZE+1)*sizeof(int));	
		memset((*board)[i].arr, 0, (SIZE+1)*sizeof(int));
	}
	return *board;
}

int destroy_board(cell_t* board) { // Goes in deep and destroys the board entirely.
	int i;
	for(i=1; i<BOARD_SIZE; ++i) {
		free(board[i].lst);
		free(board[i].arr);
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
		dst[i].val = src[i].val;
		dst[i].nal = src[i].nal;
		// ASSERT: dst[i].list and dst[i].base_array have been allocated.
		memcpy(dst[i].lst, src[i].lst, (SIZE+1)*sizeof(int));
		memcpy(dst[i].arr, src[i].arr, (SIZE+1)*sizeof(int));
	}
	return 0;
}

int** board2mat(cell_t* board) {
	int** mat = malloc(SIZE*sizeof(int*));
	int i,j;
	for(i = 0; i<SIZE; ++i) {
		mat[i] = malloc(SIZE*sizeof(int));
		for(j = 0; j<SIZE; j++) {
			mat[i][j] = board[i*SIZE + j + 1].val;
		}
	}
	return mat;
}

typedef struct {
	cell_t** arr_boards; // It is an array of pointers, hence moving them around should be easier.
	int sz_alloc; 		 // Total size allocated
	int sz;				 // Present size ... points at 1st unoccupied position.

} pool_t;
typedef pool_t workpool_t;
typedef pool_t board_store_t;


pool_t* init_pool(pool_t** ptr_pl, int in_sz){
	*ptr_pl = malloc(sizeof(pool_t));
	(*ptr_pl)->arr_boards = malloc(in_sz*sizeof(cell_t*));
	(*ptr_pl)->sz = 0;
	(*ptr_pl)->sz_alloc = in_sz;
	return (*ptr_pl); // A pointer to a pool.
}
pool_t* init_store(board_store_t** store, int in_sz) {
	*store = malloc(sizeof(board_store_t));
	(*store)->sz_alloc = in_sz;
	(*store)->sz=0;
	(*store)->arr_boards = malloc(in_sz*sizeof(cell_t*));
	int i;
	for(i=0; i<in_sz; ++i) {
		(*store)->arr_boards[i] = init_board( &((*store)->arr_boards[i]) );
		((*store)->sz)++;
	}
	return *store;
}
void destroy_pool(pool_t* pool) { // Ideally, this function would take a graveyard and operate by oushing boards onto that.
	int i;
	for(i=0; i<pool->sz; ++i) {
		destroy_board(pool->arr_boards[i]);
	}
	free(pool->arr_boards);
	free(pool);
}

// Cheaty macro functions.
#define not_full(pool) \
	(pool->sz < pool->sz_alloc)
#define not_empty(pool) \
	(pool->sz>0)

cell_t* pop_board(workpool_t* wrkpl) {
	if( not_empty(wrkpl) ) {
		wrkpl->sz--;
		return wrkpl->arr_boards[wrkpl->sz];
	} else {
		return NULL;
	}
}
cell_t* pop_mem(board_store_t* store) { // Can be called with a store=NULL.
	cell_t* retval = NULL;
	if( store!=NULL && not_empty(store) ) {
		store->sz--;
		retval =  store->arr_boards[store->sz];
	} else {
		retval = init_board(&retval);
	}
	return retval;
}

cell_t* push_board(workpool_t* pool, cell_t* board, board_store_t* store) {
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
	pool->sz++;
	copy_board(retval, board); // Go in deep and copy stuff. O(n3) operation man.
	return retval;
}

int push_mem(board_store_t* pool, cell_t* board) { // pool is 'store'
	int retval = 0;
	if( !not_full(pool) ) { // is full :P
		cell_t** oldarr = pool->arr_boards;
		pool->arr_boards = malloc(2*pool->sz_alloc*sizeof(cell_t*));
		memcpy(pool->arr_boards, oldarr, pool->sz_alloc*sizeof(cell_t*));
		pool->sz_alloc *= 2;
		free(oldarr);
		retval = 1;
	}
	pool->arr_boards[pool->sz] = board;
	pool->sz++;
	return retval; // Not a pointer, but a int.
}

// NOTE: init functions, really.
int copy_mat2board(int** mat, cell_t* board) {
	// ASSERT: board has already been allocated.
	int r,c,idx;
	for(r=0; r<SIZE; ++r) {
		for(c=0; c<SIZE; ++c) {
			idx = r*SIZE + c + 1;
			board[idx].val = mat[r][c];
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
		if(board[i].val == 0) continue; // This imposes no constraints on anybody.
		r = (i-1)/SIZE;
		c = (i-1)%SIZE;
		rp = r - (r%MINIGRIDSIZE);
		cp = c - (c%MINIGRIDSIZE);
		for(j=0; j<SIZE;++j) {
			jr = r*SIZE + j + 1;
			jc = j*SIZE + c + 1;
			jb = (rp + (j/MINIGRIDSIZE))*SIZE + (cp +(j%MINIGRIDSIZE) ) + 1;
			remove_allowed(&board[jr], board[i].val);
			remove_allowed(&board[jc], board[i].val);
			remove_allowed(&board[jb], board[i].val);
		}
	}
	return 1;
}


int get_branchon(cell_t* brd) {
	int idx, minidx, minv;
	minv = 2*SIZE;
	minidx = -1;
	for(idx=1; idx<BOARD_SIZE; ++idx) {
		if( brd[idx].val ==0) {
			if(brd[idx].nal==0) {
				return -1; // Something's wrong here.
			} else if(brd[idx].nal < minv) {
				minidx = idx;
				minv = brd[idx].nal;
			}
		}
	}
	return minidx;
}
int has_unfilled(cell_t* brd) {
	// Outputs:
	// 		1,-1 take priority. Returns whatever is found first.
	// 		
	// 		1 if there's an unfilled cell. (primarily)
	// 		0 if there are no unfilled cells
	// 		-1 if there are "some" problems... doesn't detect all of them. 
	int idx;
	for(idx=1; idx<BOARD_SIZE; ++idx) {
		if(brd[idx].val==0) { // unfilled cell.
			if(brd[idx].nal==0) {
				return -1;
			} else {
				return 1;
			}
		}
	}
	return 0;
}
// IMPORTANT: Start here.

// NOTE: Switching gears here.
// Meaty functions follow.
// int** solveSudoku(int** originalGrid);
		int copy_mat2board(int** mat, cell_t* brd);
		int restore_invariants(cell_t* brd);
// 		
// 	
int populate_workpool(workpool_t*, cell_t*, int);			
int hdfs(cell_t**, int*, board_store_t*); // board_t*, int* ... the int* has already been allocated.
//		int get_branchon(cell_t*);
// 		
// 	RETURN: SOLVED, NO_SOLN, TIME_TO_LEAVE
int heuristic_solve(cell_t**, int*);
// 			int elimination(cell_t**)
// 			int loneranger(cell_t**, int*) 	// requires int* to be initialized.
//			int has_unfilled(cell_t* brd);

workpool_t* workpool; // Global
board_store_t* thr_store; // Thread local graveyards.

cell_t* solved_board;
int solved;
cell_t* thr_board;
int thr_solved;

int** solveSudoku(int** origmat) {
	solved_board = NULL;
	solved = PARTIAL_SOLN;

	cell_t* origbrd;
	origbrd = init_board(&origbrd);
	copy_mat2board(origmat, origbrd);
	restore_invariants(origbrd);

	// Quickly run heuristic... just to worry about it.
	int* arr_lr;
	arr_lr = malloc(3*BASE_SIZE*sizeof(int));
	solved = heuristic_solve(&origbrd, arr_lr);
	free(arr_lr); // cleanup.
	if(solved == SOLVED) {
		return board2mat(origbrd);
	} else if(solved==NO_SOLN) {
		// Free memory.
		destroy_board(origbrd);
		return origmat;
	}
	// ELSE : ... i.e, heuristic didn't solve it entirely.
	printf("Go beyond it...%d\n", solved);
	// Now, to initialize workpools.
	init_pool(&workpool, INITIAL_WORKPOOL_SZ);
	populate_workpool(workpool, origbrd, NCELL_POPULATE_WORKPOOL); // TODO: Make this function, without using the graveyard.
	// origbrd can be considered a memory leak.
	// ASSERT: workpool is ready.
	#pragma omp parallel shared(solved_board, solved, workpool) private(thr_board, thr_solved, thr_store)
	{
		// Code for each thread.
		// INIT:
		thr_solved = PARTIAL_SOLN;
		thr_board  = NULL; // Pop this from workpool soon.
		thr_store  = init_store(&thr_store, INITIAL_STORE_SZ);
		int* arr_thr_lr;
		arr_thr_lr = malloc(3*BASE_SIZE*sizeof(int));
		// EXECUTION:	
		while( not_empty(workpool) ) { // There is work to be done... BEGIN while(not_empty(wrkpl))
			// loop INIT:
			#pragma omp critical (workpool_lock)
			{
				thr_board = pop_board(workpool);
			}
			if(solved==SOLVED || solved==NO_SOLN|| thr_store == NULL) { // Happens only if workpool is empty.
				thr_solved = TIME_TO_LEAVE; // Get ready to break.
				break;
			}
			// loop EXECUTION:
			thr_solved = hdfs(&thr_board, arr_thr_lr, thr_store);
			
			if(thr_solved==TIME_TO_LEAVE) { // Somebody else has found a solution.
				break;
			} else if(thr_solved==SOLVED) {
				#pragma omp critical (solved_lock)
				{ // Critical because we are altering solved/solved_board. Sure, they'll always point to a important value, but memory leaks.
					if(solved==0) {
						solved = 1;
						solved_board = thr_board;
					}
				}
				if(solved_board == thr_board) {
					break;
				}
			} else if(thr_solved == PARTIAL_SOLN) {
				printf("HDFS returned PARTIAL_SOLN ... impossibru behaviour\n");
				assert(0);
			} else { // thr_solved == NO_SOLN
				// Nothing special to do. Simply repeat.
			}
			// loop CLEANSE
			push_mem(thr_store, thr_board);
		} // END while(not_empty(wrkpl))
		// CLEANSE:
		free(arr_thr_lr); // Unidimensional array, can be solved like so.
		arr_thr_lr = NULL;

		if(thr_solved == SOLVED && solved_board==NULL) { // If I was a solution.
			solved_board = thr_board; // And don't delete this.
		} else if(thr_solved == SOLVED && solved_board!=NULL){
			// Do absolutely nothing.
		}else { // Somebody else has already got a solution.
			// Simply destroy my board... would've pushed the board into thr_store, but that has been destroyed already.
			// push_mem(thr_store, thr_board);
			thr_solved = PARTIAL_SOLN;
			thr_board  = NULL;
		}
		destroy_pool(thr_store); 
		thr_store = NULL;
		
	}
	if(solved==SOLVED) {
		destroy_pool(workpool);
		return board2mat(solved_board);
	} else {
		return origmat;
	}
}

int hdfs(cell_t** ptrbrd, int* base, board_store_t*  store) {
	// IN: ptrbrd, base have both been allocated.
	cell_t* entrystate;

	int mysolved = PARTIAL_SOLN;
	
	if(solved==1) {
		return TIME_TO_LEAVE;
	}
	entrystate = pop_mem(store); // Called within a thread.
	copy_board(entrystate, *ptrbrd);
	if(solved==1) {
		push_mem(store, entrystate);
		return TIME_TO_LEAVE;
	}

	mysolved = heuristic_solve(ptrbrd, base); // Should return SOLVED if there are no values to fill.
	
	if(mysolved == SOLVED || mysolved==NO_SOLN || mysolved==TIME_TO_LEAVE) {
		push_mem(store, entrystate);
		return mysolved;
	}
	// Continue only if it wasn't solved.
	
	int bon, nal, val, iter, nalters;
	int* alters;

	bon = get_branchon(*ptrbrd);
	if(bon<=0) {
		// Nothing to branch on but PARTIAL_SOLN...
		// Something's wrong with heuristic_solve
		printf("hs returned a solution/fully-filled but not SOLVED/NO_SOLN\n");
		assert(0);
	}
	nal = (*ptrbrd)[bon].nal;
	alters = malloc(BOARD_SIZE*sizeof(int));
	nalters = -1;
	for(iter=1; iter<=nal; ++iter) {
		val = (*ptrbrd)[bon].lst[iter];
		set_board_value_idx(*ptrbrd, bon, val, alters, &nalters);
		mysolved = hdfs(ptrbrd, base, store); // IMPORTANT: Recursive call.
		if(mysolved==SOLVED) {
			push_mem(store, entrystate);
			return SOLVED;
		} else if(mysolved==TIME_TO_LEAVE) {
			push_mem(store, entrystate);
			return TIME_TO_LEAVE;
		}
		undo_alterations(*ptrbrd, val, alters, &nalters);
	}
	free(alters);
	// None of the branches found a child. NO_SOLN exists.
	push_mem(store, *ptrbrd);
	*ptrbrd = entrystate;
	return NO_SOLN;
}

int heuristic_solve(cell_t** ptrbrd, int* base) {
	// IMPORTANT: If at any point, we realize that the board has been solved, or there is NO_SOLN, we must return
	// Also, assume that "base" has already been allocated.
	int fl_change_made = 1;
	int unfilled_cell;
	
	int idx, val, nal, viter;
	int idx2r, idx2c, idx2b;
	int lr_type, lr_titer, biter; // also needs 'val' but can reuse.
	while(fl_change_made) {

		unfilled_cell = 0;	// I want the value of unfilled_cell after last elim loop.
		
		// Elimination.
		while(fl_change_made){ // eliminate as much as possible.
			fl_change_made = 0;
			unfilled_cell = 0; // I want the value of unfilled_cell after last elimination loop.
			for(idx=1; idx<BOARD_SIZE; ++idx) {
				if((*ptrbrd)[idx].nal == 0) {
					if((*ptrbrd)[idx].val == 0) {
						return NO_SOLN;
					}
				} else if((*ptrbrd)[idx].nal == 1) {
					if((*ptrbrd)[idx].val == 0) {
						fl_change_made = 1;
						val = (*ptrbrd)[idx].lst[1];
						set_board_value_idx(*ptrbrd, idx, val, NULL, NULL);
					} else {
						printf("heuristic solve: Board invariant violated\n");
						assert(0);
						return NO_SOLN;
					}
				} else {
					unfilled_cell = 1;	
				}
			}
		}
		// ASSERT: At this point, as many changes could be made, have been made.
		if(unfilled_cell==0) { // There are no more cells to fill, don't do lone ranger.
			return SOLVED;
		}
		fl_change_made = 0;
		
		// LONE RANGER ... ideally, we'd also have a lr_unfilled_cell. But eff it.
		memset(base, 0, 3*BASE_SIZE*sizeof(int));
		// We will now use the base array.
		// x = base[type*BASE_SIZE + iter_in_type*SIZE + val] is what we'll use to detect lone rangers.
		// 		type={ row:0, col:1, box/minigrid:2 }
		// if x==0, that type has not been found yet.
		// if x==idx>0, then there is one occurence of that type, and it must be set by loneranger
		// if x==-1, then there are two or more places where val is allowed.
		// 			 the other case of value already being there is an invalid
		for(idx=1; idx<BOARD_SIZE; ++idx) {
			nal = (*ptrbrd)[idx].nal;
			idx2r = (idx-1)/SIZE;
			idx2c = (idx-1)%SIZE;
			idx2b = (idx2r/MINIGRIDSIZE)*MINIGRIDSIZE + idx2c/MINIGRIDSIZE;
			if((*ptrbrd)[idx].val==0) {
				if(nal>0) {
					for(viter=1; viter<=nal; ++viter) {
						base[0 + idx2r*SIZE + val] = (base[0 + idx2r*SIZE + val] == 0? idx:-1);
						base[BASE_SIZE + idx2c*SIZE + val] = (base[BASE_SIZE + idx2c*SIZE + val] == 0? idx:-1);
						base[2*BASE_SIZE + idx2b*SIZE + val] = (base[2*BASE_SIZE + idx2b*SIZE + val] == 0? idx:-1);
					}
				} else {
					// Here's an unfilled cell that has no allowed values.
					return NO_SOLN;				
				}
					
			}
		} // ASSERT: base has now been prepared.


		// We need to keep track of fl_change_made and unfilled_cell
		// Dealing with rows first.
		lr_type = 0;
		for(biter=1; biter<(BASE_SIZE); ++biter) {
			// TODO... set the values and also the upper bound for biter... keep track of unfilled_cell
			lr_titer = biter/SIZE;
			val = biter%SIZE;
			if(val==0) continue;
			idx = base[biter];
			if ( idx > 0) {
				if( is_allowed(&(*ptrbrd)[idx], val)) { // returns false if already set.
					fl_change_made = 1;
					set_board_value_idx(*ptrbrd, idx, val, NULL, NULL);
				}
			}
		}
		// Dealing with columns
		lr_type = 1;
		++biter;
		for(/*biter = 1 + BASE_SIZE*/; biter< 2*BASE_SIZE; ++biter) {
			lr_titer = biter/SIZE;
			val = biter%SIZE;
			if(val==0) continue;
			idx = base[biter];
			if ( idx > 0) {
				if( is_allowed(&(*ptrbrd)[idx], val)) { // returns false if already set.
					fl_change_made = 1;
					set_board_value_idx(*ptrbrd, idx, val, NULL, NULL);
				}
			}
		}
		// Dealing with boxes
		lr_type = 2;
		++biter;
		for(/*biter=1+2*BASESIZE*/; biter<3*BASE_SIZE; ++biter) {
			lr_titer = biter/SIZE;
			val = biter%SIZE;
			if(val==0) continue;
			idx = base[biter];
			if ( idx > 0) {
				if( is_allowed(&(*ptrbrd)[idx], val)) { // returns false if already set.
					fl_change_made = 1;
					set_board_value_idx(*ptrbrd, idx, val, NULL, NULL);
				}
			}
		}
		unfilled_cell = has_unfilled(*ptrbrd);
		if(unfilled_cell==1) {
			// Do nothing.
		} else if(unfilled_cell == 0) {
			return SOLVED;
		} else { // invalid solution.
			return NO_SOLN;
		}
	
		// More heuristcs can follow here.
	}

	if(unfilled_cell==0) {
		return SOLVED;
	} else {
		return PARTIAL_SOLN;
	}
}

int populate_workpool(workpool_t* wpl, cell_t* brd, int ncl) {

	int iter, viter;
	int bon,nal,val;

	cell_t** stk;
	
	cell_t* tomv2stk;
	cell_t* tobrnch;
	cell_t* tomdfy;
	int stktp;
	stktp = -1;
	stk = malloc( mypow(SIZE, ncl)*sizeof(cell_t*) );
	
	push_board(wpl, brd, NULL);
	for(iter = 0; iter<ncl; ++iter) {
		// Branch on the iter'th value
		while( not_empty(wpl) ) {
			tomv2stk = pop_board(wpl);
			
			stktp++; 
			stk[stktp] = tomv2stk;
		}

		while(stktp>-1) {
			tobrnch = stk[stktp];
			stktp--;

			bon = get_branchon(tobrnch);
			nal = tobrnch[bon].nal;
			
			// IMPORTANT: hacky optimization follows.
			for(viter=1; viter< nal; viter++) {
				val = tobrnch[bon].lst[viter];
				tomdfy = push_board(wpl, tobrnch, NULL);
				set_board_value_idx(tomdfy, bon, val, NULL, NULL);
			}
			val = tobrnch[bon].lst[nal];
			push_mem(wpl, tobrnch); // Pushes the pointer in. Now, I don't need to delete tobrnch.
			set_board_value_idx(tobrnch, bon, val, NULL, NULL);
		}
	}
	assert(stktp == -1);
	free(stk);
}