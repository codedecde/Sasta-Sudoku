#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include <assert.h>
#include "sudoku.h"
// #define SIZE 9
#define BOARD_SIZE (SIZE*SIZE + 1)

#define NO_SOLN 2
#define PARTIAL_SOLN 1
#define SOLVED 0

#define INITIAL_POOL_SIZE 4
#define WORKPOOL_SIZE_BOUND 10

#define MAX_INT 100000
struct cell_t; // r, c, allowed, n_allowed
typedef struct cell_t cell_t;
	// Functions of cell_t
	void cell_init(cell_t* this);
	int is_allowed(cell_t* this, int val);
	int add_allowed(cell_t*, int); // returns -1 on behaviour dependent on this having a non-zero value
	int remove_allowed(cell_t* this, int val); // returns -1 on behaviour dependent on this having a non-zero value
	int set_value(cell_t*, int);
	int cell_destructor(cell_t*);
	cell_t* deepcopy_board(cell_t* this);
	
	void print_cell(cell_t*);
	void print_board(cell_t*);

struct pool_t; // Thread's access this... array of board's, FIFO/LIFO/PQ behaviour to be defined.
	struct pool_node_t; // private struct for pool_t really.
typedef struct pool_t pool_t;
	typedef struct pool_node_t pool_node_t;
	
	// Functions of : pool_t
	int is_empty(pool_t* this);
	void free_pool(pool_t* this);
	cell_t* pushf(pool_t*, cell_t*);
	cell_t* pushb(pool_t*, cell_t*);
	cell_t* popf(pool_t*);
	cell_t* popb(pool_t*);
	cell_t* pop(pool_t*);

	#define THRESHOLD_FOR_BRANCH 10
	int branch(pool_t* workpool, cell_t* board, int nbranches);
	int dfs(cell_t* board);
	int heuristic_solve(cell_t* board);	
struct cell_t {
	int row;
	int col;
	int list[SIZE + 1];
	int base_array[SIZE + 1];
	int n_allowed;
	int value;
}; // This is ok.
	void cell_init(cell_t* this) {
		//this->list = (int*) malloc( (SIZE+1)*sizeof(int) ); // idx0 corresponds to zilch.
		//this->base_array = (int*) malloc( (SIZE+1)*sizeof(int) ); // idx0 is useless... value 0 means not allowed.
		memset(this->base_array, 0, (SIZE+1)*sizeof(int));
		memset(this->list, 0, (SIZE+1)*sizeof(int)); // Not guaranteed, but no problem in doing it.
		this->n_allowed = 0;
		this->value = 0;
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

			return 1; // On success?
		} else {
			this->value = val;
			assert(this->n_allowed == 0);
		}
		return -1;
	}
	
	// WHAT WOULD HAVE BEEN A GORGEOUS cell_t struct :'(

	/*
		Serves the deep copy purpose when cell_t* definition was:
		struct cell_t {
			int row;
			int col;
			int list[SIZE+1];
			int base_array[SIZE+1];
			int n_allowed;
			int value;	
		}
	 */
	cell_t* deepcopy_board(cell_t* this) {
		cell_t* other = (cell_t*)malloc(BOARD_SIZE*sizeof(cell_t));
		int i = 1;
		// print_board(this);
		for(; i<BOARD_SIZE; ++i) {
			// Copy the i'th
			other[i].row = this[i].row;
			other[i].col = this[i].col;
			other[i].n_allowed = this[i].n_allowed;
			other[i].value = this[i].value;
			memmove(other[i].list,this[i].list,(SIZE+1)*sizeof(int));
			memmove(other[i].base_array,this[i].base_array,(SIZE+1)*sizeof(int));  
		}
		return other;
	}
struct pool_node_t {
	cell_t* board_config; // malloc-ing
	pool_node_t* next;
	pool_node_t* prev;
};
struct pool_t { // A doubly linked list. Allows for FIFO/LIFO.
	pool_node_t* hd;
	pool_node_t* tl;
	int sz;
	int turn;
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
	cell_t* pushf(pool_t* this, cell_t* new_board) { 
		// ASSERT: cell_t* is a fixed size board etc.
		struct pool_node_t* new_node = (pool_node_t*) malloc(sizeof(pool_node_t));

		new_node->board_config = deepcopy_board(new_board);
		new_node->next = this->hd;
		new_node->prev = NULL;
		if(this->hd) { // this->hd is not NULL
			this->hd->prev = new_node;
		}
		this->hd = new_node;
		this->sz++;
		if(this->sz == 1) {
			this->tl = this->hd;
		}
		return new_node->board_config;
	}
	/*
	Invocation time state : [...,p,q,r,s]
	pushb(x)
	Return time state: [...,p,q,r,s,x]
	*/
	cell_t* pushb(pool_t* this, cell_t* new_board) {
		struct pool_node_t* new_node = (pool_node_t*) malloc(sizeof(pool_node_t));
		
		new_node->board_config = deepcopy_board(new_board);
		new_node->next = NULL;
		new_node->prev = this->tl;
		if(this->tl) { // this->tl is not NULL
			this->tl->next = new_node;
		}
		this->tl = new_node;
		this->sz++;
		if(this->sz == 1) {
			this->hd = this->tl;
		}
		return new_node->board_config;
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
	cell_t* pop(pool_t* this){
		assert(this);
		this->turn = (this->turn + 1) % 2;
		return (this->turn % 2) ? popb(this) : popf(this);
	}
	// END OF WHAT WOULD'VE BEEN AN AMAZING (dead)pool.
	
	// ugliest mf ever... returns -1 if number of possible branches < THRESHOLD_TO_BRANCH (macro defined.)
	int branch(pool_t* this, cell_t* board, int max_pushes) {
		// The following arrays will help me push into the workpool.
		if(this->sz > WORKPOOL_SIZE_BOUND) { // WORKPOOL TOO LARGE.
			return -1;
		}
		int * initialization_helper[SIZE+1]; // idx-i is the array (int*) of indices in board which have i allowed values
		int  initialization_helper_sizes[SIZE+1]; // idx-i is the number of indices in board which have i allowed values
		int mypushes = 0;
		memset(initialization_helper_sizes, 0, (SIZE+1)*sizeof(int));
		memset(initialization_helper, 0, (SIZE+1)*sizeof(int*));

		// Iterate through board and fill up initialization_helper.
		int idx_board;
		for(idx_board=1; idx_board < BOARD_SIZE; ++idx_board) {
			int nvalues = board[idx_board].n_allowed;
			if( nvalues == 0 ) continue; // No values are allowed on this idx_board... ensures that initialization_helper[0] is empty.

			if (initialization_helper_sizes[nvalues]==0) { // hasn't been allocated yet.
				initialization_helper[nvalues] = (int*)malloc((BOARD_SIZE)*sizeof(int)); // This is 0 indexed.
			}
			initialization_helper[nvalues][initialization_helper_sizes[nvalues]] = idx_board;
			initialization_helper_sizes[nvalues]++;
		}
		// Check if the number of possible branches is less than some value.
		
		int hma = 0; 
		int count=0;
		for( ;hma<=SIZE;++hma) {
			count += hma* initialization_helper_sizes[hma];
		}
		if( count<THRESHOLD_FOR_BRANCH ) {
			// Cleaning up after myself...
			for(hma = 0; hma<=SIZE; ++hma) {
				free(initialization_helper[hma]);
			}
			return -1; // Don't branch
		}
		// int idx = 0;
		// for (idx = 0;idx <=SIZE;idx++){
		// 	 if(initialization_helper_sizes[idx] == 0)
		// 	   	assert(initialization_helper[idx] == NULL);
		// 	 else{
		// 	   	printf("Num Empty = %d\n",idx);
		// 	   	int jj = 0;
		// 	   	for(jj=0;jj<initialization_helper_sizes[idx];jj++)
		// 	     	printf("%d %d\t",(initialization_helper[idx][jj] -1) / SIZE, (initialization_helper[idx][jj]-1) % SIZE);
		// 	   	printf("\n");
		// 	 }
		// }
		
			// Iterates initialization_helper ... values between 2 and SIZE... Note: 0,1 are also included, but they're likely to be empty.
		int idx_hma = 0; 
			// Iterates initialization_helper[hma] ... values between 0 and initialization_helper_sizes[hma]
		int idx_cell_list;
		for(hma=0; hma<=SIZE; ++hma) {
			

			for( idx_hma=0; idx_hma < initialization_helper_sizes[hma]; ++idx_hma) {
				int branch_on = initialization_helper[hma][idx_hma];

				for(idx_cell_list=1; idx_cell_list<=board[branch_on].n_allowed; ++idx_cell_list) {
					// Some deep shit.
					// critical if in 
					cell_t* newboard;
					if( omp_in_parallel() !=0 ) {
						#pragma omp critical
						{
							newboard = pushf(this, board); // NO FRIGGING CLASS >.<	
						}
					} else {
						newboard = pushf(this, board); // NO FRIGGING CLASS >.<		
					}
					
					mypushes++;
					
					set_value(&((newboard)[branch_on]), newboard[branch_on].list[idx_cell_list]);
					//print_board(newboard);
				}
				if( mypushes >= max_pushes ) break;
			}
			if( mypushes >= max_pushes ) break;
		}

		// Cleaning up after myself...
		for(hma = 0; hma<=SIZE; ++hma) {
			free(initialization_helper[hma]);
		}
		// make initialization helper etc.
		return 0;
	}

	/*
		Convert a cell_t* board to an int** board
		return: int** 
	*/
	int** cell2board(cell_t* board){
		int** new_board = malloc(SIZE * sizeof(int*));
		int i;
		for(i = 0;i < SIZE;++i){
			new_board[i] = malloc(SIZE*sizeof(int));
			int j;
			for(j = 0;j < SIZE;++j){
				int idx = SIZE*i + j + 1;
				new_board[i][j] = board[idx].value;
			}
		}
		return new_board;
	}

	void print_cell(cell_t* cell) {
		printf("Printing cell[%d,%d]=%d\n", cell->row, cell->col, cell->value);
		printf("list:");
		int i;
		for (i=0; i<=SIZE; ++i) {
			printf("%d ", cell->list[i]);
		}
		printf("\nbase_array:");
		for (i=0; i<=SIZE; ++i) {
			printf("%d ", cell->base_array[i]);
		}
		printf("\nn_allowed=%d\n", cell->n_allowed);
	}
	/*
		print_board function launches a nuke. Obviously.
	 */
	void print_board(cell_t* board) {
		int r,c;
		printf("\n\t\tprinting board\n");
		for(c=-1; c<SIZE; c++) {
			if(c==-1) printf("  ");
			else printf("%d ", c); //TODO : add %02
		} 
		printf("\n");
		for(r=0; r<SIZE; r++) {
			printf("%d ",r);
			for(c=0; c<SIZE; c++) {
				printf("%d ", board[r*SIZE + c + 1].value);
			}
			printf("\n");
		}
		printf("\n\n");
	}


/*
  Function call made from main reaches here
*/
 void debug_print(cell_t* board){
 	// Check the initializations
	 int ii;
	 for(ii = 1; ii < BOARD_SIZE;ii++){
	 	if(board[ii].n_allowed == 0){
	 		assert(board[ii].value != 0);
	 		printf("Value of (%d,%d) = %d\n",board[ii].row,board[ii].col,board[ii].value);
	 	}
	 	else{
	 		int jj;
	 		assert(board[ii].value == 0);
	 		printf("Values allowed for (%d,%d) = ",board[ii].row,board[ii].col);
	 		for(jj = 1; jj <= board[ii].n_allowed; jj++){
	 			printf("%d ",board[ii].list[jj]);
	 		}
	 		printf("\n");

	 	}
	 }
	// All initializations are correct.
 }
 void print_workpool(pool_t* workpool){
 	pool_node_t* iter = workpool->hd;
 	int ctr=0;
 	while(iter != NULL){
 		printf("ctr=%d\n",ctr++);
 		print_board(iter->board_config);
 		iter = iter->next;
 	}
 }
int** solveSudoku(int** originalGrid){
	cell_t board[BOARD_SIZE];
	int r,c;
	for(r = 0;r < SIZE;r++){
		for(c = 0; c < SIZE; c++){
			int idx = r*SIZE + c + 1;
			cell_init(&board[idx]);
			board[idx].row = r;
			board[idx].col = c;
			if(originalGrid[r][c] != 0){
				set_value( &board[idx], originalGrid[r][c] );
			} else {
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
	// board is the original problem.
	
	
	 // cell_t* newboard = deepcopy_board(board);
	 // dfs(newboard);
	 // print_board(newboard);
	 // printf("\n\n--------------------------------------------------\n\n");
	 // print_board(board);

	// Initialize pool_t
	
	
	//print_workpool(workpool);

	int error = PARTIAL_SOLN;
	error = dfs(board); // Will become a shared variable soon. TODO : Uncomment.
	if (error == SOLVED || error == NO_SOLN) {
		return cell2board(board);
	} // Simplify the problem... without parallelization man.
	

	pool_t* workpool = (pool_t*) malloc(sizeof(pool_t)); // Really just three integers.
	workpool->hd = NULL; 
	workpool->tl = NULL; 
	workpool->sz = 0;
	workpool->turn = 0;
	branch( workpool, board, INITIAL_POOL_SIZE );	//Assert whether branch occurs correctly. Also Deepcopy
	
	int solution_found;
	solution_found = 0;
	cell_t* solved_board = NULL;
	#pragma omp parallel shared(solution_found,solved_board) private(error)
	{
		error = PARTIAL_SOLN;
		cell_t* thread_board;	//Local boards
		while(!solution_found) {
			
			#pragma omp critical
			{
				thread_board = (cell_t*)pop(workpool); // TODO: Check POP
			}
			if(solution_found) break;
			error = heuristic_solve(thread_board); // error is a private variable;
			if(error == SOLVED) {
				
				// #pragma omp critical
				// {
				// 	if(solved_board != NULL)
				// 		solved_board = thread_board;
				// }

				solution_found = 1; // Doesn't need to be critical. it is set to 1 afterall.
				break;
			}
			else if(error == NO_SOLN){
				free(thread_board);
				thread_board = NULL;
				continue;
			}

			if(branch(workpool, thread_board, 5) == -1) { // This function is thread_safe inside the 
				error = dfs(thread_board);
			}
			if( error == SOLVED ) {
				// #pragma omp critical
				// {
				// 	if(solved_board != NULL)
				// 		solved_board = thread_board;
				// }
				solution_found = 1; // Doesn't need to be critical. it is set to 1 afterall.
				break;
			}
			free(thread_board);
			thread_board = NULL;		
		
		}	//End while
		
		// Got here only by breaking or by solution_found being true...
		if(error == SOLVED) { // This thread found the solution.
			// solved_board must've been set.
			#pragma omp critical
			{
				if(!solved_board)
					solved_board = deepcopy_board(thread_board);
			}
		}
		free(thread_board); 
		/*else { // my thread_board is NOT a solution.
			free(thread_board);
		}*/
	}
	return cell2board(solved_board);

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
				  	// else{ // Don't need this block. Don't keep it.
				  	// 	continue;
				  	// }
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
		int error_code = dfs(board);
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
