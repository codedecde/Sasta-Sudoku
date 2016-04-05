#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include "sudoku.h"
// #define SIZE 9
#define BOARD_SIZE (SIZE*SIZE + 1)

#define NO_SOLN 2
#define PARTIAL_SOLN 1
#define SOLVED 0

#define INITIAL_POOL_SIZE 5
#define WORKPOOL_SIZE_BOUND 10

#define MAX_INT 100000
struct cell_t; // r, c, allowed, n_allowed
typedef struct cell_t cell_t;
	// Functions of cell_t
	void cell_init(cell_t* this);
	int is_allowed(cell_t* this, int val);
	void add_allowed(cell_t*, int);
	void remove_allowed(cell_t* this, int val);
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

	#define THRESHOLD_FOR_BRANCH 10
	int branch(pool_t* workpool, cell_t* board, int nbranches);
	int dfs(cell_t* board);
	int heuristic_solve(cell_t* board);	
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
		memset(this->base_array, 0, (SIZE+1)*sizeof(int));
		memset(this->base_array, 0, (SIZE+1)*sizeof(int));
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
		} else {
			this->value = val;
		}
		return -1;
	}
	int cell_destructor(cell_t* this) {
		free(this->list);
		free(this->base_array);
	}
	// WHAT WOULD HAVE BEEN A GORGEOUS cell_t struct :'(

	/*
		Serves the deep copy purpose when cell_t* definition was:
		struct cell_t {
			int row;
			int col;
			int* list;
			int* base_array;
			int n_allowed;
			int value;	
		}
	 */
	cell_t* deepcopy_board(cell_t* this) {
		cell_t* other = (cell_t*)malloc(BOARD_SIZE*sizeof(cell_t));
		// memset(other, 0, BOARD_SIZE*sizeof(cell_t));
		int i = 1;
		// print_board(this);
		for(; i<BOARD_SIZE; ++i) {
			// Copy the i'th
			other[i].row = this[i].row;
			other[i].col = this[i].col;
			other[i].n_allowed = this[i].n_allowed;
			other[i].value = this[i].value;
			other[i].list = (int*)malloc((SIZE+1)*sizeof(int));
			other[i].base_array = (int*)malloc((SIZE+1)*sizeof(int));
			// memcpy((void*)(other[i].list), (void*)(this[i].list), (SIZE+1)*sizeof(int));
			int j;
			for(j=0; j<=SIZE; j++) {
				(other[i].list)[j] = (this[i].list)[j];
				(other[i].base_array)[j] = (this[i].base_array)[j];
			}
			// memcpy((void*)(other[i].base_array), (void*)(this[i].base_array), (SIZE+1)*sizeof(int));
			
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
	cell_t* pushf(pool_t* this, cell_t* new_board) { // ASSERT: cell_t* is a fixed size board etc.
		struct pool_node_t* new_node = (pool_node_t*) malloc(sizeof(pool_node_t));
		// new_node->board_config = (cell_t*)malloc((BOARD_SIZE)*sizeof(cell_t));
		// memcpy((void*)new_node->board_config, (void*)new_board, (BOARD_SIZE)*sizeof(cell_t));
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
		// new_node->board_config = (cell_t*)malloc((BOARD_SIZE)*sizeof(cell_t));
		// memcpy((void*)new_node->board_config, (void*)new_board, (BOARD_SIZE)*sizeof(cell_t));
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

			// Iterates initialization_helper ... values between 2 and SIZE... Note: 0,1 are also included, but they're likely to be empty.
		int idx_hma = 0; 
			// Iterates initialization_helper[hma] ... values between 0 and initialization_helper_sizes[hma]
		int idx_cell_list;
		for(hma=0; hma<=SIZE; hma++) {
			if( mypushes >= max_pushes ) break;
			for( idx_hma=0; idx_hma< initialization_helper_sizes[hma]; idx_hma++) {
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
					set_value(&((newboard)[branch_on]), board[branch_on].list[idx_cell_list]);
					print_board(newboard);
				}
			}
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
		for(r=0; r<SIZE; r++) {
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

	// Initialize pool_t
	int error = PARTIAL_SOLN;
	error = heuristic_solve(board); // Will become a shared variable soon. TODO : Uncomment.
	if (error == SOLVED || error == NO_SOLN) {
		return cell2board(board);
	} // Simplify the problem... without parallelization man.


	pool_t* workpool = (pool_t*) malloc(sizeof(pool_t)); // Really just three integers.
	workpool->hd = NULL; 
	workpool->tl = NULL; 
	workpool->sz = 0;


	branch( workpool, board, INITIAL_POOL_SIZE );
	
	int solution_found;
	solution_found = 0;
	cell_t* solved_board;
	#pragma omp parallel shared(solution_found) private(error)
	{
		cell_t* thread_board;
		while(!solution_found) {
			#pragma omp critical
			{
				// WHY CAN'T I HAVE ALIGNED BRACKETS?!
				thread_board = (cell_t*)popf(workpool); // TODO: Convert to POP
			}
			if(solution_found) break;
			error = heuristic_solve(thread_board); // error is a private variable;
			if(error == SOLVED) {
				// #pragma omp atomic
				#pragma omp critical
					solved_board = thread_board;
				solution_found = 1; // Doesn't need to be critical. it is set to 1 afterall.
				break;
			}
			if(branch(workpool, thread_board, 5) == -1) { // This function is thread_safe inside the 
				error = dfs(thread_board);
			}
			if( error == SOLVED ) {
				#pragma omp critical
					solved_board = thread_board;
				solution_found = 1; // Doesn't need to be critical. it is set to 1 afterall.
				break;
			}
			free(thread_board);
			thread_board = NULL;		
		}
		// Got here only by breaking or by solution_found being true...
		if(error == SOLVED) { // This thread found the solution.
			// solved_board must've been set.
			#pragma omp critical
			{
				solved_board = deepcopy_board(thread_board);
			}
		} else { // my thread_board is NOT a solution.
			free(thread_board);
		}
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
	Can we do better ?? ... yes... look at brute_force.c
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

	int lst_indices[BOARD_SIZE];
	memset(lst_indices,0, BOARD_SIZE*sizeof(int));
	int lst_idx;
	// lst_idx = 0; ... initialized inside forloop.
	for(i = 1; i <= num_allowed;i++){
		int val = board[min_idx].list[i];
		set_value(&board[min_idx],val);
		int jdx = 0;
		lst_idx = 0;
		for(;jdx < SIZE;jdx++){
			int jdx_r = r*SIZE + jdx + 1;
		  	int jdx_c = jdx*SIZE + c + 1;
		  	int jdx_b = ((r_prime + (jdx / MINIGRIDSIZE))*SIZE) + (c_prime + (jdx%MINIGRIDSIZE) ) + 1;
		  	if( is_allowed(&board[jdx_r],val) ) {
		  		lst_indices[lst_idx] = jdx_r;
		  		lst_idx++;
		  		remove_allowed(&board[jdx_r],val);
		  	}
		  	if( is_allowed(&board[jdx_c],val) ) {
		  		lst_indices[lst_idx] = jdx_c;
		  		lst_idx++;
		  		remove_allowed(&board[jdx_c],val);
		  	}
		  	if( is_allowed(&board[jdx_b],val) ) {
		  		lst_indices[lst_idx] = jdx_b;
		  		lst_idx++;
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
		// jdx = 0; // Need to reset iterator, bro.

		// // WRONG! ... you add them only if they were originally allowed.
		// for(;jdx < SIZE;jdx++){
		// 	int jdx_r = r*SIZE + jdx + 1;
		//   	int jdx_c = jdx*SIZE + c + 1;
		//   	int jdx_b = ((r_prime + (jdx / MINIGRIDSIZE))*SIZE) + (c_prime + (jdx%MINIGRIDSIZE) ) + 1;
		//   	add_allowed(&board[jdx_r],val);
		//   	add_allowed(&board[jdx_c],val);
		//   	add_allowed(&board[jdx_b],val);
		// }
	}
	board[min_idx].value = 0;
	board[min_idx].n_allowed = num_allowed;
	for(i = 1; i <= num_allowed;i++){
		board[min_idx].base_array[i] = base_array[i];
	}
	return NO_SOLN;
}

