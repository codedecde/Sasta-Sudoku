/*
	Solves the board using brute force. 
	returns :
 */
#define SOLUTION 0
#define PARTIAL_SOLN 1
#define NO_SOLN 2	
int brute_force(cell_t* board) {
	// Construct list of spots to fill up.
	int lst_unfilled_idxs[BOARD_SIZE]; // Unnecessarily large, but that's okay.
	int iter_lui; // iterator for lst_unfilled_indices
	int n_unfilled;

	iter_lui = 0;
	
	int iter_board;
	for(iter_board = 1; iter_board< BOARD_SIZE; ++iter_board) { // The 0th position is always null.
		if( board[iter_board].value == 0 ) {
			lst_unfilled_idxs[iter_lui] = iter_board;
			++iter_lui;
		}
	}
	n_unfilled = iter_lui;

	// Now, to do iteration.
	// To do this, maintain a lst_(cell_node_t*)
	cell_node_t* lst_cellnodes[BOARD_SIZE];

	int fl_bfcont = 1;
	int idx_lst_idxs = 0;
	/*
		In any iteration, 
	 */
	while (fl_bfcont) {

	}

}