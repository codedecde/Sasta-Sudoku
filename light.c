/*2013CS10225_2013CS10773*/
#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include <assert.h>
#include <limits.h>
#include <omp.h>	// OpenMP
#include "sudoku.h"

#define BOARD_SIZE SIZE*SIZE

int hs(int*, long*, int*);
int hdfs(int*, long*, int*);

void printboard(int* brd) {
	int r,c;
	for(r=0; r<SIZE; ++r) {
		for(c=0; c<SIZE; ++c) {
			printf("%d ", brd[r*SIZE + c]);
		}
		printf("\n");
	}
}

#define TO_IDX(r,c) (r*SIZE + c)
int** solveSudoku(int** origmat) {
	int* brd;
	long* bits;
	int* nals;

	// Initialize board, bitmasks, nals.
	brd = malloc(BOARD_SIZE*sizeof(int));
	int i;
	for(i=0; i<SIZE; ++i) {
		memmove(brd+ i*SIZE, origmat[i], SIZE*sizeof(int));
	}
	
	// ASSERT: initialization done.
	return origmat;
}