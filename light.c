/*2013CS10225_2013CS10773*/
#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include <assert.h>
#include <limits.h>
#include <omp.h>	// OpenMP
#include "sudoku.h"

#define BOARD_SIZE SIZE*SIZE
#define NO_SOLN -1
#define PARTIAL_SOLN 0
#define SOLVED 1
#define TIME_TO_LEAVE 42


int hs(int*,int*, long* );
int hdfs(int**,int**,long**);
int set_board_value(int*, int*, long*, int*,int,int);
void undo_alterations(int*, int*, long*, int*, int,int);


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
	nals = malloc(BOARD_SIZE*sizeof(int));
	bits = malloc(BOARD_SIZE*sizeof(long));
	int i;
	for(i=0; i<SIZE; ++i) {
		memmove(brd+ i*SIZE, origmat[i], SIZE*sizeof(int));
	}
	int r,c,b,idx,rp,cp;
	for(i=0 ;i < BOARD_SIZE;++i){
		r = i/SIZE;
		c = i%SIZE;
		if(origmat[r][c] != 0){
			bits[i] = 0;
			nals[i] = 0;
		}
		else{
			bits[i] = (1<<(SIZE+1)) - 2; //Everything allowed
			nals[i] = SIZE;

			rp = r - (r%MINIGRIDSIZE);
			cp = c - (c%MINIGRIDSIZE);
			for(idx = 0; idx< SIZE;++idx){
				b =origmat[r][idx] ; 
				if(b != 0){
					if( bits[i] & 1<<b ) {
						nals[i]--;
						bits[i] &= ~(1<<b);
					}
				}
				b =origmat[idx][c] ; 
				if(b != 0){
					if( bits[i] & 1<<b ) {
						nals[i]--;
						bits[i] &= ~(1<<b);
					}
				}
				b = origmat[rp + (idx / MINIGRIDSIZE)][cp + (idx%MINIGRIDSIZE)];
				if(b != 0){
					if( bits[i] & 1<<b ) {
						nals[i]--;
						bits[i] &= ~(1<<b);
					}
				}
			}
		}
	}

	int global_solved = PARTIAL_SOLN;
	// ASSERT: initialization done.
	global_solved = hdfs(&brd,&nals,&bits);
	for(r = 0; r<SIZE;++r){
		for(c = 0; c < SIZE; ++c){
			origmat[r][c] = brd[TO_IDX(r,c)];
		}
	}
	return origmat;
}

int hdfs(int** ptrbrd,int** ptrnals,long** ptrbits){
	int* brd = *ptrbrd;
	int* nals = *ptrnals;
	long* bits = *ptrbits;
	int mysolved = PARTIAL_SOLN;
	int* lbrd = malloc(BOARD_SIZE*sizeof(int));
	int* lnals = malloc(BOARD_SIZE*sizeof(int));
	long* lbits = malloc(BOARD_SIZE*sizeof(long));
	
	memmove(lbrd,brd,BOARD_SIZE*sizeof(int));
	memmove(lnals,nals,BOARD_SIZE*sizeof(int));
	memmove(lbits,bits,BOARD_SIZE*sizeof(long));

	/*TODO HS*/

	int min_idx = -1;
	int min_val = 2*SIZE;
	int loc_nal;
	int idx;
	for(idx = 0;idx<BOARD_SIZE;++idx){
		if(brd[idx] == 0 ){
			loc_nal = nals[idx];
			if(loc_nal == 0){
				free(lbrd);
				free(lnals);
				free(lbits);
				return NO_SOLN;
			}
			else if(loc_nal < min_val){
				min_val = loc_nal;
				min_idx = idx;
			}
		}
	}
	if(min_idx == -1){
		free(lbrd);
		free(lnals);
		free(lbits);
		return SOLVED;
	}

	/*Branch on min_idx */
	int val_iter;
	int* alters = malloc(4*SIZE*sizeof(int));
	int nalters = -1;
	for(val_iter = 1; val_iter <= SIZE; ++val_iter){
		if(bits[min_idx] & (1<<val_iter)){
			nalters = set_board_value(brd,nals,bits,alters,min_idx,val_iter);
			mysolved = hdfs(ptrbrd,ptrnals,ptrbits);
			if(mysolved == SOLVED){
				free(lbrd);
				free(lnals);
				free(lbits);
				return SOLVED;
			}else if(mysolved == TIME_TO_LEAVE){
				//TODO: Free memory and leave
			}
			undo_alterations(brd,nals,bits,alters,nalters,val_iter);
		}
	}
	/* NO SOLUTION*/
	free(brd);
	free(nals);
	free(bits);
	brd = lbrd;
	nals = lnals;
	bits = lbits;
	return NO_SOLN;
}

int set_board_value(int* brd, int* nals, long* bits, int* alters,int idx,int val ){
	int nalters = -1;
	int r,c,i,rp,cp;
	r = idx / SIZE;
	c = idx % SIZE;
	rp = r - (r%MINIGRIDSIZE);
	cp = c - (c%MINIGRIDSIZE);
	int b,bitval = 1<<val;
	brd[idx] = val;
	nals[idx] = 0;
	bits[idx] = 0;


	for(i = 0; i < SIZE ;++i ){
		b = TO_IDX(r,i);
		if(bits[b] & bitval){
			alters[++nalters] = b;
			bits[b] &= ~bitval;
			nals[b]--;
		}
		b = TO_IDX(i,c);

		if(bits[b] & bitval){
			alters[++nalters] = b;
			bits[b] &= ~bitval;
			nals[b]--;
		}
		
		b = TO_IDX(rp + (i / MINIGRIDSIZE),cp + (i%MINIGRIDSIZE));

		if(bits[b] & bitval){
			alters[++nalters] = b;
			bits[b] &= ~bitval;
			nals[b]--;
		}
	}

	return nalters;
}

void undo_alterations(int* brd, int* nals, long* bits, int* alters, int nalters,int val){
	int bitval = 1<<val;
	while(nalters > -1){
		bits[alters[nalters]] |= bitval;
		nals[alters[nalters]]--;
		nalters--;
	}
	return;
}
