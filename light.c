/*2013CS10225_2013CS10773*/
#include <stdio.h>	// for god knows what.
#include <stdlib.h> // for malloc?
#include <string.h> // for memcpy?
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <omp.h>	// OpenMP
#include "sudoku.h"

#define BOARD_SIZE SIZE*SIZE
#define NO_SOLN -1
#define PARTIAL_SOLN 0
#define SOLVED 1
#define TIME_TO_LEAVE 42


int hs(int*,int*, long* , int*);
int hdfs(int*,int*,long*, int*);
int set_board_value(int*, int*, long*, int*,int,int);
void undo_alterations(int*, int*, long*, int*, int,int);


void print_board(int* brd) {
	int r,c;
	for(r=0; r<SIZE; ++r) {
		for(c=0; c<SIZE; ++c) {
			printf("%d ", brd[r*SIZE + c]);
		}
		printf("\n");
	}
}

#define TO_IDX(r,c) ((r)*SIZE + (c))
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
		memmove(brd + i*SIZE, origmat[i], SIZE*sizeof(int));
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
	int* arr_lr = malloc(3*BOARD_SIZE);
	// ASSERT: initialization done.

	global_solved = hdfs(brd,nals,bits, arr_lr);
	free(arr_lr);
	for(r = 0; r<SIZE;++r){
		for(c = 0; c < SIZE; ++c){
			origmat[r][c] = brd[TO_IDX(r,c)];

		}
	}
	return origmat;
}

int hdfs(int* brd,int* nals,long* bits, int* base){
	
	int mysolved = PARTIAL_SOLN;
	int* lbrd = malloc(BOARD_SIZE*sizeof(int));
	int* lnals = malloc(BOARD_SIZE*sizeof(int));
	long* lbits = malloc(BOARD_SIZE*sizeof(long));
	
	memmove(lbrd,brd,BOARD_SIZE*sizeof(int));
	memmove(lnals,nals,BOARD_SIZE*sizeof(int));
	memmove(lbits,bits,BOARD_SIZE*sizeof(long));

	mysolved = hs(brd,nals,bits, base);

	if(mysolved == SOLVED || mysolved == TIME_TO_LEAVE){
		free(lbrd);
		free(lnals);
		free(lbits);
		return mysolved;
	}
	if(mysolved == NO_SOLN){
		memmove(brd, lbrd, BOARD_SIZE*sizeof(int));
		memmove(nals, lnals, BOARD_SIZE*sizeof(int));
		memmove(bits, lbits, BOARD_SIZE*sizeof(long));
		free(lbrd);
		free(lnals);
		free(lbits);
		return NO_SOLN;
	}

	

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
	
	// printf("mindx,v:%d,%d\n", min_idx, min_val);
	// print_board(brd);

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

	// long oldbits = bits[min_idx];
	// int oldnals = nals[min_idx];
	long oldbits = bits[min_idx];
	
	for(val_iter = 1; val_iter <= SIZE; ++val_iter){
		if(oldbits & (1<<val_iter)){
			// printf("at minidx=%d, trying val=%d\n", min_idx, val_iter);
			nalters = set_board_value(brd,nals,bits,alters,min_idx,val_iter);
			mysolved = hdfs(brd,nals,bits, base);
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
	memmove(brd, lbrd, BOARD_SIZE*sizeof(int));
	memmove(nals, lnals, BOARD_SIZE*sizeof(int));
	memmove(bits, lbits, BOARD_SIZE*sizeof(long));
	free(lbrd);
	free(lnals);
	free(lbits);
	
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
	
	if(val < 0){
		printf("Fucked\n");
		assert(0);

	}

	for(i = 0; i < SIZE ;++i ){
		b = TO_IDX(r,i);
		if((bits[b] & bitval) > 0){
			if(alters != NULL) {
				alters[++nalters] = b;
			}
				
			bits[b] &= ~bitval;
			nals[b]--;
		}
		b = TO_IDX(i,c);

		if((bits[b] & bitval) > 0){
			if(alters != NULL) {
				alters[++nalters] = b;
			}
			bits[b] &= ~bitval;
			nals[b]--;
		}
		
		b = TO_IDX((rp + (i / MINIGRIDSIZE)),(cp + (i%MINIGRIDSIZE)));
		if((bits[b] & bitval) > 0){
			if(alters != NULL) {
				alters[++nalters] = b;
			}
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
		nals[alters[nalters]]++;
		nalters--;
	}
	return;
}


int hs(int* brd,int* nals,long* bits, int* base){
	int fl_changed = 1, fl_unfilled = 0;
	int idx,val;
	int type;
	int b,r,c,bidx;
	while(fl_changed){
			fl_changed = 0;
			fl_unfilled = 0;
			for(idx = 0; idx < BOARD_SIZE;++idx){
				if(nals[idx] == 0){
					if(brd[idx] == 0){
						return NO_SOLN;
					}
				}else if(nals[idx] == 1){
					val = log2(bits[idx] & -bits[idx]);
					set_board_value(brd,nals,bits,NULL,idx,val); 
					fl_changed = 1; 
				}
				else{
					fl_unfilled = 1;
				}
			}
		}
	if(fl_unfilled == 0){
		return SOLVED;
	}
	else{
		return PARTIAL_SOLN;
	}
	/*while(fl_changed){
		// Lone ranger...
		fl_changed = 0;

		// Make base
		memset(base, 0, 3*BOARD_SIZE);
		
		for(idx=0; idx<BOARD_SIZE; ++idx) {
			for(val=1; val<=SIZE; ++val) {
				if (bits[idx] & (1<<val)) {
					assert(brd[idx] == 0); //TODO: REMOVE
					
					r = (idx)/SIZE;
					base[r*SIZE + (val-1)] = (base[r*SIZE + (val-1)] ==  0 ) ? (idx+1) : -1;
					
					r = (idx%SIZE);
					base[BOARD_SIZE + r*SIZE + (val-1)] = (base[BOARD_SIZE + r*SIZE + (val-1)] == 0 ) ? (idx + 1) : -1;
					
					r = ((idx/SIZE)/(MINIGRIDSIZE))*MINIGRIDSIZE + (r/MINIGRIDSIZE);
					
					base[2*BOARD_SIZE + r*SIZE + (val-1)] = (base[2*BOARD_SIZE+ r*SIZE + (val-1)] == 0 ) ? (idx + 1) : -1;
				}
			}
		}
		// Use base
		for(r=0; r<SIZE; ++r) {
			for(c=0; c<SIZE; ++c) {
				bidx = r*SIZE + c;
				b = base[bidx] - 1;
				
				if((b >= 0 ) && (brd[b] == 0)) {
					set_board_value(brd,nals,bits,NULL,b,c+1);
					fl_changed = 1;
				}
					
				
				bidx += BOARD_SIZE;
				b = base[bidx] - 1;
				if((b >= 0 ) && (brd[b] == 0)) {
					set_board_value(brd,nals,bits,NULL,b,c+1);
					fl_changed = 1;
				}
				
				bidx += BOARD_SIZE;
				b = base[bidx] - 1;
				if((b >= 0 ) && (brd[b] == 0)) {
					set_board_value(brd,nals,bits,NULL,b,c+1);
					fl_changed = 1;
				}
			}
		}
	}
	// Set unfilled flag.
	fl_unfilled = 0;
	for(idx=0; idx<BOARD_SIZE; ++idx) {
		if(brd[idx]==0) {
			if(nals[idx]>0) {
				fl_unfilled = 1;
				break;
			} else {
				return NO_SOLN;
			}
		}
	}
	print_board(brd);
	if(fl_unfilled==0) {
		return SOLVED;
	} else {
		return PARTIAL_SOLN;
	}*/
}
