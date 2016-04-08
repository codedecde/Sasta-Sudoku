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
#define NOT_MY_SOLN 13
#define TIME_TO_LEAVE 42


int hs(int*,int*, long* , int*);
int hdfs(int*,int*,long*, int*);
int set_board_value(int*, int*, long*, int*,int,int);
void undo_alterations(int*, int*, long*, int*, int,int);

int mypow(int base, int power) { // Because C confused me .
	int i;
	int acc;
	acc = 1;
	for(i=0; i<power; ++i) {
		acc = acc*base;
	}
	return acc;
} 

/*
	Global board:
		int* brd;
		int* nals;
		long* bits;
	Workpool
		int wpl_sz, wpl_szalloc
		int** wpl_brd
		int** wpl_nals
		long** wpl_bits

	Global Graveyard ... declared inside thread
		int grvy_sz, grvy_szalloc
		int** grvy_brd
		int** grvy_nals
		long** grvy_bits
*/
// Global state...
int* brd;
long* bits;
int* nals;
int global_solved;
int wpl_sz, wpl_szalloc;
int** wpl_brd;
int** wpl_nals;
long** wpl_bits;

int grvy_sz, grvy_szalloc;
int** grvy_brd;
int** grvy_nals;
long** grvy_bits;

void initwpl(int initsz) {
	// Alloc.
	wpl_brd = malloc(initsz*sizeof(int*));
	wpl_nals =  malloc(initsz*sizeof(int*));
	wpl_bits =  malloc(initsz*sizeof(long*));
	wpl_szalloc = initsz;
	wpl_sz = 0;
}

int pushwpl(int* brd, int* nals, long* bits, int** brdcpy, int** nalscpy, long** bitscpy) {
	if(wpl_sz >= wpl_szalloc) {
		// Grows array.
		int** oldbrds = wpl_brd;
		wpl_brd = malloc(2*wpl_szalloc*sizeof(int*));
		memcpy(wpl_brd, oldbrds, wpl_szalloc*sizeof(int*));
		free(oldbrds);

		int** oldnals = wpl_nals;
		wpl_nals = malloc(2*wpl_szalloc*sizeof(int*));
		memcpy(wpl_nals, oldnals, wpl_szalloc);
		free(oldnals);

		long** oldbits = wpl_bits;
		wpl_bits  = malloc(2*wpl_szalloc*sizeof(long*));
		memcpy(wpl_bits, oldbits, wpl_szalloc);
		free(oldbits);
		wpl_szalloc *= 2;
		return 0;
	}
	wpl_brd[wpl_sz] = malloc(BOARD_SIZE*sizeof(int));
	wpl_nals[wpl_sz] = malloc(BOARD_SIZE*sizeof(int));
	wpl_bits[wpl_sz] = malloc(BOARD_SIZE*sizeof(long));
	memcpy(wpl_brd[wpl_sz], brd, BOARD_SIZE*sizeof(int));
	memcpy(wpl_nals[wpl_sz], nals, BOARD_SIZE*sizeof(int));
	memcpy(wpl_bits[wpl_sz], bits, BOARD_SIZE*sizeof(long));
	*brdcpy = wpl_brd[wpl_sz];
	*nalscpy = wpl_nals[wpl_sz];
	*bitscpy = wpl_bits[wpl_sz];
	wpl_sz++;
	return 1;
}

int popwpl(int** brdcpy, int** nalscpy, long** bitscpy) {
	if(wpl_sz>0) {
		wpl_sz--;
		*brdcpy = wpl_brd[wpl_sz];
		*nalscpy = wpl_nals[wpl_sz];
		*bitscpy = wpl_bits[wpl_sz];
		return 1;
	}
	return 0;
}

void destroy_wpl() {
	int i;
	for(i=0; i<wpl_sz; ++i) {
		free(wpl_brd[i]);
		free(wpl_nals[i]);
		free(wpl_bits[i]);
	}
	wpl_sz = 0;
	wpl_szalloc = 0;
}

void initgrvy(int isz) {
	grvy_brd  = malloc(isz*sizeof(int*));
	grvy_nals = malloc(isz*sizeof(int*));
	grvy_bits = malloc(isz*sizeof(long*));
	grvy_szalloc = isz;
	int i;
	for(i=0; i<isz; ++i) {
		grvy_brd[i] = malloc(BOARD_SIZE*sizeof(int));
		grvy_nals[i]= malloc(BOARD_SIZE*sizeof(int));
		grvy_bits[i]= malloc(BOARD_SIZE*sizeof(long));
	}
	grvy_sz = isz;
}
void destroy_grvy() {
	int i;
	for(i=0; i<grvy_sz; ++i) {
		free(grvy_brd[i]);
		free(grvy_nals[i]);
		free(grvy_bits[i]);
	}
	grvy_sz = 0;
	grvy_szalloc = 0;
}

int pushmem(int* brd, int* nals, long* bits) {
	if(omp_in_parallel()!=0) {
	#pragma omp critical (graveyard_lock)
	{
		if(grvy_sz >= grvy_szalloc) {
			// Need to grow memory.
			int** oldbrds = grvy_brd;
			grvy_brd = malloc(2*grvy_szalloc*sizeof(int*));
			memcpy(grvy_brd, oldbrds, grvy_szalloc*sizeof(int*));
			free(oldbrds);

			int** oldnals = grvy_nals;
			grvy_nals = malloc(2*grvy_szalloc*sizeof(int*));
			memcpy(grvy_nals, oldnals, grvy_szalloc*sizeof(int*));
			free(oldnals);
		
			long** oldbits = grvy_bits;
			grvy_bits = malloc(2*grvy_szalloc*sizeof(long*));
			memcpy(grvy_bits, oldbits, grvy_szalloc*sizeof(long*));
			free(oldbits);

			grvy_szalloc*=2;
		}
		grvy_brd[grvy_sz] = brd;
		grvy_nals[grvy_sz]=nals;
		grvy_bits[grvy_sz]=bits;
		grvy_sz++;		
	}
	} else {
		if(grvy_sz >= grvy_szalloc) {
			// Need to grow memory.
			int** oldbrds = grvy_brd;
			grvy_brd = malloc(2*grvy_szalloc*sizeof(int*));
			memcpy(grvy_brd, oldbrds, grvy_szalloc*sizeof(int*));
			free(oldbrds);

			int** oldnals = grvy_nals;
			grvy_nals = malloc(2*grvy_szalloc*sizeof(int*));
			memcpy(grvy_nals, oldnals, grvy_szalloc*sizeof(int*));
			free(oldnals);
		
			long** oldbits = grvy_bits;
			grvy_bits = malloc(2*grvy_szalloc*sizeof(long*));
			memcpy(grvy_bits, oldbits, grvy_szalloc*sizeof(long*));
			free(oldbits);

			grvy_szalloc*=2;
		}
		grvy_brd[grvy_sz] = brd;
		grvy_nals[grvy_sz]=nals;
		grvy_bits[grvy_sz]=bits;
		grvy_sz++;
	}
}

int popmem(int** pbrd, int** pnals, long** pbits) {
	int ret;
	if(omp_in_parallel()!=0) {
	#pragma omp critical (graveyard_lock)
	{
		if(grvy_sz == 0) {
			// Need to allocate new memory. 
			*pbrd = malloc(BOARD_SIZE*sizeof(int*));
			*pnals= malloc(BOARD_SIZE*sizeof(int*));
			*pbits= malloc(BOARD_SIZE*sizeof(long*));
			ret = 0;
		} else {
			grvy_sz--;
			*pbrd = grvy_brd[grvy_sz];
			*pnals= grvy_nals[grvy_sz];
			*pbits= grvy_bits[grvy_sz];
			
			grvy_brd[grvy_sz] = NULL;
			grvy_nals[grvy_sz] = NULL;
			grvy_bits[grvy_sz] = NULL;
			ret = 1;
		}
	}
	} else {
		if(grvy_sz == 0) {
			// Need to allocate new memory. 
			*pbrd = malloc(BOARD_SIZE*sizeof(int*));
			*pnals= malloc(BOARD_SIZE*sizeof(int*));
			*pbits= malloc(BOARD_SIZE*sizeof(long*));
			ret = 0;
		} else {
			grvy_sz--;
			*pbrd = grvy_brd[grvy_sz];
			*pnals= grvy_nals[grvy_sz];
			*pbits= grvy_bits[grvy_sz];
			ret = 1;
		}
	}
	return ret;
}
void print_board(int* brd) {
	int r,c;
	for(r=0; r<SIZE; ++r) {
		for(c=0; c<SIZE; ++c) {
			printf("%d ", brd[r*SIZE + c]);
		}
		printf("\n");
	}
}

int getbon(int* brd, int* nals, long* bits) {
	int idx;
	int minv=2*SIZE; 
	int mini=-2;
	for(idx=0; idx<BOARD_SIZE; ++idx) {
		if(brd[idx]==0) {
			if(nals[idx]==0) {
				return NO_SOLN;
			} else if(minv>nals[idx]) {
				minv = nals[idx];
				mini = idx;
			}
		} 
	}
	return mini;
}

void populate_wpl(int ncl) {
	int lol = mypow(SIZE,ncl);
	int** stkbrd = malloc(lol*sizeof(int*)); 
	int** stknals = malloc(lol*sizeof(int*));
	long** stkbits = malloc(lol*sizeof(long*));
	
	int stksz = -1;
	int val,bon;
	int *brd2br, *nals2br, *brd2md, *nals2md, *brd2ps, *nals2ps;
	long *bits2br, *bits2md, *bits2ps;
	// pushwpl(brd, nals, bits, &nbrd, &nnals, &nbits); // lel.
	// pushmem(brd, nals, bits);
	pushwpl(brd, nals, bits, &brd2ps, &nals2ps, &bits2ps);

	for(lol=0; lol<ncl; ++lol) {
		while(wpl_sz > 0) {
			popwpl(&brd2ps, &nals2ps, &bits2ps);
			stksz++;
			stkbrd[stksz] = brd2ps;
			stknals[stksz] = nals2ps;
			stkbits[stksz] = bits2ps;
		}

		while(stksz > -1) {
			// pop from stack.
			brd2br = stkbrd[stksz];
			nals2br = stknals[stksz];
			bits2br = stkbits[stksz];
			stksz--;
			// branch
			bon = getbon(brd2br, nals2br, bits2br);
			// 		if can't branch, then move on.
			if( bon == -2 ) { // SOLVED?
				// FREE STUFF
				pushmem(brd2br, nals2br, bits2br);
			} else if(bon==-1) { // Invalid board.
				// FREE STUFF
				pushmem(brd2br, nals2br, bits2br);
			} else {
				// push on wpl
				for(val=1; val<=SIZE; ++val) {
					if( bits2br[bon] & (1<<val) ) {
						pushwpl(brd2br, nals2br, bits2br, &brd2md, &nals2md, &bits2md);
						set_board_value(brd2md, nals2md, bits2md, NULL, bon, val);
					}
				}
			}
		}
	}
}

#define TO_IDX(r,c) ((r)*SIZE + (c))
int** solveSudoku(int** origmat) {

	// brd = malloc(BOARD_SIZE*sizeof(int));
	// nals = malloc(BOARD_SIZE*sizeof(int));
	// bits = malloc(BOARD_SIZE*sizeof(long));

	// Initialize board, bitmasks, nals.
	global_solved = PARTIAL_SOLN;
	initgrvy(100); // TODO: change initial graveyard size.
	initwpl(100);
	popmem(&brd, &nals, &bits);

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

	global_solved = PARTIAL_SOLN;
	int* arr_lr = malloc(3*BOARD_SIZE*sizeof(int));
	// ASSERT: initialization done.

	global_solved = hs(brd,nals,bits, arr_lr);
	free(arr_lr);
	if(global_solved == SOLVED) {
		for(r = 0; r<SIZE;++r){ // Store brd into origmat.
			for(c = 0; c < SIZE; ++c){
				origmat[r][c] = brd[TO_IDX(r,c)];

			}
		}
	}

	// ASSERT: global_solved = PARTIAL_SOLN	
	
	populate_wpl(MINIGRIDSIZE); // Number of cells to fill up.
	// printf("wpl sz= %d\n", wpl_sz);
	#pragma omp parallel shared(wpl_sz,wpl_szalloc,wpl_bits, wpl_nals, wpl_brd, grvy_sz, grvy_szalloc, grvy_bits, grvy_nals, grvy_brd)
	{
		int* tbrd;
		int* tnals;
		long* tbits;
		int* arr_tlr;
		int didpop, loc_solved;
		arr_tlr = malloc(3*BOARD_SIZE*sizeof(int));
		loc_solved = PARTIAL_SOLN;
		while(wpl_sz>0) {
			#pragma omp critical (workpool_lock)
			{
				didpop = popwpl(&tbrd, &tnals, &tbits);
			}
			if(didpop != 1) {
				break; // Did not pop.
			}
			// printf("tid:%d got here\n", omp_get_thread_num());
			loc_solved = hdfs(tbrd, tnals, tbits, arr_tlr);
			if(loc_solved == SOLVED) {
				#pragma omp critical (global_solved_lock)
				{
					if(global_solved != SOLVED) {
						global_solved = SOLVED;
						brd = tbrd;
						nals = tnals;
						bits = tbits;
					} else {
						loc_solved = NOT_MY_SOLN;
					}
				}
				break;
			} else if(loc_solved==TIME_TO_LEAVE) {
				break;
			} else if(global_solved == SOLVED) {
				break;
			}
			#pragma omp critical (avoiding_dbf)
			{
				pushmem(tbrd, tnals, tbits);
				didpop = 0; // Prevents it from getting pushed again.
			}
		}
		free(arr_tlr);
		if( loc_solved == SOLVED || didpop != 1) {
			// Do nothing.
		} else {
			pushmem(tbrd, tnals, tbits);
		}
	}

	// Done execution.
	for(r = 0; r<SIZE;++r){
		for(c = 0; c < SIZE; ++c){
			origmat[r][c] = brd[TO_IDX(r,c)];

		}
	}
	destroy_wpl();
	destroy_grvy();
	return origmat;
}

int hdfs(int* brd,int* nals,long* bits, int* base){
	
	int mysolved = PARTIAL_SOLN;
	int* lbrd; 
	int* lnals; 
	long* lbits; 
	// lbrd = malloc(BOARD_SIZE*sizeof(int));
	// lnals = malloc(BOARD_SIZE*sizeof(int));
	// lbits = malloc(BOARD_SIZE*sizeof(long));
	popmem(&lbrd, &lnals, &lbits);

	memmove(lbrd,brd,BOARD_SIZE*sizeof(int));
	memmove(lnals,nals,BOARD_SIZE*sizeof(int));
	memmove(lbits,bits,BOARD_SIZE*sizeof(long));

	mysolved = hs(brd,nals,bits, base);

	if(mysolved == SOLVED || mysolved == TIME_TO_LEAVE){
		// free(lbrd);
		// free(lnals);
		// free(lbits);
		pushmem(lbrd, lnals, lbits);
		return mysolved;
	}
	if(mysolved == NO_SOLN){
		memmove(brd, lbrd, BOARD_SIZE*sizeof(int));
		memmove(nals, lnals, BOARD_SIZE*sizeof(int));
		memmove(bits, lbits, BOARD_SIZE*sizeof(long));
		// free(lbrd);
		// free(lnals);
		// free(lbits);
		pushmem(lbrd, lnals, lbits);
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
				pushmem(lbrd, lnals, lbits);
				return NO_SOLN;
			} else if(loc_nal < min_val){
				min_val = loc_nal;
				min_idx = idx;
			}
		}
	}
	
	if(min_idx == -1){
		pushmem(lbrd, lnals, lbits);
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
				pushmem(lbrd, lnals, lbits);
				return SOLVED;
			}else if(global_solved == SOLVED || mysolved == TIME_TO_LEAVE){
				pushmem(lbrd, lnals, lbits);
				return TIME_TO_LEAVE;
			}
			undo_alterations(brd,nals,bits,alters,nalters,val_iter);
		}
	}
	
	/* NO SOLUTION*/
	memmove(brd, lbrd, BOARD_SIZE*sizeof(int));
	memmove(nals, lnals, BOARD_SIZE*sizeof(int));
	memmove(bits, lbits, BOARD_SIZE*sizeof(long));
	pushmem(lbrd, lnals, lbits);
	
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
	int idx,val,ctr;
	int type;
	int b,r,c,bidx,count,vidx,rp,cp;

	// Elimination
	while(fl_changed){
		fl_changed = 0;
		fl_unfilled = 0;
		for(idx = 0; idx < BOARD_SIZE;++idx){
			if(nals[idx] == 0){
				if(brd[idx] == 0){
					return NO_SOLN;
				}
			}else if(nals[idx] == 1){
				val = bits[idx];
				ctr = 0;
				while( 1 ) {
					if (val & 1) {
						break;
					} else {
						ctr++;
						val = val>>1;
					}
				}
				set_board_value(brd,nals,bits,NULL,idx,ctr); 
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
	
	while(fl_changed){
		fl_changed = 0;
		for(val = 1;val <= SIZE; ++val){
		  //Row iterator
		  
			for(r = 0; r < SIZE; ++r){
				count = 0;
				vidx = -1;
				for(c = 0; c < SIZE; ++c){
					idx = TO_IDX(r,c);
						if(brd[idx] == 0){
							if(nals[idx] == 0)
								return NO_SOLN;
							else if(brd[idx] && (1<<val)){
						 	 	count++;
						  	if(count == 2)
						    	break;
						  	else
						    	vidx = idx;
							}	
						}
				}
				if(count == 1){
				  fl_changed = 1;
				  set_board_value(brd,nals,bits,NULL,vidx,val);
				}
			}

			//Column Iterator


			for(c = 0; c < SIZE; ++c){
			count = 0;
			vidx = -1;
			for(r = 0; r < SIZE; ++r){
			  idx = TO_IDX(r,c);
			  if(brd[idx] == 0){
			    if(nals[idx] == 0)
			      return NO_SOLN;
			    else if(brd[idx] && (1<<val)){
			      count++;
			      if(count == 2)
			        break;
			      else
			        vidx = idx;
			    }
			  }
			}
			if(count == 1){
			  fl_changed = 1;
			  set_board_value(brd,nals,bits,NULL,vidx,val);
			}
			}

			//BOX ITERATOR
			for(r = 0; r < SIZE; ++r){
			vidx = -1;
			count = 0;
			rp = (r/MINIGRIDSIZE) * MINIGRIDSIZE;
			cp = (r%MINIGRIDSIZE) * MINIGRIDSIZE;
			for(c = 0; c < SIZE; ++c){
			  idx = TO_IDX((rp + (c/MINIGRIDSIZE)),(cp + (c%MINIGRIDSIZE)));
			  if(brd[idx] == 0){
			    if(nals[idx] == 0)
			      return NO_SOLN;
			    else if(brd[idx] && (1<<val)){
			      count++;
			      if(count == 2)
			        break;
			      else
			        vidx = idx;
			    }
			  }
			}
			if(count == 1){
			  fl_changed = 1;
			  set_board_value(brd,nals,bits,NULL,vidx,val);
			}
			}
		}
	}

	val = SOLVED;
	for(idx = 0 ; idx < BOARD_SIZE; ++idx){
		if(brd[idx]==0 ){
			if(nals[idx]==0){
				return NO_SOLN;
			}
			else{
				val = PARTIAL_SOLN;
			}
		}
	}
	return val;

}
