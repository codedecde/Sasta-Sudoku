int hs(int* brd,int* nals,long* bits, int* base){
	int fl_changed = 1, fl_unfilled = 0;
	int idx,val;
	int type;
	int b,r,c,bidx,count,vidx,rp,cp;

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
