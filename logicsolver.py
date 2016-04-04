
f = open('tinysudoku.txt')
grid = [[0 for idx in xrange(9)] for jdx in xrange(9)]
ii = 0
for line in f:
	line = line.strip().split()
	jj = 0
	for elem in line:
		grid[ii][jj] = eval(elem)
		jj+=1
	ii+=1

allowed_set = [[set() for idx in xrange(9)] for jdx in xrange(9)]
for r in xrange(9):
	for c in xrange(9):
		if(grid[r][c] != 0):
			continue
		else:
			for idx in xrange(1,10):
				allowed_set[r][c].add(idx)
			r_prime = r - (r % 3)
			c_prime = c - (c % 3)
			for idx in xrange(9):
				if(grid[r][idx] in allowed_set[r][c]):
					allowed_set[r][c].remove(grid[r][idx])
				if(grid[idx][c] in allowed_set[r][c]):
					allowed_set[r][c].remove(grid[idx][c])
				if(grid[r_prime + (idx / 3)][c_prime + (idx%3) ] in allowed_set[r][c]):
					allowed_set[r][c].remove(grid[r_prime + (idx / 3)][c_prime + (idx%3) ])

flag = True
while(flag):
	flag = False
	for r in xrange(9):
		for c in xrange(9):
			if(len(allowed_set[r][c]) == 1):
				flag = True
				val = allowed_set[r][c].pop()
				grid[r][c] = val
				r_prime = r - (r % 3)
				c_prime = c - (c % 3)
				for idx in xrange(9):
					if(val in allowed_set[r][idx]):
						allowed_set[r][idx].remove(val)
					if(val in allowed_set[idx][c]):
						allowed_set[idx][c].remove(val)
					if(val in allowed_set[r_prime + (idx / 3)][c_prime + (idx%3) ]):
						allowed_set[r_prime + (idx / 3)][c_prime + (idx%3) ].remove(val)
for idx in xrange(9):
	for jdx in xrange(9):
		print grid[idx][jdx],
	print ''