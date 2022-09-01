import sys
help_string="usage: cdIndexCalculator -m B_n n\n\nThe provided argument n should be a positive integer greater than 1."

try:
	n=int((sys.argv[:sys.argv.index('-m')]+sys.argv[sys.argv.index('-m')+2:])[1])
except Exception as e:
	print(e)
	print(help_string)
	exit()
M=[[1 if (i!=j and i&j==i) else (-1 if (i!=j and i&j==j) else 0) for j in range(0,1<<n)]for i in range(0,1<<n)] #M[i][j] == 1 if i<j, == -1 if i>j, ==0 otherwise
ranks=[[]for i in range(0,n+1)] #ranks[i] is the list of all rank i elements
for i in range(0,1<<n):
	ranks[len([c for c in bin(i)[2:] if c=='1'])].append(i)
