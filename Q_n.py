import sys
help_string="usage: cdIndexCalculator -m Q_n n\n\nThe provided argument n should be a positive integer greater than 1."

try:
	n=int((sys.argv[:sys.argv.index('-m')]+sys.argv[sys.argv.index('-m')+2:])[1])
except:
	print(help_string)
	exit()

######################################################
#faces of the cube
#are stored as numbers 0,(3**n)-1
#order is the termiwse order on the ternary digits
#induced by the order 0<2>1 with 0 and 1 incomparable
#####################################################
def less(i,j):
	if i==j: return False
#	print("less(",i,",",j,")")
	ti=i%3 #current ternary digit of i and j
	tj=j%3
	if (ti>tj) or (ti==0 and tj==1): return False
	power=1
#	print("\tti=",ti," tj=",tj)
	for _ in range(1,n):
		power*=3
#		print("\tti=",ti," tj=",tj)
		ti=(i%(3*power)-ti)//power
		tj=(j%(3*power)-tj)//power

		if (ti>tj) or (ti==0 and tj==1): return False
	return True

def rk(i): #rank is number of ternary digits equal to 2
	ti=i%3
	ret=ti//2
	power=1
	for _ in range(1,n):
		power*=3
		ti=(i%(3*power)-ti)//power

		ret += ti//2
	return ret

M=[[0 for j in range(0,3**n)]for i in range(0,3**n)] #M[i][j] == 1 if i<j, == -1 if i>j, ==0 otherwise
#add in zerohat
for r in M: r.append(-1)
M.append([1 for i in range(0,len(M))]+[0])

#make ranks
for i in range(0,len(M)):
	for j in range(i+1,len(M)):
		if less(i,j):
			M[i][j]=1
			M[j][i]=-1
		if less(j,i):
			M[i][j]=-1
			M[j][i]=1
ranks=[[]for i in range(0,n+1)] #ranks[i] is the list of all rank i elements
for i in range(0,3**n):
	ranks[rk(i)].append(i)
ranks=[[len(M)-1]]+ranks #add zerohat
#for i in range(0,3**n):
#	for j in range(i,3**n):
#		print(i," ","<" if less(i,j) else ">" if less(j,i) else "X"," ",j)
#
#exit()
