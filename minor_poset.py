#!/usr/bin/env python3
import sys
from cdIndex import *

def inputError(input, charIndex):
	print("There was an input error at character ",charIndex," that character is a ",input[charIndex])
	print("Now exiting, running with \"debug\" as an argument will print the input")
	if 'debug' in sys.argv: print(input)
	exit()

help_string=''.join([
	"Usage: cdIndexCalculator -m minor_poset INPUT",
	"\tcdIndexCalculator -m minor_poset -f FILE",
	"\n\n",
	"The input for this module specifies a lattice by for each element of the lattice, except the zerohat, listing the elements it covers.",
	"These cover lists are comma delimited, and the lists",
	"themselves are separated by semicolons.",
	"For example, B_3 is encoded as 0;0;0;1,2;1,3;2,3;4,5,6;",
	"when labeled as below:",
	"""
	         7
	       / | \
	      4  5  6
	     | \/ \/ |
	     | /\ /\ |
	      1  2  3
	       \ | /
	         0
	"""])
########################
#input specified in file
########################
if '-f' in sys.argv[:-1]:
	try:
		f = open(sys.argv[sys.argv.index('-f')+1],'rb')
		input=[]
		for line in f:
			input.append(line[:line.find(b'#')])
		input=b''.join(input).decode('ascii')
		f.close()
	except:
		print("Error reading file ",sys.argv[sys.argv.index('-f')+1])
		exit()
################################
#input specified on command line
################################
else:
	#first argument after removing the two specifying this module
	input=[sys.argv[i] for i in range(0,len(sys.argv)) if i!= sys.argv.index('-m') and i-1!=sys.argv.index('-m')][1]

##################################################
#Extra generators besides join irreducibles to add
##################################################
if '-g' in sys.argv[:-1]:
	try:
		gens_to_add=[int(i) for i in sys.argv[sys.argv.index('-g')+1].split(',')]
	except:
		print(help_string)
		exit()
else:
	gens_to_add=[]

#################################################
#format input as given into a list of cover lists
#################################################
input=[[]]+[[int(y) for y in x.split(',')] for x in ''.join([c for c in input if c in ';,0123456789']).split(';')[:-1]]


debugPrint('input',input)

irr = [i for i in range(0,len(input)) if len(input[i])==1]
L=[[1 if i in input[j] else 0 for j in range(0,len(input))]for i in range(0,len(input))] #cover matrix for L
for i in range(0,len(L)): L[i][i]=0

transClose(L)

genL = sorted(list(set(irr + gens_to_add)))
##################################


#compute a table of all joins in L
##################################
joins = [[0 for i in range(0,len(L))]for j in range(0,len(L))]
for i in range(0,len(L)):
	for j in range(i,len(L)):
		k = join(i,j,L)
		if k == None:
			print('join of ',end='')
			print(i,end='')
			print(' and ',end='')
			print(j,end='')
			print(' does not exist.',end='')
			print('input not a lattice. Exiting')
			exit()
		joins[i][j]=k
		joins[j][i]=k

######################################################
#compute all the minors
#minors are encoded as a list whose first element
#is the minimal element and the second element is the
#list of generators
######################################################
minors = [[0,genL]]
minors_M = [[0]]
minors_ranks = [[] for i in range(0,len(genL)+1)]
minors_ranks[len(genL)].append(1) #will add a zerohat later at index 0
new = [[0,genL]]
while len(new)>0:
	old = new
	new = []
	for l in old:
		r = minors.index(l)
		for i in range(0,len(l[1])):
			minor=[l[0],l[1][:i]+l[1][i+1:]] #delete i
			if minor in minors:
				s = minors.index(minor) #save index to adjust incidence matrix
			else:
				#add the minor and add a row and column to the incidence matrix
				s = len(minors_M)
				minors_ranks[len(minor[1])].append(s+1)
				minors.append(minor)
				for x in minors_M: x.append(0)
				minors_M.append([0 for x in range(-1,s)])
			minors_M[r][s] = -1
			minors_M[s][r] = 1
			if minor not in new: new.append(minor)

			#contract i
			temp = set([joins[l[1][i]][j] for j in l[1]])
			temp.remove(l[1][i])
			minor=[l[1][i],sorted(list(temp))]
			if minor in minors:
				s = minors.index(minor)
			else:
				s = len(minors_M)
				minors_ranks[len(minor[1])].append(s+1)
				minors.append(minor)
				for x in minors_M: x.append(0)
				minors_M.append([0 for x in range(-1,s)])
			minors_M[r][s] = -1
			minors_M[s][r] = 1
			if minor not in new: new.append(minor)

#add zerohat element to minor poset
minors_ranks = [[0]]+minors_ranks
for i in range(0,len(minors_M)):
	minors_M[i] = [-1]+minors_M[i]
minors_M = [[0]+[1 for i in range(0,len(minors_M))]]+minors_M

if 'debug' in sys.argv:
	print('genL: ',end='')
	print(genL)
	print('joins: ',end='')
	print(joins)
	print('minors: ',end='')
	print(minors)
	print('\n\nminors_ranks: ',end='')
	print(minors_ranks)
	print('\n\nminors_M: ',end='')
	print(minors_M)

transClose(minors_M)
M=minors_M
ranks=minors_ranks
