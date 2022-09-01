#!/usr/bin/env python3
import sys
from cdIndex import *
import importlib

help_string="".join([
	"usage: cdIndexCalculator INPUT",
	"       cdIndexCalculator -f FILE",
	"       cdIndexCalculator -m MODULE",
	"\n\n",
	"In the first form the argument INPUT describes the poset to calculate the cd-index of,",
	" it is specified as a comma seperated list of the covering elements",
	" with these lists seperated by semicolons and the lists in each rank",
	" of the poset seperated by semicolons. See the examples distributed along",
	" with this program.",
	"\n\n",
	"The second form is used to specify the input in a file instead of passing on the command line",
	" as an argument."
	"\n\n",
	"In either of the first two forms all characters in the input except ;,0123456789",
	" are ignored. Any characters after # on the same line are ignored as well.",
	"\n\n",
	"In the third form the argument MODULE is the name of a module that when imported calculates the incidence matrix",
	" and ranks of the target poset. These two variables should be named M and ranks. See the included examples for details.",
	" The module should be named ending with \".py\" and the argument MODULE should not include the extension \".py\""
	"\n\n",
	"If an argument \"debug\" is provided the program will output more information",
	" that is useful for debugging.",
	])

if '-f' in sys.argv[:-1]:
	file=sys.argv[sys.argv.index('-f')+1]
else:
	file=None
if '-m' in sys.argv[:-1]:
	module=sys.argv[sys.argv.index('-m')+1]
	if module[-3:]=='.py': module=module[:-3]
else:
	module=None

################################
#poset is calculated by a module
################################
if module!=None:
	try:
		module=importlib.import_module(module)
	except ModuleNotFoundError:
		print("Error importing module ",module)
		exit()
	M=module.M
	ranks=module.ranks

############################
#poset is specified by input
############################
else:
	if file!=None:
		try:
			f = open(file, "rb")

			input = []
			#first strip down input
			for line in f:
				s = line[:line.find(b'#')] #if # in line cut out that and all after, if not cut off \n
				for c in s:
					if c in b';,0123456789': input.append(chr(c))
			f.close()
			input=''.join(input)

		except IOError:
			print("Error opening ",file)
			exit()
	else:
		if len(sys.argv)<2:
			print(help_string)
			exit()
		input=sys.argv[1]

	########################
	#Make incidence matrix M
	#and ranks list ranks
	#from input
	########################

	M = [] #Incidence matrix
	#turn input into a list of lists of numbers, ranks ended by an empty list
	input=[[int(x) for x in y.split(',') if x!=''] for y in input.split(';')[:-1]]
	debugPrint('input', input)

	#now find ranks
	n = 1 #1 less than number of elements
	ranks = [[0]]
	temp = []
	offset = 1
	for i in range(0,len(input)):
		if input[i] == []:
			ranks.append(temp)
			temp = []
			offset -= 1
			continue
		temp.append(i+offset)
		n += 1

	#get coatoms
	m = 0
	i = len(input)-2
	#number of coatoms is largest element of the last segment of lists in input
	m=max([max(x) for x in input[len(input)-(input[-2::-1]+[[]]).index([]):-1]])
	n += m

	debugPrint('n', n)
	debugPrint('m',m)
	ranks.append(list(range(ranks[-1][-1]+1,ranks[-1][-1]+1+m)))
	ranks.append([n])
	rank = len(ranks)

	#now make cover matrix from input
	#cover matrix
	M = [[0] + [1]*len(ranks[1]) + [0]*(n-len(ranks[1]))] #zero hat row
	r = 1 #current rank
	for i in range(0,len(input)):
		if input[i] == []:
			r += 1
			continue

		temp = [0]*(n+1)
		for x in input[i]:
			temp[x + ranks[r][-1]] = 1
		M.append(temp)
	#coatom rows
	for i in range(0, len(ranks[-2])):
		M.append([0]*(n)+[1])

	#onehat
	M.append([0]*(n+1))

	if 'debug' in sys.argv:
		print('cover matrix:')
		for r in M:
			for x in r:
				print(x,end='')
				print(" "*(3-len(str(x))),end='')
			print('')

	transClose(M)

if 'debug' in sys.argv:
	print('incidence matrix:')
	for r in M:
		for x in r:
			print(x,end='')
			print(" "*(3-len(str(x))),end='')
		print('')


debugPrint('ranks',ranks)
#############################
#calculate and print cd-index
#############################
table=[]
ab=[]
cd=cdIndex(M,ranks,table=table,ab=ab)
debugPrint('table',table)

print(cdIndexLatex(cd))

#check all intervals for Eulerian for debugging
if 'debug' in sys.argv:
	for i in range(0,len(M)):
		for j in range(0,len(M[i])):
			if M[i][j] != 1: continue
			temp = mobius(M,ranks,i,j)
			if temp != (-1)**(rankOf(ranks,j)-rankOf(ranks,i)):
				print("non Eulerian interval: [" + str(i) + ", " + str(j) + "]",end='')
				print("mobius(" + str(i) + ", " + str(j) + "): " + str(temp))

