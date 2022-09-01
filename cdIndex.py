import sys

def debugPrint(s, x=None):
	if 'debug' in sys.argv:
		print('')
		if x ==None: print(s)
		else:
			print(s + ": ",end='')
			print(x)
		print('')


##########################################
#a few general poset utilities
#
#within this module a poset is represented via matrix M
#whose entry M[i][j] is 1 if i<j, -1 if i>j and 0 otherwise.
#
#many functions also require a rank list provided whose ith entry
#is a list of the indices corresponding to rank i elements.
##########################################

def rankOf(ranks,x):
	for i in range(0, len(ranks)):
		if x in ranks[i]:
			return i
	return -1

#given a matrix encoding all the cover relations of a poset this function
#alters the input to compute the matrix representing the poset
def transClose(M):
	for i in range(0,len(M)):
		#uoi for upper order ideal
		uoi = [x for x in range(0,len(M)) if M[i][x] == 1]
		while True:
			next = [x for x in uoi]
			for x in uoi:
				for y in range(0,len(M)):
					if M[x][y] == 1 and y not in next: next.append(y)
			if uoi == next: break
			uoi = next

		for x in uoi:
			M[i][x] = 1
			M[x][i] = -1

#computes the join of i and j in M, returns -1 if it does not exist
def join(i,j,M):
	if i==j: return i
	if M[i][j] == -1: return i
	if M[i][j] == 1: return j
	m = [x for x in range(0,len(M)) if M[i][x] == 1 and M[j][x] == 1]
	for x in range(0,len(m)):
		isJoin = True
		for y in range(0,len(m)):
			if x!=y and M[m[x]][m[y]] != 1:
				isJoin = False
				break
		if isJoin: return m[x]
	return None

#computes mobius function of [i,j], useful for debugging
def mobius(M, ranks, i, j):
	if i == j: return 1
	if M[i][j] != 1: return None
	ret = 0
	for x in range(0,len(M)):
		if M[x][j] == 1 and M[x][i] == -1:
			ret -= mobius(M,ranks,i,x)
	ret -= 1 #mobius(M,rank,i,i)
	return ret


##########################################
#helper functions for polynomials with noncommutative variables
#
#monomials are stored as a list containing the coefficient and a string for the variables
#polynomials are lists of monomials
##########################################

def multPolys(p,q):
	r=[[x[0]*y[0],x[1]+y[1]] for x in p for y in q]
	#collect terms
	ret=[]
	for x in r:
		monoms=[y[1] for y in ret]
		if x[1] not in monoms:
			ret.append(x)
			continue
		ret[monoms.index(x[1])][0]+=x[0]
	return ret

def multManyPolys(P):
	ret=P[0]
	while len(P)>1:
		P=P[1:]
		ret=multPolys(ret,P[0])
	return p

def addPolys(p,q):
	ret=[x for x in p]
	for x in q:
		temp=[y[1] for y in ret]
		if x[1] in temp:
			ret[temp.index(x[1])][0]+=x[0]
		else:
			ret.append(x)
	return [x for x in ret if x[0]!=0]

def subPolyForMonom(x,p,m):
	X=[[y[0],y[1].replace(m,'*')] for y in x]
	ret=[]
	for y in X:
		q=[[y[0],'']]
		for i in range(0,len(y[1])):
			if y[1][i]=='*':
				q=multPolys(q,p)
			else: #mult by the monomial
				for j in range(0,len(q)):
					q[j][1]+=y[1][i]
		ret=addPolys(ret,q)
	return ret

##########################################
#methods to calculate flag vector table
##########################################

#fVector returns the flag f vector f_S(M)
#fVectorCalc calculates this recursively

def fVectorCalc(ranks,S,M, i, count):
	newCount = count
	if S == []: return 1
	for j in ranks[S[0]]:
		if M[i][j] == 1:
			newCount += fVectorCalc(ranks, S[1:], M, j, count)
	return newCount

def fVector(M,ranks,S):
	return fVectorCalc(ranks,S,M,ranks[0][0],0)

#returns a list whose entries are lists containing a subset of the ranks
#and the flag f and h vector entries
def makeFlagVectorsTable(M,ranks):
	table = [[[],1,1]]

	if len(ranks)<=2: return table

	#iterate over all subsets of the ranks
	for i in range(1,1<<(len(ranks)-1)-1):
		#construct the corresponding set S
		pad = 1
		elem = 1
		S = []
		while pad <= i:
			if pad&i:
				S.append(elem)

			pad <<= 1
			elem += 1
		table.append([S,fVector(M,ranks,S),0])


	#do PIE to get the flag h vectors
	for i in range(1,len(table)):
		sign = (2*(len(table[i][0])%2)) - 1 #is -1 if even number of elements, 1 if odd
		for j in range(0,i+1):
			if set(table[j][0]).issubset(table[i][0]):
				table[i][2] += sign*(2*(len(table[j][0])%2)-1)*table[j][1]
	return table

##########################################
#methods for computing the ab and cd indices
##########################################

def abIndex(table, rank):
	ab = []
	for x in table:
		u = ['a']*(rank-1)
		for s in x[0]: u[s-1] = 'b'
		ab.append([x[2],''.join(u)])

	return ab

#returns a given ab polynomial expressed in terms of c and d if possible
#if not possible returns ab
def abToCd(ab):
	if len(ab)==0: return ab
	#substitue a->c+e and b->c-e
	#where e=b-a
	#this scales by a factor of 2^n (n+1 is number of ranks)
	ce=subPolyForMonom(subPolyForMonom(ab,[[1,'c'],[1,'e']],'a'),[[1,'c'],[-1,'e']],'b')
	debugPrint('ce',ce)

	cd=subPolyForMonom(ce,[[1,'cc'],[-2,'d']],'ee')
	for m in cd:
		if 'e' in m[1]: return ab
	#divide coefficients by 2^n
	power=sum([2 if cd[0][1][i]=='d' else 1 for i in range(len(cd[0][1]))])
	return [[x[0]>>power,x[1]] for x in cd]

#returns the cd-index of the poset encoded by M
#
#optionally provide a list for table or ab and these will be
#filled with the flag vectors table and ab-index respectively

def cdIndex(M, ranks, table=None, ab=None):
	if table == None: table = []
	if ab == None: ab = []

	table += makeFlagVectorsTable(M,ranks)

	debugPrint('table', table)

	ab += abIndex(table, len(ranks)-1)

	debugPrint('ab',ab)
	cd=abToCd(ab)
	cd.sort(key=lambda x:x[1])
	debugPrint('cd', cd)
	return cd

##########################################
#formatting methods
##########################################

#returns a string representing the given polynomial in latex
def cdIndexLatex(cd):
	s = ""
	for i in range(0,len(cd)):
		if cd[i][0] == 0: continue
		if cd[i][0] == -1: s+= '-'
		elif cd[i][0] != 1: s += str(cd[i][0])
		current = ''
		power = 0
		for c in cd[i][1]:
			if current == '':
				current = c
				power = 1
				continue
			if c == current:
				power += 1
				continue
			s += current
			if power != 1: s += '^{' + str(power) + '}'
			current = c
			power = 1
		s += current
		if power != 1 and power != 0: s += '^{' + str(power) + '}'
		if power == 0 and current == "": s += '1'

		if i != len(cd)-1:
			if cd[i+1][0] >= 0: s += "+"
	if s == '': return '0'
	return s

#returns a string with latex code to display the given flag vectors table
def flagVectorsLatex(table):
	ret = "\\begin{longtable}{c|c|c}\n\t$S$&$f_S$&$h_S$\\\\\n\t\\hline\n\t\\endhead\n"
	for t in table:
		ret += "\t\{"
		ret += ','.join([str(x) for x in t[0]])
		ret = ret+"\} & " + str(t[1]) + " & " + str(t[2]) + "\\\\\n\t\\hline\n"
	return ret + "\\end{longtable}"
