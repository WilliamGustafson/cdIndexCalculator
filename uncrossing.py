from cdIndex import debugPrint,transClose
import sys

#Throughout pairings are encoded as lists of numbers, each number encodes
#a pair as two bits set. For example the pairing {{1,3},{2,4}} is encoded
#as [2**0+2**2,2**1+2**3]=[5,10]

#converts the pairing given in the input into the internal format described above
def readPairing(input):
	input = input.split(',')
	t = []
	for i in range(0,len(input)//2):
		t.append(1<<(int(input[i<<1])-1)|1<<(int(input[(i<<1)+1])-1))
	return sorted(t)

def setFormat(x):
	ret = []
	i = 1
	while x!= 0:
		if x&1: ret.append(str(i))
		x >>= 1
		i += 1
	return "{"+",".join(ret)+"}"

def pairingFormat(x):
	return "{"+",".join([setFormat(y) for y in x])+"}"


#swaps i and j in the pairing p
def swap(p,i,j):
#	return sorted([(x^((((x&(1<<i))>>i)^((x&(1<<j))>>j))<<i)^((((x&(1<<i))>>i)^((x&(1<<j))>>j))<<j)) for x in p])
	ret = []
	for arc in p:
		if (arc&(1<<i))>>i != (arc&(1<<j))>>j: #arc contains one of i and j
			ret.append(arc ^ ((1<<i)|(1<<j))) #swap i and j in the arc
		else: #arc contains both i and j or neither so don't swap i and j
			ret.append(arc)
	return sorted(ret)

#returns the number of crossings of p
def c(p):
	ret = 0
	for i in range(0,len(p)):
		xi = bin(p[i])[::-1]
		Ni = xi.find('1')
		Ei = xi.rfind('1')
		for j in range(i+1,len	(p)):
			xj = bin(p[j])[::-1]
			Nj = xj.find('1')
			Ej = xj.rfind('1')
			if (Ni - Nj > 0) == (Ei - Ej > 0) == (Nj - Ei > 0): ret += 1

	return ret

#computes the lower interval generated by the given pairing
#returns a tuple (P,ranks,M) which is the list of elements, the rank list and the incidence matrix
def lowerOrderIdeal(t):
	if c(t)==0: return [t],[[1],[0]],[[0,-1],[1,0]]

	P=[t]
	ranks = [[0]] #this is built up backwards for convenience and reversed before returning
	M=[[0]]

	num = 1 #index in to P of next element to add
	level = [t] #list of current rank to expand in next step
	leveli = [0] #indices in to P of the elements of level
	newLevel = [] #we build level for the next step during the current step here
	newLeveli = [] #indices in to P for the next step
	newRank = [] #the new rank indices to add
	while len(level) > 0:
		for i in range(0,(len(t)<<1)-1): #iterate over all pairs we can uncross
			for j in range(i+1,len(t)<<1):
				for k in range(0,len(level)): #do the uncross
					temp = swap(level[k],i,j)
					c_temp = c(temp)
					if c_temp != c(level[k])-1: continue
					if temp in P:
						M[P.index(temp)][leveli[k]]=1
						continue
					P.append(temp)
					newRank.append(num)
					if c_temp > 0: #if not minimal continue uncrossing
						newLevel.append(temp)
						newLeveli.append(num)
					num+= 1

					for x in M: x.append(0)
					M.append([0 for x in range(0,len(M[0]))])
					M[-1][leveli[k]]=1

		level = newLevel
		newLevel = []
		leveli = newLeveli
		newLeveli = []
		ranks.append(newRank)
		newRank = []

	ranks.reverse()
	ranks=[[len(M)]]+ranks
	for r in M: r.append(-1)
	M.append([1 for i in range(0,len(P))]+[0])
	transClose(M)
	return P,ranks,M

#input should be first argument discounting the two used to specify this module
input=[sys.argv[i] for i in range(0,len(sys.argv)) if i!= sys.argv.index("-m") and i-1!= sys.argv.index("-m")][1]
t = readPairing(input)
debugPrint('t',t)
P,ranks,M = lowerOrderIdeal(t)
debugPrint("P",P)
