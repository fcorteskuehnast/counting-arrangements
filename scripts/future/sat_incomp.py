# run with "python sat_incomp.py && KCBox_ExactMC x.cnf"

from itertools import *

#E = [(2,0), (2,1), (3,1), (4,1), (4,3), (5,1), (5,3), (5,4), (7,1), (7,3), (7,4), (7,5), (7,6), (8,6),]

E = [(1,0), (4,0), (4,1), (4,2), (4,3), (5,3), (6,2), (6,3), (6,5), (7,3), (7,5), (9,2), (9,3), (9,5), (9,7), (9,8), (10,0), (10,1), (10,2), (10,3), (10,5), (10,6), (10,7), (10,8), (10,9), (11,2), (11,3), (11,5), (11,7), (11,8), (11,9), (12,3), (12,5), (12,8), (14,0), (14,2), (14,3), (14,5), (14,6), (14,7), (14,8), (14,9), (14,11), (14,12), (14,13), (15,13), (16,2), (16,3), (16,5), (16,7), (16,8), (16,9), (16,12), (16,13), (16,15), (17,3), (17,8), (17,13), (17,15)]

n = max(max(e) for e in E)+1
N = range(n)

all_variables = []

all_variables += [('lt',I) for I in permutations(N,2)] 


all_variables_index = {}

num_vars = 0
for v in all_variables:
	all_variables_index[v] = num_vars
	num_vars += 1

def var(L):	return 1+all_variables_index[L]

def var_lt(*L): return var(('lt',L))


constraints = []

for u,v in permutations(N,2):
	constraints.append([-var_lt(u,v),-var_lt(v,u)])
	constraints.append([+var_lt(u,v),+var_lt(v,u)])
	
for u,v,w in permutations(N,3):
	constraints.append([-var_lt(u,v),-var_lt(v,w),+var_lt(u,w)])

for u,v in E:
    constraints.append([var_lt(u,v)])
	
output = "x.cnf"
if output:
	#args.output = "_h6_rings_"+str(n_H)+"_"+str(n_I)+"_"+str(n_J)+"_"+str(n_outer)+".in"
	print ("write cnf instance to file:",output)
	with open(output,"w") as f:
		f.write("p cnf "+str(num_vars)+" "+str(len(constraints))+"\n")

		for c in constraints:
			f.write(" ".join(str(x) for x in c)+" 0\n")