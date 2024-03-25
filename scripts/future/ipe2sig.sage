#!/usr/bin/python
from sys import argv
import xml.etree.ElementTree as ET
from itertools import combinations

import argparse



# compute orientation of point triple
def o3(p,q,r): 
	return (q[0]-p[0])*(r[1]-p[1])-(r[0]-p[0])*(q[1]-p[1])

# test whether two segments intersect
def segments_intersect(e,f):
	return (o3(e[0],e[1],f[0])*o3(e[0],e[1],f[1])<0) and (o3(f[0],f[1],e[0])*o3(f[0],f[1],e[1])<0)

# compute crossing point of two segments/lines
def segment2vec3(e):
	dx = e[1][0]-e[0][0]
	dy = e[1][1]-e[0][1]
	#assert(vector([dx,dy]) * vector([-dy,dx]) == 0) # normal
	h1 =  -e[0][0]*dy + e[0][1]*dx
	h2 =  -e[0][0]*dy + e[0][1]*dx
	assert(h1==h2)
	#assert(vector([-h1,-dy,dx])*vector([1,e[0][0],e[0][1]]) == 0)
	#assert(vector([-h1,-dy,dx])*vector([1,e[1][0],e[1][1]]) == 0)
	return vector([-h1,-dy,dx])

def intersection(e,f):
	e_v = segment2vec3(e)
	f_v = segment2vec3(f)
	p_h = e_v.cross_product(f_v) 
	return (p_h[1]/p_h[0],p_h[2]/p_h[0]) if p_h[0] != 0 else None

def relative_intersection(e,f):
	p = intersection(e,f)
	if p == None: 
		return None
	else:
		e_vec = vector(e[1])-vector(e[0])
		p_vec = vector(p)-vector(e[0])
		t = (e_vec*p_vec)/(e_vec*e_vec) # segment starts at t=0 and ends at t=1
		return t

def isparallel_float(e,f):
	p = intersection(e,f)
	return p == None or vector(p).norm() > 10^6

def dist(p,q):
	return (vector(p)-vector(q)).norm()

def slope(p,q): 
	dx = q[0]-p[0]
	dy = q[1]-p[1]
	return dy/dx if dx != 0 else +Infinity 

def offset(p,q): 
	dx = q[0]-p[0]
	dy = q[1]-p[1]
	return p[1]-dy*p[0]/dx if dx != 0 else p[0] 

def convex_position(pts):
	return len(pts) == len(Polyhedron(pts).vertices())

def signotope_completions(n,partial_signotope):
	from pysat.formula import IDPool,CNF
	vpool = IDPool()
	V = {I:vpool.id(I) for I in combinations(range(n),3)}
	print("V",V)
	cnf = CNF()
	for I,s in partial_signotope.items():
		if s != 0:
			cnf.append([s*V[I]]) # fix + and - symbols, 0 indicates free variables
	for I in combinations(range(n),4):
		II = list(combinations(I,3))
		for A,B,C in combinations(II,3):
			cnf.append([+V[A],-V[B],+V[C]])
			cnf.append([-V[A],+V[B],-V[C]])
	print(cnf.clauses)
	cnf.to_file("x.cnf")
	exit()


def ipe2sig(args):
	### parse ipe file ###
	tree = ET.parse(args.ipefile)
	root = tree.getroot()
	page = root.find('page')

	segments = []
	poly_red = None
	poly_blue = None
	for u in page.iterfind('path'):
		attr = u.attrib
		if 'matrix' in attr:
			M = [float(t) for t in attr['matrix'].split(" ")]

		color = attr['stroke'] if 'stroke' in attr else None

		lines = u.text.split("\n")

		pts = []
		is_poly = False
		for l in lines:
			if l == '': continue
			words = l.split()

			if words[-1] == 'h': 
				is_poly = True
				continue

			x,y = [RR(z) for z in words[:2]]
			if 'matrix' in attr:
				x0 = x
				y0 = y
				x = M[0]*x0+M[2]*y0+M[4]
				y = M[1]*x0+M[3]*y0+M[5]
			p = (x,y)
			pts.append(p)

		if is_poly:
			if attr["stroke"] == "red":
				assert(not poly_red) # must be unique
				poly_red = pts
			if attr["stroke"] == "blue":
				assert(not poly_blue) # must be unique
				poly_blue = pts
		else:
			if len(pts)!=2:
				print("ERROR: line with more than 2 points:",pts)
			assert(len(pts)==2)
			if pts in segments:
				print("warning: duplicate lines in ipe file!")
			else:
				segments.append(pts)

	print("segments:",segments)
	print("poly_red:",poly_red)
	print("poly_blue:",poly_blue)

	assert(poly_blue)
	assert(poly_red)

	assert(convex_position(poly_red))
	
	n = len(segments)

	segments.sort(key=lambda s: (slope(*s),offset(*s)))
	# sort segments by increasing slope
	# vertical lines have +Infinity slope and will be sorted left to right

	# orient lines from left to right, vertical bottom to top
	for i in range(n):
		p,q = segments[i]
		dx = q[0]-p[0]
		dy = q[1]-p[1]
		if dx < 0 or (dx == 0 and dy < 0):
			segments[i] = (q,p)


	# omit segments which do not intersect the shape
	relevant_segments = []
	segments_pruned = []
	for i in range(n):		
		s = segments[i]
		bound_int = []
		for j in range(len(poly_red)):
			bs = [poly_red[j-1],poly_red[j]] # boundary segment
			
			if segments_intersect(s,bs):
				t = relative_intersection(s,bs)
				p = intersection(s,bs)
				bound_int.append((t,j,p))
				assert(t > 0 and t < 1) # intersection should lie on segment
		
		if bound_int:
			assert(len(bound_int) == 2)
			relevant_segments.append(segments[i])
			segments_pruned.append([p for (t,j,p) in bound_int])

	segments = relevant_segments
	n0 = n
	n = len(segments)
	print(f"only {n} of {n0} segments intersect the given shape")


	# compute partial signotope
	partial_signotope = {}

	for i in range(n):
		s = segments[i]
		s_vec = vector(s[1])-vector(s[0])

		rest_int = []
		rest_par = []
		for j in range(i+1,n):
			t = relative_intersection(s,segments[j])
			if t == None:
				rest_par.append(j)
			else:
				t = round(t,2) #TODO: input param precision
				rest_int.append((t,j))

		rest_int.sort()

		bound_int = []
		for j in range(len(poly_red)):
			bs = [poly_red[j-1],poly_red[j]] # boundary segment
			
			if segments_intersect(s,bs):
				t = relative_intersection(s,bs)
				bound_int.append((t,j))
				assert(t > 0 and t < 1) # intersection should lie on segment
		
		assert(len(bound_int) == 2)
		l = bound_int[0][0]
		r = bound_int[1][0]
		L = [j for (t,j) in rest_int if t<l]
		M = [j for (t,j) in rest_int if t>l and t<r]
		R = [j for (t,j) in rest_int if t>r] + rest_par
		J = L+M+R

		for j,k in combinations(range(i+1,n),2):
			if j in M and k in M:
				partial_signotope[i,j,k] = 0
			else:
				partial_signotope[i,j,k] = +1 if J.index(j) < J.index(k) else -1

	print("partial_signotope",partial_signotope)
	signotope_completions(n,partial_signotope)
	exit()



	### compute internal crossings ###
	print("compute crossings")
	crossings = {}
	n = len(segments_pruned)
	for i,j in combinations(range(n),2):
		if segments_intersect(segments_pruned[i],segments_pruned[j]):
			p_inter = intersection(segments_pruned[i],segments_pruned[j])		
			
			match = 0
			for q in crossings:
				if dist(p_inter,q) < args.precision:
					crossings[q] |= {i,j}
					match += 1

			assert(match <= 1)
			if not match:
				crossings[p_inter] = {i,j}
			
	#print("crossings",crossings.values())


	cross_stat = {}
	for c in crossings.values():
		k = len(c)
		if k not in cross_stat: 
			cross_stat[k] = 0
		cross_stat[k] += 1
	print("cross_stat",cross_stat)


	print("precision:",args.precision)
	if len(crossings) >= 2:
		cr_min_dist = min(dist(p,q) for p,q in combinations(crossings,2))
		print("min distance between crossings:",cr_min_dist)


	outpath = args.ipefile+".crossings"
	print(f"write information on crossings to file {outpath}")
	with open(outpath,"w") as f:
		f.write(str(cross_stat)+"\n")


	if args.visualize:
		plt = []
		if args.visualizeoriginal:
			for i,l in enumerate(segments):
				plt.append(line(l,color='gray'))
				for p in l:
					plt.append(text(str(i),p,color="gray"))

		for i,l in enumerate(segments_pruned):
			plt.append(line(l,color='black'))
			for p in l:
				plt.append(text(str(i),vector(p)+vector([0,1]),color='gray',zorder=2))
				plt.append(point2d(p,color='gray',zorder=2))

		plt.append(polygon(poly_red,color='red',fill=False))
		if args.visualizeoriginal:
			plt.append(polygon(poly_blue,color='blue',fill=False))

		for j in range(len(poly_red)):
			s0 = poly_red[j-1] 
			plt.append(text(str(j),vector(s0)+vector([0,1]),color="red",zorder=2))
			plt.append(point2d(s0,color='red',zorder=2))

		for p,lines in crossings.items():
			plt.append(text(str(len(lines)),vector(p)+vector([0,1]),color='green',zorder=2))
			plt.append(point2d(p,color='green',zorder=2))

		outpath = args.ipefile+"."+args.visualizeformat
		print(f"write plot to file {outpath}")
		sum(plt).save(outpath,axes=False,figsize=(args.figsize,args.figsize))



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("ipefile",type=str,help="ipefile")
	parser.add_argument("--precision",default=1.,type=float,help="precision for crossings")
	parser.add_argument("--figsize",default=10,type=int,help="figsize")
	parser.add_argument("--visualize","-v",action='store_true',help="visualize")
	parser.add_argument("--visualizeformat","-vf",choices=['pdf','png'],default='pdf',help="plot output format")
	parser.add_argument("--visualizeoriginal","-vo",action='store_true',help="by default only visualized pruned segments")
	args = parser.parse_args()
	print("args",args)

	ipe2sig(args)


