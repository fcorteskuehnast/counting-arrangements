#!/usr/bin/python
from sys import argv
import xml.etree.ElementTree as ET
from itertools import combinations

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("ipefile",type=str,help="ipefile")
parser.add_argument("--precision",default=1.,type=float,help="precision for crossings")
parser.add_argument("--figsize",default=10,type=int,help="figsize")
parser.add_argument("--visualize","-v",action='store_true',help="visualize")
parser.add_argument("--visualizeformat","-vf",choices=['pdf','png'],default='pdf',help="plot output format")
parser.add_argument("--visualizeoriginal","-vo",action='store_true',help="by default only visualized pruned segments")
parser.add_argument("--debug",type=int,default=0,help="debug level")
#parser.add_argument("--check-blue",action='store_true',help="check blue reference shapes in ipe files")
args = parser.parse_args()

if args.debug: print("args",args)


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
	return (p_h[1]/p_h[0],p_h[2]/p_h[0])

def dist(p,q):
	return (vector(p)-vector(q)).norm()


def simplify_bip(bip0):
	bip = bip0[:] # copy
	i = 1
	while i < len(bip):
		if bip[i-1] == bip[i]:
			del bip[i]
			del bip[i-1]
			if i>1: i -= 1
		
		else:
			i += 1

	while bip and bip[0] == bip[-1]:
		del bip[-1]
		del bip[0]


	values = {}
	bip_simp = []
	for i in bip:
		if i not in values:
			values[i] = len(values)
		bip_simp.append(values[i])
	return bip_simp

### parse ipe file ###
assert(os.path.isfile(args.ipefile))
tree = ET.parse(args.ipefile)
root = tree.getroot()
page = root.find('page')

segments_original = []
red_polygons = []
blue_polygons = []
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

		x,y = [float(z) for z in words[:2]]
		if 'matrix' in attr:
			x0 = x
			y0 = y
			x = M[0]*x0+M[2]*y0+M[4]
			y = M[1]*x0+M[3]*y0+M[5]
		p = (x,y)
		pts.append(p)

	if is_poly:
		if attr["stroke"] == "red":
			#assert(not poly_red) # must be unique
			red_polygons.append(pts)
		if attr["stroke"] == "blue":
			#assert(not poly_blue) # must be unique
			blue_polygons.append(pts)
	else:
		if len(pts)!=2:
			print("ERROR: line with more than 2 points:",pts)
		assert(len(pts)==2)
		if pts in segments_original:
			print("warning: duplicate lines in ipe file!")
		else:
			segments_original.append(pts)

if args.debug >=2:
	print("segments_original:",segments_original)
	print("red_polygons:",red_polygons)
	print("blue_polygons:",blue_polygons)

if blue_polygons: assert(len(blue_polygons) == len(red_polygons))
assert(red_polygons)


outfile = f"{args.ipefile}.polygons"
print(f"write number of red polygons ({len(red_polygons)}) to file: {outfile}")
with open(outfile,"w") as f:
	f.write(str(len(red_polygons))+"\n")


for poly_index,poly_red in enumerate(red_polygons):
	print(f"****** handle red polygon #{poly_index} ******")

	### get bip ###
	bip = {}
	bip_pos = {}

	segments_pruned = []
	if args.debug >= 1 : print("compute bip and prune segments")
	for i in range(len(segments_original)):
		crossings_along_i = []

		for j in range(len(poly_red)):
			s = [poly_red[j-1],poly_red[j]]
			s_vec = vector(s[1])-vector(s[0])
			
			if segments_intersect(s,segments_original[i]):
				p_inter = intersection(s,segments_original[i])

				vec2 = vector(p_inter)-vector(s[0])
				t = (s_vec*vec2)/(s_vec*s_vec)
				assert(t > 0 and t < 1)
		
				for q in bip_pos.values():
					assert(dist(p_inter,q) > args.precision) # redraw red polygon or increase precision if this happens, but be careful

				assert((j+t) not in bip_pos) # redraw red polygon or increase precision if this happens, but be careful
				bip_pos[j+t] = p_inter
				crossings_along_i.append(p_inter)

		# prune
		assert(len(crossings_along_i)%2 == 0) # even
		crossings_along_i.sort()
		for p,q in combinations(crossings_along_i,2): 
			assert(p[0]<q[0] or (p[0]==q[0] and p[1]<q[1]))
		for i in range(0,len(crossings_along_i),2):
			segments_pruned.append([crossings_along_i[i],crossings_along_i[i+1]])



	if args.debug >= 2: print("compute bip")
	bip_pos = [bip_pos[t] for t in sorted(bip_pos)] # crossings ordered along polygon
	bip = [i for pos in bip_pos for i,s in enumerate(segments_pruned) if pos in s]


	if args.debug >= 1: print("bip",bip)

	for x in bip:
		assert(bip.count(x)==2) # non-convex polygons not implemented yet

	simple_bip = simplify_bip(bip)
	if args.debug >= 1: print("simplified bip",simple_bip)

	outpath = f"{args.ipefile}_poly{poly_index}.bip"
	print(f"write simplified bip to file: {outpath}")
	with open(outpath,"w") as f:
		f.write(' '.join(str(i) for i in simple_bip)+"\n")


	### compute internal crossings ###
	if args.debug >= 2: print("compute crossings")
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
	if args.debug >= 1: print("cross_stat",cross_stat)


	if args.debug >= 2: print("precision:",args.precision)
	if len(crossings) >= 2:
		cr_min_dist = min(dist(p,q) for p,q in combinations(crossings,2))
		if args.debug >= 2: print("min distance between crossings:",cr_min_dist)

	pos_min_dist = min(dist(p,q) for p,q in combinations(bip_pos,2))
	if args.debug >= 2: print("min distance between intersections on polygon:",pos_min_dist)



	outpath = f"{args.ipefile}_poly{poly_index}.crossings"
	print(f"write information on crossings to file: {outpath}")
	with open(outpath,"w") as f:
		f.write(str(cross_stat)+"\n")


	if args.visualize:
		plt = []
		if args.visualizeoriginal:
			for i,l in enumerate(segments):
				plt.append(line(l,color='gray'))
				for p in l:
					plt.append(text(str(i),p,color="gray"))

		for l in segments_pruned:
			plt.append(line(l,color='black'))

		plt.append(polygon(poly_red,color='red',fill=False))
		if args.visualizeoriginal:
			plt.append(polygon(poly_blue,color='blue',fill=False))

		for j in range(len(poly_red)):
			s0 = poly_red[j-1] 
			plt.append(text(f'p{j}',vector(s0)+vector([0,1]),color="red",zorder=2))
			plt.append(point2d(s0,color='red',zorder=2))


		multiplicities = {len(lines) for p,lines in crossings.items()}
		rainbow = sage.plot.colors.rainbow(len(multiplicities))
		rainbow.reverse()
		color = {k:rainbow[i] for i,k in enumerate(sorted(multiplicities))}

		for p,lines in crossings.items():
			multiplicity = len(lines)
			plt.append(text(str(multiplicity),vector(p)+vector([0,1]),color=color[multiplicity],zorder=2))
			plt.append(point2d(p,color=color[multiplicity],zorder=2))

		for i in range(len(bip)):
			plt.append(text(str(bip[i]),vector(bip_pos[i])+vector([0,1]),color='gray',zorder=2))
			plt.append(point2d(bip_pos[i],color='gray',zorder=2))

		outpath = f"{args.ipefile}_poly{poly_index}.{args.visualizeformat}"
		print(f"write plot to file: {outpath}")
		sum(plt).save(outpath,axes=False,figsize=(args.figsize,args.figsize))

