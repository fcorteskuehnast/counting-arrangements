#!/usr/bin/python
from sys import argv
import xml.etree.ElementTree as ET
from itertools import combinations

import argparse

import pprint
pp = pprint.PrettyPrinter()

parser = argparse.ArgumentParser()
parser.add_argument("ipefile",type=str,help="ipefile")
parser.add_argument("--precision",default=1.,type=float,help="precision for crossings")
args = parser.parse_args()
print("args",args)

def polygon_area(vertices):
    n = len(vertices) 
    area = 0
    for i in range(n):
        area += vertices[i  ][0] * vertices[i-1][1]
        area -= vertices[i-1][0] * vertices[i  ][1]
    return abs(area)/2


tree = ET.parse(args.ipefile)
root = tree.getroot()
page = root.find('page')

segments_original = []
regions = {}

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
		poly_color = attr['fill']

		if poly_color in regions.keys(): 
			regions[poly_color].append(pts)
		else:
			regions[poly_color] = [pts]
		
	else:
		if len(pts)!=2:
			print("ERROR: line with more than 2 points:",pts)
		assert(len(pts)==2)
		if pts in segments_original:
			print("warning: duplicate lines in ipe file!")
		else:
			segments_original.append(pts)

print("segments_original:",segments_original)

def map_vals(f, d):
    return dict({ key : f(val) for (key, val) in d.items()})

print("\nregions: ", regions)

areas = {}
for (color, polys) in regions.items():
	areas[color] = sum(map(polygon_area, polys))

print("\nareas:")
pp.pprint(areas)

print("\nnormalized areas:")
n_areas = map_vals(lambda x : x / ((2*64)**2 * 64), areas)
pp.pprint(n_areas)

print("\ndifferences:")


