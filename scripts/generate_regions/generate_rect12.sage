from sympy import Point, Polygon, convex_hull, Line, oo, zoo, sympify
import xml.etree.ElementTree as ET
from xml.dom import minidom
from pathlib import Path	


OUTPUT_FOLDER = "output/"

# IPE stuff

def prettify(elem):
	rough_string = ET.tostring(elem, 'utf-8')
	reparsed = minidom.parseString(rough_string)
	return reparsed.toprettyxml(indent="  ")

def pt_in_poly(p,poly):
	if p in poly.vertices or any(p in s for s in poly.sides):
		return True
	return poly.encloses_point(p)

def intersect_polygons(poly1,poly2):
	pts = list(poly1.intersection(poly2)) + list(poly1.vertices) + list(poly2.vertices)
	inside = []
	for p in pts:
		if pt_in_poly(p,poly1) and pt_in_poly(p,poly2):
			inside.append(p)
	return convex_hull(*inside)
	
	
def add_line(line, page, layer, w, h):
	square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
	intersection = square.intersection(line)
	if len(intersection) == 2:
		p1,p2 = intersection
		p1 = ( float(50*(p1[0]+w/2)), float(50*(p1[1]+h/2)) )
		p2 = ( float(50*(p2[0]+w/2)), float(50*(p2[1]+h/2)) )
		path = ET.SubElement(page, "path")
		path.set("layer", layer)
		path.text = "{} {} m {} {} l".format(p1[0],p1[1],p2[0],p2[1])
	

def add_polygon(poly, page, layer, w, h):
	square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
	
	poly = intersect_polygons(poly,square)

	path = ET.SubElement(page, "path")
	path.set("layer", layer)
	path.set("fill", "red")
	path.set("opacity","50%")
	
	txt = ""
	p = poly.vertices[0]
	txt += "{} {} m ".format( float(50*(p[0]+w/2)), float(50*(p[1]+h/2)) )
	for p in poly.vertices[1:]:
		txt += "{} {} l ".format( float(50*(p[0]+w/2)), float(50*(p[1]+h/2)) )
	txt += "h"
	
	path.text = txt


def make_IPE(bundles, polygons, area, filename, w=10, h=10):
	
	tree = ET.parse("template.ipe")
	ipe = tree.getroot()
	
	ipestyle = ipe.find("ipestyle")
	
	preamble = ET.SubElement(ipestyle, "preamble")
	preamble.text = "area = {}".format(area)
	
	layout = ET.SubElement(ipestyle, 'layout')
	layout.set("paper", "{} {}".format(50*w, 50*h))
	layout.set("origin", "0 0")
	layout.set("frame", "{} {}".format(50*w, 50*h))
	
	page = ET.SubElement(ipe, "page")
	
	square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
	
	for i in range(len(bundles)):
		layer = ET.SubElement(page, "layer")
		layer.set("name", "bundle_{}".format(i))
	
	if len(polygons):
		layer = ET.SubElement(page, "layer")
		layer.set("name", "polygons")
		
	for i,bundle in enumerate(bundles):
		for line in bundle:
			add_line(line, page, "bundle_{}".format(i), w, h)
			
	for poly in polygons:
		add_polygon(poly, page, "polygons", w, h)

	
	xml_string = prettify(ipe)
	with open(OUTPUT_FOLDER+filename, "w") as text_file:
		text_file.write(xml_string)


# Generating regions and computing areas

s1 = sympify(1)

def get_bundles_old(slopes, shifts, w=10,h=10):
	square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
	
	bundles = []
	for s in slopes:
		bundle = []
		k = 1
		while True:
			
			i = k//2*(-1)**k
			k += 1
			if s == zoo:
				line = Line(Point(i,-100), Point(i,100))
			elif abs(s) >= 1 or s==0: 
				line = Line(Point(-100,-100*s+i), Point(100,100*s+i))
			else:
				line = Line(Point(-100+i,-100*s), Point(100+i,100*s))
			
			intersection = square.intersection(line)
			if len(intersection) == 2:
				bundle.append(line)
			else:
				break
			
		bundles.append(bundle)
	return bundles

def get_bundles(slopes, shifts, w,h):
	square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
	
	bundles = []
	for i in range(len(slopes)):
		slope, shift = slopes[i], Point(shifts[i])
		bundle = []
		k = 1
		while True:
			n = k/2*(-1)**k
			k += 1
			
			if slope == zoo:
				line = Line(Point(0,-1000)+shift*n, Point(0,1000)+shift*n)
			else:
				line = Line(Point(-1000,-1000*slope)+shift*n, Point(1000,1000*slope)+shift*n)
			
			intersection = square.intersection(line)
			if len(intersection) == 2:
				bundle.append(line)
			else:
				break
			
		bundles.append(bundle)
	return bundles
	
def all_slope_symmetries(slopes):
	sym1 = []
	for slope in slopes:
		sym1.append(s1/slope)
	sym2 = []
	for slope in slopes:
		sym2.append(-slope)
	sym3 = []
	for slope in sym1:
		sym3.append(-slope)
	
	yield tuple(sorted(slopes, key=lambda x: oo if x==zoo else x))
	yield tuple(sorted(sym1, key=lambda x: oo if x==zoo else x))
	yield tuple(sorted(sym2, key=lambda x: oo if x==zoo else x))
	yield tuple(sorted(sym3, key=lambda x: oo if x==zoo else x))


def generate_regions(slopes, shifts):
	
	num_bundles = len(slopes)
	
	slope_shift = {}
	for i in range(num_bundles):
		slope_shift[slopes[i]] = shifts[i]

	mid_lines = []
	for s in slopes:
		if s == zoo:
			mid_lines.append(Line( (0,0), (0,1) ))
		else:
			mid_lines.append(Line( (0,0), (1,s) ))

	square = Polygon(Point(-100,-100), Point(100,-100), Point(100,100), Point(-100,100))
	upper_lines = []
	lower_lines = []
	basic_regions = []
	for i in range(num_bundles): #line in enumerate(mid_lines):
		p1,p2 = square.intersection(mid_lines[i])
		shift = Point(shifts[i])
		upper_lines.append(Line(p1+shift/2, p2+shift/2))
		lower_lines.append(Line(p1-shift/2, p2-shift/2))
		basic_regions.append( convex_hull(p1+shift/2, p2+shift/2, p1-shift/2, p2-shift/2) )

	print("Generating regions...")
	regions = [square]
	for i,line in enumerate(upper_lines + lower_lines):
		print("{}/{}".format(i+1, len(upper_lines) + len(lower_lines)))
		new_regions = []
		for r in regions:
			if len(r.intersection(line)):
				r1, r2 = r.cut_section(line)
				if r1 and r1.area > 0:
					new_regions.append(r1)
				if r2 and r2.area > 0:
					new_regions.append(r2)
			else:
				new_regions.append(r)
		regions = new_regions

	print("Finding bundles in regions...")
	region_slopes = []
	for i,r in enumerate(regions):
		print("{}/{}".format(i+1, len(regions)))
		r_slopes = []
		r_shifts = []
		c = r.centroid
		for i,br in enumerate(basic_regions):
			if br.encloses_point(c):
				r_slopes.append(slopes[i])
		region_slopes.append(tuple(r_slopes))
  

	print("Grouping symmetric regions...")
	agreggated_regions = {}
	areas = {}
	
	for i in range(len(regions)):
		region_bundles = region_slopes[i]
		if len(region_bundles) < 3:
			continue
		a = regions[i].area
		if a == 0:
			continue
			
		region_bundles = tuple(sorted(region_bundles, key=lambda x: oo if x==zoo else x))

		rep = None
		for sym in all_slope_symmetries(region_bundles):
			rep = sym
			if sym in agreggated_regions:
				break
		if rep in agreggated_regions:
			agreggated_regions[rep].append(regions[i])
			areas[rep] += a
		else:
			agreggated_regions[rep] = [regions[i]]
			areas[rep] = a	 
	
	extremal_lines = []
	for i in range(num_bundles):
		extremal_lines.append([lower_lines[i],upper_lines[i]])
		
	print("Making Ipe files...")
	
	Path(OUTPUT_FOLDER+"regions/").mkdir(parents=True, exist_ok=True)
	Path(OUTPUT_FOLDER+"shapes/").mkdir(parents=True, exist_ok=True)
	make_IPE(extremal_lines, [], 0, "base_construction.ipe", w=8, h=8)
	
	for i,rep in enumerate(sorted(agreggated_regions, key=lambda x:-len(x))):
		print("{}/{}".format(i+1, len(agreggated_regions)))
		make_IPE(extremal_lines, agreggated_regions[rep], areas[rep], "regions/region_{}.ipe".format(i), w=8, h=8)
		rep_shifts = [slope_shift[s] for s in rep]
		grid_bundles = get_bundles(rep,rep_shifts, w=12, h=12)
		make_IPE(grid_bundles, [], areas[rep], "shapes/grid_{}.ipe".format(i), w=12, h=12)

slopes = [0,1,2,3,-1,-2,-3, s1/2, -s1/2, s1/3, -s1/3, zoo]
shifts = [(0,1)]*7 +[(0,s1/2)]*2 +[(0,s1/3)]*2 + [(1,0)]
generate_regions(slopes, shifts)
