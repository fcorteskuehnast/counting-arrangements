from sympy import Point, Polygon, Circle, Segment, convex_hull, Line, oo, zoo, sqrt, sympify, pi
import xml.etree.ElementTree as ET
from xml.dom import minidom
from pathlib import Path    
import pprint

OUTPUT_FOLDER = "output/"

def prettify(elem):
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def pt_in_poly(p,poly):
    if p in poly.vertices or any(p in s for s in poly.sides):
        return True
    return poly.encloses_point(p)

def intersect_polygons(poly1,poly2):
    
    if isinstance(poly1, Segment):
        pts1 = list(poly1.points)
    else:
        pts1 = list(poly1.vertices)
    if isinstance(poly2, Segment):
        pts2 = list(poly2.points)
    else:
        pts2 = list(poly2.vertices)
    
    
    pts = list(poly1.intersection(poly2)) + pts1 + pts2
    inside = []
    for p in pts:
        if pt_in_poly(p,poly1) and pt_in_poly(p,poly2):
            inside.append(p)
    return convex_hull(*inside)
    
    
def add_line(line, page, layer, w, h, r):
    
    if r is not None:
        boundary = Circle(Point(0, 0), r)
    else:
        boundary = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
    intersection = boundary.intersection(line)
    
    if len(intersection) == 2:
        p1,p2 = intersection
        p1 = ( float(50*(p1[0]+w/2)), float(50*(p1[1]+h/2)) )
        p2 = ( float(50*(p2[0]+w/2)), float(50*(p2[1]+h/2)) )
        path = ET.SubElement(page, "path")
        path.set("layer", layer)
        path.text = f"{p1[0]} {p1[1]} m {p2[0]} {p2[1]} l"
    

def add_polygon(poly, page, layer, w, h):
    square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
    
    poly = intersect_polygons(poly,square)

    path = ET.SubElement(page, "path")
    path.set("layer", layer)
    path.set("fill", "red")
    path.set("opacity","50%")
    
    txt = ""
    p = poly.vertices[0]
    txt += f"{float(50*(p[0]+w/2))} {float(50*(p[1]+h/2))} m "
    for p in poly.vertices[1:]:
        txt += f"{float(50*(p[0]+w/2))} {float(50*(p[1]+h/2))} l "
    txt += "h"
    
    path.text = txt


def make_IPE(bundles, regions, regions_names, filename, w=10, h=10, r=None):
    
    tree = ET.parse("template2.ipe")
    ipe = tree.getroot()
    
    ipestyle = ipe.find("ipestyle")
    
    layout = ET.SubElement(ipestyle, 'layout')
    layout.set("paper", f"{50*w} {50*h}")
    layout.set("origin", "0 0")
    layout.set("frame", f"{50*w} {50*h}")
    
    page = ET.SubElement(ipe, "page")
    
    square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
    
    
    for i in range(len(bundles)):
        layer = ET.SubElement(page, "layer")
        layer.set("name", f"bundle_{i+1}")
        
    for i in range(len(regions)):
        layer = ET.SubElement(page, "layer")
        layer.set("name", f"{regions_names[i]}")
        
    for i,bundle in enumerate(bundles):
        for line in bundle:
            add_line(line, page, f"bundle_{i+1}", w, h, r)
    
    for i,region in enumerate(regions):
        for poly in region:
            add_polygon(poly, page, f"{regions_names[i]}", w, h)
                
    xml_string = prettify(ipe)
    with open(filename, "w") as text_file:
        text_file.write(xml_string)

        
def make_config(k, region_names, areas, filename, precision = 0.3):
    
    config_filenames = {name:f"shapes/{name}/{name}.ipe" for name in region_names}
    config_areas = {name:area for name, area in zip(region_names,areas)}
    
    config = {'k':k, 'precision': precision, }
    config['k'] = k
    config['precision'] = 0.3
    config['filenames'] = {name:f"shapes/{name}/{name}.ipe" for name in region_names}
    config['LGV3'] = {'regions': [], 'patchsize_x' : 500, 'patchsize_y' : 500, 'log_count' : 349033}
    config['areas'] = {name:float(area) for name, area in zip(region_names,areas)}
    config['exact areas'] = {name:str(area) for name, area in zip(region_names,areas)}
    
    with open(filename, "w") as file:
        file.write(pprint.pformat(config, indent = 2, sort_dicts=False))
    
s1 = sympify(1)

def get_bundles(slopes, shifts, w,h):
    square = Polygon(Point(-w/2,-h/2), Point(w/2,-h/2), Point(w/2,h/2), Point(-w/2,h/2))
    
    bundles = []
    for i in range(len(slopes)):
        slope, shift = slopes[i], Point(shifts[i])
        bundle = []
        k = 1
        while True:
            n = k//2*(-1)**k
            k += 1
            
            if slope == zoo or slope == oo:
                line = Line(Point(0,-1000)+n*shift, Point(0,1000)+n*shift)
            else:
                line = Line(Point(-1000,-1000*slope)+n*shift, Point(1000,1000*slope)+n*shift)
            
            intersection = square.intersection(line)
            if len(intersection) == 2:
                bundle.append(line)
            else:
                break
            
        bundles.append(bundle)
    return bundles
    

def fingerprint(s):
    return str((s+1e-6).n())[:8]

def list_fingerprint(L):
    return tuple(fingerprint(x) for x in sorted(L))
    
def all_slope_symmetries(slopes,k):
    for r in range(k):
        sym = [tan(arctan(x)+2*pi*r/k) for x in slopes]
        sym_pos = [x if x!=zoo else oo for x in sym]
        sym_neg = [-x if x!=oo else oo for x in sym_pos]
        yield tuple(sorted(sym_pos))
        yield tuple(sorted(sym_neg))

        
def dist_line_poly(line, poly):
    return min(line.distance(x) for x in poly.vertices)

def get_widths(slopes, shifts):
    num_bundles = len(slopes)

    square = Polygon(Point(-100,-100), Point(100,-100), Point(100,100), Point(-100,100))
    upper_lines = []
    lower_lines = []
    central_region = square
    
    mid_lines = []
    for s in slopes:
        if s == zoo or s == oo:
            mid_lines.append(Line( (0,0), (0,1) ))
        else:
            mid_lines.append(Line( (0,0), (1,s) ))
    for i in range(num_bundles):
        p1,p2 = square.intersection(mid_lines[i])        
        shift = Point(shifts[i])
        upper_lines.append(Line(p1+shift/2, p2+shift/2))
        lower_lines.append(Line(p1-shift/2, p2-shift/2))
        central_region = intersect_polygons(central_region, convex_hull(p1+shift/2, p2+shift/2, p1-shift/2, p2-shift/2) )
    
    widths = []
    for i in range(num_bundles):
        upper = upper_lines[i]
        lower = lower_lines[i]
    
        du = dist_line_poly(upper, central_region)
        dl = dist_line_poly(lower, central_region)
        dul = upper.distance(lower.points[0])
        
        widths.append((dul - du - dl)/dul)
        
    return widths


def generate_regions(slopes, shifts, num_sym, compact=False):
    
    num_bundles = len(slopes)
    
    slope_shift = {}
    for i in range(num_bundles):
        slope_shift[fingerprint(slopes[i])] = shifts[i]
    
    slopes = sorted(slopes)
    shifts = [slope_shift[fingerprint(s)] for s in slopes]
    
    if compact:
        widths = get_widths(slopes, shifts)
        multiplier = num_bundles*(num_bundles-1)/sum(widths)/(sum(widths)-1)
        print(f"Note: the bound on log(B(n)) given by ipe2bound.py will have to be multiplied by {multiplier} ~ {multiplier.n():.{7}f}")
    else:
        widths = [1]*num_bundles

    mid_lines = []
    for s in slopes:
        if s == oo:
            mid_lines.append(Line( (0,0), (0,1) ))
        else:
            mid_lines.append(Line( (0,0), (1,s) ))

    square = Polygon(Point(-100,-100), Point(100,-100), Point(100,100), Point(-100,100))
    upper_lines = []
    lower_lines = []
    basic_regions = []
    for i in range(num_bundles): #line in enumerate(mid_lines):
        p1,p2 = square.intersection(mid_lines[i])        
        shift = widths[i]*Point(shifts[i])
        upper_lines.append(Line(p1+shift/2, p2+shift/2))
        lower_lines.append(Line(p1-shift/2, p2-shift/2))
        basic_regions.append( convex_hull(p1+shift/2, p2+shift/2, p1-shift/2, p2-shift/2) )

    print("Generating regions...")
    regions = [square]
    for i,line in enumerate(upper_lines + lower_lines):
        print(f"{i+1}/{len(upper_lines) + len(lower_lines)}")
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
        print(f"{i+1}/{len(regions)}")
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
            
        region_bundles = tuple(sorted(region_bundles))

        rep = None
        for sym in all_slope_symmetries2(region_bundles, num_sym):
            for x in agreggated_regions:
                if list_fingerprint(sym) == list_fingerprint(x): 
                    rep = x
                    break
            if rep is not None:
                break
        if rep is None:
            rep = region_bundles
            
        if rep in agreggated_regions:
            agreggated_regions[rep].append(regions[i])
            areas[rep] += a
        else:
            agreggated_regions[rep] = [regions[i]]
            areas[rep] = a     
    
    extremal_lines = []
    for i in range(num_bundles):
        extremal_lines.append([lower_lines[i],upper_lines[i]])
    
    
    sorted_reps = sorted(agreggated_regions, key=lambda x:-len(x))
    grouped_reps = [[]]
    for rep in sorted_reps:
        if len(grouped_reps[-1]) == 0 or len(rep) == len(grouped_reps[-1][0]):
            grouped_reps[-1].append(rep)
        else:
            grouped_reps.append([rep])
    
    regions = []
    region_names = []
    region_areas = []
    
    for i,group in enumerate(grouped_reps):
        for j,rep in enumerate(group):
            regions.append(agreggated_regions[rep])
            region_areas.append(areas[rep])
            if len(group) > 1:
                region_names.append(f"R{len(rep)}{'abcdefghijklmnopqrstuvwxyz'[j]}")
            else:
                region_names.append(f"R{len(rep)}")
    
    print("Making base construction Ipe and config file...")
    
    Path(OUTPUT_FOLDER+"regions/").mkdir(parents=True, exist_ok=True)
    Path(OUTPUT_FOLDER+"shapes/").mkdir(parents=True, exist_ok=True)

    make_IPE(extremal_lines, regions, region_names, f"{OUTPUT_FOLDER}base_construction.ipe", w=7, h=7, r=3)
    make_config(len(slopes), region_names, region_areas, f"{OUTPUT_FOLDER}k{len(slopes)}.config", 0.3)
    
    
    print("Making other Ipe files...")
    
    for i,rep in enumerate(sorted_reps):
        print(f"{i+1}/{len(agreggated_regions)}")
        name = region_names[i]
        
        Path(f"{OUTPUT_FOLDER}shapes/{name}/").mkdir(parents=True, exist_ok=True)
        Path(f"{OUTPUT_FOLDER}regions/{name}/").mkdir(parents=True, exist_ok=True)
        
        make_IPE(extremal_lines, [agreggated_regions[rep]], [name], f"{OUTPUT_FOLDER}regions/{name}/{name}.ipe", w=8, h=8)
        
        rep_shifts = [slope_shift[fingerprint(s)] for s in rep]
        grid_bundles = get_bundles(rep,rep_shifts, w=12, h=12)
        make_IPE(grid_bundles, [], [], f"{OUTPUT_FOLDER}shapes/{name}/{name}.ipe", w=12, h=12)
  



############################################################
############################################################


slopes = [0,1,2,3,-1,-2,-3, s1/2, -s1/2, s1/3, -s1/3, oo]
shifts = [(0,1)]*7 +[(0,s1/2)]*2 +[(0,s1/3)]*2 + [(1,0)]
rotational_symmetries = 8
compact = False
OUTPUT_FOLDER = "output/"

generate_regions(slopes, shifts, rotational_symmetries, compact=compact)


###
#slopes = [      0,              oo,  sqrt(3),  -sqrt(3),  1/sqrt(3),  -1/sqrt(3)]
#shifts = [(0,1/2),  (sqrt(3)/2, 0),    (0,1),     (0,1),      (0,1),       (0,1)]
#slopes = [0,1,2,3,-1,-2,-3, s1/2, -s1/2, s1/3, -s1/3, zoo]
#shifts = [(0,1),(0,1),(0,2),(0,3), (0,1), (0,2), (0,3),  (0,1),  (0,1),    (0,1),     (0,1),    (1,0)]
###
    

