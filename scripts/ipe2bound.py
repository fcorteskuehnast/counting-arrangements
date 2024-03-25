from fractions import Fraction
from math import log2
from operator import mul
from functools import reduce
import ast
import os
import argparse
import numpy as np
import time
import multiprocessing

import os, sys
script_path = os.path.dirname(__file__)
cpp_exe = script_path+"/../cpp/prt_arr_counting"
if not os.path.isfile(cpp_exe):
    print("please compile c++ program first; checkout the README file for more information")
    exit()


import hashlib

# https://stackoverflow.com/a/44873382
'''
def sha256sum(filename):
    with open(filename, 'rb', buffering=0) as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
'''

def sha256sum(filename):
    h  = hashlib.sha256()
    b  = bytearray(128*1024)
    mv = memoryview(b)
    with open(filename, 'rb', buffering=0) as f:
        while n := f.readinto(mv):
            h.update(mv[:n])
    return h.hexdigest()

parser = argparse.ArgumentParser()
parser.add_argument("config_file",type=str,help="config file")
parser.add_argument("--only",type=str,help="only compute the specified patches/regions")
parser.add_argument("--threads",type=int,default=1,help="number of threads (for parallization)")
parser.add_argument("--max_memory",type=int,default=0,help="maximum estimated memory usage (in MB)")
parser.add_argument("--debug",type=int,default=0,help="debug level")
parser.add_argument("--precision",type=float,default=None,help="precision")
parser.add_argument("--load",action='store_false',help="by default, load previously computed values (use parameter to disable)")

parser.add_argument("--visualize","-v",action='store_true',help="visualize shapes")
parser.add_argument("--latex_table",action='store_true',help="print table in latex format")
#parser.add_argument("--use-known-values","-k",action='store_false',help="by default use known values for up to 16")
parser.add_argument("--incomplete",action='store_true',help="disable completeness checks in config file")

parser.add_argument("--ratio",action='store_true',help="compute contribution/time ratio")

parser.add_argument("--replace_region",type=str,help="replace the specified regions")
parser.add_argument("--replace_filename",type=str,help="replace filename for the specified regions")

#parser.add_argument("--redo_calculations","-r",default=False,type=bool,help="redo calculations")
args = parser.parse_args()
print("args",args)

if  args.threads <= 0:
    args.threads = 2*(1+multiprocessing.cpu_count())

assert(args.max_memory == 0 or args.max_memory > 70 + 14*args.threads) #Memory limit too small (increase limit or decrease threads)


def compute_LGV(gridsize):
    path = f"lgv_gridsize{gridsize}.txt"
    print(f"compute LGV values for gridsize={gridsize} via lgv.sage:")
    cmd = f"sage {script_path}/lgv.sage {gridsize} -o {path}"
    assert(os.system(cmd) == 0)
    with open(path) as f:
        return ast.literal_eval(f.read())


def make_cr_bip_files(shape_path):
    new_checksum = sha256sum(shape_path) 
    checksum_path = shape_path+".checksum"
    try:
        old_checksum = open(checksum_path).readline().replace("\n","")  
    except FileNotFoundError:
        old_checksum = None

    if not args.load or new_checksum != old_checksum:
        print("compute bip via ipe2bip.sage:")
        cmd = f"sage {script_path}/ipe2bip.sage {shape_path}"
        if args.precision: cmd += f" --precision {args.precision}"
        if args.visualize: cmd += " --visualize"
        if args.debug: cmd += f" --debug {args.debug}"
        assert(os.system(cmd) == 0)
        open(checksum_path,"w").write(new_checksum)
    else:
        print("ipe file did not change")


def make_count_file(bip_path, outpath, threads, max_memory, debug):
    if debug >= 1:
        print("call prt_arr_counting:") #TODO: maybe remove
    cmd = f"{cpp_exe} {bip_path} -o {outpath} -t {threads} -m {max_memory} -d {debug}"
    print("run command:",cmd)
    assert(os.system(cmd) == 0)
    
def object_from_file(path):
    if os.path.isfile(path):
        with open(path) as f:
            return ast.literal_eval(f.read())

def area_from_crs(crs):
    return crs[max(crs)] * k**2 # ... /n²

def map_vals(f, d):
    return dict({ key : f(val) for (key, val) in d.items()})

def zip_vals(d1, d2):
    assert set(d1.keys()) == set(d2.keys())
    return { key : (val1, val2) for ((key,val1),val2) in zip(d1.items(), d2.values())}

def intersect_keys(d1, d2):
    return ( {key : val for (key,val) in d1.items() if key in d2.keys()}, {key : val for (key,val) in d2.items() if key in d1.keys()} )


config = args.config_file
directory = os.path.dirname(config)
#print(directory)

assert(os.path.isfile(config))
data = object_from_file(config)
k = data['k']

if args.replace_filename:
    assert(args.replace_region)
if args.replace_region:
    assert(args.replace_region in data['filenames'])
    assert(args.replace_filename)
    data['filenames'][args.replace_region] = args.replace_filename
    print(f"*** replacing file for region {args.replace_region}: {args.replace_filename}")

if 'LGV3' in data: # use Lindström–Gessel–Viennot lemma for regions with 3 slopes
    LGV3 = data['LGV3'] 
    LGV3_regions = LGV3['regions']
else:
    LGV3 = None
    LGV3_regions = []

if 'precision' in data:
    assert(not args.precision)
    args.precision = data['precision']
IPE_regions = map_vals(lambda s : os.path.join(directory, s) , data['filenames'])
area_region = map_vals(Fraction, data['areas'])
#area_region, IPE_regions = intersect_keys(area_region, IPE_regions)


bip = {}
num_of_crossings = {}
num_of_prt_arr = {}
computing_time = {}


for p in area_region:
    if not args.incomplete:
        assert((p in IPE_regions) != (p in LGV3_regions)) # each region must either come an IPE file for dynamic programming or must be specified for the Lindström–Gessel–Viennot lemma
    else:
        assert(not (p in IPE_regions) or not (p in LGV3_regions)) # not both is given IPE file for dynamic programming or must be specified for the Lindström–Gessel–Viennot lemma


print()
print("============== DYNAMIC PROGRAMMING ==============")
for i,p in IPE_regions.items():
    if args.only and i != args.only: continue # if specified, only consider this patch
    print()
    print(80*"#")
    print(f"patch P_{i} @ {p}")
    print(80*"#")

    assert(os.path.isfile(p)) # ipe file exists

    print(f"compute and store bip, cr from ipe file {p}")
    print()
    make_cr_bip_files(p)

    num_file = p+".polygons"
    print()
    print("read number of polygons from",num_file)
    num_polygons = int(open(num_file).readline())
    assert(num_polygons >= 1)

    num_of_prt_arr[i] = 1 # []
    computing_time[i] = 0 # []
    num_of_prt_arr_list = []
    computing_time_list = []
    num_of_crossings[i] = {}

    for poly_index in range(num_polygons):
        ### get bips and counts ###
        bip_path = f"{p}_poly{poly_index}.bip"
        cr_path  = f"{p}_poly{poly_index}.crossings"

        
        print(f"\n### shape {poly_index}")
        print(f"load crossing data from file {cr_path}")
        crossings = object_from_file(cr_path)
        assert(crossings)
        #print("crossings",crossings)
        for (key,val) in crossings.items():
            if key not in num_of_crossings[i]: 
                num_of_crossings[i][key] = 0
            num_of_crossings[i][key] += val

        bip = open(bip_path).readline().replace("\n","")
        if args.debug: print(f"bip {bip}")

        summary_path = f"{bip_path}.summary"
        summary = object_from_file(summary_path)
        if not args.load or not summary or summary['bip'] != bip: 
            print(f"compute counts from scratch")
            make_count_file(bip_path, summary_path, args.threads, args.max_memory, args.debug)
            print(f"store in file {summary_path}")

        print(f"load counts from file {summary_path}")

        summary = object_from_file(summary_path)
        assert(summary) # must have been created
        assert(summary['bip'] == bip)
        if 'time' in summary:
            summary['real_time'] = summary['cpu_time'] = summary['time']
        time = summary['cpu_time']
        count = summary['count']
        #results = summary['results'] # only for statistical purposes

        num_of_prt_arr[i] *= count
        computing_time[i] += time 
        num_of_prt_arr_list.append(count)
        computing_time_list.append(time)

    print()
    if len(num_of_prt_arr_list) > 1:
        print(f"F(P_{i}) = {' * '.join(str(x) for x in num_of_prt_arr_list)} = {num_of_prt_arr[i]} (computing time: {' + '.join(str(x) for x in computing_time_list)} = {computing_time[i]}s)")
    else:
        print(f"F(P_{i}) = {num_of_prt_arr[i]} (computing time: {computing_time[i]}s)")

    print(f"num_of_crossings[{i}] = {num_of_crossings[i]}")



if LGV3:
    print()
    print("============== LGV3 ==============")
    patchsize = LGV3['patchsize']
    if args.load and 'log_count' in LGV3:
        logcount = LGV3['log_count']
        print(f"loaded logcount {logcount}")
        count = 2**logcount 
        cpu_time = 0
    else:
        results = compute_LGV(patchsize)
        assert(results['patchsize'] == patchsize)
        count = results['count']
        cpu_time = results['cpu_time']

    num_crossings = patchsize*patchsize
    print(f"loaded num_crossings {num_crossings}")

    for i in LGV3_regions:
        if args.only and i != args.only: continue # if specified, only consider this patch
        print(f"setting region {i}")
        num_of_prt_arr[i] = count
        computing_time[i] = cpu_time
        num_of_crossings[i] = {3: num_crossings}


print()
print("============== SUMMARY ==============")
for i in num_of_prt_arr:
    if log2(num_of_prt_arr[i]) < 100:
        print(f"\t{f'F(P_{i})':<12} = {num_of_prt_arr[i]:>42} (computing time: {computing_time[i]:.2f}s)")
    else:
        print(f"\t{f'log(F(P_{i}))':<12} = {log2(num_of_prt_arr[i]):>42} (computing time: {computing_time[i]:.2f}s)")


### calculate bound ###
num_of_patches_region = {}

for i in num_of_crossings:
    num_of_patches_region[i] = Fraction(area_region[i], area_from_crs(num_of_crossings[i]))

print(f"numbers of patches:")
for (i, f) in num_of_patches_region.items():
    print(f"\t{f'µ_{i}(P_{i}, n)':<18} = {f'{f} n² - O(n)':<21} = {f'{f*k**2} m² - O(m)':<18}")

lg_F_n_region = map_vals(lambda x : x[0]*log2(x[1]), zip_vals(num_of_patches_region, num_of_prt_arr))
lg_T_n_region = map_vals(lambda a : a * k/(k-1) , lg_F_n_region)

print("contributions by region:")
for (i, t) in lg_F_n_region.items():
    print(f"\tc_{i} \t≈ {t:.5f}")

#lg_F_n = log2(reduce(mul, [b**l for (b,l) in zip(num_of_prt_arr.values(), num_ofes_region.values())], 1))
lg_F_n = sum(log2(b)*l for (b,l) in zip(num_of_prt_arr.values(), num_of_patches_region.values()))
lg_T_n = (k/(k-1)) * lg_F_n

print(f"complete bound:")
print(f"\t{f'log(F_{k}(n))':<13}> {lg_F_n}")
print(f"\t{'log(B(n))':<13}> {lg_T_n}")


print()
print("============== TABLE ==============")

rows = {}
for i in num_of_prt_arr:
    rows[i] = {}
    rows[i]['region'] = i
    if i in LGV3_regions: rows[i]['region'] += "*"
    rows[i]['log2(F)'] = f"{log2(num_of_prt_arr[i]):.2f}"
    rows[i]['lim µ/n²'] = f"{num_of_patches_region[i]}"
    rows[i]['contribution'] = f"{lg_F_n_region[i]:.5f}"
    rows[i]['CPU time'] = f"{computing_time[i]:.2f}s"
    if args.ratio: rows[i]['contribution/time'] = f"{lg_F_n_region[i]/computing_time[i]:.9f}" if computing_time[i] else '???'


summary = {}
summary['region'] = '\sum' if args.latex_table else 'Σ'
summary['log2(F)'] = '-'
summary['lim µ/n²'] = '-'
summary['contribution'] = f"{sum(lg_F_n_region[i] for i in num_of_prt_arr):.5f}"
summary['CPU time'] = f"{sum(computing_time[i] for i in num_of_prt_arr):.2f}s"
if args.ratio: summary['contribution/time'] = '-'


columns = rows[min(rows)].keys()

num_of_columns = len(columns)
col_widths = {col_name: 2+max([len(col_name)] + [len(e[col_name]) for e in rows.values()]) for col_name in columns}


if args.latex_table:
    endc = "&"
    endl = "\\\\\n"
else:
    endc = ""
    endl = "\n"


line_format = (endc.join(["{:>%s}"]*num_of_columns)) % tuple(col_widths.values())+endl


header = line_format.format(*columns) # header
if args.latex_table:
    separating_line = "\\hline \n" 
else:
    separating_line = len(header)*"-"+"\n" 

table = header
table += separating_line
for row in rows.values():
    table += line_format.format(*tuple(row.values())) # row
table += separating_line
table += line_format.format(*summary.values()) # summary
print(table)

