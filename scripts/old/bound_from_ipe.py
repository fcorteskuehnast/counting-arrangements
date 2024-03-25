import ctypes
from fractions import Fraction
from math import log2
from operator import mul
from functools import reduce
import ast
import os
import argparse

#from ipe3bip.py import ipe2bip

num_prt_arr_char_p = ctypes.CDLL("../prt_arr_counting.so").num_prt_arr_char_p
num_prt_arr_char_p.restype = ctypes.c_char_p

parser = argparse.ArgumentParser()
parser.add_argument("config_file",type=str,help="config file")
#parser.add_argument("--redo_calculations","-r",default=False,type=bool,help="redo calculations")
args = parser.parse_args()
print("args",args)

def bip_from_file(path):
    with open(path, "r") as f:
        data = f.read()
        return [int(a) for a in data.split(" ")]
    
def crossings_from_file(path):
    with open(path, "r") as f:
        return ast.literal_eval(f.read())
    
def count_from_file(path):
    with open(path, "r") as f:
        return int(f.read())

def bip_cr_from_ipe(path):
    os.system(f"sage ipe2bip.sage {path}")

    #assert False, "TODO"
    bip = bip_from_file(path+".bip")
    cr = crossings_from_file(path+".crossings")

    return bip, cr

def count_from_bip(bip):
    size = len(bip)
    T = ctypes.c_char * size
    return int(num_prt_arr_char_p(T(*bip), size))

config = "6sl_small.config"
config = args.config_file

with open(config, 'r') as f:
    data = ast.literal_eval(f.read())
    k = data['k']
    filename_region = data['filename_region']
    num_of_cr_region = { i:Fraction(s) for (i,s) in data['num_of_cr_region'].items() }

'''
# bipermutation, number of crossings and number of partial arrangements should be in files w/ same basename and extensions .bip, .crossings, .count
filename_region = {
    3 : "hex_tilings/small_3sl.ipe",
    4 : "hex_tilings/small_4sl.ipe",
    5 : "hex_tilings/small_5sl.ipe",
    6 : "hex_tilings/small_6sl.ipe"
}
num_of_cr_region = {
    3 : Fraction(3, 4 * 6**2), 
    4 : Fraction(1, 4 * 6**2), 
    5 : Fraction(1, 4 * 6**2), 
    6 : Fraction(1, 4 * 6**2)
}
'''

bip_patch = {}
num_of_cr_patch = {}
num_of_prt_arr_patch = {}

for i, p in filename_region.items():
    ### get bips and counts ###
    try:
        bip_patch[i] = bip_from_file(p+".bip")
        temp = crossings_from_file(p+".crossings")
        num_of_cr_patch[i] = temp[i]
    except:
        bip_patch[i], temp = bip_cr_from_ipe(p)
        num_of_cr_patch[i] = temp[i]

    

    #### get numbers of partial arrangements ###
    try:
        num_of_prt_arr_patch[i] = count_from_file(p+".count")
    except:
        num_of_prt_arr_patch[i] = count_from_bip(bip_patch[i])
        with open(p+".count", "w") as f:
            f.write(str(num_of_prt_arr_patch[i]))

print("Numbers of partial arrangements:")
for (i, n) in num_of_prt_arr_patch.items():
    print(f"\tF(P{i}) = {n}")

### calculate bound ###
num_of_patches_region = {i : Fraction(cr_region, cr_patch) for ((i, cr_region), cr_patch) in zip(num_of_cr_region.items(), num_of_cr_patch.values())}

print(f"Numbers of patches:")
for (i, f) in num_of_patches_region.items():
    print(f"\tm{i}(P{i}, n) = {f} n² - O(n)")

lg_F_n_region = {i : l * log2(B) for ((i,l),B) in zip(num_of_patches_region.items(), num_of_prt_arr_patch.values())}
lg_T_n_region = {i : (k/(k-1)) * a for (i, a) in lg_F_n_region.items()}

print("Contributions by region:")
for (i, t) in lg_T_n_region.items():
    print(f"\tlog( F(P{i})^mi ) = {t:.4f}")

lg_F_n = log2(reduce(mul, [b**l for (b,l) in zip(num_of_prt_arr_patch.values(), num_of_patches_region.values())], 1))
lg_T_n = (k/(k-1)) * lg_F_n

print(f"Complete bound:\n\tlog(F(n)) = {lg_F_n}\n\tlog(B(n)) = {lg_T_n}")
