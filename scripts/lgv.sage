#!/usr/bin/python
import argparse
import datetime

parser = argparse.ArgumentParser()
parser.add_argument("patchsize",type=int,help="patchsize")
parser.add_argument("--outfile","-o",type=str,help="output file")
args = parser.parse_args()

time_start = datetime.datetime.now()

l = args.patchsize

M = Matrix([[binomial(2*l-abs(l-i)-abs(l-j),(2*l-abs(l-i)-abs(l-j)+3*abs(i-j))/2) for j in [1..2*l-1]] for i in [1..2*l-1]])
d = M.det()
c = RR(log(d)/log(2))
#print(M)
end_time = datetime.datetime.now()
time_diff = (end_time-time_start).total_seconds()


print(f"l={l}\tlog2(det)={c}\tratio={c/l^2}\tcpu_time={time_diff}")

if args.outfile:
	print(f"write result to file {args.outfile}")
	with open(args.outfile,"w") as f:
		f.write(str({'patchsize':l,'count':d,'cpu_time':time_diff})+"\n")
		