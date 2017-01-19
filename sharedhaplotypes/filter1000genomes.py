#!/usr/bin/env python
import sys
data=[a.strip().split(",") for a in open("HMERF_Haplotype.csv")]

header,data = data[0],data[1:]
outHMERF = open("HMERFmatch1K.csv","w")
outHMERF.write(",".join(h for h in header)+"\n")
## lets step through the positions available 
positions=dict()

for index,b in enumerate(data):
    positions[int(b[0])] = (index,b)

infile = open("impute1K.impute.legend")
header=infile.readline()
outlegend = open("match1K.legend","w")
outlegend.write(header)

seen=set()
indices=dict()
print "=============================\nSNPs present in both"
print "-----------------------------"
for linenum,a in enumerate(infile):
    ID,pos,allele0,allele1 = a.strip().split()
   
    if int(pos) in positions:
        outlegend.write(a)
        print ID,pos,allele0,allele1,"| HMERF: "," ".join(v for v in positions[int(pos)][1])
        seen.add(int(pos))
        indices[linenum] = ([allele0,allele1],positions[int(pos)][0])
   

infile.close()
outlegend.close()

missing=[]
print "=============================\nMissing in 1K"
print "-----------------------------"

for k,v in positions.iteritems():
    if k not in seen:
        missing.append(positions[k])
        print " ".join(b for b in v[1])

seen.add(179410829)
pos = [s for s in seen]
pos.sort()
print pos

for p in pos:
    if p in seen:
        outHMERF.write(",".join(v for v in positions[p][1])+"\n")

outHMERF.close()



missing=sorted(missing)
print "============================="
# now lets try to filter the data
infile = open("impute1K.impute.hap")
outfile = open("match1K","w")
for index,line in enumerate(infile):
    if index in indices:
        b = line.strip().split()
        code =   indices[index][0]
        outfile.write(" ".join(code[int(a)] for a in b)+"\n")

outfile.close()
infile.close()
