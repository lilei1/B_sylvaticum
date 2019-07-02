#!/usr/bin/env python

#Written by Li Lei 2017/12/12 
#required python3
#This is for calculate the indel length for each locus from vcf file called with GATK

import sys
print ("\t".join(["chr","length","class"])) #print header

with open(sys.argv[1]) as f:
	for line in f:
		if line.startswith('##'):
			continue
		else:
			tmp = line.strip().split('\t')
			cate = tmp[8].strip().split(';')
			cate1 = cate[3].strip().split('=')
			cls = cate1[1]
			dist = abs(int(tmp[4]) - int(tmp[3]))
			print ("\t".join([str(tmp[0]),str(dist),str(cls)]))