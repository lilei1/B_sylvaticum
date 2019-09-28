#!/usr/bin/env python
"""This is to convert the Annotation file into association file for GoaTools to do GO enrichment analysis."""

#   This script was written by Li Lei, 09-29-2019, Walnut Creek


"""
For parsing goatools output into a REViGO-compatible two-column format.
Run on command line like so:

python parseFuncOut.py <inFile> <outFile>

e.g.
python parseFuncOut.py downreg_Bstacei_vs_Bhybridum_syntOrthos_lengthNorm_sigOnly.out downreg_Bstacei_vs_Bhybridum_syntOrthos_lengthNorm_sigOnly_outList.bonferroniPvals
"""
import sys

inputList = sys.argv
annots = {}

with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                GID = tmp[2]
                #print (GID)
                GOids = []
                for e in tmp:
                    if e.startswith('GO:'):
                        L = e.split(',')
                        for x in L:
                            GOids.append(x.strip('GO:'))
                if GOids:
                        annots[GID] = ['GO:%s'%(ID) for ID in GOids]           
    
    
outFile = open('Bsylvaticum_490_v1.1.association', 'w')

for k in annots.keys():
    outFile.write('%s\t%s\n'%(k, ';'.join(annots[k])))
    
    
outFile.close()
