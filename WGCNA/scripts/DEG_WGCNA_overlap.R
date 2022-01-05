library(GeneOverlap)
#read module files
turquoise <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-turquoise.txt", header = T, 
                        stringsAsFactors = FALSE)
head(turquoise)
brown <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-brown.txt", header = T, 
                    stringsAsFactors = FALSE)
head(brown)
pink <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-pink.txt", header = T, 
                   stringsAsFactors = FALSE)
head(pink)
red <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-red.txt", header = T, 
                  stringsAsFactors = FALSE)
head(red)
blue <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-blue.txt", header = T, 
                   stringsAsFactors = FALSE)
head(blue)

yellow <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-yellow.txt", header = T, 
                     stringsAsFactors = FALSE)
head(yellow)

black <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-black.txt", header = T, 
                    stringsAsFactors = FALSE)
head(black)

green <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-green.txt", header = T, 
                    stringsAsFactors = FALSE)
head(green)

grey <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-grey.txt", header = T, 
                   stringsAsFactors = FALSE)

###Read the DGE data
drought1h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought1h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
drought2h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought2h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
drought5h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought5h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
drought10h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought10h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
drought24h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought24h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
rownames(drought24h)
#heat:
heat1h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat1h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
heat2h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat2h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
heat5h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat5h_DEG.txt", header = T, 
                        stringsAsFactors = FALSE)
heat10h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat10h_DEG.txt", header = T, 
                         stringsAsFactors = FALSE)
heat24h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat24h_DEG.txt", header = T, 
                         stringsAsFactors = FALSE)
#salt:
salt1h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt1h_DEG.txt", header = T, 
                     stringsAsFactors = FALSE)
salt2h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt2h_DEG.txt", header = T, 
                     stringsAsFactors = FALSE)
salt5h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt5h_DEG.txt", header = T, 
                     stringsAsFactors = FALSE)
salt10h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt10h_DEG.txt", header = T, 
                      stringsAsFactors = FALSE)
salt24h <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt24h_DEG.txt", header = T, 
                      stringsAsFactors = FALSE)
head(green)
#green module vs all of the DGE from each treatment
###drought1h
green.drought1h <- intersect(rownames(drought1h),green$nodeName)
green.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                         green$nodeName,
                                            genome.size=36927)
green.go.obj.drought1h <- testGeneOverlap(green.go.obj.drought1h)
print(green.go.obj.drought1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=187, e.g. Brasy1G033300.v1.1 Brasy1G074800.v1.1 Brasy1G075300.v1.1
#Union size=4721, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32206 3569
#inB    965  187
#Overlapping p-value=7.3e-11
#Odds ratio=1.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###drought2h
green.drought2h <- intersect(rownames(drought2h),green$nodeName)
green.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                         green$nodeName,
                                         genome.size=36927)
green.go.obj.drought2h <- testGeneOverlap(green.go.obj.drought2h)
print(green.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=299, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G033300.v1.1
#Union size=6850, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30077 5698
#inB    853  299
#Overlapping p-value=9.9e-18
#Odds ratio=1.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought5h
green.drought5h <- intersect(rownames(drought5h),green$nodeName)
green.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                         green$nodeName,
                                         genome.size=36927)
green.go.obj.drought5h <- testGeneOverlap(green.go.obj.drought5h)
print(green.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=580, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Union size=10694, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 26233 9542
#inB    572  580
#Overlapping p-value=4.2e-63
#Odds ratio=2.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###drought10h
green.drought10h <- intersect(rownames(drought10h),green$nodeName)
green.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                          green$nodeName,
                                          genome.size=36927)
green.go.obj.drought10h <- testGeneOverlap(green.go.obj.drought10h)
print(green.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
#  listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=609, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G011100.v1.1
#Union size=10926, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 26001 9774
#inB    543  609
#Overlapping p-value=4.4e-72
#Odds ratio=3.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###drought24h
green.drought24h <- intersect(rownames(drought24h),green$nodeName)
green.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                         green$nodeName,
                                         genome.size=36927)
green.go.obj.drought24h <- testGeneOverlap(green.go.obj.drought24h)
print(green.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=659, e.g. Brasy1G003600.v1.1 Brasy1G010400.v1.1 Brasy1G011100.v1.1
#Union size=12047, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 24880 10895
#inB    493   659
#Overlapping p-value=5.2e-76
#Odds ratio=3.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat1h
green.heat1h <- intersect(rownames(heat1h),green$nodeName)
green.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                      green$nodeName,
                                      genome.size=36927)
green.go.obj.heat1h <- testGeneOverlap(green.go.obj.heat1h)
print(green.go.obj.heat1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=946, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Union size=6605, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30322 5453
#inB    206  946
#Overlapping p-value=0e+00
#Odds ratio=25.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###heat2h
green.heat2h <- intersect(rownames(heat2h),green$nodeName)
green.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                      green$nodeName,
                                      genome.size=36927)
green.go.obj.heat2h <- testGeneOverlap(green.go.obj.heat2h)
print(green.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=835, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Union size=5997, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30930 4845
#inB    317  835
#Overlapping p-value=0e+00
#Odds ratio=16.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat5h
green.heat5h <- intersect(rownames(heat5h),green$nodeName)
green.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                      green$nodeName,
                                      genome.size=36927)
green.go.obj.heat5h <- testGeneOverlap(green.go.obj.heat5h)
print(green.go.obj.heat5h)


###heat10h
green.heat10h <- intersect(rownames(heat10h),green$nodeName)
green.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                       green$nodeName,
                                       genome.size=36927)
green.go.obj.heat10h <- testGeneOverlap(green.go.obj.heat10h)
print(green.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=477, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G011100.v1.1
#Union size=4324, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32603 3172
#inB    675  477
#Overlapping p-value=3.8e-182
#Odds ratio=7.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###heat24h
green.heat24h <- intersect(rownames(heat24h),green$nodeName)
green.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                       green$nodeName,
                                       genome.size=36927)
green.go.obj.heat24h <- testGeneOverlap(green.go.obj.heat24h)
print(green.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=659, e.g. Brasy1G003600.v1.1 Brasy1G010400.v1.1 Brasy1G011100.v1.1
#Union size=12047, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 24880 10895
#inB    493   659
#Overlapping p-value=5.2e-76
#Odds ratio=3.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###salt1h
green.salt1h <- intersect(rownames(salt1h),green$nodeName)
green.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                      green$nodeName,
                                      genome.size=36927)
green.go.obj.salt1h <- testGeneOverlap(green.go.obj.salt1h)
print(green.go.obj.salt1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=187, e.g. Brasy1G033300.v1.1 Brasy1G074800.v1.1 Brasy1G075300.v1.1
#Union size=4721, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32206 3569
#inB    965  187
#Overlapping p-value=7.3e-11
#Odds ratio=1.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###salt2h
green.salt2h <- intersect(rownames(salt2h),green$nodeName)
green.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                      green$nodeName,
                                      genome.size=36927)
green.go.obj.salt2h <- testGeneOverlap(green.go.obj.salt2h)
print(green.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=52, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G071700.v1.1
#Union size=1830, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 35097 678
#inB   1100  52
#Overlapping p-value=4e-08
#Odds ratio=2.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
green.salt5h <- intersect(rownames(salt5h),green$nodeName)
green.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                      green$nodeName,
                                      genome.size=36927)
green.go.obj.salt5h <- testGeneOverlap(green.go.obj.salt5h)
print(green.go.obj.salt5h)
#Detailed information about this GeneOverlap object:
#  listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=17, e.g. Brasy1G080300.v1.1 Brasy1G434000.v1.1 Brasy1G515500.v1.1
#Union size=1507, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 35420 355
#inB   1135  17
#Overlapping p-value=0.077
#Odds ratio=1.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
green.salt10h <- intersect(rownames(salt10h),green$nodeName)
green.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                       green$nodeName,
                                       genome.size=36927)
green.go.obj.salt10h <- testGeneOverlap(green.go.obj.salt10h)
print(green.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=53, e.g. Brasy1G080300.v1.1 Brasy1G214000.v1.1 Brasy1G290300.v1.1
#Union size=1882, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 35045 730
#inB   1099  53
#Overlapping p-value=1.5e-07
#Odds ratio=2.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
green.salt24h <- intersect(rownames(salt24h),green$nodeName)
green.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                       green$nodeName,
                                       genome.size=36927)
green.go.obj.salt24h <- testGeneOverlap(green.go.obj.salt24h)
print(green.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=88, e.g. Brasy1G033300.v1.1 Brasy1G080300.v1.1 Brasy1G127300.v1.1
#Union size=2490, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 34437 1338
#inB   1064   88
#Overlapping p-value=1.2e-09
#Odds ratio=2.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

#brown module vs all of the DGE from each treatment
###drought1h
brown.drought1h <- intersect(rownames(drought1h),brown$nodeName)
brown.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                         brown$nodeName,
                                         genome.size=36927)
brown.go.obj.drought1h <- testGeneOverlap(brown.go.obj.drought1h)
print(brown.go.obj.drought1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=184, e.g. Brasy1G028800.v1.1 Brasy1G052900.v1.1 Brasy1G053000.v1.1
#Union size=6137, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30790 3572
#inB   2381  184
#Overlapping p-value=1
#Odds ratio=0.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

########
###drought2h
brown.drought2h <- intersect(rownames(drought2h),brown$nodeName)
brown.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                         brown$nodeName,
                                         genome.size=36927)
brown.go.obj.drought2h <- testGeneOverlap(brown.go.obj.drought2h)
print(brown.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=329, e.g. Brasy1G028400.v1.1 Brasy1G028800.v1.1 Brasy1G037000.v1.1
#Union size=8233, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 28694 5668
#inB   2236  329
#Overlapping p-value=1
#Odds ratio=0.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###drought5h
brown.drought5h <- intersect(rownames(drought5h),brown$nodeName)
brown.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                         brown$nodeName,
                                         genome.size=36927)
brown.go.obj.drought5h <- testGeneOverlap(brown.go.obj.drought5h)
print(brown.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=897, e.g. Brasy1G005700.v1.1 Brasy1G007600.v1.1 Brasy1G009500.v1.1
#Union size=11790, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 25137 9225
#inB   1668  897
#Overlapping p-value=2.1e-18
#Odds ratio=1.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###drought10h
brown.drought10h <- intersect(rownames(drought10h),brown$nodeName)
brown.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                          brown$nodeName,
                                          genome.size=36927)
brown.go.obj.drought10h <- testGeneOverlap(brown.go.obj.drought10h)
print(brown.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
#  listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=879, e.g. Brasy1G005700.v1.1 Brasy1G006000.v1.1 Brasy1G007600.v1.1
#Union size=12069, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 24858 9504
#inB   1686  879
#Overlapping p-value=1e-12
#Odds ratio=1.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###drought24h
brown.drought24h <- intersect(rownames(drought24h),brown$nodeName)
brown.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                          brown$nodeName,
                                          genome.size=36927)
brown.go.obj.drought24h <- testGeneOverlap(brown.go.obj.drought24h)
print(brown.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=659, e.g. Brasy1G003600.v1.1 Brasy1G010400.v1.1 Brasy1G011100.v1.1
#Union size=12047, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 24880 10895
#inB    493   659
#Overlapping p-value=5.2e-76
#Odds ratio=3.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat1h
brown.heat1h <- intersect(rownames(heat1h),brown$nodeName)
brown.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                      brown$nodeName,
                                      genome.size=36927)
brown.go.obj.heat1h <- testGeneOverlap(brown.go.obj.heat1h)
print(brown.go.obj.heat1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=1223, e.g. Brasy1G005700.v1.1 Brasy1G008100.v1.1 Brasy1G009900.v1.1
#Union size=7741, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 29186 5176
#inB   1342 1223
#Overlapping p-value=4.5e-301
#Odds ratio=5.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###heat2h
brown.heat2h <- intersect(rownames(heat2h),brown$nodeName)
brown.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                      brown$nodeName,
                                      genome.size=36927)
brown.go.obj.heat2h <- testGeneOverlap(brown.go.obj.heat2h)
print(brown.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=1303, e.g. Brasy1G008100.v1.1 Brasy1G009500.v1.1 Brasy1G009900.v1.1
#Union size=6942, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 29985 4377
#inB   1262 1303
#Overlapping p-value=0e+00
#Odds ratio=7.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###heat5h
brown.heat5h <- intersect(rownames(heat5h),brown$nodeName)
brown.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                      brown$nodeName,
                                      genome.size=36927)
brown.go.obj.heat5h <- testGeneOverlap(brown.go.obj.heat5h)
print(brown.go.obj.heat5h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=910, e.g. Brasy1G008100.v1.1 Brasy1G009500.v1.1 Brasy1G009900.v1.1
#Union size=5304, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31623 2739
#inB   1655  910
#Overlapping p-value=6.2e-300
#Odds ratio=6.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###heat10h
brown.heat10h <- intersect(rownames(heat10h),brown$nodeName)
brown.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                       brown$nodeName,
                                       genome.size=36927)
brown.go.obj.heat10h <- testGeneOverlap(brown.go.obj.heat10h)
print(brown.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=4225, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=831, e.g. Brasy1G005800.v1.1 Brasy1G008100.v1.1 Brasy1G009500.v1.1
#Union size=5959, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30968 3394
#inB   1734  831
#Overlapping p-value=2e-193
#Odds ratio=4.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat24h
brown.heat24h <- intersect(rownames(heat24h),brown$nodeName)
brown.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                       brown$nodeName,
                                       genome.size=36927)
brown.go.obj.heat24h <- testGeneOverlap(brown.go.obj.heat24h)
print(brown.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=5793, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=1150, e.g. Brasy1G005800.v1.1 Brasy1G006000.v1.1 Brasy1G007600.v1.1
#Union size=7208, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 29719 4643
#inB   1415 1150
#Overlapping p-value=5.8e-294
#Odds ratio=5.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###salt1h
brown.salt1h <- intersect(rownames(salt1h),brown$nodeName)
brown.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                      brown$nodeName,
                                      genome.size=36927)
brown.go.obj.salt1h <- testGeneOverlap(brown.go.obj.salt1h)
print(brown.go.obj.salt1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=1648, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=47, e.g. Brasy1G035300.v1.1 Brasy1G257300.v1.1 Brasy1G440200.v1.1
#Union size=4166, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32761 1601
#inB   2518   47
#Overlapping p-value=1
#Odds ratio=0.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt2h
brown.salt2h <- intersect(rownames(salt2h),brown$nodeName)
brown.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                      brown$nodeName,
                                      genome.size=36927)
brown.go.obj.salt2h <- testGeneOverlap(brown.go.obj.salt2h)
print(brown.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=51, e.g. Brasy1G035300.v1.1 Brasy1G166100.v1.1 Brasy1G168700.v1.1
#Union size=3244, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 33683 679
#inB   2514  51
#Overlapping p-value=0.5
#Odds ratio=1.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

brown.salt5h <- intersect(rownames(salt5h),brown$nodeName)
brown.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                      brown$nodeName,
                                      genome.size=36927)
brown.go.obj.salt5h <- testGeneOverlap(brown.go.obj.salt5h)
print(brown.go.obj.salt5h)
#Detailed information about this GeneOverlap object:
#  listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=40, e.g. Brasy1G035300.v1.1 Brasy1G101000.v1.1 Brasy1G124200.v1.1
#Union size=2897, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 34030 332
#inB   2525  40
#Overlapping p-value=4.1e-03
#Odds ratio=1.6
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
brown.salt10h <- intersect(rownames(salt10h),brown$nodeName)
brown.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                       brown$nodeName,
                                       genome.size=36927)
brown.go.obj.salt10h <- testGeneOverlap(brown.go.obj.salt10h)
print(brown.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=87, e.g. Brasy1G076600.v1.1 Brasy1G101000.v1.1 Brasy1G124200.v1.1
#Union size=3261, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 33666 696
#inB   2478  87
#Overlapping p-value=1.1e-05
#Odds ratio=1.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
brown.salt24h <- intersect(rownames(salt24h),brown$nodeName)
brown.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                       brown$nodeName,
                                       genome.size=36927)
brown.go.obj.salt24h <- testGeneOverlap(brown.go.obj.salt24h)
print(brown.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=2565, e.g. Brasy1G005700.v1.1 Brasy1G005800.v1.1 Brasy1G006000.v1.1
#Intersection size=134, e.g. Brasy1G010000.v1.1 Brasy1G010600.v1.1 Brasy1G076600.v1.1
#Union size=3857, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33070 1292
#inB   2431  134
#Overlapping p-value=2.3e-04
#Odds ratio=1.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

#turquoise module vs all of the DGE from each treatment
###drought1h
turquoise.drought1h <- intersect(rownames(drought1h),turquoise$nodeName)
turquoise.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                             turquoise$nodeName,
                                             genome.size=36927)
turquoise.go.obj.drought1h <- testGeneOverlap(turquoise.go.obj.drought1h)
print(turquoise.go.obj.drought1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=2200, e.g. Brasy1G002500.v1.1 Brasy1G004300.v1.1 Brasy1G004400.v1.1
#Union size=8221, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 28706 1556
#inB   4465 2200
#Overlapping p-value=0e+00
#Odds ratio=9.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.3
########
###drought2h
turquoise.drought2h <- intersect(rownames(drought2h),turquoise$nodeName)
turquoise.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                             turquoise$nodeName,
                                             genome.size=36927)
turquoise.go.obj.drought2h <- testGeneOverlap(turquoise.go.obj.drought2h)
print(turquoise.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=3456, e.g. Brasy1G000500.v1.1 Brasy1G002500.v1.1 Brasy1G004400.v1.1
#Union size=9206, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 27721 2541
#inB   3209 3456
#Overlapping p-value=0e+00
#Odds ratio=11.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.4
###drought5h
turquoise.drought5h <- intersect(rownames(drought5h),turquoise$nodeName)
turquoise.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                             turquoise$nodeName,
                                             genome.size=36927)
turquoise.go.obj.drought5h <- testGeneOverlap(turquoise.go.obj.drought5h)
print(turquoise.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=4752, e.g. Brasy1G000500.v1.1 Brasy1G001900.v1.1 Brasy1G002400.v1.1
#Union size=12035, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 24892 5370
#inB   1913 4752
#Overlapping p-value=0e+00
#Odds ratio=11.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.4
###drought10h
turquoise.drought10h <- intersect(rownames(drought10h),turquoise$nodeName)
turquoise.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                              turquoise$nodeName,
                                              genome.size=36927)
turquoise.go.obj.drought10h <- testGeneOverlap(turquoise.go.obj.drought10h)
print(turquoise.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
# # listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=4876, e.g. Brasy1G000500.v1.1 Brasy1G002300.v1.1 Brasy1G002500.v1.1
#Union size=12172, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 24755 5507
#inB   1789 4876
#Overlapping p-value=0e+00
#Odds ratio=12.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.4
###drought24h
turquoise.drought24h <- intersect(rownames(drought24h),turquoise$nodeName)
turquoise.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                              turquoise$nodeName,
                                              genome.size=36927)
turquoise.go.obj.drought24h <- testGeneOverlap(turquoise.go.obj.drought24h)
print(turquoise.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=5054, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G002300.v1.1
#Union size=13165, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 23762 6500
#inB   1611 5054
#Overlapping p-value=0e+00
#Odds ratio=11.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.4

###heat1h
turquoise.heat1h <- intersect(rownames(heat1h),turquoise$nodeName)
turquoise.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                          turquoise$nodeName,
                                          genome.size=36927)
turquoise.go.obj.heat1h <- testGeneOverlap(turquoise.go.obj.heat1h)
print(turquoise.go.obj.heat1h)
####
# Detailed information about this GeneOverlap object:
#  listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=1551, e.g. Brasy1G002500.v1.1 Brasy1G003400.v1.1 Brasy1G004400.v1.1
#Union size=11513, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 25414 4848
#inB   5114 1551
#Overlapping p-value=3.9e-43
#Odds ratio=1.6
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###heat2h
turquoise.heat2h <- intersect(rownames(heat2h),turquoise$nodeName)
turquoise.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                          turquoise$nodeName,
                                          genome.size=36927)
turquoise.go.obj.heat2h <- testGeneOverlap(turquoise.go.obj.heat2h)
print(turquoise.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=1157, e.g. Brasy1G004400.v1.1 Brasy1G008700.v1.1 Brasy1G009600.v1.1
#Union size=11188, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 25739 4523
#inB   5508 1157
#Overlapping p-value=5.7e-07
#Odds ratio=1.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat5h
turquoise.heat5h <- intersect(rownames(heat5h),turquoise$nodeName)
turquoise.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                          turquoise$nodeName,
                                          genome.size=36927)
turquoise.go.obj.heat5h <- testGeneOverlap(turquoise.go.obj.heat5h)
print(turquoise.go.obj.heat5h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=731, e.g. Brasy1G004400.v1.1 Brasy1G008000.v1.1 Brasy1G009600.v1.1
#Union size=9583, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 27344 2918
#inB   5934  731
#Overlapping p-value=6.3e-04
#Odds ratio=1.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat10h
turquoise.heat10h <- intersect(rownames(heat10h),turquoise$nodeName)
turquoise.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                           turquoise$nodeName,
                                           genome.size=36927)
turquoise.go.obj.heat10h <- testGeneOverlap(turquoise.go.obj.heat10h)
print(turquoise.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=477, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G011100.v1.1
#Union size=4324, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32603 3172
#inB    675  477
#Overlapping p-value=3.8e-182
#Odds ratio=7.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###heat24h
turquoise.heat24h <- intersect(rownames(heat24h),turquoise$nodeName)
turquoise.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                           turquoise$nodeName,
                                           genome.size=36927)
turquoise.go.obj.heat24h <- testGeneOverlap(turquoise.go.obj.heat24h)
print(turquoise.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=659, e.g. Brasy1G003600.v1.1 Brasy1G010400.v1.1 Brasy1G011100.v1.1
#Union size=12047, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 24880 10895
#inB    493   659
#Overlapping p-value=5.2e-76
#Odds ratio=3.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###salt1h
turquoise.salt1h <- intersect(rownames(salt1h),turquoise$nodeName)
turquoise.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                          turquoise$nodeName,
                                          genome.size=36927)
turquoise.go.obj.salt1h <- testGeneOverlap(turquoise.go.obj.salt1h)
print(turquoise.go.obj.salt1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=1648, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=1102, e.g. Brasy1G005400.v1.1 Brasy1G010500.v1.1 Brasy1G020800.v1.1
#Union size=7211, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 29716  546
#inB   5563 1102
#Overlapping p-value=0e+00
#Odds ratio=10.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2
########
###salt2h
turquoise.salt2h <- intersect(rownames(salt2h),turquoise$nodeName)
turquoise.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                          turquoise$nodeName,
                                          genome.size=36927)
turquoise.go.obj.salt2h <- testGeneOverlap(turquoise.go.obj.salt2h)
print(turquoise.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=436, e.g. Brasy1G016400.v1.1 Brasy1G020800.v1.1 Brasy1G022500.v1.1
#Union size=6959, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 29968 294
#inB   6229 436
#Overlapping p-value=7.4e-142
#Odds ratio=7.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
turquoise.salt5h <- intersect(rownames(salt5h),turquoise$nodeName)
turquoise.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                          turquoise$nodeName,
                                          genome.size=36927)
turquoise.go.obj.salt5h <- testGeneOverlap(turquoise.go.obj.salt5h)
print(turquoise.go.obj.salt5h)
#Detailed information about this GeneOverlap object:
#  listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=213, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Union size=6824, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 30103 159
#inB   6452 213
#Overlapping p-value=1e-64
#Odds ratio=6.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
turquoise.salt10h <- intersect(rownames(salt10h),turquoise$nodeName)
turquoise.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                           turquoise$nodeName,
                                           genome.size=36927)
turquoise.go.obj.salt10h <- testGeneOverlap(turquoise.go.obj.salt10h)
print(turquoise.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=1152, e.g. Brasy1G003600.v1.1 Brasy1G003800.v1.1 Brasy1G010400.v1.1
#Intersection size=53, e.g. Brasy1G080300.v1.1 Brasy1G214000.v1.1 Brasy1G290300.v1.1
#Union size=1882, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 35045 730
#inB   1099  53
#Overlapping p-value=1.5e-07
#Odds ratio=2.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
turquoise.salt24h <- intersect(rownames(salt24h),turquoise$nodeName)
turquoise.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                           turquoise$nodeName,
                                           genome.size=36927)
turquoise.go.obj.salt24h <- testGeneOverlap(turquoise.go.obj.salt24h)
print(turquoise.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=6665, e.g. Brasy1G000500.v1.1 Brasy1G000800.v1.1 Brasy1G001900.v1.1
#Intersection size=633, e.g. Brasy1G002500.v1.1 Brasy1G008000.v1.1 Brasy1G011500.v1.1
#Union size=7458, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 29469 793
#inB   6032 633
#Overlapping p-value=5e-122
#Odds ratio=3.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

#pink module vs all of the DGE from each treatment
###drought1h
pink.drought1h <- intersect(rownames(drought1h),pink$nodeName)
pink.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                        pink$nodeName,
                                        genome.size=36927)
pink.go.obj.drought1h <- testGeneOverlap(pink.go.obj.drought1h)
print(pink.go.obj.drought1h)
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=15, e.g. Brasy1G105900.v1.1 Brasy1G122900.v1.1 Brasy1G141700.v1.1
#Union size=3781, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33146 3741
#inB     25   15
#Overlapping p-value=4.2e-06
#Odds ratio=5.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###drought2h
pink.drought2h <- intersect(rownames(drought2h),pink$nodeName)
pink.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                        pink$nodeName,
                                        genome.size=36927)
pink.go.obj.drought2h <- testGeneOverlap(pink.go.obj.drought2h)
print(pink.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=21, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Union size=6016, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30911 5976
#inB     19   21
#Overlapping p-value=1.4e-07
#Odds ratio=5.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought5h
pink.drought5h <- intersect(rownames(drought5h),pink$nodeName)
pink.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                        pink$nodeName,
                                        genome.size=36927)
pink.go.obj.drought5h <- testGeneOverlap(pink.go.obj.drought5h)
print(pink.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=20, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Union size=10142, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 26785 10102
#inB     20    20
#Overlapping p-value=2e-03
#Odds ratio=2.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###drought10h
pink.drought10h <- intersect(rownames(drought10h),pink$nodeName)
pink.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                         pink$nodeName,
                                         genome.size=36927)
pink.go.obj.drought10h <- testGeneOverlap(pink.go.obj.drought10h)
print(pink.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
#  listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=6, e.g. Brasy1G141700.v1.1 Brasy1G485600.v1.1 Brasy2G365200.v1.1
#Union size=10417, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 26510 10377
#inB     34     6
#Overlapping p-value=0.98
#Odds ratio=0.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought24h
pink.drought24h <- intersect(rownames(drought24h),pink$nodeName)
pink.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                         pink$nodeName,
                                         genome.size=36927)
pink.go.obj.drought24h <- testGeneOverlap(pink.go.obj.drought24h)
print(pink.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=16, e.g. Brasy1G105900.v1.1 Brasy1G122900.v1.1 Brasy1G485600.v1.1
#Union size=11578, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 25349 11538
#inB     24    16
#Overlapping p-value=0.15
#Odds ratio=1.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat1h
pink.heat1h <- intersect(rownames(heat1h),pink$nodeName)
pink.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                     pink$nodeName,
                                     genome.size=36927)
pink.go.obj.heat1h <- testGeneOverlap(pink.go.obj.heat1h)
print(pink.go.obj.heat1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=11, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy2G365200.v1.1
#Union size=6428, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30499 6388
#inB     29   11
#Overlapping p-value=0.074
#Odds ratio=1.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###heat2h
pink.heat2h <- intersect(rownames(heat2h),pink$nodeName)
pink.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                     pink$nodeName,
                                     genome.size=36927)
pink.go.obj.heat2h <- testGeneOverlap(pink.go.obj.heat2h)
print(pink.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=7, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Union size=5713, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31214 5673
#inB     33    7
#Overlapping p-value=0.42
#Odds ratio=1.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat5h
pink.heat5h <- intersect(rownames(heat5h),pink$nodeName)
pink.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                     pink$nodeName,
                                     genome.size=36927)
pink.go.obj.heat5h <- testGeneOverlap(pink.go.obj.heat5h)
print(pink.go.obj.heat5h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=4, e.g. Brasy1G122900.v1.1 Brasy3G162900.v1.1 Brasy3G220200.v1.1
#Union size=3685, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33242 3645
#inB     36    4
#Overlapping p-value=0.57
#Odds ratio=1.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat10h
pink.heat10h <- intersect(rownames(heat10h),pink$nodeName)
pink.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                      pink$nodeName,
                                      genome.size=36927)
pink.go.obj.heat10h <- testGeneOverlap(pink.go.obj.heat10h)
print(pink.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=4225, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=6, e.g. Brasy1G485600.v1.1 Brasy3G140000.v1.1 Brasy4G127400.v1.1
#Union size=4259, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32668 4219
#inB     34    6
#Overlapping p-value=0.31
#Odds ratio=1.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###heat24h
pink.heat24h <- intersect(rownames(heat24h),pink$nodeName)
pink.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                      pink$nodeName,
                                      genome.size=36927)
pink.go.obj.heat24h <- testGeneOverlap(pink.go.obj.heat24h)
print(pink.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=5793, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=13, e.g. Brasy1G105900.v1.1 Brasy1G122900.v1.1 Brasy1G485600.v1.1
#Union size=5820, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31107 5780
#inB     27   13
#Overlapping p-value=6.3e-03
#Odds ratio=2.6
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###salt1h
pink.salt1h <- intersect(rownames(salt1h),pink$nodeName)
pink.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                     pink$nodeName,
                                     genome.size=36927)
pink.go.obj.salt1h <- testGeneOverlap(pink.go.obj.salt1h)
print(pink.go.obj.salt1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=1648, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=1, e.g. Brasy8G117400.v1.1
#Union size=1687, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 35240 1647
#inB     39    1
#Overlapping p-value=0.84
#Odds ratio=0.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt2h
pink.salt2h <- intersect(rownames(salt2h),pink$nodeName)
pink.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                     pink$nodeName,
                                     genome.size=36927)
pink.go.obj.salt2h <- testGeneOverlap(pink.go.obj.salt2h)
print(pink.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=1, e.g. Brasy1G122900.v1.1
#Union size=769, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36158 729
#inB     39   1
#Overlapping p-value=0.55
#Odds ratio=1.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
pink.salt5h <- intersect(rownames(salt5h),pink$nodeName)
pink.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                     pink$nodeName,
                                     genome.size=36927)
pink.go.obj.salt5h <- testGeneOverlap(pink.go.obj.salt5h)
print(pink.go.obj.salt5h)
#Detailed information about this GeneOverlap object:
#  listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=0, e.g. 
#Union size=412, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36515 372
#inB     40   0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
pink.salt10h <- intersect(rownames(salt10h),pink$nodeName)
pink.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                      pink$nodeName,
                                      genome.size=36927)
pink.go.obj.salt10h <- testGeneOverlap(pink.go.obj.salt10h)
print(pink.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=0, e.g. 
#Union size=823, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36104 783
#inB     40   0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
pink.salt24h <- intersect(rownames(salt24h),pink$nodeName)
pink.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                      pink$nodeName,
                                      genome.size=36927)
pink.go.obj.salt24h <- testGeneOverlap(pink.go.obj.salt24h)
print(pink.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=40, e.g. Brasy1G039900.v1.1 Brasy1G105900.v1.1 Brasy1G122900.v1.1
#Intersection size=4, e.g. Brasy1G485600.v1.1 Brasy2G246100.v1.1 Brasy6G065500.v1.1
#Union size=1462, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 35465 1422
#inB     36    4
#Overlapping p-value=0.067
#Odds ratio=2.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
#blue module vs all of the DGE from each treatment
###drought1h
blue.drought1h <- intersect(rownames(drought1h),blue$nodeName)
blue.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                        blue$nodeName,
                                        genome.size=36927)
blue.go.obj.drought1h <- testGeneOverlap(blue.go.obj.drought1h)
print(blue.go.obj.drought1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=945, e.g. Brasy1G002800.v1.1 Brasy1G004000.v1.1 Brasy1G005900.v1.1
#Union size=6442, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30485 2811
#inB   2686  945
#Overlapping p-value=2.2e-186
#Odds ratio=3.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
########
###drought2h
blue.drought2h <- intersect(rownames(drought2h),blue$nodeName)
blue.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                        blue$nodeName,
                                        genome.size=36927)
blue.go.obj.drought2h <- testGeneOverlap(blue.go.obj.drought2h)
print(blue.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=1450, e.g. Brasy1G001000.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Union size=8178, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 28749 4547
#inB   2181 1450
#Overlapping p-value=3.3e-291
#Odds ratio=4.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2
###drought5h
blue.drought5h <- intersect(rownames(drought5h),blue$nodeName)
blue.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                        blue$nodeName,
                                        genome.size=36927)
blue.go.obj.drought5h <- testGeneOverlap(blue.go.obj.drought5h)
print(blue.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=2483, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Union size=11270, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 25657 7639
#inB   1148 2483
#Overlapping p-value=0e+00
#Odds ratio=7.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###drought10h
blue.drought10h <- intersect(rownames(drought10h),blue$nodeName)
blue.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                         blue$nodeName,
                                         genome.size=36927)
blue.go.obj.drought10h <- testGeneOverlap(blue.go.obj.drought10h)
print(blue.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
#  listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=2549, e.g. Brasy1G000400.v1.1 Brasy1G001000.v1.1 Brasy1G001800.v1.1
#Union size=11465, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 25462 7834
#inB   1082 2549
#Overlapping p-value=0e+00
#Odds ratio=7.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2
###drought24h
blue.drought24h <- intersect(rownames(drought24h),blue$nodeName)
blue.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                         blue$nodeName,
                                         genome.size=36927)
blue.go.obj.drought24h <- testGeneOverlap(blue.go.obj.drought24h)
print(blue.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=2712, e.g. Brasy1G000400.v1.1 Brasy1G001000.v1.1 Brasy1G002800.v1.1
#Union size=12473, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 24454 8842
#inB    919 2712
#Overlapping p-value=0e+00
#Odds ratio=8.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###heat1h
blue.heat1h <- intersect(rownames(heat1h),blue$nodeName)
blue.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                     blue$nodeName,
                                     genome.size=36927)
blue.go.obj.heat1h <- testGeneOverlap(blue.go.obj.heat1h)
print(blue.go.obj.heat1h)
####
#Detailed information about this GeneOverlap object:
# listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=1772, e.g. Brasy1G000400.v1.1 Brasy1G001000.v1.1 Brasy1G004000.v1.1
#Union size=8258, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 28669 4627
#inB   1859 1772
#Overlapping p-value=0e+00
#Odds ratio=5.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2
###heat2h
blue.heat2h <- intersect(rownames(heat2h),blue$nodeName)
blue.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                     blue$nodeName,
                                     genome.size=36927)
blue.go.obj.heat2h <- testGeneOverlap(blue.go.obj.heat2h)
print(blue.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=1318, e.g. Brasy1G000400.v1.1 Brasy1G001800.v1.1 Brasy1G005000.v1.1
#Union size=7993, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 28934 4362
#inB   2313 1318
#Overlapping p-value=8.3e-239
#Odds ratio=3.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###heat5h
blue.heat5h <- intersect(rownames(heat5h),blue$nodeName)
blue.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                     blue$nodeName,
                                     genome.size=36927)
blue.go.obj.heat5h <- testGeneOverlap(blue.go.obj.heat5h)
print(blue.go.obj.heat5h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=601, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G005000.v1.1
#Union size=6679, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30248 3048
#inB   3030  601
#Overlapping p-value=5.7e-40
#Odds ratio=2.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat10h
blue.heat10h <- intersect(rownames(heat10h),blue$nodeName)
blue.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                      blue$nodeName,
                                      genome.size=36927)
blue.go.obj.heat10h <- testGeneOverlap(blue.go.obj.heat10h)
print(blue.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=4225, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=809, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G005000.v1.1
#Union size=7047, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 29880 3416
#inB   2822  809
#Overlapping p-value=1.1e-86
#Odds ratio=2.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###heat24h
blue.heat24h <- intersect(rownames(heat24h),blue$nodeName)
blue.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                      blue$nodeName,
                                      genome.size=36927)
blue.go.obj.heat24h <- testGeneOverlap(blue.go.obj.heat24h)
print(blue.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=5793, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=1007, e.g. Brasy1G001000.v1.1 Brasy1G002800.v1.1 Brasy1G005000.v1.1
#Union size=8417, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 28510 4786
#inB   2624 1007
#Overlapping p-value=2e-85
#Odds ratio=2.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###salt1h
blue.salt1h <- intersect(rownames(salt1h),blue$nodeName)
blue.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                     blue$nodeName,
                                     genome.size=36927)
blue.go.obj.salt1h <- testGeneOverlap(blue.go.obj.salt1h)
print(blue.go.obj.salt1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=1648, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=354, e.g. Brasy1G005900.v1.1 Brasy1G026700.v1.1 Brasy1G032300.v1.1
#Union size=4925, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32002 1294
#inB   3277  354
#Overlapping p-value=4e-47
#Odds ratio=2.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
########
###salt2h
blue.salt2h <- intersect(rownames(salt2h),blue$nodeName)
blue.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                     blue$nodeName,
                                     genome.size=36927)
blue.go.obj.salt2h <- testGeneOverlap(blue.go.obj.salt2h)
print(blue.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=121, e.g. Brasy1G156700.v1.1 Brasy1G188300.v1.1 Brasy1G211700.v1.1
#Union size=4240, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 32687 609
#inB   3510 121
#Overlapping p-value=7.7e-09
#Odds ratio=1.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

blue.salt5h <- intersect(rownames(salt5h),blue$nodeName)
blue.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                     blue$nodeName,
                                     genome.size=36927)
blue.go.obj.salt5h <- testGeneOverlap(blue.go.obj.salt5h)
print(blue.go.obj.salt5h)
#Detailed information about this GeneOverlap object:
#  listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=50, e.g. Brasy1G124300.v1.1 Brasy1G214600.v1.1 Brasy1G235500.v1.1
#Union size=3953, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 32974 322
#inB   3581  50
#Overlapping p-value=0.015
#Odds ratio=1.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
blue.salt10h <- intersect(rownames(salt10h),blue$nodeName)
blue.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                      blue$nodeName,
                                      genome.size=36927)
blue.go.obj.salt10h <- testGeneOverlap(blue.go.obj.salt10h)
print(blue.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=155, e.g. Brasy1G004000.v1.1 Brasy1G033400.v1.1 Brasy1G091300.v1.1
#Union size=4259, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 32668 628
#inB   3476 155
#Overlapping p-value=1.7e-17
#Odds ratio=2.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
blue.salt24h <- intersect(rownames(salt24h),blue$nodeName)
blue.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                      blue$nodeName,
                                      genome.size=36927)
blue.go.obj.salt24h <- testGeneOverlap(blue.go.obj.salt24h)
print(blue.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=3631, e.g. Brasy1G000400.v1.1 Brasy1G000900.v1.1 Brasy1G001000.v1.1
#Intersection size=278, e.g. Brasy1G007700.v1.1 Brasy1G033400.v1.1 Brasy1G042400.v1.1
#Union size=4779, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32148 1148
#inB   3353  278
#Overlapping p-value=1.7e-29
#Odds ratio=2.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

#yellow module vs all of the DGE from each treatment
###drought1h
yellow.drought1h <- intersect(rownames(drought1h),yellow$nodeName)
yellow.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                          yellow$nodeName,
                                          genome.size=36927)
yellow.go.obj.drought1h <- testGeneOverlap(yellow.go.obj.drought1h)
print(yellow.go.obj.drought1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=100, e.g. Brasy1G028700.v1.1 Brasy1G030900.v1.1 Brasy1G042200.v1.1
#Union size=5284, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31643 3656
#inB   1528  100
#Overlapping p-value=1
#Odds ratio=0.6
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###drought2h
yellow.drought2h <- intersect(rownames(drought2h),yellow$nodeName)
yellow.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                          yellow$nodeName,
                                          genome.size=36927)
yellow.go.obj.drought2h <- testGeneOverlap(yellow.go.obj.drought2h)
print(yellow.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=214, e.g. Brasy1G028700.v1.1 Brasy1G030900.v1.1 Brasy1G039100.v1.1
#Union size=7411, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 29516 5783
#inB   1414  214
#Overlapping p-value=1
#Odds ratio=0.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought5h
yellow.drought5h <- intersect(rownames(drought5h),yellow$nodeName)
yellow.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                          yellow$nodeName,
                                          genome.size=36927)
yellow.go.obj.drought5h <- testGeneOverlap(yellow.go.obj.drought5h)
print(yellow.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=422, e.g. Brasy1G021200.v1.1 Brasy1G036400.v1.1 Brasy1G039600.v1.1
#Union size=11328, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 25599 9700
#inB   1206  422
#Overlapping p-value=0.92
#Odds ratio=0.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought10h
yellow.drought10h <- intersect(rownames(drought10h),yellow$nodeName)
yellow.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                           yellow$nodeName,
                                           genome.size=36927)
yellow.go.obj.drought10h <- testGeneOverlap(yellow.go.obj.drought10h)
print(yellow.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
#  listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=354, e.g. Brasy1G028700.v1.1 Brasy1G065000.v1.1 Brasy1G065700.v1.1
#Union size=11657, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 25270 10029
#inB   1274   354
#Overlapping p-value=1
#Odds ratio=0.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought24h
yellow.drought24h <- intersect(rownames(drought24h),yellow$nodeName)
yellow.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                           yellow$nodeName,
                                           genome.size=36927)
yellow.go.obj.drought24h <- testGeneOverlap(yellow.go.obj.drought24h)
print(yellow.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=697, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G021200.v1.1
#Union size=12485, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 24442 10857
#inB    931   697
#Overlapping p-value=9.9e-24
#Odds ratio=1.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat1h
yellow.heat1h <- intersect(rownames(heat1h),yellow$nodeName)
yellow.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                       yellow$nodeName,
                                       genome.size=36927)
yellow.go.obj.heat1h <- testGeneOverlap(yellow.go.obj.heat1h)
print(yellow.go.obj.heat1h)
#Detailed information about this GeneOverlap object:
#  listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=446, e.g. Brasy1G000200.v1.1 Brasy1G044900.v1.1 Brasy1G069100.v1.1
#Union size=7581, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 29346 5953
#inB   1182  446
#Overlapping p-value=3e-25
#Odds ratio=1.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1
###heat2h
yellow.heat2h <- intersect(rownames(heat2h),yellow$nodeName)
yellow.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                       yellow$nodeName,
                                       genome.size=36927)
yellow.go.obj.heat2h <- testGeneOverlap(yellow.go.obj.heat2h)
print(yellow.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=742, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G028700.v1.1
#Union size=6566, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30361 4938
##inB    886  742
#Overlapping p-value=3.2e-194
#Odds ratio=5.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

###heat5h
yellow.heat5h <- intersect(rownames(heat5h),yellow$nodeName)
yellow.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                       yellow$nodeName,
                                       genome.size=36927)
yellow.go.obj.heat5h <- testGeneOverlap(yellow.go.obj.heat5h)
print(yellow.go.obj.heat5h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=727, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G030900.v1.1
#Union size=4550, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32377 2922
#inB    901  727
#Overlapping p-value=4.5e-311
#Odds ratio=8.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2

###heat10h
yellow.heat10h <- intersect(rownames(heat10h),yellow$nodeName)
yellow.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                        yellow$nodeName,
                                        genome.size=36927)
yellow.go.obj.heat10h <- testGeneOverlap(yellow.go.obj.heat10h)
print(yellow.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=4225, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=890, e.g. Brasy1G000100.v1.1 Brasy1G028700.v1.1 Brasy1G030900.v1.1
#Union size=4963, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31964 3335
#inB    738  890
#Overlapping p-value=0e+00
#Odds ratio=11.6
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2
###heat24h
yellow.heat24h <- intersect(rownames(heat24h),yellow$nodeName)
yellow.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                        yellow$nodeName,
                                        genome.size=36927)
yellow.go.obj.heat24h <- testGeneOverlap(yellow.go.obj.heat24h)
print(yellow.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=5793, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=1204, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Union size=6217, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30710 4589
#inB    424 1204
#Overlapping p-value=0e+00
#Odds ratio=19.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.2
###salt1h
yellow.salt1h <- intersect(rownames(salt1h),yellow$nodeName)
yellow.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                       yellow$nodeName,
                                       genome.size=36927)
yellow.go.obj.salt1h <- testGeneOverlap(yellow.go.obj.salt1h)
print(yellow.go.obj.salt1h)
#Detailed information about this GeneOverlap object:
#  listA size=1648, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
##Intersection size=35, e.g. Brasy1G039100.v1.1 Brasy1G049300.v1.1 Brasy1G088300.v1.1
#Union size=3241, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33686 1613
#inB   1593   35
#Overlapping p-value=1
#Odds ratio=0.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###salt2h
yellow.salt2h <- intersect(rownames(salt2h),yellow$nodeName)
yellow.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                       yellow$nodeName,
                                       genome.size=36927)
yellow.go.obj.salt2h <- testGeneOverlap(yellow.go.obj.salt2h)
print(yellow.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=48, e.g. Brasy1G030900.v1.1 Brasy1G039100.v1.1 Brasy1G086000.v1.1
#Union size=2310, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 34617 682
#inB   1580  48
#Overlapping p-value=4.1e-03
#Odds ratio=1.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

yellow.salt5h <- intersect(rownames(salt5h),yellow$nodeName)
yellow.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                       yellow$nodeName,
                                       genome.size=36927)
yellow.go.obj.salt5h <- testGeneOverlap(yellow.go.obj.salt5h)
print(yellow.go.obj.salt5h)
#Detailed information about this GeneOverlap object:
#  listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=42, e.g. Brasy1G030900.v1.1 Brasy1G042200.v1.1 Brasy1G104600.v1.1
#Union size=1958, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 34969 330
#inB   1586  42
#Overlapping p-value=2.9e-08
#Odds ratio=2.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
yellow.salt10h <- intersect(rownames(salt10h),yellow$nodeName)
yellow.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                        yellow$nodeName,
                                        genome.size=36927)
yellow.go.obj.salt10h <- testGeneOverlap(yellow.go.obj.salt10h)
print(yellow.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=92, e.g. Brasy1G030900.v1.1 Brasy1G042200.v1.1 Brasy1G078300.v1.1
#Union size=2319, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 34608 691
#inB   1536  92
#Overlapping p-value=1.5e-17
#Odds ratio=3.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
yellow.salt24h <- intersect(rownames(salt24h),yellow$nodeName)
yellow.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                        yellow$nodeName,
                                        genome.size=36927)
yellow.go.obj.salt24h <- testGeneOverlap(yellow.go.obj.salt24h)
print(yellow.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=1628, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Intersection size=223, e.g. Brasy1G039100.v1.1 Brasy1G042200.v1.1 Brasy1G063600.v1.1
#Union size=2831, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 34096 1203
#inB   1405  223
#Overlapping p-value=2.6e-63
#Odds ratio=4.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.1

#black module vs all of the DGE from each treatment
###drought1h
black.drought1h <- intersect(rownames(drought1h),black$nodeName)
black.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                         black$nodeName,
                                         genome.size=36927)
black.go.obj.drought1h <- testGeneOverlap(black.go.obj.drought1h)
print(black.go.obj.drought1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=1, e.g. Brasy8G291900.v1.1
#Union size=3829, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33098 3755
#inB     73    1
#Overlapping p-value=1
#Odds ratio=0.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###drought2h
black.drought2h <- intersect(rownames(drought2h),black$nodeName)
black.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                         black$nodeName,
                                         genome.size=36927)
black.go.obj.drought2h <- testGeneOverlap(black.go.obj.drought2h)
print(black.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=8, e.g. Brasy1G321700.v1.1 Brasy1G321800.v1.1 Brasy1G322100.v1.1
#Union size=6063, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30864 5989
#inB     66    8
#Overlapping p-value=0.93
#Odds ratio=0.6
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought5h
black.drought5h <- intersect(rownames(drought5h),black$nodeName)
black.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                         black$nodeName,
                                         genome.size=36927)
black.go.obj.drought5h <- testGeneOverlap(black.go.obj.drought5h)
print(black.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=21, e.g. Brasy1G320900.v1.1 Brasy1G321500.v1.1 Brasy1G321700.v1.1
#Union size=10175, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 26752 10101
#inB     53    21
#Overlapping p-value=0.47
#Odds ratio=1.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###drought10h
black.drought10h <- intersect(rownames(drought10h),black$nodeName)
black.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                          black$nodeName,
                                          genome.size=36927)
black.go.obj.drought10h <- testGeneOverlap(black.go.obj.drought10h)
print(black.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
#  listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=4, e.g. Brasy1G321800.v1.1 Brasy8G297300.v1.1 BrasyJ073900.v1.1
#Union size=10453, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 26474 10379
#inB     70     4
#Overlapping p-value=1
#Odds ratio=0.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought24h
black.drought24h <- intersect(rownames(drought24h),black$nodeName)
black.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                          black$nodeName,
                                          genome.size=36927)
black.go.obj.drought24h <- testGeneOverlap(black.go.obj.drought24h)
print(black.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=28, e.g. Brasy1G320900.v1.1 Brasy1G321500.v1.1 Brasy1G321700.v1.1
#Union size=11600, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 25327 11526
#inB     46    28
#Overlapping p-value=0.14
#Odds ratio=1.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat1h
black.heat1h <- intersect(rownames(heat1h),black$nodeName)
black.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                      black$nodeName,
                                      genome.size=36927)
black.go.obj.heat1h <- testGeneOverlap(black.go.obj.heat1h)
print(black.go.obj.heat1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=1, e.g. Brasy1G321800.v1.1
#Union size=6472, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30455 6398
#inB     73    1
#Overlapping p-value=1
#Odds ratio=0.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat2h
black.heat2h <- intersect(rownames(heat2h),black$nodeName)
black.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                      black$nodeName,
                                      genome.size=36927)
black.go.obj.heat2h <- testGeneOverlap(black.go.obj.heat2h)
print(black.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=1, e.g. Brasy1G253000.v1.1
#Union size=5753, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31174 5679
#inB     73    1
#Overlapping p-value=1
#Odds ratio=0.1
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat5h
black.heat5h <- intersect(rownames(heat5h),black$nodeName)
black.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                      black$nodeName,
                                      genome.size=36927)
black.go.obj.heat5h <- testGeneOverlap(black.go.obj.heat5h)
print(black.go.obj.heat5h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=3723, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33204 3649
#inB     74    0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat10h
black.heat10h <- intersect(rownames(heat10h),black$nodeName)
black.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                       black$nodeName,
                                       genome.size=36927)
black.go.obj.heat10h <- testGeneOverlap(black.go.obj.heat10h)
print(black.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=4225, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=4299, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32628 4225
#inB     74    0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###heat24h
black.heat24h <- intersect(rownames(heat24h),black$nodeName)
black.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                       black$nodeName,
                                       genome.size=36927)
black.go.obj.heat24h <- testGeneOverlap(black.go.obj.heat24h)
print(black.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=5793, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=5867, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31060 5793
#inB     74    0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###salt1h
black.salt1h <- intersect(rownames(salt1h),black$nodeName)
black.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                      black$nodeName,
                                      genome.size=36927)
black.go.obj.salt1h <- testGeneOverlap(black.go.obj.salt1h)
print(black.go.obj.salt1h)
#Detailed information about this GeneOverlap object:
#  listA size=1648, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=1722, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 35205 1648
#inB     74    0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###salt2h
black.salt2h <- intersect(rownames(salt2h),black$nodeName)
black.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                      black$nodeName,
                                      genome.size=36927)
black.go.obj.salt2h <- testGeneOverlap(black.go.obj.salt2h)
print(black.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=804, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36123 730
#inB     74   0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
black.salt5h <- intersect(rownames(salt5h),black$nodeName)
black.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                      black$nodeName,
                                      genome.size=36927)
black.go.obj.salt5h <- testGeneOverlap(black.go.obj.salt5h)
print(black.go.obj.salt5h)
#Detailed information about this GeneOverlap object:
#  listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=446, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36481 372
#inB     74   0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
black.salt10h <- intersect(rownames(salt10h),black$nodeName)
black.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                       black$nodeName,
                                       genome.size=36927)
black.go.obj.salt10h <- testGeneOverlap(black.go.obj.salt10h)
print(black.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=857, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36070 783
#inB     74   0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
black.salt24h <- intersect(rownames(salt24h),black$nodeName)
black.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                       black$nodeName,
                                       genome.size=36927)
black.go.obj.salt24h <- testGeneOverlap(black.go.obj.salt24h)
print(black.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=74, e.g. Brasy1G247700.v1.1 Brasy1G247900.v1.1 Brasy1G248200.v1.1
#Intersection size=0, e.g. 
#Union size=1500, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 35427 1426
#inB     74    0
#Overlapping p-value=1
#Odds ratio=0.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

#grey module vs all of the DGE from each treatment
###drought1h
grey.drought1h <- intersect(rownames(drought1h),grey$nodeName)
grey.go.obj.drought1h <- newGeneOverlap(rownames(drought1h),
                                        grey$nodeName,
                                        genome.size=36927)
grey.go.obj.drought1h <- testGeneOverlap(grey.go.obj.drought1h)
print(grey.go.obj.drought1h)
####
#Detailed information about this GeneOverlap object:
#  listA size=3756, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=4, e.g. Brasy2G073300.v1.1 Brasy2G437100.v1.1 BrasyJ017700.v1.1
#Union size=3805, e.g. Brasy1G002500.v1.1 Brasy1G002800.v1.1 Brasy1G004000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33122 3752
#inB     49    4
#Overlapping p-value=0.8
#Odds ratio=0.7
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###drought2h
grey.drought2h <- intersect(rownames(drought2h),grey$nodeName)
grey.go.obj.drought2h <- newGeneOverlap(rownames(drought2h),
                                        grey$nodeName,
                                        genome.size=36927)
grey.go.obj.drought2h <- testGeneOverlap(grey.go.obj.drought2h)
print(grey.go.obj.drought2h)
#Detailed information about this GeneOverlap object:
#  listA size=5997, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=3, e.g. Brasy2G073300.v1.1 Brasy2G258200.v1.1 Brasy2G437100.v1.1
#Union size=6047, e.g. Brasy1G000500.v1.1 Brasy1G001000.v1.1 Brasy1G002500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30880 5994
#inB     50    3
#Overlapping p-value=0.99
#Odds ratio=0.3
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought5h
grey.drought5h <- intersect(rownames(drought5h),grey$nodeName)
grey.go.obj.drought5h <- newGeneOverlap(rownames(drought5h),
                                        grey$nodeName,
                                        genome.size=36927)
grey.go.obj.drought5h <- testGeneOverlap(grey.go.obj.drought5h)
print(grey.go.obj.drought5h)
#Detailed information about this GeneOverlap object:
#  listA size=10122, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=4, e.g. Brasy2G050400.v1.1 Brasy2G146300.v1.1 Brasy2G258200.v1.1
#Union size=10171, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 26756 10118
#inB     49     4
#Overlapping p-value=1
#Odds ratio=0.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought10h
grey.drought10h <- intersect(rownames(drought10h),grey$nodeName)
grey.go.obj.drought10h <- newGeneOverlap(rownames(drought10h),
                                         grey$nodeName,
                                         genome.size=36927)
grey.go.obj.drought10h <- testGeneOverlap(grey.go.obj.drought10h)
print(grey.go.obj.drought10h)
#Detailed information about this GeneOverlap object:
#  listA size=10383, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=8, e.g. Brasy1G197200.v1.1 Brasy2G258200.v1.1 Brasy9G136400.v1.1
#Union size=10428, e.g. Brasy1G000400.v1.1 Brasy1G000500.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 26499 10375
#inB     45     8
#Overlapping p-value=0.99
#Odds ratio=0.5
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###drought24h
grey.drought24h <- intersect(rownames(drought24h),grey$nodeName)
grey.go.obj.drought24h <- newGeneOverlap(rownames(drought24h),
                                         grey$nodeName,
                                         genome.size=36927)
grey.go.obj.drought24h <- testGeneOverlap(grey.go.obj.drought24h)
print(grey.go.obj.drought24h)
#Detailed information about this GeneOverlap object:
#  listA size=11554, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=5, e.g. Brasy2G049600.v1.1 Brasy5G446800.v1.1 BrasyJ017700.v1.1
#Union size=11602, e.g. Brasy1G000100.v1.1 Brasy1G000300.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA   inA
#notB 25325 11549
#inB     48     5
#Overlapping p-value=1
#Odds ratio=0.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat1h
grey.heat1h <- intersect(rownames(heat1h),grey$nodeName)
grey.go.obj.heat1h <- newGeneOverlap(rownames(heat1h),
                                     grey$nodeName,
                                     genome.size=36927)
grey.go.obj.heat1h <- testGeneOverlap(grey.go.obj.heat1h)
print(grey.go.obj.heat1h)
#Detailed information about this GeneOverlap object:
#  listA size=6399, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=2, e.g. Brasy5G446800.v1.1 BrasyJ017700.v1.1
#Union size=6450, e.g. Brasy1G000200.v1.1 Brasy1G000400.v1.1 Brasy1G001000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 30477 6397
#inB     51    2
#Overlapping p-value=1
#Odds ratio=0.2
###heat2h
grey.heat2h <- intersect(rownames(heat2h),grey$nodeName)
grey.go.obj.heat2h <- newGeneOverlap(rownames(heat2h),
                                     grey$nodeName,
                                     genome.size=36927)
grey.go.obj.heat2h <- testGeneOverlap(grey.go.obj.heat2h)
print(grey.go.obj.heat2h)
#Detailed information about this GeneOverlap object:
#  listA size=5680, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=8, e.g. Brasy2G023200.v1.1 Brasy2G050900.v1.1 Brasy2G258200.v1.1
#Union size=5725, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31202 5672
#inB     45    8
#Overlapping p-value=0.58
#Odds ratio=1.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat5h
grey.heat5h <- intersect(rownames(heat5h),grey$nodeName)
grey.go.obj.heat5h <- newGeneOverlap(rownames(heat5h),
                                     grey$nodeName,
                                     genome.size=36927)
grey.go.obj.heat5h <- testGeneOverlap(grey.go.obj.heat5h)
print(grey.go.obj.heat5h)
#Detailed information about this GeneOverlap object:
#  listA size=3649, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=1, e.g. Brasy2G258200.v1.1
#Union size=3701, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000400.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 33226 3648
#inB     52    1
#Overlapping p-value=1
#Odds ratio=0.2
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###heat10h
grey.heat10h <- intersect(rownames(heat10h),grey$nodeName)
grey.go.obj.heat10h <- newGeneOverlap(rownames(heat10h),
                                      grey$nodeName,
                                      genome.size=36927)
grey.go.obj.heat10h <- testGeneOverlap(grey.go.obj.heat10h)
print(grey.go.obj.heat10h)
#Detailed information about this GeneOverlap object:
#  listA size=4225, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=8, e.g. Brasy1G197200.v1.1 Brasy2G073300.v1.1 Brasy2G258200.v1.1
#Union size=4270, e.g. Brasy1G000100.v1.1 Brasy1G000400.v1.1 Brasy1G000900.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 32657 4217
#inB     45    8
#Overlapping p-value=0.26
#Odds ratio=1.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###heat24h
grey.heat24h <- intersect(rownames(heat24h),grey$nodeName)
grey.go.obj.heat24h <- newGeneOverlap(rownames(heat24h),
                                      grey$nodeName,
                                      genome.size=36927)
grey.go.obj.heat24h <- testGeneOverlap(grey.go.obj.heat24h)
print(grey.go.obj.heat24h)
#Detailed information about this GeneOverlap object:
#  listA size=5793, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=8, e.g. Brasy2G073300.v1.1 Brasy5G446800.v1.1 Brasy9G133400.v1.1
#Union size=5838, e.g. Brasy1G000100.v1.1 Brasy1G000200.v1.1 Brasy1G000300.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 31089 5785
#inB     45    8
#Overlapping p-value=0.61
#Odds ratio=1.0
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0

###salt1h
grey.salt1h <- intersect(rownames(salt1h),grey$nodeName)
grey.go.obj.salt1h <- newGeneOverlap(rownames(salt1h),
                                     grey$nodeName,
                                     genome.size=36927)
grey.go.obj.salt1h <- testGeneOverlap(grey.go.obj.salt1h)
print(grey.go.obj.salt1h)
#Detailed information about this GeneOverlap object:
#  listA size=1648, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=2, e.g. Brasy2G073300.v1.1 BrasyJ017700.v1.1
#Union size=1699, e.g. Brasy1G005400.v1.1 Brasy1G005900.v1.1 Brasy1G010500.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 35228 1646
#inB     51    2
#Overlapping p-value=0.69
#Odds ratio=0.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
########
###salt2h
grey.salt2h <- intersect(rownames(salt2h),grey$nodeName)
grey.go.obj.salt2h <- newGeneOverlap(rownames(salt2h),
                                     grey$nodeName,
                                     genome.size=36927)
grey.go.obj.salt2h <- testGeneOverlap(grey.go.obj.salt2h)
print(grey.go.obj.salt2h)
#Detailed information about this GeneOverlap object:
#  listA size=730, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=6, e.g. Brasy4G147000.v1.1 Brasy9G136400.v1.1 BrasyJ013700.v1.1
#Union size=777, e.g. Brasy1G003600.v1.1 Brasy1G011100.v1.1 Brasy1G016400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36150 724
#inB     47   6
#Overlapping p-value=6.1e-04
#Odds ratio=6.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
grey.salt5h <- intersect(rownames(salt5h),grey$nodeName)
grey.go.obj.salt5h <- newGeneOverlap(rownames(salt5h),
                                     grey$nodeName,
                                     genome.size=36927)
grey.go.obj.salt5h <- testGeneOverlap(grey.go.obj.salt5h)
print(grey.go.obj.salt5h)

#Detailed information about this GeneOverlap object:
#listA size=372, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=1, e.g. Brasy2G258200.v1.1
#Union size=424, e.g. Brasy1G011600.v1.1 Brasy1G019200.v1.1 Brasy1G026400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36503 371
#inB     52   1
#Overlapping p-value=0.42
#Odds ratio=1.9
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt10h
grey.salt10h <- intersect(rownames(salt10h),grey$nodeName)
grey.go.obj.salt10h <- newGeneOverlap(rownames(salt10h),
                                      grey$nodeName,
                                      genome.size=36927)
grey.go.obj.salt10h <- testGeneOverlap(grey.go.obj.salt10h)
print(grey.go.obj.salt10h)
#Detailed information about this GeneOverlap object:
#  listA size=783, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=4, e.g. Brasy1G197200.v1.1 Brasy1G249700.v1.1 Brasy2G073300.v1.1
#Union size=832, e.g. Brasy1G002500.v1.1 Brasy1G004000.v1.1 Brasy1G004400.v1.1
#Genome size=36927
# Contingency Table:
#notA inA
#notB 36095 779
#inB     49   4
#Overlapping p-value=0.026
#Odds ratio=3.8
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
###salt24h
grey.salt24h <- intersect(rownames(salt24h),grey$nodeName)
grey.go.obj.salt24h <- newGeneOverlap(rownames(salt24h),
                                      grey$nodeName,
                                      genome.size=36927)
grey.go.obj.salt24h <- testGeneOverlap(grey.go.obj.salt24h)
print(grey.go.obj.salt24h)
#Detailed information about this GeneOverlap object:
#  listA size=1426, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#listB size=53, e.g. Brasy1G025300.v1.1 Brasy1G045300.v1.1 Brasy1G197200.v1.1
#Intersection size=8, e.g. Brasy2G073300.v1.1 Brasy2G487000.v1.1 Brasy8G032700.v1.1
#Union size=1471, e.g. Brasy1G002500.v1.1 Brasy1G007700.v1.1 Brasy1G008000.v1.1
#Genome size=36927
# Contingency Table:
#notA  inA
#notB 35456 1418
#inB     45    8
#Overlapping p-value=9.2e-04
#Odds ratio=4.4
#Overlap tested using Fisher's exact test (alternative=greater)
#Jaccard Index=0.0
