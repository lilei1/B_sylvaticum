# Bsylvaticum

steps for denova annotate the LTR:

### 1) Run LTR-finder (Download from [here](https://github.com/oushujun/LTR_FINDER_parallel))

The unmasked genome is in this directory:
`/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa`

```
perl LTR_FINDER_parallel -seq /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa -size 100000000 -threads 36 -time 30000 -v -harvest_out


perl LTR_FINDER_parallel -seq /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa -size 5000000 -threads 36 -time 1500 -v -harvest_out
```



### 2) Run LTR-harvest

Use this version `gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz`

```
./gt suffixerator   -db /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa   -indexname Bsylvaticum_490_v1.0.fa   -tis -suf -lcp -des -ssp -sds -dna

./gt ltrharvest -index Bsylvaticum_490_v1.0.fa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 >Bsylvaticum_490_v1.0.harvest.scn
```

### 3) Install the LTR_retriever

```
module load python/3.6-anaconda-4.4

conda create -n LTR_retriever

source activate LTR_retriever

conda install -c conda-forge perl perl-text-soundex

conda install -c bioconda cd-hit

conda install -c bioconda/label/cf201901 repeatmasker
```

### 4) Run LTR_retriever

```
./LTR_retriever -genome /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bsylvaticum_490_v1.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Bsylvaticum_490_v1.0.fa.finder.combine.scn -v -threads 10

./LTR_retriever -genome /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bsylvaticum_490_v1.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Bsylvaticum_490_v1.0.fa.finder.combine.scn -v -threads 10
```

Results: 

---

##########################
### LTR_retriever v2.6 ###
##########################

Contributors: Shujun Ou, Ning Jiang

For LTR_retriever, please cite:

Ou S and Jiang N (2018). LTR_retriever: A Highly Accurate and Sensitive Program for Identification of Long Terminal Repeat Retrotransposons. Plant Physiol. 176(2): 1410-1422.

For LAI, please cite:

Ou S, Chen J, Jiang N (2018). Assessing genome assembly quality using the LTR Assembly Index (LAI). Nucleic Acids Res. 2018;46(21):e126.

Parameters: -genome /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bsylvaticum_490_v1.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Bsylvaticum_490_v1.0.fa.finder.combine.scn -v -threads 10


Mon Jul  8 15:45:37 PDT 2019 Dependency checking: All passed!
Mon Jul  8 15:46:16 PDT 2019 LTR_retriever is starting from the Init step.
Mon Jul  8 15:46:17 PDT 2019 The longest sequence ID in the genome contains 16 characters, which is longer than the limit (15)
Trying to reformat seq IDs...
Attempt 1...
Mon Jul  8 15:46:22 PDT 2019 Seq ID conversion successful!

Mon Jul  8 15:46:22 PDT 2019 Start to convert inputs...
Total candidates: 3560
Total uniq candidates: 3560

Mon Jul  8 15:46:27 PDT 2019 Module 1: Start to clean up candidates...
Sequences with 10 missing bp or 0.8 missing data rate will be discarded.
Sequences containing tandem repeats will be discarded.

Mon Jul  8 15:48:10 PDT 2019 3430 clean candidates remained

Mon Jul  8 15:48:10 PDT 2019 Modules 2-5: Start to analyze the structure of candidates...
The terminal motif, TSD, boundary, orientation, age, and superfamily will be identified in this step.

Mon Jul  8 15:53:17 PDT 2019 Intact LTR-RT found: 1553

Mon Jul  8 15:55:18 PDT 2019 Module 6: Start to analyze truncated LTR-RTs...
Truncated LTR-RTs without the intact version will be retained in the LTR-RT library.
Use -notrunc if you don't want to keep them.

Mon Jul  8 15:55:18 PDT 2019 372 truncated LTR-RTs found
Mon Jul  8 15:56:32 PDT 2019 46 truncated LTR sequences have added to the library

Mon Jul  8 15:56:32 PDT 2019 Module 5: Start to remove DNA TE and LINE transposases, and remove plant protein sequences...
Total library sequences: 1463
Mon Jul  8 16:08:23 PDT 2019 Retained clean sequence: 1462

Mon Jul  8 16:08:23 PDT 2019 Sequence clustering for /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.ltrTE ...
Mon Jul  8 16:08:23 PDT 2019 Unique lib sequence: 1461

Mon Jul  8 16:09:37 PDT 2019 Module 6: Start to remove nested insertions in internal regions...
Mon Jul  8 16:11:15 PDT 2019 Raw internal region size (bit): 4786419
Clean internal region size (bit): 3132263

Mon Jul  8 16:11:15 PDT 2019 Sequence number of the redundant LTR-RT library: 4705
The redundant LTR-RT library size (bit): 13546223

Mon Jul  8 16:11:15 PDT 2019 Module 8: Start to make non-redundant library...

Mon Jul  8 16:12:26 PDT 2019 Final LTR-RT library entries: 1360
Final LTR-RT library size (bit): 4329935

Mon Jul  8 16:12:26 PDT 2019 Total intact LTR-RTs found: 1553
Total intact non-TGCA LTR-RTs found: 63

Mon Jul  8 16:12:30 PDT 2019 Start to annotate whole-genome LTR-RTs...
Use -noanno if you don't want whole-genome LTR-RT annotation.


######################################
### LTR Assembly Index (LAI) beta3.1 ###
######################################

Developer: Shujun Ou

Please cite:

Ou S., Chen J. and Jiang N. (2018). Assessing genome assembly quality using the LTR Assembly Index (LAI). Nucleic Acids Res. gky730: https://doi.org/10.1093/nar/gky730

Parameters: -genome /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod -intact /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.pass.list -all /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.out -t 10 -q -blast /global/homes/l/llei2019/.conda/envs/LTR_retriever/bin/


Mon Jul  8 17:16:48 PDT 2019 Dependency checking: Passed!
Mon Jul  8 17:16:48 PDT 2019 Calculation of LAI will be based on the whole genome.
Please use the -mono parameter if your genome is a recent ployploid, for high identity between homeologues will overcorrect raw LAI scores.
Mon Jul  8 17:16:48 PDT 2019 Estimate the identity of LTR sequences in the genome: quick mode
Mon Jul  8 17:22:11 PDT 2019 The identity of LTR sequences: 96.4916337767503%

【Warning】 The identity jumps above the safe limit. Instead, identity of 96% will be used for LAI adjustment.

Mon Jul  8 17:22:11 PDT 2019 Calculate LAI:

Done!

Mon Jul  8 17:22:18 PDT 2019 Result file: /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.out.LAI

You may use either raw_LAI or LAI for intraspecific comparison
but please use ONLY LAI for interspecific comparison

Mon Jul  8 17:22:18 PDT 2019 All analyses were finished!

##############################
####### Result files #########
##############################

Table output for intact LTR-RTs (detailed info)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.pass.list (All LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.nmtf.pass.list (Non-TGCA LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.pass.list.gff3 (GFF3 format for intact LTR-RTs)

LTR-RT library
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.LTRlib.redundant.fa (All LTR-RTs with redundancy)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.LTRlib.fa (All non-redundant LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.nmtf.LTRlib.fa (Non-TGCA LTR-RTs)

Whole-genome LTR-RT annotation by the non-redundant library
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.out.gff (GFF format)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.out.fam.size.list (LTR family summary)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.out.superfam.size.list (LTR superfamily summary)

LTR Assembly Index (LAI)
/global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.out.LAI

---


### 4) I can sumerize TE and repetitive 

```
./extract_family.pl /global/cscratch1/sd/llei2019/genomes/Bsylvaticum_490_v1.0.fa.mod.out.fam.size.list >Bsylvaticum_490_v1.ltr.prop
```