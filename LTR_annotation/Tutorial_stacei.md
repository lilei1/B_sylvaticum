# Bstacei

steps for denova annotate the LTR:

### 1) Run LTR-finder (Download from [here](https://github.com/oushujun/LTR_FINDER_parallel))

The unmasked genome is in this directory:
`/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bstacei_316_v1.0.fa`

```
 perl LTR_FINDER_parallel -seq /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa -size 5000000 -threads 36 -time 1500 -v -harvest_out

```



### 2) Run LTR-harvest

Use this version `gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz`

```
./gt suffixerator -db /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa -indexname Bstacei_316_v1.0.fa -tis -suf -lcp -des -ssp -sds -dna

./gt ltrharvest -index Bstacei_316_v1.0.fa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 >Bstacei_316_v1.0.harvest.scn



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
./LTR_retriever -genome /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bstacei_316_v1.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Bstacei_316_v1.0.fa.finder.combine.scn -v -threads 15

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

Parameters: -genome /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bstacei_316_v1.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Bstacei_316_v1.0.fa.finder.combine.scn -v -threads 15


Fri Jul 12 16:16:34 PDT 2019 Dependency checking: All passed!
Fri Jul 12 16:16:55 PDT 2019 LTR_retriever is starting from the Init step.
Fri Jul 12 16:16:55 PDT 2019 Start to convert inputs...
Total candidates: 790
Total uniq candidates: 790

Fri Jul 12 16:16:57 PDT 2019 Module 1: Start to clean up candidates...
Sequences with 10 missing bp or 0.8 missing data rate will be discarded.
Sequences containing tandem repeats will be discarded.

Fri Jul 12 16:17:05 PDT 2019 599 clean candidates remained

Fri Jul 12 16:17:05 PDT 2019 Modules 2-5: Start to analyze the structure of candidates...
The terminal motif, TSD, boundary, orientation, age, and superfamily will be identified in this step.

Fri Jul 12 16:17:51 PDT 2019 Intact LTR-RT found: 127

Fri Jul 12 16:17:54 PDT 2019 Module 6: Start to analyze truncated LTR-RTs...
Truncated LTR-RTs without the intact version will be retained in the LTR-RT library.
Use -notrunc if you don't want to keep them.

Fri Jul 12 16:17:54 PDT 2019 85 truncated LTR-RTs found
Fri Jul 12 16:18:18 PDT 2019 60 truncated LTR sequences have added to the library

Fri Jul 12 16:18:18 PDT 2019 Module 5: Start to remove DNA TE and LINE transposases, and remove plant protein sequences...
Total library sequences: 278
Fri Jul 12 16:19:24 PDT 2019 Retained clean sequence: 278

Fri Jul 12 16:19:24 PDT 2019 Sequence clustering for /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.ltrTE ...
Fri Jul 12 16:19:24 PDT 2019 Unique lib sequence: 275

Fri Jul 12 16:19:25 PDT 2019 Module 6: Start to remove nested insertions in internal regions...
Fri Jul 12 16:19:43 PDT 2019 Raw internal region size (bit): 694460
Clean internal region size (bit): 683727

Fri Jul 12 16:19:43 PDT 2019 Sequence number of the redundant LTR-RT library: 441
The redundant LTR-RT library size (bit): 834863

Fri Jul 12 16:19:43 PDT 2019 Module 8: Start to make non-redundant library...

Fri Jul 12 16:19:45 PDT 2019 Final LTR-RT library entries: 274
Final LTR-RT library size (bit): 751617

Fri Jul 12 16:19:45 PDT 2019 Total intact LTR-RTs found: 127
Total intact non-TGCA LTR-RTs found: 12

Fri Jul 12 16:19:46 PDT 2019 Start to annotate whole-genome LTR-RTs...
Use -noanno if you don't want whole-genome LTR-RT annotation.


######################################
### LTR Assembly Index (LAI) beta3.1 ###
######################################

Developer: Shujun Ou

Please cite:

Ou S., Chen J. and Jiang N. (2018). Assessing genome assembly quality using the LTR Assembly Index (LAI). Nucleic Acids Res. gky730: https://doi.org/10.1093/nar/gky730

Parameters: -genome /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa -intact /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.pass.list -all /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.out -t 15 -q -blast /global/homes/l/llei2019/.conda/envs/LTR_retriever/bin/


Fri Jul 12 16:30:56 PDT 2019 Dependency checking: Passed!
Fri Jul 12 16:30:56 PDT 2019 Calculation of LAI will be based on the whole genome.
Please use the -mono parameter if your genome is a recent ployploid, for high identity between homeologues will overcorrect raw LAI scores.
Fri Jul 12 16:30:56 PDT 2019 Estimate the identity of LTR sequences in the genome: quick mode
Fri Jul 12 16:31:24 PDT 2019 The identity of LTR sequences: 94.0940248661505%
Fri Jul 12 16:31:24 PDT 2019 Calculate LAI:

Done!

Fri Jul 12 16:31:28 PDT 2019 Result file: /global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.out.LAI

You may use either raw_LAI or LAI for intraspecific comparison
but please use ONLY LAI for interspecific comparison

Fri Jul 12 16:31:28 PDT 2019 All analyses were finished!

##############################
####### Result files #########
##############################

Table output for intact LTR-RTs (detailed info)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.pass.list (All LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.nmtf.pass.list (Non-TGCA LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.pass.list.gff3 (GFF3 format for intact LTR-RTs)

LTR-RT library
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.LTRlib.redundant.fa (All LTR-RTs with redundancy)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.LTRlib.fa (All non-redundant LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.nmtf.LTRlib.fa (Non-TGCA LTR-RTs)

Whole-genome LTR-RT annotation by the non-redundant library
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.out.gff (GFF format)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.out.fam.size.list (LTR family summary)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.out.superfam.size.list (LTR superfamily summary)

LTR Assembly Index (LAI)
/global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.out.LAI

(LTR_retriever) 

---


### 4) I can sumerize TE and repetitive 

```
./extract_family.pl /global/cscratch1/sd/llei2019/genomes/Bstacei//global/cscratch1/sd/llei2019/genomes/Bstacei/Bstacei_316_v1.0.fa.out.fam.size.list >Bstacei_316_v1.0.fa.ltr.prop

```