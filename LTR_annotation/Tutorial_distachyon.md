# Bdistachyon

steps for denova annotate the LTR:

### 1) Run LTR-finder (Download from [here](https://github.com/oushujun/LTR_FINDER_parallel))

The unmasked genome is in this directory:
`/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa`

```
perl LTR_FINDER_parallel -seq /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa -size 5000000 -threads 36 -time 1500 -v -harvest_out

```



### 2) Run LTR-harvest

Use this version `gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz`

```
./gt suffixerator -db /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa -indexname Bdistachyon_314_v3.0.fa -tis -suf -lcp -des -ssp -sds -dna

./gt ltrharvest -index Bdistachyon_314_v3.0.fa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 >Bdistachyon_314_v3.0.harvest.scn


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
./LTR_retriever -genome /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bdistachyon_314_v3.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Bdistachyon_314_v3.0.fa.finder.combine.scn -v -threads 10

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

Parameters: -genome /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bdistachyon_314_v3.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Bdistachyon_314_v3.0.fa.finder.combine.scn -v -threads 10


Fri Jul 12 14:59:07 PDT 2019 Dependency checking: All passed!
Fri Jul 12 14:59:31 PDT 2019 LTR_retriever is starting from the Init step.
Fri Jul 12 14:59:31 PDT 2019 The longest sequence ID in the genome contains 38 characters, which is longer than the limit (15)
Trying to reformat seq IDs...
Attempt 1...
Fri Jul 12 14:59:34 PDT 2019 Seq ID conversion successful!

Fri Jul 12 14:59:34 PDT 2019 Start to convert inputs...
Total candidates: 1221
Total uniq candidates: 1221

Fri Jul 12 14:59:36 PDT 2019 Module 1: Start to clean up candidates...
Sequences with 10 missing bp or 0.8 missing data rate will be discarded.
Sequences containing tandem repeats will be discarded.

Fri Jul 12 14:59:56 PDT 2019 1179 clean candidates remained

Fri Jul 12 14:59:56 PDT 2019 Modules 2-5: Start to analyze the structure of candidates...
The terminal motif, TSD, boundary, orientation, age, and superfamily will be identified in this step.

Fri Jul 12 15:01:13 PDT 2019 Intact LTR-RT found: 313

Fri Jul 12 15:01:23 PDT 2019 Module 6: Start to analyze truncated LTR-RTs...
Truncated LTR-RTs without the intact version will be retained in the LTR-RT library.
Use -notrunc if you don't want to keep them.

Fri Jul 12 15:01:23 PDT 2019 239 truncated LTR-RTs found
Fri Jul 12 15:02:03 PDT 2019 66 truncated LTR sequences have added to the library

Fri Jul 12 15:02:03 PDT 2019 Module 5: Start to remove DNA TE and LINE transposases, and remove plant protein sequences...
Total library sequences: 458
Fri Jul 12 15:05:10 PDT 2019 Retained clean sequence: 457

Fri Jul 12 15:05:10 PDT 2019 Sequence clustering for /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.ltrTE ...
Fri Jul 12 15:05:10 PDT 2019 Unique lib sequence: 456

Fri Jul 12 15:05:15 PDT 2019 Module 6: Start to remove nested insertions in internal regions...
Fri Jul 12 15:05:44 PDT 2019 Raw internal region size (bit): 1453749
Clean internal region size (bit): 1310542

Fri Jul 12 15:05:53 PDT 2019 Sequence number of the redundant LTR-RT library: 1005
The redundant LTR-RT library size (bit): 2392180

Fri Jul 12 15:05:53 PDT 2019 Module 8: Start to make non-redundant library...

Fri Jul 12 15:05:58 PDT 2019 Final LTR-RT library entries: 438
Final LTR-RT library size (bit): 1442091

Fri Jul 12 15:05:59 PDT 2019 Total intact LTR-RTs found: 313
Total intact non-TGCA LTR-RTs found: 26

Fri Jul 12 15:06:01 PDT 2019 Start to annotate whole-genome LTR-RTs...
Use -noanno if you don't want whole-genome LTR-RT annotation.


######################################
### LTR Assembly Index (LAI) beta3.1 ###
######################################

Developer: Shujun Ou

Please cite:

Ou S., Chen J. and Jiang N. (2018). Assessing genome assembly quality using the LTR Assembly Index (LAI). Nucleic Acids Res. gky730: https://doi.org/10.1093/nar/gky730

Parameters: -genome /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod -intact /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.pass.list -all /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.out -t 10 -q -blast /global/homes/l/llei2019/.conda/envs/LTR_retriever/bin/


Fri Jul 12 15:27:09 PDT 2019 Dependency checking: Passed!
Fri Jul 12 15:27:09 PDT 2019 Calculation of LAI will be based on the whole genome.
Please use the -mono parameter if your genome is a recent ployploid, for high identity between homeologues will overcorrect raw LAI scores.
Fri Jul 12 15:27:09 PDT 2019 Estimate the identity of LTR sequences in the genome: quick mode
Fri Jul 12 15:28:53 PDT 2019 The identity of LTR sequences: 93.2001230152184%
Fri Jul 12 15:28:53 PDT 2019 Calculate LAI:

Done!

Fri Jul 12 15:28:59 PDT 2019 Result file: /global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.out.LAI

You may use either raw_LAI or LAI for intraspecific comparison
but please use ONLY LAI for interspecific comparison

Fri Jul 12 15:28:59 PDT 2019 All analyses were finished!

##############################
####### Result files #########
##############################

Table output for intact LTR-RTs (detailed info)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.pass.list (All LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.nmtf.pass.list (Non-TGCA LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.pass.list.gff3 (GFF3 format for intact LTR-RTs)

LTR-RT library
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.LTRlib.redundant.fa (All LTR-RTs with redundancy)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.LTRlib.fa (All non-redundant LTR-RTs)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.nmtf.LTRlib.fa (Non-TGCA LTR-RTs)

Whole-genome LTR-RT annotation by the non-redundant library
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.out.gff (GFF format)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.out.fam.size.list (LTR family summary)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.out.superfam.size.list (LTR superfamily summary)

LTR Assembly Index (LAI)
/global/cscratch1/sd/llei2019/genomes/Bdistachyon/Bdistachyon_314_v3.0.fa.mod.out.LAI

(LTR_retriever) 

---


### 4) I can sumerize TE and repetitive 

```
./extract_family.pl /global/cscratch1/sd/llei2019/genomes/Bdistachyon_314_v3.0.fa.mod.out.fam.size.list >Bdistachyon_314_v3.0.fa.ltr.prop

```