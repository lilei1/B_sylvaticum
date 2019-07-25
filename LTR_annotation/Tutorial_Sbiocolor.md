# Sbiocolor

steps for denova annotate the LTR:

### 1) Run LTR-finder (Download from [here](https://github.com/oushujun/LTR_FINDER_parallel))

The unmasked genome is in this directory:
`/global/cscratch1/sd/llei2019/genomes/Sbicolor/Sbicolor_454_v3.0.1.fa`

```
perl LTR_FINDER_parallel -seq /global/cscratch1/sd/llei2019/genomes/Sbicolor/Sbicolor_454_v3.0.1.fa -size 5000000 -threads 36 -time 1500 -v -harvest_out

```



### 2) Run LTR-harvest

Use this version `gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz`

```
./gt suffixerator -db /global/cscratch1/sd/llei2019/genomes/Sbicolor/Sbicolor_454_v3.0.1.fa -indexname Sbicolor_454_v3.0.1.fa -tis -suf -lcp -des -ssp -sds -dna

./gt ltrharvest -index Sbicolor_454_v3.0.1.fa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 >Sbicolor_454_v3.0.1.fa.harvest.scn


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
./LTR_retriever -genome /global/cscratch1/sd/llei2019/genomes/Sbicolor/Sbicolor_454_v3.0.1.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bstacei_316_v1.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Sbicolor_454_v3.0.1.fa.finder.combine.scn -v -threads 15


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

Parameters: -genome /global/cscratch1/sd/llei2019/genomes/Sbicolor/Sbicolor_454_v3.0.1.fa -inharvest /global/u2/l/llei2019/softwares/gt-1.5.10-Linux_x86_64-64bit-complete/bin/Bstacei_316_v1.0.harvest.scn -infinder /global/u2/l/llei2019/softwares/LTR_FINDER_parallel/Sbicolor_454_v3.0.1.fa.finder.combine.scn -v -threads 15


Mon Jul 15 12:49:59 PDT 2019 Dependency checking: All passed!
Mon Jul 15 12:50:23 PDT 2019 LTR_retriever is starting from the Init step.
Mon Jul 15 12:50:25 PDT 2019 Start to convert inputs...
Total candidates: 790
Total uniq candidates: 790

Mon Jul 15 12:50:31 PDT 2019 Module 1: Start to clean up candidates...
Sequences with 10 missing bp or 0.8 missing data rate will be discarded.
Sequences containing tandem repeats will be discarded.

Mon Jul 15 12:50:44 PDT 2019 753 clean candidates remained

Mon Jul 15 12:50:44 PDT 2019 Modules 2-5: Start to analyze the structure of candidates...
The terminal motif, TSD, boundary, orientation, age, and superfamily will be identified in this step.

Mon Jul 15 12:51:36 PDT 2019 Intact LTR-RT found: 0

Mon Jul 15 12:51:36 PDT 2019 No LTR-RT was found in your data.

Mon Jul 15 12:51:36 PDT 2019 All analyses were finished!
(LTR_retriever) 

---


### 4) No new LTR annotated!! I did not adapte the files