# Objective: To find the intersection of the orthologs in different intersections


## Inputs: John Lovell create a bunch of the files for orthougs and I need to organized them.

```
#creat a file list
find $PWD -type f -name "*.csv" >file.list

#count the number of orthologs
for i in $(cat file.list); do  echo $(wc -l $i); done|tr ' ' '\t'>count_gid.txt

#creat the correct file format
paste <(cut -f 2 count_gid.txt|cut -d '/' -f 6|cut -d '.' -f1) <(cut -f 1 count_gid.txt) >gene_ortholog.count

paste <(cut -f 1 gene_ortholog.count |cut -d '_' -f 1) <(cut -f 1 gene_ortholog.count |cut -d '_' -f 2) <(cut -f 1 gene_ortholog.count |cut -d '_' -f 3) <(cut -f 1 gene_ortholog.count |cut -d '_' -f 4) <(cut -f 1 gene_ortholog.count |cut -d '_' -f 5) <(cut -f 2 gene_ortholog.count)>gene_ortholog.count.matrix
```

## Mask the multiple versus 1 or 1 versus multiple as 1, and if absent as 0:

```
./mask_count.pl gene_ortholog.count.matrix>masked_gene_ortholog.count.matrix
```

## Sort the file based on the categories:

```
 sort -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n masked_gene_ortholog.count.matrix >sorted_masked_gene_ortholog.count.matrix
```

## Sum the number of the files
```
./sum_count.pl sorted_masked_gene_ortholog.count.matrix >sorted_masked_gene_ortholog.categories.txt

```

## expand the file and format it as R plot file

```
./expand_matrix.pl sorted_masked_gene_ortholog.categories.txt >sorted_masked_gene_ortholog.matrix
```

## R plot with UpSetR
The R code is in the scripts and the plot is in the results file.