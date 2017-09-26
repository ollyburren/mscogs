# MSCOGS - **M**ulti **S**NP **c**apture-HIC omnibus **g**ene **s**core. 
 
This attempts to port the previous software [COGS](https://github.com/ollyburren/CHIGP) from R programming language to python. There are two main objectives:

1. Get better at python
2. Allow COGS to cope with models consisting of multiple causal variants. Such as those generated by [GUESSFM](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005272).
 
## Step 1

Here we will develop software to interface with GUESSFM R objects and retrieve the relevant information. Hope to complete this activity by the end of day 1

## Step 2

Here we will integrate functional capture Hi-C data and other annotations to compute a genescore and highest scoring PIR,cSNP and promoter. Hope to complete this activity by the end of day 2/3

## Step 3

Compute tree based Bayes Factors and analyse to ascertain the most likely set of genes/tissue to be prioritised for a give trait.

## Step 4

Extend the above to deal with single variant posterior probabilities and hence apply to selected phenotypes for which public GWAS summary statistics are available.

Further information on the methods employed is available from: 

- Burren,O.S. et al. (2017) Chromosome contacts in activated T cells identify autoimmune disease candidate genes. Genome Biol., 18, 165.
- Javierre,B.M. et al. (2016) Lineage-Specific Genome Architecture Links Enhancers and Non-coding Disease Variants to Target Gene Promoters. Cell, 167, 1369–1384.e19.
 
