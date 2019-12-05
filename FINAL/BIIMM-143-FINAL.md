BIMM 143 FINAL
================
Michael Nguyen
12/5/2019

## R Markdown

r

``` r
library(bio3d)
library(BiocManager)
```

    ## Bioconductor version 3.10 (BiocManager 1.30.9), ?BiocManager::install for
    ##   help

``` r
seqs <- read.fasta("unk-sequences.fa")
seqs
```

    ##          1        .         .         .         .         .         .         70 
    ## wt       MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQ
    ## mutant   MTEYKLVVVGAVGVGKSALTINLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQ
    ##          *********** *********^************************************************ 
    ##          1        .         .         .         .         .         .         70 
    ## 
    ##         71        .         .         .         .         .         .         140 
    ## wt       YMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIP
    ## mutant   YMRSGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQVQDLARSYGIP
    ##          ***^******************************************************* ********** 
    ##         71        .         .         .         .         .         .         140 
    ## 
    ##        141        .         .         .         .        189 
    ## wt       FIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGCVKIKKCIIM
    ## mutant   FIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGCVKIKKCIIM
    ##          ************************************************* 
    ##        141        .         .         .         .        189 
    ## 
    ## Call:
    ##   read.fasta(file = "unk-sequences.fa")
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   2 sequence rows; 189 position columns (189 non-gap, 0 gap) 
    ## 
    ## + attr: id, ali, call

``` r
ide <- conserv(seqs$ali, method = "identity")
mutant.sites <- which(ide < 1)   # 1 being conserved!!
mutant.sites
```

    ## [1]  12  22  74 130
