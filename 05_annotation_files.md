05\_annotation\_files
================

``` r
knitr::opts_chunk$set(echo = TRUE)
library(biomaRt)
```

Retrieve GO terms from Ensembl
------------------------------

``` r
mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=90)

go_ids<-getBM(attributes = c("name_1006","go_id","namespace_1003"), mart = mart)
gene_ids<-getBM(attributes = c("ensembl_gene_id","external_gene_name","description","gene_biotype","name_1006","go_id"), mart = mart)

save(gene_ids,go_ids,file = "./raw_data/annotation.rds")
```

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] biomaRt_2.36.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.18         AnnotationDbi_1.42.1 knitr_1.20          
    ##  [4] magrittr_1.5         hms_0.4.2            progress_1.2.0      
    ##  [7] IRanges_2.14.10      BiocGenerics_0.26.0  bit_1.1-14          
    ## [10] R6_2.2.2             rlang_0.2.1          httr_1.3.1          
    ## [13] stringr_1.3.1        blob_1.1.1           tools_3.5.0         
    ## [16] parallel_3.5.0       packrat_0.4.9-3      Biobase_2.40.0      
    ## [19] DBI_1.0.0            htmltools_0.3.6      assertthat_0.2.0    
    ## [22] yaml_2.2.0           bit64_0.9-7          rprojroot_1.3-2     
    ## [25] digest_0.6.15        crayon_1.3.4         S4Vectors_0.18.3    
    ## [28] bitops_1.0-6         curl_3.2             RCurl_1.95-4.11     
    ## [31] memoise_1.1.0        evaluate_0.11        RSQLite_2.1.1       
    ## [34] rmarkdown_1.10       stringi_1.2.4        compiler_3.5.0      
    ## [37] prettyunits_1.0.2    backports_1.1.2      stats4_3.5.0        
    ## [40] XML_3.98-1.12        pkgconfig_2.0.1
