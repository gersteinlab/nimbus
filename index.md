# NIMBUS
### A Negative Binomial Regression based Integrative Method for Mutation Burden Analysis
Identifying highly mutated regions is a key way that scientists can use on population scale sequencing to discover key genomic regions associated with complex diseases such as cancer. Nevertheless, it is challenging to identify such regions because severe mutation rate heterogeneity, across different genome regions of the same individual and also across different individuals, gives rise to highly over-dispersed counts of mutations. Moreover, it is known that part of this heterogeneity relates to confounding genomic features, such as replication timing and chromatin organization. Here, we address these issues with a Negative binomial regression based Integrative Method for mutation Burden analysis (NIMBus). This approach uses a Gamma-Poisson mixture model to capture the mutation rate heterogeneity across different individuals and thus models the over dispersed mutation counts by a negative binomial distribution. Furthermore, it regresses the mutation counts against 381 features extracted from REMC and ENCODE to accurately estimate the local background mutation rate. This framework can be readily extended to accommodate additional genomic features in the future. NIMBus was used to analyze 649 whole-genome cancer sequences. It successfully controlled P value inflation and identified well-known coding and noncoding drivers, such as TP53 and the TERT promoter. The NIMBus model and results are available as an online resource (github.gersteinlab.org/nimbus).

#### Full Results and README
[Download All Results (81MB)](http://files.gersteinlab.org/public-docs/2019/06.15/all-results.zip)

[README file](http://files.gersteinlab.org/public-docs/2019/06.15/README.txt)
#### Results for Noncoding Regions
[P-Values for Individual Cancers and Annotations](http://files.gersteinlab.org/public-docs/2019/06.15/noncoding-individual.zip) 

[P-Values for Fisher's Method](http://files.gersteinlab.org/public-docs/2019/06.15/noncoding-fisher.zip)

[P-Values for Promoter-TF burdening](http://files.gersteinlab.org/public-docs/2019/06.15/promoter-TF-burden.zip)
#### Results for Coding Regions 
[P-Values for Fisher's Method](http://files.gersteinlab.org/public-docs/2019/06.15/coding-fisher.zip)
#### Model Covariates
The following covariate files may be used with the NIMBus source code on GitHub to generate mutation rate models.

[Covariates (7.8MB)](http://files.gersteinlab.org/public-docs/2019/06.15/chrY.rm.cov.merge.hg19.1mb.txt)

[ENCODE histone mark covariates (30 principle components) (1.3MB)](http://files.gersteinlab.org/public-docs/2019/06.15/1Mb_pca_covariates_top30)
