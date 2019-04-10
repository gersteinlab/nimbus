# NIMBUS
### A Negative Binomial Regression based Integrative Method for Mutation Burden Analysis
Identifying highly mutated regions is a key way that scientists can use on population scale sequencing to discover key genomic regions associated with complex diseases such as cancer. Nevertheless, it is challenging to identify such regions because severe mutation rate heterogeneity, across different genome regions of the same individual and also across different individuals, gives rise to highly over-dispersed counts of mutations. Moreover, it is known that part of this heterogeneity relates to confounding genomic features, such as replication timing and chromatin organization. Here, we address these issues with a Negative binomial regression based Integrative Method for mutation Burden analysis (NIMBus). This approach uses a Gamma-Poisson mixture model to capture the mutation rate heterogeneity across different individuals and thus models the over dispersed mutation counts by a negative binomial distribution. Furthermore, it regresses the mutation counts against 381 features extracted from REMC and ENCODE to accurately estimate the local background mutation rate. This framework can be readily extended to accommodate additional genomic features in the future. NIMBus was used to analyze 649 whole-genome cancer sequences. It successfully controlled P value inflation and identified well-known coding and noncoding drivers, such as TP53 and the TERT promoter. We make NIMBus available and release our results as an online resource (nimbus.gersteinlab.org).
  
### [Results](https://github.com/gersteinlab/nimbus)
#### Results for Noncoding Regions
[P-Values for Individual Cancers and Annotations](https://github.com/gersteinlab/nimbus/tree/master/noncoding-individual) 

[P-Values for Fisher's Method](https://github.com/gersteinlab/nimbus/tree/master/noncoding-fisher)
#### Results for Coding Regions 
[P-Values for Fisher's Method](https://github.com/gersteinlab/nimbus/tree/master/coding-fisher)

