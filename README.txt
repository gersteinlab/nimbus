### NIMBus README ###

NIMBus is implemented across two main scripts:

1) covariate.regression.R
2) pval.calculation.R (with convolution.methods.nb.pois.R as a dependency)

The first script, covariate.regression.R, uses a set of covariates and mutation rate information to generate a negative bionomial mutation rate model.

covariates.regression.R requires as inputs:

— A list of all 64 trimer combinations (3mer.txt)
— Covariate values for each bin of the genome. In our analysis, the genome was binned into 1MB regions, for 2521 total bins. (files.gersteinlab.org/public-docs/2019/06.15/chrY.rm.cov.merge.hg19.1mb.txt)
— The number variants affecting a particular trimer in each bin of the genome (e.g ./CaType.AAA.matrix.txt shows simulated variant information for variants affecting the AAA trimer in a 5-patient cohort. 64 such files are necessary as input.)
— The total number of trimers of each type in each bin in the genome (./trimer_offset/$trimer$.1Mb.txt provides this information for a 1Mb bin size)

covariates.regression.R provides as output:

— Mu parameter files for simulating a Poisson distribution (./$trimer$.$CaType$.poisson.mu.txt).
— Mu and theta parameter files for simulating a negative binomial distribution (./$trimer$.$CaType$.nb.theta.txt ./$trimer$.$CaType$.mu.nb.txt)

pval.calculation.R requires as inputs:


— The poisson, mu, and theta parameter files generated from pval.calculation.R. These parameter files must include the parameters for all trimers in a single tab separated file. Each row in these files corresponds to a separate bin. Single-line/bin examples are shown in (CaType.nb.theta.txt, CaType.mu.nb.txt, CaType.poisson.mu.txt).
— A test region file (CaType.input.txt).

The test region file contains the following columns:

(1) test region name (e.g., gene).
(2) total region length in nucleotides.
(3) total number of variants affecting this test region across the disease cohort.
(4) bin name(s). If the test region spans multiple bins, they are listed with comma separation. These bin names must be included in the parameter files (mu, theta, mu poission).
(5) length of the test region overlapping the bin(s). If the test region overlaps multiple bins, the length of overlap for each bin is listed with comma separation.
(6) the number of variants affecting the bin(s). If the test region overlaps multiple bins, the number of variants affecting each bin is listed with comma separation.
(7) the number of each the 64 trimers in the test region(s) with comma separation. If the test region overlaps multiple bins, these are listed with colon separation.
(8) the number of variants affecting each the 64 trimers in the test region(s) with comma separation. If the test region overlaps multiple bins, these are listed with colon separation.

An example test region input is given a CaType.input.txt. This file contains two simulated test regions. The first overlaps just a single bin. The second test region overlaps multiple bins.
