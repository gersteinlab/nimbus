### NIMBus README ###

NIMBus is implemented across two main scripts:

1) covariate.regression.R
2) pval.calculation.R (with pval_cal_functions.R and load_data_functions.R as dependencies)

The first script, covariate.regression.R, uses a set of covariates and mutation rate information to generate a negative bionomial mutation rate model.

covariates.regression.R requires as inputs:

— A list of all 64 trimer combinations (3mer.txt)
— Covariate values for each bin of the genome. In our analysis, the genome was binned into 1MB regions, for 2521 total bins. (files.gersteinlab.org/public-docs/2019/06.15/chrY.rm.cov.merge.hg19.1mb.txt)
— The number variants affecting a particular trimer in each bin of the genome (e.g ./CaType.AAA.matrix.txt shows simulated variant information for variants affecting the AAA trimer in a 5-patient cohort. 64 such files are necessary as input.)
— The total number of trimers of each type in each bin in the genome (./trimer_offset/$trimer$.1Mb.txt provides this information for a 1Mb bin size)

covariates.regression.R provides as output:

— Mu parameter files for simulating a Poisson distribution (./$trimer$.$CaType$.poisson.mu.txt).
— Mu and theta parameter files for simulating a negative binomial distribution (./$trimer$.$CaType$.nb.theta.txt ./$trimer$.$CaType$.mu.nb.txt)

	The poisson, mu, and theta parameter files generated from pval.calculation.R. These parameter files must include the parameters for all trimers in a single tab separated file. Each row in these files corresponds to a separate bin. Single-line/bin examples are shown in (CaType.nb.theta.txt, CaType.mu.nb.txt, CaType.poisson.mu.txt).

	Parameter files must be formatted to be used by the pval.calculation.R script. The example formatted files are shown in the formatted.CaType.mu.nb.txt, formatted.CaType.nb.theta.txt, and formatted.CaType.poisson.mu.txt files. 

pval.calculation.R requires as inputs:
- The results from the covariates.regression.R script

- A CaType.kmer.varct.txt file
	This file must be tab delimited. Column 1: chromosome number, Column 2: start of region, Column 3: end of region, Columns 4: associated gene, Column 5-68: the occurence of the particular mutated trimer in that region. Trimers are referenced in alphabetical order. 

- A raw.kmer.count.txt file
	This file must be tab delimited. Column 1: gene, Column 2: length of region, Column 3-66: the occurence of each trimer in that region. Trimers are referenced in alphabetical order

- A closest file 
	This file must be tab delimited. This file contains the closest 1m genome bins to the test regions. 
	i.e. chr1	14000	15000	gene:1	chr1	1000000	2000000	chr1_2

The CaType.kmer.varct.txt file, raw.kmer.count.txt file, and the closest file must be located in the same folder. In each file, they should share the first four columns corresponding to the test region sites. 

Example files are shown on the Github as CaType.kmer.varct.txt, raw.kmer.count.txt, and closest.1m.union.txt files. 

The test regions folder contains the relevant test regions we used during our analysis. These regions were the intersection of the DHS hotspots and promoters that were also intersected with the transcription factor binding sites (TFBSs). 
	Column 4: transcription factor 
	Column 5: the TF binding site 
	Column 6: gene 

The covariates folder contains the 1M covariate matrix, 1M row names, and the top 30 covariates from the PCA analysis. 
