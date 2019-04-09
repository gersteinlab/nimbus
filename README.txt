Files are separated into three folders: (1) coding-fisher, (2) noncoding-fisher, (3) noncoding-individual. 

Files in the /coding-fisher/ folder contain at least six columns: (1) gene, (2) len, (3) ct, (4) mu, (5) sigma, (6) region.
The len column indicates the length of the gene. The ct column is the number of variants affecting the test region (e.g. gene 
of interest). Mu and sigma are parameters from the negative binomial model which correspond to the mean and degree of dispersion
of mutation rate, respectively. Region is the bin of the genome where the test region is located. The genome was binned into 1 
megabase bins in order to correct the mutation rate for various genomic features. The remaining columns in the files contain 
the pvalue. 

Files in the /noncoding-fisher/ and /noncoding-individual/ folders have five columns: (1) chr, (2) start, (3) stops, (4) name, 
and (p). Chr means _____.  The name corresonds to the noncoding annotation file used. 

The /coding-fisher/all-cancers folder contains (1) a file with the pvalues for all cancers and (2) cancer-specific file 
containing only the significant genes and pvalues. 

Similarly, the noncoding-fisher folder contains region-specific files in both original and significant folders. Files containing only 
significant pvalues are located within the /noncoding-fisher/significant folder. 

Lastly, in the noncoding-individual folder, within each folder exists (1) a file containing all pvalues for all cancers and 
(2) significant pvalues file filtered by region and by cancer.  
