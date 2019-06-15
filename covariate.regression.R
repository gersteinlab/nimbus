#!/usr/bin/env Rscript
library("base")
library("AER")
library("foreign")
library("ggplot2")
library("MASS")
args = commandArgs(trailingOnly=TRUE)
cancer <- args[1]
trimer <- args[2]

# e.g.,
# cancer <- "Breast-AdenoCa"
# trimer <- "AAA"

dir <- "./regression"
cancer_dir <- "./trimer_variants_matrix"
trimer_list <- read.table("./3mer.txt",as.is=T)$V1
# top 30 pcs, after running pca
pcs_cov <- read.table("./1Mb_pca_covariates_top30", as.is = TRUE, header = F)
# off set - keeps track of the number of times each trimer appears
kmer_length <- read.table(paste("./", trimer,".1Mb.txt", sep = ""), as.is = TRUE, header = F)$V1

# get the offset
len.offset <- kmer_length

# get the variant count
file = paste(cancer_dir, "/", cancer, "_3mer_matrix/", cancer, ".", trimer, ".matrix.txt" ,sep="")
y <- read.table(file, as.is=T)
y <- rowSums(y)

######## MODEL 1: POISSON ########
######## 1 best feature from all catogories of feature, #######
poisson.1 <- glm(y ~ .+offset(log(len.offset)), data=data.frame(pcs_cov), family = poisson())

poisson_mu <- (poisson.1$fitted.values / len.offset)
disptest_results <- dispersiontest(poisson.1)$p.value

######## MODEL 2: Negative binomial ########
nb.1 <- glm.nb(y ~ .+offset(log(len.offset)), data=data.frame(pcs_cov))

nb_mu <- (nb.1$fitted.values / len.offset)
nb_theta <- (nb.1$theta)

cmd <- "mkdir -p "
data_dir <- "./regression/results/"
system(paste(cmd, data_dir, cancer, sep = ""))

write.table(poisson_mu, file = paste(data_dir, cancer, "/", trimer, ".", cancer, ".poisson.mu.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(disptest_results, file = paste(data_dir, cancer, "/", trimer, ".", cancer, ".disp.pval.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(nb_mu, file = paste(data_dir, cancer, "/", trimer, ".", cancer, ".nb.mu.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(nb_theta, file = paste(data_dir, cancer, "/", trimer, ".", cancer, ".nb.theta.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

print(cancer)
print(trimer)