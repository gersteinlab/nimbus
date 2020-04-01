#!/usr/bin/env Rscript
rm(list=ls())
setwd("./nimbus")
library(data.table)

print("before loading files")
args=(commandArgs(TRUE))
if(length(args)==0){
    stop("No arguments supplied.")
}

cancer.type <- args[1]
prefix.input <- args[2]
print(paste("cancer.type =",cancer.type))
print(paste("prefix.input =",prefix.input))

#prefix.v <- "gene.annotation.txt"
prefix.v <- args[2]
if(!is.element(prefix.input,prefix.v)){
	stop("prefix.input is invalid")
}

resolution.train = "1m"
if(!is.element(prefix.input, prefix.v)){
  stop("prefix.input is not valid");
}
flag.merge.input <- 1

work.dir <- args[3]

source("./load.data.functions.R")

#/*============ read input files ==========*/
path.v <- common_file_define(prefix.input,resolution.train, work.dir, cancer.type)
file.v <- file_name(prefix.input,path.v,resolution.train,cancer.type)
dtest.prepare <- data.prepare.for.testing(file.v,flag.merge.input)

############## only for testing purpose ##########
sanity_check_of_read_data(dtest.prepare)

#/*============ pval calculation ==========*/
source("./pval.cal.functions.R")
data.pval <- pval_calculation_all_genes(dtest.prepare,flag.merge.input,file.v)
save.results.for.pval(path.v, prefix.input,cancer.type,data.pval)
