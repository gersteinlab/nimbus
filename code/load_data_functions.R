#!/usr/bin/env Rscript

common_file_define <- function(prefix.input,resolution.train, work.dir, cancer.type){
  file.kmer <- "./trimer.txt"
  file.cancer.list <- "./cancer_list.txt"
  
  #/*======= path.train has different data structure of the local code =======*/
  
  # all directories must be followed by a a forward slash (i.e. ..../)

  # this is the directory with the poisson, mu, and theta, parameter values outputted from the covariate.regresion.R script
  # formatted into "AML.mu.nb.txt" type formats
  path.train=paste(work.dir, "/covariate_regression_output/", prefix.input, "/", cancer.type, "/formatted/", sep="")

  # this is the directory containing a KmerCt file for each type of variant. In our case, we have a file for each 
  # type of cancer. The first four columns are taken directly from the cell type file. The fifth column
  # corresponds to the number of variants affecting each the 64 trimers in the test region(s) with comma separation.
  # this directory should also contain the closest.1m.union.txt or closest.1m.intersect.txt file 
  # should also have the raw.kmer.count.txt file  
  
  
  path.varCt = paste(work.dir, "/KmerCt/",prefix.input,"/",sep="")
  path.kmerCt = paste(work.dir, "/KmerCt/", prefix.input, "/",sep="")

  path.save=paste(work.dir, "/pval/" , prefix.input,"/",sep="")
  
  if (!file.exists(path.save)){
    dir.create(path.save)
  }
  if(!file.exists(path.varCt)){
    print(path.varCt)
    stop("common_file_define: path.varCt does not exist")
  }
  if(!file.exists(path.kmerCt)){
    print(path.kmerCt)
    stop("common_file_define: path.kmerCt does not exist")
  }
  if(!file.exists(path.train)){
    print(path.train)
    stop("common_file_define: path.train does not exist")
  }
  if(!file.exists(path.save)){
    print(path.save)
    stop("common_file_define: path.save does not exist")
  }
  return(list(file.kmer=file.kmer, file.cancer.list=file.cancer.list, path.train=path.train, 
              path.varCt=path.varCt, path.kmerCt=path.kmerCt,path.save=path.save))
}

file_name <- function(prefix.input,path.v,resolution.train,cancer.type){
  if(nchar(cancer.type) < 1){
    stop("common_file_define: invalid cancer.type")
  }
  if(nchar(resolution.train) < 1){
    stop("common_file_define: invalid resolution.train")
  }
  if(nchar(prefix.input) < 1){
    stop("common_file_define: invalid prefix.input")
  }
  #/*==== find the closest bin mu ===*/
  f.test.closest <- paste(path.v$path.varCt,"closest.",resolution.train,".",prefix.input,sep="")
  if(! file.exists(f.test.closest)){
    print(f.test.closest)
    stop("common_file_define: invalid f.test.closest")
  }
  #/*==== find pois mu est  ===*/
  f.train.pois.mu <- paste(path.v$path.train,cancer.type,".poisson.mu.txt",sep="")
  if(! file.exists(f.train.pois.mu)){
    print(f.train.pois.mu)
    stop("common_file_define: invalid f.train.pois.mu")
  }
  #/*==== find ng mu est ===*/
  f.train.nb.mu <- paste(path.v$path.train,cancer.type,".mu.nb.txt",sep="")
  if(! file.exists(f.train.nb.mu)){
    print(f.train.nb.mu)
    stop("common_file_define: invalid f.train.nb.mu")
  }
  
  #/*==== find ng other parameters ===*/
  f.train.nb.theta <- paste(path.v$path.train,cancer.type,".nb.theta.txt",sep="")
  if(! file.exists(f.train.nb.theta)){
    print(f.train.nb.theta)
    stop("common_file_define: invalid f.train.nb.theta")
  } 
  
  #/*==== find raw kmer count ===*/
  f.test.kmer.raw <- paste(path.v$path.kmerCt,"/","raw.kmer.count.txt",sep="")
  if(! file.exists(f.test.kmer.raw)){
    print(f.test.kmer.raw)
    stop("common_file_define: invalid f.test.kmer.raw")
  }
  
  #/*==== find var count ===*/
  f.test.var.ct <- paste(path.v$path.varCt,"/",cancer.type,".kmer.varct.txt",sep="")
  if(! file.exists(f.test.var.ct)){
    print(f.test.var.ct)
    stop("common_file_define: invalid f.test.var.ct")
  }
  
  #/*==== find all kmer counts ===*/
  f.org.kmer <- "/ysm-gpfs/pi/gerstein/yc774/project/nimbus/input/3mer.txt"
  if(! file.exists(f.org.kmer)){
    print(f.org.kmer)
    stop("common_file_define: invalid f.org.kmer")
  }
  
  #/*==== find all cancer types ===*/
  f.org.cancerType <- "/ysm-gpfs/pi/gerstein/yc774/project/nimbus/input/cancer_list.txt"
  if(! file.exists(f.org.cancerType)){
    print(f.org.cancerType)
    stop("common_file_define: invalid f.org.cancerType")
  }
  
  file.v <- list(f.test.closest=f.test.closest, f.train.pois.mu=f.train.pois.mu, 
                 f.train.nb.mu=f.train.nb.mu,f.train.nb.theta=f.train.nb.theta,
                 f.test.kmer.raw=f.test.kmer.raw, f.test.var.ct=f.test.var.ct, f.org.kmer=f.org.kmer, f.org.cancerType=f.org.cancerType)
                 
  return(file.v)
}
separate.aggregate.at.gene.level <- function(id.unique,d.org,kmer.v){
  ####/*get the sum of var/kmer/sample counts for each gene, if necessary*/  
  index.val <- match(kmer.v, colnames(d.org))
  if(sum(is.na(index.val))==0){
    d.mat <- d.org[,index.val]
    colnames(d.mat) <- kmer.v
  }else{
    d.mat <- d.org[,-c(1:4)]
  }
  dtmp <- data.frame(id.unique,d.mat)
  dtmp.agg <- aggregate(. ~ id.unique, data = dtmp, sum)
  return(dtmp.agg)
}


data.prepare.for.testing <- function(file.v,flag.merge.input){
  library(data.table)
  #/*====== read all kmers ======*/
  d.kmerList <- read.table(file.v$f.org.kmer,header = F, sep="\t",as.is = T)
  kmer.v <- d.kmerList$V1
  #/*====== read all cancer types ======*/
  d.cancerList <- read.table(file.v$f.org.cancerType,header = F, sep="\t",as.is = T)
  cancerType.v <- d.cancerList$V1
  
  #/*====== read raw kmer count data ======*/
  d.test.kmer <- fread(file.v$f.test.kmer.raw,header = F, sep="\t")
  d.test.kmer <- data.frame(d.test.kmer)
  colnames(d.test.kmer) <- c("name","size",kmer.v)

  #/*====== read variant kmer count data ======*/
  d.test.var <- fread(file.v$f.test.var.ct,header = F, sep="\t")
  d.test.var <- data.frame(d.test.var)
  d.test.var <- d.test.var[,1:68]
  print(dim(d.test.var))
  print(head(d.test.var))
  colnames(d.test.var) <- c("chr","start","stop","name",kmer.v)

  #/*====== read closest bin ======*/
  d.test.closest <- fread(file.v$f.test.closest,header = F, sep="\t")
  d.test.closest <- data.frame(d.test.closest)
  colnames(d.test.closest) <- c("chr","start","stop","name","bin_target","bin_start","bin_stop","bin_name")
  print("d.test.closest: ")
  print(dim(d.test.closest))
  print(head(d.test.closest))

  #/*====== read pois estimates ======*/
  d.train.pois.mu <- fread(file.v$f.train.pois.mu,header = F, sep="\t")
  d.train.pois.mu <- data.frame(d.train.pois.mu)
  colnames(d.train.pois.mu) <-c("chr","start","stop","name",kmer.v)  
  print("d.train.pois.mu")
  print(dim(d.train.pois.mu))
  print(head(d.train.pois.mu))


  #/*====== read NB mu estimates ======*/
  d.train.nb.mu <- fread(file.v$f.train.nb.mu,header = F, sep="\t")
  d.train.nb.mu <- data.frame(d.train.nb.mu)
  colnames(d.train.nb.mu) <-c("chr","start","stop","name",kmer.v)
  print("d.train.nb.mu: ")
  print(dim(d.train.nb.mu))
  print(head(d.train.nb.mu))
  
  #/*====== read NB theta estimates ======*/
  d.train.nb.theta <- fread(file.v$f.train.nb.theta,header = F, sep="\t")
  d.train.nb.theta <- data.frame(d.train.nb.theta)
  colnames(d.train.nb.theta) <- c("kmer", "theta ")
  print("nb theta:")
  print(dim(d.train.nb.theta))
  print(head(d.train.nb.theta))

  
  #/*====== matching all the data to test ======*/
  name.used <- d.test.kmer$name
  d.test.kmer <- d.test.kmer[match(name.used,d.test.kmer$name),]  
  d.test.var <- d.test.var[match(name.used,d.test.var$name),]  
  d.test.closest <- d.test.closest[match(name.used,d.test.closest$name),] 
  print("d.test.closest: ")
  print(dim(d.test.closest))
  print(head(d.test.closest))
  
  
  if(grepl("ENH.geneMerge",prefix.input)){ 
    index.val <- grepl("\\.",name.used)
    d.test.kmer <- d.test.kmer[!index.val,]
    d.test.var <- d.test.var[!index.val,]
    d.test.closest <- d.test.closest[!index.val,]
    name.used <- d.test.kmer$name
  }
  
  if(flag.merge.input==1){
    tmp <- strsplit(name.used, split = ":")
    geneName <- sapply(tmp,function(x) x[1])
  }else if(flag.merge.input==0){
    geneName <- name.used
  }
  
  #/*====== remove chromosome x y, and genes with NA as names*/
  index.rm <- ( ! is.element(d.test.var$chr, paste("chr",1:22,sep="")) )
  index.rm <- (index.rm | (geneName=="NA" ))
  index.rm <- (index.rm | (geneName=="." ))
  d.test.kmer <- d.test.kmer[!index.rm,]
  d.test.var <- d.test.var[!index.rm,]
  d.test.closest <- d.test.closest[!index.rm,]
  
  name.used <- d.test.kmer$name
  if(flag.merge.input==1){
    tmp <- strsplit(name.used, split = ":")
    geneName <- sapply(tmp,function(x) x[1])
  }else if(flag.merge.input==0){
    geneName <- name.used
  }
  
  target.bin.name.used = d.test.closest$bin_name

  id.unique <- paste(geneName,target.bin.name.used,sep = ":")
  
  d.test.kmer.merge <- separate.aggregate.at.gene.level(id.unique,d.test.kmer,kmer.v)
  d.test.var.merge <- separate.aggregate.at.gene.level(id.unique,d.test.var,kmer.v)
  tmp <- strsplit(as.character(d.test.kmer.merge$id.unique), split = ":")

  bin_target_merge <- sapply(tmp,function(x) x[2])

  d.reform.mu.pois.merge <- d.train.pois.mu[match(bin_target_merge,d.train.pois.mu$name),]

  d.reform.mu.nb.merge <- d.train.nb.mu[match(bin_target_merge,d.train.nb.mu$name),]
  bin_target_merge[is.element(bin_target_merge, d.train.pois.mu$name)][1]
  
  dtest.prepare <- list(kmer.v=kmer.v, cancerType.v=cancerType.v,
                        d.test.kmer=d.test.kmer, d.test.var=d.test.var,#d.test.sample=d.test.sample,
                        d.test.closest=d.test.closest,d.train.pois.mu=d.train.pois.mu,
                        d.train.nb.mu=d.train.nb.mu,d.train.nb.theta=d.train.nb.theta,
                        d.reform.mu.pois.merge=d.reform.mu.pois.merge,
                        d.reform.mu.nb.merge=d.reform.mu.nb.merge,
                        d.test.kmer.merge=d.test.kmer.merge,d.test.var.merge=d.test.var.merge,
                        #d.test.sample.merge=d.test.sample.merge,
                        bin_target_merge=bin_target_merge,id.unique=id.unique)
  return(dtest.prepare)
}

sanity_check_of_read_data <- function(dtest.prepare){
  print("===== d.test.kmer.merge =====");
  print(dim(dtest.prepare$d.test.kmer.merge))
  print(dim(dtest.prepare$d.test.kmer))
  print(dtest.prepare$d.test.kmer.merge[1:3,1:6])
  print(dtest.prepare$d.test.kmer[1:3,1:6])
  
  print("===== d.test.var.merge =====");
  print(dtest.prepare$d.test.var.merge[1:3,1:6])
  print(dtest.prepare$d.test.var[1:3,1:6])
  print(dim(dtest.prepare$d.test.var.merge))
  print(dim(dtest.prepare$d.test.var))
  
  print("===== d.train mu estimates =====");
  print(dtest.prepare$d.reform.mu.pois.merge[1:3,1:6])
  print(dtest.prepare$d.reform.mu.nb.merge[1:3,1:6])
  print(dim(dtest.prepare$d.reform.mu.pois.merge))
  print(dim(dtest.prepare$d.reform.mu.nb.merge))
  return(1)
}


save.results.for.pval <- function(path.v, prefix.input,cancer.type,data.pval){
  file.pval <- paste(path.v$path.save,"/pval.",cancer.type,".",prefix.input,".txt",sep = "")
  write.table(data.pval, file = file.pval, append = F, quote = F, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)  
}
