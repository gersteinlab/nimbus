mu2p <- function(mu.val,sigma.val){
  #/*======== prepare for the Negative binomial regression =========*/
  alfa.val <- 1.0/sigma.val
  p.val <- 1.0/(1.0+sigma.val*mu.val)
  return(list(alfa.val=alfa.val, p.val=p.val))
}


pcal.ng.many.plus.pois <- function(n.obs, alfa.v, p.v, lambda.pois){
  try(if(sum(p.v<=0) > 0) stop("pcal.ng.many: invalid p.v input"))
  if(n.obs == 0){
    p.value = 1
  }else{
    x1 <- seq(0,n.obs)
    p.calculate <- dnbinom(x1, size=alfa.v[1], prob=p.v[1], log = FALSE)
    for(i in 2:length(alfa.v)){
      p.calculate.new <- lapply(x1, pcal.ng.two, alfa=alfa.v[i],p=p.v[i],p.org = p.calculate)
      p.calculate <- unlist(p.calculate.new)
      print(paste(i,p.calculate))
      if(sum(is.na(p.calculate))>0){
        stop("pcal.ng.many.plus.pois: NAN generated", p.calculate)
      }
    }
    #/*---------- add poisson ----------*/
    p.calculate.new <- lapply(x1, pcal.convolute.pois, lambda.pois=lambda.pois,p.org = p.calculate)
    p.value <- 1.0 - sum(unlist(p.calculate.new))
    if(p.value<0){ #/*for numerical calculation*/
      p.value = 0
    }
  }
  return (p.value)
}

pcal.ng.two <- function(n.obs, alfa, p, p.org){
  x2 <- n.obs - seq(0,n.obs)
  pval <- dnbinom(x2, size=alfa, prob=p, log = FALSE)
  return(sum(p.org[1:length(x2)]*pval))
}

pcal.convolute.pois <- function(n.obs, lambda.pois,p.org){
  x2 <- n.obs - seq(0,n.obs)
  pval <- dpois(x2, lambda = lambda.pois, log = FALSE)
  return(sum(p.org[1:length(x2)]*pval))
}

pcal.ng.many <- function(n.obs, alfa.v, p.v){
  try(if(sum(p.v<=0) > 0) stop("pcal.ng.many: invalid p.v input"))
  if(n.obs == 0){
    p.value = 1
  }else{
   if(length(alfa.v)==1){
     p.value <- pnbinom(x1, size=alfa.v[1], prob=p.v[1], log = FALSE,lower.tail=F)
   } else{
      x1 <- seq(0,n.obs)
      p.calculate <- dnbinom(x1, size=alfa.v[1], prob=p.v[1], log = FALSE)
      for(i in 2:length(alfa.v)){
        p.calculate.new <- lapply(x1, pcal.ng.two, alfa=alfa.v[i],p=p.v[i],p.org = p.calculate)
        p.calculate <- unlist(p.calculate.new)
        if(sum(is.na(p.calculate))>0){
          stop("pcal.ng.many: NAN generated,p.calculate=",p.calculate )
        }
      }
      p.value <- (1 - sum(p.calculate))
      if(p.value<0){ #/*for numerical calculation*/
        p.value = 0
      }
   }
  }
  return (p.value)
}

pval_calculation_per_gene <- function(gene.org, dtest.prepare, geneName, index.pois, index.nb){
  index.tmp = which(geneName==gene.org)
  #print("----------------------------------------------------")
  #print(index.tmp)
  var.mat <- dtest.prepare$d.test.var.merge[index.tmp,-1]
  kmer.mat <- dtest.prepare$d.test.kmer.merge[index.tmp,-1]
  mu.pois.mat <- dtest.prepare$d.reform.mu.pois.merge[index.tmp,-c(1:4)]
  mu.nb.mat <- dtest.prepare$d.reform.mu.nb.merge[index.tmp,-c(1:4)]
  theta.v <- dtest.prepare$d.train.nb.theta$theta
  
  n.obs = sum(var.mat)
  if(n.obs==0){
    stop("pval_calculation_per_gene: no mutation observed but sent to pval calculation, double check the data\n");
  }
    mu.nb <- as.vector(t(mu.nb.mat[,index.nb]))    
    theta.nb <- rep(x=theta.v[index.nb], nrow(mu.nb.mat[,index.nb]))
    pos.nb <- as.vector(t(kmer.mat[,index.nb]))
    alfa.nb <- theta.nb
    p.v <- 1.0/(1.0+mu.nb/theta.nb)
    alfa.v <- alfa.nb*pos.nb

    p.obs <- pcal.ng.many(n.obs, alfa.v, p.v) 
  if(p.obs<0){
    p.obs = 0
  }
  print(p.obs)
  return(p.obs)
}

pval_calculation_all_genes <- function(dtest.prepare,flag.merge.input,file.v){
  id.unique <- dtest.prepare$d.test.kmer.merge$id.unique
  tmp <- strsplit(as.character(id.unique), split = ":")
  geneName <- sapply(tmp,function(x) x[1])
  geneName.uniq <- unique(geneName)

  #/*====== calculate merged annotation length ======*/
  len.mat <- dtest.prepare$d.test.kmer.merge[,-1]
  len.tmp <- aggregate(. ~ geneName, data = len.mat, sum)
  len.val <- apply(len.tmp[,-1], 1, sum)
  len.val <- len.val[match(geneName.uniq, len.tmp$geneName)]

  #/*====== calculate number of samples with mutation ======*/
  #sampleCt <- function(x){return(sum(x>0))}
  #sample.mat <- dtest.prepare$d.test.sample.merge[,-1]
  #sample.tmp <- aggregate(. ~ geneName, data = sample.mat, sum)
  #sample.val <- apply(sample.tmp[,-1], 1, sampleCt)
  #sample.val <- sample.val[match(geneName.uniq, sample.tmp$geneName)]
  
  #/*====== note that it is possible for a test to have only ppois or nb background ======*/
  #index.pois <- which(dtest.prepare$d.train.nb.theta$pval_rou>0.05)
  #index.nb <- which(dtest.prepare$d.train.nb.theta$pval_rou<=0.05)
  index.nb <- c(1:64)
  index.pois <- 0
  
  #/*====== variant and background count ======*/
  p.val <- rep(1,length(geneName.uniq))
  #/*====== select elments with mutation 0 variants and their pval  = 1, no need to calculate by default ======*/
  var.mat <- data.frame(geneName,dtest.prepare$d.test.var.merge[,-1])
  var.mat.merge <- aggregate(. ~ geneName, data = var.mat, sum )
  var.sum <- apply(as.matrix(var.mat.merge[,-1]),1,sum)
  var.sum <- var.sum[match(geneName.uniq, var.mat.merge$geneName)]
  index.pval.calculate <- which(var.sum>0 & len.val>0)
  print(index.pval.calculate)
  print(length(index.pval.calculate))

  #/*====== calculate pval for genes with elements ======*/
  for(i in 1:length(index.pval.calculate)){
    if(i%%100==1){
      print(paste("------ pval ",i,"out of",length(index.pval.calculate),"calculations"))
    }
    i.org <- index.pval.calculate[i]
    gene.org <- geneName.uniq[i.org]
    p.obs = pval_calculation_per_gene(gene.org, dtest.prepare, geneName, index.pois, index.nb)
    print("pause")
    p.val[i.org] = p.obs
  }
  data.pval <- data.frame(geneName.uniq,len.val,var.sum,var.sum/len.val,#sample.val, rep(ncol(sample.tmp),length(p.val))
                          p.val)
  colnames(data.pval) <- c("name","len","ct","den","p")    
  data.pval$padj <- p.adjust(data.pval$p)
  data.pval <- data.pval[order(data.pval$p),]
  return(data.pval)
}