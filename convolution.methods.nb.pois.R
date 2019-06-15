mu2p <- function(mu.val,sigma.val){
  alfa.val <- 1.0/sigma.val
  p.val <- 1.0/(1.0+sigma.val*mu.val)
  return(list(alfa.val=alfa.val, p.val=p.val))
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
    x1 <- seq(0,n.obs)
    p.calculate <- dnbinom(x1, size=alfa.v[1], prob=p.v[1], log = FALSE)
    for(i in 2:length(alfa.v)){
      p.calculate.new <- lapply(x1, pcal.ng.two, alfa=alfa.v[i],p=p.v[i],p.org = p.calculate)
      p.calculate <- unlist(p.calculate.new)
      # print(paste(i,p.calculate))
      if(sum(is.na(p.calculate))>0){
        stop("pcal.ng.many: NAN generated,p.calculate=",p.calculate )
      }
    }
    p.value <- (1 - sum(p.calculate))
    if(p.value<0){ #/*for numerical calculation*/
      p.value = 0
    }
  }
  return (p.value)
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
      # print(paste(i,p.calculate))
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
