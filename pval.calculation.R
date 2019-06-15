rm(list=ls())
setwd("~")

# e.g.,
# cellType = "HepG2"

ca_type = "Lymph-BNHL"
path_data="~/"
source(paste(path_data,"convolution.methods.nb.pois.R",sep=""))
#/*============ load data ============*/
library(data.table)
name = paste(path_data,ca_type,"input.txt",sep="")
dct <- fread(name,header=F,sep="\t")
dct <- data.frame(dct)
colnames(dct) <- c("name","binCt","len","var_ct","bin_str","len_str","var_ct_sum_str","kmer_str","var_ct_str")

name = paste(path_data,ca_type,".poisson.mu.txt",sep="")
dmu.pois <- fread(name,header=F,sep="\t")
dmu.pois <- data.frame(dmu.pois)
colnames(dmu.pois) <- c('bin','AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT')


name = paste(path_data,ca_type,".nb.mu.txt",sep="")
dmu.nb <- fread(name,header=F,sep="\t")
dmu.nb <- data.frame(dmu.nb)
colnames(dmu.nb) <- c('bin','AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT')


name = paste(path_data,ca_type,".nb.theta.txt",sep="")
dtheta.nb <- fread(name,header=F,sep="\t")
dtheta.nb <- data.frame(dtheta.nb)
colnames(dtheta.nb) <- c('bin','AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT')


name = paste(path_data,ca_type,".disp.pval.txt",sep="")
d.disp.para <- fread(name,header=T,sep="\t")
d.disp.para <- data.frame(d.disp.para)
print(d.disp.para)
d.disp.para$p.adj <- p.adjust(d.disp.para$pval)
print(d.disp.para$p.adj)
flag.nb <- which(d.disp.para$p.adj<0.05)
flag.pois <- which(d.disp.para$p.adj>=0.05)
print('done')
table(dct$var_ct)
#/*============ calculate pvalue ============*/
n.annotation <- dim(dct)[1]
pval.anno <- rep(1,n.annotation)
index.calculate <- which(dct$var_ct>0)
print('done2')
for(i in 1:length(index.calculate)){
  if(i%%500==0){
    print(paste(i,"lines out of",length(index.calculate),"completed"));
  }
  k = index.calculate[i]
  ct.val = dct$var_ct[k]
  if(ct.val>0){
    bin.val = as.character(unlist(strsplit(dct$bin_str[k],split=",")))
    if(length(bin.val)<=6){
      flag.nb.total = rep(0,0)
      if(length(bin.val)==1){ #/*most cases there is just one bin*/
        kmer.v <- as.numeric(unlist(strsplit(dct$kmer_str[k],split=",")))
        mu.pois.v <- dmu.pois[dmu.pois$bin==bin.val,flag.pois+1]
        mu.nb.v  <- dmu.nb[dmu.nb$bin==bin.val,flag.nb+1]
	mu.nb.v = as.numeric(unlist(mu.nb.v))
        theta.v <- dtheta.nb[dtheta.nb$bin==bin.val,flag.nb+1]
	theta.v = as.numeric(unlist(theta.v))
	mu.pois.v = as.numeric(unlist(mu.pois.v))
        lambda.pois = sum(mu.pois.v*kmer.v[flag.pois])
        flag.nb.total = flag.nb
      } else { #/*========== genes with >1 bin =========*/
        nbin = length(bin.val)
        kmer.v.str <- unlist(strsplit(dct$kmer_str[k],split = ":"))
        kmer.v <- rep(0,0)
        mu.pois.v <- rep(0,0)
        mu.nb.v <- rep(0,0)
        theta.v <- rep(0,0)
        lambda.pois = 0;
        for(jj in 1:length(bin.val)){
          kmer.v <- c(kmer.v, as.numeric(unlist(strsplit(kmer.v.str[jj],split = ","))))
          mu.pois.v <- c(mu.pois.v, as.numeric(dmu.pois[dmu.pois$bin==bin.val[jj],flag.pois+1]))
          mu.nb.v <- c(mu.nb.v, as.numeric(dmu.nb[dmu.nb$bin==bin.val[jj],flag.nb+1]))
          theta.v <- c(theta.v, as.numeric(dtheta.nb[dtheta.nb$bin==bin.val[jj],flag.nb+1]))
          lambda.pois = lambda.pois + sum(mu.pois.v*kmer.v[flag.pois])
          flag.nb.total=c(flag.nb.total, flag.nb+64*(jj-1))
        }
      }
      sigma.ng.used <- 1.0/theta.v
      ng.re <- mu2p(as.numeric(mu.nb.v), as.numeric(sigma.ng.used))

      alfa.used <- as.numeric(ng.re$alfa.val*kmer.v[flag.nb.total])

      p.val.used <- as.numeric(ng.re$p.val)

      p.value = 1; #/*default value*/
      
      if(sum(mu.nb.v*kmer.v[flag.nb])==0){#/*pois only case*/
        if(lambda.pois>0){
          p.value = ppois(ct.val-1, lambda=lambda.pois, lower.tail=F)
        }
      }else if(lambda.pois==0){ #/*this is the only NG model*/
        if(sum(alfa.used)>0){
          p.value = pcal.ng.many(ct.val,alfa.used, p.val.used)
        }
      }else{ #/*this is the many nb and pois case*/
        p.value = pcal.ng.many.plus.pois(ct.val, alfa.used, p.val.used, lambda.pois)
      }
      pval.anno[k] = p.value
      # stop("find a non-zero one")
    }
  }
}

sum(p.adjust(pval.anno)<0.05)
dtmp <- data.frame(dct$name, dct$var_ct, dct$len, pval.anno,p.adjust(pval.anno))
colnames(dtmp) <- c("name","ct","len","pval","p.adj")
dtmp <- dtmp[order(dtmp$pval),]
dtmp$den = dtmp$ct/dtmp$len
dtmp[p.adjust(dtmp$pval)<0.1,]
name = paste(path_data,ca_type,".pval.estimation.txt",sep="")
write.table(dtmp, file = name, append = FALSE, quote = F, sep = "\t", row.names = F,col.names = TRUE)
