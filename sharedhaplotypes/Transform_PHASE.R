#!/usr/bin/Rscript
##################################################################
Splitmarkers <- function(str,SNP)
  {
    str <- paste(str)
     if (str=="") {
       if (SNP) return (c("?","?"))
       else return(c("-1","-1"))
     }
    if (grepl('/',str)==TRUE) return(unlist(strsplit(str,"/")))
    else return(unlist(strsplit(str,"")))
  }
##################################################################
## test for SNPs by determining whether the first character of the data
## is a base
locusSNP <- function(col) {
  issnp <- function(yy) {
    if (substr(yy,1,1)  %in% c("A","C","T","G")) return (TRUE)
    return (FALSE)
    }
  return(sapply(col,issnp)) 
}
####################################################################

ped <- read.csv("phase1_samples_integrated_20101123.ped",sep="\t",header=TRUE,quote="")
use <- ped[,3]=="0" & ped[,4]=="0"
eth = c("CEU","GBR")
useB = ped$Population %in% eth
ped = ped[use & useB,]
allSampleNames <- scan("impute1K.impute.hap.indv",what=character(0))
keep = allSampleNames %in% paste(ped[,2])
UsedSampleNames <- allSampleNames[keep]


hmerf <- read.csv("HMERFmatch1K.csv")
f <- file("phase1K.inp","w")
known = file("phase1K.known","w")
n <- ncol(hmerf) -2
markers=nrow(hmerf)

SNPs = locusSNP(hmerf[,3])       ## the most complete data
cat(n+length(UsedSampleNames),markers,sep="\n",file=f)
cat("P ",hmerf[,1],"\n",file=f)
loc = rep("M",markers)
loc[SNPs] <- "S"
cat(loc,"\n",file=f)
for (ind in 1:n) {
  col=ind+2
  cat(colnames(hmerf)[col],"\n",file=f)
  vals <- mapply(Splitmarkers,hmerf[,col],SNPs)
  cat(vals[1,],"\n",file=f)
  cat(vals[2,],"\n",file=f)
  cat(rep("*",ncol(vals)),"\n",file=known)
}

raw <- t(as.matrix(read.table("match1K")))
## now add the mutation - first add the column
raw <- raw[,c(1:14,1,15:ncol(raw))]
raw[,15] <- "A"
PhasedLine = paste(rep("0",ncol(raw)),collapse=" ")


raw <- apply(raw,1,paste,collapse=" ")
raw <- matrix(raw,ncol=2,byrow=T)
raw <- raw[keep,]

for (i in 1:length(UsedSampleNames)) {
  cat(UsedSampleNames[i],"\n",file=f)
  cat(raw[i,],sep="\n",file=f)
  cat(PhasedLine,"\n",file=known)
}

close(f)
close(known)

