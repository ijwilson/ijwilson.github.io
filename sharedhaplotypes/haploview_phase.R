#!/usr/bin/Rscript

## output in hapmap PHASE output.   We want to be able to use the 1000 genomes
## data to get a nice picture of the phasing information to furhter inform these data

## this script outputs data in the haps format with an associated legend file.

ped <- read.csv("phase1_samples_integrated_20101123.ped",sep="\t",header=TRUE,quote="")
NoParents <- ped[,3]=="0" & ped[,4]=="0"
RequiredPopulations = c("CEU","GBR")
ped = ped[NoParents &  ped$Population %in% RequiredPopulations,]
SampleNames <- scan("impute1K.impute.hap.indv",what=character(0))
keepSamples = SampleNames %in% ped[,2]
doubleKeepSamples <- rep(keepSamples,rep(2,length(keepSamples)))
UsedSampleNames <- SampleNames[keepSamples]
USN2 <- rep(UsedSampleNames,rep(2,length(UsedSampleNames)))
families <- paste("fam",USN2,sep="_")
## the data file
raw2 <- t(as.matrix(read.table("match1K")))[doubleKeepSamples,]
raw2[raw2=="A"] <- 1
raw2[raw2=="C"] <- 2
raw2[raw2=="G"] <- 3
raw2[raw2=="T"] <- 4

raw2b <- cbind(families,USN2,raw2)
write.table(raw2b,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,file="match1K.haps")

## the legend file
m1kdata <- read.table("match1K.legend",header=TRUE)
write.table(m1kdata[,1:2],col.names=FALSE,row.names=FALSE,quote=FALSE,file="match1K.info",sep="\t")
