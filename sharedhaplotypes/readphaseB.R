## script to read PHASE pairs output

readpairs <- function(filename)
{
  conn = file(filename,"r")
  res <- readLines(conn)
  inds = grep("IND",res)
  indnames <- substring(res[inds],6)
  n <- diff(c(inds,length(res)+1))-1
  res <- gsub("\\s","",res[-inds])
  namecol <- rep(indnames,n)
  data <- matrix(unlist(strsplit(res,",")),ncol=3,byrow=T)
  r = data.frame(cbind(namecol,data))
  r[,4] <- as.numeric(paste(r[,4]))
  close(conn)
  r
}


splitter <- function(pairdata,wh=1) {
  sp <- function(x) {
    c(substr(x,15,15),substr(x,1,7),substr(x,8,8),substr(x,8,11),substr(x,12,13)
      ,paste(substr(x,14,14),substr(x,16,19),sep="")
             ,substr(x,20,22),substr(x,23,24))
  }
  t(sapply(paste(pairdata[,1+wh]),sp))
}
## Now using a vim editing mode


splitter2 <- function(pairdata,wh=1) {
  sp <- function(x) {
    c(substr(x,15,15),substr(x,1,9),paste(c(substr(x,10,14),substr(x,16,18)),sep=""),substr(x,19,24))
  }
  t(sapply(paste(pairdata[,1+wh]),sp))
}



splitter3 <- function(pairdata,wh=1,left=11,right=18) {
  sp3 <- function(x) {
    c(HMERF=substr(x,15,15),left=substr(x,1,left-1),mid=paste(c(substr(x,left,14),substr(x,16,right)),collapse=""),right=substr(x,right+1,24))
  }
  t(sapply(paste(pairdata[,1+wh]),sp3))
}


splitter4 <- function(pairdata,wh=1,ll=8,left=11,right=18) {
  sp4 <- function(x) {
    c(HMERF=substr(x,15,15),ll=substr(x,1,ll-1),left=substr(x,ll,left-1),mid=paste(c(substr(x,left,14),substr(x,16,right)),collapse=""),right=substr(x,right+1,24))
  }
  t(sapply(paste(pairdata[,1+wh]),sp4))
}


filenames <- c("phase1K.out_pairs","phase1K.out_pairs")
a <- lapply(filenames,readpairs)
w <- lapply(a,function(x) c(min(grep("HG",x[,1])),nrow(x)))
#b <- lapply(a,function(x) cbind(x,splitter(x,1),splitter(x,2)))
b <- lapply(a,function(x) cbind(x,splitter4(x,1,9,12,18),splitter4(x,2,9,12,18)))
b[[1]][1:3,]
bcase <- mapply(function(x,y) x[1:(y[1]-1),],b,w,SIMPLIFY=FALSE)
casefreq(bcase[[1]],3,nsplit=4)
casefreq(bcase[[1]],2,nsplit=4)

bcontrol <- mapply(function(x,y){return(x[y[1]:y[2],])},b,w,SIMPLIFY=FALSE)
tb1 = table(c(paste(bcontrol[[1]][,7]),paste(bcontrol[[1]][,11])))
table(c(paste(bcontrol[[1]][,7]),paste(bcontrol[[1]][,12])))

b[[1]][1:10,]
table(c(paste(bcontrol[[1]][,6]),paste(bcontrol[[1]][,14])))

casefreq <- function(bbb,cl=1,nsplit=3) {
  fam <- paste(bbb[,1])
  hap <- paste(bbb[,5+cl])
  hap[bbb[,5]=="A"] <- paste(bbb[bbb[,5]=="A",6+nsplit+cl])
  tapply(bbb[,4],data.frame(fam=fam,hap=hap),sum)
}

casefreq2 <- function(bbb,cl=1) {
  fam <- paste(bbb[,1])
  hap1 <- paste(bbb[,5+cl])
  hap2 <- paste(bbb[,13+cl])
  
  tapply(bbb[,4],data.frame(hap1=hap1,hap2=hap2,fam=fam),sum)
}




casefreq(bcase[[1]],1)
casefreq(bcase[[1]],3)
casefreq(bcase[[1]],4)
casefreq(bcase[[1]],5)

casefreq(bcase[[1]],6)
casefreq(bcase[[1]],7)


write.csv(bcase[[1]],file="bcase1.csv")
write.csv(bcase[[2]],file="bcase2.csv")
write.csv(bcontrol[[1]],file="bcontrol1.csv")
write.csv(bcontrol[[2]],file="bcontrol2.csv")



tapply(bcase[[1]][,4],list(paste(bcase[[1]][,1]),bcase[[1]][,6]),sum)



r <- mapply(function(x,y) allfamilymatches(x[1:(y-1),]),a,w)
allfamilymatches()



a = readpairs("../PHASE/phaseall1.outshort_pairs")
haps <- read.table("../PHASE/phase.out_freqs",header=T)
removemut <- matrix(unlist(strsplit(paste(haps[,2]),"")),nrow=nrow(haps),byrow=TRUE)
diseasemut <- removemut[,15]
newhap <- apply(removemut[,-15],1,paste,collapse="")


allfamilymatches <- function(xx)  {
         fams <- names(table(paste(xx[,1])))
         nf = length(fams)
         for (ii in 1:(nf-1)) {
           fam1=fams[ii]
           cat("Trying disease haplotype from ",fam1,"\n")
               for (jj in (ii+1):nf) {
                 fam2=fams[jj]
                 cat(fam2)
                     
                     famhaps <- c(paste(xx[xx[,1]==fam1,2]),paste(xx[xx[,1]==fam1,3]))
                      fp = c(paste(xx[xx[,1]==fam1,4]),paste(xx[xx[,1]==fam1,4]))
                #     famhaps <- famhaps[substr(famhaps,15,15)=="A"]
                #      fp = fp[substr(famhaps,15,15)=="A"]
                       notfamhaps <- c(paste(xx[xx[,1]==fam2,2]),paste(xx[xx[,1]==fam2,3]))
                      nfp <-  c(paste(xx[xx[,1]==fam2,4]),paste(xx[xx[,1]==fam2,4]))
              #       notfamhaps <- notfamhaps[substr(notfamhaps,15,15)=="A"]
               #      nfp <- nfp[substr(notfamhaps,15,15)=="A"]
                     m <- match(famhaps,notfamhaps)
                     if (sum(!is.na(m))==0) {
                       cat("  No match","\n")  
                 } else {
                   cat("  Matches\n")
                   u=which(!is.na(m))
                   mm <- m[u]
                   res <- cbind(rep(fam2,length(u)),famhaps[u],notfamhaps[mm],fp[u],nfp[mm])
                   write.table(res,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" & ")
                  cat("\n")
                 }
           }
     }
}
allfamilymatches(a[[1]][1:300,])


## all data
raw <- as.matrix(read.table("../data/1000match"))
raw <- t(raw)
raw <- raw[keep,]

b = apply(raw,1,paste,collapse="")

a = apply(raw,1,paste,collapse="")

raw <- as.matrix(read.table("../data/1000match"))
raw <- t(raw)
raw <- raw[keep,]

b = apply(raw,1,paste,collapse="")
u <- substr(paste(b),15,23)=="ATGTGTCAC"
u2 <-  paste(b)=="ACACATACGGATACATGTGTCAC"    # X34, X40 and Calgary
u3 <- paste(b) == "ACACATAGAGGTATATGTGTCAC"
u4 <- paste(b) == "ACACATACGGGTATATGTGTCAC"
u5 <- paste(b) == "ACACATACGGGTACACGGGTCGT"
sum(u5)/length(u5)
u6 <- paste(b) == "ACACATACGGGTACGCGGGTCGT"
u7 <- paste(a) == "ACACATACAAGCATGTGTGTCAC"
sum(u7)/length(u5)
vals <- c("GTACCCGCGGGTACATGTGTCAC","GTACCCGCGGGTATGTGTGTCAC","GTACCCGCGGGTATGTGTGTCAC","GTACCCGCGGGTACATGTGTCAC"
          ,"GTACCCGCGGGTATATGTGTCAC","GTACCCGCGGGTATATGTGTCAC")
for (v in vals) {
  print(sum(paste(b)==v)/length(u5))
}
u8 <- paste(b) == "GTACCCGCGGGTATGATGTGTCAC"
sum(u8)/length(u5)

