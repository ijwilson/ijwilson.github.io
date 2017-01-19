## function to read PHASE pairs output
readpairs <- function(filename) {
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
  return(r)
}
# function to split pairs data by a single positions
splitter <- function(pairdata,splits=c(1)) {
  nsplit <- length(splits)+1
  l <- nchar(paste(pairdata[1,3]))
  splitvecstart <- c(1,splits)
  splitvecend <- c(splits,l)
  pd <- rep(paste(pairdata[,2]),rep(nsplit,nrow(pairdata)))
  s1 <- substring(pd,splitvecstart,splitvecend)
  pd1 <- rep(paste(pairdata[,3]),rep(nsplit,nrow(pairdata)))
  s1 <- substring(pd,splitvecstart,splitvecend)
  data.frame(pairdata,hap1 = matrix(s1,ncol=nsplit,byrow=T),hap2 = matrix(s1,ncol=nsplit,byrow=T))
}

pairs <- readpairs("phase1K.out_pairs")
splitter(pairs,c(4,8,11))

w = min(grep("HG",pairs[,1]))

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

raw <- as.matrix(read.table("match1K"))
raw <- t(raw)

b = apply(raw,1,paste,collapse="")
u <- substr(paste(b),15,23)=="ATGTGTCAC"

print(sum(u))
