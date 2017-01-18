
sperm_count <-
  c(80.295, 77.191, 169.566, 19.369, 136.779, 12.938, 12.514, 44.975, 
    117.89, 22.329, 22.113, 146.105, 79.137, 95.176, 36.811, 61.604, 
    89.553, 88.86, 35.404, 99.811, 43.441, 77.926, 62.545, 41.976, 
    102.642, 151.866, 20.169, 25.997, 17.705, 24.782, 37.28, 89.019, 
    55.345, 241.77, 198.505, 21.919, 37.32, 60.079, 82.141, 93.074, 
    104.943, 162.588, 118.786)
oocyte_count <-
  c(1249200L, 1393650L, 1031100L)

## Hypothesis Test
##
## This function tests the hypothesis that the frequency of the 
## paternal haplotype in the child is at the same relative frequency
## as the relative numbers of mtDNA copies transmitted in the sperm
## and oocyte, against an alternative hypothesis that the frequency 
## is lower.
##
##  bootHT returns the p-value.

bootHT <- function(mother,child,nm,nc,n_boot=100000) {
  mo <- mean(oocyte_count)
  sdo <- sd(oocyte_count)
  s.cn <- replicate(n_boot,mean(sample(sperm_count,43,replace=TRUE)))
  o.cn <- replicate(n_boot,mean(rnorm(3,mo,sdo)))
  ratio <- s.cn/o.cn
  noise.m <- rpois(n_boot,nm*mother/nm)
  noise.c <- rpois(n_boot,nc*mother/nm)
  child.fromfather <- rbinom(n_boot,nc,ratio)
  d <- (child.fromfather+noise.c)/nc - noise.m/nm
  obs <- child/nc - mother/nm
  return(rank(c(obs,d))[1]/n_boot)
}


PaternalTransmission0 <- function(mother,child,nm,nc,n_boot=100000) {
  noise.m <- rpois(n_boot,nm*mother/nm)
  noise.c <- rpois(n_boot,nc*mother/nm)
  d <- noise.c/nc - noise.m/nm
  obs <- child/nc - mother/nm
  return(1- rank(c(obs,d))[1]/n_boot)  ## upper tail needed
}


## tests done in manuscript.  

if (FALSE) {  ## do not run
  ## Extract data for A1
  trioHaplogroups <- read.csv("TrioHaplogroups.csv")
  A1 <- trioHaplogroups[trioHaplogroups$trio==1 & trioHaplogroups$motif=="A",]
  father <- A1$haplotype[A1$individual=="father"]
  bootHT(11,6,
         nm=  sum(A1$count[A1$individual=="mother"]),
         nc=sum(A1$count[A1$individual=="child"]),1E6)
  PaternalTransmission0(11,6,
         nm=  sum(A1$count[A1$individual=="mother"]),
         nc=sum(A1$count[A1$individual=="child"]),1E6)
  #B2
  B2 <- trioHaplogroups[trioHaplogroups$trio==2 & trioHaplogroups$motif=="B",]
  bootHT(25,16,
         nm=  sum(B2$count[B2$individual=="mother"]),
         nc=sum(B2$count[B2$individual=="child"]),1E6)
  PaternalTransmission0(25,16,
         nm=  sum(B2$count[B2$individual=="mother"]),
         nc=sum(B2$count[B2$individual=="child"]),1E6)
  #C3
  C3 <- trioHaplogroups[trioHaplogroups$trio==3 & trioHaplogroups$motif=="C",]
  bootHT(5,5,
         nm=  sum(C3$count[C3$individual=="mother"]),
         nc=sum(C3$count[C3$individual=="child"]),1E6)
  PaternalTransmission0(5,5,
         nm=  sum(C3$count[C3$individual=="mother"]),
         nc=sum(C3$count[C3$individual=="child"]),1E6)
  
  #D4
  D4 <- trioHaplogroups[trioHaplogroups$trio==4 & trioHaplogroups$motif=="D",]
  bootHT(2,4,
         nm=  sum(D4$count[D4$individual=="mother"]),
         nc=sum(D4$count[D4$individual=="child"]),1E6)
  PaternalTransmission0(2,4,
         nm=  sum(D4$count[D4$individual=="mother"]),
         nc=sum(D4$count[D4$individual=="child"]),1E6)
}



powerFunc <- function(coverage, het, freq, reps = 1e+05, p = 0.05, show_plot = FALSE) {
  if (require(VGAM) == FALSE) {
    stop("This function requires the R package VGAM.\n"
         , "Install it using \n>install.packages(\"VGAM\")")
  }
  ## need to find the probability that the different between a poisson
  ## with a mean coverage*(het+freq) and a poisson with mean coverage*het
  ## is in the bottom 5% of the difference between two poissons.
  
  ## do this by simulation We simulate two sets of Poisson differences.
  ## For a one-tailed test the critical value is the lower tail so we get
  ## the lower tail of H_0 and see what proportion of
  
  d_h0 <- rskellam(reps, coverage * (het + freq), coverage * het)  ## null
  d_h1 <- rskellam(reps, coverage * het, coverage * het)  ## alternative  
  critical_value <- quantile(d_h0, probs = c(p))  ## the critical value
  if (show_plot) {
    plot(density(d_h0), xlim = range(density(c(d_h0, d_h1))$x))
    lines(density(d_h1), col = "blue")
    lines(c(critical_value, critical_value), c(0, max(density(d_h0)$y/3)), 
          col = "red")
  }
  sum(d_h1 <= critical_value)/reps
}

##---------------------------------------------------------------------------------------------


powerFuncB <- function(coverage, het, freq, reps = 1e+05, p = 0.05, show_plot = FALSE) {
  if (require(VGAM) == FALSE) {
    stop("This function requires the R package VGAM.\n"
         , "Install it using \n>install.packages(\"VGAM\")")
  }
  ## need to find the probability that the different between a poisson
  ## with a mean coverage*(het+freq) and a poisson with mean coverage*het
  ## is in the bottom 5% of the difference between two poissons.
  
  ## 
  d_h0 <- rskellam(reps, coverage * het, coverage * het)  ## null distribution
  d_h1 <- rskellam(reps, coverage * (het + freq), coverage * het)  ## H_1
  critical_value <- quantile(d_h0, probs = c(1 - p))  ## the critical value
  if (show_plot) {
    plot(density(d_h0), xlim = range(density(c(d_h0, d_h1))$x))
    lines(density(d_h1), col = "blue")
    lines(c(critical_value, critical_value), c(0, max(density(d_h0)$y/3)), 
          col = "red")
  }
  sum(d_h1 > critical_value)/reps
}

if (FALSE) {
  # example code
  myCoverage <- c(10000, 20000, 50000, 1e+05, 2e+05, 5e+05, 1e+06)
  pow <- sapply(myCoverage, powerFunc, het = 0.005, freq = 0.005)
  # note that this is equivalent but longer winded
  pow <- numeric(length(myCoverage))
  for (i in 1:length(myCoverage)) {
    pow[i] <- powerFunc(myCoverage[i], het = 0.001, freq = 1e-04)
  }
  
  plot(myCoverage, pow, log = "x", ylim = c(0, 1))
  powB <- sapply(myCoverage, powerFuncB, het = 0.005, freq = 0.005)
  points(myCoverage, powB, col = "red")
}
