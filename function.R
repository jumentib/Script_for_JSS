#' simulator_ewas : function to simulate DNA methylation data for EWAS
#'
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param freq : (vector) mean methylation values (if NULL, set randomly)
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.variance : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param sigma :	standard deviation of residual errors
#' @param sd.A :	standard deviation for effect sizes
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#'
#' @return Y : matrix of methylation beta values
#' @return Ynorm : pnorm(Y)
#' @return X : exposure
#' @return A : effect sizes exposure
#' @return causal : set of CpGs associated with the exposure
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#'
#' @details
#'
#' This function is used to simulate datasets for EWAS.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X and U and prop.variance
#' (intensity of the confounders or correlation between X and U).
#' Then this matrix is used to simulate via normal laws X and U.
#' Thereafter, the effect sizes of X (A) and U (V) are calculated
#' using mean parameters of effect sizes (meanA) and standard deviations (sdA and sdV).
#' Note that the effect sizes of X are calculated only for causal mediators with X.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' To finish the methylation matrix is calculated thanks to the formula : Y = V*U + A*X + Z
#' @examples
#' # Simulate data :
#' simu <- simulator_ewas(100, 500, 5)
#' @export
simulator_ewas <- function (n = 100,
                            p = 500,
                            K = 5,
                            freq = NULL,
                            prop.causal = 0.010,
                            prop.variance = 0.4,
                            sigma = 1,
                            sd.A = 0.2,
                            mean.A = 1.0,
                            sd.U = 1.0,
                            sd.V = 1.0)
{
  causal <- sample.int(p, prop.causal * p)
  
  x.nb = length(causal)
  
  if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site
  
  cs <- runif(K, min = -1, max = 1)
  
  theta <- sqrt( prop.variance /sum((cs/sd.U)^2) )
  
  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)
  
  Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))
  
  Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))
  
  UX <- MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma)
  U <- UX[, 1:K, drop = FALSE]   # confounders
  X <- UX[, K + 1, drop = FALSE] # outcome
  
  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))
  
  A <- matrix(0, p, 1)
  A[causal, 1] <- rnorm(x.nb, mean.A, sd.A)
  
  Epsilon <- apply(matrix(rep(0,p),nrow = 1), 2, function(x) rnorm(n,x,sigma))
  
  Z = U %*% t(V) + X %*% t(A) + Epsilon
  
  M = matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z
  
  N = M
  M = pnorm(M)
  
  return(list(Y = N,
              X = X,
              A = A,
              U = U,
              V = V,
              Z = Z,
              K = K,
              Ynorm = M,
              causal = causal))
}

# function for Root Mean Square Error calculation
RMSE <- function(beta, true.beta) {
  return(sqrt(mean((true.beta - beta)^2)))
}

# function use for Generative simulation for F1 score, precision and recall 
#' F1 : Calcul F1 score for a simulation (on effect sizes)
#'
#' @param beta : effect from a regression
#' @param causal : true causal variable (SNP or CpG)
#' @param nb.hit : number of first hit for calculate the F1 score
#'
#' @return the F1 score, the precision and the recall
#' 
#' @export
F1 <- function(beta, causal, nb.hit = 10) {
  abs.beta <- abs(beta)
  p <- length(beta)
  # create list for F1
  first.hit <- order(abs.beta, decreasing = T)[1:nb.hit]
  last.hit <- order(abs.beta, decreasing = T)[(nb.hit + 1):p]
  # F1 score
  tp <- sum(first.hit %in% causal)
  fp <- sum(!(first.hit %in% causal))
  fn <- sum(causal %in% last.hit)
  
  preci <- tp / (tp + fp)
  recal <- tp / (tp + fn)
  
  f1 <- 2 * (preci * recal) / (preci + recal)
  
  preci <- ifelse(is.na(preci), 0, preci)
  recal <- ifelse(is.na(recal), 0, recal)
  f1 <- ifelse(is.na(f1), 0, f1)
  
  return(list(f1 = f1, precision = preci, recall = recal))
}

# F1 score, precision and recall for empirical simulation
#' F1.LD2 : Calcul F1 score for a simulation (on effect sizes), take in account the linkage desequilibrum
#'
#' @param beta : effect from a regression
#' @param causal : true causal variable (SNP or CpG)
#' @param nb.hit : number of first hit for calculate the F1 score
#' @param poschr : Genetic map: positions on chromosome 
#' @param dist.max : see details section
#' @param cor.min : see details section
#'
#' @return the F1 score, the precision and the recall
#' 
#' @details 
#' 
#' Candidate SNPs that are not within a 10 Kb (dist.max) window centered on a causal SNP 
#' (and with a correlation of less than 0.2 (cor.min) with this same SNP) 
#' are considered false discoveries.
#'
#' @export
F1.LD2 <- function(beta, geno, poschr, causal, nb.hit = 50, dist.max = 10000, cor.min = 0.2) {
  
  neighbour.p = 50
  pos <- sort(order(abs(beta), decreasing = T)[1:nb.hit])
  pN <- NULL
  for (j in 1:length(pos)) {
    sel <- pos[j]
    p.sel <- poschr[sel]
    
    dis <- abs(p.sel - poschr)
    
    colnames(geno) <- 1:ncol(geno)
    
    g <- geno[, (dis <= dist.max)]
    
    if (is.null(ncol(g))) {
      pN <- c(pN, sel)
    }
    
    else {
      c <- rep(NA, ncol(g))
      for (i in 1:ncol(g)) {
        c[i] <- cor(geno[, sel], g[, i])^2
      }
      
      g <- g[, c >= cor.min]
      if (is.vector(g)) {
        pN <- c(pN, sel)
      }
      else {
        pN <- c(pN, as.numeric(colnames(g)))
      }
    }
  }
  
  
  NbP <- NULL
  for (j in 1:length(pos)) {
    pos1 <- pos[j]
    pos2 <- pos[j + 1]
    
    pos50 <- pos1:(pos1 + neighbour.p)
    
    if ((pos1 %in% pos50) & (pos2 %in% pos50)) {
      NbP <- c(NbP, 0)
    } 
    else {NbP <- c(NbP, 1)}
  }
  # 
  np <- sum(NbP) - 1
  
  # 
  tp <- sum(causal %in% pN)
  # 
  fp <- np - tp
  fp <- ifelse(fp < 0, 0, fp)
  # 
  fn <- sum(!(causal %in% pN))  
  
  preci <- tp / (tp + fp)
  recal <- tp / (tp + fn)
  
  f1 <- 2 * (preci * recal) / (preci + recal)
  
  preci <- ifelse(is.na(preci), 0, preci)
  recal <- ifelse(is.na(recal), 0, recal)
  f1 <- ifelse(is.na(f1), 0, f1)
  
  return(list(f1 = f1, precision = preci, recall = recal, tp = tp, fp = fp, fn = fn, np = np))
}

##### GEMMA ######

# Write the phenotype data to a file in the format used by GEMMA. Each
# line of the file contains one phenotype observation.

#' write.gemma.pheno
#'
#' @param pheno phenotype
#' @return nothing
#'
#' @export
write.gemma.pheno <- function (file, pheno) {
  phenotype <- 1:nrow(pheno)
  y <- pheno[phenotype] # j'ai retiré une "[]"
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y,file,quote = FALSE,row.names = FALSE,col.names = FALSE)
}


# Write the covariate data to a file in the format used by GEMMA. Each
# line corresponds to a sample. We must include an additional
# covariate for the intercept.

#' write.gemma.covariates
#'
#' @param covariates covariable
#' @return nothing
#'
#' @export
write.gemma.covariates <- function (file, covariates, pheno) {
  if (length(covariates) == 0) {
    write.table(data.frame(rep(1,nrow(pheno))),file,sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)
  } else {
    round.col <- function (x) {
      if (is.numeric(x))
        round(x,digits = 6)
      else
        x
    }
    write.table(cbind(1,data.frame(lapply(pheno[covariates],round.col))),
                file,sep = " ",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
  }
}

# Write the SNP information to a space-delimited text file in the
# format used by GEMMA. This file contains one line per SNP, with
# three columns: (1) SNP label, (2) base-pair position, (3)
# chromosome.

#' write.gemma.map
#'
#' @param map map
#' @return nothing
#'
#' @export
write.gemma.map <- function (file, map)
  write.table(map,file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

# Store the mean genotypes as a space-delimited text file in the
# format used by GEMMA, in which we have one row per SNP, and one
# column per sample. The first three column give the SNP label, and
# the two alleles.

#' write.gemma.geno
#'
#' @param geno genotype
#' @return nothing
#'
#' @export
write.gemma.geno <- function (file, geno, map) {
  geno <- t(geno)
  geno <- as.data.frame(geno,check.names = FALSE)
  geno <- round(geno,digits = 3)
  geno <- cbind(map,geno)
  write.table(geno,file,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
}


#' gemma2 : GEMMA programme in R (better version of gemma)
#'
#' @param covar : covariable et NULL marche
#' @param pheno : phenotype d'interet (doit etre une matrice)
#' @param geno : genotype --> one row per SNP, and one column per sample.
#' @param map : This file contains one line per SNP, with
#' @param three columns: (1) SNP label, (2) base-pair position, (3) chromosome.
#' @param gemma.exe : the executable gemma
#' @param opt.mat : option for the calcul of the related.matrix, "1" calculates the centered relatedness # matrix while "2" calculates the standardized relatedness matrix
#' @param met.reg : linear mixed model of Bayesian sparse linear mixed model.
#' @param lmm : if lmm choose option : where the  "1" performs Wald test "2" performs likelihood ratio test, "3" performs score test, "4" performs all the three tests.
#' @param bslmm : if bslmm choose option : "1" fits a standard linear BSLMM "2" fits a ridge regression/GBLUP, "3" fits a probit BSLMM.
#' @param pmin : for bslmm : initial no zero proportion of beta (pi)
#' @param pmax : for bslmm : max no zero proportion of beta (pi)
#' @param opt.bs : option for bslmm : w and s default as in Zeng 2017 (10000 for both)
#' @param smax : Number of variants with sparse effects (~ number of major effect loci). Default to 300 like in gemma
#' @return : a gemma result
#'
#' @details
#'
#' See paper of GEMMA
#'
#' @export
gemma2 <- function(geno, pheno, map, covar = NULL, gemma.exe = "./gemma.macosx",
                   opt.mat = c(1,2), met.reg = c("lmm", "bslmm"),
                   lmm = c(1,2,3,4), bslmm = c(1,2,3), pmin = 0.01, pmax = 0.05,
                   opt.bs = "-w  10000 -s  10000", smax = 300) {
  
  system("mkdir gemma")
  
  write.gemma.pheno("gemma/pheno.txt",pheno)
  write.gemma.covariates("gemma/covariates.txt",covar,pheno)
  write.gemma.geno("gemma/geno.txt",geno,map)
  write.gemma.map("gemma/map.txt",map)
  
  # related.matrix
  system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt -c gemma/covariates.txt -gk",
               opt.mat,"-o kin"),
         ignore.stdout = F)
  
  # lmm
  if (met.reg == "lmm") {
    print(met.reg)
    system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt -c gemma/covariates.txt",
                 " -a gemma/map.txt -k output/kin.cXX.txt -lmm", lmm,"-o lmm"),
           ignore.stdout = F)
    res <- read.table("output/lmm.assoc.txt", header = T)
  }
  # bslmm
  else {
    print(met.reg)
    # opt.bs <- "-w  10000 -s  10000" # Zeng 2017
    pmin <- log10(pmin)
    pmax <- log10(pmax)
    system(paste(gemma.exe,"-g gemma/geno.txt -p gemma/pheno.txt",
                 " -a gemma/map.txt -k output/kin.cXX.txt -bslmm", bslmm," -o bslmm",
                 opt.bs, "-pmin", pmin, "-pmax", pmax, "-smax", smax),
           ignore.stdout = F)
    res <- read.table("output/bslmm.param.txt", header = T)
  }
  return(res)
}
