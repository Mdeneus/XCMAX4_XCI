
###########################################################################

 #compute  XCI status based on XCMAX4 function 

###########################################################################


 XCMAX4_XCI <- function(data) {

  Phen <- as.matrix(data[, 1])
  Snp <- as.matrix(data[, 2])
  X <- as.matrix(data[, -c(1, 2)])
  Sex <- as.matrix(data["Sex"])
  n <- length(Phen)

  # Indicator variables
  Ind_AA <- rep(0, n)
  Ind_Aa <- rep(0, n)
  Ind_AO <- rep(0, n)
  Ind_AA[Sex == 1 & Snp == 2] <- 1
  Ind_Aa[Sex == 1 & Snp == 1] <- 1
  Ind_AO[Sex == 0 & Snp == 2] <- 1

  # Logistic regression coefficients
  esta <- summary(glm(Phen ~ X, family = binomial(link = "logit")))$coefficients[, 1]

  # Helper functions
  estf <- function(esta, X) {
    1 / (1 + exp(-X %*% esta))
  }
  infora <- function(estpen, X) {
    l <- dim(X)[2]
    Ia <- matrix(0, nrow = l, ncol = l)
    for (i in 1:l) {
      for (j in 1:l) {
        Ia[i, j] <- sum(X[, i] * X[, j] * (1 - estpen) * estpen)
      }
    }
    return(Ia)
  }
  inforb <- function(estpen, z1, z2, X, Ind_AO, Ind_Aa, Ind_AA) {
    Snp <- 2 * Ind_AA + z1 * Ind_Aa + z2 * Ind_AO
    Ib <- sum(Snp * Snp * (1 - estpen) * estpen)
    return(Ib)
  }
  inforba <- function(estpen, z1, z2, X, Ind_AO, Ind_Aa, Ind_AA) {
    l <- dim(X)[2]
    Snp <- 2 * Ind_AA + z1 * Ind_Aa + z2 * Ind_AO
    Iba <- NULL
    for (i in 1:l) {
      Iba[i] <- sum(X[, i] * Snp * (1 - estpen) * estpen)
    }
    return(Iba)
  }

  # Calculate estimated probabilities
  X <- cbind(rep(1, n), X)
  estpen <- estf(esta, X)

  # Compute variance components
  Ia <- infora(estpen, X)
  Ib11 <- inforb(estpen, 1, 1, X, Ind_AO, Ind_Aa, Ind_AA)
  Ib02 <- inforb(estpen, 0, 2, X, Ind_AO, Ind_Aa, Ind_AA)
  Ib12 <- inforb(estpen, 1, 2, X, Ind_AO, Ind_Aa, Ind_AA)
  Ib22 <- inforb(estpen, 2, 2, X, Ind_AO, Ind_Aa, Ind_AA)
  Iba11 <- inforba(estpen, 1, 1, X, Ind_AO, Ind_Aa, Ind_AA)
  Iba02 <- inforba(estpen, 0, 2, X, Ind_AO, Ind_Aa, Ind_AA)
  Iba12 <- inforba(estpen, 1, 2, X, Ind_AO, Ind_Aa, Ind_AA)
  Iba22 <- inforba(estpen, 2, 2, X, Ind_AO, Ind_Aa, Ind_AA)

  stest <- function(estpen, Phen, z1, z2, Ind_AO, Ind_Aa, Ind_AA, Ib, Iba, Ia) {
    Snp <- 2 * Ind_AA + z1 * Ind_Aa + z2 * Ind_AO
    score <- sum(Snp * (Phen - estpen))
    variance <- Ib - t(as.matrix(Iba)) %*% solve(Ia) %*% as.matrix(Iba)
    re <- score / sqrt(variance)
    return(re)
  }

  # Compute test statistics
  s11 <- stest(estpen, Phen, 1, 1, Ind_AO, Ind_Aa, Ind_AA, Ib11, Iba11, Ia)
  s02 <- stest(estpen, Phen, 0, 2, Ind_AO, Ind_Aa, Ind_AA, Ib02, Iba02, Ia)
  s12 <- stest(estpen, Phen, 1, 2, Ind_AO, Ind_Aa, Ind_AA, Ib12, Iba12, Ia)
  s22 <- stest(estpen, Phen, 2, 2, Ind_AO, Ind_Aa, Ind_AA, Ib22, Iba22, Ia)

  # Determine zmax1
  zmax1 <- max(abs(s11), abs(s02), abs(s12), abs(s22))

  classify_xci <- function(z1, z2) {
    if (z1 == 1 && z2 == 1) {
      return("XCI-E (Escape)")
    } else if (z1 == 1 && z2 == 2) {
      return("XCI-R (Random)")
    } else if (z1 == 0 && z2 == 2) {
      return("XCI-s (SN)")  # Semi-null
    } else if (z1 == 2 && z2 == 2) {
      return("XCI-s (SR)")  # Semi-recurrent
    } else {
      return("Undefined")
    }
  }


  if (zmax1 == abs(s11)) {
    xci_type <- classify_xci(1, 1)
  } else if (zmax1 == abs(s02)) {
    xci_type <- classify_xci(0, 2)
  } else if (zmax1 == abs(s12)) {
    xci_type <- classify_xci(1, 2)
  } else if (zmax1 == abs(s22)) {
    xci_type <- classify_xci(2, 2)
  } else {
    xci_type <- "Undefined"
  }

  # Output results
  return(list(
    "statistic" = zmax1,
    "xci_type" = xci_type
  ))
}

##########################
