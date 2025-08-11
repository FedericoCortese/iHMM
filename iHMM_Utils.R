
iHmmNormalSampleBeam <- function(Y, hypers, numb=100, nums=1000, numi=5, S0) {
  # Y: numeric vector (length T)
  # hypers: list with either alpha0/gamma or alpha0_a,alpha0_b,gamma_a,gamma_b, and mu_0, sigma2_0, sigma2
  # numb: burn-in iterations
  # nums: number of posterior samples to keep
  # numi: thinning (iterations between kept samples)
  # S0: integer vector of initial states (length T)
  
  Tlen <- length(Y)
  
  sample <- list()
  
  if(is.null(S0)){
    S0=sample(1:5,Ylen)
  }
  
  sample$S <- as.integer(S0)
  sample$K <- max(sample$S)
  
  niter <- numb + (nums - 1L) * numi
  
  S_out <- list()
  stats <- list(
    K       = numeric(niter),
    alpha0  = numeric(niter),
    gamma   = numeric(niter),
    jml     = numeric(niter),  # kept for parity with MATLAB (unused here)
    jll     = numeric(niter),
    trellis = numeric(niter)
  )
  
  # Initialize hypers; resample a few times as our initial guess might be off.
  if (!is.null(hypers$alpha0)) {
    sample$alpha0 <- hypers$alpha0
  } else {
    # MATLAB: gamrnd(shape, scale). Here: rgamma(shape=., rate=.)
    sample$alpha0 <- rgamma(1L, shape = hypers$alpha0_a, rate = hypers$alpha0_b)
  }
  if (!is.null(hypers$gamma)) {
    sample$gamma <- hypers$gamma
  } else {
    sample$gamma <- rgamma(1L, shape = hypers$gamma_a, rate = hypers$gamma_b)
  }
  
  for (i in 1:5) {
    sample$Beta <- rep(1 / (sample$K + 1L), sample$K + 1L)
    tmp <- iHmmHyperSample(sample$S, sample$Beta, sample$alpha0, sample$gamma, hypers, 20L)
    sample$Beta   <- tmp$Beta
    sample$alpha0 <- tmp$alpha0
    sample$gamma  <- tmp$gamma
  }
  
  # Sample the emission and transition probabilities.
  sample$Mus <- SampleNormalMeans_R(sample$S, Y, sample$K, hypers$sigma2, hypers$mu_0, hypers$sigma2_0)
  sample$Pi  <- SampleTransitionMatrix(sample$S, sample$alpha0 * sample$Beta)
  # Remove the (K+1)-th row as in MATLAB
  sample$Pi  <- sample$Pi[-(sample$K + 1L), , drop = FALSE]
  
  iter <- 1L
  cat(sprintf("Iteration 0: K = %d, alpha0 = %f, gamma = %f.\n", sample$K, sample$alpha0, sample$gamma))
  
  while (iter <= niter) {
    
    # Reset trellis size count
    stats$trellis[iter] <- 0
    
    # Sample auxiliary variables u
    u <- numeric(Tlen)
    for (t in 1:Tlen) {
      if (t == 1L) {
        u[t] <- runif(1L) * sample$Pi[1L, sample$S[t]]
      } else {
        u[t] <- runif(1L) * sample$Pi[sample$S[t - 1L], sample$S[t]]
      }
    }
    
    # Extend Pi and Mus (and break Beta stick) while needed
    while (max(sample$Pi[, ncol(sample$Pi)]) > min(u)) {
      
      pl <- ncol(sample$Pi)
      bl <- length(sample$Beta)
      
      # Safety check: #columns == length(Beta)
      stopifnot(bl == pl)
      
      # Add a row to transition matrix (Dirichlet with alpha0 * Beta)
      sample$Pi <- rbind(sample$Pi, dirichlet_sample(sample$alpha0 * sample$Beta))
      
      # Add a mean for the new state
      sample$Mus <- c(sample$Mus, rnorm(1L, mean = hypers$mu_0, sd = sqrt(hypers$sigma2_0)))
      
      # Break Beta stick (split the last mass)
      be <- sample$Beta[length(sample$Beta)]
      bg <- rbeta(1L, 1, sample$gamma)
      sample$Beta[length(sample$Beta)] <- bg * be
      sample$Beta <- c(sample$Beta, (1 - bg) * be)
      
      # Split the last column of Pi according to pg
      pe <- sample$Pi[, ncol(sample$Pi)]
      a  <- rep(sample$alpha0 * sample$Beta[length(sample$Beta) - 1L], bl)
      b  <- sample$alpha0 * (1 - sum(sample$Beta[1:(length(sample$Beta) - 1L)]))
      pg <- betarnd_vec(a, b)  # vectorized Beta draws
      
      if (is.nan(sum(pg))) {   # approximation when a or b are really small
        p <- a / (a + b)
        pg <- rbinom(length(p), size = 1L, prob = p)
      }
      
      # Update existing last column and append a new split column
      sample$Pi[, pl] <- pg * pe
      sample$Pi <- cbind(sample$Pi, (1 - pg) * pe)
    }
    
    sample$K <- nrow(sample$Pi)
    
    # Safety checks
    stopifnot(sample$K == (length(sample$Beta) - 1L))
    stopifnot(sample$K == length(sample$Mus))
    
    # Forward dynamic program on truncated trellis
    dyn_prog <- matrix(0, nrow = sample$K, ncol = Tlen)
    
    dyn_prog[, 1L] <- as.numeric(sample$Pi[1L, 1:sample$K] > u[1L])
    stats$trellis[iter] <- stats$trellis[iter] + sum(dyn_prog[, 1L])
    
    for (k in 1:sample$K) {
      dyn_prog[k, 1L] <- exp(-0.5 * (Y[1L] - sample$Mus[k]) * (Y[1L] - sample$Mus[k]) / hypers$sigma2) * dyn_prog[k, 1L]
    }
    dyn_prog[, 1L] <- dyn_prog[, 1L] / sum(dyn_prog[, 1L])
    
    for (t in 2:Tlen) {
      A <- (sample$Pi[1:sample$K, 1:sample$K, drop = FALSE] > u[t])
      dyn_prog[, t] <- t(A) %*% dyn_prog[, t - 1L]
      stats$trellis[iter] <- stats$trellis[iter] + sum(A)
      
      for (k in 1:sample$K) {
        dyn_prog[k, t] <- exp(-0.5 * (Y[t] - sample$Mus[k]) * (Y[t] - sample$Mus[k]) / hypers$sigma2) * dyn_prog[k, t]
      }
      dyn_prog[, t] <- dyn_prog[, t] / sum(dyn_prog[, t])
    }
    
    # Backtrack to sample a path
    if (sum(dyn_prog[, Tlen]) != 0 && is.finite(sum(dyn_prog[, Tlen]))) {
      sample$S[Tlen] <- sample_categorical(dyn_prog[, Tlen])
      
      for (t in (Tlen - 1L):1L) {
        r <- dyn_prog[, t] * as.numeric(sample$Pi[, sample$S[t + 1L]] > u[t + 1L])
        r <- r / sum(r)
        sample$S[t] <- sample_categorical(r)
      }
      
      # Cleanup: remove unused states
      zind <- setdiff(seq_len(sample$K), unique(sample$S))
      zind <- sort(zind, decreasing = TRUE)
      if (length(zind) > 0) {
        for (i_rm in zind) {
          sample$Beta[length(sample$Beta)] <- sample$Beta[length(sample$Beta)] + sample$Beta[i_rm]
          sample$Beta <- sample$Beta[-i_rm]
          sample$Pi   <- sample$Pi[, -i_rm, drop = FALSE]
          sample$Pi   <- sample$Pi[-i_rm, , drop = FALSE]
          sample$Mus  <- sample$Mus[-i_rm]
          sample$S[sample$S > i_rm] <- sample$S[sample$S > i_rm] - 1L
        }
        sample$K <- nrow(sample$Pi)
      }
      
      # Resample Beta (and optionally alpha0, gamma)
      tmp <- iHmmHyperSample(sample$S, sample$Beta, sample$alpha0, sample$gamma, hypers, 20L)
      sample$Beta   <- tmp$Beta
      sample$alpha0 <- tmp$alpha0
      sample$gamma  <- tmp$gamma
      
      # Resample emissions
      sample$Mus <- SampleNormalMeans_R(sample$S, Y, sample$K, hypers$sigma2, hypers$mu_0, hypers$sigma2_0)
      
      # Resample transitions
      sample$Pi <- SampleTransitionMatrix(sample$S, sample$alpha0 * sample$Beta)
      sample$Pi <- sample$Pi[-(sample$K + 1L), , drop = FALSE]
      
      # Safety checks
      stopifnot(nrow(sample$Pi) == sample$K)
      stopifnot(ncol(sample$Pi) == sample$K + 1L)
      stopifnot(sample$K == (length(sample$Beta) - 1L))
      stopifnot(min(sample$Pi) >= 0)
      stopifnot(sample$K == max(sample$S))
      
      # Stats + logging
      stats$alpha0[iter] <- sample$alpha0
      stats$gamma[iter]  <- sample$gamma
      stats$K[iter]      <- sample$K
      
      stats$jll[iter] <- iHmmNormalJointLogLikelihood(
        sample$S, Y, sample$Beta, sample$alpha0, hypers$mu_0, hypers$sigma2_0, hypers$sigma2
      )
      
      cat(sprintf("Iteration: %d: K = %d, alpha0 = %f, gamma = %f, JL = %f.\n",
                  iter, sample$K, sample$alpha0, sample$gamma, stats$jll[iter]))
      
      # Save a sample (post burn-in, every numi)
      if (iter >= numb && ((iter - numb) %% numi == 0L)) {
        # Deep copy (to be safe)
        S_out[[length(S_out) + 1L]] <- list(
          S     = sample$S,
          K     = sample$K,
          Beta  = sample$Beta,
          Pi    = sample$Pi,
          Mus   = sample$Mus,
          alpha0 = sample$alpha0,
          gamma  = sample$gamma
        )
      }
      
      iter <- iter + 1L
      
    } else {
      cat("Wasted computation as there were no paths through the iHMM.\n")
      # (loop continues without incrementing iter)
    }
  }
  
  list(S = S_out, stats = stats)
}

# ---- Helpers ----

# Dirichlet sampler (returns a probability vector)
dirichlet_sample <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

# Vectorized Beta random draws with (possibly) vector shape1 and scalar/vec shape2
betarnd_vec <- function(a, b) {
  # a: numeric vector; b: scalar or vector of same length
  if (length(b) == 1L) b <- rep(b, length(a))
  rbeta(length(a), shape1 = a, shape2 = b)
}

# Categorical sampler: returns index in 1..length(p)
sample_categorical <- function(p) {
  p <- as.numeric(p)
  p <- p / sum(p)
  # Equivalent to: 1 + sum(runif(1) > cumsum(p))
  which.max(runif(1L) <= cumsum(p))
}

# SampleNormalMeans (R version)
SampleNormalMeans_R <- function(S, Y, K, sigma2, mu_0, sigma2_0) {
  Mus <- numeric(K)
  for (k in 1:K) {
    Ys <- Y[S == k]
    N  <- length(Ys)
    if (N == 0L) {
      Mus[k] <- rnorm(1L, mean = mu_0, sd = sqrt(sigma2_0))
    } else {
      Mu_ml <- sum(Ys) / N
      Mu_n  <- (sigma2 * mu_0 + N * sigma2_0 * Mu_ml) / (N * sigma2_0 + sigma2)
      S2_n  <- 1 / (1 / sigma2_0 + N / sigma2)
      Mus[k] <- rnorm(1L, mean = Mu_n, sd = sqrt(S2_n))
    }
  }
  Mus
}

# ------ iHmmHyperSample ------
iHmmHyperSample <- function(S, ibeta, ialpha0, igamma, hypers, numi) {
  K <- length(ibeta) - 1L
  Tlen <- length(S)
  
  # N: transition counts (K+1 x K+1 so we can index safely)
  N <- matrix(0, nrow = K, ncol = K)
  # Count "start" -> S[1]
  N[1L, S[1L]] <- N[1L, S[1L]] + 1L
  for (t in 2:Tlen) {
    N[S[t - 1L], S[t]] <- N[S[t - 1L], S[t]] + 1L
  }
  
  # M: number of tables (K x K)
  M <- matrix(0, nrow = K, ncol = K)
  for (j in 1:K) {
    for (k in 1:K) {
      if (N[j, k] > 0) {
        for (ell in 1:N[j, k]) {
          p <- (ialpha0 * ibeta[k]) / (ialpha0 * ibeta[k] + ell - 1)
          M[j, k] <- M[j, k] + rbinom(1L, size = 1L, prob = p)
        }
      }
    }
  }
  
  # Resample Beta ~ Dirichlet([m_.k]_{k=1..K}, gamma)
  m_dotk <- colSums(M)
  ibeta  <- dirichlet_sample(c(m_dotk, igamma))
  
  # Resample alpha0 (Escobar & West style for HDP rows)
  if (!is.null(hypers$alpha0)) {
    ialpha0 <- hypers$alpha0
  } else {
    for (iter in 1:numi) {
      Nj <- rowSums(N)                          # length K
      w  <- rbeta(K, shape1 = ialpha0 + 1, shape2 = Nj)
      p  <- Nj / ialpha0
      p  <- p / (p + 1)
      s  <- rbinom(K, size = 1L, prob = p)
      shape <- hypers$alpha0_a + sum(M) - sum(s)
      rate  <- hypers$alpha0_b - sum(log(w))
      ialpha0 <- rgamma(1L, shape = shape, rate = rate)
    }
  }
  
  # Resample gamma (Escobar & West 1995)
  if (!is.null(hypers$gamma)) {
    igamma <- hypers$gamma
  } else {
    k <- length(ibeta)        # == K+1
    m <- sum(M)
    for (iter in 1:numi) {
      mu    <- rbeta(1L, shape1 = igamma + 1, shape2 = m)
      pi_mu <- 1 / (1 + (m * (hypers$gamma_b - log(mu))) / (hypers$gamma_a + k - 1))
      shape <- if (runif(1L) < pi_mu) hypers$gamma_a + k else hypers$gamma_a + k - 1
      rate  <- hypers$gamma_b - log(mu)
      igamma <- rgamma(1L, shape = shape, rate = rate)
    }
  }
  
  list(Beta = ibeta, alpha0 = ialpha0, gamma = igamma, N = N, M = M)
}

# ---- SampleTransitionMatrix ------

SampleTransitionMatrix <- function(S, H) {
  K <- length(H)
  Tlen <- length(S)
  
  N <- matrix(0, nrow = K, ncol = K)
  for (t in 2:Tlen) {
    N[S[t - 1L], S[t]] <- N[S[t - 1L], S[t]] + 1L
  }
  
  Pi <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    Pi[k, ] <- dirichlet_sample(N[k, ] + H)
  }
  Pi
}

# ----- iHmmNormalJointLogLikelihood----
iHmmNormalJointLogLikelihood <- function(S, Y, Beta, alpha0, mu_0, sigma2_0, sigma2) {
  K <- max(S)
  Tlen <- length(S)
  
  N <- matrix(0, nrow = K, ncol = K)
  E <- matrix(0, nrow = K, ncol = 3)
  
  # Initial contributions
  N[1L, S[1L]] <- N[1L, S[1L]] + 1L
  E[S[1L], 1] <- E[S[1L], 1] + Y[1L]
  E[S[1L], 2] <- E[S[1L], 2] + 1
  E[S[1L], 3] <- E[S[1L], 3] + Y[1L]^2
  
  for (t in 2:Tlen) {
    N[S[t - 1L], S[t]] <- N[S[t - 1L], S[t]] + 1L
    E[S[t], 1] <- E[S[t], 1] + Y[t]
    E[S[t], 2] <- E[S[t], 2] + 1
    E[S[t], 3] <- E[S[t], 3] + Y[t]^2
  }
  
  logp <- 0
  
  for (k in 1:K) {
    R  <- c(N[k, ], 0) + alpha0 * Beta
    ab <- alpha0 * Beta
    nz <- which(R != 0)
    
    # Transition part: log Dirichlet-multinomial row likelihood
    logp <- logp +
      lgamma(alpha0) -
      lgamma(sum(c(N[k, ], 0)) + alpha0) +
      sum(lgamma(R[nz])) -
      sum(lgamma(ab[nz]))
    
    # Emission part: Normal-Normal conjugacy
    sigma2_n <- 1 / (1 / sigma2_0 + E[k, 2] / sigma2)
    mu_n     <- (mu_0 / sigma2_0 + E[k, 1] / sigma2) * sigma2_n
    
    logp <- logp +
      0.5 * log(sigma2_n) -
      0.5 * (log(sigma2_0) + (E[k, 2] / 2) * log(2 * pi * sigma2)) -
      0.5 * (E[k, 3] / sigma2 + mu_0^2 / sigma2_0 - mu_n^2 / sigma2_n)
  }
  
  as.numeric(logp)
}



# Simulation --------------------------------------------------------------

sim_data_stud_t=function(seed=123,
                         TT,
                         P,
                         Ktrue=3,
                         mu=1.5,
                         rho=0,
                         nu=4,
                         pers=.95){
  
  
  MU=seq(mu, -mu, length.out=Ktrue)
  
  # Markov chain simulation
  x <- numeric(TT)
  Q <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(Q)=rep(pers,Ktrue)
  init <- rep(1/Ktrue,Ktrue)
  set.seed(seed)
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for(i in 2:TT){
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Continuous variables simulation
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, TT, P * Ktrue)
  SimData = matrix(0, TT, P)
  
  set.seed(seed)
  for(k in 1:Ktrue){
    # u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    u = mvtnorm::rmvt(TT, sigma = (nu-2)*Sigma/nu, df = nu, delta = rep(MU[k],P))
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
  }
  
  SimData=data.frame(SimData)
  
  return(list(
    SimData=SimData,
    mchain=x,
    TT=TT,
    P=P,
    K=Ktrue,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
}
