rcura <- function(n, beta, z = NULL, q = NULL, alpha, sigma, C, model, dist) {
  
  
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != round(n)) {
    stop("n debe ser un numero entero positivo")
  }
  if (!is.numeric(beta) || length(beta) < 1) {
    stop("beta debe ser un vector numerico de largo >= 1")
  }
  if (!is.null(z)) {
    z <- as.matrix(z)
    if (nrow(z) != length(beta)) stop("z debe tener length(beta) filas (covariables en filas, individuos en columnas)")
    if (ncol(z) != n) stop("z debe tener n columnas")
  }
  if (!model %in% c("bernoulli", "poisson", "binneg", "geometrica", "polilog")) {
    stop("model debe ser 'bernoulli', 'poisson', 'binneg', 'geometrica' o 'polilog'")
  }
  if (!dist %in% c("weibull", "gamma", "bs")) {
    stop("dist debe ser 'weibull', 'gamma' o 'bs' (Birnbaum-Saunders)")
  }
  if (model %in% c("binneg", "polilog") && (is.null(q) || q <= 0)) {
    stop("para model='", model, "', q debe ser positivo y no puede ser NULL")
  }
  if (alpha <= 0) stop("alpha debe ser positivo")
  if (sigma <= 0) stop("sigma debe ser positivo")
  if (C <= 0) stop("C (censura) debe ser positivo")
  
  
  theta_link <- function(z, beta, model) {
    eta <- as.vector(crossprod(z, matrix(beta, ncol = 1)))
    switch(model,
           "bernoulli"  = exp(eta),
           "poisson"    = exp(eta),
           "binneg"     = stats::plogis(eta),
           "geometrica" = stats::plogis(eta),
           "polilog"    = stats::plogis(eta),
           stop("model no reconocido en theta_link")
    )
  }
  
  rpl <- function(n, theta, alpha_par) {
    polyloga <- function(theta, alpha) {
      exp(log(theta) + log1p(VGAM::lerch(x = theta, s = alpha, v = 1) - 1))
    }
    flag <- 0
    maximo <- 2000
    while (flag == 0) {
      m <- 1:maximo
      prob <- theta^m / (m^alpha_par * polyloga(theta, alpha_par))
      if (1 - sum(prob) < 1e-3) flag <- 1
      maximo <- maximo + 5000
    }
    sample(m, n, replace = TRUE, prob = prob)
  }
  
  sim_M <- function(theta, model, q) {
    n <- length(theta)
    switch(model,
           "bernoulli"  = stats::rbinom(n, size = 1, prob = theta / (1 + theta)),
           "poisson"    = stats::rpois(n, lambda = theta),
           "binneg"     = stats::rnbinom(n, size = q, prob = 1 - theta),
           "geometrica" = stats::rnbinom(n, size = 1, prob = 1 - theta),
           "polilog"    = {
             M <- integer(n)
             theta_unique <- unique(theta)
             for (th in theta_unique) {
               idx <- which(theta == th)
               m_vals <- rpl(length(idx), th, q)
               M[idx] <- m_vals - 1
             }
             M
           },
           stop("model no reconocido en sim_M")
    )
  }
  
  sim_W <- function(n, dist, alpha, sigma) {
    switch(dist,
           "weibull" = stats::rweibull(n, shape = alpha, scale = sigma),
           "gamma"   = stats::rgamma(n, shape = alpha, scale = sigma),
           "bs"      = VGAM::rbisa(n, shape = alpha, scale = sigma),
           stop("dist no reconocido en sim_W")
    )
  }
  
  
  r <- length(beta)
  if (is.null(z)) {
    z <- matrix(1, nrow = r, ncol = n)
    if (r > 1) {
      for (i in 2:r) z[i, ] <- stats::rbinom(n, 1, 0.5)
    }
  } else {
    z <- as.matrix(z)
  }
  
  theta <- theta_link(z, beta, model)
  M <- sim_M(theta, model, q)
  
  Tt <- rep(Inf, n)
  no_curados <- which(M > 0)
  
  if (length(no_curados) > 0) {
    for (i in no_curados) {
      W <- sim_W(M[i], dist, alpha, sigma)
      Tt[i] <- min(W)
    }
  }
  
  t_obs <- pmin(Tt, C)
  delta <- as.numeric(Tt <= C)
  
  out <- data.frame(t = t_obs, delta = delta)
  attr(out, "z") <- z
  attr(out, "theta") <- theta
  attr(out, "M") <- M
  attr(out, "p0_empirico") <- mean(M == 0)
  attr(out, "cens_pct") <- mean(delta == 0)
  out
}
