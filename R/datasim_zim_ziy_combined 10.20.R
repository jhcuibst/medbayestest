

# data sim M
dat.sim.zim <- function(
    x.d,               # Predictors for log link, if missing, then "x "is generated based on "px" with level "0", "1" etc.
    x.z,               # Predictors for logit link
    coef.d = NULL,     # coefficients of log model, including intercept
    coef.z = NULL,     # coefficients of logit model, including intercept
    dist.m = c("ZINB", "ZIP", "NB", "Poi", "ZILN"),
    sd = NULL,
    theta = NULL     # shape parameter of NB distribution, Var = u + (u^2)/theta
)
{
  # ------------------------------------------------------------------------------------------------
  ## Simulate exposure variable
  # ------------------------------------------------------------------------------------------------
  if (missing(px) & missing(x.d)) {
    warning("Note: a bianry distribution and p=0.5 will be used by default")
    px = c(0.5, 0.5)                       # binary exposure with p=0.5 by default
    x.d <- c(rep(0, n/2), rep(1, n/2))     # generate a binary exposure with p=0.5 by default
  }
  if (!missing(px) & missing(x.d)){
    x.d <- list()
    for (i in 1:length(px)){
      x.d[[i]] <- rep((i-1), n*px[i])
    }
    x.d <- unlist(x.d)
  }

  x.d <- as.matrix(cbind(1, x.d) )      # combine all predictors in M model (both logit and log models)
  coef.d <- as.matrix(coef.d)

  if (!(ncol(x.d) == length(coef.d))){
    stop("Number of coefficients and predictors not match in log link")
  }

  # ------------------------------------------------------------------------------------------------
  ## Simulate Covariates
  # ------------------------------------------------------------------------------------------------
  if (missing(x.z) & is.null(p.zero)) x.z <- x.d
  if (missing(x.z) & !is.null(p.zero)) x.z <- 0
  if (missing(coef.z)) coef.z = NULL

  if (is.null(coef.z)){
    x.z = 0
  } else {
    x.z <- x.z
    coef.z <- as.matrix(coef.z)
    if (!(ncol(x.z) == length(coef.z))){
      stop("Number of coefficients and predictors not match in logit link")
    }
  }

  coef.d <- as.matrix(coef.d)
  coef.z <- as.matrix(coef.z)

  # ------------------------------------------------------------------------------------------------
  ## Simulate Mediator
  # ------------------------------------------------------------------------------------------------

  if (missing(p.zero)) p.zero = NULL
  if (missing(py)) py = NULL

  varx = x.d[,2]
  x.lvl = nlevels(as.factor(varx))

  if (!is.null(p.zero)){
    if (length(p.zero) == 1 & x.lvl > 1){
      p.zero[1:x.lvl] = p.zero[1]
    } else if(length(p.zero) != x.lvl & length(p.zero) >1 ) {
      stop("Number of 'p.zero' not match 'X levels' ")
    }
  }

  ## Mediator Simulation Function
  msim <- function(x.lvl){

    if (is.null(coef.z) & is.null(p.zero)){
      x.z = 0
    } else if (!is.null(coef.z) & is.null(p.zero)){
      x.z <- x.z[varx == unique(varx)[x.lvl], ]
    } else if (is.null(coef.z) & !is.null(p.zero)){
      x.z = 0
    }
    x.d <- x.d[varx == unique(varx)[x.lvl], ]

    n = n*px[x.lvl]
    ind <- sort(rep(1:n.ind, n))

    b <- rep(NA, n.ind)                  # random effect: b ~ N(0, tau.z^2)
    for (j in 1:n.ind) b[j] <- rnorm(1, 0, tau.z)

    ## calculate the excess-zero proportion given each X-level

    if (is.null(p.zero)){
      etaz <- b[ind] + x.z %*% coef.z
      m.normal <- rnorm(n, -etaz, 1.6)
      p.zero = exp(g0+gx*x.d[,2][x.lvl]+coef.z[-(1:2)]*x.d[,-(1:2)])/(1+exp(g0+gx*x.d[,2][x.lvl]+coef.z[-(1:2)]*x.d[,-(1:2)]))
      p.zero = round(mean(p.zero),2)
      quantiles <- quantile(m.normal, p.zero)
      m.z <-  as.numeric( factor(cut(m.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1
    } else{
      p.zero = p.zero[x.lvl]
      m.z = rep(1, length(x.d[,1]))
    }

    ## simulate M given zero proportion
    b <- rep(NA, n.ind)                  # random effect: b ~ N(0, tau.c^2)
    for (j in 1:n.ind) b[j] <- rnorm(1, 0, tau.d)

    etad <- b[ind] + x.d %*% coef.d
    etad <- mean(etad)
    if (dist.m %in% c("ZINB", "NB") )  m.nb <- rnbinom(n, mu = exp(etad), size = theta)
    if (dist.m %in% c("ZIP", "Poi") )  m.nb <- rpois(n, lambda =  exp(etad))

    mnb <- ifelse(m.z == 0, 0, m.nb)

    return(mnb)
  }

  ls <- list()
  for(i in 1:x.lvl) ls[[i]] <- msim(x.lvl = i)
  comb_m <- c(ls[[1]], ls[[2]])
  dfm <- cbind(int = 1, x = x.d[1:n,2], mnb = comb_m, c = x.d[1:n,-(1:2)])
  colnames(dfm) <- c("int", "x", "mnb", colnames(cov_var))
  return(dfm)
}

# dat sim Y ----
## If not ZIM ----
dat.sim.ziy <- function(
    x.d,               # Predictors for log link, if missing, then "x "is generated based on "px" with level "0", "1" etc.
    x.z,               # Predictors for logit link
    coef.d = NULL,     # coefficients of log model, including intercept
    coef.z = NULL,     # coefficients of logit model, including intercept
    dist.y = c("ZINB", "ZIP", "NB", "Poi", "ZILN"),
    sd = NULL,
    theta = NULL     # shape parameter of NB distribution, Var = u + (u^2)/theta
)
{
  # set.seed(seed)
  if (!(ncol(x.d) == length(coef.d))){
    stop("Number of coefficients and predictors not match in log link. Coefficients should include the intercept.")
  }

  # ------------------------------------------------------------------------------------------------
  ## Simulate Covariates
  # ------------------------------------------------------------------------------------------------
  if (missing(x.z)) x.z <- x.d
  if (missing(coef.z)) coef.z = NULL

  coef.d <- as.matrix(coef.d)
  coef.z <- as.matrix(coef.z)

  # ------------------------------------------------------------------------------------------------
  ## Simulate Outcome
  # ------------------------------------------------------------------------------------------------

  etaz <- x.z %*% as.matrix(coef.z)
  m.normal <- rnorm(n, -etaz, 1.6)
  p.zero = exp(etaz)/ (1+exp(etaz))
  p.zero = round(mean(p.zero),2)
  quantiles <- quantile(m.normal, p.zero)
  m.z <-  as.numeric( factor(cut(m.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1

  etad <- x.d %*% coef.d
  if (dist.y %in% c("ZINB", "NB") )  m.nb <- rnbinom(n, mu = exp(etad), size = theta)
  if (dist.y %in% c("ZIP", "Poi") )  m.nb <- rpois(n, lambda =  exp(etad))
  if (dist.y %in% c("ZILN") )        m.nb <- exp(rnorm(n, mean = etad, sd = sd))
  # if (dist %in% c("ZILN") )        m.nb <- rlnorm(n, meanlog = etad, sdlog = exp(sd))

  zi.y = ifelse(m.z == 0, 0, m.nb)

  df.y <- data.frame(x = x, m = m, c = c, zi.y = zi.y)
  df.y$ind.y = ifelse(df.y$zi.y == 0, 0, 1)

  return(df.y)

}

# dat sim Y ----
## If ZIM simulated from dat.sim.zim ----
dat.sim.ziy <- function(
    x.d,               # Predictors for log link, if missing, then "x "is generated based on "px" with level "0", "1" etc.
    x.z,               # Predictors for logit link
    coef.d = NULL,     # coefficients of log model, including intercept
    coef.z = NULL,     # coefficients of logit model, including intercept
    dist.y = c("ZINB", "ZIP", "NB", "Poi", "ZILN"),
    sd = NULL,
    theta = NULL     # shape parameter of NB distribution, Var = u + (u^2)/theta
)
{
  # set.seed(seed)
  if (!(ncol(dfm) == length(coef.d))){
    stop("Number of coefficients and predictors not match in log link. Coefficients should include the intercept.")
  }

  # ------------------------------------------------------------------------------------------------
  ## Simulate Covariates
  # ------------------------------------------------------------------------------------------------
  x.d <- as.matrix(dfm)
  if (missing(x.z)) x.z <- x.d
  if (missing(coef.z)) coef.z = NULL

  coef.d <- as.matrix(coef.d)
  coef.z <- as.matrix(coef.z)

  # ------------------------------------------------------------------------------------------------
  ## Simulate Outcome
  # ------------------------------------------------------------------------------------------------

  etaz <- x.z %*% as.matrix(coef.z)
  m.normal <- rnorm(n, -etaz, 1.6)
  p.zero = exp(etaz)/ (1+exp(etaz))
  p.zero = round(mean(p.zero),2)
  quantiles <- quantile(m.normal, p.zero)
  m.z <-  as.numeric( factor(cut(m.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1

  etad <- x.d %*% coef.d
  if (dist.y %in% c("ZINB", "NB") )  m.nb <- rnbinom(n, mu = exp(etad), size = theta)
  if (dist.y %in% c("ZIP", "Poi") )  m.nb <- rpois(n, lambda =  exp(etad))
  if (dist.y %in% c("ZILN") )        m.nb <- exp(rnorm(n, mean = etad, sd = sd))
  # if (dist %in% c("ZILN") )        m.nb <- rlnorm(n, meanlog = etad, sdlog = exp(sd))
  zi.y = ifelse(m.z == 0, 0, m.nb)

  if (ydist == "gaussian"){
    y = rnorm(n, mean = dfm %*% coef.y, sd = sd)
  }else if (ydist == "binom"){
    if ( is.null(py) ){
      py0 <- expit(dfm[x==0,] %*% coef.y)
      py1 <- expit(dfm[x==1,] %*% coef.y)
    } else if ( !(is.null(py)) ){
      py0 <- py[1]
      py1 <- py[2]
    }
    p0 = mean(py0);     p1 = mean(py1)
    eta_0 <- dfm[x==0,] %*% coef.y
    y.normal <- rnorm(n/2, -eta_0, 1.6)
    quantiles <- quantile(y.normal, p0)
    y0 <- 2- as.numeric( factor(cut(y.normal, breaks = c(-Inf, quantiles, Inf))) )

    eta_1 <- dfm[x==1,] %*% coef.y
    y.normal <- rnorm(n/2, -eta_1, 1.6)
    quantiles <- quantile(y.normal, p1)
    y1 <- 2- as.numeric( factor(cut(y.normal, breaks = c(-Inf, quantiles, Inf))) )
    y <- c(y0, y1)
  }

  df.y <- data.frame(cbind(dfm, y))

  return(df.y)
}










