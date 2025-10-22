
# Sim function ----
dat.sim.zim <- function(
    x.d,               # Predictors for log link, if missing, then "x "is generated based on "px" with level "0", "1" etc.
    x.z,        # Predictors for logit link
    coef.d = NULL,     # coefficients of log model, including intercept
    coef.z = NULL,     # coefficients of logit model, including intercept
    dist.m = c("ZINB", "ZIP", "NB", "Poi", "ZILN"),
    dist.y,
    sd = NULL,
    theta = NULL,     # shape parameter of NB distribution, Var = u + (u^2)/theta
    coef.y = 0
)
{
  # ------------------------------------------------------------------------------------------------
  ## Simulate exposure variable
  # ------------------------------------------------------------------------------------------------
  px = mean(x.d[,1])
  x.d <- as.matrix(cbind(1, x.d) )      # combine all predictors in M model (both logit and log models)
  coef.d <- as.matrix(coef.d)

  if (!(ncol(x.d) == length(coef.d))){
    stop("Number of coefficients and predictors not match in log link")
  }

  # ------------------------------------------------------------------------------------------------
  ## Simulate Covariates
  # ------------------------------------------------------------------------------------------------
  if (missing(x.z)){
    x.z <- x.d
  }else{
      x.z = as.matrix(cbind(1, x.z))
    }

  if (missing(coef.z)) coef.z = NULL
  if (!(ncol(x.z) == length(coef.z))){
      stop("Number of coefficients and predictors not match in logit link")
    }

  coef.d <- as.matrix(coef.d)
  coef.z <- as.matrix(coef.z)

  # ------------------------------------------------------------------------------------------------
  ## Simulate Mediator
  # ------------------------------------------------------------------------------------------------
  varx = x.d[,2]
  x.lvl = nlevels(as.factor(varx))

  ## Mediator Simulation Function
  msim <- function(x.lvl){

    x.d <- x.d[varx == unique(varx)[x.lvl], ]
    if (is.null(coef.z)){
      x.z = 0
    } else if (!is.null(coef.z)){
      x.z <- x.z[varx == unique(varx)[x.lvl], ]
    }

    ## calculate the excess-zero proportion given each X-level
    etaz <-  x.z %*% coef.z
    m.normal <- rnorm(n, -etaz, 1.6)
    p.zero = exp(etaz) / (1 + exp(etaz))
    p.zero = round(mean(p.zero),2)
    quantiles <- quantile(m.normal, p.zero)
    m.z <-  as.numeric( factor(cut(m.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1

    etad <- x.d %*% coef.d
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
  # colnames(dfm) <- c("int", "x", "mnb", colnames(cov_var))

  ## -----------------------------------------------------------------------------------------------
  ## Simulate Outcome
  ## -----------------------------------------------------------------------------------------------
  xm <- c(dfm[,"x"]*dfm[,"mnb"])
  dfm <- data.frame(cbind(dfm, xm)) %>%
    mutate(im = ifelse(mnb == 0, 0, 1)) %>% mutate(xi = x*im) %>%
    as.matrix()
  expit <- function(x) exp(x)/(1 + exp(x))

  if (!(length(coef.y) == ncol(dfm))){
    stop("Coefficients in y model does not match number of predictors")
  }

  if (missing(dist.y)) dist.y = "gaussian"
  if (missing(sd) & dist.y == "gaussian") sd = 1

  if (dist.y == "gaussian"){
    y = rnorm(n, mean = dfm %*% coef.y, sd = sd)
  }else if (dist.y == "binom"){

      py0 <- expit(dfm[x==0,] %*% coef.y)
      py1 <- expit(dfm[x==1,] %*% coef.y)

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
  dat <- data.frame(cbind(dfm, y))

  return(dat)
}

