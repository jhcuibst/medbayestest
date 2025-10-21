
# devtools::install_github("jhcuibst/medbayestest")

library(medbayestest)

medbayes_sens <- function(object, rho = 0.3, rho01 = 0, rho0Y = 0, ...){

  if (!inherits(object$params, "params") )
    stop("Input must be a 'mediation_result' object from mediation_analysis().")
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("This function requires R package MASS. Please make sure it's loaded.")



  # Reuse existing components
  model.m <- object$params$model.m
  model.y <- object$params$model.y
  treat <- object$params$treat
  mediator <- object$params$mediator
  ind_mediator <- object$params$ind_mediator
  outcome <- object$params$outcome
  control.value <- object$params$control.value
  treat.value <- object$params$treat.value

  sens.out <- matrix(nrow = length(rho), ncol = 3, 3)

  for(i in (1:length(rho))){
    rho_i = rho[i]
    rho01 = rho01
    rho0Y = rho0Y

    results <- run_sensitivity(model.m = model.m, model.y = model.y,
                               treat = treat, mediator = mediator, ind_mediator = ind_mediator, outcome = outcome,
                               control.value = control.value, treat.value = treat.value,
                               rho01 = rho01, rho0Y = rho0Y, rho = rho_i )
    val <- results$effects.rr[term,1]
    p <- results$effects.rr[term,2]
    sens.out[i,1] <- rho_i
    sens.out[i,2] <- val
    sens.out[i,3] = p
  }
  sens.out <- data.frame(sens.out)
  return(sens.out)
}

run_sensitivity <- function(model.m = model.m, model.y = model.y,
                            treat = treat, mediator = mediator, ind_mediator = ind_mediator, outcome = outcome,
                            control.value = control.value, treat.value = treat.value,
                            rho01 = 0, rho0Y = 0, rho = 0){
#
#   if (!inherits(object$params, "params") )
#     stop("Input must be a 'mediation_result' object from mediation_analysis().")

    value = c(control.value, treat.value)
    datm <- model.m$data
    daty <- model.y$data

    if(rho01 != 0 & rho0Y !=0){
      Sigma <- matrix(c(1,    rho01,  rho0Y,
                        rho01,    1,  rho,
                        rho0Y, rho,   1  ), 3, 3)
      E <- mvrnorm(4000, mu = c(0,0,0), Sigma = Sigma)
      E0 <- E[,1]
      E1 <- E[,2]
      EY <- E[,3]
    }

    if(rho01 == 0 & rho0Y !=0){
      Sigma <- matrix(c(1,    0,   rho0Y,
                        0,    1,  rho,
                        rho0Y,rho,1), 3, 3)
      E <- mvrnorm(4000, mu = c(0,0,0), Sigma = Sigma)
      E0 <- E[,1]
      E1 <- E[,2]
      EY <- E[,3]
    }

    if(rho01 == 0 & rho0Y ==0){
      Sigma <- matrix(c(1,  rho,
                        rho,1), 2, 2)
      E <- mvrnorm(4000, mu = c(0,0), Sigma = Sigma)
      E0 = 0
      E1 <- E[,1]
      EY <- E[,2]
    }

    # Mediator ----
    ef_m = as_draws_df(model.m)
    dpar.m = c("mu", "hu")
    ## ZI part ----
    predict.ind_ms <- array(0, dim = c(2, ndraws(model.m), 1))
    dat.new = datm
    if(!is.null(ind_mediator)){
      predict.ind_ms <- array(NA, dim = c(2, ndraws(model.m), nrow(dat.new)))

      for(i in 1:length(value)){
        dat.new[, treat] = value[i]
        predict.ind_ms[i,,] = 0
        # --------------------------------------------------------------------------------------------------#
        predict.ind_ms[i,,] = plogis(1- posterior_linpred(model.m, newdata = dat.new, dpar = dpar.m[2]) + E0)
      }
    }else{
      predict.ind_ms <- array(0, dim = c(2, ndraws(model.m), 1))
    }

    ## Count Part ----
    predict.ms <- array(NA, dim = c(2, ndraws(model.m), nrow(dat.new)))

    for(i in 1:length(value)){
      dat.new[, treat] = value[i]
      predict.ms[i,,] = 0
      predict.ms[i,,] = exp(posterior_linpred(model.m, newdata = dat.new, dpar = NULL) + E1)
    }

    # Outcome ----
    daty = model.y$data
    ef_y = as_draws_df(model.y)

    predict.y.cov.mu <- array(NA, dim = c(2, ndraws(model.y), nrow(daty)) )
    for(i in 1:length(value)){
      dat.y.temp <- daty
      dat.y.temp[, treat] <- value[i]
      dat.y.temp[, mediator] = 0
      predict.y.cov.mu[i,,] = posterior_linpred(model.y, newdata = dat.y.temp)
    }

    coef.mediator = paste("b_", mediator, sep = "")
    bm = as.matrix(ef_y)[,coef.mediator]

    ## check indicator ----

    if( grepl("zero", family(model.y)$family) ){
      coef.mediator.zi = paste("b_zi_", mediator, sep = "")
      b_zi_m = as.matrix(ef_y)[,coef.mediator.zi]
    }
    if( grepl("hurdle", family(model.y)$family) ){
      coef.mediator.zi = paste("b_hu_", mediator, sep = "")
      b_zi_m = as.matrix(ef_y)[,coef.mediator.zi]
    }

    if(!is.null(ind_mediator)){
      coef.ind_mediator = paste("b_", ind_mediator, sep = "")
      b.ind_m = as.matrix(ef_y)[,coef.ind_mediator]

      if( grepl("zero", family(model.y)$family) ){
        coef.ind_mediator.zi = paste("b_zi_", ind_mediator, sep = "")
        b_zi_indm = as.matrix(ef_y)[,coef.ind_mediator.zi]
      }
      if( grepl("hurdle", family(model.y)$family) ){
        coef.ind_mediator.zi = paste("b_hu_", ind_mediator, sep = "")
        b_zi_indm = as.matrix(ef_y)[,coef.ind_mediator.zi]
      }
    }else{
      b.ind_m = 0
      b_zi_indm = 0
    }

    ## check interactions ----
    INT.xm  <- c(paste("b_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)),
                 paste("b_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )
    INT.xIm <- c(paste("b_", treat, ":", ind_mediator, sep="") %in% colnames(as.matrix(model.y)) ,
                 paste("b_", ind_mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )

    INT.xm_hu <- c(paste("b_hu_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)) ,
                   paste("b_hu_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )
    INT.xm_zi <- c(paste("b_zi_", treat, ":", mediator, sep="") %in% colnames(as.matrix(model.y)) ,
                   paste("b_zi_", mediator, ":", treat, sep="") %in% colnames(as.matrix(model.y))  )

    coef.intxm = c(paste("b_", treat, ":", mediator, sep=""), paste("b_", mediator, ":", treat, sep="") )
    coef.intxm = coef.intxm[INT.xm]

    coef.intxIm = c(paste("b_", treat, ":", ind_mediator, sep=""), paste("b_", ind_mediator, ":", treat, sep=""))
    coef.intxIm = coef.intxIm[INT.xIm]

    int_of_xm <- any(coef.intxm %in% colnames(as.matrix(model.y)))
    int_of_xIm <- any(coef.intxIm %in% colnames(as.matrix(model.y)))

    if(int_of_xm) {
      bxm = as.matrix(ef_y)[,coef.intxm]
    } else(
      bxm = 0
    )

    if(int_of_xIm ) {
      bxIm = as.matrix(ef_y)[,coef.intxIm]
    } else {
      bxIm = 0
    }

    outcome.linpred.mu = array(NA, dim = c(2, 2, 2, ndraws(model.y), nrow(daty)))

    ## Calculation function ----
    calc_linpred.mu <- function(i, j, r, bm, bxm, int_of_xm, int_of_xIm, b.ind_m ){
      if(int_of_xm){
        int_xm = 1
      } else {
        int_xm=0
        bxm=0
      }
      if(int_of_xIm){
        int_xIm = 1
      } else {
        int_xIm=0
        bxIm=0}

      if(is.null(ind_mediator)){
        as.numeric(
          as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*int_xm*value[i]*predict.ms[j,,]) )  +
          predict.y.cov.mu[i,,]
      }else{
        as.numeric(
          as.matrix(predict.ms[j,,]*bm) + as.matrix(bxm*int_xm*value[i]*predict.ms[j,,]) +
            as.matrix(b.ind_m*predict.ind_ms[r,,] + bxIm*int_xIm*value[i]*predict.ind_ms[r,,]) ) +
          predict.y.cov.mu[i,,]
      }
    }

    for(i in 1:length(value)){
      for(j in 1:2){
        for(r in 1:2){
          outcome.linpred.mu[i,j,r,,] =  calc_linpred.mu(i,j,r, bm, bxm, int_of_xm, int_of_xIm, b.ind_m ) + 3*EY
        }
      }
    }


    # Calculate results ----
    y_link = family(model.y)$link

    zi.outcome = grepl("zero", family(model.y)$family) | grepl("hurdle", family(model.y)$family)
    if(y_link == "identity" & !(zi.outcome))  outcome.pred = outcome.linpred.mu
    if(y_link == "logit"& !(zi.outcome))      outcome.pred  = exp(outcome.linpred.mu)/(1+exp(outcome.linpred.mu))

      if(y_link == "logit") {
        res.rd = cal.sens.rd(outcome.pred)
        res.rr = cal.sens.rr(outcome.pred)
        rst <- list(res.rd, res.rr)
      }else{
        res.rd = cal.sens.rd(outcome.pred)
      }
    rst <- list(effects.rd = res.rd, effects.rr = res.rr)
}

cal.sens.rr <- function(outcome.pred)
{
  indirect_control = outcome.pred[1,2,2,,] / outcome.pred[1,1,2,,]
  indirect_treated = outcome.pred[2,2,2,,] / outcome.pred[2,1,2,,]

  indirect = (indirect_control + indirect_treated)/2

  # **************************************************
  indirect_Im = outcome.pred[2,1,2,,] / outcome.pred[2,1,1,,]
  indirect = indirect * indirect_Im

  # **************************************************

  res = rbind(
    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<1), mean(indirect_control>1))),

    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<1), mean(indirect_treated>1))),

    # *******************************************************
    c(mean(indirect_Im), median(indirect_Im), sd(indirect_Im),
      quantile(indirect_Im, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im<1), mean(indirect_Im>1)))
    # *******************************************************
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=4)

  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  res = res[,c(2,4, 6)]
  res
}

cal.sens.rd <- function(outcome.pred)
{
  indirect_control = outcome.pred[1,2,2,,] - outcome.pred[1,1,2,,]
  indirect_treated = outcome.pred[2,2,2,,] - outcome.pred[2,1,2,,]

  indirect = (indirect_control + indirect_treated)/2

  # **************************************************
  indirect_Im = outcome.pred[2,1,2,,] - outcome.pred[2,1,1,,]
  indirect = indirect + indirect_Im

  # **************************************************

  res = rbind(
    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_control<1), mean(indirect_control>1))),

    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated),
      quantile(indirect_treated, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_treated<1), mean(indirect_treated>1))),

    # *******************************************************
    c(mean(indirect_Im), median(indirect_Im), sd(indirect_Im),
      quantile(indirect_Im, probs=c(0.025,0.975), na.rm = T),
      2*min(mean(indirect_Im<1), mean(indirect_Im>1)))
    # *******************************************************
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=4)

  rownames(res) = c("Indirect_control", "Indirect_treated", "Indirect_Indicator")
  colnames(res) = c("Mean", "Median", "SD", "l-95% CI", "u-95% CI", "Bayes_p")
  res = res[,c(2,6)]
  res
}











