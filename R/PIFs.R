# Functions to calculate potential impact fractions ----

# Potential impact fraction calculation; different methods possible
Calculate_single_PIF <- function(bmi.mean, bmi.sd,
                          height.mean,
                          wtChange.mean,
                          bmi.target.min, bmi.target.max,
                          rr.per5bmi,
                          bmi.minRisk,
                          pif.type = "relative risk shift"){

  # BASELINE ----

  #BMI data
  df.est <- data.frame(bmi = seq(10, 50, 1))
  df.est$bmi.cat <- cut(df.est$bmi, breaks = c(-Inf, 25, 30, 35, 40, Inf),
    labels = c("normal", "overwt", "obese I", "obese II", "obese III"),
    right = FALSE)

  #highest bmi
  bmi.max <- tail(df.est$bmi, 1)

  #log mean and standard deviation
  bmi.log.sd <- sqrt(log(bmi.sd ^ 2 + exp(2 * log(bmi.mean))) -
      2 * log(bmi.mean))
  bmi.log.mean <- log(bmi.mean) - .5 * bmi.log.sd^2

  #cumulative density
  df.est$cumDens <- ifelse(df.est$bmi == bmi.max, 1,
    plnorm(df.est$bmi, bmi.log.mean, bmi.log.sd))

  #Density & product of density & average BMI
  df.est[, c('dens', 'C.B')] <- NA
  for (i in 1:(nrow(df.est)-1)){
    df.est$dens[i] <- df.est$cumDens[i+1] - df.est$cumDens[i]
    df.est$C.B[i] <- ((df.est$bmi[i] + df.est$bmi[i+1])/2) * df.est$dens[i]
  }

  #mean weight within BMI cat
  df.est$wt <- ifelse(df.est$bmi == bmi.max, NA,
    (df.est$bmi + .5) * (height.mean ^ 2))

  # INTERVENTION ----

  #mean weight following intervention
  df.est$wt2 <- ifelse(df.est$bmi < bmi.target.min | df.est$bmi > bmi.target.max,
    df.est$wt, df.est$wt - wtChange.mean)

  #mean BMI following intervention
  df.est$bmi2 <- df.est$wt2/(height.mean^2) - .5
  df.est$bmi2[df.est$bmi == bmi.max] <- df.est$bmi2[df.est$bmi == bmi.max - 1] + 1

  #product of density & average BMI following interventio
  df.est$C.B2 <- NA
  for (i in 1:(nrow(df.est)-1)){
    if (df.est$bmi[i] < (bmi.target.min-1) | df.est$bmi[i] > bmi.target.max){
      df.est$C.B2[i] <- ((df.est$bmi[i] + df.est$bmi[i+1])/2) * df.est$dens[i]
    } else {
      df.est$C.B2[i] <- ((df.est$bmi2[i] + df.est$bmi2[i+1])/2) * df.est$dens[i]
    }
  }

  #density in weight cats following intervention
  for (cat in levels(df.est$bmi.cat)){
    temp.max <- max(df.est$bmi[df.est$bmi.cat == cat])
    df.est[, cat] <- ifelse(df.est$bmi2 < temp.max, df.est$dens,
      pmax(0, (temp.max+1 - df.est$bmi2)*df.est$dens))
  }

  # DATA SET UP FOR PIF CALC ----

  #create data frame by BMI cat
  df.ref <- data.frame(bmi.cat = levels(df.est$bmi.cat))

  #BMI cut-points
  df.ref$bmi.cut
  for (i in 1:nrow(df.ref)){
    df.ref$bmi.cut[i] <- max(df.est$bmi[df.est$bmi.cat == df.ref$bmi.cat[i]])
    df.ref$bmi.cut[i] <- if (i < nrow(df.ref)) df.ref$bmi.cut[i]+1 else df.ref$bmi.cut[i]
  }

  #distribution pop by BMI cat
  df.ref$prop <- NA
  for (i in 1:nrow(df.ref)){
    df.ref$prop[i] <- df.est$cumDens[df.est$bmi == df.ref$bmi.cut[i]]
    df.ref$prop[i] <- if (i == 1) df.ref$prop[i] else df.ref$prop[i] - sum(df.ref$prop[1:(i-1)])
  }

  #mean BMI within each category
  df.ref$bmi <- NA
  for (i in 1:nrow(df.ref)){
    df.ref$bmi[i] <- sum(df.est$C.B[df.est$bmi.cat == df.ref$bmi.cat[i]], na.rm = TRUE) / df.ref$prop[i]
  }

  #(normalised) relative risk for each BMI category
  df.ref$rr <- pmax(1, rr.per5bmi^((df.ref$bmi - bmi.minRisk)/5))
  df.ref$rr <- df.ref$rr / df.ref$rr[1]

  # CALCULATE PIF, SELECTED SCENARIO ----

  if (pif.type == "proportion shift"){

    #new distribution across categories
    df.ref$prop2 <- NA
    for (i in 1:(nrow(df.ref)-1)){
      df.ref$prop2[i] <- sum(df.est[, as.character(df.ref$bmi.cat[i])], na.rm = T)
      df.ref$prop2[i] <- if (i == 1) df.ref$prop2[i] else df.ref$prop2[i] - sum(df.ref$prop2[1:(i-1)])
    }
    df.ref$prop2[nrow(df.ref)] <- 1 - sum(df.ref$prop2, na.rm = T)

    #calculating PIF
    pif <- (df.ref$prop %*% df.ref$rr - df.ref$prop2 %*% df.ref$rr)/(df.ref$prop %*% df.ref$rr)

  } else if (pif.type == "relative risk shift"){

    #new mean bmi within cats
    df.ref$bmi2 <- tapply(df.est$C.B2, df.est$bmi.cat, sum, na.rm = T)/df.ref$prop

    #new relative risk (normalised)
    df.ref$rr2 <- pmax(1, rr.per5bmi^((df.ref$bmi2 - bmi.minRisk)/5))
    df.ref$rr2 <- df.ref$rr2 / df.ref$rr2[1]

    #calculating PIF
    pif <- (df.ref$prop %*% df.ref$rr - df.ref$prop %*% df.ref$rr2)/(df.ref$prop %*% df.ref$rr)

  } else{

    #simulated BMI data
    bmi.sim <- data.frame("bmi" = rlnorm(10000, bmi.log.mean, bmi.log.sd))

    #define relative risk function
    f.rr <- function(X, theta){theta ^ abs(X - bmi.minRisk)}

    #convert change in weight to change in BMI
    bmiChange.mean <- wtChange.mean / height.mean^2

    #define counterfactual distribution of BMI
    f.cft <- function(X){
      X2 <- X
      X2$bmi <- ifelse(X2$bmi < bmi.target.min | X2$bmi >= bmi.target.max,
        X2$bmi, X2$bmi - bmiChange.mean)
      return(X2$bmi)
    }

    #calculating PIF
    pif <- pif(X = bmi.sim, thetahat = rr.per5bmi, rr = f.rr,
      cft = f.cft, check_rr = FALSE)

  }

  # RETURN PIF ----
  pif <- as.numeric(pif)
  pif

}


# PIFs for each year of follow-up ----
calculate_PIFs_byYear <- function(wt.loss.mntnd = FALSE,
                                  mean.wt.loss.mantnd = NULL,
                                  timeH = NULL,
                                  list.data = NULL,
                                  trt.delay = NULL,
                                  yrs.to.bl.regain = NULL,
                                  wt.loss.yr1 = NULL,
                                  bmi.min = NULL,
                                  bmi.max = NULL,
                                  bmi.min.risk = NULL){

  # Make data available in present environment ----
  list2env(list.data, envir = environment())

  # Set up PIF list ----
  pif.list <- list()

  # Outcomes for PIF estimation - diseases & mortality
  outcomes <- append(disease.names, "mortality")

  # Calculations if no weight loss is maintained ----
  if (!wt.loss.mntnd){

    for (t in 1:timeH){ # for each year of modelling

      # set up diseases by age and sex
      pif.list[[t]] <- rrData

      # Calculate PIFs for each year
      if (t < trt.delay + 1 | t >= trt.delay + yrs.to.bl.regain){

        # If pre-trt effect or after weight regained, pif = 0
        pif.list[[t]][, outcomes] <- 0

      } else if (t == trt.delay + 1) {

        # If first year of treatment effect
        for (d in outcomes){ # for each disease
          for (r in 1:nrow(rrData)){ # for each age and sex group
            pif.list[[t]][r, d] <- Calculate_single_PIF(bmi.mean = rrData$mean[r],
                                                        bmi.sd = rrData$sd[r],
                                                        height.mean = rrData$height[r],
                                                        wtChange.mean = wt.loss.yr1,
                                                        bmi.target.min = bmi.min,
                                                        bmi.target.max = bmi.max,
                                                        rr.per5bmi = as.numeric(rrData[r, d]),
                                                        bmi.minRisk = bmi.min.risk,
                                                        pif.type = "relative risk shift")
            }
        }
      } else {

        # After first year of effect & until weight is regained
        pif.list[[t]] <- pif.list[[trt.delay + 1]]
        pif.list[[t]][, outcomes] <- pif.list[[t]][, outcomes] -
          (t - (trt.delay + 1)) * (pif.list[[t]][, outcomes] / (yrs.to.bl.regain-1))

      }
    }
  } else {

  # Calculations if some weight loss is maintained ----

    # Order years to allow for easier calculations in wt maintenance phase
    for (t in c(trt.delay+1, trt.delay + yrs.to.bl.regain,
      (1:timeH)[-c(trt.delay+1, trt.delay + yrs.to.bl.regain)])){

      # Set up diseases by age and sex
      pif.list[[t]] <- rrData

      # PIFs before treatment effect materialises
      if (t < trt.delay + 1){
        pif.list[[t]][, outcomes] <- 0

      # PIF for first year of treatment effect
      } else if (t == trt.delay + 1){
        for (d in outcomes){ # for each disease
          for (r in 1:nrow(rrData)){ # for each age and sex group
            pif.list[[t]][r, d] <- Calculate_single_PIF(bmi.mean = rrData$mean[r],
                                                        bmi.sd = rrData$sd[r],
                                                        height.mean = rrData$height[r],
                                                        wtChange.mean = wt.loss.yr1,
                                                        bmi.target.min = bmi.min,
                                                        bmi.target.max = bmi.max,
                                                        rr.per5bmi = as.numeric(rrData[r, d]),
                                                        bmi.minRisk = bmi.min.risk,
                                                        pif.type = "relative risk shift")
          }
        }
    # PIF for period when weight maintenance starts
      } else if (t == yrs.to.bl.regain + trt.delay) {
        for (d in outcomes){ # for each disease
          for (r in 1:nrow(rrData)){ # for each age and sex group
            pif.list[[t]][r, d] <- Calculate_single_PIF(bmi.mean = rrData$mean[r],
                                                        bmi.sd = rrData$sd[r],
                                                        height.mean = rrData$height[r],
                                                        wtChange.mean = mean.wt.loss.mantnd,
                                                        bmi.target.min = bmi.min,
                                                        bmi.target.max = bmi.max,
                                                        rr.per5bmi = as.numeric(rrData[r, d]),
                                                        bmi.minRisk = bmi.min.risk,
                                                        pif.type = "relative risk shift")
          }
        }
    # PIFs between start trt effect and start of weight mainteance
      } else if (t > trt.delay + 1 & t < yrs.to.bl.regain + trt.delay){

        pif.list[[t]] <- pif.list[[trt.delay + 1]]
        pif.list[[t]][, outcomes] <- pif.list[[t]][, outcomes] -
          (t - (trt.delay + 1)) * ((pif.list[[t]][, outcomes] -
              pif.list[[yrs.to.bl.regain + trt.delay]][, outcomes]) / (yrs.to.bl.regain-1))

    # PIFs in weight maintenance period
      } else {
          pif.list[[t]] <- pif.list[[yrs.to.bl.regain + trt.delay]]
    }
  }
}

  #Return PIF list
  pif.list

}

