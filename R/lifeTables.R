# Life tables and associated functions ----

# Standard lifetable ----
Produce_lifetable <- function(dataT = NULL,
                              timeH = NULL,
                              intervention = FALSE,
                              pifs = NULL,
                              targeted.pop.analysis = FALSE,
                              targeted.pop.data = NULL,
                              dis.lt.c = NULL,
                              dis.lt.t = NULL,
                              bmi.direct = NULL){

  # Set-up mortality data ----

  # For whole population
  lifeTable <- dataT

  # Adjust if analysis of targeted population by BMI
  if (targeted.pop.analysis){
    lifeTable <- select(dataT, -mortalityRate) %>%
      left_join(., targeted.pop.data[["targetedMortality"]], by = c("age", "sex"))
  }

  # Adjust mortality rates if intervention arm
  if (intervention){
    # If model mortality through diseases only
    if (!bmi.direct){
      for (d in disease.names){
        lifeTable$mortalityRate <- lifeTable$mortalityRate +
          (dis.lt.t[[d]]$bx - dis.lt.c[[d]]$bx)
      }
    } else {
      # If model mortality directly as a function of BMI
      lifeTable$pif <- NA
      for (i in 1:nrow(lifeTable)){
        if (lifeTable$age[i] == 100){
          lifeTable$pif[i] <- 1
        } else if (lifeTable$year[i] < timeH){
          pif.temp <- pifs[[lifeTable$year[i]+1]]
          lifeTable$pif[i] <- 1 - as.numeric(pif.temp[pif.temp$sex == lifeTable$sex[i] &
              pif.temp$ageGrp == lifeTable$ageGrp[i], "mortality"])
        } else{
          lifeTable$pif[i] <- 1
        }
      }
      lifeTable$mortalityRate <- lifeTable$mortalityRate * lifeTable$pif
    }
  }

  # Calculate life expectancy ----

  #generate new variables
  lifeTable[, c("Qx", "lx", "Dx", "Lx")] <- NA

  #probability of dying between age x and x+1
  lifeTable$Qx = ifelse(lifeTable$age < 100, 1 - exp(-lifeTable$mortalityRate), 1)

  #Number survivors at those at risk at start model &
  # number dying between ages x and x+1
  for (i in 1:nrow(lifeTable)){

    if (i == 1){ #starting number in cohort for men
      lifeTable$lx[i] <- 1
      lifeTable$Dx[i] <- lifeTable$lx[i] * lifeTable$Qx[i]
    } else if (i == nrow(lifeTable)/2 + 1){ #starting number in cohort for men
      lifeTable$lx[i] <- 1
      lifeTable$Dx[i] <- lifeTable$lx[i] * lifeTable$Qx[i]
    } else{ #the rest of the table
      lifeTable$lx[i] <- lifeTable$lx[i - 1] - lifeTable$Dx[i - 1]
      lifeTable$Dx[i] <- lifeTable$lx[i] * lifeTable$Qx[i]
    }
  }

  #Years lived by cohort to age x+.5
  for (i in 1:nrow(lifeTable)){
    lifeTable$Lx[i] <- ifelse(lifeTable$year[i] < timeH & lifeTable$age[i] < 100,
      (lifeTable$lx[i] + lifeTable$lx[i+1])/2, NA)
  }

  # Return lifetable
  lifeTable

}

# Disease lifetable ----
Produce_disease_lifetable <- function(dataT = NULL,
  disease = NULL,
  timeH = NULL,
  intervention = FALSE,
  pifs = NULL,
  second.trt = NULL,
  targeted.pop.analysis = FALSE,
  targeted.pop.data = NULL,
  trends = NULL,
  list.data = NULL){

  # load primetime data into environment
  list2env(list.data, envir = environment())
  
  # Select appropriate incidence & prevalence data
  if (targeted.pop.analysis){
    incidenceData <- targeted.pop.data[["targetedRates"]]
    prevData <- targeted.pop.data[["targetedPrev"]]
  } else {
    incidenceData <- baselineIncidenceRates
    prevData <- baselinePrevalence
  }

  # Add in incidence and case fatality rates ----
  lifeTable <- dataT %>%
    left_join(., select_(incidenceData, "age", "sex", ix = "disease"),
      by = c("age", "sex")) %>%
    left_join(., select_(baselineCaseFatality, "age", "sex", fx = "disease"),
      by = c("age", "sex"))

  # adjust for trends over 20years if selected ----
  if (trends){

    # estimate trend over time
    trends <- dataT %>%
      left_join(., select_(trendsIncidence, "age", "sex", incidence.data = disease),
        by = c("age", "sex")) %>%
      group_by(sex) %>%
      mutate(incidence.trend = ifelse(year <= 20, exp(incidence.data * year),
        exp(incidence.data * 20))) %>%
      left_join(., select_(trendsCaseFatality, "age", "sex", caseFatality.data = "disease"),
        by = c("age", "sex")) %>%
      group_by(sex) %>%
      mutate(caseFatality.trend = ifelse(year <= 20, exp(caseFatality.data * year),
        exp(caseFatality.data * 20))) %>%
      select(age, sex, incidence.trend, caseFatality.trend)

    # Adjust incidence and case-fatality rates
    lifeTable <- lifeTable %>%
      left_join(., trends, by = c("age", "sex")) %>%
      mutate(ix = ix * incidence.trend) %>%
      mutate(fx = fx * caseFatality.trend) %>%
      select(-incidence.trend, -caseFatality.trend)
  }

  # Adjust incidence rate in intervention arm ----
  if (intervention){
    lifeTable$pif <- NA
    for (i in 1:nrow(lifeTable)){
      if (lifeTable$age[i] == 100){
        lifeTable$pif[i] <- 1
      } else if (lifeTable$year[i] < timeH){
        pif.temp <- pifs[[lifeTable$year[i]+1]]
        lifeTable$pif[i] <- 1 - as.numeric(pif.temp[pif.temp$sex == lifeTable$sex[i] &
            pif.temp$ageGrp == lifeTable$ageGrp[i], disease])
      } else{
        lifeTable$pif[i] <- 1
      }
    }
    lifeTable$ix <- lifeTable$ix * lifeTable$pif
  }

  # Progression of hypothetical cohort ----

  # create columns for cohort number (healthy, diseased, dead)
  lifeTable[, c("Sx", "Cx", "Dx", "PYx", "cx", "bx")] <- NA

  # populate with data
  for (i in 1:nrow(lifeTable)){ #loop bc interdependent. ./2 cos gender.

    if (i == 1){ #starting number in cohort for men
      lifeTable$Cx[i] <- 10^3 *
        prevData[prevData$age == lifeTable$age[i] &
            prevData$sex == "male", disease]
      lifeTable$Sx[i] <- 10^3 - lifeTable$Cx[i]
      lifeTable$Dx[i] <- 0
    } else if (i == nrow(lifeTable)/2 + 1){ #starting number in cohort for women
      lifeTable$Cx[i] <- 10^3 *
        prevData[prevData$age == lifeTable$age[i] &
            prevData$sex == "female", disease]
      lifeTable$Sx[i] <- 10^3 - lifeTable$Cx[i]
      lifeTable$Dx[i] <- 0
    } else{ #the rest of the table
      lifeTable$Sx[i] <- lifeTable$Sx[i-1] * exp(-lifeTable$ix[i-1])
      lifeTable$Cx[i] <-
        lifeTable$Sx[i-1] * (1 - exp(-lifeTable$ix[i-1])) +
        lifeTable$Cx[i-1] * exp(-lifeTable$fx[i-1])
      lifeTable$Dx[i] <- lifeTable$Dx[i-1] +
        lifeTable$Cx[i-1] * (1- exp(-lifeTable$fx[i-1]))

    }
  }

  # num alive at start of period
  lifeTable <- mutate(lifeTable, num.alive = Sx + Cx)

  # Person-years lived between x and x+1, prevalence and mortality rates
  for (i in 1:nrow(lifeTable)){
    if (lifeTable$age[i] == 100 | lifeTable$year[i] >= timeH){
      lifeTable[i, c("PYx", "cx", "bx")] <- 0
    } else {

      #Years of life lived between x and x+1
      lifeTable$PYx[i] <- sum(lifeTable$num.alive[i:(i+1)])/2

      #Prevalence rate
      lifeTable$cx[i] <- (sum(lifeTable$Cx[i:(i+1)])/2) /
        lifeTable$PYx[i]

      #Mortality rate
      lifeTable$bx[i] <- (lifeTable$Dx[i+1] - lifeTable$Dx[i]) /
        lifeTable$PYx[i]

    }
  }

  # Clean up and return ----
  lifeTable <- select(lifeTable, sex,age, ageGrp, ix, cx, bx)
  lifeTable

}

# Estimate disease costs  ----
Estimate_disease_costs <- function(intervention = FALSE,
  disease = NULL,
  lifeTab = NULL,
  disease.lifeTabs = NULL,
  list.data = NULL){

  # Make data available in present environment ----
  list2env(list.data, envir = environment())

  # Select and combine lifetables ----
  disData <- disease.lifeTabs[[disease]] %>%
    left_join(., select(lifeTab, age, sex, Lx),
      by = c("sex", "age"))


  # Disease incidence ----
  disData[, paste0(disease, ".ix")] <- disData$ix * disData$Lx

  # Healthcare costs ----
  disData[,paste0("nhs.cost.", disease)] <- disData$cx * disData$Lx *
    costs.hc[which(costs.hc$disease == disease), "unitCost"]

  # Social care costs ----

  # NB - calc differs for ihd & stroke b/c diff costs inc & prev

  # Create new column name
  tempColName <- paste0("sc.cost.", disease)

  #ihd & stroke
  if (disease %in% c("ihd", "stroke")){

    tempData <- Estimate_socialCare_costs(costs.formalCare, c("nonDiseased", paste0(disease, ".i")),
      paste0(tempColName, ".i"))
    tempData <- Estimate_socialCare_costs(tempData, c("nonDiseased", paste0(disease, ".p")),
      paste0(tempColName, ".p"))
    tempData <- select_(tempData, "age", "sex", paste0(tempColName,".i"), paste0(tempColName,".p"))
    disData <- left_join(disData, tempData, by = c("age", "sex"))
    disData[tempColName] <- ((disData$ix * disData[paste0(tempColName, ".i")]) +
        (disData$cx - disData$ix) *
        disData[paste0(tempColName, ".p")] ) * disData$Lx
    disData[paste(tempColName, c("i", "p"), sep = ".")] <- NULL

    #other conditions
  } else {

    tempData <- Estimate_socialCare_costs(costs.formalCare,
      c("nonDiseased", disease), tempColName)
    tempData <- select(tempData, "age", "sex", tempColName)
    disData <- left_join(disData, tempData, by = c("age", "sex"))
    disData[tempColName] <- disData[tempColName] * disData$cx * disData$Lx

  }

  # Return output ----
  disData

}

#program for disease-specific social care costs
Estimate_socialCare_costs <- function(df, list_of_cols, new_col) {
  df %>%
    mutate_(.dots = ~Reduce(`-`, .[list_of_cols])) %>%
    setNames(c(names(df), new_col))
}

# Combine results from overall and disease specific lifetables to generate aggregate outputs per age group from model. Allows results per year up to time horizon or for a partic year to be produced.
Generate_outcomes <- function(lifeTab = NULL,
  disease.lifeTabs.c = NULL,
  disease.lifeTabs.t = NULL,
  list.data = NULL,
  intervention = FALSE,
  age.start = NULL,
  dr.health = NULL,
  dr.costs = NULL,
  timeH = NULL,
  cost.trt = 0,
  wt.loss.yr1 = NULL,
  yrs.to.regain = NULL,
  wt.mntnd = NULL,
  bmi.direct = TRUE,
  popData = NULL){

  # Make data available in present environment ----
  list2env(list.data, envir = environment())

  # Select disease life-tables based on intervention status ----
  disease.lifeTabs <- if (intervention) disease.lifeTabs.t else disease.lifeTabs.c


  # DALYS ----

  # merge in daly weight by age and sex for gen pop
  lifeTab <- lifeTab %>%
    left_join(., select(dalyWeights, age, sex, wx = genPop),
      by = c("age", "sex"))

  # adjust if intervention arm
  if (intervention){
    lifeTab <- left_join(lifeTab,
      select(dalyWeights, age, sex, disease.names),
      by = c("age", "sex"))
    for (d in disease.names){
      lifeTab$wx <- lifeTab$wx + pull(lifeTab, d) *
        (disease.lifeTabs.t[[d]]$cx - disease.lifeTabs.c[[d]]$cx)
    }
    lifeTab[, disease.names] <- NULL
  }

  # calculate DALYs
  lifeTab$Lwx <- lifeTab$Lx * (1 - lifeTab$wx)

  # QALYs ----

  # merge in utility weights by age and sex
  lifeTab <- left_join(lifeTab,
    select(baselineQol, age, sex, ux = utility),
    by = c("age", "sex"))

  # adjust if intervention arm
  if (intervention){
    for (d in disease.names){
      if (d %in% c("ihd", "stroke")){
        lifeTab$ux <- lifeTab$ux +
          utilityDecDisease[utilityDecDisease$disease == d, "utilityDec.inc"] *
          (disease.lifeTabs.t[[d]]$ix - disease.lifeTabs.c[[d]]$ix) +
          ((disease.lifeTabs.t[[d]]$cx - disease.lifeTabs.c[[d]]$cx) -
              (disease.lifeTabs.t[[d]]$ix - disease.lifeTabs.c[[d]]$ix)) *
          utilityDecDisease[utilityDecDisease$disease == d, "utilityDec.prev"]
      } else{
        lifeTab$ux <- lifeTab$ux +
          (disease.lifeTabs.t[[d]]$cx - disease.lifeTabs.c[[d]]$cx) *
          utilityDecDisease[utilityDecDisease$disease == d, "utilityDec.prev"]
      }
    }
  }

  # Direct effect of weight on EQ-5D...

  # Independent QALY gain per kg weight loss
  qaly.gain.per.kg.wt.loss <- 0.0028

  # Add in if intervention
  if (intervention & bmi.direct){
    # Weight loss per year
    if (is.null(wt.mntnd)){
      wt.per.yr <- wt.loss.yr1 - (wt.loss.yr1 / yrs.to.regain) * 0:(yrs.to.regain-1)
      wt.per.yr <- append(wt.per.yr, rep(0, length = timeH - yrs.to.regain + 1))
    } else {
      wt.per.yr <- wt.loss.yr1 - (wt.loss.yr1 / yrs.to.regain) * 0:(yrs.to.regain-1)
      wt.per.yr <- append(wt.per.yr, c(rep(wt.mntnd, length = timeH - yrs.to.regain), 0))
    }
    wt.per.yr <- wt.per.yr[1:(nrow(lifeTab)/2)]

    # add in proportion targeted
    tempD <- popData %>%
      filter(ageMedian == lifeTab$age[1]) %>%
      select(sex, propTarget)
    lifeTab <- left_join(lifeTab, tempD, by = "sex")

    # Adjust utility
    lifeTab <- lifeTab %>%
      mutate(ux = ux + propTarget * wt.per.yr * qaly.gain.per.kg.wt.loss) %>%
      select(-propTarget)
  }

  # problem: these reductions should only apply to those targeted not whole pop.

  #calculations need to be adjusted for older ages.

  # calculate QALYs
  lifeTab$Lux <- lifeTab$Lx * lifeTab$ux

  # Healthcare costs, unrelated ----
  lifeTab <- lifeTab %>%
    left_join(., costs.hc.other, by = c("age", "sex")) %>%
    mutate(nhs.cost.other = cost * Lx) %>%
    select(-cost)

  # Social care costs, unrelated ----
  lifeTab <- lifeTab %>%
    left_join(., select(costs.formalCare, age, sex, nonDiseased),
      by = c("age", "sex")) %>%
    mutate(sc.cost.other = -nonDiseased * Lx) %>%
    select(-nonDiseased)

  # Add in disease specific costs (and incidence rates) ----
  for (disease in disease.names){

    #disease costs
    tempData <- disease.lifeTabs[[disease]] %>%
      select(age, sex, contains(disease, vars = names(disease.lifeTabs[[disease]])))

    #merge in to lifetable
    lifeTab <- lifeTab %>%
      left_join(., tempData, by = c("age", "sex"))

  }

  #combined costs
  lifeTab$nhs.cost.total <- rowSums(lifeTab[, contains("nhs.cost", vars = names(lifeTab))])
  lifeTab$nhs.cost.cancerAll <- rowSums(lifeTab[, contains("nhs.cost.cancer", vars = names(lifeTab))])
  lifeTab$sc.cost.total <- rowSums(lifeTab[, contains("sc.cost", vars = names(lifeTab))])
  lifeTab$nhs.cost.disease <- lifeTab$nhs.cost.total - lifeTab$nhs.cost.other
  lifeTab$sc.cost.disease <- lifeTab$sc.cost.total - lifeTab$sc.cost.other

  #combined cancer incidence
  tempNames <- names(lifeTab)[contains("cancer", vars = names(lifeTab))]
  tempNames <- tempNames[contains(".ix", vars = tempNames)]
  lifeTab$cancer.ix <- rowSums(lifeTab[, tempNames])

  # Treatment costs ----
  if (!intervention){
    lifeTab$trt.cost <- 0
  } else {
    lifeTab <- lifeTab %>%
      group_by(sex) %>%
      mutate(trt.cost = ifelse(row_number() == 1, cost.trt, 0))
  }

  # Discounting health & costs ----

  #Health outcomes
  tempNames <- c("Lx", "Lwx", "Lux")
  lifeTab <- lifeTab %>%
    group_by(sex) %>%
    mutate_at(tempNames, funs(. / ((1 + dr.health/100)^(row_number()-1))))

  #Costs
  tempNames <- names(lifeTab)[contains("cost", vars = names(lifeTab))]
  lifeTab <- lifeTab %>%
    group_by(sex) %>%
    mutate_at(tempNames, funs(. / ((1 + dr.costs/100)^(row_number()-1))))

  # Calculate cumulative outcomes over time ----

  #define outcome variables
  outcomes <- c("Lx", "Lux", "Lwx",
                "nhs.cost.total", "nhs.cost.other","nhs.cost.disease",
                "sc.cost.total", "sc.cost.other", "sc.cost.disease",
                "trt.cost",
                names(lifeTab)[contains(".ix", vars = names(lifeTab))])

  #discounted cumulative outcomes over time
  lifeTab <- lifeTab %>%
    filter(age < 100) %>%
    select(sex, age, outcomes) %>%
    group_by(sex) %>%
    mutate(year = row_number(), age = age.start) %>%
    filter(year <= timeH) %>%
    mutate_at(outcomes, cumsum) %>%
    select(sex, age, year, outcomes)

  #add in year 0 results (for graphing)
  year0 <- lifeTab %>%
    group_by(sex) %>%
    filter(row_number() == 1) %>%
    mutate_at(outcomes, funs(replace(., is.numeric(.), 0))) %>%
    mutate(year = 0)

  #carry observations forward beyond age 100 to time-horizon.
  if (age.start + timeH > 100){ 
    year.high <- max(lifeTab$year)
    locf <- lifeTab %>%
      group_by(sex) %>%
      filter(row_number() == n()) %>%
      splitstackshape::expandRows(count = timeH - year.high, count.is.col = FALSE) %>%
      group_by(sex) %>%
      mutate(year = year + row_number())
    lifeTab <- bind_rows(lifeTab, locf)
  }
  lifeTab <- bind_rows(lifeTab, year0) %>%
    ungroup()

  # Return data ----
  lifeTab

}

# Summarise data function ----
Summarise_ouput <- function(popData = NULL,
                            outData = NULL){

  # Define model outcomes
  outcomes <- names(outData$control)[sapply(outData$control, is.numeric)][-c(1:2)]

  # Number of treatments modelled
  num.trts <- length(contains("trt", vars = names(outData)))

  # Compare modelled treatments to control group
  out.list <- outData
  for (n in 1:num.trts){
    out.list[[paste0("diff_", n, "c")]] <- outData$control
    out.list[[paste0("diff_", n, "c")]][, outcomes] <- outData[[paste0("trt", n)]][, outcomes] - outData$control[, outcomes]
  }

  # Add in population data
  for (n in 1:length(out.list)){
    out.list[[n]] <- out.list[[n]] %>%
      left_join(., popData, by = c("sex", "age" = "ageMedian"))
  }

  # Convert proportion targeted to 1 for absolute outcomes
  for (n in c("targeted", "control", paste0("trt", 1:num.trts))){
    out.list[[n]] <- mutate(out.list[[n]], propTarget = 1)
  }

  # Adjust differences for prop pop targeted, and clean up
  for (n in 1:length(out.list)){
    out.list[[n]] <- out.list[[n]] %>%
      mutate_at(outcomes[-which(outcomes == "trt.cost")], funs(. / propTarget)) %>%
      select(sex, age, ageGrp, year, outcomes)
  }

  # Return data
  out.list

  }


# ESTIMATING RATE FOR TARGETED BMI POPULATION ----
Calculate_diseaseData_targetedPop <- function(list.data = NULL,
                                              bmi.max = NULL,
                                              bmi.min = NULL,
                                              bmi.minRisk = NULL,
                                              age.min = NULL){

  # Make data available in present environment ----
  list2env(list.data, envir = environment())

  # BASELINE INCIDENCE RATES ----
  # Set up data to store rate multipliers
  rateMultipliers <- rrData %>%
    select(sex, ageGrp, disease.names, mortality) %>%
    mutate_at(c(disease.names, "mortality"), funs(replace(., values = 1)))

  # Identify lowest age in category
  rateMultipliers$tempAge <- as.numeric(str_sub(str_split(rateMultipliers$ageGrp, ",", simplify = TRUE)[, 1], 2, -1))
  
  # For each age-sex group
  for (g in 1:nrow(rateMultipliers)){

    # keep as 1 if age less than first estimation age
    if (rateMultipliers$tempAge[g] >= age.min){
    
      # BMI data
      df.est <- data.frame(bmi = seq(10, max(50, bmi.max+1), 1))
      df.est$bmi.cat <- cut(df.est$bmi,
        breaks = c(-Inf, min(25, bmi.min-1), bmi.min, bmi.max, Inf),
        labels = c("normal", "belowTarget", "targeted", "aboveTarget"),
        right = FALSE)
  
      # Highest bmi
      max.bmi <- tail(df.est$bmi, 1)
  
      # Log mean and standard deviation
      bmi.log.sd <- sqrt(log(rrData$sd[g] ^ 2 + exp(2 * log(rrData$mean[g]))) -
          2 * log(rrData$mean[g]))
      bmi.log.mean <- log(rrData$mean[g]) - .5 * bmi.log.sd^2
  
      # Cumulative density
      df.est$cumDens <- ifelse(df.est$bmi == max.bmi, 1,
        plnorm(df.est$bmi, bmi.log.mean, bmi.log.sd))
  
      # Density & product of density & average BMI
      df.est[, c('dens', 'C.B')] <- NA
      for (i in 1:(nrow(df.est)-1)){
        df.est$dens[i] <- df.est$cumDens[i+1] - df.est$cumDens[i]
        df.est$C.B[i] <- ((df.est$bmi[i] + df.est$bmi[i+1])/2) * df.est$dens[i]
      }
  
      # Create data frame by BMI cat
      df.ref <- data.frame(bmi.cat = levels(df.est$bmi.cat))
  
      # BMI cut-points
      df.ref$bmi.cut
      for (i in 1:nrow(df.ref)){
        df.ref$bmi.cut[i] <- max(df.est$bmi[df.est$bmi.cat == df.ref$bmi.cat[i]])
        df.ref$bmi.cut[i] <- if (i < nrow(df.ref)) df.ref$bmi.cut[i]+1 else df.ref$bmi.cut[i]
      }
  
      # Distribution pop by BMI cat
      df.ref$prop <- NA
      for (i in 1:nrow(df.ref)){
        df.ref$prop[i] <- df.est$cumDens[df.est$bmi == df.ref$bmi.cut[i]]
        df.ref$prop[i] <- if (i == 1) df.ref$prop[i] else df.ref$prop[i] - sum(df.ref$prop[1:(i-1)])
      }
  
      # Mean BMI within each category
      df.ref$bmi <- NA
      for (i in 1:nrow(df.ref)){
        df.ref$bmi[i] <- sum(df.est$C.B[df.est$bmi.cat == df.ref$bmi.cat[i]], na.rm = TRUE) / df.ref$prop[i]
      }
  
      # Number of people by BMI
      df.ref$n <- df.ref$prop * 1
  
      # For each disease
      for (d in c(disease.names, "mortality")){
  
        # Reset data
        df.ref[, c("rr", "events", "rate")] <- NA
  
        # Normalised relative risks by BMI
        df.ref$rr <- as.numeric(rrData[g, d]) ^ (abs(df.ref$bmi - bmi.minRisk) / 5)
        df.ref$rr <- df.ref$rr / df.ref$rr[1]
  
        # Number of events  & rate by BMI
        df.ref$events[1] <- 1 / (1 + (df.ref$n[-1] %*% df.ref$rr[-1])/
            df.ref$n[1])
        df.ref$rate[1] <- df.ref$events[1] / df.ref$n[1]
        for (i in (1:nrow(df.ref))[-1]){
          df.ref$rate[i] <- df.ref$rr[i] * df.ref$rate[1]
          df.ref$events[i] <- df.ref$rate[i] * df.ref$n[i]
        }
  
        # Calculate rate multiplier and add to data frames
        rateMultipliers[g, d] <- df.ref$rate[df.ref$bmi.cat == "targeted"]
  
      }
    }
  }

  # Combine data on baseline rates and rate multipliers
  targetedRates <- baselineIncidenceRates %>%
    mutate(ageGrp = as.character(cut(age, breaks = seq(0, 100, 5), right = FALSE)),
      ageGrp = ifelse(age == 100, "[95,100)", ageGrp)) %>%
    left_join(., rateMultipliers, by = c("sex", "ageGrp"))

  # Calculate rates in targeted pop for each disease
  for (d in disease.names){
    targetedRates[, d] <- targetedRates[, paste0(d, ".x")] * targetedRates[, paste0(d, ".y")]
  }

  # Select data to retain
  targetedRates <- select(targetedRates, sex, age, disease.names)

  # MORTALITY RATE ----

  targetedMortality <- totMortalityRates %>%
    mutate(ageGrp = as.character(cut(age, breaks = seq(0, 100, 5), right = FALSE)),
      ageGrp = ifelse(age == 100, "[95,100)", ageGrp)) %>%
    left_join(., rateMultipliers, by = c("sex", "ageGrp")) %>%
    mutate(mortalityRate = mortalityRate * mortality) %>%
    select(sex, age, mortalityRate)

  # PREVALENCE RATES ----

  # Data frame to store output
  targetedPrev <- baselinePrevalence

  # For each disease
  for (d in disease.names){

    # Combine data on prevalence, incidence, and case-fatality
    calcPrev <- select_(targetedRates, "sex", "age", i = d) %>%
      left_join(., select_(baselineCaseFatality, "sex", "age", f = d),
        by = c("sex", "age")) %>%
      left_join(., select_(baselinePrevalence, "sex", "age", d),
        by = c("sex", "age")) %>%
      mutate(r = 0,
        I = i + r + f,
        q = sqrt(i^2 + 2*i*r - 2*i*f + r^2 + 2*f*r + f^2),
        w = exp(-(I+q)/2),
        v = exp(-(I-q)/2))

    # Calculate baseline prevalences in targeted population
    calcPrev[, c("S", "C", "D", "PY", "c", "b")] <- NA
    for (i in 1:nrow(calcPrev)){
      if (calcPrev$q[i]==0 | calcPrev$age[i]==min(calcPrev$age)){
        calcPrev[i, "S"] <- 1; calcPrev[i, c("C", "D")] <- 0
      } else {
        calcPrev$S[i] <- (2*(calcPrev$v[i]-calcPrev$w[i])*(calcPrev$S[i-1]*(calcPrev$f[i]+calcPrev$r[i]) +
            calcPrev$C[i-1]*calcPrev$r[i]) + calcPrev$S[i-1]*(calcPrev$v[i]*(calcPrev$q[i]-
                calcPrev$I[i]) + calcPrev$w[i]*(calcPrev$q[i]+calcPrev$I[i]))) / (2*calcPrev$q[i])
        calcPrev$C[i] <- -((calcPrev$v[i]-calcPrev$w[i])*(2*((calcPrev$f[i]+calcPrev$r[i])*(calcPrev$S[i-1]+
            calcPrev$C[i-1])-calcPrev$I[i]*calcPrev$S[i-1]) - calcPrev$C[i-1]*calcPrev$I[i]) -
            calcPrev$C[i-1]*calcPrev$q[i]*(calcPrev$v[i]+calcPrev$w[i])) / (2*calcPrev$q[i])
        calcPrev$D[i] <- (((calcPrev$v[i]-calcPrev$w[i])*(2*calcPrev$f[i]*calcPrev$C[i-1] -
            calcPrev$I[i]*(calcPrev$S[i-1]+calcPrev$C[i-1]))) - (calcPrev$q[i]*
                (calcPrev$S[i-1]+calcPrev$C[i-1])*(calcPrev$v[i]+calcPrev$w[i])) +
            (2*calcPrev$q[i]*(calcPrev$S[i-1]+calcPrev$C[i-1]+calcPrev$D[i-1]))) / (2*calcPrev$q[i])

      }
    }
    for (i in 1:nrow(calcPrev)){
      if (calcPrev$age[i] != 100) {
        calcPrev$PY[i] <- .5 * (calcPrev$S[i]+calcPrev$C[i]+calcPrev$S[i+1]+calcPrev$C[i+1])
        calcPrev$c[i] <- .5 * ((calcPrev$C[i]+calcPrev$C[i+1]) / calcPrev$PY[i])
        calcPrev$b[i] <- (calcPrev$D[i+1] - calcPrev$D[i]) / (calcPrev$PY[i])
      }
    }

    # Replace output in targeted group baseline prevalence data
    targetedPrev[, d] <- calcPrev$c

  }

  # COMBINE DATA ----
  out.list <- list("targetedRates" = targetedRates,
                   "targetedMortality" = targetedMortality,
                   "targetedPrev" = targetedPrev)
  out.list

}


