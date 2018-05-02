
#Function to prepare data for PRIMEtime model.
# As default, does for deterministic analysis; but can also do PSA.
Manipulate_data <- function(psa = FALSE, singleCostMultiplier = FALSE, diab.sg = NULL){

  # Set up output list ----
  data.list <- list()

  # Single cost multiplier
  costMultiplier <- f_randomCosts(1)

  # NHS disease costs ----
  data.list$costs.hc <- if (psa){
    if (singleCostMultiplier){
      mutate(data.costs.hc, unitCost = unitCost * costMultiplier)
    } else {
      mutate_at(data.costs.hc, vars(unitCost), funs(f_randomCosts))
    }
  } else data.costs.hc

  # NHS non-disease costs ----
  data.list$costs.hc.other <- if (psa){
    if (singleCostMultiplier){
      mutate(data.costs.hc.other, cost = cost * costMultiplier)
    } else {
      mutate_at(data.costs.hc.other,
        vars(cost), funs(f_randomCosts))
    }
  } else data.costs.hc.other

  # Social care costs ----
  if (psa){
    if (singleCostMultiplier){
      data.list$costs.formalCare <- mutate_at(data.costs.formalCare,
        vars(-sex, -age), funs(. * costMultiplier))
    } else {
      data.list$costs.formalCare <- data.costs.formalCare
      data.list$costs.formalCare[, -c(1:2)] <- lapply(data.list$costs.formalCare[, -c(1:2)],
        function(x) f_randomCosts(x))
    }
  } else data.list$costs.formalCare <- data.costs.formalCare

  # BMI DATA ----
  
  data.list$bmiData <- if (diab.sg == "All adults"){
    data.bmi
  } else if (diab.sg == "Diabetics"){
    data.bmi.byDiab %>%
      filter(diabetes == TRUE) %>%
      select(-diabetes)
  } else {
    data.bmi.byDiab %>%
      filter(diabetes == FALSE) %>%
      select(-diabetes)
  }
  
  # RELATIVE RISK DATA ----

  # RRs for BMI on disease risks
  
  # copy data
  rrData <- data.rr

  # random draw from log-normal distribution
  if (psa) {
    log.sd <- sqrt(log(rrData$rr_se ^ 2 + exp(2 * log(rrData$rr_mean))) -
        2 * log(rrData$rr_mean))
    log.mean <- log(rrData$rr_mean) - .5 * log.sd^2
    rrData$rr.out <- rlnorm(length(log.mean), log.mean, log.sd)
  } else rrData$rr.out <- rrData$rr_mean

  #create frequency indicator by age
  rrData$freq <- round((rrData$ageEnd - rrData$ageStart) / 5)

  #expand rows
  rrData <- splitstackshape::expandRows(rrData, "freq")

  #create mean age within each row
  rrData <- rrData %>%
    group_by(disease, gender, ageStart) %>%
    mutate(seq = row_number() - 1,
      ageMean = (ageStart + 2) + 5 * seq)

  #create frequency indicator by sex
  rrData$freq <- ifelse(rrData$gender == "both", 2, 1)

  #expand rows
  rrData <- splitstackshape::expandRows(rrData, "freq")

  #row per male female
  rrData <- rrData %>%
    group_by(disease, gender, ageMean) %>%
    mutate(seq = row_number())
  rrData$sex = ifelse(rrData$gender == "both", ifelse(rrData$seq == 1, "male", "female"), rrData$gender)

  #group data
  rrData <- rrData %>%
    ungroup() %>%
    select(disease, ageMean, sex, rr.out) %>%
    mutate(ageGrp = as.character(cut(ageMean, breaks = seq(0, 100, 5), right = FALSE))) %>%
    filter(!is.na(ageGrp)) %>%
    select(-ageMean) %>%
    spread(disease, rr.out, fill = 1)

  #merge RR data into BMI data
  data.list$rrData <- right_join(rrData, data.list$bmiData, by = c("sex", "ageGrp"))

  # DISEASE NAMES ----
  data.list$disease.names <- if (diab.sg == "Diabetics") {
    data.disease.names[data.disease.names != "diabetes"]
  } else {
    data.disease.names
  }
  
  # BASELINE DISEASE RATES ----
  
  # Generate RRs for effect diabetes on disease; random draw from log-normal distribution
  rr.diab.inc <- data.diab.rr.incidence
  rr.diab.m <- data.diab.rr.mortality
  if (psa) {
      
      # incidence
      log.sd <- sqrt(log(data.diab.rr.incidence$rr_se ^ 2 + 
          exp(2 * log(data.diab.rr.incidence$rr_mean))) - 
          2 * log(data.diab.rr.incidence$rr_mean))
      log.mean <- log(data.diab.rr.incidence$rr_mean) - .5 * log.sd^2
      rr.diab.inc$rr.out <- rlnorm(length(log.mean), log.mean, log.sd)

      # mortality
      log.sd <- sqrt(log(data.diab.rr.mortality$rr_se ^ 2 + 
          exp(2 * log(data.diab.rr.mortality$rr_mean))) - 
          2 * log(data.diab.rr.mortality$rr_mean))
      log.mean <- log(data.diab.rr.mortality$rr_mean) - .5 * log.sd^2
      rr.diab.m$rr.out <- rlnorm(length(log.mean), log.mean, log.sd)
    } else {
      
    rr.diab.inc$rr.out <- rr.diab.inc$rr_mean
    rr.diab.m$rr.out <- rr.diab.m$rr_mean
  }
  rr.diab.data <- list("rr.diab" = rr.diab.inc, "rr.mortality" = rr.diab.m)
  
  # Calculate adjusted disease incidence, prevalence, and mortality rates
  diab.disease.data <- generate_diseaseRates_byDiabetes(relative.risks = rr.diab.data)
  
  # Select disease data according to diabetes status
  if (diab.sg == "All adults"){
    data.list$baselineIncidenceRates <- data.incidenceRates
    data.list$baselinePrevalence  <- data.prevalence
    data.list$totMortalityRates <- data.mortalityRates
  } else if (diab.sg == "Diabetics"){
    data.list$baselineIncidenceRates <- diab.disease.data$diab.incidenceRates %>%
      filter(has.diabetes == TRUE) %>%
      select(-has.diabetes)
    data.list$baselinePrevalence  <- diab.disease.data$diab.prevalence %>%
      filter(has.diabetes == TRUE) %>%
      select(-has.diabetes)
    data.list$totMortalityRates <- diab.disease.data$diab.mortalityRates %>%
      filter(has.diabetes == TRUE) %>%
      select(-has.diabetes)
  } else {
    data.list$baselineIncidenceRates <- diab.disease.data$diab.incidenceRates %>%
      filter(has.diabetes == FALSE) %>%
      select(-has.diabetes)
    data.list$baselinePrevalence  <- diab.disease.data$diab.prevalence %>%
      filter(has.diabetes == FALSE) %>%
      select(-has.diabetes)
    data.list$totMortalityRates <- diab.disease.data$diab.mortalityRates %>%
      filter(has.diabetes == FALSE) %>%
      select(-has.diabetes)
  }
  
  # unaffected by diabetes status
  data.list$baselineCaseFatality <- data.caseFatality
  data.list$trendsCaseFatality <- data.caseFatalityTrends
  data.list$trendsIncidence <- data.incidenceTrends
  
  # DALY WEIGHTS ----
  data.list$dalyWeights <- data.dalyWt
  
  # QUALITY OF LIFE DATA ----

  #mean utilities by age categories
  utilityAge <- if (psa){
    mutate(data.utilityAge,
      uAge = rnorm(n(), meanU, seU)) %>%
      select(ageBand, uAge, meanAge)
  } else {
    rename(data.utilityAge,
      uAge = meanU) %>%
      select(ageBand, uAge, meanAge)
  }

  #decrements per year of age (compared to mean in cat) and for men vs women
  uDec.age <- ifelse(psa, rnorm(1, -0.0002747, 0.000165), -0.0002747)
  uDec.male <- ifelse(psa, rnorm(1, 0.0010046, 0.0006241), 0.0010046)

  #utility decrements by incident & prevalent cases
  utilityDecDisease <- if (psa){
    data.utilityDecs %>%
      mutate(utilityDec.inc = rnorm(n(), utilityDec.inc.m, utilityDec.inc.se),
        utilityDec.prev = rnorm(n(), utilityDec.prev.m, utilityDec.prev.se)) %>%
      select(disease, utilityDec.inc, utilityDec.prev)
  } else {
    data.utilityDecs %>%
      rename_at(vars(contains(".m")),
        funs(str_replace(., ".m", ""))) %>%
      select(disease, utilityDec.inc, utilityDec.prev)
  }
  data.list$utilityDecDisease <- utilityDecDisease

  #data for baseline QoL calc, create age groups
  baselineQol <- data.pop %>%
    mutate(ageBand = as.character(cut(age, breaks = c(-Inf, seq(10, 80, 10), Inf),
      right = FALSE,
      labels = c(paste(seq(0, 70, 10), seq(9, 79, 10), sep = "-"), "80+"))))

  #add in mean QoL by 10-yr age
  baselineQol <- left_join(baselineQol, utilityAge, by = "ageBand")

  #calculate pre-disease utility (i.e. age and sex only)
  baselineQol <- baselineQol %>%
    mutate(utility = uAge + (age - meanAge) * uDec.age,
      utility = ifelse(sex == "male", utility + uDec.male, utility),
      utility = ifelse(utility > 1, 1, utility))

  # remove unneeded ages
  baselineQol <- baselineQol %>%
    filter(age >= min(data.list$baselinePrevalence$age))
  
  #calc decrements across diseases for prevalent & incident cases by age and sex
  baselineQol$prev.dec <- as.matrix(data.list$baselinePrevalence[, data.disease.names]) %*%
    utilityDecDisease$utilityDec.prev[utilityDecDisease$disease %in% data.disease.names]
  baselineQol$inc.dec <- as.matrix(data.list$baselineIncidenceRates[, data.disease.names]) %*%
    utilityDecDisease$utilityDec.inc[utilityDecDisease$disease %in% data.disease.names]

  #calculate final utility by age and sex
  baselineQol$utility <- rowSums(baselineQol[, c("utility", "prev.dec", "inc.dec")])

  #restrict to required output
  baselineQol <- baselineQol %>%
    select(sex, age, utility) %>%
    filter(age >= 20)

  #save output
  data.list$baselineQol <- baselineQol

  # Return data list ----
  data.list

}


# Function for creating new population distribution
# based on selected population characteristics
Define_targeted_population <- function(min.bmi = NULL, max.bmi = NULL, data.list = NULL){

  # Counts by sex and age (in 5-year bands)
  data.pop <- data.pop %>%
    filter(age >= min(data.list$baselineIncidenceRates$age)) %>%
    mutate(ageGrp = as.character(cut(age, seq(0, 100, 5), right = FALSE)),
      ageGrp = ifelse(age==100, "[95,100)", ageGrp)) %>%
    group_by(sex, ageGrp) %>%
    summarise(count = sum(count),
              ageMedian = floor(median(age)))

  # Proportion of population in defined BMI range
  data.pop <- left_join(data.pop, data.list$bmiData, by = c("sex", "ageGrp")) %>%
      mutate(ln.sd = sqrt(log(sd*sd + exp(2 * log(mean))) - 2 * log(mean)),
             ln.mean = log(mean) - .5 * ln.sd^2,
             propTarget = plnorm(max.bmi, ln.mean, ln.sd) - plnorm(min.bmi, ln.mean, ln.sd))

  # Clean up pop data to be returned
  data.pop <- select(data.pop, sex, ageMedian, ageGrp, count, propTarget)

  # Return population data
  data.pop

}


#function for converting class of empty disease to numeric
# and replacing NA's with zeros
f.convertClassEmptyDis <- function(data){
  outData <- data
  outData[sapply(outData, is.logical)] <- lapply(outData[sapply(outData, is.logical)], as.numeric)
  outData[sapply(outData, is.numeric)] <- lapply(outData[sapply(outData, is.numeric)],
    function(x) ifelse(is.na(x), 0, x))
  return(outData)
}

#Function for calculating random gamma draw from vector of costs
f_randomCosts <- function(input, formalCare = FALSE){
  m <- if (formalCare) -input else input
  s <- m * .1
  a <- (m / s) ^ 2
  b <- 1 / ((s ^ 2) / m)
  out <- ifelse(m == 0, 0, rgamma(length(m), a, b))
  out <- if (formalCare) -out else out
  return(out)
}


# Functioning for calculate disease rates by diabetes subgroup
generate_diseaseRates_byDiabetes <- function(relative.risks = NULL){
  
  # manipulate relative risk data ----
  
  # Add RR data to environment
  list2env(relative.risks, envir = environment())
  
  # expand relative risks for mortality by age
  rr.mortality$freq <- rr.mortality$end - rr.mortality$start + 1
  sum(rr.mortality$freq) == 101
  rr.mortality <- splitstackshape::expandRows(rr.mortality, "freq")
  rr.mortality$age <- 0:(nrow(rr.mortality)-1)
  rr.mortality <- select(rr.mortality, age, rr = rr.out)
  
  # select outcome diab on disease incidence
  rr.diab <- select(rr.diab, disease, sex, rr = rr.out)

  # estimate incidence and mortality rates by diabetes status ====
  
  # add mortality rate to incidence data
  temp.incidence <- left_join(data.incidenceRates, data.mortalityRates, by = c("sex", "age"))
  
  # create new incidence rates for each outcome
  new.incidence.data <- NULL
  for (d in c("ihd", "stroke", "mortalityRate")){
    
    # extract data and add diabetes prevalence 
    tempData <- select(temp.incidence, sex, age, rate = d) %>%
      filter(age >= 20) %>%
      left_join(., select(data.prevalence, sex, age, diabetes), by = c("sex", "age")) 
    
    # add in relative risk data  
    if (d != "mortalityRate"){ #FIX ERROR HERE,  "mortalityRates"
      tempData <- tempData %>%
        left_join(., filter(rr.diab, disease == d), by = "sex") #FIX ERROR HERE, disease = d
    } else {
      tempData <- tempData %>%
        left_join(., rr.mortality, by = "age")
    }
    
    # adjust data
    tempData$N <- 1000
    tempData$Nd <- tempData$N * tempData$diabetes
    tempData$Nh <- tempData$N - tempData$Nd
    tempData$E <- tempData$N * tempData$rate
    tempData$Rh <- tempData$E / (tempData$Nd * tempData$rr + tempData$Nh)
    tempData$Rd <- tempData$Rh * tempData$rr
    tempData$Ed <- tempData$Nd * tempData$Rd
    tempData$Eh <- tempData$Nh * tempData$Rh
    tempData$rate.T <- tempData$Ed / tempData$Nd
    tempData$rate.F <- tempData$Eh / tempData$Nh
    
    # tidy data
    tempData <- tempData %>%
      select(sex, age, rate.T, rate.F) %>%
      gather(key = has.diabetes, value = disease, rate.T, rate.F) %>%
      mutate(has.diabetes = has.diabetes == "rate.T") %>%
      rename_at("disease", funs(paste(d)))
    
    # add data to list
    if (is.null(new.incidence.data)){
      new.incidence.data <- tempData
    } else {
      new.incidence.data <- left_join(new.incidence.data, tempData, 
        by = c("sex", "age", "has.diabetes"))
    }
  }
  
  # extract mortality data
  new.mortality.data <- select(new.incidence.data, sex, age, has.diabetes, mortalityRate)
  
  #incorporate new incidence data into old, and adjust diabetes.  
  new.incidence.data <- new.incidence.data %>%
    left_join(., select(data.incidenceRates, -ihd, -stroke), by = c("sex", "age")) %>%
    mutate(diabetes = ifelse(has.diabetes, 0, diabetes)) %>%
    select(-mortalityRate)
  
  # Estimate disease prevalence ----
  
  # set-up output data
  new.prevalence.rates <- bind_rows(data.prevalence, data.prevalence) %>%
    mutate(has.diabetes = ifelse(row_number() <= n() / 2, TRUE, FALSE)) %>%
    filter(age >= 20)
  
  # Estimate for each disease
  for (d in c("ihd", "stroke")){
    
    # whether has diabetes
    for (diab in c(TRUE, FALSE)){
      
      # Combine data on prevalence, incidence, and case-fatality
      calcPrev <- select(new.incidence.data, sex, age, has.diabetes, i = d) %>%
        filter(has.diabetes == diab) %>%
        left_join(., select(data.caseFatality, sex, age, f = d),
          by = c("sex", "age")) %>%
        mutate(r = 0,
          I = i + r + f,
          q = sqrt(i^2 + 2*i*r - 2*i*f + r^2 + 2*f*r + f^2),
          w = exp(-(I+q)/2),
          v = exp(-(I-q)/2))
      
      # Calculate baseline prevalences in targeted population
      calcPrev[, c("S", "C", "D", "PY", "c", "b")] <- NA
      for (i in 1:nrow(calcPrev)){
        if (calcPrev$q[i]==0 | calcPrev$age[i]==20){
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
      new.prevalence.rates[new.prevalence.rates$has.diabetes == diab, d] <- calcPrev$c
    }
  }
  
  # changes to diabetes prevalence by status
  new.prevalence.rates <- mutate(new.prevalence.rates, 
    diabetes = ifelse(has.diabetes, 1, 0))
  
  # Data to extract ----
  diab.disease.data <- list(
    "diab.incidenceRates" = new.incidence.data,
    "diab.prevalence" = new.prevalence.rates,
    "diab.mortalityRates" = new.mortality.data
  )
  
  
}
