
#Function to prepare data for PRIMEtime model.
# As default, does for deterministic analysis; but can also do PSA.
Manipulate_data <- function(psa = FALSE, singleCostMultiplier = FALSE){

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

  # RELATIVE RISK DATA ----

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
  data.list$rrData <- right_join(rrData, data.bmi, by = c("sex", "ageGrp"))

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

  # remove unnecesary ages
  baselineQol <- baselineQol %>%
    filter(age >= min(baselinePrevalence$age))
  
  #calc decrements across diseases for prevalent & incident cases by age and sex
  baselineQol$prev.dec <- as.matrix(baselinePrevalence[, disease.names]) %*%
    utilityDecDisease$utilityDec.prev[utilityDecDisease$disease %in% disease.names]
  baselineQol$inc.dec <- as.matrix(baselineIncidenceRates[, disease.names]) %*%
    utilityDecDisease$utilityDec.inc[utilityDecDisease$disease %in% disease.names]

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
Define_targeted_population <- function(min.bmi = NULL, max.bmi = NULL){

  # Counts by sex and age (in 5-year bands)
  data.pop <- data.pop %>%
    filter(age >= min(baselineIncidenceRates$age)) %>%
    mutate(ageGrp = as.character(cut(age, seq(0, 100, 5), right = FALSE)),
      ageGrp = ifelse(age==100, "[95,100)", ageGrp)) %>%
    group_by(sex, ageGrp) %>%
    summarise(count = sum(count),
              ageMedian = floor(median(age)))

  # Proportion of population in defined BMI range
  data.pop <- left_join(data.pop, data.bmi, by = c("sex", "ageGrp")) %>%
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

