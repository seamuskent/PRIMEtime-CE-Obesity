# Data file for generating data for primetimeCE package
# last updated: 26.04.2018

#clear workspace
rm(list =ls())

# Global settings ----
options(stringsAsFactors = FALSE)

# LOAD DATA ----
Dir.path <- "J:/obesity_modelling/PRIME/PRIMEtime-CE-Obesity/data-raw/"

# Original PRIMEtime data
data.pop <- read.csv(paste0(Dir.path, "n_pop_by_age_1yr_31102017.csv"))
data.mortalityRates <- read.csv(paste0(Dir.path, "totMortalityRate_byAgeSex_31102017.csv"))
data.incidenceRates <- read.csv(paste0(Dir.path, "diseaseIncidenceRates_byAgeSex_31102017.csv"))
data.caseFatality <- read.csv(paste0(Dir.path, "caseFatalityRates_byAgeSex_31102017.csv"))
data.prevalence <- read.csv(paste0(Dir.path, "baselinePrevalence_byAgeSex_31102017.csv"))
data.incidenceTrends <- read.csv(paste0(Dir.path, "trends_incidence_byDiseaseAgeSex_01112017.csv"))
data.caseFatalityTrends <- read.csv(paste0(Dir.path, "trends_caseFatality_byDiseaseAgeSex_01112017.csv"))
data.dalyWt <- read.csv(paste0(Dir.path, "dalyWeights_byAgeSexDisease_31102017.csv"))
data.utilityDecs <- read.csv(paste0(Dir.path, "utilityDecrements_disease_01112017.csv"))
data.costs.hc <- read.csv(paste0(Dir.path, "costs_disease_01112017.csv"))
data.costs.hc.other <- read.csv(paste0(Dir.path, "costs_otherNonDisease_01112017.csv"))
data.costs.formalCare <- read.csv(paste0(Dir.path, "societalcosts_ageSexDisease_02112017.csv"))
data.rr <- read.csv(paste0(Dir.path, "relativeRisk_byAgeSex_06022018.csv"))
data.utilityAge <- read.csv(paste0(Dir.path, "utilities_by_ageCat_02032018.csv"))
data.disease.names <- c("ihd", "stroke", "diabetes", "cancerBreast", "cancerColorectum",
                   "cancerLiver", "cancerKidney", "cancerPancreas")

# BMI smoothed from HSE 2014
data.bmi <- read.table(paste0(Dir.path, "bmi_by_sex_and_age5yr_HSE2014_smoothed_24042018.txt"), header = TRUE)

# BMI by diabetes from HSE 2014
data.bmi.byDiab <- read.table(paste0(Dir.path, "bmi_by_sex__age5yr_and_diab_HSE2014_25042018.txt"), header = TRUE)

# Disease data by diabetes
#load(paste0(Dir.path, "disease_data_by_diabetes_25042018.Rd"))
#data.incidenceRates.byDiab <- diab.data$incidence.by.diabetes
#data.prevalence.byDiab <- diab.data$prevalence.by.diabetes
#data.mortalityRates.byDiab <- diab.data$mortality.by.diabetes

# DETERMINISTIC DATA ----

# Disease incidence, prev, & case-fatality
data.incidenceRates <- f.convertClassEmptyDis(data.incidenceRates)
data.caseFatality <- f.convertClassEmptyDis(data.caseFatality)
data.prevalence <- f.convertClassEmptyDis(data.prevalence)

# Incidence & case-fatality trends
data.incidenceTrends <- select(data.pop, age, sex) %>%
  mutate(ageBand = as.character(cut(age, breaks = c(0, 35, 65, Inf), right = FALSE))) %>%
  left_join(., data.incidenceTrends, by = c("ageBand", "sex"))
data.caseFatalityTrends <- select(data.pop, age, sex) %>%
  mutate(ageBand = as.character(cut(age, breaks = c(0, 35, 65, Inf), right = FALSE))) %>%
  left_join(., mutate(data.caseFatalityTrends, diabetes = ihd), by = c("ageBand", "sex"))

# Inflate costs ====

# cost multiplier
cm <- 302.3 / 290.5

# Disease-related healthcare costs
data.costs.hc$unitCost <- data.costs.hc$unitCost * cm

# NHS unrelated costs
data.costs.hc.other$cost <- data.costs.hc.other$cost * cm

# formal care costs
tempNames <- names(data.costs.formalCare)[-c(1:2)]
data.costs.formalCare <- data.costs.formalCare %>%
  mutate_at(tempNames, funs(. * cm))

# RR for on disease rates ====

data.diab.rr.incidence <- read.csv(paste0(Dir.path, "relativeRisk_diabetes_byAgeSex_22122017.csv"))
data.diab.rr.mortality <- read.csv(paste0(Dir.path, "relativeRisk_diabetes_mortality_byAgeSex_25042018.csv"))

# SAVE AS R OBJECTS ----
usethis::use_data_raw()
usethis::use_data(data.pop, data.bmi,
                  data.mortalityRates, data.incidenceRates,
                  data.caseFatality, data.prevalence, data.incidenceTrends,
                  data.caseFatalityTrends, data.dalyWt,
                  data.costs.hc, data.costs.hc.other, data.costs.formalCare,
                  data.rr, data.utilityAge, data.utilityDecs,
                  data.disease.names, data.bmi.byDiab,
                  data.incidenceRates.byDiab, data.prevalence.byDiab,
                  data.mortalityRates.byDiab,
                  data.diab.rr.incidence, data.diab.rr.mortality,
  overwrite = TRUE, internal = FALSE)
