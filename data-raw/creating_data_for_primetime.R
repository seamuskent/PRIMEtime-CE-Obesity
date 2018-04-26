# Data file for generating data for primetimeCE package
# last updated: 26.04.2018

#clear workspace
rm(list =ls())

# Global settings ----
options(stringsAsFactors = FALSE)

# LOAD DATA ----
Dir.path <- "J:/obesity_modelling/PRIME/primetime_pkg/primetimeCE/data-raw/"
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
data.bmi <- read.table("J:/obesity_modelling/PRIME/input_data/bmi_by_sex_and_age5yr_HSE2014_smoothed_24042018.txt", header = TRUE)

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


# SAVE AS R OBJECTS ----
usethis::use_data_raw()
usethis::use_data(data.pop, data.bmi,
                  data.mortalityRates, data.incidenceRates,
                  data.caseFatality, data.prevalence, data.incidenceTrends,
                  data.caseFatalityTrends, data.dalyWt,
                  data.costs.hc, data.costs.hc.other, data.costs.formalCare,
                  data.rr, data.utilityAge, data.utilityDecs,
                  data.disease.names,
  overwrite = TRUE, internal = FALSE)
