#' Estimation of the PRIMEtime-CE Obesity model
#'
#' \code{primetime} generates estimates of disease incidence, life-years, quality-adjusted life-years, and health and social care costs over the specified time-horizon up to age 100 years by sex and age group (in 5-year bands).
#'
#'
#' @param time.horizon A positive integer. Outcomes are estimated up to age 100 years.
#' @param age.range A two element numeric vector giving the lowest and highest ages over which to estimate outcomes.
#' @param bmi.target.min A positive integer giving the lowest BMI group at whom the intervention(s) is targeted.
#' @param bmi.target.max A positive integer giving the highest BMI group at whom the intervention(s) is targeted. A value of \emph{Inf} means no upper bound.
#' @param cost.of.treatment A numeric vector with length equal to the number of active interventions giving intervention costs in year 1. 
#' @param mean.wt.loss.yr1 A numeric vector with length equal to the number of active interventions to assess. Values are mean weight loss in kilograms 1-year following intervention.
#' @param se.diff A non-negative integer giving the standard error of the difference in effect between two treatment options. Only required if \emph{deterministic = FALSE}.
#' @param time.to.trt.effect A non-negative integer giving the delay between initiation of treatment and health effects.
#' @param time.to.weight.regain A positive integer giving the years till weight returns to baseline, or in the case of some long-term residual weight difference, the time till the weight difference stabilises.
#' @param some.wt.loss.maintained A logical value indicating whether there is some residual weight difference between treatment and the control groups after weight regain.
#' @param mean.wt.loss.maintained A numeric vector with length equal to the number of active interventions to assess. Values are mean weight loss in kilograms. Only required in \emph{some.wt.loss.maintained = TRUE}.
#' @param discount.rate.health,discount.rate.cost Non-negative values giving the annual discount rate for health benefits and costs. A value of 3.5 means a discount rate of 3.5\%.
#' @param deterministic Whether to run a deterministic (default) or probabilistic analysis.
#' @param whole.uk.population Whether to model outcomes over the whole UK population of selected age (default) or a user-specified population.
#' @param population.characeristics A five-element numeric vector reporting in order: the proportion of men, mean age in men, standard deviation of age in men, mean age in women, and standard deviation of age in women. Only required when \emph{whole.uk.population = TRUE}. 
#' @param subgroup Indicates whether to model outcomes across all adults ("All adults"; default) or in subgroup of individuals with type-2 diabetes ("Diabetics"") or without ("Non-diabetics").  
#' @return Mean per-person health and economic outcomes in groups defined by sex and age (in 5-year bands) for each year of analysis, for the control and treatment groups, and differences between treatment and control groups. The output can be accessed by other functions to produce standardised outputs (e.g. summary tables), or can be used directly by the user for bespoke analyses.
#' @examples
#' primetime(time.horizon = 25,age.range = c(20, 80),cost.of.treatment = c(34.06, 1001.58),bmi.target.min = 30, bmi.target.max = 60,mean.wt.loss.yr1 = c(3.1, 10.3),time.to.trt.effect = 1,time.to.weight.regain = 5,some.wt.loss.maintained = FALSE,mean.wt.loss.maintained = NULL,discount.rate.health = 3.5,discount.rate.cost = 3.5,whole.uk.population = TRUE)
#' @export
primetime <- function(time.horizon = NULL,
  age.range = c(20, 100),
                      cost.of.treatment = NULL,
                      bmi.target.min = NULL,
                      bmi.target.max = NULL,
                      mean.wt.loss.yr1 = NULL,
                      se.diff = NULL,
                      time.to.trt.effect = NULL,
                      time.to.weight.regain = NULL,
                      some.wt.loss.maintained = FALSE,
                      mean.wt.loss.maintained = NULL,
                      discount.rate.health = NULL,
                      discount.rate.cost = NULL,
                      deterministic = TRUE,
                      bmi.directEffect.onMortality = TRUE,
                      bmi.directEffect.onQoL = FALSE,
                      model.incidence.trends = TRUE,
                      whole.uk.population = TRUE,
                      population.characteristics = NULL,
                      subgroup = "All adults") {

  # Extract arguments - called and default ----
  arguments <- mget(names(formals()),sys.frame(sys.nframe()))

  # ERROR MANAGEMENT ----

  # All numeric variables are numeric
  for (v in c("time.horizon", "age.range", "bmi.target.min", "bmi.target.max",
    "mean.wt.loss.yr1", "time.to.trt.effect", "time.to.weight.regain",
    "discount.rate.health", "discount.rate.cost")){
    if (!is.numeric(get(v))){
      stop(v, " must be numeric.", call. = FALSE)
    }
  }


  # Time horizon
  if (time.horizon < 0 | time.horizon %% 1 != 0){
    stop("time.horizon must be a positive integer.", call. = FALSE)
  }

  # Acceptable age ranges
  if (any(!all(age.range %% 5 == 0), age.range[1] < 20, age.range[2] > 100,
    age.range[2] < age.range[1])){
    stop("Both ages must be multiples of 5, and be between 20 and 100 inclusive.", call. = FALSE)
  }

  # Non-contradictory number of interventions
  if (length(mean.wt.loss.yr1) != length(mean.wt.loss.maintained) & some.wt.loss.maintained){
    stop("mean.wt.loss.yr1 and mean.wt.loss.maintained must have the same length,
      equal to the number of interventions to be modelled.", call. = FALSE)
  }

  # Minimum BMI for targetting interventions
  if (bmi.target.min < 25 | !is.numeric(bmi.target.min)){
    stop("bmi.target.bmi must be 25 kg/m2 or greater.", call. = FALSE)
  }

  if (bmi.target.max > 57) bmi.target.max <- 57

  # Consistency of max and min bmi targets
  if (bmi.target.max <= bmi.target.min){
    stop("bmi.target.max must be more than bmi.target.min", call. = FALSE)
  }

  # Time to trt effect and weight regain must be positive integers
  if (time.to.trt.effect < 0 | time.to.trt.effect %% 1 != 0){
    stop("time.to.trt.effect must be a positive integer", call. = FALSE)
  }

  # must specify a population
  if (!whole.uk.population & is.null(population.characteristics)){
    stop("You must specify demographic characteristics of the population.")
  }
  
  # non-negative discount rates
  if (discount.rate.health < 0 | discount.rate.cost < 0){
    stop("Discount rates must be non-negative.")
  }
  
  # population characteristics
  if (!whole.uk.population & length(population.characteristics) != 5){
    stop("population.characteristics must be a numeric vector with length 5, containing 
the following information in order: proportion of men, mean age in men, standard deviation 
of age in men, mean age in women, standard deviation of age in women.")
  }
 
  # Standard errror of treatment effect
  if (!deterministic & (is.null(se.diff) | se.diff < 0)){
    stop("Please specifiy a non-negative numeric value for SE of the treatment effect. if you do not wish to model uncertainty in the treatment effect, please specify as 0.")
  }
  
  # DATA MANIPULATION ----

  # Data that is potentially probabilistic
  primetime.data <- Manipulate_data(psa = !deterministic, diab.sg = subgroup)

  # Generate proportion individuals targeted by age and sex group
  pop.data <- Define_targeted_population(min.bmi = bmi.target.min,
    max.bmi = bmi.target.max, data.list = primetime.data)

  # Minimum risk BMI
  bmi.min.risk <- 21

  # Repeat mean weight loss by arm
  bc.wt.loss <- mean.wt.loss.yr1

  # Random treatment effect
  if (!deterministic){
    if (length(mean.wt.loss.yr1) == 1){
      mean.wt.loss.yr1 <- rnorm(1, mean.wt.loss.yr1, se.diff)
    } else {
    mean.wt.loss.yr1[2] <- mean.wt.loss.yr1[1] +
      rnorm(1, mean.wt.loss.yr1[2] - mean.wt.loss.yr1[1], se.diff)
    }
  }
  
  # POTENTIAL IMPACT FRACTIONS ----

  pif.list <- list()
  for (n in 1:length(mean.wt.loss.yr1)){
    pif.list[[n]] <- calculate_PIFs_byYear(wt.loss.mntnd = some.wt.loss.maintained,
      mean.wt.loss.mantnd = mean.wt.loss.maintained[n],
      timeH = time.horizon,
      list.data = primetime.data,
      trt.delay = time.to.trt.effect,
      yrs.to.bl.regain = time.to.weight.regain,
      wt.loss.yr1 = mean.wt.loss.yr1[n],
      bmi.min = bmi.target.min,
      bmi.max = bmi.target.max,
      bmi.min.risk = bmi.min.risk,
      min.age = age.range[1])
  }

  # ESTIMATION BY AGE GROUP ----

  # Register parallelisation
  no_cores <- parallel::detectCores() - 1
  doParallel::registerDoParallel(no_cores)

  # Start analysis
  out.list <- foreach::foreach(age = seq(age.range[1], age.range[2] - 5, 5),
    .packages = "primetimeCE") %dopar% {

      # Initialise store of output for each type
      temp.list <- list()

      #Index age for estimation
      age.index <- age + 2

      # Create data template for lifetables for age group
      dataTemplate <- data.pop %>%
        select(age, sex) %>%
        left_join(., primetime.data$totMortalityRates, by = c("age", "sex")) %>%
        filter(age >= age.index) %>%
        mutate(ageGrp = as.character(cut(age, breaks = seq(20, 100, 5), right = FALSE))) %>%
        group_by(sex) %>%
        mutate(year = row_number() - 1) %>%
        filter(year <= min(time.horizon, 100 - age.index))

      # TARGETED POPULATION STATUS QUO ESTIMATION ----

      # Calculate disease data for targeted population
      diseaseData.targetedPop <- Calculate_diseaseData_targetedPop(list.data = primetime.data,
        bmi.max = bmi.target.max,
        bmi.min = bmi.target.min,
        bmi.minRisk = bmi.min.risk,
        age.min = age.range[1])

      # Standard life-table
      lifeTable <- Produce_lifetable(dataT = dataTemplate,
        timeH = time.horizon,
        targeted.pop.analysis = TRUE,
        targeted.pop.data = diseaseData.targetedPop)

      # Disease life-tables
      disease.lifetables.targeted <- list()
      for (d in primetime.data$disease.names){
        disease.lifetables.targeted[[d]] <- Produce_disease_lifetable(dataT = dataTemplate,
          disease = d,
          timeH = time.horizon,
          targeted.pop.analysis = TRUE,
          targeted.pop.data = diseaseData.targetedPop,
          trends = model.incidence.trends,
          list.data = primetime.data)
      }

      # Estimate costs for each disease
      for (d in primetime.data$disease.names){
        disease.lifetables.targeted[[d]] <- Estimate_disease_costs(disease = d,
          lifeTab = lifeTable,
          disease.lifeTabs = disease.lifetables.targeted,
          list.data = primetime.data)
      }

      # Generate output
      out.data <- Generate_outcomes(lifeTab = lifeTable,
        disease.lifeTabs.c = disease.lifetables.targeted,
        disease.lifeTabs.t = NULL,
        list.data = primetime.data,
        intervention = FALSE,
        age.start = age.index,
        dr.health = discount.rate.health,
        dr.costs = discount.rate.cost,
        timeH = time.horizon)

      # Combine output with that from other age groups
      temp.list$targeted <- out.data

      # WHOLE POPULATION STATUS QUO ESTIMATION ----

      # Standard life-table
      lifeTable <- Produce_lifetable(dataT = dataTemplate, timeH = time.horizon)

      # Disease life-tables
      disease.lifetables.control <- list()
      for (d in primetime.data$disease.names){
        disease.lifetables.control[[d]] <- Produce_disease_lifetable(dataT = dataTemplate,
          disease = d,
          timeH = time.horizon,
          trends = model.incidence.trends,
          list.data = primetime.data)
      }

      # Estimate costs for each disease
      for (d in primetime.data$disease.names){
        disease.lifetables.control[[d]] <- Estimate_disease_costs(disease = d,
          lifeTab = lifeTable,
          disease.lifeTabs = disease.lifetables.control,
          list.data = primetime.data)
      }

      # Generate output
      out.data <- Generate_outcomes(lifeTab = lifeTable,
        disease.lifeTabs.c = disease.lifetables.control,
        disease.lifeTabs.t = NULL,
        list.data = primetime.data,
        intervention = FALSE,
        age.start = age.index,
        dr.health = discount.rate.health,
        dr.costs = discount.rate.cost,
        timeH = time.horizon)

      # Combine output with that from other age groups
      temp.list$control <- out.data

      # TREATMENT GROUP ESTIMATION ----

      for (n in 1:length(mean.wt.loss.yr1)){ #for however many interventions there are

        # Disease life-tables
        disease.lifetables.trt <- list()
        for (d in primetime.data$disease.names){
          disease.lifetables.trt[[d]] <- Produce_disease_lifetable(dataT = dataTemplate,
            disease = d,
            intervention = TRUE,
            pifs = pif.list[[n]],
            timeH = time.horizon,
            trends = model.incidence.trends,
            list.data = primetime.data)
        }

        # Standard life-table
        lifeTable2 <- Produce_lifetable(dataT = dataTemplate,
          timeH = time.horizon,
          intervention = TRUE,
          pifs = pif.list[[n]],
          dis.lt.c = disease.lifetables.control,
          dis.lt.t = disease.lifetables.trt,
          bmi.direct = bmi.directEffect.onMortality)

        # Estimate costs for each disease
        for (d in primetime.data$disease.names){
          disease.lifetables.trt[[d]] <- Estimate_disease_costs(disease = d,
            lifeTab = lifeTable2,
            disease.lifeTabs = disease.lifetables.trt,
            list.data = primetime.data)
        }

        # Generate output
        out.data <- Generate_outcomes(lifeTab = lifeTable2,
          disease.lifeTabs.c = disease.lifetables.control,
          disease.lifeTabs.t = disease.lifetables.trt,
          list.data = primetime.data,
          intervention = TRUE,
          age.start = age.index,
          dr.health = discount.rate.health,
          dr.costs = discount.rate.cost,
          timeH = time.horizon,
          cost.trt = cost.of.treatment[n],
          wt.loss.yr1 = bc.wt.loss[n],
          yrs.to.regain = time.to.weight.regain,
          wt.mntnd = if (some.wt.loss.maintained) mean.wt.loss.maintained[n] else NULL,
          bmi.direct = bmi.directEffect.onQoL,
          popData = pop.data)

        # store data
        temp.list[[paste0("trt", n)]] <- out.data

      }

      # output to save per age group
      temp.list

    }

  # End parallelisation
  doParallel::stopImplicitCluster()

  # combine into list by intervention
  output.list <- list()
  for (i in 1:length(out.list)){
    if (i == 1){
      for (n in names(out.list[[1]])){
        output.list[[n]] <- out.list[[i]][[n]]
      }
    } else {
      for (n in names(out.list[[1]])){
        output.list[[n]] <- bind_rows(output.list[[n]], out.list[[i]][[n]])
      }
    }
  }

  # SUMMARISE RESULTS ----
  output.list <- Summarise_ouput(popData = pop.data,
    outData = output.list)

  # RETURN DATA ----
  output.list <- append(list("arguments" = arguments), output.list)
  output.list

  }


