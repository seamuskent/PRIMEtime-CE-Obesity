# Data summaries ----

# Create overall summary table ----
summary_table <- function(
  model.results = NULL,
  comparator = "control",
  active.intervention = "trt1",
  produce.results.table = TRUE,
  extended.summary = FALSE,
  timeH = NULL,
  costs.to.include = "disease-related nhs costs",
  present.only.overall.results = FALSE){ #other options: "unrelated nhs costs", "disease-related nhs and social care", "all nhs and social care"

  # Check that interventions to assess exist
  if (!exists(comparator, where = model.results)){
    stop("There is no comparator called ", comparator)
  }
  if (!exists(active.intervention, where = model.results)){
    stop("There is no intervention called ", active.intervention)
  }

  # Number in cohort
  number.in.cohort <- 1000

  # Defining difference of interest ----
  # If comparator is with 'no intervention'
  if (comparator == "control"){
    comparison <- c(paste0("diff_", str_sub(active.intervention, -1), "c"))
    # If comparing two active treatments
  } else {
    comparison <- c(paste0("diff_", str_sub(active.intervention, -1), str_sub(comparator, -1)))
    model.results[[comparison]] <- model.results[["control"]]
  }


  # Extract PRIMEtime arguments ----

  # Extract all arguments
  arguments <- model.results$arguments

  # Define population to apply model to
  whole.uk.population <- arguments$whole.uk.population
  population.characteristics <- arguments$population.characteristics

  # Time horizon for presentation of results
  timeH <- if(is.null(timeH)) arguments$time.horizon else timeH
  if (timeH > arguments$time.horizon){
    stop("timeH must be no greater than the time horizon specified for model estimation.")
  }

  # Define outcome variables ----
  outcomes <- names(model.results$control)[sapply(model.results$control, is.numeric)][-c(1:2)]

  # Create differences between arms, and absolute values by arm ----

  # Analysis versus no intervention...

  if (comparator == "control"){

    # Set targeted results as control
    model.results[["control"]] <- model.results[["targeted"]]

    # create absolute outcomes for treatment group
    model.results[[active.intervention]] <- model.results[["control"]]
    model.results[[active.intervention]][, outcomes] <- model.results[["control"]][, outcomes] +
      model.results[[comparison]][, outcomes]

    # list components to keep
    toKeep <- c(comparator, active.intervention, comparison)
    model.results <- model.results[names(model.results) %in% toKeep]

  } else {

    # Analysis comparing two interventions...

    # Difference between two treatments
    model.results[[comparison]] <- model.results[["targeted"]]
    model.results[[comparison]][, outcomes] <-
      model.results[[paste0("diff_", str_sub(active.intervention, -1), "c")]][, outcomes] -
      model.results[[paste0("diff_", str_sub(comparator, -1), "c")]][, outcomes]

    # create absolute outcomes for comparator
    model.results[[comparator]] <- model.results[["targeted"]]
    model.results[[comparator]][, outcomes] <- model.results[["targeted"]][, outcomes] +
      model.results[[paste0("diff_", str_sub(comparator, -1), "c")]][, outcomes]

    # create absolute outcomes for main treatment
    model.results[[active.intervention]] <- model.results[["targeted"]]
    model.results[[active.intervention]][, outcomes] <- model.results[["targeted"]][, outcomes] +
      model.results[[paste0("diff_", str_sub(active.intervention, -1), "c")]][, outcomes]

    # remove irrelevant components of list
    toKeep <- c(comparator, active.intervention, comparison)
    model.results <- model.results[names(model.results) %in% toKeep]
  }

  # Define population of interest ----

  # Whole UK Population (default)
  if (whole.uk.population){
    pop <- Define_targeted_population(min.bmi = arguments$bmi.target.min,
      max.bmi = arguments$bmi.target.max)
    pop <- mutate(pop, count = count * propTarget) %>%
      select(sex, ageGrp, count)
  } else {

    # User-defined population
    pc <- population.characteristics
    pop <- data.pop %>%
      mutate(count = number.in.cohort * ifelse(sex == "male",
        pc[1] * (pnorm(age, pc[2], pc[3]) - pnorm(age-1, pc[2], pc[3])),
        (1 - pc[1]) * (pnorm(age, pc[4], pc[5]) - pnorm(age-1, pc[4], pc[5]))),
        count = count * number.in.cohort / sum(count)) %>%
      mutate(ageGrp = as.character(cut(age, breaks = seq(0, 100, 5),
        right = FALSE)),
        ageGrp = ifelse(age == 100, "[95,100)", ageGrp)) %>%
      select(sex, ageGrp, count)
  }

  # Summarise across age bands ----

  # Sum results over fifteen year age bands
  results1 <- list()
  for (c in names(model.results)){
    results1[[c]] <- model.results[[c]] %>%
      left_join(., pop, by = c("sex", "ageGrp")) %>%
      mutate_at(outcomes, funs(. * count)) %>%
      mutate(ageGrp = as.character(cut(age, breaks = c(0, 20, 35, 50, 65, 80, Inf),
        labels = c("0-19", "20-34", "35-49", "50-64", "65-79", "80+"),
        right = FALSE))) %>%
      group_by(sex, ageGrp, year) %>%
      summarise_if(is.numeric, sum) %>%
      mutate_at(outcomes, funs(. / count)) %>%
      select(-age)
  }

  # Sum results across all ages, separately by sex
  results2 <- list()
  for (c in names(model.results)){
    results2[[c]] <- model.results[[c]] %>%
      left_join(., pop, by = c("sex", "ageGrp")) %>%
      mutate_at(outcomes, funs(. * count)) %>%
      group_by(year) %>%
      summarise_if(is.numeric, sum) %>%
      mutate_at(outcomes, funs(. / count)) %>%
      select(-age) %>%
      mutate(ageGrp = "All", sex = "All")
  }

  # Sum results across all ages by sex
  results3 <- list()
  for (c in names(model.results)){
    results3[[c]] <- model.results[[c]] %>%
      left_join(., pop, by = c("sex", "ageGrp")) %>%
      mutate_at(outcomes, funs(. * count)) %>%
      group_by(sex, year) %>%
      summarise_if(is.numeric, sum) %>%
      mutate_at(outcomes, funs(. / count)) %>%
      select(-age) %>%
      mutate(ageGrp = "All")
  }


  # Combine agg results with those by age and sex
  results <- list()
  for (c in names(model.results)){
    results[[c]] <- bind_rows(results2[[c]], results1[[c]], results3[[c]])
  }

  # Define total costs ----
  for (n in names(results)){
    if (costs.to.include == "disease-related nhs costs"){
      results[[n]] <- mutate(results[[n]],
        cost.total = trt.cost + nhs.cost.disease)
    } else if (costs.to.include == "unrelated nhs costs"){
      results[[n]] <- mutate(results[[n]],
        cost.total = trt.cost + nhs.cost.total)
    } else if (costs.to.include == "disease-related nhs and social care"){
      results[[n]] <- mutate(results[[n]],
        cost.total = trt.cost + nhs.cost.disease + sc.cost.disease)
    } else {
      results[[n]] <- mutate(results[[n]],
        cost.total = trt.cost + nhs.cost.total + sc.cost.total)
    }
  }

  # Create summary table ----
  if (produce.results.table){

    # Simple summary...
    if (!extended.summary){

      # variables to report (in order of presentation)
      report.vars <- c("Lx", "Lux", "cost.total")

      # create output template
      out.template <- data.frame(matrix(NA, nrow = length(report.vars) + 1, ncol = 4))
      names(out.template) <- c("outcome", comparator, active.intervention, comparison)
      out.template$outcome <- c("Life-years", "QALYs", "Total costs", "ICER")

      # empty outcome list
      out.list <- list()

      # Create gender and age groups
      sexL <- unique(results[[1]]$sex)
      ageL <- unique(results[[1]]$ageGrp)

      # Populate by sex & age groups
      for (s in sexL){
        for (a in ageL){

          if (!(s == "All" & a != "All")){
            # extract data
            out <- out.template
            for (m in names(out)[-1]){
              res.temp <- filter(results[[m]], sex == s, ageGrp == a)
              for (r in 1:length(report.vars)){
                out[r, m] <- res.temp %>%
                  filter(year == timeH) %>%
                  select(report.vars[r]) %>%
                  .[[report.vars[r]]]
              }
            }

            # calculate ICER
            out[out$outcome == "ICER", 4] <- out[3, 4] / out[2, 4]

            # clean up data
#            out[1:2, 2:3] <- round(out[1:2, 2:3], 2)
#            out[1:2, 4] <- round(out[1:2, 4], 3)
#            out[3, 2:3] <- round(out[3, 2:3])
#            out[3, 4] <- round(out[3, 4])
#            out[4, 4] <- round(out[4, 4])

            # store results in list
            out.list[[paste0("sex: ", s, "; age group: ", a)]] <- out
          }
        }
      }

    } else {

      # Extended summary...

      # cost vars to report
      cost.vars <- c("cost.total", "trt.cost")
      if (costs.to.include == "disease-related nhs costs"){
        cost.vars <- c(cost.vars, "nhs.cost.disease")
      } else if (costs.to.include == "unrelated nhs costs"){
        cost.vars <- c(cost.vars, "nhs.cost.total", "nhs.cost.disease", "nhs.cost.other")
      } else if (costs.to.include == "disease-related nhs and social care"){
        cost.vars <- c(cost.vars, "nhs.cost.disease", "sc.cost.disease")
      } else {
        cost.vars <- c(cost.vars, "nhs.cost.total", "nhs.cost.disease", "nhs.cost.other",
                       "sc.cost.total", "sc.cost.disease", "sc.cost.other")
      }

      # variables to report (in order of presentation)
      report.vars <- c("Lx", "Lux", cost.vars, outcomes[str_detect(outcomes, "ix")], "ICER")

      # create output template
      out.template <- data.frame(matrix(NA, nrow = length(report.vars), ncol = 4))
      names(out.template) <- c("outcome", comparator, active.intervention, comparison)
      out.template$outcome <- report.vars
      out.template$cat <- c(rep("health", 2), rep("cost", length(cost.vars)),
                            rep("ix", length(outcomes[str_detect(outcomes, "ix")])), "icer")

      # empty outcome list
      out.list <- list()

      # Create gender and age groups
      sexL <- unique(results[[1]]$sex)
      ageL <- unique(results[[1]]$ageGrp)

      # Populate by sex & age groups
      for (s in sexL){
        for (a in ageL){

          if (!(s == "All" & a != "All")){
            # extract data
            out <- out.template
            for (m in names(out)[c(-1,-5)]){
              res.temp <- filter(results[[m]], sex == s, ageGrp == a)
              for (r in 1:(length(report.vars)-1)){
                out[r, m] <- res.temp %>%
                  filter(year == timeH) %>%
                  select(report.vars[r]) %>%
                  .[[report.vars[r]]]
              }
            }

            # calculate ICER
            out[out$outcome == "ICER", 4] <- out[out$outcome == "cost.total", 4] / out[out$outcome == "Lux", 4]

            # clean up data
            results.table <- out.template
            results.table[out$cat=="health", 2:4] <- round(out[out$cat=="health", 2:4], 3)
            results.table[out$cat=="cost", 2:4] <- round(out[out$cat=="cost", 2:4])
            results.table[out$cat=="ix", 2:4] <- round(out[out$cat=="ix", 2:4] * 10^5)
            results.table[out$cat=="icer", 4] <- round(out[out$cat=="icer", 4])
            results.table$cat <- NULL

            # store results in list
            out.list[[paste0("sex: ", s, "; age group: ", a)]] <- out #results.table

          }
        }
      }
    }
  }

  # Data to return ----
  if (present.only.overall.results) out.list <- out.list[[1]]
  if (produce.results.table) out.list else results
}


# Summary including PSA ----
summary_psa <- function(model = NULL, psa.out = NULL){

  # Generate summary for deterministic results
  det.sum <- summary_table(
    model.results = model,
    comparator = "trt1",
    active.intervention = "trt2",
    produce.results.table = TRUE,
    extended.summary = TRUE)

  # For each group calculate standard deviation of probabilistic means ----

  # Output template
  out.template <- data.frame(matrix(NA, nrow = length(psa.out), ncol = nrow(det.sum[[1]])))
  names(out.template) <- det.sum[[1]][,1]

  # Generate outputs
  for (n in names(det.sum)){

    # Record output for each PSA iteration
    out.temp <- out.template
    for (i in 1:length(psa.out)){
      tempD <- psa.out[[i]][[n]]
      out.temp[i, ] <- tempD$diff_21
    }

    # Calculate standard deviations & CIs
    det.sum[[n]]$sd <- apply(out.temp, 2, sd, na.rm = T)
    det.sum[[n]]$lci <- det.sum[[n]]$diff_21 - qnorm(.975) * det.sum[[n]]$sd
    det.sum[[n]]$uci <- det.sum[[n]]$diff_21 + qnorm(.975) * det.sum[[n]]$sd
  }

  # Return output
  det.sum

}

