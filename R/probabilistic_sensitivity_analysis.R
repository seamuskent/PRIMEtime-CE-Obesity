#' Probabilisitic sensitivity analysis
#'
#' \code{primetime_psa} Estimates the PRIMEtime model \emph{n} times representing the uncertainty in selected input parameters using Monte Carlo Simulation.
#'
#'
#' @inheritParams summary_table
#' @return Returns the same output as \code{primetime} but for each Monte Carlo simulation. The output can be accessed by other functions to produce standardised outputs (e.g. summary tables), or can be used directly by the user for bespoke analyses.
#' @examples
#' primetime_psa(n = 1000, model = primetime.results)
#' @export

primetime_psa <- function(n = 1000, model.results = NULL, costs.to.include = "disease-related nhs costs", all.years = FALSE){

  # Extract arguments from model
  inputs <- model.results$arguments

  # change from determinstic to probabilsitic
  inputs$deterministic <- FALSE

  # Loop over iterations
  out.list <- foreach::foreach(i = 1:n,
    .packages = "primetimeCE") %do% {

      # Report iteration number
      print(paste0("Iteration number: ", i))

      # Specify probabilistic model
      tempModel <- rlang::invoke(primetime, inputs)

      # Detailed summary
      if (all.years){
        tempSum <- list()
        for (y in c(1, seq(5, inputs$time.horizon, 5))){
          tempSum[[y]] <- summary_table(
            model.results = tempModel,
            comparator = "trt1",
            active.intervention = "trt2",
            nicely.presented.results = FALSE,
            extended.summary = TRUE,
            costs.to.include = costs.to.include,
            timeH = y)
        }
      } else {
        tempSum <- summary_table(
          model.results = tempModel,
          comparator = "trt1",
          active.intervention = "trt2",
          nicely.presented.results = FALSE,
          extended.summary = TRUE,
          costs.to.include = costs.to.include)
      }
      
      # Define output to save
      tempSum

    }

  #return output
  out.list

}
