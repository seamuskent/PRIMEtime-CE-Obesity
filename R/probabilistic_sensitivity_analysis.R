

primetime_psa <- function(n = 1000, model = NULL){

  # Extract arguments from model
  inputs <- model$arguments

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
      tempSum <- summary_table(
        model.results = tempModel,
        comparator = "trt1",
        active.intervention = "trt2",
        nicely.presented.results = FALSE,
        extended.summary = TRUE
      )

      # Define output to save
      tempSum

    }

  #return output
  out.list

}
