## Script used to run sleepdata analysis
## Written by: Jasper Ginn
##
## For more information about the analysis, please visit
##   <https://github.com/JasperHG90/sleepsimR-documentation>
##
## Version: 0.01b

library(sleepsimR)
library(argparser, quietly=TRUE)
library(logger)

# Set up a parser
args <- arg_parser("Run a single chain of the sleep data analysis used in my thesis. For more information and examples, see <https://github.com/JasperHG90/sleepsimR-documentation>.")
# Add arguments
args <- add_argument(args, "iterations", help="Number of MCMC iterations that are used to sample the posterior distribution of the parameters.", type = "numeric", default = NULL)
args <- add_argument(args, "burn_in", help="Number of samples that will be discarded (burn-in samples) at the beginning of the chain.", type = "numeric", default = NULL)
args <- add_argument(args, "variables", help="1 to 4 character names of the variables to be used in the analysis. Three values are expected. Accepted variable names are: 'EEG_Fpz_Cz_mean_theta', 'EEG_Fpz_Cz_mean_beta', 'EOG_min_beta', 'EOG_median_theta'", nargs=3, default = c("EEG_Fpz_Cz_mean_theta", "EOG_min_beta", "EOG_median_theta"))
# Parse
argv <- parse_args(args)

# Set up the logger
log_info("Application is starting up ...")

# Main function
main <- function(iterations = argv$iterations, burn_in = argv$burn_in, variables = argv$variables) {
  # Read sleep data
  sleepdat <- readRDS("app/sleep_data_subset.rds")
  # Load summary statistics
  sumstats <- readRDS("app/summary_statistics.rds")
  total_var <- readRDS("app/total_variance.rds")
  # Round values
  sumstats$mmvar <- round(sumstats$mmvar, 2)
  total_var$tvar <- round(total_var$tvar, 2)
  # Make unique id
  uid <- uuid::UUIDgenerate()
  # Iterations / burn_in > 0
  # Check if passed
  if(iterations <= 0) {
    stop("Number of iterations cannot be less than 0.")
  }
  if(burn_in <= 0) {
    stop("Number of burn-in samples cannot be less than 0.")
  }
  if(burn_in >= iterations) {
    stop("Number of iterations must be larger than the number of burn-in samples.")
  }
  if(all("EEG_Fpz_Cz_mean_theta", "EEG_Fpz_Cz_mean_beta") %in% variables) {
    stop("Cannot use both 'EEG_Fpz_Cz_mean_theta' and 'EEG_Fpz_Cz_mean_beta' in analysis.")
  }
  if(!all(variables %in% c("EEG_Fpz_Cz_mean_theta", "EOG_min_beta", "EOG_median_theta", "EEG_Fpz_Cz_mean_beta"))) {
    stop('Variables must match three of "EEG_Fpz_Cz_mean_theta", "EOG_min_beta", "EOG_median_theta", "EEG_Fpz_Cz_mean_beta" exactly.')
  }
  # Force ordering of variables
  if('EEG_Fpz_Cz_mean_theta' %in% variables) {
    order <- c(
      "EEG_Fpz_Cz_mean_theta" = 1,
      "EOG_min_beta" = 2,
      "EOG_median_theta" = 3
    )
  } else {
    order <- c(
      "EEG_Fpz_Cz_mean_beta" = 1,
      "EOG_min_beta" = 2,
      "EOG_median_theta" = 3
    )
  }
  variables <- c("EOG_min_beta", "EOG_median_theta", "EEG_Fpz_Cz_mean_beta")
  variables <- variables[order(match(variables,names(order)))]
  # Model properties
  mprop = list(
    "m" = 3,
    "n_dep" = 3
  )
  # Subset summary stats and variances
  sss <- sumstats[sumstats$variable %in% variables,]
  tvv <- total_var[total_var$variable %in% variables,]
  # Make list of hyperprior means
  em_prior_means <- list(
    sss$mmvar[1:3],
    sss$mmvar[4:6],
    sss$mmvar[7:9]
  )
  # Draw seed and set
  seed <- sample.int(10000000, 1)
  set.seed(seed)
  # Initial values
  ## TPM gamma
  diag_value <- runif(1, 0.6, 0.8)
  gam <- diag(diag_value, mprop$m)
  gam[lower.tri(gam) | upper.tri(gam)] <- (1-diag_value) / 2
  # Check if all 1
  if(!all(apply(gam, 1, sum) == 1)) {
    stop("Initial values for the TPM do not exactly sum to 1 on all states.")
  }
  # For the 3 continuous emission distributions
  start_EM <- list(
    # Gamma
    start_gamma,
    #EEG_Fpz_Cz_max_gamma
    matrix(c( sss$mmvar[1] + runif(1, -0.05, 0.05), total_var$tvar[1]+ runif(1, -0.05, 0.05),
              sss$mmvar[2] + runif(1, -0.05, 0.05), total_var$tvar[2]+ runif(1, -0.05, 0.05),
              sss$mmvar[3] + runif(1, -0.05, 0.05), total_var$tvar[3] + runif(1, -0.05, 0.05)),
           nrow=mprop$m,
           ncol=2,
           byrow = TRUE),
    # EOG_median_theta
    matrix(c( sss$mmvar[4]+ runif(1, -0.05, 0.05), total_var$tvar[4]+ runif(1, -0.05, 0.05),
              sss$mmvar[5]+ runif(1, -0.05, 0.05), total_var$tvar[5]+ runif(1, -0.05, 0.05),
              sss$mmvar[6]+ runif(1, -0.05, 0.05), total_var$tvar[6]+ runif(1, -0.05, 0.05)),
           nrow=mprop$m,
           ncol=2,
           byrow = TRUE),
    # EOG_min_beta
    matrix(c( sss$mmvar[7]+ runif(1, -0.05, 0.05), total_var$tvar[7]+ runif(1, -0.05, 0.05),
              sss$mmvar[8]+ runif(1, -0.05, 0.05), total_var$tvar[8]+ runif(1, -0.05, 0.05),
              sss$mmvar[9]+ runif(1, -0.05, 0.05), total_var$tvar[9]+ runif(1, -0.05, 0.05)),
           nrow=mprop$m,
           ncol=2,
           byrow = TRUE)
  )
  # Remove sleep states from data
  sleepdat_subs <- sleepdat
  sleepdat_subs$sleep_state <- NULL
  sleepdat_subs <- sleepdat_subs[,c("id", variables)]
  # Run model
  mod <- sleepsimR::run_mHMM(as.matrix(sleepdat_subs),
                             start_EM,
                             mprop,
                             em_prior_means,
                             seed,
                             mcmc_iterations=iterations,
                             mcmc_burn_in=burn_in,
                             show_progress = FALSE,
                             order_data = FALSE)
  # Make results
  results <- list(
    "unique_id" = uid,
    "seed_initial_values" = seed,
    "initial_values" = start_EM,
    "prior_values" = em_prior_means,
    "variables" = variables,
    "iterations" = iterations,
    "burn_in" = burn_in,
    "model" = mod
  )
  # Save model
  saveRDS(mod, paste0("/var/sleepsimr_sleepdata_analysis/model_", uid, ".rds"))
}



