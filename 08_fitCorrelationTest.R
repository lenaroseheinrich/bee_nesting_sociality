# #===============================================================================
# # Correlated evolution tests
# #===============================================================================
# # Tests for correlated evolution between bee sociality and nesting strategy
# # using corHMM fitCorrelationTest across binary and multi-state encodings.
# # Outputs model tables and fitted models for comparison.
# #===============================================================================
# # Setup
# #-------------------------------------------------------------------------------
# rm(list = ls())

# repo_dir <- "/home/lenarh/repos/bee_nesting_sociality"
# data_dir <- "/home/lenarh/data/bee_nesting_sociality"

# curated_data_dir <- file.path(data_dir, "curated_data")
# results_dir <- file.path(data_dir, "results", "corHMM", "correlation_tests")
# dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# library(corHMM)

# source(file.path(repo_dir, "00_utility_functions.R"))

# #-------------------------------------------------------------------------------
# # Load data
# #-------------------------------------------------------------------------------
# traits <- read.csv(
#   file.path(curated_data_dir, "bee_traits_clean.csv")
# )

# phy <- read.tree(
#   file.path(curated_data_dir, "bee_tree_pruned.tre")
# )

# dat <- traits

# #-------------------------------------------------------------------------------
# # Trait encodings to test
# #-------------------------------------------------------------------------------
# sociality_scoring <- c("sociality_binary", "sociality_three_states")
# nest_scoring <- c("nest_binary", "nest_three_states")

# #-------------------------------------------------------------------------------
# # Correlated evolution tests
# #-------------------------------------------------------------------------------
# for (sociality_index in 1:length(sociality_scoring)) {
#   for (nest_index in 1:length(nest_scoring)) {
#     dat1 <- dat[, c(
#       "tips",
#       sociality_scoring[sociality_index],
#       nest_scoring[nest_index]
#     )]

#     corHMM_fits <- corHMM:::fitCorrelationTest(phy, dat1)
#     corHMM_tbl <- corHMM:::getModelTable(corHMM_fits)

#     write.csv(
#       corHMM_tbl,
#       file = file.path(
#         results_dir,
#         paste0(
#           "corHMM_fits_cortest_table_",
#           sociality_scoring[sociality_index],
#           "_",
#           nest_scoring[nest_index],
#           ".csv"
#         )
#       )
#     )

#     save(
#       corHMM_fits,
#       file = file.path(
#         results_dir,
#         paste0(
#           "corHMM_fits_cortest_",
#           sociality_scoring[sociality_index],
#           "_",
#           nest_scoring[nest_index],
#           ".Rsave"
#         )
#       )
#     )
#   }
# }

# #-------------------------------------------------------------------------------
# # Likelihood ratio test
# #-------------------------------------------------------------------------------
# teststat <- -2 * (corHMM_tbl$lnLik[2] - corHMM_tbl$lnLik[4])
# p.val <- pchisq(teststat, df = 8, lower.tail = FALSE)

#===============================================================================
# Correlated evolution tests
#===============================================================================
# Tests for correlated evolution between bee sociality and nesting strategy
# using corHMM fitCorrelationTest across alternative binary and three-state
# scoring methods.
#
# The three-state variables retain parasites as a separate state.
#
# Outputs:
#   1. A model-comparison table for each scoring-method pair
#   2. An Rsave containing the fitted models for each scoring-method pair
#   3. A summary table of likelihood-ratio tests comparing the hidden-Markov
#      independent and hidden-Markov correlated models
#===============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

set.seed(1234)

data_dir <- "/home/lenarh/data/bee_nesting_sociality"

curated_data_dir <- file.path(
  data_dir,
  "curated_data"
)

results_dir <- file.path(
  data_dir,
  "results",
  "corHMM",
  "correlation_tests"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

library(corHMM)

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
traits <- read.csv(
  file.path(
    curated_data_dir,
    "bee_traits_clean.csv"
  )
)

phy <- ape::read.tree(
  file.path(
    curated_data_dir,
    "bee_tree_pruned.tre"
  )
)

#-------------------------------------------------------------------------------
# Trait scoring methods to test
#-------------------------------------------------------------------------------
sociality_scoring <- c(
  "sociality_binary",
  "sociality_three_states"
)

nest_scoring <- c(
  "nest_binary",
  "nest_three_states"
)

# Check that all required columns are present
required_columns <- c(
  "tips",
  sociality_scoring,
  nest_scoring
)

missing_columns <- setdiff(
  required_columns,
  names(traits)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the trait dataset: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Check for duplicated species names
if (anyDuplicated(traits$tips) > 0) {
  stop(
    "The trait dataset contains duplicated values in the 'tips' column."
  )
}

#-------------------------------------------------------------------------------
# Correlated evolution tests
#-------------------------------------------------------------------------------
lrt_results_list <- list()
result_index <- 1

for (sociality_variable in sociality_scoring) {
  for (nest_variable in nest_scoring) {
    message(
      "\nRunning correlation test for: ",
      sociality_variable,
      " + ",
      nest_variable
    )

    # Select species names and the two trait-scoring columns
    dat1 <- traits[, c(
      "tips",
      sociality_variable,
      nest_variable
    )]

    # Fit independent and correlated models, with and without hidden rates
    corHMM_fits <- corHMM:::fitCorrelationTest(
      phy,
      dat1
    )

    # Create model-comparison table
    corHMM_tbl <- corHMM:::getModelTable(
      corHMM_fits
    )

    # Add model names as an explicit column
    corHMM_tbl <- data.frame(
      model = names(corHMM_fits),
      corHMM_tbl,
      row.names = NULL,
      check.names = FALSE
    )

    #---------------------------------------------------------------------------
    # Save model table
    #---------------------------------------------------------------------------
    table_filename <- paste0(
      "corHMM_fits_cortest_table_",
      sociality_variable,
      "_",
      nest_variable,
      ".csv"
    )

    write.csv(
      corHMM_tbl,
      file = file.path(
        results_dir,
        table_filename
      ),
      row.names = FALSE
    )

    #---------------------------------------------------------------------------
    # Save fitted models
    #---------------------------------------------------------------------------
    fits_filename <- paste0(
      "corHMM_fits_cortest_",
      sociality_variable,
      "_",
      nest_variable,
      ".Rsave"
    )

    save(
      corHMM_fits,
      file = file.path(
        results_dir,
        fits_filename
      )
    )

    #---------------------------------------------------------------------------
    # Likelihood-ratio test
    #
    # Compare:
    #   hidden-Markov independent model
    # versus
    #   hidden-Markov correlated model
    #---------------------------------------------------------------------------

    hidden_independent_row <- grep(
      "hidden.*independent",
      corHMM_tbl$model,
      ignore.case = TRUE
    )

    hidden_correlated_row <- grep(
      "hidden.*correlated",
      corHMM_tbl$model,
      ignore.case = TRUE
    )

    # Confirm that exactly one of each model was found
    if (
      length(hidden_independent_row) != 1 ||
        length(hidden_correlated_row) != 1
    ) {
      stop(
        paste0(
          "Could not uniquely identify the hidden independent and hidden ",
          "correlated models for ",
          sociality_variable,
          " + ",
          nest_variable,
          ". Model names were: ",
          paste(corHMM_tbl$model, collapse = ", ")
        )
      )
    }

    independent_lnLik <- corHMM_tbl$lnLik[
      hidden_independent_row
    ]

    correlated_lnLik <- corHMM_tbl$lnLik[
      hidden_correlated_row
    ]

    independent_np <- corHMM_tbl$np[
      hidden_independent_row
    ]

    correlated_np <- corHMM_tbl$np[
      hidden_correlated_row
    ]

    # Likelihood-ratio statistic
    teststat <- 2 *
      (correlated_lnLik -
        independent_lnLik)

    # Difference in number of estimated parameters
    df_lrt <- correlated_np -
      independent_np

    if (df_lrt <= 0) {
      stop(
        paste0(
          "The correlated model does not have more parameters than the ",
          "independent model for ",
          sociality_variable,
          " + ",
          nest_variable,
          "."
        )
      )
    }

    if (teststat < 0) {
      warning(
        paste0(
          "The correlated model has a lower likelihood than the independent ",
          "model for ",
          sociality_variable,
          " + ",
          nest_variable,
          ". Check model convergence."
        )
      )
    }

    p_value <- pchisq(
      teststat,
      df = df_lrt,
      lower.tail = FALSE
    )

    # Store LRT results
    lrt_results_list[[result_index]] <- data.frame(
      sociality_scoring = sociality_variable,
      nest_scoring = nest_variable,
      independent_model = corHMM_tbl$model[
        hidden_independent_row
      ],
      correlated_model = corHMM_tbl$model[
        hidden_correlated_row
      ],
      independent_np = independent_np,
      correlated_np = correlated_np,
      independent_lnLik = independent_lnLik,
      correlated_lnLik = correlated_lnLik,
      LR_statistic = teststat,
      df = df_lrt,
      p_value = p_value
    )

    result_index <- result_index + 1
  }
}

#-------------------------------------------------------------------------------
# Combine and save likelihood-ratio test results
#-------------------------------------------------------------------------------
lrt_results <- do.call(
  rbind,
  lrt_results_list
)

rownames(lrt_results) <- NULL

write.csv(
  lrt_results,
  file = file.path(
    results_dir,
    "corHMM_correlation_test_LRT_results.csv"
  ),
  row.names = FALSE
)

save(
  lrt_results,
  file = file.path(
    results_dir,
    "corHMM_correlation_test_LRT_results.Rsave"
  )
)

# Display final results
print(lrt_results)
