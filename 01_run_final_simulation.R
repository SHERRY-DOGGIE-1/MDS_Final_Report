# ============================================================
# 01_run_final_simulation.R
# MDS Final Report: unified direct-coordinate Monte Carlo study
#
# Run this script from the PROJECT ROOT.
#
# Formal spatial patterns:
#   1. CSR
#   2. Beta centre-high: X,Y ~ Beta(4,4)
#   3. Beta right-high:  X ~ Beta(5,1.5), Y ~ Uniform(0,1)
#   4. Beta left-high:   X ~ Beta(1.5,5), Y ~ Uniform(0,1)
#
# Common settings:
#   nominal marker total = 105000
#   marker-dose SD       = 3112
#   ratios               = 1, 2, 5, 10
#   stopping counts      = 333, 800
#
# Each generated point pattern is reused for both stopping counts.
# ============================================================


# ------------------------------------------------------------
# 0. Clean session and reproducibility
# ------------------------------------------------------------

rm(list = ls())
graphics.off()
options(stringsAsFactors = FALSE)

seed_main <- 123L


# ------------------------------------------------------------
# 1. Choose pilot or full run
# ------------------------------------------------------------

# Start with "pilot". After checking all outputs, change to "full".
run_mode <- "full"

n_repetitions <- if (run_mode == "pilot") {
  20L
} else {
  400L
}

# FALSE means completed scenario-ratio files are skipped on rerun.
overwrite_existing_settings <- FALSE


# ------------------------------------------------------------
# 2. Output folders
# ------------------------------------------------------------

output_root <- "outputs"
run_tag <- paste0(run_mode, "_R", n_repetitions)

run_output_dir <- file.path(output_root, run_tag)
raw_setting_dir <- file.path(run_output_dir, "raw_settings")

dir.create(output_root, showWarnings = FALSE)
dir.create(run_output_dir, showWarnings = FALSE)
dir.create(raw_setting_dir, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Formal experimental settings
# ------------------------------------------------------------

V <- 1

marker_total_nominal <- 105000L
marker_dose_sd <- 3112
marker_relative_sd <- marker_dose_sd / marker_total_nominal

ratio_values <- c(1, 2, 5, 10)
x_stop_values <- c(333L, 800L)

y_low <- 0.83
y_high <- 0.88

scenarios <- c(
  "CSR",
  "Beta centre-high",
  "Beta right-high",
  "Beta left-high"
)

scenario_parameters <- data.frame(
  scenario = scenarios,
  x_distribution = c(
    "Uniform(0,1)",
    "Beta(4,4)",
    "Beta(5,1.5)",
    "Beta(1.5,5)"
  ),
  y_distribution = c(
    "Uniform(0,1)",
    "Beta(4,4)",
    "Uniform(0,1)",
    "Uniform(0,1)"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  scenario_parameters,
  file.path(run_output_dir, "scenario_parameters.csv"),
  row.names = FALSE
)

parameter_table <- data.frame(
  parameter = c(
    "Monte Carlo repetitions per scenario-ratio setting",
    "Nominal marker total",
    "Marker-dose standard deviation",
    "Marker-dose proportional standard deviation",
    "Target-to-marker ratios",
    "Stopping counts",
    "Transect lower y-bound",
    "Transect upper y-bound",
    "Sample volume"
  ),
  value = c(
    as.character(n_repetitions),
    as.character(marker_total_nominal),
    as.character(marker_dose_sd),
    format(marker_relative_sd, digits = 8),
    paste(ratio_values, collapse = ", "),
    paste(x_stop_values, collapse = ", "),
    as.character(y_low),
    as.character(y_high),
    as.character(V)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  parameter_table,
  file.path(run_output_dir, "simulation_parameters.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Coordinate generator
# ------------------------------------------------------------

generate_xy <- function(n, scenario) {
  
  n <- as.integer(n)
  
  if (n < 1L) {
    stop("n must be a positive integer.")
  }
  
  if (scenario == "CSR") {
    
    x <- runif(n)
    y <- runif(n)
    
  } else if (scenario == "Beta centre-high") {
    
    # Controlled separable centre concentration.
    x <- rbeta(n, shape1 = 4, shape2 = 4)
    y <- rbeta(n, shape1 = 4, shape2 = 4)
    
  } else if (scenario == "Beta right-high") {
    
    x <- rbeta(n, shape1 = 5, shape2 = 1.5)
    y <- runif(n)
    
  } else if (scenario == "Beta left-high") {
    
    x <- rbeta(n, shape1 = 1.5, shape2 = 5)
    y <- runif(n)
    
  } else {
    
    stop("Unknown scenario: ", scenario)
  }
  
  list(x = x, y = y)
}


# ------------------------------------------------------------
# 5. Safe filename helper
# ------------------------------------------------------------

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}


# ------------------------------------------------------------
# 6. Run one scenario-ratio setting
# ------------------------------------------------------------

run_one_setting <- function(
    scenario,
    ratio_value,
    n_repetitions,
    setting_seed
) {
  
  set.seed(setting_seed)
  
  n_target <- as.integer(
    round(ratio_value * marker_total_nominal)
  )
  
  c_true <- n_target / V
  
  n_stops <- length(x_stop_values)
  n_rows <- n_repetitions * n_stops
  
  results <- data.frame(
    scenario = rep(scenario, n_rows),
    ratio_value = rep(ratio_value, n_rows),
    setting_seed = rep(setting_seed, n_rows),
    repetition = rep(seq_len(n_repetitions), each = n_stops),
    x_stop = rep(x_stop_values, times = n_repetitions),
    n_target_total = rep(n_target, n_rows),
    n_marker_actual = rep(NA_integer_, n_rows),
    n_target_band = rep(NA_integer_, n_rows),
    x_end = rep(NA_real_, n_rows),
    x_count = rep(NA_integer_, n_rows),
    n_count = rep(NA_integer_, n_rows),
    c_true = rep(c_true, n_rows),
    c_hat = rep(NA_real_, n_rows),
    sigma_L = rep(NA_real_, n_rows),
    valid = rep(FALSE, n_rows),
    invalid_reason = rep(NA_character_, n_rows),
    stringsAsFactors = FALSE
  )
  
  setting_start <- Sys.time()
  
  for (m in seq_len(n_repetitions)) {
    
    if (
      m == 1L ||
      m %% 25L == 0L ||
      m == n_repetitions
    ) {
      message(
        "  ", scenario,
        " | ratio = ", ratio_value,
        " | repetition ", m, "/", n_repetitions
      )
    }
    
    # Realised marker total for this simulated preparation.
    n_marker_actual <- as.integer(
      round(
        rnorm(
          1,
          mean = marker_total_nominal,
          sd = marker_dose_sd
        )
      )
    )
    
    n_marker_actual <- max(n_marker_actual, 1L)
    
    # Directly generate all target and marker coordinates.
    target <- generate_xy(n_target, scenario)
    marker <- generate_xy(n_marker_actual, scenario)
    
    # Restrict to the maximal horizontal transect band.
    target_in_band <- (
      target$y >= y_low &
        target$y <= y_high
    )
    
    marker_in_band <- (
      marker$y >= y_low &
        marker$y <= y_high
    )
    
    band_target_x <- target$x[target_in_band]
    band_marker_x <- marker$x[marker_in_band]
    
    n_target_band <- length(band_target_x)
    
    feasible_stops <- x_stop_values[
      x_stop_values <= n_target_band
    ]
    
    if (length(feasible_stops) > 0L) {
      
      # Partial sorting is sufficient because only selected order
      # statistics are required.
      ordered_target_x <- sort.int(
        band_target_x,
        partial = feasible_stops
      )
      
    } else {
      
      ordered_target_x <- numeric(0)
    }
    
    for (j in seq_along(x_stop_values)) {
      
      x_stop <- x_stop_values[j]
      row_id <- (m - 1L) * n_stops + j
      
      results$n_marker_actual[row_id] <- n_marker_actual
      results$n_target_band[row_id] <- n_target_band
      
      # Invalid case 1: the stopping count cannot be reached.
      if (n_target_band < x_stop) {
        
        results$invalid_reason[row_id] <- "not_enough_targets"
        next
      }
      
      x_end <- ordered_target_x[x_stop]
      x_count <- x_stop
      
      # Direct count of marker points in the final rectangle.
      n_count <- sum(band_marker_x <= x_end)
      
      results$x_end[row_id] <- x_end
      results$x_count[row_id] <- x_count
      results$n_count[row_id] <- n_count
      
      # Invalid case 2: concentration estimate is undefined.
      if (n_count == 0L) {
        
        results$invalid_reason[row_id] <- "zero_marker_count"
        next
      }
      
      c_hat <- (
        x_count *
          marker_total_nominal /
          (n_count * V)
      )
      
      sigma_L <- 100 * sqrt(
        marker_relative_sd^2 +
          1 / x_count +
          1 / n_count
      )
      
      results$c_hat[row_id] <- c_hat
      results$sigma_L[row_id] <- sigma_L
      results$valid[row_id] <- TRUE
    }
    
    rm(
      target,
      marker,
      target_in_band,
      marker_in_band,
      band_target_x,
      band_marker_x,
      ordered_target_x
    )
    
    if (m %% 25L == 0L) {
      invisible(gc(verbose = FALSE))
    }
  }
  
  elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      setting_start,
      units = "secs"
    )
  )
  
  attr(results, "elapsed_seconds") <- elapsed_seconds
  results
}


# ------------------------------------------------------------
# 7. Run all scenario-ratio settings with checkpoints
# ------------------------------------------------------------

run_grid <- expand.grid(
  scenario_index = seq_along(scenarios),
  ratio_index = seq_along(ratio_values),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

timing_log <- data.frame(
  scenario = character(0),
  ratio_value = numeric(0),
  n_repetitions = integer(0),
  elapsed_seconds = numeric(0),
  output_file = character(0),
  stringsAsFactors = FALSE
)

for (g in seq_len(nrow(run_grid))) {
  
  scenario_index <- run_grid$scenario_index[g]
  ratio_index <- run_grid$ratio_index[g]
  
  scenario <- scenarios[scenario_index]
  ratio_value <- ratio_values[ratio_index]
  
  setting_seed <- as.integer(
    seed_main +
      scenario_index * 100000L +
      ratio_index * 1000L
  )
  
  setting_stub <- paste0(
    safe_name(scenario),
    "_ratio_",
    ratio_value
  )
  
  setting_rds <- file.path(
    raw_setting_dir,
    paste0(setting_stub, ".rds")
  )
  
  setting_csv <- file.path(
    raw_setting_dir,
    paste0(setting_stub, ".csv")
  )
  
  if (
    file.exists(setting_rds) &&
    !overwrite_existing_settings
  ) {
    
    message(
      "\nSkipping completed setting: ",
      scenario,
      " | ratio = ",
      ratio_value
    )
    
    existing <- readRDS(setting_rds)
    elapsed_seconds <- attr(existing, "elapsed_seconds")
    
    if (is.null(elapsed_seconds)) {
      elapsed_seconds <- NA_real_
    }
    
    rm(existing)
    
  } else {
    
    message(
      "\nStarting setting: ",
      scenario,
      " | ratio = ",
      ratio_value
    )
    
    setting_result <- run_one_setting(
      scenario = scenario,
      ratio_value = ratio_value,
      n_repetitions = n_repetitions,
      setting_seed = setting_seed
    )
    
    elapsed_seconds <- attr(
      setting_result,
      "elapsed_seconds"
    )
    
    saveRDS(
      setting_result,
      setting_rds
    )
    
    write.csv(
      setting_result,
      setting_csv,
      row.names = FALSE
    )
    
    rm(setting_result)
    invisible(gc(verbose = FALSE))
  }
  
  timing_log <- rbind(
    timing_log,
    data.frame(
      scenario = scenario,
      ratio_value = ratio_value,
      n_repetitions = n_repetitions,
      elapsed_seconds = elapsed_seconds,
      output_file = setting_rds,
      stringsAsFactors = FALSE
    )
  )
  
  write.csv(
    timing_log,
    file.path(run_output_dir, "timing_log.csv"),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 8. Combine all raw setting outputs
# ------------------------------------------------------------

setting_files <- list.files(
  raw_setting_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

if (length(setting_files) == 0L) {
  stop("No setting files were found.")
}

raw_results <- do.call(
  rbind,
  lapply(setting_files, readRDS)
)

row.names(raw_results) <- NULL

raw_results$scenario <- factor(
  raw_results$scenario,
  levels = scenarios
)

raw_results <- raw_results[
  order(
    raw_results$scenario,
    raw_results$ratio_value,
    raw_results$repetition,
    raw_results$x_stop
  ),
]

raw_results$scenario <- as.character(
  raw_results$scenario
)

write.csv(
  raw_results,
  file.path(run_output_dir, "all_raw_results.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Summarise one scenario-ratio-stopping-count group
# ------------------------------------------------------------

summarise_group <- function(df) {
  
  valid_df <- df[df$valid, , drop = FALSE]
  
  m_planned <- nrow(df)
  m_eff <- nrow(valid_df)
  
  failed_stop <- sum(
    df$invalid_reason == "not_enough_targets",
    na.rm = TRUE
  )
  
  zero_marker <- sum(
    df$invalid_reason == "zero_marker_count",
    na.rm = TRUE
  )
  
  success_rate <- m_eff / m_planned
  
  scenario <- df$scenario[1]
  ratio_value <- df$ratio_value[1]
  x_stop <- df$x_stop[1]
  n_target_total <- df$n_target_total[1]
  c_true <- df$c_true[1]
  
  mean_marker_total_actual <- mean(
    df$n_marker_actual,
    na.rm = TRUE
  )
  
  mean_target_band <- mean(
    df$n_target_band,
    na.rm = TRUE
  )
  
  if (m_eff <= 1L) {
    
    return(
      data.frame(
        scenario = scenario,
        ratio_value = ratio_value,
        x_stop = x_stop,
        n_target_total = n_target_total,
        c_true = c_true,
        m_planned = m_planned,
        m_eff = m_eff,
        success_rate = success_rate,
        failed_stop = failed_stop,
        zero_marker = zero_marker,
        mean_marker_total_actual = mean_marker_total_actual,
        mean_target_band = mean_target_band,
        mean_x_end = NA_real_,
        mean_n = NA_real_,
        mean_c_hat = NA_real_,
        c_jensen = NA_real_,
        rel_bias_raw = NA_real_,
        rel_bias_jensen = NA_real_,
        abs_rel_bias_jensen = NA_real_,
        mean_sigma = NA_real_,
        sigma_exact = NA_real_,
        signed_gap = NA_real_,
        absolute_gap = NA_real_,
        relative_gap = NA_real_,
        mc_se_mean_sigma = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  mean_x_end <- mean(valid_df$x_end)
  mean_n <- mean(valid_df$n_count)
  mean_c_hat <- mean(valid_df$c_hat)
  
  c_jensen <- (
    x_stop *
      marker_total_nominal /
      (mean_n * V)
  )
  
  rel_bias_raw <- (
    mean_c_hat - c_true
  ) / c_true
  
  rel_bias_jensen <- (
    c_jensen - c_true
  ) / c_true
  
  abs_rel_bias_jensen <- abs(
    rel_bias_jensen
  )
  
  mean_sigma <- mean(
    valid_df$sigma_L
  )
  
  sigma_exact <- 100 * sqrt(
    sum(
      (valid_df$c_hat - c_true)^2
    ) /
      (m_eff - 1L)
  ) / c_true
  
  signed_gap <- mean_sigma - sigma_exact
  absolute_gap <- abs(signed_gap)
  relative_gap <- absolute_gap / sigma_exact
  
  mc_se_mean_sigma <- (
    sd(valid_df$sigma_L) /
      sqrt(m_eff)
  )
  
  data.frame(
    scenario = scenario,
    ratio_value = ratio_value,
    x_stop = x_stop,
    n_target_total = n_target_total,
    c_true = c_true,
    m_planned = m_planned,
    m_eff = m_eff,
    success_rate = success_rate,
    failed_stop = failed_stop,
    zero_marker = zero_marker,
    mean_marker_total_actual = mean_marker_total_actual,
    mean_target_band = mean_target_band,
    mean_x_end = mean_x_end,
    mean_n = mean_n,
    mean_c_hat = mean_c_hat,
    c_jensen = c_jensen,
    rel_bias_raw = rel_bias_raw,
    rel_bias_jensen = rel_bias_jensen,
    abs_rel_bias_jensen = abs_rel_bias_jensen,
    mean_sigma = mean_sigma,
    sigma_exact = sigma_exact,
    signed_gap = signed_gap,
    absolute_gap = absolute_gap,
    relative_gap = relative_gap,
    mc_se_mean_sigma = mc_se_mean_sigma,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 10. Create final summary table
# ------------------------------------------------------------

split_key <- interaction(
  raw_results$scenario,
  raw_results$ratio_value,
  raw_results$x_stop,
  drop = TRUE,
  lex.order = TRUE
)

summary_results <- do.call(
  rbind,
  lapply(
    split(raw_results, split_key),
    summarise_group
  )
)

row.names(summary_results) <- NULL

summary_results$scenario <- factor(
  summary_results$scenario,
  levels = scenarios
)

summary_results <- summary_results[
  order(
    summary_results$x_stop,
    summary_results$scenario,
    summary_results$ratio_value
  ),
]

summary_results$scenario <- as.character(
  summary_results$scenario
)

write.csv(
  summary_results,
  file.path(run_output_dir, "summary_results.csv"),
  row.names = FALSE
)

print(summary_results)


# ------------------------------------------------------------
# 11. Runtime summary
# ------------------------------------------------------------

total_elapsed <- sum(
  timing_log$elapsed_seconds,
  na.rm = TRUE
)

cat("\n==================================================\n")
cat("Simulation completed.\n")
cat("Mode:", run_mode, "\n")
cat("Repetitions per scenario-ratio setting:",
    n_repetitions, "\n")
cat("Output folder:", run_output_dir, "\n")
cat("Measured setting time (seconds):",
    round(total_elapsed, 1), "\n")

if (run_mode == "pilot" && n_repetitions > 0L) {
  
  estimated_full_seconds <- (
    total_elapsed *
      400 /
      n_repetitions
  )
  
  cat(
    "Simple linear estimate for R = 400:",
    round(estimated_full_seconds / 3600, 2),
    "hours\n"
  )
  
  cat(
    "This is a rough estimate only.\n"
  )
}

cat("==================================================\n")