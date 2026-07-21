# ============================================================
# 03_make_spatial_diagnostics.R
# Create formal spatial-pattern, transect and CSR diagnostics
#
# Run this script from the PROJECT ROOT.
# It uses one reproducible direct-coordinate example per scenario.
# It does NOT run the full Monte Carlo experiment.
# ============================================================


# ------------------------------------------------------------
# 0. Clean session and set output folder
# ------------------------------------------------------------

rm(list = ls())
graphics.off()
options(stringsAsFactors = FALSE)

seed_diagnostics <- 987L
set.seed(seed_diagnostics)

diagnostic_dir <- file.path(
  "outputs",
  "spatial_diagnostics"
)

dir.create(
  diagnostic_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Formal settings used in the main experiment
# ------------------------------------------------------------

marker_total_nominal <- 105000L

example_ratio <- 2
example_n_target <- as.integer(
  example_ratio * marker_total_nominal
)
example_n_marker <- marker_total_nominal

y_low <- 0.83
y_high <- 0.88

# Use one formal stopping count for the explanatory figure.
x_stop_example <- 333L

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
  file.path(
    diagnostic_dir,
    "scenario_parameters.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 2. Coordinate generator
# ------------------------------------------------------------

generate_xy <- function(n, scenario) {
  
  if (scenario == "CSR") {
    
    x <- runif(n)
    y <- runif(n)
    
  } else if (scenario == "Beta centre-high") {
    
    x <- rbeta(n, 4, 4)
    y <- rbeta(n, 4, 4)
    
  } else if (scenario == "Beta right-high") {
    
    x <- rbeta(n, 5, 1.5)
    y <- runif(n)
    
  } else if (scenario == "Beta left-high") {
    
    x <- rbeta(n, 1.5, 5)
    y <- runif(n)
    
  } else {
    
    stop("Unknown scenario: ", scenario)
  }
  
  data.frame(
    x = x,
    y = y
  )
}


# ------------------------------------------------------------
# 3. Generate one reproducible full point pattern per scenario
# ------------------------------------------------------------

simulation_list <- vector(
  "list",
  length(scenarios)
)

names(simulation_list) <- scenarios

for (i in seq_along(scenarios)) {
  
  scenario <- scenarios[i]
  
  # Separate deterministic seed for each scenario.
  set.seed(seed_diagnostics + i * 1000L)
  
  target <- generate_xy(
    example_n_target,
    scenario
  )
  
  marker <- generate_xy(
    example_n_marker,
    scenario
  )
  
  simulation_list[[scenario]] <- list(
    scenario = scenario,
    target = target,
    marker = marker
  )
}


# ------------------------------------------------------------
# 4. Save one figure in both PNG and PDF formats
# ------------------------------------------------------------

save_both_formats <- function(
    filename_stub,
    draw_function,
    png_width = 1800,
    png_height = 1500,
    png_res = 180,
    pdf_width = 10,
    pdf_height = 8
) {
  
  png(
    file.path(
      diagnostic_dir,
      paste0(filename_stub, ".png")
    ),
    width = png_width,
    height = png_height,
    res = png_res
  )
  
  draw_function()
  dev.off()
  
  pdf(
    file.path(
      diagnostic_dir,
      paste0(filename_stub, ".pdf")
    ),
    width = pdf_width,
    height = pdf_height
  )
  
  draw_function()
  dev.off()
}


# ------------------------------------------------------------
# 5. Figure A: formal spatial-pattern illustration
# ------------------------------------------------------------

draw_spatial_patterns <- function() {
  
  max_target_display <- 1500L
  max_marker_display <- 750L
  
  par(
    mfrow = c(2, 2),
    mar = c(4, 4, 3.5, 1),
    oma = c(0, 0, 2.5, 0)
  )
  
  for (i in seq_along(scenarios)) {
    
    scenario <- scenarios[i]
    sim <- simulation_list[[scenario]]
    
    set.seed(seed_diagnostics + i * 10L)
    
    target_id <- sample(
      seq_len(nrow(sim$target)),
      min(max_target_display, nrow(sim$target))
    )
    
    marker_id <- sample(
      seq_len(nrow(sim$marker)),
      min(max_marker_display, nrow(sim$marker))
    )
    
    plot(
      sim$target$x[target_id],
      sim$target$y[target_id],
      xlim = c(0, 1),
      ylim = c(0, 1),
      pch = 16,
      cex = 0.48,
      col = rgb(0, 0, 1, 0.35),
      xlab = "x",
      ylab = "y",
      main = scenario
    )
    
    points(
      sim$marker$x[marker_id],
      sim$marker$y[marker_id],
      pch = 16,
      cex = 0.48,
      col = rgb(1, 0, 0, 0.35)
    )
    
    abline(
      h = c(y_low, y_high),
      lty = 2,
      col = "black"
    )
    
    legend(
      "topright",
      legend = c(
        "Target subset",
        "Marker subset",
        "Transect band"
      ),
      col = c(
        rgb(0, 0, 1, 0.35),
        rgb(1, 0, 0, 0.35),
        "black"
      ),
      pch = c(16, 16, NA),
      lty = c(NA, NA, 2),
      bty = "n",
      cex = 0.76
    )
  }
  
  mtext(
    "Direct-coordinate spatial patterns used in the simulation",
    outer = TRUE,
    font = 2,
    cex = 1.18
  )
}

save_both_formats(
  filename_stub = "fig_spatial_patterns",
  draw_function = draw_spatial_patterns
)


# ------------------------------------------------------------
# 6. Figure B: formal linear-transect examples
# ------------------------------------------------------------

transect_summary <- data.frame(
  scenario = character(0),
  x_stop = integer(0),
  x_end = numeric(0),
  marker_count = integer(0),
  target_count_in_band = integer(0),
  stringsAsFactors = FALSE
)

draw_linear_transects <- function() {
  
  max_target_display <- 1500L
  max_marker_display <- 750L
  
  par(
    mfrow = c(2, 2),
    mar = c(4, 4, 4.3, 1),
    oma = c(0, 0, 2.6, 0)
  )
  
  local_summary <- data.frame(
    scenario = character(0),
    x_stop = integer(0),
    x_end = numeric(0),
    marker_count = integer(0),
    target_count_in_band = integer(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(scenarios)) {
    
    scenario <- scenarios[i]
    sim <- simulation_list[[scenario]]
    
    target_in_band <- (
      sim$target$y >= y_low &
        sim$target$y <= y_high
    )
    
    marker_in_band <- (
      sim$marker$y >= y_low &
        sim$marker$y <= y_high
    )
    
    band_target_x <- sim$target$x[target_in_band]
    band_marker_x <- sim$marker$x[marker_in_band]
    
    if (length(band_target_x) < x_stop_example) {
      
      plot.new()
      title(
        main = paste0(
          scenario,
          "\nStopping count not reached"
        )
      )
      
      next
    }
    
    ordered_target_x <- sort.int(
      band_target_x,
      partial = x_stop_example
    )
    
    x_end <- ordered_target_x[x_stop_example]
    marker_count <- sum(band_marker_x <= x_end)
    
    local_summary <- rbind(
      local_summary,
      data.frame(
        scenario = scenario,
        x_stop = x_stop_example,
        x_end = x_end,
        marker_count = marker_count,
        target_count_in_band =
          length(band_target_x),
        stringsAsFactors = FALSE
      )
    )
    
    set.seed(seed_diagnostics + i * 20L)
    
    target_id <- sample(
      seq_len(nrow(sim$target)),
      min(max_target_display, nrow(sim$target))
    )
    
    marker_id <- sample(
      seq_len(nrow(sim$marker)),
      min(max_marker_display, nrow(sim$marker))
    )
    
    plot(
      sim$target$x[target_id],
      sim$target$y[target_id],
      xlim = c(0, 1),
      ylim = c(0, 1),
      pch = 16,
      cex = 0.46,
      col = rgb(0, 0, 1, 0.32),
      xlab = "x",
      ylab = "y",
      main = paste0(
        scenario,
        "\nx_stop = ", x_stop_example,
        ", x_end = ", round(x_end, 3),
        ", n = ", marker_count
      )
    )
    
    points(
      sim$marker$x[marker_id],
      sim$marker$y[marker_id],
      pch = 16,
      cex = 0.46,
      col = rgb(1, 0, 0, 0.32)
    )
    
    rect(
      0,
      y_low,
      x_end,
      y_high,
      border = "black",
      col = rgb(1, 0, 0, 0.18),
      lwd = 2
    )
    
    abline(
      h = c(y_low, y_high),
      lty = 2,
      col = "black"
    )
    
    legend(
      "topright",
      legend = c(
        "Target subset",
        "Marker subset",
        "Counting rectangle"
      ),
      col = c(
        rgb(0, 0, 1, 0.32),
        rgb(1, 0, 0, 0.32),
        "black"
      ),
      pch = c(16, 16, NA),
      lty = c(NA, NA, 1),
      bty = "n",
      cex = 0.72
    )
  }
  
  mtext(
    "Linear transect with a shared target stopping rule",
    outer = TRUE,
    font = 2,
    cex = 1.18
  )
  
  assign(
    "transect_summary",
    local_summary,
    envir = .GlobalEnv
  )
}

save_both_formats(
  filename_stub = "fig_linear_transect_examples",
  draw_function = draw_linear_transects
)

write.csv(
  transect_summary,
  file.path(
    diagnostic_dir,
    "transect_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. CSR theoretical pairwise-distance CDF on 0 <= t <= 1
# ------------------------------------------------------------

H_csr_unit_square <- function(t) {
  
  ifelse(
    t < 0,
    0,
    ifelse(
      t <= 1,
      pi * t^2 -
        (8 / 3) * t^3 +
        (1 / 2) * t^4,
      NA_real_
    )
  )
}


# ------------------------------------------------------------
# 8. Figure C: descriptive CSR discrepancy diagnostic
# ------------------------------------------------------------

diagnostic_subset_size <- 1000L

csr_gap_summary <- data.frame(
  scenario = character(0),
  subset_size = integer(0),
  csr_gap = numeric(0),
  stringsAsFactors = FALSE
)

draw_csr_diagnostic <- function() {
  
  par(
    mfrow = c(2, 2),
    mar = c(4, 4, 4.2, 1),
    oma = c(0, 0, 2.6, 0)
  )
  
  t_values <- seq(
    0,
    1,
    length.out = 500
  )
  
  theoretical_values <- H_csr_unit_square(
    t_values
  )
  
  local_gap_summary <- data.frame(
    scenario = character(0),
    subset_size = integer(0),
    csr_gap = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(scenarios)) {
    
    scenario <- scenarios[i]
    target <- simulation_list[[scenario]]$target
    
    set.seed(seed_diagnostics + i * 30L)
    
    subset_size <- min(
      diagnostic_subset_size,
      nrow(target)
    )
    
    id <- sample(
      seq_len(nrow(target)),
      subset_size
    )
    
    sampled_points <- target[
      id,
      c("x", "y")
    ]
    
    pairwise_distances <- as.vector(
      dist(sampled_points)
    )
    
    empirical_cdf <- ecdf(
      pairwise_distances
    )
    
    empirical_values <- empirical_cdf(
      t_values
    )
    
    csr_gap <- max(
      abs(
        empirical_values -
          theoretical_values
      ),
      na.rm = TRUE
    )
    
    local_gap_summary <- rbind(
      local_gap_summary,
      data.frame(
        scenario = scenario,
        subset_size = subset_size,
        csr_gap = csr_gap,
        stringsAsFactors = FALSE
      )
    )
    
    plot(
      t_values,
      empirical_values,
      type = "l",
      xlim = c(0, 1),
      ylim = c(0, 1),
      col = "blue",
      lwd = 2,
      xlab = "Inter-event distance t",
      ylab = "H(t)",
      main = paste0(
        scenario,
        "\nCSR gap = ",
        round(csr_gap, 3)
      )
    )
    
    lines(
      t_values,
      theoretical_values,
      col = "red",
      lwd = 2,
      lty = 2
    )
    
    legend(
      "bottomright",
      legend = c(
        "Empirical CDF",
        "Diggle Eqn (2.2)"
      ),
      col = c("blue", "red"),
      lwd = 2,
      lty = c(1, 2),
      bty = "n",
      cex = 0.78
    )
  }
  
  mtext(
    "CSR diagnostic using target subset: inter-event distances",
    outer = TRUE,
    font = 2,
    cex = 1.18
  )
  
  assign(
    "csr_gap_summary",
    local_gap_summary,
    envir = .GlobalEnv
  )
}

save_both_formats(
  filename_stub = "fig_csr_diagnostic_inter_event",
  draw_function = draw_csr_diagnostic
)

write.csv(
  csr_gap_summary,
  file.path(
    diagnostic_dir,
    "csr_gap_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Finish
# ------------------------------------------------------------

cat(
  "\nSpatial diagnostics saved in:\n",
  diagnostic_dir,
  "\n\nImportant reporting note:\n",
  "Only displayed subsets are plotted for legibility. ",
  "Stopping positions and marker counts use the full simulated ",
  "point patterns.\n",
  "The CSR discrepancy is descriptive, not a formal p-value.\n"
)