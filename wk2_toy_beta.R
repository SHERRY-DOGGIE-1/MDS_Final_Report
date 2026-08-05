# ============================================================
# Part B toy baseline: Scenario B with shared Beta patterns
# No Stockmarr uncertainty yet
#
# Aim:
#   1. Generate CSR baseline and three Beta spatial patterns
#   2. Targets and markers share the same spatial pattern
#   3. Run CSR diagnostic from new_single_run
# ============================================================


# ------------------------------------------------------------
# 0. Basic settings
# ------------------------------------------------------------

set.seed(123)

out_dir <- "partB_beta_toy_figures"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Fixed total numbers for toy demonstration.
# This is conditional on the total number of points.
# For a fully Poisson version, set use_poisson_total <- TRUE.
use_poisson_total <- FALSE

mu_target <- 12000
mu_marker <- 6000

get_total_count <- function(mu) {
  if (use_poisson_total) {
    return(rpois(1, mu))
  } else {
    return(mu)
  }
}

# Plot only a subset, otherwise the figure is too crowded
plot_target_size <- 1200
plot_marker_size <- 600

# CSR diagnostic subset size
check_size <- 700

# Linear transect setting, same idea as Part A
y_low <- 0.83
y_high <- 0.88
x_stop <- 100


# ------------------------------------------------------------
# 1. Generate points
# ------------------------------------------------------------

generate_xy <- function(n, scenario) {
  
  if (scenario == "CSR") {
    x <- runif(n)
    y <- runif(n)
  }
  
  if (scenario == "Beta centre-high") {
    # X,Y ~ Beta(3,3)
    # Points concentrate near the centre of the slide
    x <- rbeta(n, 3, 3)
    y <- rbeta(n, 3, 3)
  }
  
  if (scenario == "Beta right-high") {
    # X ~ Beta(3,1.5), Y ~ Uniform(0,1)
    # Points concentrate toward the right side
    x <- rbeta(n, 3, 1.5)
    y <- runif(n)
  }
  
  if (scenario == "Beta left-high") {
    # X ~ Beta(1.5,3), Y ~ Uniform(0,1)
    # Points concentrate toward the left side
    x <- rbeta(n, 1.5, 3)
    y <- runif(n)
  }
  
  data.frame(x = x, y = y)
}


generate_one_scenario <- function(scenario) {
  
  n_target <- get_total_count(mu_target)
  n_marker <- get_total_count(mu_marker)
  
  target <- generate_xy(n_target, scenario)
  marker <- generate_xy(n_marker, scenario)
  
  list(
    scenario = scenario,
    target = target,
    marker = marker,
    n_target = n_target,
    n_marker = n_marker
  )
}


scenarios <- c(
  "CSR",
  "Beta centre-high",
  "Beta right-high",
  "Beta left-high"
)

sim_list <- lapply(scenarios, generate_one_scenario)
names(sim_list) <- scenarios


# ------------------------------------------------------------
# 2. Plot point patterns for group meeting
# ------------------------------------------------------------

plot_one_pattern <- function(sim_obj) {
  
  target <- sim_obj$target
  marker <- sim_obj$marker
  
  target_id <- sample(seq_len(nrow(target)),
                      min(plot_target_size, nrow(target)))
  marker_id <- sample(seq_len(nrow(marker)),
                      min(plot_marker_size, nrow(marker)))
  
  plot(target$x[target_id], target$y[target_id],
       xlim = c(0, 1), ylim = c(0, 1),
       pch = 16,
       cex = 0.55,
       col = rgb(0, 0, 1, 0.45),
       xlab = "x",
       ylab = "y",
       main = sim_obj$scenario)
  
  points(marker$x[marker_id], marker$y[marker_id],
         pch = 16,
         cex = 0.55,
         col = rgb(1, 0, 0, 0.45))
  
  legend("topright",
         legend = c("Target", "Marker"),
         col = c(rgb(0, 0, 1, 0.45), rgb(1, 0, 0, 0.45)),
         pch = 16,
         bty = "n",
         cex = 0.9)
}


png(file.path(out_dir, "fig1_beta_point_patterns.png"),
    width = 1800, height = 1500, res = 180)

par(mfrow = c(2, 2),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0))

for (sc in scenarios) {
  plot_one_pattern(sim_list[[sc]])
}

mtext("Scenario B toy examples: targets and markers share the same spatial pattern",
      outer = TRUE,
      cex = 1.2,
      font = 2)

dev.off()


# ------------------------------------------------------------
# 3. CSR diagnostic: inter-event distances
# Same idea as new_single_run
# ------------------------------------------------------------

H_theoretical <- function(t) {
  ifelse(
    t <= 1,
    pi * t^2 - (8 / 3) * t^3 + t^4 / 2,
    1 / 3 -
      2 * t^2 +
      4 * sqrt(t^2 - 1) * (2 * t^2 + 1) / 3 +
      2 * t^2 * asin(1 / t) -
      (t^4 / 2) * (2 * t^(-2) - 1)
  )
}


csr_inter_event_diagnostic <- function(target_points,
                                       scenario_name,
                                       check_size = 700) {
  
  id <- sample(seq_len(nrow(target_points)),
               min(check_size, nrow(target_points)))
  
  check_x <- target_points$x[id]
  check_y <- target_points$y[id]
  
  inter_event_dist <- as.vector(dist(cbind(check_x, check_y)))
  edf_fun <- ecdf(inter_event_dist)
  
  t_vals <- seq(0, sqrt(2), length.out = 500)
  h_vals <- H_theoretical(t_vals)
  edf_vals <- edf_fun(t_vals)
  
  # A simple visual discrepancy number.
  # Not a formal p-value.
  csr_gap <- max(abs(edf_vals - h_vals), na.rm = TRUE)
  
  plot(edf_fun,
       main = paste0(scenario_name, "\nCSR gap = ", round(csr_gap, 3)),
       xlab = "Inter-event distance t",
       ylab = "H(t)",
       col = "blue",
       lwd = 2)
  
  lines(t_vals, h_vals,
        col = "red",
        lwd = 2,
        lty = 2)
  
  legend("bottomright",
         legend = c("Empirical CDF", "CSR theoretical curve"),
         col = c("blue", "red"),
         lwd = 2,
         lty = c(1, 2),
         bty = "n",
         cex = 0.85)
  
  return(csr_gap)
}


png(file.path(out_dir, "fig2_csr_diagnostic_inter_event.png"),
    width = 1800, height = 1500, res = 180)

par(mfrow = c(2, 2),
    mar = c(4, 4, 3.5, 1),
    oma = c(0, 0, 2, 0))

csr_gap_vec <- numeric(length(scenarios))
names(csr_gap_vec) <- scenarios

for (sc in scenarios) {
  csr_gap_vec[sc] <- csr_inter_event_diagnostic(
    target_points = sim_list[[sc]]$target,
    scenario_name = sc,
    check_size = check_size
  )
}

mtext("CSR diagnostic using target subset: inter-event distances",
      outer = TRUE,
      cex = 1.2,
      font = 2)

dev.off()


# ------------------------------------------------------------
# 4. linear transect visualization
# This is only for PPT illustration
# ------------------------------------------------------------

plot_linear_transect <- function(sim_obj,
                                 y_low = 0.83,
                                 y_high = 0.88,
                                 x_stop = 100) {
  
  target <- sim_obj$target
  marker <- sim_obj$marker
  
  target_in_band <- (target$y >= y_low & target$y <= y_high)
  band_target_x <- target$x[target_in_band]
  
  if (length(band_target_x) < x_stop) {
    plot.new()
    title(main = paste0(sim_obj$scenario, "\nNot enough targets in band"))
    return(NA)
  }
  
  x_end <- sort(band_target_x, partial = x_stop)[x_stop]
  
  target_id <- sample(seq_len(nrow(target)),
                      min(plot_target_size, nrow(target)))
  marker_id <- sample(seq_len(nrow(marker)),
                      min(plot_marker_size, nrow(marker)))
  
  plot(target$x[target_id], target$y[target_id],
       xlim = c(0, 1), ylim = c(0, 1),
       pch = 16,
       cex = 0.5,
       col = rgb(0, 0, 1, 0.35),
       xlab = "x",
       ylab = "y",
       main = paste0(sim_obj$scenario, "\n",
                     "x_stop = ", x_stop, ", x_end = ", round(x_end, 3)))
  
  points(marker$x[marker_id], marker$y[marker_id],
         pch = 16,
         cex = 0.5,
         col = rgb(1, 0, 0, 0.35))
  
  rect(0, y_low, x_end, y_high,
       border = "black",
       col = rgb(1, 0, 0, 0.20),
       lwd = 2)
  
  abline(h = c(y_low, y_high),
         lty = 2,
         col = "black")
  
  legend("topright",
         legend = c("Target subset", "Marker subset", "Counting rectangle"),
         col = c(rgb(0, 0, 1, 0.35), rgb(1, 0, 0, 0.35), "black"),
         pch = c(16, 16, NA),
         lty = c(NA, NA, 1),
         bty = "n",
         cex = 0.8)
  
  return(x_end)
}


png(file.path(out_dir, "fig4_linear_transect_examples.png"),
    width = 1800, height = 1500, res = 180)

par(mfrow = c(2, 2),
    mar = c(4, 4, 3.5, 1),
    oma = c(0, 0, 2, 0))

x_end_vec <- numeric(length(scenarios))
names(x_end_vec) <- scenarios

for (sc in scenarios) {
  x_end_vec[sc] <- plot_linear_transect(
    sim_obj = sim_list[[sc]],
    y_low = y_low,
    y_high = y_high,
    x_stop = x_stop
  )
}

mtext("Linear transect toy examples under CSR and shared Beta patterns",
      outer = TRUE,
      cex = 1.2,
      font = 2)

dev.off()



cat("\nFigures saved in folder:", out_dir, "\n")
cat("1. fig1_beta_point_patterns.png\n")
cat("2. fig2_csr_diagnostic_inter_event.png\n")
cat("3. fig3_csr_diagnostic_nearest_neighbour.png\n")
cat("4. fig4_linear_transect_examples.png\n")
