# ============================================================
# 02_make_report_figures.R
# Create Final Report result figures from summary_results.csv
#
# Run this script from the PROJECT ROOT after Script 01.
# It does NOT rerun the Monte Carlo simulation.
# ============================================================


# ------------------------------------------------------------
# 0. Clean session
# ------------------------------------------------------------

rm(list = ls())
graphics.off()
options(stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 1. Select which simulation output to plot
# ------------------------------------------------------------

# Use "pilot_R20" while testing.
# Change to "full_R400" after the formal simulation is complete.
run_tag <- "full_R400"

input_file <- file.path(
  "outputs",
  run_tag,
  "summary_results.csv"
)

if (!file.exists(input_file)) {
  stop(
    "Cannot find: ", input_file,
    "\nRun 01_run_final_simulation.R first, ",
    "or correct run_tag."
  )
}

figure_dir <- file.path(
  "outputs",
  run_tag,
  "report_figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Read and check results
# ------------------------------------------------------------

summary_results <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "scenario",
  "ratio_value",
  "x_stop",
  "mean_n",
  "mean_sigma",
  "sigma_exact",
  "signed_gap",
  "relative_gap",
  "abs_rel_bias_jensen"
)

missing_columns <- setdiff(
  required_columns,
  names(summary_results)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing columns in summary_results.csv: ",
    paste(missing_columns, collapse = ", ")
  )
}

scenarios <- c(
  "CSR",
  "Beta centre-high",
  "Beta right-high",
  "Beta left-high"
)

ratio_values <- sort(
  unique(summary_results$ratio_value)
)

x_stop_values <- sort(
  unique(summary_results$x_stop)
)

summary_results$scenario <- factor(
  summary_results$scenario,
  levels = scenarios
)


# ------------------------------------------------------------
# 3. Consistent plotting settings
# ------------------------------------------------------------

scenario_colours <- c(
  "CSR" = "black",
  "Beta centre-high" = "blue",
  "Beta right-high" = "red",
  "Beta left-high" = "darkgreen"
)

scenario_symbols <- c(
  "CSR" = 19,
  "Beta centre-high" = 17,
  "Beta right-high" = 15,
  "Beta left-high" = 18
)

stop_symbols <- c(
  "333" = 16,
  "800" = 17
)


# ------------------------------------------------------------
# 4. Save one figure in both PNG and PDF formats
# ------------------------------------------------------------

save_both_formats <- function(
    filename_stub,
    draw_function,
    png_width = 1800,
    png_height = 850,
    png_res = 180,
    pdf_width = 10,
    pdf_height = 5.5
) {
  
  png(
    file.path(
      figure_dir,
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
      figure_dir,
      paste0(filename_stub, ".pdf")
    ),
    width = pdf_width,
    height = pdf_height
  )
  
  draw_function()
  dev.off()
}


# ------------------------------------------------------------
# 5. Two-panel metric plot by stopping count
# ------------------------------------------------------------

draw_metric_by_stop <- function(
    data,
    metric,
    y_label,
    plot_title,
    add_zero_line = FALSE
) {
  
  par(
    mfrow = c(1, length(x_stop_values)),
    mar = c(4.8, 5.0, 4.0, 1.0),
    oma = c(0, 0, 2.4, 0)
  )
  
  for (x_stop in x_stop_values) {
    
    tmp_stop <- data[
      data$x_stop == x_stop,
      ,
      drop = FALSE
    ]
    
    y_values <- tmp_stop[[metric]]
    y_range <- range(y_values, na.rm = TRUE)
    
    if (
      !all(is.finite(y_range)) ||
      diff(y_range) == 0
    ) {
      y_range <- c(0, 1)
    }
    
    y_padding <- 0.06 * diff(y_range)
    
    plot(
      ratio_values,
      rep(NA_real_, length(ratio_values)),
      type = "n",
      xlim = range(ratio_values),
      ylim = y_range + c(-y_padding, y_padding),
      xaxt = "n",
      xlab = "Target-to-marker ratio",
      ylab = y_label,
      main = paste0("x_stop = ", x_stop)
    )
    
    axis(
      1,
      at = ratio_values,
      labels = ratio_values
    )
    
    if (add_zero_line) {
      abline(
        h = 0,
        lty = 2,
        col = "grey50"
      )
    }
    
    for (scenario in scenarios) {
      
      tmp <- tmp_stop[
        tmp_stop$scenario == scenario,
        ,
        drop = FALSE
      ]
      
      tmp <- tmp[
        order(tmp$ratio_value),
        ,
        drop = FALSE
      ]
      
      lines(
        tmp$ratio_value,
        tmp[[metric]],
        type = "b",
        pch = scenario_symbols[scenario],
        col = scenario_colours[scenario],
        lwd = 2
      )
    }
    
    legend(
      "topleft",
      legend = scenarios,
      col = scenario_colours[scenarios],
      pch = scenario_symbols[scenarios],
      lwd = 2,
      bty = "n",
      cex = 0.80
    )
  }
  
  mtext(
    plot_title,
    outer = TRUE,
    font = 2,
    cex = 1.12
  )
}


# ------------------------------------------------------------
# 6. Main report figure: relative precision gap
# ------------------------------------------------------------

save_both_formats(
  filename_stub = "fig_relative_precision_gap",
  draw_function = function() {
    draw_metric_by_stop(
      data = summary_results,
      metric = "relative_gap",
      y_label = expression(
        "|" * bar(sigma)[L] - sigma[exact] * "|" /
          sigma[exact]
      ),
      plot_title =
        "Relative error of the formula-based precision estimate"
    )
  }
)


# ------------------------------------------------------------
# 7. Signed gap: direction of formula error
# ------------------------------------------------------------

save_both_formats(
  filename_stub = "fig_signed_precision_gap",
  draw_function = function() {
    draw_metric_by_stop(
      data = summary_results,
      metric = "signed_gap",
      y_label = expression(
        bar(sigma)[L] - sigma[exact]
      ),
      plot_title =
        "Signed precision gap (negative values indicate underestimation)",
      add_zero_line = TRUE
    )
  }
)


# ------------------------------------------------------------
# 8. Mean counted marker number
# ------------------------------------------------------------

save_both_formats(
  filename_stub = "fig_mean_counted_markers",
  draw_function = function() {
    draw_metric_by_stop(
      data = summary_results,
      metric = "mean_n",
      y_label = expression(
        "Mean counted marker number  " * bar(n)
      ),
      plot_title =
        "Mean marker count in the final counting rectangle"
    )
  }
)


# ------------------------------------------------------------
# 9. Concentration bias
# ------------------------------------------------------------

save_both_formats(
  filename_stub = "fig_concentration_bias",
  draw_function = function() {
    draw_metric_by_stop(
      data = summary_results,
      metric = "abs_rel_bias_jensen",
      y_label = expression(
        "|" * hat(c)[Jensen] - c[true] * "|" /
          c[true]
      ),
      plot_title =
        "Absolute relative bias of the ratio-of-means concentration summary"
    )
  }
)


# ------------------------------------------------------------
# 10. Pooled diagnostic: relative gap versus mean marker count
# ------------------------------------------------------------

# Keep FALSE unless the report text provides a precise,
# properly cited justification for using n = 100 as a reference.
show_reference_n100 <- FALSE

draw_gap_vs_nbar <- function() {
  
  data <- summary_results
  
  par(
    mar = c(4.8, 5.2, 3.5, 1.2)
  )
  
  x_range <- range(
    data$mean_n,
    na.rm = TRUE
  )
  
  y_range <- range(
    data$relative_gap,
    na.rm = TRUE
  )
  
  plot(
    data$mean_n,
    data$relative_gap,
    log = "x",
    type = "n",
    xlim = x_range,
    ylim = y_range,
    xlab = expression(
      "Mean counted marker number  " * bar(n)
    ),
    ylab = expression(
      "|" * bar(sigma)[L] - sigma[exact] * "|" /
        sigma[exact]
    ),
    main =
      "Relative precision gap versus mean counted marker number"
  )
  
  for (scenario in scenarios) {
    
    for (x_stop in x_stop_values) {
      
      sub <- data[
        data$scenario == scenario &
          data$x_stop == x_stop,
        ,
        drop = FALSE
      ]
      
      points(
        sub$mean_n,
        sub$relative_gap,
        pch = stop_symbols[as.character(x_stop)],
        col = scenario_colours[scenario],
        cex = 1.35
      )
    }
  }
  
  if (show_reference_n100) {
    
    abline(
      v = 100,
      lty = 3,
      col = "grey50"
    )
    
    text(
      100,
      max(data$relative_gap, na.rm = TRUE),
      expression(bar(n) == 100),
      pos = 4,
      col = "grey40",
      cex = 0.8
    )
  }
  
  legend(
    "topright",
    legend = scenarios,
    col = scenario_colours[scenarios],
    pch = 15,
    bty = "n",
    title = "Spatial pattern",
    cex = 0.82
  )
  
  legend(
    "bottomleft",
    legend = paste0("x_stop = ", x_stop_values),
    col = "grey25",
    pch = stop_symbols[
      as.character(x_stop_values)
    ],
    bty = "n",
    title = "Stopping count",
    cex = 0.82
  )
}

save_both_formats(
  filename_stub = "fig_gap_vs_nbar_pooled",
  draw_function = draw_gap_vs_nbar,
  png_width = 1500,
  png_height = 1100,
  pdf_width = 8,
  pdf_height = 6
)


# ------------------------------------------------------------
# 11. Export concise tables for report drafting
# ------------------------------------------------------------

main_table_columns <- c(
  "scenario",
  "ratio_value",
  "x_stop",
  "m_eff",
  "success_rate",
  "mean_n",
  "mean_sigma",
  "sigma_exact",
  "signed_gap",
  "relative_gap",
  "abs_rel_bias_jensen"
)

main_results_table <- summary_results[
  ,
  main_table_columns,
  drop = FALSE
]

write.csv(
  main_results_table,
  file.path(
    figure_dir,
    "table_main_results.csv"
  ),
  row.names = FALSE
)

cat(
  "\nReport figures and table saved in:\n",
  figure_dir,
  "\n"
)