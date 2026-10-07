# ============================================================
# 01. Implementation of the distribution functions (Subbab 5.2)
#
# Evaluates every distribution function at a few points and draws the
# PDF, survival, and hazard functions for parameter sets that show the
# hazard shapes discussed in Subbab 2.1.4: decreasing (beta <= 1),
# increasing (large alpha), and modified unimodal (small alpha, beta > 1).
#
# Output: results/01_distribution/, figures/01_distribution/
# ============================================================

source("scripts/helpers.R")

result_dir <- output_dir("results", "01_distribution")
fig_dir <- output_dir("figures", "01_distribution")

# ------------------------------------------------------------
# 1. Example values of each function (alpha = 1, beta = 2, lambda = 0.5)
# ------------------------------------------------------------
alpha <- 1
beta <- 2
lambda <- 0.5
x <- c(0.5, 1, 1.5, 2, 3)

example_values <- data.frame(
  x = x,
  pdf = dnmeweibull(x, alpha, beta, lambda),
  cdf = pnmeweibull(x, alpha, beta, lambda),
  survival = snmeweibull(x, alpha, beta, lambda),
  hazard = hnmeweibull(x, alpha, beta, lambda)
)
save_table(example_values, file.path(result_dir, "example_values.csv"))

p <- c(0.10, 0.25, 0.50, 0.75, 0.90)
example_quantiles <- data.frame(
  p = p,
  quantile = qnmeweibull(p, alpha, beta, lambda),
  check_cdf_of_quantile = pnmeweibull(qnmeweibull(p, alpha, beta, lambda),
                                      alpha, beta, lambda)
)
save_table(example_quantiles, file.path(result_dir, "example_quantiles.csv"))

set.seed(2026)
example_random <- data.frame(
  index = 1:5,
  value = rnmeweibull(5, alpha, beta, lambda)
)
save_table(example_random, file.path(result_dir, "example_random.csv"))

# ------------------------------------------------------------
# 2. Shapes of the PDF, survival, and hazard functions
# ------------------------------------------------------------
shapes <- data.frame(
  label = c("alpha = 0,5; beta = 0,8; lambda = 1 (menurun)",
            "alpha = 5; beta = 2; lambda = 1 (meningkat)",
            "alpha = 0,05; beta = 2; lambda = 1 (unimodal termodifikasi)",
            "alpha = 1; beta = 1; lambda = 1 (menurun ke lambda)"),
  alpha = c(0.5, 5, 0.05, 1),
  beta = c(0.8, 2, 2, 1),
  lambda = c(1, 1, 1, 1)
)
cols <- c("#1b6ca8", "#d1495b", "#2e8b57", "#8c564b")
grid_x <- seq(0.01, 3, length.out = 500)

# One panel per parameter set so that every shape stays visible.
draw_family <- function(fun, ylab, file) {
  open_png(file.path(fig_dir, file), width = 2000, height = 1500)
  graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  for (i in seq_len(nrow(shapes))) {
    values <- fun(grid_x, shapes$alpha[i], shapes$beta[i], shapes$lambda[i])
    graphics::plot(grid_x, values, type = "l", lwd = 2, col = cols[i],
                   xlab = "x", ylab = ylab,
                   ylim = c(0, max(values[is.finite(values)])),
                   main = shapes$label[i], cex.main = 0.85)
  }
  grDevices::dev.off()
}

draw_family(dnmeweibull, "f(x)", "pdf_shapes.png")
draw_family(snmeweibull, "S(x)", "survival_shapes.png")
draw_family(hnmeweibull, "h(x)", "hazard_shapes.png")

# Hazard of the modified unimodal case on its own, where the
# increase-decrease-increase pattern is visible.
open_png(file.path(fig_dir, "hazard_modified_unimodal.png"))
graphics::plot(grid_x, hnmeweibull(grid_x, 0.05, 2, 1), type = "l", lwd = 2,
               col = cols[3], xlab = "x", ylab = "h(x)")
grDevices::dev.off()

# ------------------------------------------------------------
# 3. Random sample against the theoretical density
# ------------------------------------------------------------
set.seed(2026)
sample_x <- rnmeweibull(10000, alpha, beta, lambda)
open_png(file.path(fig_dir, "random_sample_histogram.png"))
graphics::hist(sample_x, breaks = 60, freq = FALSE, col = "grey90",
               border = "grey60", main = "", xlab = "x", ylab = "Kepadatan")
graphics::curve(dnmeweibull(x, alpha, beta, lambda), add = TRUE, lwd = 2,
                col = cols[1])
grDevices::dev.off()

cat("\nResults written to", result_dir, "and", fig_dir, "\n")
