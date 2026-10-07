# ============================================================
# Evaluasi survival dan hazard: NME-Weibull vs Weibull
# Data GBSG2: waktu sampai kekambuhan atau kematian (hari),
# atau sampai pengamatan terakhir untuk pasien tersensor.
# Jalankan dari folder utama proyek nmeweibull (berisi R/ dan data/).
# Fungsi distribusi diambil langsung dari kode package Anda.
# ============================================================

if (!file.exists("R/distribution.R")) {
  stop("Buka proyek nmeweibull.Rproj terlebih dahulu.")
}
source("R/distribution.R")
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Pasang package survival: install.packages('survival')")
}

# Folder terpisah agar hasil skrip lama tetap tersedia.
result_dir <- "results/"
figure_dir <- "figures/"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Baca dan periksa data
# status = 1: kekambuhan atau kematian.
# status = 0: masih hidup tanpa kekambuhan pada pengamatan terakhir
#             (tersensor kanan).
# Analisis mengasumsikan sensor noninformatif dan observasi independen.
# ------------------------------------------------------------

dat <- read.csv("data/gbsg2_recurrence_free_time.csv")
if (!all(c("time", "status") %in% names(dat))) {
  stop("Data harus memiliki kolom time dan status.")
}
if (!is.numeric(dat$time) || !is.numeric(dat$status)) {
  stop("Kolom time dan status harus berupa angka.")
}
if (any(!is.finite(dat$time) | dat$time <= 0) ||
    any(!is.finite(dat$status) | !(dat$status %in% c(0, 1)))) {
  stop("Periksa data: time harus positif dan finite; status hanya 0 atau 1.")
}
time <- dat$time
status <- dat$status
n <- length(time)
if (n < 3 || sum(status) == 0 || length(unique(time)) < 2) {
  stop("Data terlalu sedikit/tidak bervariasi, atau tidak memiliki kejadian.")
}

# ------------------------------------------------------------
# 2. Ringkasan data
# Median waktu teramati mencakup waktu kejadian dan waktu sensor;
# angka ini bukan estimasi median survival.
# ------------------------------------------------------------

descriptive_stats <- data.frame(
  n = n,
  events = sum(status == 1),
  censored = sum(status == 0),
  censoring_percent = 100 * mean(status == 0),
  minimum_observed_days = min(time),
  median_observed_days = median(time),
  maximum_observed_days = max(time)
)
write.csv(descriptive_stats, file.path(result_dir, "ringkasan_data.csv"),
          row.names = FALSE)
print(descriptive_stats)

# ------------------------------------------------------------
# 3. Estimasi NME-Weibull dengan likelihood tersensor
# L = hasil kali f(t) untuk kejadian dan S(t) untuk sensor.
# fitnmeweibull() yang sekarang belum menerima status sensor.
# Waktu diskalakan saat optimisasi, lalu lambda dikembalikan ke hari.
# ------------------------------------------------------------

fit_nme_censored <- function(time, status) {
  time_scale <- median(time)
  x <- time / time_scale
  event <- status == 1

  neg_loglik <- function(log_par) {
    par <- exp(log_par)
    log_f <- dnmeweibull(x[event], par[1], par[2], par[3], log = TRUE)
    log_s <- snmeweibull(x[!event], par[1], par[2], par[3], log = TRUE)
    # Koreksi densitas ke satuan hari; survival tidak perlu koreksi.
    ll <- c(log_f - log(time_scale), log_s)
    if (any(!is.finite(ll))) return(1e100)
    -sum(ll)
  }

  # Beberapa nilai awal sederhana untuk mengurangi ketergantungan start.
  starts <- expand.grid(alpha = c(0.01, 0.1, 1, 10),
                        beta = c(0.7, 1.5, 3))
  lower <- c(alpha = 1e-6, beta = 0.05, lambda = 1e-10)
  upper <- c(alpha = 1e4, beta = 20, lambda = 100)
  best_fit <- NULL

  for (i in seq_len(nrow(starts))) {
    a0 <- starts$alpha[i]
    start <- c(a0, starts$beta[i], a0 / (a0 + 2))
    fit <- try(stats::optim(
      par = log(start), fn = neg_loglik, method = "L-BFGS-B",
      lower = log(lower), upper = log(upper),
      control = list(maxit = 5000, factr = 1e5)
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    if (fit$convergence != 0 || !is.finite(fit$value) || fit$value >= 1e99) next
    if (is.null(best_fit) || fit$value < best_fit$value) best_fit <- fit
  }
  if (is.null(best_fit)) stop("Estimasi NME-Weibull tidak konvergen.")

  estimate <- setNames(exp(best_fit$par), c("alpha", "beta", "lambda"))
  # Batas berlaku pada skala waktu internal, sebelum konversi lambda.
  if (any(abs(best_fit$par - log(lower)) < 0.01) ||
      any(abs(best_fit$par - log(upper)) < 0.01)) {
    warning("Estimasi mendekati batas parameter; interpretasikan dengan hati-hati.")
  }
  estimate["lambda"] <- estimate["lambda"] / time_scale^estimate["beta"]
  list(estimate = estimate, logLik = -best_fit$value,
       AIC = 6 + 2 * best_fit$value,
       BIC = 3 * log(length(time)) + 2 * best_fit$value)
}

fit_nme <- fit_nme_censored(time, status)
alpha_hat <- unname(fit_nme$estimate["alpha"])
beta_hat <- unname(fit_nme$estimate["beta"])
lambda_hat <- unname(fit_nme$estimate["lambda"])

# ------------------------------------------------------------
# 4. Estimasi Weibull dan perbandingan model
# survreg menangani sensor secara langsung, tanpa kovariat (~ 1).
# Parameterisasi pembanding: S(t) = exp(-lambda * t^beta).
# Pada survreg: beta = 1/scale; lambda = exp(-intercept * beta).
# ------------------------------------------------------------

fit_weibull <- survival::survreg(
  survival::Surv(time, status) ~ 1, dist = "weibull"
)
weibull_beta <- 1 / fit_weibull$scale
weibull_scale <- exp(unname(stats::coef(fit_weibull)[1]))
weibull_lambda <- weibull_scale^(-weibull_beta)
model_comparison <- data.frame(
  model = c("NME-Weibull", "Weibull"),
  alpha = c(alpha_hat, NA_real_),
  beta = c(beta_hat, weibull_beta),
  lambda = c(lambda_hat, weibull_lambda),
  logLik = c(fit_nme$logLik, as.numeric(stats::logLik(fit_weibull))),
  AIC = c(fit_nme$AIC, stats::AIC(fit_weibull)),
  BIC = c(fit_nme$BIC, stats::BIC(fit_weibull))
)
write.csv(model_comparison, file.path(result_dir, "perbandingan_model.csv"),
          row.names = FALSE)
print(model_comparison)

# ------------------------------------------------------------
# 5. Survival dan hazard pada waktu tertentu
# Seperti p = 0.90, 0.95, 0.99 pada evaluasi ukuran risiko,
# di sini titik evaluasi adalah 1, 2, 3, 4, dan 5 tahun (dalam hari).
# S(t) adalah peluang; h(t) adalah laju per hari, bukan peluang.
# ------------------------------------------------------------

evaluation_days <- c(365, 730, 1095, 1460, 1825)
# Tidak melakukan ekstrapolasi di luar waktu follow-up maksimum.
evaluation_days <- evaluation_days[evaluation_days <= max(time)]
if (length(evaluation_days) == 0) {
  stop("Ubah evaluation_days agar berada dalam rentang follow-up data.")
}
evaluation_table <- data.frame(
  time_days = evaluation_days,
  n_risk = vapply(evaluation_days, function(t) sum(time >= t), numeric(1)),
  S_NME = snmeweibull(evaluation_days, alpha_hat, beta_hat, lambda_hat),
  S_Weibull = stats::pweibull(evaluation_days, shape = weibull_beta,
                            scale = weibull_scale, lower.tail = FALSE),
  h_NME_per_day = hnmeweibull(evaluation_days, alpha_hat, beta_hat, lambda_hat),
  h_Weibull_per_day = weibull_beta * weibull_lambda *
    evaluation_days^(weibull_beta - 1)
)
write.csv(evaluation_table, file.path(result_dir, "survival_hazard.csv"),
          row.names = FALSE)
print(evaluation_table)

# ------------------------------------------------------------
# 6. Grafik survival dan hazard dari kedua model
# S(t): peluang masih hidup tanpa kekambuhan melewati waktu t.
# h(t): laju kekambuhan atau kematian, bersyarat masih bebas kejadian.
# ------------------------------------------------------------

grid_time <- seq(min(time), max(time), length.out = 300)
S_nme <- snmeweibull(grid_time, alpha_hat, beta_hat, lambda_hat)
S_weibull <- stats::pweibull(grid_time, shape = weibull_beta,
                           scale = weibull_scale, lower.tail = FALSE)
h_nme <- hnmeweibull(grid_time, alpha_hat, beta_hat, lambda_hat)
h_weibull <- weibull_beta * weibull_lambda * grid_time^(weibull_beta - 1)

png(file.path(figure_dir, "survival.png"), width = 1000, height = 700, res = 150)
plot(c(0, grid_time), c(1, S_nme), type = "l", col = "blue", lwd = 2,
     xlab = "Waktu (hari)", ylab = "S(t)", ylim = c(0, 1),
     main = "Peluang masih hidup tanpa kekambuhan")
lines(c(0, grid_time), c(1, S_weibull), col = "red", lwd = 2, lty = 2)
legend("topright", c("NME-Weibull", "Weibull"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
dev.off()

# Mulai dari waktu positif karena hazard di t = 0 dapat tak hingga.
png(file.path(figure_dir, "hazard.png"), width = 1000, height = 700, res = 150)
plot(grid_time, h_nme, type = "l", col = "blue", lwd = 2,
     ylim = range(c(0, h_nme, h_weibull)),
     xlab = "Waktu (hari)", ylab = "h(t) per hari",
     main = "Hazard kekambuhan atau kematian")
lines(grid_time, h_weibull, col = "red", lwd = 2, lty = 2)
legend("bottomright", c("NME-Weibull", "Weibull"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
dev.off()

cat("\nSelesai. Tabel:", result_dir, "\nGrafik:", figure_dir, "\n")
# Cara membaca hasil:
# - AIC/BIC lebih kecil menunjukkan kecocokan relatif yang lebih baik.
# - Grafik membandingkan pola survival dan hazard kedua model.
# - Grafik hazard hanya menggambarkan bentuk hazard model terestimasi.
# - Ekor kurva dengan sedikit subjek berisiko harus dibaca hati-hati.
# - Evaluasi pada data fitting ini bukan validasi prediksi di data baru.
# Rujukan API:
# https://stat.ethz.ch/R-manual/R-patched/library/survival/html/survreg.html
