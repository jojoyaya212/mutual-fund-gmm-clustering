# --- packages ---
install.packages("MCMCpack")
library(MCMCpack)

##### FUND MODEL (150 funds clustered into 3 states) #####

# -------------------------------------------------------------------
# DATA: use FIRST 150 FUNDS (first 150 columns) as individual items
# Each column is a fund; rows are monthly returns
# -------------------------------------------------------------------
mf_data <- read.csv("MF_monthly_returns_full_len.csv")

# Keep only the first 150 funds (columns). Coerce to numeric matrix.
fund_returns <- as.matrix(mf_data[, 1:250, drop = FALSE])
storage.mode(fund_returns) <- "double"

# Number of time periods and number of funds
T_time  <- nrow(fund_returns)  # months per fund
T_fund  <- ncol(fund_returns)  # number of funds (150)

# For initialization and grids we’ll still use a pooled vector of returns,
# but NOTE: s_fund is now a length-150 vector assigning EACH FUND to {1,2,3}.
y_fund <- as.numeric(fund_returns)  # pooled returns for init & plotting
K <- 2

# --- priors (unchanged) ---
m0 <- 0
v0 <- 10^2
v0_inv <- 1 / v0
alpha <- rep(1, K)
v_prior <- 0.01
s_prior <- 0.01

# --- initialization (uses pooled returns for rough scale/location) ---
set.seed(123)
mu_fund <- rnorm(K, mean(y_fund, na.rm = TRUE), sd(y_fund, na.rm = TRUE))
sigma2_fund <- rep(var(y_fund, na.rm = TRUE), K)
pi_k_fund <- rep(1/K, K)

# IMPORTANT CHANGE: s_fund now assigns EACH FUND to a state (length = T_fund = 150)
s_fund <- sample(1:K, T_fund, replace = TRUE)

# --- MCMC store (unchanged shapes) ---
n_iter <- 10000
burn_in <- 5000
mu_store_fund <- matrix(NA, n_iter, K)
sigma2_store_fund <- matrix(NA, n_iter, K)
pi_store_fund <- matrix(NA, n_iter, K)
s_store_fund <- matrix(NA, n_iter, T_fund)

# --- helper to get all returns belonging to component j ---
# Given current s_fund (fund -> state), collect all returns from funds in state j
gather_component_returns <- function(j, fund_returns, s_vec) {
  if (!any(s_vec == j)) return(numeric(0))
  # cbind the columns whose funds belong to state j, then vectorize
  as.numeric(fund_returns[, s_vec == j, drop = FALSE])
}

# --- Gibbs sampler ---
for (iter in 1:n_iter) {
  
  # --- 1) Update (mu, sigma2) for each component using funds currently in that component ---
  for (j in 1:K) {
    y_j <- gather_component_returns(j, fund_returns, s_fund)
    n_j <- length(y_j)
    
    # Conjugate Normal-Inverse-Gamma updates (unchanged form)
    v_post <- 1 / (n_j / sigma2_fund[j] + v0_inv)
    m_post <- if (n_j > 0) v_post * (sum(y_j) / sigma2_fund[j] + m0 * v0_inv) else m0
    mu_fund[j] <- rnorm(1, m_post, sqrt(v_post))
    
    shape <- v_prior / 2 + n_j / 2
    rate  <- s_prior / 2 + if (n_j > 0) sum((y_j - mu_fund[j])^2) / 2 else 0
    sigma2_fund[j] <- 1 / rgamma(1, shape = shape, rate = max(rate, 1e-8))
  }
  
  # --- 2) Update s_fund (state of EACH FUND) using its entire return series ---
  # Use log-likelihood to avoid underflow: sum_t log N(y_it | mu_j, sigma2_j) + log pi_j
  for (t in 1:T_fund) {
    y_it <- fund_returns[, t]  # all time points for fund t
    ll <- sapply(1:K, function(j) {
      sum(dnorm(y_it, mean = mu_fund[j], sd = sqrt(sigma2_fund[j]), log = TRUE)) + log(pi_k_fund[j])
    })
    # stabilize & softmax
    ll <- ll - max(ll)
    probs <- exp(ll)
    probs <- probs / sum(probs)
    s_fund[t] <- sample(1:K, 1, prob = probs)
  }
  
  # --- 3) Update mixture weights pi ---
  n_k <- sapply(1:K, function(j) sum(s_fund == j))
  pi_k_fund <- as.numeric(rdirichlet(1, alpha + n_k))
  
  # --- store ---
  mu_store_fund[iter, ] <- mu_fund
  sigma2_store_fund[iter, ] <- sigma2_fund
  pi_store_fund[iter, ] <- pi_k_fund
  s_store_fund[iter, ] <- s_fund
}

# --- posterior means (unchanged) ---
mu_post_fund <- colMeans(mu_store_fund[(burn_in + 1):n_iter, ])
sigma2_post_fund <- colMeans(sigma2_store_fund[(burn_in + 1):n_iter, ])
pi_post_fund <- colMeans(pi_store_fund[(burn_in + 1):n_iter, ])
cat("FUND Posterior Means (mu):", mu_post_fund, "\n")
cat("FUND Posterior Sigmas (sigma):", sigma2_post_fund, "\n")
cat("FUND Posterior Pis (pi):", pi_post_fund, "\n")

# ----- FUND trace plot and density (unchanged) -----
par(mfrow = c(K, 1))
for (i in 1:K) {
  plot(mu_store_fund[,i], type = "l", col = "blue", ylab = paste0("mu", i, " (fund)"),
       main = paste0("Trace plot of mu", i, " (fund)"))
  abline(v = burn_in, col = "red", lty = 2)
}

mu_post_draws_fund <- mu_store_fund[(burn_in + 1):n_iter, ]
sigma2_post_draws_fund <- sigma2_store_fund[(burn_in + 1):n_iter, ]
pi_post_draws_fund <- pi_store_fund[(burn_in + 1):n_iter, ]

summary_stats <- function(draws) {
  t(apply(draws, 2, function(x) {
    c(
      Mean = mean(x),
      Median = median(x),
      St.Dev = sd(x),
      Lower = quantile(x, 0.025),
      Upper = quantile(x, 0.975)
    )
  }))
}

summary_df_fund <- rbind(
  summary_stats(mu_post_draws_fund),
  summary_stats(sigma2_post_draws_fund),
  summary_stats(pi_post_draws_fund)
)
rownames(summary_df_fund) <- c(paste0("mu", 1:K, "_fund"),
                               paste0("sigma2_", 1:K, "_fund"),
                               paste0("pi", 1:K, "_fund"))
summary_df_fund <- as.data.frame(summary_df_fund)
summary_df_fund$`95% Interval` <- paste0("(", round(summary_df_fund$Lower, 4), ", ",
                                         round(summary_df_fund$Upper, 4), ")")
summary_df_fund <- summary_df_fund[, c("Mean", "Median", "St.Dev", "95% Interval")]
print(summary_df_fund)

par(mfrow = c(K, 3))
for (i in 1:K) {
  plot(density(mu_post_draws_fund[,i]), main=paste0("Posterior of mu", i, " (fund)"),
       xlab=paste0("mu", i, "_fund"))
  plot(density(sigma2_post_draws_fund[,i]), main=paste0("Posterior of sigma2_", i, " (fund)"),
       xlab=paste0("sigma2_", i, "_fund"))
  plot(density(pi_post_draws_fund[,i]), main=paste0("Posterior of pi", i, " (fund)"),
       xlab=paste0("pi", i, "_fund"))
}

# --- posterior predictive density over returns (unchanged form) ---
# We construct predictive density for a future return y_{T+1} from the mixture
y_grid_fund <- seq(from = min(y_fund, na.rm = TRUE) - 0.05,
                   to   = max(y_fund, na.rm = TRUE) + 0.05, length.out = 200)
predictive_density_fund <- numeric(length(y_grid_fund))
M <- nrow(mu_post_draws_fund)
for (m in 1:M) {
  for (k in 1:K) {
    predictive_density_fund <- predictive_density_fund +
      pi_post_draws_fund[m, k] *
      dnorm(y_grid_fund, mu_post_draws_fund[m, k], sqrt(sigma2_post_draws_fund[m, k]))
  }
}
predictive_density_fund <- predictive_density_fund / M

# Set normal plot size and readable font (unchanged)
par(mfrow = c(1, 1))
par(mar = c(5, 5, 4, 2) + 0.1)
par(cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)

# Plot posterior predictive density for fund (unchanged)
plot(y_grid_fund, predictive_density_fund, type = "l", col = "blue", lwd = 2,
     main = expression("Posterior Predictive Density of " ~ y[T+1] ~ "(fund)"),
     xlab = expression(y[T+1]), ylab = "Density")

# --- OPTIONAL: a quick look at final fund-to-state assignment (kept separate, no deletion of your outputs) ---
cat("\nFinal fund-to-state assignment counts:\n")
print(table(s_fund))


# --- Fund-level state assignment details ---

# Use column names if available, otherwise create Fund_1 ... Fund_T
fund_names <- colnames(fund_returns)
if (is.null(fund_names)) fund_names <- paste0("Fund_", seq_len(T_fund))

# 1) Final iteration assignment (already in s_fund)
final_assign <- s_fund

# 2) Posterior mode assignment per fund (most frequent state after burn-in)
s_post <- s_store_fund[(burn_in + 1):n_iter, , drop = FALSE]  # draws × funds
mode_state <- apply(s_post, 2, function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
})

# 3) Posterior probability for each state (proportion of post–burn-in draws)
post_probs <- sapply(1:K, function(j) colMeans(s_post == j))
colnames(post_probs) <- paste0("P_state", 1:K)

# Combine into a tidy data.frame
assign_df <- data.frame(
  Fund = fund_names,
  FinalState = final_assign,
  ModeState = mode_state,
  post_probs,
  check.names = FALSE
)

# Show a quick preview
cat("\nHead of fund-to-state mapping (final vs posterior mode and probabilities):\n")
print(head(assign_df, 10))

# Optional: sort by ModeState for a clean look
assign_df_sorted <- assign_df[order(assign_df$ModeState, assign_df$FinalState), ]

# Optional: write to CSV for inspection
# write.csv(assign_df_sorted, file = "fund_state_assignments.csv", row.names = FALSE)
# cat("\nSaved mapping to fund_state_assignments.csv\n")

cat("\nFunds by posterior-mode state:\n")
for (j in 1:K) {
  cat(sprintf("\nState %d:\n", j))
  print(fund_names[which(mode_state == j)])
}

# --- Three-lane plot: fixed margins, visible counts (axis 4), highest→lowest ranking ---

use_mode <- TRUE
state_raw <- if (use_mode) mode_state else s_fund   # length T_fund, values in {1,2,3}

# Posterior μ summaries
mu_mean <- colMeans(mu_post_draws_fund)
mu_ci   <- t(apply(mu_post_draws_fund, 2, quantile, probs = c(0.025, 0.975)))

# Rank states by DESC μ (highest first)
ord <- order(mu_mean, decreasing = TRUE)         # e.g., c(2,3,1)
rank_map <- match(1:K, ord)                      # original state -> rank position 1..K
state_ranked <- rank_map[state_raw]              # each fund's rank lane (1..K)

# Reordered μ and CI by rank
mu_mean_r <- mu_mean[ord]
mu_ci_r   <- mu_ci[ord, , drop = FALSE]

# Colors per rank (top→bottom)
lane_cols <- c("#d95f02", "#1b9e77", "#7570b3")

x_pos <- seq_len(T_fund)

# Big margins; use oma for extra room; put y-label with mtext
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1,1), mar = c(5.5, 10.5, 4, 18), oma = c(0, 0, 0, 0), mgp = c(3, 0.7, 0))

# We want rank 1 at the TOP. Map rank r -> y position:
y_for_rank <- rev(seq_len(K))   # for K=3: c(3,2,1) so rank1→y=3, rank2→y=2, rank3→y=1

plot(NA,
     xlim = c(1, T_fund), ylim = c(0.5, K + 0.5),
     xlab = "Fund index", ylab = "",
     yaxt = "n", xaxs = "i",
     main = if (use_mode)
       "Funds by Posterior-Mode State (Ranked by Mean Return)"
     else
       "Funds by Final-Iteration State (Ranked by Mean Return)")

# Left y-axis labels: "Rank r (original S = ...)"
axis(2, at = y_for_rank,
     labels = paste0("Rank ", 1:K, "  (original S = ", ord, ")"),
     las = 1)

# Y-axis title placed with mtext so it's never clipped
#mtext("State (ranked by mean return, highest at top)", side = 2, line = 5.5)

# Baselines
abline(h = y_for_rank, lwd = 2)

# Fund ticks per ranked lane
tick_half <- 0.12
for (r in 1:K) {
  idx <- which(state_ranked == r)
  if (length(idx)) {
    y0 <- y_for_rank[r]
    segments(x0 = x_pos[idx], y0 = y0 - tick_half,
             x1 = x_pos[idx], y1 = y0 + tick_half, lwd = 1.1, col = lane_cols[r])
  }
}

# Light vertical guides
abline(v = seq(20, T_fund, by = 20), col = "gray85", lty = 3)

# Counts per lane
counts_by_rank <- sapply(1:K, function(r) sum(state_ranked == r))

# Draw counts on a RIGHT-HAND AXIS so they always show
axis(4, at = y_for_rank, labels = paste0("n = ", counts_by_rank), las = 1, tick = FALSE)

# (Optional, if you want colored counts instead of black, uncomment this block)
# usr <- par("usr")
# text(usr[2] - 2, y_for_rank, paste0("n = ", counts_by_rank),
#      adj = 0, cex = 0.95, col = lane_cols)

# μ annotations OUTSIDE the plot (as percentages)
usr <- par("usr")
gap  <- 0.09 * (usr[2] - usr[1])
xdot <- usr[2] + 2.05 * gap
xtxt <- usr[2] + 2.35 * gap

par(xpd = NA)  # allow drawing in right margin
for (r in 1:K) {
  y0 <- y_for_rank[r]
  points(xdot, y0, pch = 19, col = lane_cols[r])
  line1 <- sprintf("μ = %.2f%%", 100 * mu_mean_r[r])
  line2 <- sprintf("[%.2f%%, %.2f%%]", 100 * mu_ci_r[r,1], 100 * mu_ci_r[r,2])
  text(xtxt, y0 + 0.14, line1, adj = 0, cex = 0.95, col = lane_cols[r])
  text(xtxt, y0 - 0.14, line2, adj = 0, cex = 0.95, col = lane_cols[r])
}
par(xpd = FALSE)

# Console check: ranked order really is highest → lowest
summ_tab <- data.frame(
  Rank            = 1:K,
  OriginalState   = ord,
  Count           = counts_by_rank,
  Mu_Mean_pct     = round(100 * mu_mean_r, 3),
  Mu_Lower95_pct  = round(100 * mu_ci_r[,1], 3),
  Mu_Upper95_pct  = round(100 * mu_ci_r[,2], 3)
)
cat("\nRanked state summary (highest μ at top):\n")
print(summ_tab, row.names = FALSE)

par(old_par)

