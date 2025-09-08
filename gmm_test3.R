library(MCMCpack)

##### FUND MODEL (First 10 Funds Combined) #####
mf_data <- read.csv("MF_monthly_returns_full_len.csv")
mf_data <- mf_data[, 1:10]  # Use only the first 10 funds
y_fund <- as.numeric(unlist(mf_data))  # Flatten into one vector
y_fund <- y_fund[!is.na(y_fund)]       # Remove any NA values
T_fund <- length(y_fund)
K <- 2

m0 <- 0
v0 <- 10^2
v0_inv <- 1 / v0
alpha <- rep(1, K)
v_prior <- 0.01
s_prior <- 0.01

set.seed(123)
mu_fund <- rnorm(K, mean(y_fund), sd(y_fund))
sigma2_fund <- rep(var(y_fund), K)
pi_k_fund <- rep(1/K, K)
s_fund <- sample(1:K, T_fund, replace = TRUE)

n_iter <- 10000
burn_in <- 5000
mu_store_fund <- matrix(NA, n_iter, K)
sigma2_store_fund <- matrix(NA, n_iter, K)
pi_store_fund <- matrix(NA, n_iter, K)
s_store_fund <- matrix(NA, n_iter, T_fund)

for (iter in 1:n_iter) {
  for (j in 1:K) {
    idx <- which(s_fund == j)
    n_j <- length(idx)
    y_j <- y_fund[idx]
    v_post <- 1 / (n_j / sigma2_fund[j] + v0_inv)
    m_post <- v_post * (sum(y_j) / sigma2_fund[j] + m0 * v0_inv)
    mu_fund[j] <- rnorm(1, m_post, sqrt(v_post))
    
    shape <- v_prior / 2 + n_j / 2
    rate <- s_prior / 2 + sum((y_j - mu_fund[j])^2) / 2
    sigma2_fund[j] <- 1 / rgamma(1, shape = shape, rate = max(rate, 1e-8))
  }
  
  for (t in 1:T_fund) {
    probs <- sapply(1:K, function(j) {
      dnorm(y_fund[t], mu_fund[j], sqrt(sigma2_fund[j])) * pi_k_fund[j]
    })
    probs <- probs / sum(probs)
    s_fund[t] <- sample(1:K, 1, prob = probs)
  }
  
  n_k <- sapply(1:K, function(j) sum(s_fund == j))
  pi_k_fund <- as.numeric(rdirichlet(1, alpha + n_k))
  
  mu_store_fund[iter, ] <- mu_fund
  sigma2_store_fund[iter, ] <- sigma2_fund
  pi_store_fund[iter, ] <- pi_k_fund
  s_store_fund[iter, ] <- s_fund
}


mu_post_fund <- colMeans(mu_store_fund[(burn_in + 1):n_iter, ])
sigma2_post_fund <- colMeans(sigma2_store_fund[(burn_in + 1):n_iter, ])
pi_post_fund <- colMeans(pi_store_fund[(burn_in + 1):n_iter, ])
cat("FUND Posterior Means (mu):", mu_post_fund, "\n")
cat("FUND Posterior Sigmas (sigma):", sigma2_post_fund, "\n")
cat("FUND Posterior Pis (pi):", pi_post_fund, "\n")

# FUND trace plot and density
par(mfrow = c(2, 1))
plot(mu_store_fund[,1], type = "l", col = "blue", ylab = "mu1 (fund)", main = "Trace plot of mu1 (fund)")
abline(v = burn_in, col = "red", lty = 2)
plot(mu_store_fund[,2], type = "l", col = "blue", ylab = "mu2 (fund)", main = "Trace plot of mu2 (fund)")
abline(v = burn_in, col = "red", lty = 2)

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
rownames(summary_df_fund) <- c("mu1_fund", "mu2_fund", "sigma2_1_fund", "sigma2_2_fund", "pi1_fund", "pi2_fund")
summary_df_fund <- as.data.frame(summary_df_fund)
summary_df_fund$`95% Interval` <- paste0("(", round(summary_df_fund$Lower, 4), ", ", round(summary_df_fund$Upper, 4), ")")
summary_df_fund <- summary_df_fund[, c("Mean", "Median", "St.Dev", "95% Interval")]
print(summary_df_fund)

par(mfrow = c(2, 3))
plot(density(mu_post_draws_fund[,1]), main="Posterior of mu1 (fund)", xlab="mu1_fund")
plot(density(mu_post_draws_fund[,2]), main="Posterior of mu2 (fund)", xlab="mu2_fund")
plot(density(sigma2_post_draws_fund[,1]), main="Posterior of sigma2_1 (fund)", xlab="sigma2_1_fund")
plot(density(sigma2_post_draws_fund[,2]), main="Posterior of sigma2_2 (fund)", xlab="sigma2_2_fund")
plot(density(pi_post_draws_fund[,1]), main="Posterior of pi1 (fund)", xlab="pi1_fund")
plot(density(pi_post_draws_fund[,2]), main="Posterior of pi2 (fund)", xlab="pi2_fund")

y_grid_fund <- seq(from = min(y_fund) - 0.05, to = max(y_fund) + 0.05, length.out = 200)
predictive_density_fund <- numeric(length(y_grid_fund))
M <- nrow(mu_post_draws_fund)
for (m in 1:M) {
  for (k in 1:K) {
    predictive_density_fund <- predictive_density_fund +
      pi_post_draws_fund[m, k] * dnorm(y_grid_fund, mu_post_draws_fund[m, k], sqrt(sigma2_post_draws_fund[m, k]))
  }
}
predictive_density_fund <- predictive_density_fund / M

# Set normal plot size and readable font
par(mfrow = c(1, 1))                      # Single plot layout
par(mar = c(5, 5, 4, 2) + 0.1)            # Margins: bottom, left, top, right
par(cex.lab = 1.2, cex.axis = 1.1,        # Label and axis size
    cex.main = 1.3)                       # Main title size

# Plot posterior predictive density for fund
plot(y_grid_fund, predictive_density_fund, type = "l", col = "blue", lwd = 2,
     main = expression("Posterior Predictive Density of " ~ y[T+1] ~ "(fund)"),
     xlab = expression(y[T+1]), ylab = "Density")

predictive_mean <- mean(rowSums(pi_post_draws_fund * mu_post_draws_fund))
cat("Posterior predictive mean of y_{T+1} (fund):", predictive_mean, "\n")

###############################################################################################
# Use posterior means to compute state assignment
posterior_probs <- matrix(NA, nrow = T_fund, ncol = K)

for (k in 1:K) {
  posterior_probs[, k] <- pi_post_fund[k] * dnorm(y_fund, mean = mu_post_fund[k], sd = sqrt(sigma2_post_fund[k]))
}

# Normalize to sum to 1
posterior_probs <- posterior_probs / rowSums(posterior_probs)

# Assign each y_t to the most probable state
state_assignments <- apply(posterior_probs, 1, which.max)

# Optional: Reshape assignments back into original fund matrix structure
original_shape <- dim(mf_data)
state_matrix <- matrix(state_assignments, nrow = nrow(mf_data), ncol = ncol(mf_data))

# View a sample
head(state_matrix)
# Install package if not already installed
install.packages("writexl")
library(writexl)

# Step 1: Assign states as previously calculated
posterior_probs <- matrix(NA, nrow = T_fund, ncol = K)
for (k in 1:K) {
  posterior_probs[, k] <- pi_post_fund[k] * dnorm(y_fund, mean = mu_post_fund[k], sd = sqrt(sigma2_post_fund[k]))
}
posterior_probs <- posterior_probs / rowSums(posterior_probs)
state_assignments <- apply(posterior_probs, 1, which.max)

# Step 2: Reshape to original fund-by-time format
state_matrix <- matrix(state_assignments, nrow = nrow(mf_data), ncol = ncol(mf_data))
colnames(state_matrix) <- colnames(mf_data)
rownames(state_matrix) <- rownames(mf_data)

# Step 3: Export to Excel
write_xlsx(as.data.frame(state_matrix), "fund_state_assignments.xlsx")
