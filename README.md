# Bayesian Gaussian Mixture Model for Mutual Fund Clustering (10 Funds)

## 📌 Overview
This project implements a **Bayesian Gaussian Mixture Model (GMM)** with Gibbs sampling to analyze the return dynamics of the first 10 mutual funds. Each monthly return is modeled as coming from one of two latent states, capturing different mean–variance regimes. Posterior inference provides parameter distributions, state assignments, and predictive densities.

---

## 📊 Key Findings
- **Two-state clustering** effectively separates high-return/low-volatility and low-return/high-volatility regimes.  
- Posterior summaries show clear separation of state means:  
  - μ₁ ≈ 0.8% (σ² ≈ 0.00094)  
  - μ₂ ≈ -2.4% (σ² ≈ 0.0089)  
- Mixture weights: π₁ ≈ 0.77, π₂ ≈ 0.23, indicating most returns fall in the stable state.  
- State 2 strongly aligns with **financial crises**, including:  
  - Global Financial Crisis (2008–09)  
  - European Debt Crisis (2011)  
  - Oil price crash & China slowdown (2016)  
  - Fed tightening & Trade War (2018)  
  - COVID-19 crash (2020)  
  - Inflation surge & Ukraine war (2022).  
- Posterior predictive density validates robustness of the clustering.  

---

## 🧪 Methods
1. **Data**  
   - Input: `MF_monthly_returns_full_len.csv` (first 10 funds).  
   - Converted to pooled vector of monthly returns.  

2. **Model**  
   - Gaussian mixture with **K=2 latent states**.  
   - Priors: Normal-Inverse-Gamma for (μ, σ²); Dirichlet for mixture weights.  

3. **Inference**  
   - Gibbs sampler with 10,000 iterations (5,000 burn-in).  
   - Posterior draws for μ, σ², π stored for analysis.  
   - State assignment per return via posterior probabilities.  

4. **Outputs**  
   - Trace plots and posterior densities of μ, σ², and π.  
   - Predictive density for new returns.  
   - State assignment matrix reshaped to fund × time, exported to Excel (`fund_state_assignments.xlsx`).  
   - Crisis alignment visualization in **GMM output.pdf**.  

---

## 📂 Repository Structure
- **gmm_test3.R** – R script implementing Bayesian GMM for first 10 funds.  
- **GMM_output.pdf** – Output file with parameter estimates, posterior plots, and crisis mapping.  
- **fund_state_assignments.xlsx** – Exported Excel file with fund × time state assignments.  

---

## ⚙️ Requirements
This project uses **R** with the following packages:
```r
install.packages(c("MCMCpack", "writexl"))
```

---

## 🚀 Usage
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/mutual-fund-gmm-10funds.git
   cd mutual-fund-gmm-10funds
   ```
2. Open **fund_gmm_10funds.R** in RStudio or R.  
3. Run the script to:  
   - Perform Gibbs sampling and posterior inference.  
   - Generate posterior summaries and plots.  
   - Export fund × time state assignments to Excel.  
4. Review results in **GMM_output.pdf**.  

---

## 📑 References
- Bayesian Gaussian Mixture Models – *MCMCpack* R Package  
- Historical financial crises (2008–2023)  
