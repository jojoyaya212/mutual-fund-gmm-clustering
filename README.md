# Bayesian Gaussian Mixture Model for Mutual Fund Clustering

## 📌 Overview
This project applies a **Bayesian Gaussian Mixture Model (GMM)** with Gibbs sampling to cluster mutual funds based on their monthly returns. Using MCMC methods, each fund is assigned to a latent state representing return-risk dynamics. The analysis identifies structural breaks and crisis periods in financial markets【72†source】.

---

## 📊 Key Findings
- Funds are clustered into **latent states** that capture return and volatility patterns.  
- Posterior inference provides credible intervals for means, variances, and mixture weights.  
- **Crisis detection**: State 2 aligns with major crises such as:  
  - Global Financial Crisis (2008–09)  
  - European Debt Crisis (2011)  
  - Oil Price Crash & China Slowdown (2016)  
  - Fed Tightening & Trade War (2018)  
  - COVID-19 Crash (2020)  
  - Inflation Surge & Ukraine War (2022)【72†source】.  
- Posterior predictive densities validate clustering robustness.  

---

## 🧪 Methods
1. **Data**  
   - Input: Mutual fund monthly returns (`MF_monthly_returns_full_len.csv`).  
   - Subset: First 150 funds used for clustering.  

2. **Model**  
   - Gaussian mixture with **K=2 latent states** (extendable).  
   - Priors: Normal-Inverse-Gamma for (μ, σ²), Dirichlet for mixture weights.  

3. **Inference**  
   - Gibbs sampling (10,000 iterations, 5,000 burn-in).  
   - Posterior means, variances, and state probabilities computed.  
   - Fund-level assignments by posterior mode and probability.  

4. **Outputs**  
   - Trace plots and density plots for μ, σ², and π.  
   - Posterior predictive density for new fund returns.  
   - Fund-to-state mapping with posterior probabilities.  
   - Visualizations of fund state distributions across crises.  

---

## 📂 Repository Structure
- **fund_gmm_model.R** – R script implementing Bayesian GMM via Gibbs sampling.  
- **GMM_output.pdf** – Report with parameter estimates, plots, and crisis mapping.  

---

## ⚙️ Requirements
This project uses **R** with the following packages:
```r
install.packages(c("MCMCpack"))
```

---

## 🚀 Usage
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/mutual-fund-gmm-clustering.git
   cd mutual-fund-gmm-clustering
   ```
2. Open **fund_gmm_model.R** in RStudio or R.  
3. Run the script to perform Gibbs sampling and clustering.  
4. Review results:  
   - Fund-level assignments (`assign_df`)  
   - Crisis alignment in **GMM_output.pdf**  

---

## 📑 References
- Bayesian Gaussian Mixture Models – MCMCpack R Package  
- Crisis periods identified from financial history (2008–2023)【72†source】  
