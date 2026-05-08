# Dupire Equation and Local Volatility Calibration (CEV & Gatheral)
 
## Overview
 
This project implements the **numerical resolution of the Dupire equation** using the **Crank-Nicolson implicit finite difference scheme**, and calibrates **local volatility** models using the **Levenberg-Marquardt algorithm**.
 
Two local volatility parametrizations are studied:
- The **CEV model** (Constant Elasticity of Variance): $\sigma_{locale}(K) = \beta_1 / K^{\beta_2}$
- The **Gatheral model**: $\sigma_{locale}(K) = b\left(\rho(K - m) + \sqrt{(K-m)^2 + a^2}\right)$
The implementation is written in **Python** and relies on tridiagonal linear system solvers and iterative gradient-based optimization.
 
---
 
## Project Structure
 
```
TP4/
├── TP4.py                    # Python implementation of all computations
└── README.md                 # Project documentation
```
 
---
 
# Theoretical Background
 
## Dupire Equation
 
The **Dupire equation** gives the price of a European option $V(T, K)$ as a function of strike $K$ and maturity $T$, for a fixed initial asset price $S_0$ and time $t_0$:
 
$$\frac{\partial V}{\partial T} + rK\frac{\partial V}{\partial K} - \frac{1}{2}\sigma^2_{locale}(K,T) K^2 \frac{\partial^2 V}{\partial K^2} = 0$$
 
with boundary and initial conditions:
 
$$V(T=0, K) = \max(S_0 - K, 0)$$
$$V(T, K=0) = S_0, \qquad V(T, K=K_{max}) = 0$$
 
---
 
## Crank-Nicolson Discretization
 
The spatial and temporal grids are:
 
$$K_i \in [0, K_{max}], \quad \Delta K = \frac{K_{max}}{N+1}, \qquad T^n \in [0, T_{max}], \quad \Delta t = \frac{T_{max}}{M+1}$$
 
The Crank-Nicolson scheme leads to the tridiagonal system at each time step:
 
$$B_i V^{n+1}_{i-1} + D_i V^{n+1}_i + A_i V^{n+1}_{i+1} = C^n_i$$
 
with coefficients:
 
$$A_i = \frac{\Delta t}{4}\left(\frac{K_i}{\Delta K}r - \sigma^2_{loc}\left(\frac{K_i}{\Delta K}\right)^2\right)$$
 
$$B_i = -\frac{\Delta t}{4}\left(\frac{K_i}{\Delta K}r + \sigma^2_{loc}\left(\frac{K_i}{\Delta K}\right)^2\right)$$
 
$$D_i = 1 + \frac{\Delta t}{2}\sigma^2_{loc}\left(\frac{K_i}{\Delta K}\right)^2$$
 
At each time step the system is solved via **Thomas algorithm** (LU decomposition of the tridiagonal matrix).
 
---
 
## Numerical Parameters
 
| Parameter | Value |
|-----------|-------|
| $K_{max}$ | 20 |
| $S_0$ | 10 |
| $r$ | 0.1 |
| $T_{max}$ | 0.5 |
| $N$ | 199 |
| $M$ | 49 |
 
---
 
# Implemented Functions
 
### Dupire Pricing
- `sigma_locale(h, i, beta1, beta2)` — Computes the local volatility at grid node $i$:
$$\sigma(K_i) = \beta_1 / K_i^{\beta_2} + h$$
- `Prix_Dupire(h, beta1, beta2)` — Solves the Dupire PDE via Crank-Nicolson and returns the full price matrix $V^n_i$.
---
 
### Vega
- `Vega_Dupire(h, beta1, beta2)` — Computes the numerical Vega via finite difference:
$$\text{Vega}(T, K, \sigma_l) = \frac{V(T, K, \sigma_l + h) - V(T, K, \sigma_l)}{h}, \quad h = 0.01$$
---
 
### Calibration Utilities
- `Prix_Dupire_Utiles(beta1, beta2)` — Extracts option prices at the market strikes $K_p$ from the full matrix, using:
$$i_p = \frac{K_p}{\Delta K} + 1$$
- `Vega_Utiles(beta1, beta2)` — Extracts Vega values at market strikes $K_p$.
---
 
### Levenberg-Marquardt Calibration
- `levenberg_marquardt_CEV(V_market, K_market)` — Calibrates $(\beta_1, \beta_2)$ in the CEV model by minimizing:
$$I(\beta_1, \beta_2) = \sum_{p=1}^{P} \omega_p \left|V^{marche}(T_p, K_p) - V^{dupire}(T_p, K_p, \sigma_{locale}(\beta_1, \beta_2))\right|^2$$
- `levenberg_marquardt_Gatheral(V_market, K_market)` — Calibrates $(a, m)$ in the Gatheral model for fixed $b = 0.05$, $\rho = 0.1$.
---
 
# Part I — Dupire Equation: Numerical Resolution
 
## Two Cases for Local Volatility
 
**Case a) — Constant volatility:**
 
$$\sigma_l = 0.3$$
 
**Case b) — CEV parametrization:**
 
$$\sigma(K_i) = \frac{\beta_1}{K_i^{\beta_2}}, \quad \beta_1 = 1, \quad \beta_2 = 1$$
 
## Visualizations
 
- 2D plots of $V(T, K)$ for maturities $T$, $T/2$ and $T=0$
- 3D surface $V(T, K)$
- 2D and 3D Vega surfaces for $\sigma_l = 0.3$ and $\sigma_l = 1/K$
---
 
# Part II — CEV Local Volatility Calibration
 
## Model
 
$$\sigma_{locale}(S, t) = \frac{\beta_1}{S^{\beta_2}}$$
 
## Market Data ($T = 0.5$)
 
| $K_p$ | 7 | 7.5 | 8 | 8.5 | 9 | 9.5 | 10 | 10.5 | 11 | 11.5 | 12 | 12.5 | 13 | 13.5 | 14 |
|--------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| $V^{marche}_p$ | 3.3634 | 2.9092 | 2.4703 | 2.0536 | 1.6666 | 1.3167 | 1.0100 | 0.7504 | 0.5389 | 0.3733 | 0.2491 | 0.1599 | 0.0986 | 0.0584 | 0.0332 |
 
## Jacobian
 
The Jacobian matrix $J \in \mathbb{R}^{P \times 2}$ has rows:
 
$$\frac{\partial r_p}{\partial \beta_1} = -\text{Vega}(K_p, T) \cdot \frac{1}{K_p^{\beta_2}}$$
 
$$\frac{\partial r_p}{\partial \beta_2} = \text{Vega}(K_p, T) \cdot \frac{\beta_1 \ln(K_p)}{K_p^{\beta_2}}$$
 
## Algorithm
 
**Initialization:** $\beta_1^{(0)} = 1$, $\beta_2^{(0)} = 1$, $\varepsilon = 10^{-5}$, $\lambda = 0.001$
 
**Update rule at iteration $k$:**
 
$$M = J^T J + \lambda I$$
$$d^{(k)} = -M^{-1} J^T r^{(k)}$$
$$\beta^{(k+1)} = \beta^{(k)} + d^{(k)}$$
 
**Stopping criterion:** $\sqrt{(d_1^{(k)})^2 + (d_2^{(k)})^2} < \varepsilon$
 
---
 
# Part III — Gatheral Model Calibration
 
## Model
 
$$\sigma_{locale}(K) = b\left(\rho(K - m) + \sqrt{(K - m)^2 + a^2}\right)$$
 
with fixed $b = 0.05$, $\rho = 0.1$, and calibrated parameters $a = \beta_1$, $m = \beta_2$.
 
## Market Data ($T = 0.5$)
 
| $K_p$ | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|--------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| $V^{marche}_p$ | 5.2705 | 4.3783 | 3.5510 | 2.8138 | 2.1833 | 1.6651 | 1.2541 | 0.9374 | 0.6983 | 0.5195 | 0.3851 | 0.2817 | 0.1987 | 0.1277 |
 
## Visualizations
 
- 3D Dupire surface for initial parameters $a = 5$, $m = 5$
- 3D Vega surface for $a = 5$, $m = 5$
- 2D plot of calibrated local volatility $K \mapsto \sigma_{locale}(K)$
- 3D Dupire surface after calibration
- 3D Vega surface after calibration
---
 
# Libraries Used
 
```
numpy
matplotlib
scipy
```
 
Install dependencies with:
 
```bash
pip install numpy matplotlib scipy
```
 
---
 
# Usage
 
Run the main script:
 
```bash
python TP4.py
```
 
The script will:
1. Solve the Dupire PDE via Crank-Nicolson for constant and CEV volatility (Part I)
2. Visualize option price surfaces in 2D and 3D
3. Compute and visualize Vega surfaces
4. Calibrate CEV parameters $(\beta_1, \beta_2)$ via Levenberg-Marquardt (Part II)
5. Calibrate Gatheral parameters $(a, m)$ via Levenberg-Marquardt (Part III)
6. Plot calibrated local volatility curves and surfaces
---
 
# Key Concepts Covered
 
- Dupire forward PDE for option pricing
- Crank-Nicolson implicit finite difference scheme
- Tridiagonal system resolution (Thomas algorithm)
- Local volatility models (CEV, Gatheral)
- Vega computation by finite difference
- Levenberg-Marquardt nonlinear least squares
- Jacobian computation for gradient-based calibration
- Ill-posed inverse problems and Tikhonov regularization
---
 
# Disclaimer
 
This repository contains the **code and methodology** used in an academic laboratory session at CY Tech (Université Paris-Cergy). The original written report has been removed due to **academic and data usage restrictions**.
 
