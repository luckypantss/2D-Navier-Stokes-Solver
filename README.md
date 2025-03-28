# 🌀 2D Navier-Stokes Solver

This project implements a 2D fluid dynamics solver from scratch using the Finite Difference Method (FDM) and Explicit Euler integration.  
It includes turbulence modeling via the Smagorinsky LES model, pressure smoothing using the Artificial Compressibility Method (ACM), and parameter tuning via Bayesian Optimization.

---

## 🔧 Features

- Lid-driven cavity flow simulation
- Finite Difference discretization
- Artificial Compressibility Method (ACM) for pressure-velocity coupling
- Smagorinsky LES turbulence model
- Bayesian optimization for timestep and smoothing parameter
- GPU-ready option using MATLAB's `gpuArray`

---

## 📊 Results

> Add your images here after upload!

- ✅ Pressure contours over time
- ✅ 2D velocity field visualization
- ✅ Reynolds number evolution

---

## 🧠 Why it Matters

While not a machine learning project directly, this simulation demonstrates:
- Numerical modeling skills relevant to scientific ML
- Data generation pipelines for surrogate modeling or PINNs
- Optimization workflows applicable to hyperparameter tuning in ML

---

## 📘 Related Thesis

This project is based on my Bachelor's Thesis:  
📄 **“Enabling Interactive Fluid Simulations”** — [Download it here](YOUR_LINK_HERE)  
> (256+ downloads)

---

## 🖥️ How to Run

1. Open `src/main.m` in MATLAB
2. Adjust parameters (Reynolds number, grid size, etc.)
3. Run the script to generate quiver plots and pressure contours

---

## 📁 Project Structure

