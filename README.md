# Physics of Behavior — DPG 2026 Tutorial

**Tutorial session at the [DPG Spring Meeting Dresden 2026](https://dresden26.dpg-tagungen.de/)**

Given by [Antonio Carlos Costa, PI INSERM, Paris Brain Institute](https://antonioccosta.github.io/)

---

## Overview

How can ideas from statistical mechanics and dynamical systems help us understand animal behavior? This tutorial introduces a physics-inspired, data-driven framework for extracting long-timescale structure from behavioral time series.

Starting from stochastic dynamics in a double-well potential and building up to real *C. elegans* data, participants will learn how **transfer operator methods** can identify metastable behavioral states, slow modes, and coarse-grained dynamics.

---

## 🚀 Run in Google Colab (no installation required)

The easiest way to follow along is to open the notebooks directly in Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AntonioCCosta/DPG_PhysicsBehavior_2026)

Each notebook automatically installs all dependencies and clones the repository on first run.

---

## 📓 Tutorial Notebooks

The tutorial is structured as four progressive notebooks:

### [Notebook 1 — Theory: Stochastic dynamics on a double well](https://colab.research.google.com/github/AntonioCCosta/DPG_PhysicsBehavior_2026/blob/main/notebooks/notebook_1_theory_dw.ipynb)

An introduction to the theoretical foundations. We model behavior as a stochastic dynamical system and use a 1D double-well potential as a minimal example. Topics covered:
- Simulating Langevin dynamics (Euler–Maruyama scheme) 
- Steady-state distributions from the Fokker-Planck equation
- Timescale separation and Kramers escape rates
- Semi-analytical computation of the **transfer operator** (Fokker-Planck discretized via the Chang-Cooper scheme) and its eigenspectrum
- Identifying metastable states from slow eigenvectors; optimizing the coherence of metastable sets

---

### [Notebook 2 — Data-driven transfer operator on a double well](https://colab.research.google.com/github/AntonioCCosta/DPG_PhysicsBehavior_2026/blob/main/notebooks/notebook_2_toy_dw.ipynb)

Moving from theory to data. We approximate the transfer operator directly from simulated trajectories using the **Ulam-Galerkin method**, without access to the underlying equations.

---

### [Notebook 3 — Extending to a chaotic system: the Lorenz attractor](https://colab.research.google.com/github/AntonioCCosta/DPG_PhysicsBehavior_2026/blob/main/notebooks/notebook_3_Lorenz.ipynb)

Generalizing the framework beyond equilibrium dynamics. We apply the same approach to the Lorenz system — a canonical example of chaos — and tackle the realistic scenario of **partial observations**. Topics covered:
- Integrating the Lorenz equations and identifying coherent (almost-invariant) sets via the transfer operator
- Locating the unstable periodic orbit
- Handling non-Markovian observations: restoring Markovianity via **Takens delay embedding**
- Optimizing the delay embedding dimension using entropy rate and predictive information

---

### [Notebook 4 — Real data: *C. elegans* foraging behavior](https://colab.research.google.com/github/AntonioCCosta/DPG_PhysicsBehavior_2026/blob/main/notebooks/notebook_4_worm.ipynb)

Applying the full pipeline to real behavioral recordings. We work with body posture time series of freely crawling *C. elegans* and recover a multiscale picture of behavior.

---

## 📁 Repository Structure

```
DPG_PhysicsBehavior_2026/
├── notebooks/          # Tutorial Jupyter notebooks (run in Colab or locally)
│   ├── notebook_1_theory_dw.ipynb
│   ├── notebook_2_toy_dw.ipynb
│   ├── notebook_3_Lorenz.ipynb
│   └── notebook_4_worm.ipynb
├── data/               # Pre-computed simulation results and worm behavioral data
├── utils/              # Analysis utilities
│   ├── operator_calculations.py   # Transfer operator estimation and eigendecomposition
│   ├── delay_embedding.py         # Takens delay embedding (trajectory matrix)
│   ├── partition_methods.py       # k-means state space partitioning
│   └── coarse_graining.py         # Optimal metastable partition from slow eigenvectors
├── requirements.txt
└── README.md
```

---

## 🖥️ Running Locally

**1. Clone the repository**
```bash
git clone https://github.com/AntonioCCosta/DPG_PhysicsBehavior_2026.git
cd DPG_PhysicsBehavior_2026
```

**2. Install dependencies**
```bash
pip install -r requirements.txt
```

**3. Launch Jupyter**
```bash
jupyter notebook
```

Python 3.8+ recommended.

---

## 📚 References

The tutorial covers methods and concepts from the following works:

**Transfer operators & metastable dynamics**
- Risken. "Fokker-planck equation." In The Fokker-Planck equation: methods of solution and applications, pp. 63-95. Berlin, Heidelberg: Springer Berlin Heidelberg, 1989.
- Hänggi, Talkner, Borkovec. "Reaction-rate theory: fifty years after Kramers." *Reviews of modern physics* 62, no. 2 (1990): 251.
- Chang, Cooper. "A practical difference scheme for Fokker-Planck equations." *Journal of Computational Physics* 6, no. 1 (1970): 1-16.
- Froyland, Dellnitz. "Detecting and locating near-optimal almost-invariant sets and cycles." *SIAM Journal on Scientific Computing* 24, no. 6 (2003): 1839-1863.
- Froyland, Padberg. "Almost-invariant sets and invariant manifolds—connecting probabilistic and geometric descriptions of coherent structures in flows." *Physica D: Nonlinear Phenomena* 238, no. 16 (2009): 1507-1523.
- Bollt, Santitissadeekorn, eds. "Applied and computational measurable dynamics". Society for Industrial and Applied Mathematics, 2013.

**Markov state models & data-driven approaches**
- Bowman, Pande, Noé, eds. "An introduction to Markov state models and their application to long timescale molecular simulation". Springer Science & Business Media, 2013.
- Coifman, Lafon. "Diffusion maps." *Applied and computational harmonic analysis* 21, no. 1 (2006): 5-30.

**Delay embedding & predictability**
- Sauer, Yorke, Casdagli. "Embedology." *Journal of statistical Physics*, 65(3) (1991) pp.579-616.
- Gaspard, Wang. "Noise, chaos, and (ε, τ)-entropy per unit time." *Physics reports* 235, no. 6 (1993): 291-343.

**Unstable periodic orbits**
- Barrio, Dena, Tucker. "A database of rigorous and high-precision periodic orbits of the Lorenz model." *Computer Physics Communications* 194 (2015): 76-83.

**C. elegans behavioral analysis**
- Stephens, Johnson-Kerner, Bialek, Ryu. "Dimensionality and dynamics in the behavior of *C. elegans*". *PLoS Computational Biology*, 4(4), e1000028 (2008).
- Costa A.C., Ahamed, Jordan, Stephens. "Maximally predictive states: from partial observations to long timescales" *Chaos* 2023
- Costa A.C., Ahamed, Jordan, Stephens. "A Markovian dynamics for Caenorhabditis elegans behavior across scales." *PNAS* 121, no. 32 (2024): e2318805121.

**Further applications of these techniques
- Costa A.C., Sridhar, Wyart, Vergassola. "Fluctuating landscapes and heavy tails in animal behavior." *PRX life* 2, no. 2 (2024): 023001.
- Sridhar, Vergassola, Marques, Orger, Costa A.C., Wyart. "Uncovering multiscale structure in the variability of larval zebrafish navigation." *PNAS* 121, no. 47 (2024): e2410254121.


---

## 📬 Contact

Questions? Reach out:

- 📧 antonioccosta.phys@gmail.com
- 🌐 [antonioccosta.github.io](https://antonioccosta.github.io/)

---

## 📄 License

Released under [CC0 1.0 Universal](LICENSE) (Public Domain). Free to use, adapt, and redistribute without restriction.
