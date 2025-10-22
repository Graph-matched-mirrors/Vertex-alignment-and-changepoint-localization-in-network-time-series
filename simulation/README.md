Scripts are organized by **model** (London, Atlanta) and by **functionality** (simulation, analysis, visualization).  Many functions (mostly in utility_functions.r) builds upon tools and methods from the **Euclidean Mirror Graph Matching** project, available at:

[TianyiChen97/Euclidean_mirror_graph_matching](https://github.com/TianyiChen97/Euclidean-mirrors-and-first-order-changepoints-in-network-time-series)

---

### 1. `utility_functions.r`
Core library of helper routines shared across all scripts:
- Implements **Random Dot Product Graph (RDPG)** sampling (`rdpg.sample`) and **Adjacency Spectral Embedding (ASE)**.
- Provides **Procrustes** (`procrustes2`), **distance estimators** (`getD`, `getDW1`, `getDW2`), and **Multidimensional Scaling (MDS)** operations (`doMDS`).
- Implements required calculations for both **London** and **Atlanta** models:
  - `trueLondondMV`, `trueW1squareLondon` for the London changepoint model.
  - `trueAtlantadmv`, `trueshuffleAtlantadmv` for the Atlanta latent process model.
- Defines error functions such as `pairederrorinshuffling`, and shuffling/perturbation utilities.
- Serves as support functions for simulation scripts (`MCSim_London.R`, `MCSim_Atlanta.R`).

---

### 2. `London_model_simulation_functions.R`
Functions for **London model simulations**:
- Experiments to estimate changepoint detection accuracy under different perturbations.
- Functions include `simulatenetworkchangepoint`, `analyzenetworkchangepoint`, and `analyzeshufflingerrors`.
- Computes **L∞ errors** across multiple network distance metrics: average degree, W₁, W₂, and dMV.
- Generates plots summarizing estimation accuracy across network sizes and shuffling ratios.

---

### 3. `MCSim_London.R`
Full Monte Carlo simulation for the **London Model**:
- Uses functions including embedding, MDS, ISOmap, shuffling to simulate London Model
- Produces analyses and visualisation for metrics including MDS/ISOMAP.
- Comprehensive simulation for different parameters (e.g. change probabilities `p`, `q`, number of nodes `n`).

---

### 4. `MCSim_Atlanta.R`
Full Monte Carlo simulation for the **Atlanta Model**:
- Same purpose as above but for Atlanta model incorporating parameters `(p, q, m, δ, t*)`.
