Scripts are organized by **model** (London, Atlanta) and by **functionality** (simulation, analysis, visualization).

---

### 1. `utility_functions.r`
Core library of helper routines shared across all scripts:
- Implements **Random Dot Product Graph (RDPG)** sampling (`rdpg.sample`) and **Adjacency Spectral Embedding (ASE)**.
- Provides **Procrustes** (`procrustes2`), **distance estimators** (`getD`, `getDW1`, `getDW2`), and **Multidimensional Scaling (MDS)** operations (`doMDS`).
- Implements required calculations for both true shuddled/dMV/W1 distance of **London** and **Atlanta** models:
  - `trueLondondMV`, `trueW1squareLondon` for the London model.
  - `trueAtlantadmv`, `trueshuffleAtlantadmv` for the Atlanta model.
- Defines error functions such as `pairederrorinshuffling`, and shuffling/perturbation utilities.
- Serves as support functions for simulation scripts (`MCSim_London.R`, `MCSim_Atlanta.R`).
- multiple functions adopted from the other [paper](https://github.com/TianyiChen97/Euclidean-mirrors-and-first-order-changepoints-in-network-time-series).
---

### 2. `London_model_simulation_functions.R`
Functions for **London model simulations**:
- Experiments to estimate changepoint detection accuracy under different perturbations.
- Functions include `simulatenetworkchangepoint`, `analyzenetworkchangepoint`, and `analyzeshufflingerrors`.
- Computes **$l_{\infty}$ errors** across multiple network distance metrics: average degree, W₁, W₂, and dMV.
- Generates plots summarizing estimation accuracy across network sizes and shuffling ratios.

---

### 3. `MCSim_London.R`
Full Monte Carlo simulation for the **London Model**:
- Uses functions including embedding, MDS, ISOmap, shuffling to simulate London Model for changepoint localization task.
- Output simulation results for different metrics including MDS/ISOMAP.
- Simulation for different parameters (e.g. change probabilities `p`, `q`, number of nodes `n`).

---

### 4. `MCSim_Atlanta.R`
Full Monte Carlo simulation for the **Atlanta Model**:
- Same purpose as above but for Atlanta model incorporating parameters `p, q, m, $\delta$, $t^*$`.

---

### 5. `summary_London-tables-and-plots-_NEW.R`/`summary_Atlanta-tables-and-plots-_NEW.R`
- Input simulation results stored in 'simulation results folder...' and generate tables and plots (Table ?? and Figure ?? ) to summarize simulation results.
