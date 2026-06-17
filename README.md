# ChronoCCI 

**ChronoCCI** is a specialized single-cell bioinformatics tool designed to align cell-cell communication (CCI) networks along branching cell-fate trajectories. By combining the Markov chain destiny mapping of **CellRank** with the consensus ligand-receptor databases of **LIANA**, it reconstructs the chronological cascade of signaling waves and automatically prioritizes molecular gatekeepers governing lineage choices.

## Features
* **Automated Trajectory Alignment:** Integrates with CellRank's `GPCCA` to track signals across continuous cell lineages.
* **Cell-Type Relay Mapping:** Tracks the donor-to-acceptor communication channels active at the exact peak of each molecular wave.
* **Bifurcation Screening:** Utilizes a chronological Mean Absolute Error (MAE) scoring formula to scan for asymmetric signaling patterns between diverging paths.

## Installation
You can install the package directly from GitHub using `pip`:

```bash
pip install git+https://github.com/ronnaug/chronocci
```

## Quick Start
Run a comprehensive, end-to-end joint trajectory and interaction analysis with just a few lines of code:

```python
import scanpy as sc
import chronocci as ch

# 1. Load your preprocessed AnnData object
adata = sc.datasets.pbmc68k_reduced()

# 2. Compute the chronological multi-lineage signaling matrices
df_results = ch.run_joint_chronological_cci_pipeline(
    adata=adata,
    cluster_key="louvain",
    root_cell_type="CD14+ Monocytes",
    species="human"
)

# 3. Plot the continuous dot-relay cascade timeline map
ch.plot_cell_type_relay_timeline(
    adata=adata,
    df_final_sorted=df_results,
    cluster_key="louvain",
    top_n=12,
    output_name="my_blood_study"
)

# 4. Screen for asymmetric switches and plot the bifurcation map
df_asymmetric = ch.snoop_py (
    adata=adata,
    df_final_sorted=df_results,
    lineage_A="CD14+ Monocytes",
    lineage_B="CD34+ Dendritic"
)

ch.plot_signaling_bifurcation(
    adata=adata,
    df_final_sorted=df_asymmetric,
    lineage_A="CD14+ Monocytes",
    lineage_B="CD34+ Dendritic",
    top_n_pairs=4,
    output_name="my_blood_study"
)
```
