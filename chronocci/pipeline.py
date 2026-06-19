import cellrank as cr
import liana as li
import numpy as np
import pandas as pd
import scanpy as sc

def run_joint_chronological_cci_pipeline(
    adata: sc.AnnData,
    cluster_key: str,
    root_cell_type: str = None,
    species: str = "human",
    top_n_per_lineage: int = 4,
    n_lineages: int = None  
):
    """Executes the master chronological CellRank + LIANA processing pipeline with non-negative expression floors."""
    print(f"LAUNCHING JOINT CHRONOLOGICAL CCI PIPELINE...")
    
    if "iroot" in adata.uns: del adata.uns["iroot"]
    if "diffmap" not in adata.obsm: sc.tl.diffmap(adata)

    if root_cell_type is None:
        diffmap_coords = adata.obsm["X_diffmap"][:, 1:3]
        center_of_mass = np.mean(diffmap_coords, axis=0)
        distances_to_center = np.linalg.norm(diffmap_coords - center_of_mass, axis=1)
        root_cell_type = adata.obs[cluster_key].values[np.argmin(distances_to_center)]
        print(f"Automatically identified root: '{root_cell_type}'")

    root_cells_indices = np.where(adata.obs[cluster_key] == root_cell_type)
    adata.uns["iroot"] = int(root_cells_indices[0][0])

    sc.tl.dpt(adata)
    adata.obs["timeline_time"] = adata.obs["dpt_pseudotime"].astype(float)

    pk = cr.kernels.PseudotimeKernel(adata, time_key="timeline_time")
    pk.compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    combined_kernel = 0.8 * pk + 0.2 * ck

    estimator = cr.estimators.GPCCA(combined_kernel)
    estimator.compute_schur()

    states_to_find = n_lineages if n_lineages is not None else (2, 6)
    
    estimator.compute_macrostates(n_states=states_to_find,  n_cells=None)
    
    estimator.predict_terminal_states()
    estimator.compute_fate_probabilities(solver="direct", check_sum_tol=1e-2)

    adata.obsm["lineages_fwd"] = estimator.fate_probabilities
    available_lineages = list(estimator.fate_probabilities.names)

    resource_name = "consensus" if species.lower() == "human" else "mouseconsensus"
    lr_db = li.resource.select_resource(resource_name)
    lr_db["ligand_upper"] = lr_db["ligand"].str.upper()
    lr_db["receptor_upper"] = lr_db["receptor"].str.upper()
    adata_genes_upper = {g.upper(): g for g in adata.var_names}
    
    valid_pairs = lr_db[
        lr_db["ligand_upper"].isin(adata_genes_upper.keys()) & 
        lr_db["receptor_upper"].isin(adata_genes_upper.keys())
    ].copy()

    time_grid = np.linspace(0, 1, 30)
    model = cr.models.GAM(adata)
    lineage_dfs = []

    for lineage in available_lineages:
        print(f"Modeling temporal profiles for lineage: {lineage}...")
        gene_trends = {}
        genes_to_fit = set(valid_pairs["ligand_upper"]).union(set(valid_pairs["receptor_upper"]))
        
        for gene in genes_to_fit:
            real_gene_name = adata_genes_upper[gene]
            try:
                model.fit(real_gene_name, lineage=lineage, time_key="timeline_time")
                _, y_pred, _ = model.predict(x_test=time_grid)
                # FIX 1: Clip negative extrapolated values from GAM model to absolute zero
                gene_trends[gene] = np.clip(y_pred, 0, None)
            except:
                fate_weights = np.array(adata.obsm["lineages_fwd"][lineage].X).flatten()
                sorted_idx = np.argsort(adata.obs["timeline_time"].values)
                t_sort = adata.obs["timeline_time"].values[sorted_idx]
                g_sort = (adata[:, real_gene_name].X.toarray().flatten() if hasattr(adata[:, real_gene_name].X, "toarray") else adata[:, real_gene_name].X.flatten())[sorted_idx]
                w_sort = fate_weights[sorted_idx]
                
                y_pred = []
                for t in time_grid:
                    mask = (t_sort >= t - 0.1) & (t_sort <= t + 0.1)
                    if np.sum(mask) > 0 and np.sum(w_sort[mask]) > 0:
                        y_pred.append(np.average(g_sort[mask], weights=w_sort[mask]))
                    else:
                        y_pred.append(0.0)
                # FIX 2: Safeguard fallback window tracking values from negative raw values
                gene_trends[gene] = np.clip(np.array(y_pred), 0, None)

        df_trends = pd.DataFrame(gene_trends, index=time_grid)

        cci_results = []
        for _, row in valid_pairs.iterrows():
            lig, rec = row["ligand_upper"], row["receptor_upper"]
            if lig in df_trends.columns and rec in df_trends.columns:
                score_vector = df_trends[lig].values * df_trends[rec].values
                max_val = np.max(score_vector)
                if max_val <= 1e-9: continue # Discard mathematically dead combinations
                
                cci_results.append({
                    "interaction_pair": f"{row['ligand']}_{row['receptor']}",
                    "ligand": row["ligand"], "receptor": row["receptor"],
                    "lineage": lineage, "max_signal_raw": max_val,
                    "trajectory": score_vector,
                    "peak_pseudotime": time_grid[np.argmax(score_vector)]
                })

        df_lineage_cci = pd.DataFrame(cci_results)
        if not df_lineage_cci.empty:
            max_val_raw = df_lineage_cci["max_signal_raw"]
            df_lineage_cci["z_score_lineage"] = (max_val_raw - max_val_raw.mean()) / (max_val_raw.std() + 1e-6)
            df_lineage_cci = df_lineage_cci.sort_values(by="z_score_lineage", ascending=False).head(top_n_per_lineage)
            lineage_dfs.append(df_lineage_cci)

    df_global_pool = pd.concat(lineage_dfs, ignore_index=True)
    df_global_pool = df_global_pool.sort_values(by="z_score_lineage", ascending=False)
    df_unique_pairs = df_global_pool.drop_duplicates(subset=["interaction_pair"]).copy()

    norm_trajectories = []
    for _, row in df_unique_pairs.iterrows():
        vec = row["trajectory"]
        norm_trajectories.append(vec / np.mean(vec) if np.mean(vec) > 0 else vec)
    df_unique_pairs["trajectory_norm"] = norm_trajectories

    df_final_sorted = df_unique_pairs.sort_values(by="peak_pseudotime").reset_index(drop=True)
    return df_final_sorted



def snoop_py (
    adata: sc.AnnData,
    df_final_sorted: pd.DataFrame,
    lineage_A: str,
    lineage_B: str
) -> pd.DataFrame:
    """Ranks and filters ligand-receptor pairs based on chronological MAE asymmetry with non-negative floors."""
    print(f"Screening for highly asymmetric CCI markers:\n{lineage_A} vs {lineage_B}...")

    time_grid = np.linspace(0, 1, 30)
    model = cr.models.GAM(adata)
    adata_genes_upper = {g.upper(): g for g in adata.var_names}
    bifurcation_scores = []
    
    for idx, row in df_final_sorted.iterrows():
        lig, rec = row["ligand"], row["receptor"]
        trends = {}
        for l in [lineage_A, lineage_B]:
            gene_trends_local = {}
            for gene_sym in [lig, rec]:
                real_g = adata_genes_upper[gene_sym.upper()]
                try:
                    model.fit(real_g, lineage=l, time_key="timeline_time")
                    _, y_pred, _ = model.predict(x_test=time_grid)
                    # FIX 1: Enforce zero floor on snoop's internal GAM predictions
                    gene_trends_local[gene_sym.upper()] = np.clip(y_pred, 0, None)
                except:
                    w = np.array(adata.obsm["lineages_fwd"][l].X).flatten()
                    s_idx = np.argsort(adata.obs["timeline_time"].values)
                    t_s = adata.obs["timeline_time"].values[s_idx]
                    g_s = (adata[:, real_g].X.toarray().flatten() if hasattr(adata[:, real_g].X, "toarray") else adata[:, real_g].X.flatten())[s_idx]
                    
                    y_fallback = [np.mean(g_s[(t_s >= t-0.1) & (t_s <= t+0.1) & (w[s_idx] > 0)]) if np.sum((t_s >= t-0.1) & (t_s <= t+0.1) & (w[s_idx] > 0)) > 0 else 0.0 for t in time_grid]
                    # FIX 2: Enforce zero floor on snoop's fallback mechanism
                    gene_trends_local[gene_sym.upper()] = np.clip(np.array(y_fallback), 0, None)
                    
            trends[l] = gene_trends_local[lig.upper()] * gene_trends_local[rec.upper()]
            
        vec_A, vec_B = trends[lineage_A], trends[lineage_B]
        max_total = max(np.max(vec_A), np.max(vec_B))
        if max_total > 0:
            bif_score = np.mean(np.abs((vec_A / max_total) - (vec_B / max_total)))
        else:
            bif_score = 0.0
        bifurcation_scores.append(bif_score)
        
    df_bif = df_final_sorted.copy()
    df_bif["bifurcation_score"] = bifurcation_scores
    df_bif = df_bif.sort_values(by="bifurcation_score", ascending=False).reset_index(drop=True)
    return df_bif

