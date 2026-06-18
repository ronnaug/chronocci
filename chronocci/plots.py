import colorsys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import cellrank as cr

def plot_cell_type_relay_timeline(
    adata: sc.AnnData, 
    df_final_sorted: pd.DataFrame, 
    cluster_key: str, 
    top_n: int = 15, 
    figsize: tuple = (15, 9),
    output_name: str = "study"
):
    """Renders the chronological dot-relay timeline map."""
    if df_final_sorted.empty: return
    top_df = df_final_sorted.head(top_n).copy()
    
    donors, acceptors = [], []
    adata_genes_upper = {g.upper(): g for g in adata.var_names}

    for _, row in top_df.iterrows():
        lig_g = adata_genes_upper[row["ligand"].upper()]
        rec_g = adata_genes_upper[row["receptor"].upper()]
        peak_t = row["peak_pseudotime"]
        
        l_exp = adata[:, lig_g].X.toarray().flatten() if hasattr(adata[:, lig_g].X, "toarray") else adata[:, lig_g].X.flatten()
        r_exp = adata[:, rec_g].X.toarray().flatten() if hasattr(adata[:, rec_g].X, "toarray") else adata[:, rec_g].X.flatten()
        time_mask = (adata.obs["timeline_time"].values >= peak_t - 0.15) & (adata.obs["timeline_time"].values <= peak_t + 0.15)
        
        if np.sum(time_mask) > 0:
            df_peak = pd.DataFrame({"cluster": adata.obs[cluster_key].values[time_mask], "l": l_exp[time_mask], "r": r_exp[time_mask]})
            donors.append(df_peak.groupby("cluster", observed=True)["l"].mean().idxmax())
            acceptors.append(df_peak.groupby("cluster", observed=True)["r"].mean().idxmax())
        else:
            donors.append("Unknown"); acceptors.append("Unknown")
            
    top_df["donor_cell"], top_df["acceptor_cell"] = donors, acceptors
    top_df["communication_axis"] = top_df["donor_cell"].astype(str) + " ➔ " + top_df["acceptor_cell"].astype(str)
    top_df = top_df.sort_values(by="peak_pseudotime", ascending=True).reset_index(drop=True)

    unique_lineages = top_df["lineage"].unique()
    lineage_colors_dict = {lin: colorsys.hls_to_rgb(i/len(unique_lineages), 0.70, 0.55) for i, lin in enumerate(unique_lineages)}
    point_colors = [lineage_colors_dict[lin] for lin in top_df["lineage"]]

    plt.rcParams["font.family"] = "sans-serif"
    fig = plt.figure(figsize=figsize, dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[0.80, 0.20], wspace=0.05)
    ax, ax_leg = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    ax_leg.axis("off")
    
    #ax.hlines(y=top_df.index, xmin=-0.05, xmax=1.05, colors="#f5f5f5", zorder=1)
    size_factor = 4.0
    scatter = ax.scatter(top_df["peak_pseudotime"], top_df.index, s=top_df["max_signal_raw"]*size_factor, c=point_colors, edgecolors="#222222", linewidths=0.8, zorder=2, alpha=0.95)

    for i, row in top_df.iterrows():
        ax.text(row["peak_pseudotime"] + 0.015, i + 0.1, row["interaction_pair"], fontsize=9.5, weight="bold", color="#1a252f", va="bottom", ha="left")

    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.75, len(top_df) - 0.25)
    ax.set_yticks(top_df.index); ax.set_yticklabels(top_df["communication_axis"], fontsize=11, weight="bold", color="#2c3e50")
    ax.set_xlabel("Global Chronological Pseudotime", fontsize=11, labelpad=10)
    ax.set_title("Chronological Cascade of Inter-Lineage Cell Communications", fontsize=13, weight="bold", pad=20)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.spines["left"].set_color("#cccccc"); ax.spines["bottom"].set_color("#cccccc")

    color_lines = [plt.Line2D([], [], marker="o", color="w", markerfacecolor=lineage_colors_dict[lin], markeredgecolor="#222222", markersize=10, alpha=0.95) for lin in unique_lineages]
    legend_color = ax_leg.legend(color_lines, unique_lineages, title="Biological Lineage", title_fontsize=9.5, fontsize=9, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.75))
    ax_leg.add_artist(legend_color)

    legend_vals = np.round(np.linspace(top_df["max_signal_raw"].min(), top_df["max_signal_raw"].max(), 5), 1)
    size_lines = [plt.Line2D([], [], marker="o", color="w", markerfacecolor="#7f8c8d", markeredgecolor="#222222", markersize=np.sqrt(val * size_factor) * 1.5, alpha=0.7) for val in legend_vals]
    legend_size = ax_leg.legend(size_lines, [f"{val}" for val in legend_vals], title="Interaction Score", title_fontsize=9.5, fontsize=9, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.40))
    ax_leg.add_artist(legend_size)

    plt.subplots_adjust(left=0.20, right=0.95, top=0.90, bottom=0.12)
    plt.savefig(f"CCI_Relay_{output_name}.pdf", bbox_inches="tight", dpi=300)
    plt.show()


def plot_signaling_bifurcation(
    adata: sc.AnnData,
    df_final_sorted: pd.DataFrame,
    lineage_A: str,
    lineage_B: str,
    top_n_pairs: int = 6,
    figsize: tuple = (12, 8),
    output_name: str = "study"
):
    """Renders the mirrored Signaling Bifurcation Area Map with premium pastel palette."""
    if df_final_sorted.empty: return
    top_df = df_final_sorted.head(min(top_n_pairs, len(df_final_sorted))).copy()
    
    time_grid = np.linspace(0, 1, 30)
    model = cr.models.GAM(adata)
    adata_genes_upper = {g.upper(): g for g in adata.var_names}
    
    fig, axes = plt.subplots(len(top_df), 1, figsize=figsize, sharex=True, dpi=150)
    if len(top_df) == 1: axes = [axes]
        
    plt.rcParams["font.family"] = "sans-serif"
    
    # Эстетичная пастельная гамма
    color_A, color_B = "#A0C4DF", "#FBC4B6"
    edge_A, edge_B = "#78A2CC", "#EAA294"

    for idx, (df_idx, row_data) in enumerate(top_df.iterrows()):
        ax = axes[idx]
        pair_name = row_data["interaction_pair"]
        lig, rec = row_data["ligand"], row_data["receptor"]
        trends = {}
        
        for l in [lineage_A, lineage_B]:
            gene_trends_local = {}
            for gene_sym in [lig, rec]:
                real_g = adata_genes_upper[gene_sym.upper()]
                try:
                    model.fit(real_g, lineage=l, time_key="timeline_time")
                    _, y_pred, _ = model.predict(x_test=time_grid)
                    gene_trends_local[gene_sym.upper()] = np.clip(y_pred, 0, None)
                except:
                    w = np.array(adata.obsm["lineages_fwd"][l].X).flatten()
                    s_idx = np.argsort(adata.obs["timeline_time"].values)
                    t_s = adata.obs["timeline_time"].values[s_idx]
                    g_s = (adata[:, real_g].X.toarray().flatten() if hasattr(adata[:, real_g].X, "toarray") else adata[:, real_g].X.flatten())[s_idx]
                    
                    y_fallback = [np.mean(g_s[(t_s >= t-0.1) & (t_s <= t+0.1) & (w[s_idx] > 0)]) if np.sum((t_s >= t-0.1) & (t_s <= t+0.1) & (w[s_idx] > 0)) > 0 else 0.0 for t in time_grid]
                    gene_trends_local[gene_sym.upper()] = np.clip(np.array(y_fallback), 0, None)
                    
            trends[l] = gene_trends_local[lig.upper()] * gene_trends_local[rec.upper()]
            
        # Рассчитываем нормировочный максимум ПОСЛЕ того, как собрали оба тренда
        max_total = max(np.max(trends[lineage_A]), np.max(trends[lineage_B])) or 1.0
        
        # Отрисовка пути А с аккуратной полупрозрачной обводкой
        ax.fill_between(time_grid, trends[lineage_A]/max_total, 0, 
                        color=color_A, edgecolor=edge_A, linewidth=1.0, alpha=0.85, 
                        label=f"Path to {lineage_A}" if idx==0 else "")
                        
        # Отрисовка пути Б с аккуратной полупрозрачной обводкой (зеркально вниз)
        ax.fill_between(time_grid, -trends[lineage_B]/max_total, 0, 
                        color=color_B, edgecolor=edge_B, linewidth=1.0, alpha=0.85, 
                        label=f"Path to {lineage_B}" if idx==0 else "")
        
        # Аккуратная разделительная пунктирная ось
        ax.axhline(0, color="#b0b0b0", linestyle="--", linewidth=0.8, alpha=0.7)
        
        ax.set_ylabel("Intensity", fontsize=8)
        ax.set_title(f"Bifurcation of Dynamics: {pair_name}", fontsize=10, weight="bold", pad=5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#cccccc")
        ax.spines["bottom"].set_color("#cccccc")
        ax.get_yaxis().set_ticks([])

    axes[-1].set_xlim(0, 1)
    axes[-1].set_xlabel("Global Chronological Pseudotime", fontsize=11, labelpad=10)
    fig.suptitle(f"CCI Fate Decision Map: {lineage_A} vs {lineage_B}", fontsize=13, weight="bold", y=0.98)
    fig.legend(loc="upper right", bbox_to_anchor=(0.95, 0.98), frameon=False, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f"CCI_Bifurcation_{lineage_A}_vs_{lineage_B}_{output_name}.pdf", bbox_inches="tight", dpi=300)
    
    if matplotlib.get_backend().lower() != 'agg':
        plt.show()
    else:
        plt.close(fig)


def plot_signaling_streamgraph(
    adata: sc.AnnData,
    df_final_sorted: pd.DataFrame,
    cluster_key: str,
    time_grid: np.ndarray,
    figsize: tuple = (16, 9),
    output_name: str = "study"
):
    """
    Renders an integrative Streamgraph visualization showing the chronological cascade 
    of multi-lineage cell-cell interactions over pseudotime.
    """
    if df_final_sorted.empty:
        print("Warning: df_final_sorted is empty. Skipping streamgraph plot.")
        return df_final_sorted

    # Match genes regardless of lowercase/uppercase format in your dataframe
    adata_genes_upper = {g.upper(): g for g in adata.var_names}
    joint_labels = []
    
    print("Analyzing paracrine interactions across cell types at trajectory peaks...")
    
    for _, row in df_final_sorted.iterrows():
        lig_g = adata_genes_upper.get(row["ligand"].upper())
        rec_g = adata_genes_upper.get(row["receptor"].upper())
        peak_t = row["peak_pseudotime"]
        
        # Guard against missing genes in the object
        if not lig_g or not rec_g:
            joint_labels.append(f"{row['ligand']} → {row['receptor']} [Lineage: {row['lineage']}] (Gene Missing)")
            continue
            
        l_exp = adata[:, lig_g].X.toarray().flatten() if hasattr(adata[:, lig_g].X, "toarray") else adata[:, lig_g].X.flatten()
        r_exp = adata[:, rec_g].X.toarray().flatten() if hasattr(adata[:, rec_g].X, "toarray") else adata[:, rec_g].X.flatten()
        
        time_mask = (adata.obs["timeline_time"].values >= peak_t - 0.15) & (adata.obs["timeline_time"].values <= peak_t + 0.15)
        
        if np.sum(time_mask) > 0:
            df_peak_cells = pd.DataFrame({
                "cluster": adata.obs[cluster_key].values[time_mask],
                "l": l_exp[time_mask],
                "r": r_exp[time_mask]
            })
            donor = df_peak_cells.groupby("cluster", observed=True)["l"].mean().idxmax()
            acceptor = df_peak_cells.groupby("cluster", observed=True)["r"].mean().idxmax()
        else:
            donor, acceptor = "Unknown", "Unknown"
            
        joint_labels.append(f"{row['ligand']} ({donor}) → {row['receptor']} ({acceptor}) [Lineage: {row['lineage']}]")

    # Reconstruct the trajectories matrix into a 2D array
    trajectories_matrix = np.array(df_final_sorted["trajectory_norm"].tolist())

    # Generate aesthetic pastel colors dynamically
    n_colors = len(df_final_sorted)
    colors = [colorsys.hls_to_rgb(i / n_colors, 0.70, 0.55) for i in range(n_colors)]

    # Draw the streamgraph canvas
    plt.rcParams["font.family"] = "sans-serif"
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    
    ax.stackplot(
        time_grid, 
        trajectories_matrix, 
        labels=joint_labels, 
        baseline="wiggle", 
        colors=colors,            
        alpha=0.95,               
        edgecolor="white", 
        linewidth=0.4
    )
    
    # Styling and Labels
    ax.set_xlim(np.min(time_grid), np.max(time_grid))
    ax.set_xlabel("Global Chronological Pseudotime (Tissue Differentiation)", fontsize=11, labelpad=10)
    ax.set_ylabel("Relative Interaction Intensity (Stream Width)", fontsize=11, labelpad=10)
    ax.set_title("Unified Chronological Cascade of Multi-Lineage Interactions", fontsize=13, weight="bold", pad=15)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.get_yaxis().set_ticks([])
    
    # Legend Placement outside the chart
    ax.legend(
        loc="center left", 
        bbox_to_anchor=(1.02, 0.5), 
        title="Global Interactions Cascade\n[Donor → Acceptor (Trajectory)]", 
        frameon=False, 
        fontsize=9
    )
    
    plt.tight_layout()
    
    # Save Outputs
    stream_filename = f"CCI_Streamgraph_{output_name}.pdf"
    plt.savefig(stream_filename, bbox_inches="tight", dpi=300)
    print(f"Integrative Streamgraph successfully saved to '{stream_filename}'!")
    
    # Avoid execution halts on servers without graphical interfaces (like Github Actions)
    if mpl.get_backend().lower() != 'agg':
        plt.show()
    else:
        plt.close(fig)

    return df_final_sorted
