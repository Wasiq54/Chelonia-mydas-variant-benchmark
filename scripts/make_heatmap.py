# =========================
# Clustered Jaccard Heatmaps (Normal + Abnormal)
# Based on YOUR results
# =========================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage


# -------------------------
# Step 1: Your Normal dataset
# -------------------------
normal_data = """
ToolA,ToolB,Jaccard
BCF_normal_SNPs,BCF_normal_SNPs,1.000000
BCF_normal_SNPs,Deepvariant_normal_SNPs,0.930868
BCF_normal_SNPs,Freebayes_normal_SNPs,0.803898
BCF_normal_SNPs,Gatk_normal_SNPs,0.937860
BCF_normal_SNPs,Varscan_normal_SNPs,0.943323

Deepvariant_normal_SNPs,BCF_normal_SNPs,0.930868
Deepvariant_normal_SNPs,Deepvariant_normal_SNPs,1.000000
Deepvariant_normal_SNPs,Freebayes_normal_SNPs,0.786252
Deepvariant_normal_SNPs,Gatk_normal_SNPs,0.918006
Deepvariant_normal_SNPs,Varscan_normal_SNPs,0.928128

Freebayes_normal_SNPs,BCF_normal_SNPs,0.803898
Freebayes_normal_SNPs,Deepvariant_normal_SNPs,0.786252
Freebayes_normal_SNPs,Freebayes_normal_SNPs,1.000000
Freebayes_normal_SNPs,Gatk_normal_SNPs,0.809775
Freebayes_normal_SNPs,Varscan_normal_SNPs,0.804201

Gatk_normal_SNPs,BCF_normal_SNPs,0.937860
Gatk_normal_SNPs,Deepvariant_normal_SNPs,0.918006
Gatk_normal_SNPs,Freebayes_normal_SNPs,0.809775
Gatk_normal_SNPs,Gatk_normal_SNPs,1.000000
Gatk_normal_SNPs,Varscan_normal_SNPs,0.920129

Varscan_normal_SNPs,BCF_normal_SNPs,0.943323
Varscan_normal_SNPs,Deepvariant_normal_SNPs,0.928128
Varscan_normal_SNPs,Freebayes_normal_SNPs,0.804201
Varscan_normal_SNPs,Gatk_normal_SNPs,0.920129
Varscan_normal_SNPs,Varscan_normal_SNPs,1.000000
"""


# -------------------------
# Step 2: Your Abnormal dataset
# -------------------------
abnormal_data = """
ToolA,ToolB,Jaccard
BCF_abnormal_SNPs,BCF_abnormal_SNPs,1.000000
BCF_abnormal_SNPs,Deepvariant_abnormal_SNPs,0.931014
BCF_abnormal_SNPs,Freebayes_abnormal_SNPs,0.804477
BCF_abnormal_SNPs,Gatk_abnormal_SNPs,0.937287
BCF_abnormal_SNPs,Varscan_abnormal_SNPs,0.944160

Deepvariant_abnormal_SNPs,BCF_abnormal_SNPs,0.931014
Deepvariant_abnormal_SNPs,Deepvariant_abnormal_SNPs,1.000000
Deepvariant_abnormal_SNPs,Freebayes_abnormal_SNPs,0.787334
Deepvariant_abnormal_SNPs,Gatk_abnormal_SNPs,0.918906
Deepvariant_abnormal_SNPs,Varscan_abnormal_SNPs,0.929124

Freebayes_abnormal_SNPs,BCF_abnormal_SNPs,0.804477
Freebayes_abnormal_SNPs,Deepvariant_abnormal_SNPs,0.787334
Freebayes_abnormal_SNPs,Freebayes_abnormal_SNPs,1.000000
Freebayes_abnormal_SNPs,Gatk_abnormal_SNPs,0.809711
Freebayes_abnormal_SNPs,Varscan_abnormal_SNPs,0.804360

Gatk_abnormal_SNPs,BCF_abnormal_SNPs,0.937287
Gatk_abnormal_SNPs,Deepvariant_abnormal_SNPs,0.918906
Gatk_abnormal_SNPs,Freebayes_abnormal_SNPs,0.809711
Gatk_abnormal_SNPs,Gatk_abnormal_SNPs,1.000000
Gatk_abnormal_SNPs,Varscan_abnormal_SNPs,0.919892

Varscan_abnormal_SNPs,BCF_abnormal_SNPs,0.944160
Varscan_abnormal_SNPs,Deepvariant_abnormal_SNPs,0.929124
Varscan_abnormal_SNPs,Freebayes_abnormal_SNPs,0.804360
Varscan_abnormal_SNPs,Gatk_abnormal_SNPs,0.919892
Varscan_abnormal_SNPs,Varscan_abnormal_SNPs,1.000000
"""


# -------------------------
# Step 3: Prepare matrix
# -------------------------
def prepare_matrix(data_str: str) -> pd.DataFrame:
    df = pd.read_csv(StringIO(data_str.strip()))
    mat = df.pivot(index="ToolA", columns="ToolB", values="Jaccard")

    # Ensure full square matrix (same tools in rows/cols)
    tools = sorted(set(df["ToolA"]).union(set(df["ToolB"])))
    mat = mat.reindex(index=tools, columns=tools)

    # Fill missing values safely (should not happen if complete)
    mat = mat.fillna(0.0)

    # Force diagonal = 1 (safety)
    for t in tools:
        mat.loc[t, t] = 1.0

    return mat


# -------------------------
# Step 4: Clustered heatmap function
# -------------------------
def plot_clustered_heatmap(matrix: pd.DataFrame, title: str, output_filename: str):
    # Distance matrix = 1 - similarity
    dist = 1.0 - matrix.values

    # Make it perfectly symmetric (avoid tiny float issues)
    dist = (dist + dist.T) / 2.0

    # Diagonal must be 0 for distance
    for i in range(dist.shape[0]):
        dist[i, i] = 0.0

    # Convert to condensed format for linkage
    condensed = squareform(dist, checks=False)

    # Hierarchical clustering
    Z = linkage(condensed, method="average")

    # Plot
    g = sns.clustermap(
        matrix,
        row_linkage=Z,
        col_linkage=Z,
        cmap="viridis",
        linewidths=0.5,
        annot=True,
        fmt=".3f",
        figsize=(11, 9),
        cbar_pos=(0.02, 0.82, 0.03, 0.16)
    )

    g.fig.suptitle(title, y=1.02, fontsize=16, fontweight="bold")
    g.savefig(output_filename, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"✅ Saved: {output_filename}")


# -------------------------
# Step 5: Run for Normal and Abnormal
# -------------------------
normal_matrix = prepare_matrix(normal_data)
abnormal_matrix = prepare_matrix(abnormal_data)

plot_clustered_heatmap(normal_matrix, "Normal Dataset (Jaccard Similarity)", "normal_heatmap_clustered.png")
plot_clustered_heatmap(abnormal_matrix, "Abnormal Dataset (Jaccard Similarity)", "abnormal_heatmap_clustered.png")
