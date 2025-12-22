# ============================================================
# Clustered Jaccard Heatmaps (Normal + Abnormal)
# - Keeps original layout/area
# - Adds spacing between labels and heatmap (pad)
# - Shows figures in Jupyter output
# - Also saves PNGs
# ============================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage


# -------------------------
# Font sizes (adjust if needed)

#BASE_FONT = 14          # General text size
#TICK_FONT = 13          # X/Y tick labels
#ANNOT_FONT = 13         # Numbers inside heatmap cells
#TITLE_FONT = 18         # Title size
#CBAR_TICK_FONT = 12     # Colorbar tick labels
# -------------------------


BASE_FONT = 15
TICK_FONT = 16
ANNOT_FONT = 16
TITLE_FONT = 16
CBAR_TICK_FONT = 12

plt.rcParams.update({
    "font.size": BASE_FONT,
    "xtick.labelsize": TICK_FONT,
    "ytick.labelsize": TICK_FONT,
})
sns.set_theme(style="white")


# -------------------------
# Normal dataset
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
# Abnormal dataset
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
# Convert text -> matrix
# -------------------------
def prepare_matrix(data_str: str) -> pd.DataFrame:
    df = pd.read_csv(StringIO(data_str.strip()))
    mat = df.pivot(index="ToolA", columns="ToolB", values="Jaccard")

    tools = sorted(set(df["ToolA"]).union(set(df["ToolB"])))
    mat = mat.reindex(index=tools, columns=tools).fillna(0.0)

    # Ensure diagonal is 1
    for t in tools:
        mat.loc[t, t] = 1.0
    return mat


# -------------------------
# Short labels mapping
# -------------------------
def shorten_labels(matrix: pd.DataFrame) -> pd.DataFrame:
    label_map = {
        "BCF_normal_SNPs": "BCF",
        "Deepvariant_normal_SNPs": "DeepVariant",
        "Freebayes_normal_SNPs": "FreeBayes",
        "Gatk_normal_SNPs": "GATK",
        "Varscan_normal_SNPs": "VarScan",
        "BCF_abnormal_SNPs": "BCF",
        "Deepvariant_abnormal_SNPs": "DeepVariant",
        "Freebayes_abnormal_SNPs": "FreeBayes",
        "Gatk_abnormal_SNPs": "GATK",
        "Varscan_abnormal_SNPs": "VarScan",
    }
    return matrix.rename(index=label_map, columns=label_map)


# -------------------------
# Plot function (keeps old area; adds label spacing; shows output)
# -------------------------
def plot_clustered_heatmap(matrix: pd.DataFrame, title: str, output_filename: str):
    # distance = 1 - similarity
    dist = 1.0 - matrix.values
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)

    Z = linkage(squareform(dist, checks=False), method="average")

    # NOTE: Keeping your original figure size and cbar position
    g = sns.clustermap(
        matrix,
        row_linkage=Z,
        col_linkage=Z,
        cmap="viridis",
        linewidths=0.5,
        annot=True,
        fmt=".3f",
        annot_kws={"size": ANNOT_FONT},
        figsize=(9, 7),                      # ✅ old area
        cbar_pos=(0.02, 0.82, 0.03, 0.16)     # ✅ old area
    )

    # Title
    g.fig.suptitle(title, y=1.02, fontsize=TITLE_FONT, fontweight="bold")

    # ✅ Fix X labels (rotate + anchor + extra distance)
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(),
        rotation=45,
        ha="right",
        rotation_mode="anchor"
    )
    g.ax_heatmap.tick_params(axis="x", labelsize=TICK_FONT, pad=10)  # <-- distance from heatmap

    # ✅ Fix Y labels (keep horizontal + extra distance)
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(),
        rotation=0,
        va="center"
    )
    g.ax_heatmap.tick_params(axis="y", labelsize=TICK_FONT, pad=12)  # <-- distance from heatmap

    # Colorbar tick font
    if g.cax is not None:
        g.cax.tick_params(labelsize=CBAR_TICK_FONT)

    # Save (with a little padding so labels never cut)
    g.savefig(output_filename, dpi=300, bbox_inches="tight", pad_inches=0.25)

    # ✅ Show in Jupyter output
    plt.show()

    print(f"✅ Saved: {output_filename}")


# -------------------------
# Run
# -------------------------
normal_matrix = shorten_labels(prepare_matrix(normal_data))
abnormal_matrix = shorten_labels(prepare_matrix(abnormal_data))

plot_clustered_heatmap(
    normal_matrix,
    "(A)",
    "normal_heatmap_clustered_shortlabels.png"
)

plot_clustered_heatmap(
    abnormal_matrix,
    "(B)",
    "abnormal_heatmap_clustered_shortlabels.png"
)
