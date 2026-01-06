# to run on the dataset of scores on our designed binders (i.e. no labelling, some sort of unsupervised approach). This would allow us to understand better our dataset's structure, and see if we have some clusters.
# outline: clean dataset, standardize, run PCA, plot in 3D space and see scores contribution, try umap and clustering? 


# 1--------------df loading and cleaning: ---------------------------
import pandas as pd
import numpy as np

# load your merged metrics dataframe
df = pd.read_csv("metrics.csv")

# keep a copy of identifiers / metadata
non_numeric_cols = df.select_dtypes(exclude=[np.number]).columns
df_meta = df[non_numeric_cols] #should include: ID, 

# numeric data only
X = df.select_dtypes(include=[np.number])
# look up the "weird" cols: ...pDockQ_min, _pDockQ_max, pDockQ2_min, pDockQ2_max, everything with ipSAE
import seaborn as sns
#_ = sns.pairplot(X.loc[:, ~X.columns.str.contains("pDockQ_|pDockQ2|ipSAE|LIS",  case=False)])
#_ = sns.pairplot(X.loc[:, ~X.columns.str.contains("ipSAE",  case=False)])
#plt.savefig("./ipSAE_cols.png")
#plt.close()

# remove those cols for now:
# af3binary_bin_has_clash, for now : ...pDockQ_min, _pDockQ_max, pDockQ2_min, pDockQ2_max, everything with ipSAE
X = X.loc[:, ~df.columns.str.contains("pDockQ_|pDockQ2|ipSAE|LIS",  case=False)]


# inspect:
X.head

#Sanity check 
#assert X.shape[1] > 1, "PCA requires at least 2 numeric variables"

#  2-------------Standardize the data----------------------------------
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 3----------run PCA-------------
from sklearn.decomposition import PCA

pca = PCA(n_components=3)
X_pca = pca.fit_transform(X_scaled)

explained_var = pca.explained_variance_ratio_

for i, v in enumerate(explained_var, start=1):
    print(f"PC{i}: {v:.3f}")
  
# 4 -------------3D PCA scatter plot--------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")

ax.scatter(
    X_pca[:, 0],
    X_pca[:, 1],
    X_pca[:, 2],
    alpha=0.7
)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")

plt.tight_layout()
plt.savefig("./output/analyses/3D_PCA.png")
plt.close()


# 5 ----------metrics biplot----------------:
loadings = pca.components_.T
feature_names = X.columns
# maybe filter the loadings to keep only the most important, ie, longest arrows (not the most correlated)
"""pc_x = 0  # PC1
pc_y = 1  # PC2

arrow_lengths = np.sqrt(
    loadings[:, pc_x]**2 + loadings[:, pc_y]**2
)
top_idx = np.argsort(arrow_lengths)[-2:]

plt.figure(figsize=(8, 8))
plt.scatter(X_pca[:, pc_x], X_pca[:, pc_y], alpha=0.4)

scale = 3  # facteur visuel, inchangé mathématiquement

for i in top_idx:
    plt.arrow(
        0, 0,
        loadings[i, pc_x] * scale,
        loadings[i, pc_y] * scale,
        head_width=0.05,
        alpha=0.9
    )
    plt.text(
        loadings[i, pc_x] * scale * 1.1,
        loadings[i, pc_y] * scale * 1.1,
        feature_names[i],
        fontsize=10,
        weight="bold"
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.axhline(0, color="grey", linewidth=0.8)
plt.axvline(0, color="grey", linewidth=0.8)
plt.grid(alpha=0.3)

plt.savefig("./biplot_PC1_PC2_top2.png")
plt.close()

"""

"""
plt.figure(figsize=(8, 8))
plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.4)

for i, feature in enumerate(feature_names):
    plt.arrow(
        0, 0,
        loadings[i, 0] * 3,
        loadings[i, 1] * 3,
        head_width=0.05,
        alpha=0.8
    )
    plt.text(
        loadings[i, 0] * 3.2,
        loadings[i, 1] * 3.2,
        feature,
        fontsize=9
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.axhline(0)
plt.axvline(0)
plt.grid()
plt.savefig("./output/analyses/biplot.png")
plt.close()
"""


def plot_biplot(
    X_pca,
    loadings,
    feature_names,
    pc_x,
    pc_y,
    output_path,
    scale_arrows=3,
    n_top_features=2
):
    """
    pc_x, pc_y: indices PCA (0-based)
    """

    # contribution des variables au plan (PCx, PCy)
    contrib = loadings[:, pc_x]**2 + loadings[:, pc_y]**2

    # indices des variables les plus contributives
    top_idx = np.argsort(contrib)[-n_top_features:]

    plt.figure(figsize=(8, 8))
    plt.scatter(
        X_pca[:, pc_x],
        X_pca[:, pc_y],
        alpha=0.4
    )

    for i in top_idx:
        plt.arrow(
            0, 0,
            loadings[i, pc_x] * scale_arrows,
            loadings[i, pc_y] * scale_arrows,
            head_width=0.05,
            alpha=0.9
        )
        plt.text(
            loadings[i, pc_x] * scale_arrows * 1.1,
            loadings[i, pc_y] * scale_arrows * 1.1,
            feature_names[i],
            fontsize=10,
            weight="bold"
        )

    plt.xlabel(f"PC{pc_x+1}")
    plt.ylabel(f"PC{pc_y+1}")
    plt.axhline(0, color="grey", linewidth=0.8)
    plt.axvline(0, color="grey", linewidth=0.8)
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()



# 5 ---------- biplots PC1/PC2 et PC2/PC3 ----------------

plot_biplot(
    X_pca=X_pca,
    loadings=loadings,
    feature_names=feature_names,
    pc_x=0,
    pc_y=1,
    output_path="./biplot_PC1_PC2.png"
)

plot_biplot(
    X_pca=X_pca,
    loadings=loadings,
    feature_names=feature_names,
    pc_x=1,
    pc_y=2,
    output_path="./biplot_PC2_PC3.png"
)


assert "af3ternary_ter_iptm" in df.columns
assert "af3binary_bin_iptm" in df.columns
highlight_mask = (
    (df["af3ternary_ter_iptm"] >= 0.7) &
    (df["af3binary_bin_iptm"] <= 0.5)
)

### ------------------------ threshold biplot------------
def plot_biplot(
    X_pca,
    loadings,
    feature_names,
    pc_x,
    pc_y,
    output_path,
    scale_arrows=3,
    n_top_features=2,
    highlight_mask=None,
    highlight_label="Highlighted",
    highlight_color="red"
):
    """
    pc_x, pc_y: indices PCA (0-based)
    highlight_mask: boolean array of shape (n_samples,)
    """

    # contribution des variables au plan (PCx, PCy)
    contrib = loadings[:, pc_x]**2 + loadings[:, pc_y]**2

    # indices des variables les plus contributives
    top_idx = np.argsort(contrib)[-n_top_features:]

    plt.figure(figsize=(8, 8))

    # points non sélectionnés
    if highlight_mask is None:
        plt.scatter(
            X_pca[:, pc_x],
            X_pca[:, pc_y],
            alpha=0.4,
            s=25
        )
    else:
        highlight_mask = np.asarray(highlight_mask, dtype=bool)

        plt.scatter(
            X_pca[~highlight_mask, pc_x],
            X_pca[~highlight_mask, pc_y],
            alpha=0.4,
            s=25,
            label="Other samples"
        )

        plt.scatter(
            X_pca[highlight_mask, pc_x],
            X_pca[highlight_mask, pc_y],
            alpha=0.9,
            s=40,
            color=highlight_color,
            label=highlight_label
        )

    # flèches (top contributive variables)
    for i in top_idx:
        plt.arrow(
            0, 0,
            loadings[i, pc_x] * scale_arrows,
            loadings[i, pc_y] * scale_arrows,
            head_width=0.05,
            alpha=0.9
        )
        plt.text(
            loadings[i, pc_x] * scale_arrows * 1.1,
            loadings[i, pc_y] * scale_arrows * 1.1,
            feature_names[i],
            fontsize=10,
            weight="bold"
        )

    plt.xlabel(f"PC{pc_x+1}")
    plt.ylabel(f"PC{pc_y+1}")
    plt.axhline(0, color="grey", linewidth=0.8)
    plt.axvline(0, color="grey", linewidth=0.8)
    plt.grid(alpha=0.3)

    if highlight_mask is not None:
        plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
highlight_mask = (
    (df["af3ternary_ter_iptm"] >= 0.7) &
    (df["af3binary_bin_iptm"] <= 0.5)
)
plot_biplot(
    X_pca=X_pca,
    loadings=loadings,
    feature_names=feature_names,
    pc_x=0,
    pc_y=1,
    output_path="./biplot_PC1_PC2_thresholded.png",
    highlight_mask=highlight_mask,
    highlight_label="ter_iptm ≥ 0.7 & bin_iptm ≤ 0.5"
)


### ----------------------------------------------------------heatmap des corrélations entre variables----------------------------------------------------------
corr = X.corr()
plt.figure(figsize=(12, 10))

sns.heatmap(
    corr,
    cmap="vlag",
    center=0,
    square=True,
    linewidths=0.5,
    cbar_kws={"label": "Pearson correlation"}
)

plt.xticks(rotation=90)
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig("./correlation_heatmap.png", dpi=300)
plt.close()

# heatmap avec clusters:
sns.clustermap(
    corr,
    cmap="vlag",
    center=0,
    linewidths=0.5,
    figsize=(14, 14),
    cbar_kws={"label": "Pearson correlation"}
)

plt.savefig("./correlation_clustermap.png", dpi=300)
plt.close()



### -----------------------to plot a scatter plot of pc1pc2 with a variable colorgradient

pc1 = X_pca[:, 0]
pc2 = X_pca[:, 1]

color_var = df["af3binary_bin_iptm"].values
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))

sc = plt.scatter(
    pc1,
    pc2,
    c=color_var,
    cmap="viridis",
    s=30,
    alpha=0.8
)

plt.xlabel("PC1")
plt.ylabel("PC2")

cbar = plt.colorbar(sc)
cbar.set_label("af3binary_bin_iptm")

plt.axhline(0, color="grey", linewidth=0.5)
plt.axvline(0, color="grey", linewidth=0.5)

plt.tight_layout()
plt.savefig("./PCA_colored_af3_binary_iptm.png", dpi=300)
plt.close()

