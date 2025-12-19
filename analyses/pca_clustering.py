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
_ = sns.pairplot(X.loc[:, ~df.columns.str.contains("pDockQ_|pDockQ2|ipSAE|LIS",  case=False)])
plt.savefig("./output/analyses/weird_cols.png")
plt.close()

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


