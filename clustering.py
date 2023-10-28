import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as spc


filename = 'filtered_normalized_counts.csv'  # correct filename here
df = pd.read_csv(filename, index_col=0)

transposed = np.log2(df.transpose())
corr = transposed.corr(method="spearman")  # be careful, could take much time
corr.to_csv('gene_pairwise_correlation.csv')

# create linkage matrix with cluster distances
pdist = spc.distance.pdist(corr)
linkage = spc.linkage(pdist, method="complete")
linkage = pd.DataFrame(linkage)
linkage.columns = ['cluster1', 'cluster2', 'size', 'elements']
linkage.to_csv('correlation_based_linkage.csv.gz', compression="gzip", index=False)

# can further use cut_tree function to cut linkage into specified number of clusters

# plot expression dynamics in each cluster
cluster_mapping = {
    1: ["gene1", "gene2", "gene3"],
    2: ["gene4", "gene5"],
}  # insert proper mapping here (keys are cluster names, values - gene arrays)
num_clusters = 10
v_size, h_size = 2, 5  # dimentions of the plot

log2_df = np.log2(df)
f, axes = plt.subplots(v_size, h_size, figsize=(30, 10))
clusters = list(cluster_mapping.keys())

for i, cluster in enumerate(sorted(clusters)):
    # each cluster as a separate subplot
    current_df = log2_df.loc[set(cluster_mapping[cluster]).intersection(log2_df.index)]
    ax = axes[i // h_size][i - h_size * (i // h_size)]
    for name in current_df.index:
        ax.plot(current_df.columns.tolist(), current_df.loc[name].tolist(), c="grey", alpha=0.1)

    # plot median trajectory
    ax.plot(current_df.columns.tolist(), current_df.median().tolist(), c="red")

    # add stage for the lower subplots only
    if i // h_size == v_size - 1:
        stages = current_df.columns.tolist()
        max_lables = 18
        ratio = current_df.shape[1] // max_lables
        if ratio > 1:
            stages = [name if ind % ratio == 0 else "" for ind, name in enumerate(stages)]
        ax.set_xticklabels(stages, rotation=90, fontsize=18)
    else:
        ax.set_xticklabels([])
    # add names to the subplots
    ax.title.set_text(f"cluster {cluster}")
    ax.title.set_size(18)

    # add number of genes to each subplot
    ymin, ymax = ax.get_ylim()
    y_text = ymax - (ymax - ymin) / 7
    ax.text(0, y_text, f"{current_df.shape[0]} genes", fontsize=18)

plt.show()
