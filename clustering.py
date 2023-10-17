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
