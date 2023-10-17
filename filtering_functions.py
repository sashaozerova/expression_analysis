import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


filename = 'normalized_counts.csv'  # correct filename here
df = pd.read_csv(filename, index_col=0)

# Plot expression values disctribution for each gene
expression_threshold = 1  # chosen threshold for log-transformed expression values here
sns.displot(np.log2(df).mean(axis=1).tolist(), kde=True)
plt.axvline(x=expression_threshold, color='red')
plt.plot()

# plot cumulative distribution
sns.displot(np.log2(df).mean(axis=1).tolist(), kind="ecdf")
plt.axvline(x=expression_threshold, color='red')
plt.plot()

expr_filtered_df = df[df.mean(axis=1) >= 2**expression_threshold]

# Plot STD values distribution for each gene. Chosen by peak density value
kde = sns.kdeplot(np.log2(expr_filtered_df).std(axis=1).tolist())
x_kde_values = kde.lines[0].get_xdata() # Get the x data of the distribution
y_kde_values = kde.lines[0].get_ydata() # Get the y data of the distribution
maxid = np.argmax(y_kde_values) # The id of the peak (maximum of y data)

std_threshold = x_kde_values[maxid]
sns.displot(np.log2(expr_filtered_df).std(axis=1).tolist(), kde=True)
plt.axvline(x=std_threshold, color='red')
plt.plot()

sns.displot(np.log2(expr_filtered_df).std(axis=1).tolist(), kind="ecdf")
plt.axvline(x=std_threshold, color='red')
plt.plot()

std_expr_filtered_df = expr_filtered_df[np.log2(expr_filtered_df).std(axis=1) >= std_threshold]
std_expr_filtered_df.to_csv('filtered_normalized_counts.csv')