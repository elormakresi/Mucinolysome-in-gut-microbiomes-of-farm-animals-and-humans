import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm

# Read input data and labels
data = open("final_metrics.tsv").readlines()[1:]
f_list = [line.strip() for line in open("doc_5_order_cper.txt")]
g_list = [line.strip() for line in open("coh_5_order_cper.txt")]


print (data, g_list, f_list)

# Metrics info for plotting
metrics = [
    {"index": 3, "bounds": [0.0, 0.7, 0.8, 0.9, 1.0], "annot_thresh": 0.7, "fmt": ".2f", "title": "ipTM Score"},
    {"index": 4, "bounds": [0.0, 0.8, 0.9, 1.0], "annot_thresh": 0.8, "fmt": ".2f", "title": "pTM Score"},
    {"index": 5, "bounds": [0.0, 0.7, 0.8, 0.9, 1.0], "annot_thresh": 0.7, "fmt": ".2f", "title": "Ranking Score"},
    {"index": 6, "bounds": [0, 50, 60, 70, 80, 90, 100], "annot_thresh": 50, "fmt": ".0f", "title": "# Contacts"},
    {"index": 7, "bounds": [0, 80, 90, 100], "annot_thresh": 80, "fmt": ".2f", "title": "Avg IF pLDDT"},
    {"index": 8, "bounds": [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], "annot_thresh": 0.3, "fmt": ".2f", "title": "pDockQ Score"},
]

# Initialize plot
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 9))
axes = axes.flatten()

# Generate each subplot
for i, metric in enumerate(metrics):
    matrix = pd.DataFrame(index=g_list, columns=f_list)

    for d in data:
        fields = d.strip().split("\t")
        combo = fields[0]
        value = float(fields[metric["index"]])
        if "-" in combo:
            g, f = combo.split("-", 1)
            print (value, combo, )
            if f in f_list and g in g_list:
                matrix.at[g, f] = value

    matrix = matrix.astype(float)
    print (matrix)
    #statistics of the matrix
    print (matrix)
    # Flatten and drop NaNs
    values = matrix.values.flatten()
    values = values[~np.isnan(values)]

    # Compute statistics
    mean_val = np.mean(values)
    median_val = np.median(values)
    min_val = np.min(values)
    max_val = np.max(values)

    # Apply threshold
    thresh = metric["annot_thresh"]
    above_thresh = values[values >= thresh]
    count_above = len(above_thresh)
    total_count = len(values)
    percent_above = 100 * count_above / total_count if total_count > 0 else 0

    # Print report
    print(f"\n=== {metric['title']} ===")
    print(f"Total Values: {total_count}")
    print(f"Mean: {mean_val:.2f}")
    print(f"Median: {median_val:.2f}")
    print(f"Min: {min_val:.2f}")
    print(f"Max: {max_val:.2f}")
    print(f"Above Threshold ({thresh}): {count_above} ({percent_above:.1f}%)")
    
    #end of the statistics of the matrix

    bounds = metric["bounds"]
    cmap = ListedColormap(['#ffffff'] + sns.color_palette("viridis", len(bounds) - 2))
    norm = BoundaryNorm(boundaries=bounds, ncolors=len(bounds) - 1)
    annot_matrix = matrix.applymap(lambda x: f"{x:{metric['fmt']}}" if x >= metric["annot_thresh"] else "")

    ax = axes[i]
    sns.heatmap(
        matrix,
        annot=annot_matrix,
        fmt="",
        cmap=cmap,
        norm=norm,
        ax=ax,
        cbar=True,
        linewidths=0.25,
        linecolor='gray',
        annot_kws={"size": 6},
        cbar_kws={"ticks": bounds}
    )

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=5, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=5, rotation=0)
    ax.set_title(metric["title"], fontsize=8)
    ax.set_xlabel("Cohesins", fontsize=6)
    ax.set_ylabel("Dockerins", fontsize=6)

    # Format colorbar
    cbar = ax.collections[0].colorbar
    cbar.outline.set_edgecolor('gray')
    cbar.outline.set_linewidth(0.25)
    cbar.ax.tick_params(labelsize=6)

plt.suptitle("Coh-Doc Interactions in $\\it{C. perfringens}$", fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.96])
for ax in axes:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.7)
        spine.set_edgecolor("gray")

labels = ['A', 'B', 'C', 'D', 'E', 'F']

for label, ax in zip(labels, axes):
    ax.text(
        -0.15, 1.05, label, transform=ax.transAxes,
        fontsize=12, fontweight='bold', va='top', ha='right'
    )
plt.savefig("Cper_full_seq_heatmaps.jpg", dpi=300)
plt.savefig("Cper_full_seq_heatmaps.svg", dpi=300)
plt.savefig("Cper_full_seq_heatmaps.pdf", dpi=300)
# plt.show()
