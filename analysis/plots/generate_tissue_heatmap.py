import numpy as np
import collections
import matplotlib.pyplot as plt
import matplotlib
plt.switch_backend('agg')
import os
import sys
import json
sys.path.append('../')
from utils import utils
from collections import defaultdict
from mpl_toolkits.axes_grid1 import AxesGrid

empirical_percents = [5, 1, 0.05]
sbs = ["0.01", "0.05"]
methods = ["tejaas_perm", "matrixeqtl"]

tissue_file = "/usr/users/fsimone/trans-eqtl-pipeline/analysis/plots/tissues.txt"
tissues, descriptions = utils.read_tissues(tissue_file)

json_file = "../gtex_metadata.json"
with open(json_file) as instream:
    gtex_meta = json.load(instream)
my_colors = []
for d in descriptions:
    my_colors.append("#"+gtex_meta[d.replace(" ", "_")]["colorHex"])


# datasets = ['gtex-wb', 'gtex-sse', 'gtex-thy', 'gtex-aa', 'gtex-bca', 'gtex-bce']
datasets = ["gtex-"+t for t in tissues]
print("Working with:", datasets)
N = len(datasets)
matrix_dict = dict(zip(datasets, np.arange(N)))


def load_plot_data(auc_file, N, tissue_dict):
    auc_array = np.zeros((N,N))
    if os.path.exists(auc_file):
        with open(auc_file) as instream:
            for line in instream:
                t1  = line.split("\t")[0]
                t2  = line.split("\t")[1]
                auc = np.float64(line.split("\t")[2])
                auc1k = np.float64(line.split("\t")[3].rstrip())
                i = matrix_dict.get(t1, None)
                j = matrix_dict.get(t2, None)
                if i is None:
                    print("Oops, {:s} is not in the data".format(t1))
                    continue
                if j is None:
                    print("Oops, {:s} is not in the data".format(t2))
                    continue
                if i > j:
                    i,j = j,i
                auc_array[i,j] = auc
            #auc_dict[method] = auc_array
            # auc_dict[mkey] = auc_array
    return auc_array

data_dir = "/cbscratch/franco/trans-eqtl/analysis/data"
auc_dict = dict()
for ep in empirical_percents:
    for method in methods:
        if method == "tejaas_perm":
            for sb in sbs:
                method_sb = method + "_" + sb
                mkey = method_sb + "_ep" + str(ep)
                auc_dir = os.path.join(data_dir, method_sb, str(ep))
                auc_file = os.path.join(auc_dir, method_sb+".auc.txt")
                print(auc_file)
                auc_array = load_plot_data(auc_file, N, matrix_dict)
                auc_dict[mkey] = auc_array
        else:
            mkey = method + "_ep" + str(ep)
            auc_dir = os.path.join(data_dir, method_sb, str(ep))
            auc_file = os.path.join(auc_dir, method+".auc.txt")    
            print(auc_file)
            auc_array = load_plot_data(auc_file, N, matrix_dict)
            auc_dict[mkey] = auc_array

mykeys = list(auc_dict.keys())
print("Loaded:", mykeys)
# combiplots = [x for x in itertools.combinations(mykeys,2)]
# for m1, m2 in combiplots:
# m1 = mykeys[0]
# m2 = mykeys[1]
m1 = "tejaas_perm_0.01_ep5"
m2 = "tejaas_perm_0.05_ep5"

outfile = m1+"_"+m2+"_crossval_map_"+str(ep)+".png"
method_array = auc_dict[m1].T + auc_dict[m2] + np.identity(N)
method_array[range(N), range(N)] = (method_array.min() - 0.01)

fig = plt.figure(figsize = (15, 15))
ax = fig.add_subplot(111)

# orig_cmap = matplotlib.cm.seismic
# shifted_cmap = shiftedColorMap(orig_cmap, start=0.4, midpoint=0.7, stop=0.8, name='shifted')

cmap = plt.cm.RdBu_r
norm = plt.Normalize(method_array.min(), method_array.max())
rgba = cmap(norm(method_array))
rgba[range(N), range(N), :3] = 0.7,0.7,0.7
im = ax.imshow(rgba)

# Add the colorbar using a fake (not shown) image
im = ax.imshow(method_array, visible=False, cmap="RdBu_r")

# Create colorbar
#cbar = ax.figure.colorbar(im, ax=ax)

# plt.colorbar(im,fraction=0.046, pad=0.04)
cbar = ax.figure.colorbar(im,fraction=0.046, pad=0.04)
cbar.ax.set_ylabel("AUC", rotation=-90, va="bottom")




# delete axis lines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y',which='both',left=False)
ax.tick_params(axis='x',which='both',bottom=False)

ax.set_xticks(np.arange(N))
ax.set_yticks(np.arange(N))
ax.set_xticklabels(datasets, fontsize=14, color="black")
ax.set_yticklabels(["  " for x in datasets], fontsize=14, color="black")
bbox = dict(boxstyle="round", alpha=0.9, ec="white", fc="white")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", va="center", fontsize=12,
        rotation_mode="anchor", bbox=bbox)
plt.setp(ax.get_yticklabels(), bbox=bbox)

# ax.tick_params(axis='x', colors='red')
# my_colors = ['#FFF000', 'b', 'r', 'r', 'g', 'k']
for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), my_colors):
    ticklabel.set_backgroundcolor(tickcolor)
for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), my_colors):
    ticklabel.set_backgroundcolor(tickcolor)

# for i in range(N):
#     for j in range(N):
#         text = ax.text(j, i, "{:.3f}".format(method_array[i, j]).replace("{:.3f}".format(method_array.min()), ""),
#                        ha="center", va="center", color="w", fontsize=14)

# ax.set_title("Tissue-Tissue validation AUC", fontsize=14)
fig.tight_layout()
plt.savefig(outfile, bbox_inches='tight')
plt.show()