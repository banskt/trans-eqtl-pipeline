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

# code by paul-h
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

method1 = "tejaas_perm"
method2 = "matrixeqtl"
outfile = method1+"_"+method2+"_crossval_map.png"

tissue_file = "/usr/users/fsimone/trans-eqtl-pipeline/main/tissues.txt"
tissues, descriptions = utils.read_tissues(tissue_file)

json_file = "../gtex_metadata.json"
with open(json_file) as instream:
    gtex_meta = json.load(instream)


my_colors = []
for d in descriptions:
    my_colors.append("#"+gtex_meta[d.replace(" ", "_")]["colorHex"])


datasets = ["gtex-"+t for t in tissues]
print("Working with:", datasets)
N = len(datasets)
matrix_dict = dict(zip(datasets, np.arange(N)))

methods = ['tejaas_perm', 'matrixeqtl' ]#, 'matrixeqtl_fdr', 'matrixeqtl_rand', 'tejaas_rand_perm']

auc_dir = "/cbscratch/franco/trans-eqtl/analysis/data"

auc_dict = dict()
for method in methods:
    # auc_file = method+".testauc.txt"
    auc_file = os.path.join(auc_dir, method+".auc.txt")
    print(auc_file)
    auc_array = np.zeros((N,N))
    if os.path.exists(auc_file):
        with open(auc_file) as instream:
            for line in instream:
                t1  = line.split("\t")[0]
                t2  = line.split("\t")[1]
                auc = np.float64(line.split("\t")[2].rstrip())
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
            auc_dict[method] = auc_array

print("Loaded:", auc_dict.keys())

method_array = auc_dict[method1].T + auc_dict[method2] + np.identity(N)
method_array[range(N), range(N)] = (method_array.min() - 0.01)
fig = plt.figure(figsize = (12, 12))
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
cbar = ax.figure.colorbar(im, ax=ax)
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
ax.set_yticklabels(datasets, fontsize=14, color="black")
bbox = dict(boxstyle="round", alpha=0.9, ec="white", fc="white")
plt.setp(ax.get_xticklabels(), rotation=90, ha="right", va="center",
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

ax.set_title("Tissue-Tissue validation AUC", fontsize=14)
fig.tight_layout()
plt.savefig(outfile, bbox_inches='tight')
plt.show()