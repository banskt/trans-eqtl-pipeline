import numpy as np
import collections
import matplotlib.pyplot as plt
import matplotlib
plt.switch_backend('agg')
import os
import sys
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
tissues, descriprtions = utils.read_tissues(tissue_file)

datasets = ["gtex-"+t for t in tissues]
print("Working with:", datasets)
N = len(datasets)
matrix_dict = dict(zip(datasets, np.arange(N)))

methods = ['tejaas_perm', 'tejaas_rand_perm', 'matrixeqtl', 'matrixeqtl_fdr', 'matrixeqtl_rand']



auc_dict = dict()
for method in methods:
    auc_file = method+".testauc.txt"
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

method_array = auc_dict[method1].T + auc_dict[method2] + np.identity(N)*0.4
fig = plt.figure(figsize = (8, 8))
ax = fig.add_subplot(111)

# orig_cmap = matplotlib.cm.seismic
# shifted_cmap = shiftedColorMap(orig_cmap, start=0.4, midpoint=0.7, stop=0.8, name='shifted')

im = ax.imshow(method_array, cmap="hot_r")
# Create colorbar
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel("AUC", rotation=-90, va="bottom")



ax.set_xticks(np.arange(N))
ax.set_yticks(np.arange(N))
ax.set_xticklabels(datasets, fontsize=14)
ax.set_yticklabels(datasets, fontsize=14)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
for i in range(N):
    for j in range(N):
        text = ax.text(j, i, "{:.3f}".format(method_array[i, j]).replace("0.500", ""),
                       ha="center", va="center", color="w", fontsize=14)

ax.set_title("Tissue-Tissue validation AUC", fontsize=14)
fig.tight_layout()
plt.savefig(outfile, bbox_inches='tight')
plt.show()