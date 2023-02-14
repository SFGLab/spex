# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
# download test data
# this file is 145 Mb, and may take a few seconds to download
import cooltools
import cooler

clr = cooler.Cooler(f'ChIA_PIPE_GM12878_RNAPOL2_4DNFIXG4BEEL.mcool::resolutions/100000')
chr = "chr17"
chromstarts = []
for i in clr.chromnames:
    print(f'{i} : {clr.extent(i)}')
    chromstarts.append(clr.extent(i)[0])

from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')

def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

f, ax = plt.subplots(
    figsize=(7,6))



im = ax.matshow(
    clr.matrix(balance=False).fetch(chr),
    vmax=5,
    extent=(0,clr.chromsizes[chr], clr.chromsizes[chr], 0),
    cmap="Reds"
)
plt.colorbar(im, ax=ax ,fraction=0.046, pad=0.04, label='raw counts')
ax.set_title(chr, y=1.08)
ax.set_ylabel('position, Mb')
format_ticks(ax)

plt.savefig("heatmap.svg", dpi=400)