import numpy as np
from matplotlib import pyplot as plt
import math
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

def f(x, factor):
    return np.exp(-factor*abs(x))

x_v = np.linspace(-10, 10, 1000)

i = 0
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])


ax.patch.set_edgecolor('black')  

ax.patch.set_linewidth(3)  

for factor in [1, 0.5, 0.2, 0.1]:
    plt.plot(x_v, [f(x, factor) for x in x_v], color='red', alpha=math.sqrt(factor))

plt.savefig("plot.png", dpi=400)