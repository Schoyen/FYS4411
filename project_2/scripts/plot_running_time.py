import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import os

sns.set(color_codes=True)

file_path = os.path.join("..", "dat")
filename = "iter_run_times_{0}.dat"

fig_path = os.path.join(file_path, "figures")
fig_name = "iter_run_times_{0}.pdf"

n = 6
filename = os.path.join(file_path, filename.format(n))
fig_name = os.path.join(fig_path, fig_name.format(n))

df = pd.read_csv(filename, sep=" *, *", engine="python")
df.index = df["R"]
del df["R"]

df.plot()
plt.xlabel(r"$R$", fontsize=12)
plt.ylabel(r"$t$ [s]", fontsize=12)
plt.title(
        r"Time spent per CCD iteration for $N = {0}$ particles".format(n),
        fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(fig_name)
plt.show()

zoom_fig_name = "iter_run_times_zoom_{0}.pdf"
zoom_fig_name = os.path.join(fig_path, zoom_fig_name.format(n))

del df["CoupledClusterDoubles"]

df.plot()
plt.xlabel(r"$R$", fontsize=12)
plt.ylabel(r"$t$ [s]", fontsize=12)
plt.title(
        r"Sparse and optimized CCD for $N = {0}$ particles".format(n),
        fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(zoom_fig_name)
plt.show()
