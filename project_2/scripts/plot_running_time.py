import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import os

sns.set(color_codes=True)

file_path = os.path.join("..", "dat")
filename = "iter_run_times_{0}.dat"

fig_path = os.path.join(file_path, "figures")
fig_name = "iter_run_times_{0}.pdf"

n = 2
filename = os.path.join(file_path, filename.format(n))
fig_name = os.path.join(fig_path, fig_name.format(n))

df = pd.read_csv(filename, sep=" *, *", engine="python")
df.index = df["R"]
del df["R"]

df.plot()
plt.xlabel(r"$R$")
plt.ylabel(r"$t$ [s]")
plt.title(r"Time spent per CCD iteration")
plt.savefig(fig_name)
plt.show()
