import pickle
import tqdm

l = [90, 110, 132, 156]
filename = "coulomb_{0}.pkl"
out_file = "coulomb_{0}.dat"

data = []
for _l in l:
    with open(filename.format(_l), "rb") as f:
        data.append(pickle.load(f))

for i, _l in enumerate(l):
    print (_l)
    with open(out_file.format(_l), "w") as f:
        dat = data[i]
        _p, _q, _r, _s = dat.coords
        for p, q, r, s in tqdm.tqdm(zip(_p, _q, _r, _s)):
            f.write("{0} {1} {2} {3} {4}\n".format(p, q, r, s, dat[p, q, r, s]))
