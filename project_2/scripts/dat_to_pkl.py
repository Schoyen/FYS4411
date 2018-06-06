import numpy as np
import sparse
import os
import pickle

file_path = os.path.join("..", "dat")
read_file = os.path.join(file_path, "coulomb_{0}.dat")
write_file = os.path.join(file_path, "coulomb_{0}.pkl")

l = 90
read_file = read_file.format(l)
write_file = write_file.format(l)

numpy_from_file = np.loadtxt(read_file)

data = numpy_from_file[:, 4]
coords = np.transpose(numpy_from_file[:, 0:4])

sparse_to_file = sparse.COO(coords, data, shape=(l//2, l//2, l//2, l//2))

pickle.dump(sparse_to_file, open(write_file, "wb"))