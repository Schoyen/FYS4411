import numpy as np
import sparse

def transform_one_body_elements(h, c):
    _h = np.einsum("ip, jq, ij -> pq", c.conj(), c, h)

    _h[np.abs(_h) < 1e-8] = 0

    return sparse.COO.from_numpy(_h)

def transform_two_body_elements(u, c):
    _u = np.einsum("pi, qj, kr, ls, ijkl -> pqrs", c.conj(), c.conj(), c, c, u)

    _u[np.abs(_u) < 1e-8] = 0

    return sparse.COO.from_numpy(_u)
