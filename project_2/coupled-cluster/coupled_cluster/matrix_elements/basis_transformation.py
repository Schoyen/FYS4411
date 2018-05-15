import numpy as np
import sparse

def transform_one_body_elements(h, c):
    _h = np.einsum("pi, jq, ij -> pq", c.conj(), c, h)

    _h[np.abs(_h) < 1e-8] = 0

    return sparse.COO.from_numpy(_h)
