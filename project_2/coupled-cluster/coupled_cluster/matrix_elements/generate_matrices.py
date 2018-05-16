import sparse
import pickle
import numba
from .coulomb_interface import (
        get_coulomb_element, get_energy, _get_antisymmetrized_elements
)
from .index_map import get_indices_nm

ORBITAL_INTEGRALS = None

@numba.njit(cache=True)
def spin_delta(p, q):
    return not ((p & 0x1) ^ (q & 0x1))

def get_coulomb_elements(l: int, filename="") -> sparse.COO:
    global ORBITAL_INTEGRALS

    ORBITAL_INTEGRALS = sparse.DOK((l//2, l//2, l//2, l//2))

    indices = {p: get_indices_nm(p) for p in range(l//2)}

    for p in range(l//2):
        for q in range(l//2):
            for r in range(l//2):
                for s in range(l//2):
                    element = get_coulomb_element(
                            *indices[p],
                            *indices[q],
                            *indices[r],
                            *indices[s]
                    )

                    if abs(element) < 1e-8:
                        continue

                    ORBITAL_INTEGRALS[p, q, r, s] = element

    ORBITAL_INTEGRALS = ORBITAL_INTEGRALS.to_coo()
    if filename:
        with open(filename, "w") as f:
            pickle.dump(ORBITAL_INTEGRALS, f)

    return ORBITAL_INTEGRALS

def get_antisymmetrized_elements(l: int, oi=ORBITAL_INTEGRALS,
        filename="") -> sparse.COO:

    if oi is None:
        oi = get_coulomb_elements(l)

    u = sparse.COO(
            *_get_antisymmetrized_elements(l, oi.todense()), shape=(l, l, l, l))

    if filename:
        with open(filename, "w") as f:
            pickle.dump(u, f)

    return u

def get_one_body_elements(l: int) -> sparse.COO:
    h = sparse.DOK((l//2, l//2))

    for p in range(l//2):
        h[p, p] = get_energy(*get_indices_nm(p))

    return h.to_coo()

def add_spin_to_one_body_elements(hi, l):
    h = sparse.DOK((l, l))

    for p in range(l):
        for q in range(l):
            h[p, q] = spin_delta(p, q) * hi[p//2, q//2]

    return h.to_coo()

def get_one_body_elements_spin(l: int) -> sparse.COO:
    h = sparse.DOK((l, l))
    _h = get_one_body_elements(l)

    for p in range(l):
        h[p, p] = _h[p//2, p//2]

    return h.to_coo()
