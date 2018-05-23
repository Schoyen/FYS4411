import sparse
import pickle
import os
import numba
from .coulomb_interface import (
        get_coulomb_element, get_energy, _get_antisymmetrized_elements,
        _get_coulomb_elements
)
from .index_map import get_indices_nm

ORBITAL_INTEGRALS = None

@numba.njit(cache=True)
def spin_delta(p, q):
    return not ((p & 0x1) ^ (q & 0x1))

def get_coulomb_elements(l: int, filename=None, tol=1e-8) -> sparse.COO:
    global ORBITAL_INTEGRALS

    if filename is not None:
        try:
            with open(filename, "rb") as f:
                ORBITAL_INTEGRALS = pickle.load(f)
                return ORBITAL_INTEGRALS
        except FileNotFoundError:
            pass

    indices = {p: get_indices_nm(p) for p in range(l//2)}
    ORBITAL_INTEGRALS = sparse.COO(
            *_get_coulomb_elements(l, indices, tol=tol),
            shape=(l//2, l//2, l//2, l//2))

    if filename is not None:
        with open(filename, "wb") as f:
            pickle.dump(ORBITAL_INTEGRALS, f)

    return ORBITAL_INTEGRALS

def get_antisymmetrized_elements(l: int, oi=ORBITAL_INTEGRALS,
        filename=None, tol=1e-8) -> sparse.COO:

    if oi is None:
        oi = get_coulomb_elements(l, filename=filename, tol=tol)

    u = sparse.COO(
            *_get_antisymmetrized_elements(l, oi.todense(), tol=tol),
            shape=(l, l, l, l))

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
