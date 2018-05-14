import sparse
import pickle
import numba
from .coulomb_interface import get_coulomb_element, get_energy
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

def get_antisymmetrized_elements(l: int, filename="") -> sparse.COO:
    global ORBITAL_INTEGRALS

    u = sparse.DOK((l, l, l, l))

    if ORBITAL_INTEGRALS is None:
        ORBITAL_INTEGRALS = get_coulomb_elements(l)

    _oi = ORBITAL_INTEGRALS

    for p in range(l):
        for q in range(l):
            for r in range(l):
                for s in range(l):
                    u_pqrs = spin_delta(p, r) * spin_delta(q, s) \
                            * _oi[p//2, q//2, r//2, s//2]
                    u_pqsr = spin_delta(p, s) * spin_delta(q, r) \
                            * _oi[p//2, q//2, s//2, r//2]

                    u_as = u_pqrs - u_pqsr
                    if abs(u_as) < 1e-8:
                        continue

                    u[p, q, r, s] = u_as

    u = u.to_coo()
    if filename:
        with open(filename, "w") as f:
            pickle.dump(u, f)

    return u

def get_one_body_elements(l: int) -> sparse.COO:
    h = sparse.DOK((l, l))

    for p in range(l):
        h[p, p] = get_energy(*get_indices_nm(p//2))

    return h.to_coo()
