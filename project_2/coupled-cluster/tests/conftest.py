import pytest
import numpy as np
import sparse

n = 2
l = 12
filename = "tests/dat/coulomb.dat"

def spin_delta(p, q):
    return not ((p & 0x1) ^ (q & 0x1))

def pytest_namespace():
    return {"n": n, "l": l, "orbital_integrals": orbital_integrals, "u": u}

def get_file_orbital_integrals():
    orbital_integrals = sparse.DOK((l//2, l//2, l//2, l//2))

    with open(filename, "r") as f:
        for line in f:
            line = line.split()

            if not line:
                continue

            p, q, r, s, val = line

            if abs(float(val)) < 1e-8:
                continue

            orbital_integrals[int(p), int(q), int(r), int(s)] = float(val)

    return orbital_integrals.to_coo()

orbital_integrals = get_file_orbital_integrals()

def get_file_antisymmetrized_integrals():
    u = sparse.DOK((l, l, l, l))

    for p in range(l):
        for q in range(l):
            for r in range(l):
                for s in range(l):
                    u_pqrs = spin_delta(p, r) * spin_delta(q, s) \
                            * orbital_integrals[p//2, q//2, r//2, s//2]
                    u_pqsr = spin_delta(p, s) * spin_delta(q, r) \
                            * orbital_integrals[p//2, q//2, s//2, r//2]

                    u_as = u_pqrs - u_pqsr
                    if abs(u_as) < 1e-8:
                        continue

                    u[p, q, r, s] = u_as

    u = u.to_coo()
    return u

u = get_file_antisymmetrized_integrals()

@pytest.fixture
def index_map():
    i_map = {
        0: (0, 0),
        1: (0, -1),
        2: (0, 1),
        3: (0, -2),
        4: (1, 0),
        5: (0, 2)
    }

    return i_map
