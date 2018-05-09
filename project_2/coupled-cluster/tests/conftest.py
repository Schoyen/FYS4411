import pytest
import numpy as np
import sparse

n = 2
l = 12
filename = "tests/dat/coulomb.dat"

def pytest_namespace():
    return {"n": n, "l": l, "orbital_integrals": orbital_integrals}

def get_orbital_integrals():
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

orbital_integrals = get_orbital_integrals()
