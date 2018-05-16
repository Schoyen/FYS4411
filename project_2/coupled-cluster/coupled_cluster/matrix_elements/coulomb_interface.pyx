import numpy as np
cimport numpy as np

np.import_array()

def get_energy(int n, int m):
    return 2*n + abs(m) + 1

def get_coulomb_element(int n_p, int m_p, int n_q, int m_q, int n_r, int m_r,
        int n_s, int m_s):

    return coulomb_ho(n_p, m_p, n_q, m_q, n_r, m_r, n_s, m_s)

cdef int spin_delta(int p, int q):
    return not ((p & 0x1) ^ (q & 0x1))

def _get_antisymmetrized_elements(int l, np.ndarray[double, ndim=4] oi,
        double tol=1e-8):
    cdef list data, indices
    cdef int p, q, r, s
    cdef double u_pqrs, u_pqsr, u_as

    data = []
    indices = [[], [], [], []]

    for p in range(l):
        for q in range(l):
            for r in range(l):
                for s in range(l):
                    u_pqrs = spin_delta(p, r) * spin_delta(q, s) \
                            * oi[p//2, q//2, r//2, s//2]
                    u_pqsr = spin_delta(p, s) * spin_delta(q, r) \
                            * oi[p//2, q//2, s//2, r//2]

                    u_as = u_pqrs - u_pqsr

                    if abs(u_as) < tol:
                        continue

                    indices[0].append(p)
                    indices[1].append(q)
                    indices[2].append(r)
                    indices[3].append(s)
                    data.append(u_as)

    return indices, data
