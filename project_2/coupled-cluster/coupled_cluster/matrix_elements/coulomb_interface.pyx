import numpy as np
cimport numpy as np

np.import_array()

def get_energy(int n, int m):
    return 2*n + abs(m) + 1

def get_coulomb_element(int n_p, int m_p, int n_q, int m_q, int n_r, int m_r,
        int n_s, int m_s):

    return coulomb_ho(n_p, m_p, n_q, m_q, n_r, m_r, n_s, m_s)
