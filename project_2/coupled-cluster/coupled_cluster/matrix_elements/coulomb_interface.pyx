import numpy as np
cimport numpy as np

np.import_array()

def get_coulomb_element(p, q, r, s):
    n_p, m_p = p
    n_q, m_q = q
    n_r, m_r = r
    n_s, m_s = s

    return coulomb_ho(n_p, m_p, n_q, m_q, n_r, m_r, n_s, m_s)
