import numba
import numpy as np



@numba.njit(nogil=True, parallel=True)
def compute_chi_abcd(chi_abcd, t, u, n_size, m_size):
    for a in numba.prange(m_size):
        a_virt = a + n_size
        for b in range(a, m_size):
            b_virt = b + n_size
            for c in range(m_size):
                c_virt = c + n_size
                for d in range(c, m_size):
                    d_virt = d + n_size

                    val = 0
                    for m in range(n_size):
                        for n in range(n_size):
                            val += t[a, b, m, n] * u[m, n, c_virt, d_virt]

                    val = 0.25 * val + 0.5 * u[a_virt, b_virt, c_virt, d_virt]

                    chi_abcd[a, b, c, d] = val
                    chi_abcd[b, a, c, d] = -val
                    chi_abcd[a, b, d, c] = -val
                    chi_abcd[b, a, d, c] = val

@numba.njit(nogil=True, parallel=True)
def compute_chi_ad(chi_ad, t, u, n_size, m_size):
    for a in numba.prange(m_size):
        for d in range(m_size):
            d_virt = d + n_size

            val = 0
            for c in range(m_size):
                c_virt = c + n_size
                for n in range(n_size):
                    for m in range(n_size):
                        val += t[a, c, n, m] * u[n, m, c_virt, d_virt]

            chi_ad[a, d] = 0.5 * val

@numba.njit(nogil=True, parallel=True)
def compute_chi_bmjc(chi_bmjc, t, u, n_size, m_size):
    for b in numba.prange(m_size):
        b_virt = b + n_size
        for m in range(n_size):
            for j in range(n_size):
                for c in range(m_size):
                    c_virt = c + n_size

                    val = 0
                    for d in range(m_size):
                        d_virt = d + n_size
                        for n in range(n_size):
                            val += t[b, d, j, n] * u[m, n, c_virt, d_virt]

                    chi_bmjc[b, m, j, c] = 0.5 * val + u[b_virt, m, j, c_virt]

@numba.njit(nogil=True, parallel=True)
def compute_chi_nj(chi_nj, t, u, n_size, m_size):
    for n in numba.prange(n_size):
        for j in range(n_size):

            val = 0
            for c in range(m_size):
                c_virt = c + n_size
                for d in range(m_size):
                    d_virt = d + n_size
                    for m in range(n_size):
                        val += t[c, d, j, m] * u[m, n, c_virt, d_virt]

            chi_nj[n, j] = 0.5 * val

@numba.njit(nogil=True, parallel=True)
def compute_chi_abcd_contraction(term, t, chi_abcd, n_size, m_size):
    for a in numba.prange(m_size):
        for b in range(a, m_size):
            for i in range(n_size):
                for j in range(i, n_size):

                    val = 0
                    for c in range(m_size):
                        for d in range(m_size):
                            val += t[c, d, i, j] * chi_abcd[a, b, c, d]

                    term[a, b, i, j] = val
                    term[b, a, i, j] = -val
                    term[a, b, j, i] = -val
                    term[b, a, j, i] = val

@numba.njit(nogil=True, parallel=True)
def compute_chi_bmjc_contraction(term, t, chi_bmjc, n_size, m_size):
    for a in numba.prange(m_size):
        for b in range(a, m_size):
            for i in range(n_size):
                for j in range(i, n_size):

                    val = 0
                    for c in range(m_size):
                        for m in range(n_size):
                            val += t[a, c, i, m] * chi_bmjc[b, m, j, c]
                            val -= t[b, c, i, m] * chi_bmjc[a, m, j, c]
                            val -= t[a, c, j, m] * chi_bmjc[b, m, i, c]
                            val += t[b, c, j, m] * chi_bmjc[a, m, i, c]

                    term[a, b, i, j] = val
                    term[b, a, i, j] = -val
                    term[a, b, j, i] = -val
                    term[b, a, j, i] = val

@numba.njit(nogil=True, parallel=True)
def compute_t_u_contraction(term, t, u, n_size, m_size):
    for a in numba.prange(m_size):
        for b in range(a, m_size):
            for i in range(n_size):
                for j in range(i, n_size):

                    val = 0
                    for m in range(n_size):
                        for n in range(n_size):
                            val += t[a, b, m, n] * u[m, n, i, j]

                    term[a, b, i, j] = 0.5 * val
                    term[b, a, i, j] = -0.5 * val
                    term[a, b, j, i] = -0.5 * val
                    term[b, a, j, i] = 0.5 * val
