import cython
#from cython.parallel import prange

import numpy as np
cimport numpy as np

from libc.math cimport fabs

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def amplitude_scaling_two_body(
        np.ndarray[double, ndim=4] t,
        np.ndarray[double, ndim=2] h, int m, int n, double tol=1e-10):
    cdef int a, b, i, j
    cdef double divisor, val

    for a in range(m):
        for b in range(a, m):
            for i in range(n):
                for j in range(i, n):
                    divisor = h[i, i] + h[j, j] \
                            - h[a + n, a + n] - h[b + n, b + n]

                    if fabs(divisor) < tol:
                        continue

                    val = t[a, b, i, j] / divisor

                    t[a, b, i, j] = val
                    t[a, b, j, i] = -val
                    t[b, a, i, j] = -val
                    t[b, a, j, i] = val

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def amplitude_scaling_two_body_sparse(
        np.ndarray[unsigned char, ndim=2] indices,
        np.ndarray[double, ndim=1] data,
        np.ndarray[double, ndim=2] h, int n, double tol=1e-10):
    cdef int a, b, i, j, index, length
    cdef double divisor
    cdef np.ndarray[unsigned char, ndim=1] a_arr, b_arr, i_arr, j_arr

    a_arr = indices[0]
    b_arr = indices[1]
    i_arr = indices[2]
    j_arr = indices[3]

    length = len(data)

    for index in range(length):
        a = a_arr[index]
        b = b_arr[index]
        i = i_arr[index]
        j = j_arr[index]

        divisor = h[i, i] + h[j, j] - h[n + a, n + a] - h[n + b, n + b]

        if fabs(divisor) < tol:
            continue

        data[index] = data[index]/divisor
