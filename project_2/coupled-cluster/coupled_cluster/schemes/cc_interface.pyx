import numpy as np
cimport numpy as np

def _amplitude_scaling_two_body(np.ndarray[int, ndim=2, mode="c"] indices,
        np.ndarray[double, ndim=1, mode="c"] data,
        np.ndarray[double, ndim=2, mode="c"] h, int n, double tol=1e-10):
    cdef int a, b, i, j, index
    cdef double divisor

    for index in range(len(data)):
        a, b, i, j = indices[index]

        divisor = h[i, i] + h[j, j] - h[n + a, n + a] - h[n + b, n + b]

        if abs(divisor) < tol:
            continue

        data[index] = data[index]/divisor
