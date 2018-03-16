'''
This is supposed to find the "true" variance
by way of the blocking method.
Code is adapted from
'''

import pickle
import numpy as np
import time

def block(x):
    n = len(x)
    d = np.log2(n)
    s = np.zeros(int(d))
    gamma = np.zeros(int(d))
    mu = np.mean(x)
    start_time = time.time()

    # Auto-covariance and variances for each blocking transformation
    for i in range(int(d)):

        # n changes in length
        n = len(x)

        # Autocovariance of x
        gamma[i] = (1/n) * np.sum((x[0:n-1] - mu)*(x[1:n] - mu))

        # Variance of x
        s[i] = np.var(x)

        # We might get a situaion where the array is not easily split in two equal sizes
        x_1 = x[0::2] # Extracting all numbers at odd positions
        x_2 = x[1::2] # Numbers at even positions
        # If length is not equal, remove highest number
        if (len(x_1) > len(x_2)):
            x_1 = x_1[:-1]
        elif (len(x_2) > len(x_1)):
            x_2 = x_2[:-1]

        # Blocking transformation
        x = 0.5*(x_1 + x_2)

    # Test observator from theorem (chi^2-distributed)
    factor_1 = (gamma/s)**2
    factor_2 = 2**np.arange(1, d+1)
    # Do the same length check again
    if (len(factor_1) > len(factor_2)):
        factor_1 = factor_1[:-1]
    elif (len(factor_2) > len(factor_1)):
        factor_2 = factor_2[:-1]

    M = (np.cumsum((factor_1 * factor_2[::-1])[::-1]))[::-1]

    # Test percentiles
    q = np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272,\
            16.811894, 18.475307, 20.090235, 21.665994, 23.209251,\
            24.724970, 26.216967, 27.688250, 29.141238, 30.577914,\
            31.999927, 33.408664, 34.805306, 36.190869, 37.566235,\
            38.932173, 40.289360, 41.638398, 42.979820, 44.314105,\
            45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # The actual Chi squared test - should we have stopped blocking?
    for k in range(0, int(d)):
        if (M[k] < q[k]):
            break
        if (k >= d-1):
            print("More data is needed!")
    
    result = s[k] / 2**(d-k)

    print("Runtime: {:12.9f} seconds".format(time.time() - start_time))
    print("Mean: {:8.5f}, Iterations: {}, STD. {:8.5f}".format(mu, k, result**.5))

    return result

x = pickle.load(open("local_energies.p", "rb"))
block(x)
