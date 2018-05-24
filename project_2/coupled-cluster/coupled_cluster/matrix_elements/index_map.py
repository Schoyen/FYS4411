import numpy as np

P_LIST = []
NM_DICT = {}
SHELL_ARR = None

def generate_index_map(num_shells: int) -> None:
    global P_LIST, NM_DICT, SHELL_ARR

    P_LIST = []
    NM_DICT = {}
    SHELL_ARR = np.cumsum(range(1, num_shells + 1)) * 2

    counter = 0
    decrease = False

    for shell in range(num_shells + 1):
        n_max = shell // 2
        n = 0

        for m in range(-shell, shell + 1, 2):
            P_LIST.append((n, m))
            NM_DICT[(n, m)] = counter
            counter += 1

            if not decrease and n < n_max:
                n += 1
            elif not decrease and n == n_max:
                decrease = True
                if shell % 2 == 0:
                    n -= 1
            elif decrease:
                n -= 1
            else:
                raise Exception((
                    "Invalid value: decrease = {0}, n = {1}, " \
                            + "n_max = {2}").format(decrease, n, n_max))

def get_index_p(n: int, m: int) -> int:
    try:
        return NM_DICT[(n, m)]
    except KeyError:
        raise KeyError((
            "Index map does not contain a p-value for " + \
                    "(n = {0}, m = {1})").format(n, m))

def get_indices_nm(p: int) -> int:
    try:
        return P_LIST[p]
    except IndexError:
        raise IndexError(
            "Index map does not contain a (n, m)-tuple for p = {0}".format(p))
