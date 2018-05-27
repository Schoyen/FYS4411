import numpy as np

class IndexMap:
    p_list = []
    nm_dict = {}
    shell_arr = None

def generate_index_map(num_shells: int) -> None:
    IndexMap.p_list = []
    IndexMap.nm_dict = {}
    IndexMap.shell_arr = np.cumsum(range(1, num_shells + 1)) * 2

    counter = 0

    for shell in range(num_shells):
        n_max = shell // 2
        n = 0
        decrease = False

        for m in range(-shell, shell + 1, 2):
            IndexMap.p_list.append((n, m))
            IndexMap.nm_dict[(n, m)] = counter
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
        return IndexMap.nm_dict[(n, m)]
    except KeyError:
        raise KeyError((
            "Index map does not contain a p-value for " + \
                    "(n = {0}, m = {1})").format(n, m))

def get_indices_nm(p: int) -> int:
    try:
        return IndexMap.p_list[p]
    except IndexError:
        raise IndexError(
            "Index map does not contain a (n, m)-tuple for p = {0}".format(p))
