import pandas as pd
from .coulomb_interface import get_energy

INDEX_MAP = None

def _generate_index_map(n_range: int, m_range: int) -> None:
    global INDEX_MAP

    INDEX_MAP = pd.DataFrame(
            index=[(n, m) for n in range(n_range)
                for m in range(-m_range, m_range + 1)],
            columns=["energy", "n", "m", "p"])
    INDEX_MAP.n = list(map(lambda x: x[0], INDEX_MAP.index))
    INDEX_MAP.m = list(map(lambda x: x[1], INDEX_MAP.index))
    INDEX_MAP.energy = list(map(lambda x: get_energy(*x), INDEX_MAP.index))
    INDEX_MAP = INDEX_MAP.sort_values(by=["energy", "m"])
    INDEX_MAP.p = list(range(len(INDEX_MAP)))

def get_index_p(n: int, m: int) -> int:
    if INDEX_MAP is None:
        _generate_index_map(n + 1, m + 1)

    n_mask = INDEX_MAP.n == n
    if not any(n_mask):
        _generate_index_map(n + 1, m + 1)
        n_mask = INDEX_MAP.n == n

    m_mask = INDEX_MAP.m == m
    if not any(m_mask):
        _generate_index_map(n + 1, m + 1)
        m_mask = INDEX_MAP.m == m

    return int(INDEX_MAP[n_mask & m_mask].p)

def get_indices_nm(p: int) -> int:
    if INDEX_MAP is None:
        _generate_index_map(p//2 + 1, p + 1)

    mask = INDEX_MAP.p == p
    if not any(mask):
        _generate_index_map(p//2 + 1, p + 1)
        mask = INDEX_MAP.p == p

    return int(INDEX_MAP[mask].n), int(INDEX_MAP[mask].m)
