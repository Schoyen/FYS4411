from coupled_cluster.matrix_elements.index_map import (
        get_index_p, get_indices_nm
)

def test_index_map(index_map):
    for p in index_map:
        assert index_map[p] == get_indices_nm(p)

    for p, (n, m) in index_map.items():
        assert p == get_index_p(n, m)
