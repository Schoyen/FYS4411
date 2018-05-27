import pytest

from coupled_cluster.matrix_elements.index_map import (
        get_index_p, get_indices_nm, generate_index_map, IndexMap
)

def test_index_map(index_map):
    num_shells = pytest.large_num_shells
    generate_index_map(num_shells)

    print (num_shells)
    print (IndexMap.p_list)
    print (IndexMap.nm_dict)
    print (IndexMap.shell_arr)

    for p in index_map:
        assert index_map[p] == get_indices_nm(p)

    for p, (n, m) in index_map.items():
        assert p == get_index_p(n, m)

    assert False
