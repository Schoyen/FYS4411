import pandas as pd
from .coulomb_interface import get_energy

def _generate_index_map(num_shells):
    global INDEX_MAP
    INDEX_MAP = pd.DataFrame(index=range(l), columns=["energy", "n", "m"])
