from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

source_files = []


include_dirs = [
        "include",
        np.get_include()
]


libraries = []

extensions = []

setup(
    name="Bosonic vmc",
    version="0.0.1",
    ext_modules=cythonize(extensions)
)
