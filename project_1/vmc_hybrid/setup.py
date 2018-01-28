from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

source_files = [
    "vmc/particles.pyx"
]


include_dirs = [
        "include",
        np.get_include()
]


libraries = []

extensions = [
    Extension(
        name="vmc_hybrid",
        sources=source_files,
        language="c",
        include_dirs=include_dirs,
        libraries=libraries
    )
]

setup(
    name="Bosonic vmc",
    version="0.0.1",
    ext_modules=cythonize(extensions)
)
