from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

source_files = [
    "vmc/metropolis_sampling.pyx",
    "src/wavefunction.c",
    "src/metropolis_sampling.c",
    "src/harmonic_oscillator.c"
]


include_dirs = [
        "include",
        np.get_include()
]


libraries = []

extensions = [
    Extension(
        name="metropolis_sampling",
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
