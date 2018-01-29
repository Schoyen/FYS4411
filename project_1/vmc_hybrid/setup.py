from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

source_files = [
    "vmc/particles.pyx",
    "src/metropolis_sampling.c",
    #"src/bosonic_hard_sphere.c",
    "src/one_dimensional_ho.c"
]


include_dirs = [
        "include",
        np.get_include()
]


libraries = []

extensions = [
    Extension(
        name="particles",
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
