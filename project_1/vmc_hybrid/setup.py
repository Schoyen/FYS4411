from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

particles_source_files = [
    "vmc/particles.pyx",
    "src/metropolis_sampling.c",
    #"src/bosonic_hard_sphere.c",
    "src/one_dimensional_ho.c"
]

wavefunction_source_files = [
    "vmc/wavefunction.pyx"
]


include_dirs = [
        "include",
        np.get_include()
]


libraries = []

extensions = [
    Extension(
        name="particles",
        sources=particles_source_files,
        language="c",
        include_dirs=include_dirs,
        libraries=libraries
    ),
    Extension(
        name="wavefunction",
        sources=wavefunction_source_files,
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
