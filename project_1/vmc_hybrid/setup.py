from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np
import os
import platform


os.environ["CFLAGS"] = "-g -std=c++11"

if platform.system() == "Darwin":
    os.environ["CFLAGS"] += "-stdlib=libc++"

base_path = ["vmc"]
source_path = base_path + ["src"]

config_path = source_path + ["config"]
wavefunctions_path = source_path + ["wavefunctions"]
hamiltonians_path = source_path + ["hamiltonians"]
math_path = source_path + ["math"]
solvers_path = source_path + ["solvers"]


source_files = [
        os.path.join(*base_path, "interface.pyx"),
        os.path.join(*wavefunctions_path, "wavefunction.cpp"),
        os.path.join(*wavefunctions_path, "simple_gaussian.cpp"),
        os.path.join(*hamiltonians_path, "harmonic_oscillator.cpp"),
        os.path.join(*solvers_path, "monte_carlo_method.cpp"),
        os.path.join(*solvers_path, "metropolis_algorithm.cpp")
]


include_dirs = [
        os.path.join(*config_path),
        os.path.join(*wavefunctions_path),
        os.path.join(*hamiltonians_path),
        os.path.join(*math_path),
        os.path.join(*solvers_path),
        np.get_include()
]


libraries = []

extensions = [
    Extension(
        name="vmc.interface",
        sources=source_files,
        language="c++",
        include_dirs=include_dirs,
        libraries=libraries
    )
]

setup(
    name="Bosonic vmc",
    version="0.0.1",
    ext_modules=cythonize(extensions)
)
