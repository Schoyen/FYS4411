from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np
import os

os.environ["CFLAGS"] = "-g -std=c++11 -stdlib=libc++"

base_path = ["vmc"]
source_path = base_path + ["src"]

config_path = source_path + ["config"]
wavefunctions_path = source_path + ["wavefunctions"]
hamiltonians_path = source_path + ["hamiltonians"]
math_path = source_path + ["math"]
solvers_path = source_path + ["solvers"]


source_files = [
        os.path.join(*base_path, "interface.pyx")
]


include_dirs = [
        os.path.join(*config_path),
#        wavefunctions_path,
#        hamiltonians_path,
#        math_path,
#        solvers_path,
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
