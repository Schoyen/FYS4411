from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np
import os
import glob
import platform

os.environ["CFLAGS"] = "-std=c++11"

if platform.system() == "Darwin":
    os.environ["CFLAGS"] += " -stdlib=libc++"

base_path = ["coupled_cluster"]
interface_path = base_path + ["matrix_elements"]
source_path =  interface_path + ["src"]

source_files = [
    *glob.glob(os.path.join(*interface_path, "*.pyx")),
    *glob.glob(os.path.join(*source_path, "*.cpp"))
]

include_dirs = [
    os.path.join(*source_path),
    np.get_include()
]

libraries = []

define_macros = []

undef_macros = []

extensions = [
    Extension(
        name="coupled_cluster.matrix_elements.coulomb_interface",
        sources=source_files,
        language="c++",
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=define_macros,
        undef_macros=undef_macros
    )
]

setup(
    name="Coupled Cluster",
    version="0.0.1",
    packages=["coupled_cluster"],
    ext_modules=cythonize(extensions)
)
