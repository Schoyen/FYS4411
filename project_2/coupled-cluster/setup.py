from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

import numpy as np
import os
import glob
import platform

os.environ["CFLAGS"] = "-std=c++11"

if platform.system() == "Darwin":
    os.environ["CFLAGS"] += " -stdlib=libc++"

base_path = ["coupled_cluster"]

coulomb_interface_path = base_path + ["matrix_elements"]
coulomb_source_path =  coulomb_interface_path + ["src"]

cc_interface_path = base_path + ["schemes"]

source_files = {
    "coulomb_interface": [
        *glob.glob(os.path.join(*coulomb_interface_path, "*.pyx")),
        *glob.glob(os.path.join(*coulomb_source_path, "*.cpp"))
    ],
    "cc_interface": [
        *glob.glob(os.path.join(*cc_interface_path, "*.pyx"))
    ]
}

include_dirs = {
    "coulomb_interface": [
        os.path.join(*coulomb_source_path),
        np.get_include()
    ],
    "cc_interface": [
        np.get_include()
    ]
}

libraries = []

define_macros = []

undef_macros = []

extensions = [
    Extension(
        name="coupled_cluster.matrix_elements.coulomb_interface",
        sources=source_files["coulomb_interface"],
        language="c++",
        include_dirs=include_dirs["coulomb_interface"],
        libraries=libraries,
        define_macros=define_macros,
        undef_macros=undef_macros
    ),
    Extension(
        name="coupled_cluster.schemes.cc_interface",
        sources=source_files["cc_interface"],
        language="c",
        include_dirs=include_dirs["cc_interface"],
        libraries=libraries,
        define_macros=define_macros,
        undef_macros=undef_macros,
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"]
    )
]

setup(
    name="Coupled Cluster",
    version="0.0.1",
    packages=find_packages(),
    ext_modules=cythonize(extensions)
)
