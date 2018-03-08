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

base_path = ["vmc"]
source_path = base_path + ["src"]

source_files = [
        os.path.join(*base_path, "interface.pyx"),
        *glob.glob(os.path.join(*source_path, "*", "*.cpp"))
]


include_dirs = [
        *glob.glob(os.path.join(*source_path, "*")),
        np.get_include()
]


libraries = []

define_macros=[]

undef_macros=[
        "NDEBUG"
]

extensions = [
    Extension(
        name="vmc.interface",
        sources=source_files,
        language="c++",
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=define_macros,
        undef_macros=undef_macros
    )
]

setup(
    name="Bosonic vmc",
    version="0.0.1",
    packages=["vmc"],
    ext_modules=cythonize(extensions)
)
