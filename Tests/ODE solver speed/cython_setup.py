""" Generic setup.py to compile each .pyx module. Build using
python cython_setup.py build_ext --inplace
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

extensions = cythonize([
    Extension(
        "*",
        sources=["*.pyx"],
        include_dirs=[np.get_include()],
        define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        extra_compile_args = ["-ffast-math", "-O3", "-Wno-unreachable-code", "-Wunused-function"],
        libraries=["m"]
    ),
], compiler_directives={'language_level' : "3"})

setup(
    ext_modules = extensions,
)
