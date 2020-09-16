""" Generic setup.py to compile each .pyx module. Build using
python cython_setup.py build_ext --inplace
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

extensions = cythonize([
    Extension(
        "*",
        sources=["*.pyx"],
        include_dirs=[numpy.get_include()],
        define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        extra_compile_args = ["-ffast-math", "-O3", "-Wno-unreachable-code", "-Wunused-function"],
        libraries=["m"],
    ),
])

setup(
    ext_modules = extensions,
)
