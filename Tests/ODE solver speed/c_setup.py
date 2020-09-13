""" Generic setup.py to compile c module. Build using
python cython_setup.py build_ext --inplace
"""

from distutils.core import setup, Extension

setup(name             = "fast_c_functions",
      ext_modules      = [Extension('c_functions', ['c_functions.h'],
                                    extra_compile_args=["-fastmath", "-O3", "-march=native"])],
)
