""" Generic setup.py to compile c module. Build using
python c_setup.py build_ext --inplace
"""

from distutils.core import setup, Extension
import numpy as np

def configure():
    ext_modules = [Extension('c_lib/c_functions',
                    sources = ['c_functions.c'],
                    extra_compile_args=["-ffast-math", "-O3"])]

    setup(
            include_dirs = [np.get_include()],
            ext_modules = ext_modules
          )

if __name__ == '__main__':
    configure()
