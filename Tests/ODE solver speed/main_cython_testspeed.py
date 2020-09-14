# Run me with:
# python3 cython_setyp.py build_ext --inplace && python3 main_cython_testspeed.py

# See https://github.com/rochus/cython-odesolver-nD-shootout

from main_cython_testspeed import main
from main_cython_testspeed import testtiming_cython
from main_cython_testspeed import testtiming

if __name__ == "__main__":
    testtiming()
