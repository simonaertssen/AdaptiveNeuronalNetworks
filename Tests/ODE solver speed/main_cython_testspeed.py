# Run me with:
# python3 cython_setyp.py build_ext --inplace && python3 main_cython_testspeed.py

from main_cython_testspeed import main
from main_cython_testspeed import testtiming_cython
from main_cython_testspeed import testtiming

if __name__ == "__main__":
    main()
