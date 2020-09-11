from numba import njit

@njit(fastmath=True)
def factorial(n):
    if n==1 or n==0:
        return 1
    else:
        return n * factorial(n-1)

@njit(fastmath=True)
def a_n(n):
    return (2**n)*(factorial(n)**2)/factorial(2*n)
