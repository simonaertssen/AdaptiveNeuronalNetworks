cpdef factorial(int n) except -2:
    if n==1 or n==0:
        return 1
    else:
        return n * factorial(n-1)

cpdef cython_a_n(int n) except -1:
    return (2**n)*(factorial(n)**2)/factorial(2*n)

# def factorial(n):
#     if n==1 or n==0:
#         return 1
#     else:
#         return n * factorial(n-1)
#
# def cython_a_n(n):
#     return (2**n)*(factorial(n)**2)/factorial(2*n)
