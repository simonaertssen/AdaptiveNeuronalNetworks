#include <stdio.h>
#include <math.h>

double factorial(int n){
  if (n==1 || n==0) return 1;
  else return n * factorial(n-1);
}

double c_a_n(int n){
  return pow(2,n)*pow(factorial(n),2)/factorial(2*n);
}


double DOPRI54_step()


double pytorch_pulse(theta):
    return pow(1 - cos(theta), 2)
