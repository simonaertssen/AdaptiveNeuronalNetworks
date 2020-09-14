#include <stdio.h>
#include <math.h>

double sumtest(double *theta, int length){
  int i = 0;
  double sum = 0;
  for (i = 0; i < length; ++i){
    sum += *theta;
    theta++;
  }
  return sum;
}
