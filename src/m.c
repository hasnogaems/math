#include <float.h>
#include <math.h>
#include <s21_math.h>
#include <stdio.h>
#include <stdlib.h>
double exp_hypergeom(double x) {
    int n = 0;
    double term = x;
    double result = 0.0;
    while (term > 1e-15) {  // условие остановки
        result += term;
        n++;
        term *= x / n;
    }
    return result;
}
int main() {
  for (double k = -2000; k < 1000; k += 1) {
     double a = s21_exp(k);
//     if(k<-20){
// a = exp_hypergeom(k);
//     }
    double b = exp(k);
    //printf("k=%.7f\n", k);
    if (fabs(a - b) > 1e-8) {
        printf("k=%.7f\n", k);
      //printf("s21=%.7f\n sm=%.7f   k=%.7f pres=%.7f\n", a, b, k, a - b);
      // break;
    }
  }
}
