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
  int check = 1;
  for (double k = -400; k <= 400; k += 1.7) {
    for (double g = -400; g < 400; g += 1) {
      double a = s21_pow(k, g);
      double b = pow(k, g);
      if (fabsl(a - b) > 1e-6) {
        printf("k=%.7f\ng=%.7f\n21=%.7f\nm=%.7f\n", k, g, a, b);
        check = 0;
        break;
      }
      if (!check) {
        break;
      }
    }
  }
}