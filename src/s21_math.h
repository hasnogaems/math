#ifndef S21_MATH_H
#define S21_MATH_H
#include <stdbool.h>
#define s21_INFINITY (1.0 / 0.0)
#define s21_NEGINFINITY (-1.0 / 0.0)
#include <float.h>
#include <stdio.h>

// 1
extern const double s21_M_E;
extern const double s21_M_PI;

int s21_abs(int x);

// 2
long double s21_acos(double x);

// 3
long double s21_asin(double x);

// 4
long double s21_atan(double x);

// 5
long double s21_ceil(double x);

// 6
long double s21_cos(double x);

// 7
long double s21_exp(double x);

// 8
long double s21_fabs(double x);

// 9
long double s21_floor(double x);

// 10
long double s21_fmod(double x, double y);

// 11
long double s21_log(double x);

// 12
long double s21_pow(double base, double exp);

// 13
long double s21_sin(double x);

// 14
long double s21_sqrt(double x);

// 15
long double s21_tan(double x);


int s21_isnan(long double x);
int s21_isinf(long double x);
long double s21_sqrtx(double base, double exp) ;
#endif  // MATHLIB_H
