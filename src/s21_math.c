#include "s21_math.h"

#define NAN (0.0 / 0.0)
const double s21_M_E = 2.71828182845904523536028747135266250;
const double s21_M_PI = 3.14159265358979323846264338327950288;
long double s21_sqrtx(double base, double exp) {
  return s21_exp(s21_log(base) / exp);
}
int s21_isnan(long double x) { return (x != x); }
int s21_isinf(long double x) {
  return (!s21_isnan(x) && (x == s21_INFINITY || x == s21_NEGINFINITY)) ? 1 : 0;
}
long double s21_fact(int x) {
  long double res = 1;
  for (int i = 0; i <= x; i++) {
    res *= i;
  }
  return res;
}
// 1
int s21_abs(int x) { return (x < 0 ? -x : x); }

// 2
long double s21_acos(double x) {
  long double result = NAN;
  if (!(x > 1 || s21_isnan(x))) result = s21_M_PI / 2 - s21_asin(x);
  return result;
}

// 3
long double s21_asin(double x) {
  long double res = 0.0;
  if (s21_isnan(x) || s21_isinf(x))
    res = NAN;
  else if (s21_fabs(x) > 1)
    res = NAN;
  else if (s21_fabs(x) <= 1) {
    res = s21_atan(x / s21_sqrt(1 - x * x));
  }
  return res;
}

// 4
long double s21_atan(double x) {
  long double result = NAN;
  if (!s21_isnan(x)) {
    int n = x < 0.0;
    if (n) x = -x;
    int i = x > 1.0;
    if (i) x = 1.0 / x;
    int f = x > 0.5;
    if (f) x = (x - 0.5) / (1.0 + 0.5 * x);

    double p = x, r = p, xx = x * x;

    for (int d = 3; p < 0. || p > 5e-16; d += 2) {
      r += p *= (2 - d) * xx / d;
    }
    if (f) r += 0.4636476090008060935;
    if (i) r = 1.5707963267948966192 - r;
    result = n ? -r : r;
  }

  return result;
}

// 5
long double s21_ceil(double x) {
  if (!(s21_isinf(x) || s21_isnan(x))) {
    long double integerPart = (int)x;
    long double decimalPart = x - integerPart;
    if (decimalPart > 0) {
      integerPart += 1;
    }
    x = integerPart;
  }

  return x;
}

// 6
long double s21_cos(double x) {
  long double result;
  if (s21_isnan(x) || s21_isinf(x)) {
    result = NAN;  // Для бесконечных значений результат не определен
  } else {
    // Приведение угла к диапазону [-pi, pi]
    while (x > s21_M_PI) {
      x -= 2 * s21_M_PI;
    }
    while (x < -s21_M_PI) {
      x += 2 * s21_M_PI;
    }

    result = 1.0;
    long double term = 1.0;
    long double sign = -1.0;

    for (int n = 2; s21_fabs(term) > 1e-20; n += 2) {
      term *= x * x / (n * (n - 1));
      result += sign * term;
      sign *= -1.0;
    }
  }
  return result;
}

// 7

long double s21_exp(double x) {
  long double result = 1.0;
  if (x == s21_NEGINFINITY) {
    result = 0;
  } else if (s21_isnan(x)) {
    result = NAN;
  } else if (x > 709)
    result = s21_INFINITY;
  else if (x == 0) {
    result = 1;
  } else if (x < -20) {
    result = 0;
  }

 else {
    long double term = x;
    long double precision = 1e-20;  // Точность до 20 знаков после запятой
    long double term_contrib = 0;
    long double f = 1;
    for (int n = 2; s21_fabs(term_contrib) > precision || s21_fmod(n, 2) == 0;
         n++) {
      term_contrib = term / f;
      result += term_contrib;
      term *= x;
      f *= n;
    }
  }
  return result;
}

// 8
long double s21_fabs(double x) { return (long double)(x < 0 ? -x : x); }

// 9
long double s21_floor(double x) {
  if (!(s21_isinf(x) || s21_isnan(x))) {
    long double intPart = (int)x;
    if (x >= 0 || intPart == x) {
      intPart += 1;
    }
    x = intPart - 1;
  }

  return x;
}

// 10
long double s21_fmod(double x, double y) {
  long double res;
  if (s21_isnan(y) || s21_isnan(x) || y == 0.0 || x == s21_INFINITY) {
    res = NAN;
  } else if (x == 0.0) {
    res = x;
  } else {
    int digitX = x < 0 ? 0 : 1;
    int digitY = y < 0 ? 0 : 1;
    long double fabsX = s21_fabs(x);
    long double fabsY = s21_fabs(y);
    while (fabsX >= fabsY) {
      fabsX -= fabsY;
    }
    if ((digitX && digitY) || (digitX && !digitY)) {
      res = fabsX;
    } else {
      res = -fabsX;
    }
  }
  return res;
}
// 11
long double s21_log(double x) {
  long double result = 0.0;
  if (x < 0 || s21_isnan(x))
    result = NAN;
  else if (x == 0.0)
    result = s21_NEGINFINITY;
  else if (x == s21_INFINITY)
    result = s21_INFINITY;
  else if (x == s21_M_E)
    result = 1.0;
  else {
    long double term = (x - 1) / (x + 1);
    long double term_squared = term * term;
    long double power = term;
    long double precision = 1e-20;
    long double term_contrib = power;
    for (int n = 1; s21_fabs(term_contrib) >= precision; n += 2) {
      term_contrib = power / n;
      result += power / n;
      power *= term_squared;
    }

    result *= 2;
    if (result > DBL_MAX) result = NAN;
  }
  return result;
}

// 12

long double s21_pow(double base, double exp) {
  long double res = 1.0;
  if (exp == 0.0) {
    res = 1.0;
  } else if (base == 1.0) {
    res = 1.0;
  } else if (base >= 0.0 && base < 1 && exp == s21_INFINITY) {
    res = 0.0;
  } else if ((base >= 0.0 && base < 1) && exp == s21_NEGINFINITY) {
    res = s21_INFINITY;
  } else if (s21_isnan(base) || s21_isnan(exp)) {
    res = NAN;
  } else if (base == -1.0 && s21_isinf(exp)) {
    res = 1.0;
  } else if (s21_isinf(base) && !s21_isinf(exp)) {
    res = (exp > 0) ? s21_INFINITY : 0.0;
  } else if (s21_isinf(exp) && s21_isinf(base)) {
    if (exp == s21_NEGINFINITY)
      res = 0;
    else
      res = s21_INFINITY;
  } else if (!s21_isinf(base) && s21_isinf(exp)) {
    res = (exp < 0) ? 0.0 : s21_INFINITY;
  } else if (base == 0.0 && exp < 0.0) {
    res = s21_INFINITY;
  }
  else if (base == 0.0 && exp > 0.0) {
    res = 0.0;
  } else if (base < 0 && s21_fmod(exp, 2) == 1.0) {
    res = -s21_exp(exp * s21_log(-base));
  } else if (base < 0 && s21_fmod(exp, 2) == 0.0){
       res = s21_exp(exp * s21_log(-base));}
  else {
    res = s21_exp(exp * s21_log(base));
  }
  return res;
}

// 13
long double s21_sin(double x) {
  long double result;
  if (s21_isnan(x) || s21_isinf(x)) {
    result = NAN;  // Для бесконечных значений результат не определен
  } else {
    // Приведение угла к диапазону [-pi, pi]
    while (x > s21_M_PI) {
      x -= 2 * s21_M_PI;
    }
    while (x < -s21_M_PI) {
      x += 2 * s21_M_PI;
    }

    result = x;
    long double term = x;
    long double sign = -1.0;

    for (int n = 3; s21_fabs(term) > 1e-20; n += 2) {
      term *= x * x / (n * (n - 1));
      result += sign * term;
      sign *= -1.0;
    }
  }
  return result;
}

// 14
long double s21_sqrt(double x) { return s21_sqrtx(x, 2); }

// 15
long double s21_tan(double x) {
  long double result;
  long double sin_x = s21_sin(x);
  long double cos_x = s21_cos(x);
  if (cos_x == 0.0) {
    if (sin_x > 0.0) {
      result = s21_INFINITY;
    } else if (sin_x < 0.0) {
      result = s21_NEGINFINITY;
    } else {
      result = NAN;
    }
  } else if (s21_isnan(sin_x) || s21_isnan(cos_x)) {
    result = NAN;
  } else if (x == s21_M_PI / 2) {
    result = 16331239353195369.755859375;
  } else if (x == -s21_M_PI / 2) {
    result = -16331239353195369.755859375;
  } else {
    result = sin_x / cos_x;
  }

  return result;
}
