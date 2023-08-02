# Utility functions

"""Hill function"""
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

"""
Relative exponential function.
Returns one when x is close to zero
"""
exprel(x) = ifelse(isapprox(x, 0), one(x/expm1(x)), x/expm1(x))

"""Logistic function"""
expit(x) = 1 / (1 + exp(-x))
