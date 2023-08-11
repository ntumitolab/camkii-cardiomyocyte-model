# Constants
const R = 8314      # [J/kmol*K]
const Frdy = 96485  # [C/mol]
const Temp = 310    # [K] 310 K (37 C) for BT / 295 K (22 C) for RT
const FoRT = Frdy / R / Temp
const VT = inv(FoRT) # Thermal temparature (mV)
const Qpow = (Temp - 310) / 10 # Temp factor

# Utility functions

"""Hill function"""
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

"""
Relative exponential function.
Returns one when x is close to zero
"""
exprel(x) = ifelse(isapprox(x, 0), one(x / expm1(x)), x / expm1(x))

"""Logistic function"""
expit(x) = 1 / (1 + exp(-x))

"""Nernst potential"""
nernst(xo, xi) = VT * log(xo / xi)
nernst(xo, xi, z) = nernst(xo, xi) / z
