# Utility functions

"Hill function"
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)
