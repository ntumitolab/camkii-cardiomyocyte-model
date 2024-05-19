#===
Constants
===#

# Units
const second = 1           # second
const minute = 60second    # minute
const ms = 1e-3second      # millisecond
const Hz = inv(second)     # Herz
const kHz = 1e3Hz          # kilohertz
const metre = 1            # meter
const cm = 0.01metre       # centimeter
const cm² = cm^2           # square centimeter
const μm = 1E-6metre       # Micrometer
const mL = cm^3            # milliliter = cubic centimeter
const Liter = 1e3mL        # liter
const μL = 1E-6Liter
const pL = 1E-12Liter
const mM = 1
const Molar = 1000mM       # molar (1000 since the SI units is mM)
const μM = 1E-3mM          # micromolar
const nM = 1E-6mM          # nanomolar
const Amp = 1              # ampere
const mA = 1E-3Amp         # milliampere
const μA = 1E-6Amp         # micropampere
const Volt = 1             # volt
const mV = 1E-3Volt        # millivolt
const mS = mA / Volt       # milliseimens
const Farad = Amp * second / Volt
const μF = 1E-6Farad
const T₀ = 310.0           # Default temp (37C)
const Faraday = 96485.0    # Faraday constant (columb / mol)
const RGAS = 8.314            # Ideal gas constant
const VT = RGAS * T₀ / Faraday      # Thermal voltage (@37C), around 26.7 mV
const iVT = inv(VT)        # Reciprocal of thermal voltage (@37C)
# Utility functions

"""
Regular Hill/MM function
"""
hil(x, k = one(x)) = x / (x + k)
hil(x, k, n) = hil(NaNMath.pow(x, n), NaNMath.pow(k, n))

"""
Repressive Hill/MM function
"""
hilr(x, k = one(x)) = hil(k, x)
hilr(x, k, n) = hil(k, x, n)

"""
Logistic sigmoid function.
"""
expit(x) = hilr(exp(-x))

"""
    exprel(x, em1 = expm1(x))
"""
exprel(x, em1 = expm1(x)) = x / em1

"""Nernst potential"""
nernst(x_out, x_in) = VT * NaNMath.log(x_out / x_in)
nernst(x_out, x_in, z) = nernst(x_out, x_in) / z

"GHK flux equation"
ghk(px, x_i, x_o, zvfrt, ezvfrtm1 = expm1(zvfrt), z = 1) = px * z * F * ((ezvfrtm1 + 1) * x_i - x_o) * exprel(zvfrt, ezvfrtm1)

"GHK flux equation from voltage across the membrane"
function ghkVm(px, vm, x_i, x_o, z = 1)
    zvfrt = z * vm * iVT
    em1 = expm1(zvfrt)
    return ghk(px, x_i, x_o, zvfrt, em1, z)
end

"Accumulate reaction rates in a reaction network"
function add_rate!(rates, v, substrates, products)
    for s in substrates
        rates[s] -= v
    end
    for p in products
        rates[p] += v
    end
    return rates
end
