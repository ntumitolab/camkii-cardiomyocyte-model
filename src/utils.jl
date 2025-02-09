# Units
const ms = 1                    # millisecond
const μM = 1                    # micromolar = 1
const mV = 1                    # millivolt
const Kelvin = 1                # temperature unit Kelvin
const mmol = 1                  # millimole = 1
const second = 1000ms           # second is the SI unit
const minute = 60second         # minute
const Hz = inv(second)          # Herz
const kHz = inv(ms)             # kilohertz
const m³ = mmol / μM            # cubic meter
const metre = cbrt(m³)          # meter
const cm = 0.01metre            # centimeter
const cm² = cm^2                # square centimeter
const μm = 1e-6metre            # micrometer
const mL = cm^3                 # milliliter = cubic centimeter
const Liter = 1000mL            # liter
const μL = μm^3                 # microliter
const pL = 1e-12Liter           # picoliter
const mol = 1000mmol            # mole
const mM = 1000μM               # mM is the SI unit
const Molar = 1000mM            # molarity is used in equilibrium constants
const nM = 0.001μM              # nanomolar
const μFcm⁻² = 1                # area capacitance (μF/cm²)
const μF = μFcm⁻² * cm²         # microfarad
const Farad = 1e6μF             # Farad
const μAμF = mV * inv(ms)       # common current density
const μA = μAμF * μF            # micropampere
const μAcm⁻² = μAμF * μFcm⁻²    # real current density
const Ampere = 1e6μA            # electric current unit Ampere
const Columb = Ampere * second  # electric charge unit Columb
const Volt = 1000mV             # electric potential unit Volt
const Joule = Columb * Volt     # energy unit Joule
const Seimens = Ampere / Volt   # conductance unit
const milliseimens = 0.001Seimens # milliseimens
const mScm⁻² = milliseimens / cm²
const mSμF = μAμF / mV          # conductance density
const Faraday = 96485Columb / mol # Faraday constant (columb / mol)
const T₀ = 310Kelvin            # Default temp (37C)
const RGAS = 8.314Joule / Kelvin / mol # Ideal gas constant (J/K⋅mol)
const VT = RGAS * T₀ / Faraday  # Thermal voltage (@37C), 26.7 mV
const iVT = inv(VT)             # Reciprocal of thermal voltage (0.037 per mV)

"""
Regular Hill/MM function
"""
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n::Int) = hil(x^n, k^n)
hil(x, k, n) = hil(NaNMath.pow(x, n), NaNMath.pow(k, n))

"""
Repressive Hill/MM function
"""
hilr(x, k=one(x)) = hil(k, x)
hilr(x, k, n) = hil(k, x, n)

"""
Logistic sigmoid function.
"""
expit(x, a=1, b=1) = a / (b + exp(-x))

"""
    exprel(x, em1 = expm1(x))
"""
exprel(x, em1=expm1(x)) = x / em1

"""Nernst potential"""
nernst(x_out, x_in) = VT * NaNMath.log(x_out / x_in)
nernst(x_out, x_in, z) = nernst(x_out, x_in) / z

"GHK flux equation"
ghk(px, x_i, x_o, zvfrt, ezvfrtm1=expm1(zvfrt), z=1) = px * z * Faraday * ((ezvfrtm1 + 1) * x_i - x_o) * exprel(zvfrt, ezvfrtm1)

"GHK flux equation from voltage across the membrane"
function ghkVm(px, vm, x_i, x_o, z=1)
    zvfrt = z * vm * iVT
    em1 = expm1(zvfrt)
    return ghk(px, x_i, x_o, zvfrt, em1, z)
end

"Accumulate chemical reaction rates into a look-up table"
function add_raw_rate!(lut, rate, substrates, products)
    for s in substrates
        lut[s] -= rate
    end
    for p in products
        lut[p] += rate
    end
    return lut
end

"Accumulate chemical reaction rates with law of mass action into a look-up table"
function add_rate!(lut, kf, substrates, kb, products)
    rate = prod(substrates; init=kf) - prod(products; init=kb)
    return add_raw_rate!(lut, rate, substrates, products)
end
