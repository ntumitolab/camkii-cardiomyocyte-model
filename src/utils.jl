# Units
const second = 1           # second
const minute = 60second    # minute
const ms = second//1000    # millisecond
const Hz = 1 // second     # Herz
const kHz = 1000Hz         # kilohertz
const metre = 1            # meter
const cm = metre//100      # centimeter
const cm² = cm^2           # square centimeter
const μm = metre//10^6     # Micrometer
const mL = cm^3            # milliliter = cubic centimeter
const Liter = 1000mL       # liter
const μL = Liter//10^6
const pL = Liter//10^12
const mM = 1
const Molar = 1000mM       # molar (1000 since the SI units is mM)
const μM = mM//10^3        # micromolar
const nM = mM//10^6        # nanomolar
const Amp = 1              # ampere
const mA = Amp//10^3         # milliampere
const μA = Amp//10^6         # micropampere
const Volt = 1               # volt
const mV = Volt//10^3        # millivolt
const mS = mA // Volt        # milliseimens
const Farad = Amp * second // Volt
const μF = Farad//10^6
const T₀ = 310                  # Default temp (37C)
const Faraday = 96485           # Faraday constant (columb / mol)
const RGAS = 8314//1000         # Ideal gas constant (J/K⋅mol)
const VT = RGAS * T₀ // Faraday # Thermal voltage (@37C), about 26.7 mV
const iVT = inv(VT)             # Reciprocal of thermal voltage
const μAμF = μA // μF           # Common unit for current density, normalized by capacitance
const mSμF = ms // μF           # Common unit for conductance, normalized by capacitance

"""
Regular Hill/MM function
"""
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(NaNMath.pow(x, n), NaNMath.pow(k, n))

"""
Repressive Hill/MM function
"""
hilr(x, k=one(x)) = hil(k, x)
hilr(x, k, n) = hil(k, x, n)

"""
Logistic sigmoid function.
"""
expit(x) = inv(1 + exp(-x))

"""
    exprel(x, em1 = expm1(x))
"""
exprel(x, em1=expm1(x)) = x / em1

"""Nernst potential"""
nernst(x_out, x_in) = VT * NaNMath.log(x_out / x_in)
nernst(x_out, x_in, z) = nernst(x_out, x_in) // z

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
