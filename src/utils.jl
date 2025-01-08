# Units
const ms = 1               # millisecond
const second = 1000ms      # second is the SI unit
const minute = 60second    # minute
const Hz = inv(second)     # herz
const kHz = inv(second)    # kilohertz
const metre = 1            # meter
const cm = 0.01metre       # centimeter
const cm² = cm^2           # square centimeter
const μm = metre / 10^6    # Micrometer
const mL = cm^3            # milliliter = cubic centimeter
const Liter = 1000mL       # liter
const μL = μm^3            #
const pL = Liter / 10^12   # picoliter
const mmol = 1
const mol = 1000mmol
const μM = mmol/metre^3    # micromolar
const mM = 1000μM          # mM is the SI unit
const Molar = 1000mM       # Molar is used in equilibrium constants
const nM = μM / 10^3       # nanomolar
const Ampere = 1           # current unit Ampere
const mA = Ampere / 10^3   # milliampere
const μA = Ampere / 10^6   # micropampere
const Joule = 10^6         # energy unit Joule
const Kelvin = 1           # temperature unit Kelvin
const Columb = Ampere * second # unit of electric charge
const Volt = Joule / Columb # voltage
const mV = Volt / 10^3     # millivolt
const milliseimens = Ampere / Volt / 10^3 # milliseimens
const Farad = Columb / Volt
const μF = Farad / 10^6
const T₀ = 310Kelvin           # Default temp (37C)
const Faraday = 96485Columb / mol # Faraday constant (columb / mol)
const RGAS = 8.314Joule/Kelvin/mol # Ideal gas constant (J/K⋅mol)
const VT = RGAS * T₀ / Faraday # Thermal voltage (@37C), about 26.7 mV
const iVT = inv(VT)            # Reciprocal of therm al voltage
const μAμF = μA / μF           # Common unit for current density, normalized by capacitance
const mSμF = milliseimens / μF # Common unit for conductance, normalized by capacitance

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
