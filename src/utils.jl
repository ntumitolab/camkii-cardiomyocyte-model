# Units
const second = 1                # SI time unit
const mol = 1                   # SI quantity unit
const Volt = 1                  # SI electric potential unit
const Kelvin = 1                # SI temperature unit
const metre = 1                 # SI length unit
const ms = second//1000         # millisecond
const m³ = metre^3              # volume unit
const mM = mol//m³              # concentration unit millimolar (is actually SI)
const mV = Volt // 1000         # millivolt
const mmol = mol//1000          # millimole
const minute = 60second         # minute
const Hz = 1//second            # Herz (once per second)
const kHz = 1000Hz              # kilohertz
const MHz = 1000kHz             # megahertz
const cm = metre //100          # centimeter
const cm² = cm^2                # square centimeter
const μm = metre / 10^6         # micrometer
const mL = cm^3                 # milliliter = cubic centimeter
const litre = 1000mL            # liter
const μL = μm^3                 # microliter
const pL = 1e-12litre           # picoliter
const Molar = 1000mM            # molarity is used in equilibrium constants
const μM = mM//1000             # micromolar
const nM = μM//1000             # nanomolar
const Ampere = 1                # electric current SI unit
const Columb = Ampere * second  # electric charge unit
const Farad = Columb // Volt    # capacitance unit
const μF = Farad // 10^6        # microfarad
const μFcm⁻² = μF // cm²        # area capacitance (μF/cm²)
const μAμF⁻¹ = Ampere // Farad          # common current density
const μAcm⁻² = Ampere // 10^6 // cm²    # current area density
const Joule = Columb * Volt          # energy SI unit
const Seimens = Ampere // Volt       # conductance unit
const milliseimens = Seimens // 1000 # milliseimens
const mScm⁻² = milliseimens // cm²      # conductance area density
const mSμF⁻¹ = milliseimens // μF       # conductance density
const Faraday = 96485Columb // mol      # Faraday constant (columb / mol)
const T₀ = 310Kelvin                    # Default temp (37C)
const RGAS = 8.314Joule / Kelvin / mol  # Ideal gas constant (J/K⋅mol)
const VT = RGAS * T₀ / Faraday  # Thermal voltage (@37C), 26.7 mV
const iVT = inv(VT)             # Reciprocal of thermal voltage (0.037 per mV)

"""
Boltzmann factor for voltage
"""
_bf(vm) = exp(iVT * vm)

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

    expit(x[, a=1, b=1]) = a / (b + exp(-x))
"""
expit(x, a=one(x), b=one(x)) = a / (b + exp(-x))

"""
Relative exponential function

    exprel(x, em1 = expm1(x))
"""
exprel(x, em1=expm1(x)) = x / em1

"""Nernst potential"""
nernst(x_out, x_in) = VT * NaNMath.log(x_out / x_in)
nernst(x_out, x_in, z) = nernst(x_out, x_in) / z

"""
GHK flux equation

    ghk(px, vm, x_i, x_o, z = 1)

https://en.wikipedia.org/wiki/Goldman%E2%80%93Hodgkin%E2%80%93Katz_flux_equation

- px: the permeability of the membrane for ion x measured in meter per second
- vm: the transmembrane potential in volts
- x_i: the intracellular concentration of ion (mM)
- x_o: the extracellular concentration of ion (mM)
- z: the valence of ion x
"""
function ghk(px, vm, x_i, x_o, z=1)
    zvfrt = z * vm * iVT
    return px * z * Faraday * exprel(zvfrt) * (exp(zvfrt) * x_i - x_o)
end
const ghkVm = ghk

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
