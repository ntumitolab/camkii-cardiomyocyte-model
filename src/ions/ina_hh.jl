using ModelingToolkit

# Sodium channel (HH model) dynamics
_am(x) = 0.32 / 0.1 * exprel(-0.1 * (x + 47.13))
_bm(x) = 0.08 * exp(-x / 11)
_ahhi(x) = 0
_ahlo(x) = 0.135 * exp((x + 80) / (-6.8))
_bhhi(x) = 0.66 / 0.13 * expit((x + 10.66) / 11.1)
_bhlo(x) = 1.1 * 2.56 * exp(0.079 * (x - 2)) + 3.1e5 * exp(0.35 * (x - 2))
_ajhi(x) = 0
_ajlo(x) = (-1.2714e5 * exp(0.2444 * x) - 3.474e-5 * exp(-0.04391 * x)) * (x + 37.78) * expit(-0.311 * (x + 79.23))
_bjhi(x) = 0.3 * exp(-2.535e-7 * x) * expit(0.1 * (x + 32))
_bjlo(x) = 0.1212 * exp(-0.01052 * x) * expit(0.1378 * (x + 40.14))
_ah(x) = ifelse(x >= -40, _ahhi(x), _ahlo(x))
_bh(x) = ifelse(x >= -40, _bhhi(x), _bhlo(x))
_aj(x) = ifelse(x >= -40, _ajhi(x), _ajlo(x))
_bj(x) = ifelse(x >= -40, _bjhi(x), _bjlo(x))
_hlss(x) = expit(-(x + 91) / 6.1)

function ina_hh_eqs(inashift = 0)
    @parameters tauhl = 600 #ms
    @variables t vm(t) na_m(t) na_h(t) na_j(t) nal_h(t)
    D = Differential(t)
    v = vm(t) - inashift
    eqs = [
        D(na_m) ~ _am(v) * (1 - na_m) + _bm(v)
        D(na_h) ~ _ah(v) * (1 - na_h) + _bh(v)
        D(na_j) ~ _aj(v) * (1 - na_j) + _bj(v)
        D(nal_h) ~ (_hlss(vm) - nal_h) / tauhl
    ]
    return eqs
end

_i_na(gna, m, h, j, vm, ena) = gna * m^3 * h * j * (vm - ena)
_i_nal(gna, m, h, vm, ena) = gna * m^3 * h * (vm - ena)
