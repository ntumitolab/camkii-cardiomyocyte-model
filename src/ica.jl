# Calcium currents
"Calcium flux scaled by phosphorylated LCC"
function get_ICa_scalep(LCCb_PKAp=0)
    @parameters begin
        ICa_scale = 0.95 # 5.25
        fracLCCbp0 = 0.250657 # Derived quantity - (LCCbp(baseline)/LCCbtot)
        fracLCCbpISO = 0.525870 # Derived quantity - (LCCbp(ISO)/LCCbtot)
    end
    a_favail = (1.56 - 1) / (fracLCCbpISO / fracLCCbp0 - 1) # fracLCCbp ISO (x1.56 o.1 ISO)
    favail = (1 - a_favail) + a_favail * (LCCb_PKAp / fracLCCbp0)  # Test (max x2.52 100# phosph)
    return ICa_scale * favail
end

"Na-Ca exchanger"
# TODO: too much, might be unit problem
function get_ncx_sys(nai, cai, nao, cao, vm, ICa_scale=1; name=:ncxsys)
    @parameters begin
        fNaCa = 1
        kNaCa = 2.2680e-016 * μA / cm^2 / μM^4
        dNaCa = 1e-16 / μM^4
        gamma = 0.5
    end
    @variables INaCa(t)
    eqs = INaCa ~ ICa_scale * kNaCa * ((exp(iVT * gamma * vm) * nai^3 * cao - exp(iVT * (gamma - 1) * vm) * cai * nao^3 * fNaCa) / (1 + dNaCa * (nao^3 * cai * fNaCa + nai^3 * cao)))
    return ODESystem(eqs, t; name)
end

"L-type calcium current (LCC)"
function get_lcc_sys(cai, cao, vm, ICa_scale=1; name=:lccsys)
    @parameters begin
        GCaL = 1.3125e-4 * 0.8 * 0.6 / ms
        taufca = 10ms
    end
    @variables begin
        ICaL(t)
        i_d(t) = 0.00033
        i_f(t) = 0.99869
        i_fca(t) = 0.9911
        dinf(t)
        taud(t)
        finf(t)
        tauf(t)
        fcainf(t)
    end

    V = vm * Volt / mV # Convert voltage to mV

    alphad = 1.4 * expit((V + 35) / 13) + 0.25
    betad = 1.4 * expit(-(V + 5) / 5)
    gammad = expit((V - 50) / 20)
    alphafca = hilr(cai, 0.000325mM * 1.5, 8)
    betafca = 0.1 * expit(-(cai - 0.0005mM) / 0.0001mM)
    gammafca = 0.2 * expit(-(cai - 0.00075mM) / 0.0008mM)
    kfca = 1 - (fcainf > i_fca) * (V > -60)

    eqs = [
        ICaL ~ ICa_scale * i_d * i_f * i_fca * ghkVm(GCaL, vm, cai, 0.341 * cao, 2),
        D(i_d) ~ (dinf - i_d) /taud,
        D(i_f) ~ (finf - i_f) /tauf,
        D(i_fca) ~ kfca/taufca * (fcainf - i_fca),
        dinf ~ expit((V + 11.1) / 7.2),
        taud ~ (alphad * betad + gammad) * ms,
        finf ~ expit(-(V + 23.3) / 5.4),
        tauf ~ (1125 * exp(-(V + 27)^2 / 240) + 165 * expit((V - 25) / 10) + 120) * ms,
        fcainf ~ (alphafca + betafca + gammafca + 0.23) / 1.46,
    ]
    return ODESystem(eqs, t; name)
end

"T-type calcium current (TCC)"
function get_tcc_sys(vm, E_Ca; name=:tccsys)
    @parameters gCaT = 0.2mS / cm^2
    @variables begin
        ICaT(t)
        i_b(t) = 0.00305
        i_g(t) = 0.61179
        binf(t)
        taub(t)
        ginf(t)
        taug(t)
    end

    V = vm * Volt / mV # Convert voltage to mV

    eqs = [
        binf ~ expit((V + 37.49098) / 5.40634),
        taub ~ (0.6 + 5.4 * expit(-(V + 100) * 0.03)) * ms,
        ginf ~ expit(-(V + 66) / 6),
        taug ~ (1 + 40 * expit(-(V + 65) * 0.08)) * ms,
        D(i_b) ~ (binf - i_b)/taub,
        D(i_g) ~ (ginf - i_g)/taug,
        ICaT ~ gCaT * i_b * i_g * (vm - E_Ca + 106.5mV),
    ]
    return ODESystem(eqs, t; name)
end
