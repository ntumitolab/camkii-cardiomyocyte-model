# Plama membrane calcium currents
function get_ica_sys(nai, cai, nao, cao, vm; Acap=4π * (10μm)^2, Cm=1μF // cm^2, LCCb_PKAp=0, name=:icasys)
    @parameters begin
        ICa_scale0 = 0.95 # or 5.25
        fracLCCbp0 = 0.250657 # Derived quantity - (LCCbp(baseline)/LCCbtot)
        fracLCCbpISO = 0.525870 # Derived quantity - (LCCbp(ISO)/LCCbtot)
        fNaCa = 1
        kNaCa = 2.268e-016μAμF / μM^4
        dNaCa = 1e-16 / μM^4
        gNaCa = 0.5
        GCaL = 6.3e-5 * ((metre//10)^3 // ms // Farad)
        taufca = 10ms
        gCaT = 0.2mSμF
        gCab = 0.0008mSμF
    end

    @variables begin
        ICa_scale(t)
        JCa_SL(t)
        E_Ca(t)
        INaCa(t)
        ICaL(t)
        i_d(t) = 0.00033
        i_f(t) = 0.99869
        i_fca(t) = 0.9911
        dinf(t)
        taud(t)
        finf(t)
        tauf(t)
        fcainf(t)
        ICaT(t)
        i_b(t) = 0.00305
        i_g(t) = 0.61179
        binf(t)
        taub(t)
        ginf(t)
        taug(t)
        ICab(t)
    end

    # Calcium flux scaled by LCC (phosphorylated@beta subunit)
    a_favail = (1.56 - 1) / (fracLCCbpISO / fracLCCbp0 - 1)         # fracLCCbp ISO (x1.56 o.1 ISO)
    favail = (1 - a_favail) + a_favail * (LCCb_PKAp / fracLCCbp0)   # Test (max x2.52 100# phosph)

    # Na-Ca exchanger (NCX)
    a = nai^3 * cao
    b = nao^3 * cai * fNaCa
    # L-type calcium channel (LCC)
    V = vm * Volt / mV # Convert voltage to mV
    alphad = 1.4 * expit((V + 35) / 13) + 0.25
    betad = 1.4 * expit(-(V + 5) / 5)
    gammad = expit((V - 50) / 20)
    alphafca = hilr(cai, 0.000325mM * 1.5, 8)
    betafca = 0.1 * expit(-(cai - 0.0005mM) / 0.0001mM)
    gammafca = 0.2 * expit(-(cai - 0.00075mM) / 0.0008mM)
    kfca = 1 - (fcainf > i_fca) * (V > -60)

    return ODESystem([
            ICa_scale ~ ICa_scale0 * favail,
            E_Ca ~ nernst(cao, cai, 2),
            JCa_SL ~ (2 * INaCa - ICaL - ICaT - ICab) * (Acap * Cm // Faraday),
            INaCa ~ ICa_scale * kNaCa * (exp(iVT * gNaCa * vm) * a - exp(iVT * (gNaCa - 1) * vm) * b) / (1 + dNaCa * (a + b)),
            ICaL ~ ICa_scale * i_d * i_f * i_fca * ghkVm(GCaL, vm, cai, 0.341 * cao, 2),
            dinf ~ expit((V + 11.1) / 7.2),
            taud ~ (alphad * betad + gammad) * ms,
            finf ~ expit(-(V + 23.3) / 5.4),
            tauf ~ (1125 * exp(-(V + 27)^2 / 240) + 165 * expit((V - 25) / 10) + 120) * ms,
            fcainf ~ (alphafca + betafca + gammafca + 0.23) / 1.46,
            D(i_d) ~ (dinf - i_d) / taud,
            D(i_f) ~ (finf - i_f) / tauf,
            D(i_fca) ~ kfca / taufca * (fcainf - i_fca),
            binf ~ expit((V + 37.49098) / 5.40634),
            taub ~ (0.6 + 5.4 * expit(-(V + 100) * 0.03)) * ms,
            ginf ~ expit(-(V + 66) / 6),
            taug ~ (1 + 40 * expit(-(V + 65) * 0.08)) * ms,
            D(i_b) ~ (binf - i_b) / taub,
            D(i_g) ~ (ginf - i_g) / taug,
            ICaT ~ gCaT * i_b * i_g * (vm - E_Ca + 106.5mV),
            ICab ~ gCab * (vm - E_Ca),
        ], t; name)
end
