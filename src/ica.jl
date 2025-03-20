"Plama membrane calcium currents"
function get_ica_sys(nai, cai, nao, cao, vm; LCCb_PKAp=0, name=:icasys)
    @parameters begin
        ICa_scale0 = 0.95 # or 5.25
        fracLCCbp0 = 0.250657   # Derived quantity - (LCCbp(baseline)/LCCbtot)
        fracLCCbpISO = 0.525870 # Derived quantity - (LCCbp(ISO)/LCCbtot)
        fNaCa = 1
        kNaCa = 2.268e-4μAμF⁻¹ / mM^4
        dNaCa = 1e-4 / mM^4
        gNaCa = 0.5
        GCaL = 6.3e-5 * (metre^3 / second / Farad)
        taufca = 10ms
        gCaT = 0.2mSμF⁻¹
        gCab = 0.0008mSμF⁻¹
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

    # Calcium flux scaled by LCC (phosphorylated beta subunit)
    a_favail = (1.56 - 1) / (fracLCCbpISO / fracLCCbp0 - 1)         # fracLCCbp ISO (x1.56 0.1um ISO)
    favail = (1 - a_favail) + a_favail * (LCCb_PKAp / fracLCCbp0)   # Test (max x2.52 100# phosph)

    # Na-Ca exchanger (NCX)
    a = nai^3 * cao
    b = nao^3 * cai * fNaCa
    inaca = ICa_scale * kNaCa * (_bf(gNaCa * vm) * a - _bf((gNaCa - 1) * vm) * b) / (1 + dNaCa * (a + b))
    # L-type calcium channel (LCC)
    V = vm  # Now we use mV as the base unit
    alphad = 1.4 * expit((V + 35mV) * inv(13mV)) + 0.25
    betad = 1.4 * expit((V + 5mV) * inv(-5mV))
    gammad = expit((V - 50mV) * inv(20mV))
    alphafca = hilr(cai, 0.4875μM, 8)
    betafca = 0.1 * expit(-(cai - 0.5μM) * inv(0.1μM))
    gammafca = 0.2 * expit(-(cai - 0.75μM) * inv(0.8μM))
    kfca = 1 - (fcainf > i_fca) * (V > -60mV)

    return ODESystem([
            ICa_scale ~ ICa_scale0 * favail,
            E_Ca ~ nernst(cao, cai, 2),
            INaCa ~ inaca,
            ICaL ~ ICa_scale * i_d * i_f * i_fca * ghkVm(GCaL, vm, cai, 0.341 * cao, 2),
            dinf ~ expit((V + 11.1mV) * inv(7.2mV)),
            taud ~ (alphad * betad + gammad) * 1ms,
            finf ~ expit((V + 23.3mV) * inv(-5.4mV)),
            tauf ~ (1125 * exp(-(V + 27mV)^2 * inv(240mV^2)) + 165 * expit((V - 25mV) * inv(10mV)) + 120) * ms,
            fcainf ~ (alphafca + betafca + gammafca + 0.23) / 1.46,
            D(i_d) ~ (dinf - i_d) / taud,
            D(i_f) ~ (finf - i_f) / tauf,
            D(i_fca) ~ kfca / taufca * (fcainf - i_fca),
            binf ~ expit((V + 37.49098mV) * inv(5.40634mV)),
            taub ~ (0.6 + 5.4 * expit(-(V + 100mV) * 0.03/mV)) * ms,
            ginf ~ expit(-(V + 66mV) * inv(6mV)),
            taug ~ (1 + 40 * expit(-(V + 65mV) * inv(12.5mV))) * ms,
            D(i_b) ~ (binf - i_b) / taub,
            D(i_g) ~ (ginf - i_g) / taug,
            ICaT ~ gCaT * i_b * i_g * (vm - E_Ca + 106.5mV),
            ICab ~ gCab * (vm - E_Ca),
        ], t; name)
end
