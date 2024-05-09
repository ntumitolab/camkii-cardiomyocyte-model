# Electrophysiology of neonatal rat CMC model
using ModelingToolkit

"Calcium buffer term"
function beta_cai(Ca;
    TnI_PKAp=0,
    fracTnIpo=0.062698,
    TrpnTotal=35μM,
    CmdnTotal=50μM,
    KmCmdn=2.38μM,
    KmTrpn=0.5μM
)
    fPKA_TnI = (1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIpo)) # Max effect +61#
    KmTrpnNew = KmTrpn / fPKA_TnI
    return inv(1 + TrpnTotal * KmTrpnNew / (Ca + KmTrpnNew)^2 + CmdnTotal * KmCmdn / (Ca + KmCmdn)^2)
end

function build_neonatal_ecc_eqs(;
    LCCb_PKAp = 0,  # Fraction of LCC phosphorylated by PKA
)
    @parameters begin

        Ca_o = 1796μM
        Na_o = 154578μM
        K_o = 5366μM
        Mg_i = 1000μM
    end

    # PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph)
    ICa_scale = 0.95 # 5.25
    fracLCCbp0 = 0.250657 # Derived quantity - (LCCbp(baseline)/LCCbtot)
    fracLCCbpISO = 0.525870 # Derived quantity - (LCCbp(ISO)/LCCbtot)
    a_favail = (1.56 - 1) / (fracLCCbpISO / fracLCCbp0 - 1) # fracLCCbp ISO (x1.56 o.1 ISO)
    favail = (1 - a_favail) + a_favail * (LCCb_PKAp / fracLCCbp0)  # Test (max x2.52 100# phosph)
    ICa_scalep = ICa_scale * favail

    @variables begin
        t
        Na_i(t)
        K_i(t) = 135000μM
        Cai_sub_SL(t)
        E_Na(t)
        E_K(t)
        E_Ca(t)
        vm(t)
    end

    D = Differential(t)

    # NCX
    @parameters begin
        fNaCa = 1
        kNaCa = 2.2680e-016 * μA/cm^2 / μM^4
        dNaCa = 1e-16 / μM^4
        gamma = 0.5
    end
    @variables INaCa(t)

    # LCC
    @parameters GCaL = 1.3125e-4 * 0.8 * 0.6/ms
    @variables ICaL(t) i_d(t) i_f(t) i_fca(t)

    dinf = expit((vm + 11.1mV) / 7.2mV)
    alphad = 1.4 * expit((vm + 35mV) / 13mV) + 0.25
    betad = 1.4 * expit(-(vm + 5mV) / 5mV)
    gammad = expit((vm - 50mV) / 20mV)
    taud = (alphad * betad + gammad) * ms
    finf = expit(-(vm + 23.3mV) / 5.4mV)
    tauf = (1125 * exp(-(V + 27mV)^2 / 240mV^2) + 165 * expit((vm-25mV)/10mV) + 120) * ms



    # Nernst potentials
    eqs = [
        E_Na ~ nernst(Na_o, Na_i),
        E_K ~ nernst(K_o, K_i),
        E_Ca ~ nernst(Ca_o, Cai_sub_SL, 2),
        INaCa ~ ICa_scalep * kNaCa * ( (exp(iVT * gamma * vm) * Na_i^3 * Ca_o - exp(iVT * (gamma - 1) * vm) * Cai_sub_SL * Na_o^3 * fNaCa)/ (1 + dNaCa * (Na_o^3 * Cai_sub_SL * fNaCa + Na_i^3 * Ca_o))),
        ICaL ~ ICa_scalep * i_d * i_f * i_fca * ghkVm(GCaL, vm, Cai_sub_SL, 0.341 * Ca_o, 2)
    ]

    return eqs
end
