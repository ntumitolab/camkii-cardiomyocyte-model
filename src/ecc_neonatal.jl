# Electrophysiology of neonatal rat CMC model
# Model of ECC of rat neonatal ventricular myocyte 2009
# Code & model: Topi Korhonen, University of Oulu (topi.korhonen@oulu.fi)
#
# PLEASE MENTION THE FOLLOWING REFERENCE WHEN USING THIS CODE OR PART OF IT: Korhonen et al. "Model of excitation-contraction coupling of rat neonatal ventricular myocytes" Biophys J. 2009, Feb; 96(3):1189-1209
function get_nak_sys(na_i, na_o, k_o, vm; name=:naksys)
    @parameters begin
        INaKmax = 2.7μAμF
        KmNaiNaK = 18.6mM
        nNaK = 3.2
        KmKoNaK = 1.5mM
    end

    @variables INaK(t)
    sigma = 1 / 7 * expm1(na_o / 67.3mM)
    fNaK = inv(1 + 0.1245 * exp(-0.1vm * iVT) + 0.0365 * sigma * exp(-vm * iVT))
    fKo = hil(k_o, KmKoNaK)
    fNai = hil(na_i, KmNaiNaK, nNaK)
    return ODESystem([INaK ~ INaKmax * fNaK * fKo * fNai], t; name)
end

function build_neonatal_ecc_sys(;
    rSR_true=6μm,
    rSL_true=10.5μm,
    dx=0.1μm,
    name=:neonataleccsys,
    simplify=true,
    reduce_iso=false,
)
    @parameters begin
        ca_o = 1796μM
        na_o = 154578μM
        k_o = 5366μM
        mg_i = 1000μM
        ROS = 0μM
        ISO = 0μM
        ATP = 5mM
        Istim = 0
        # cell geometry
        Cm = 1μF / cm^2
        Acap = 4π * rSL_true^2
        Vmyo = 4 / 3 * π * (rSL_true^3 - rSR_true^3) # 3.944 pL
        ACAP_F = Acap * Cm / Faraday
        V_sub_SR = 4 / 3 * pi * ((rSR_true + dx)^3 - (rSR_true)^3) # 0.046 pL
        V_sub_SL = 4 / 3 * pi * (rSL_true^3 - (rSL_true - dx)^3)   # 0.137 pL
    end

    @variables begin
        na_i(t) = 2.243mM
        k_i(t) = 152.586mM
        vm(t) = -67.3mV
        JCa_SL(t)
        JCa_SR(t)
    end

    barsys = if reduce_iso
        get_bar_sys_reduced(ISO)
    else
        get_bar_sys(ATP, ISO)
    end
    @unpack LCCa_PKAp, LCCb_PKAp, fracPLBp, TnI_PKAp, IKUR_PKAp = barsys
    capdesys = get_ca_pde_sys(; JCa_SR, JCa_SL, TnI_PKAp, rSR_true, rSL_true, dx, V_sub_SL)
    @unpack Cai_sub_SL, Cai_sub_SR, Cai_mean = capdesys
    camkiisys = get_camkii_sys(Cai_mean; ROS)
    icasys = get_ica_sys(na_i, Cai_sub_SL, na_o, ca_o, vm; Acap, Cm, LCCb_PKAp)
    @unpack INaCa, ICaL, ICaT, ICab = icasys
    inasys = get_ina_sys(na_i, na_o, vm)
    @unpack INa, INab = inasys
    iksys = get_ik_sys(k_i, k_o, na_i, na_o, vm; IKUR_PKAp)
    @unpack IK1, Ito, IKs, IKr, IfNa, IfK, If = iksys
    sersys = get_ser_sys(Cai_sub_SR; fracPLBp, V_sub_SR)
    naksys = get_nak_sys(na_i, na_o, k_o, vm)
    @unpack INaK = naksys

    eqs = [
        D(vm) ~ -(INab + INaCa + ICaL + ICaT + If + Ito + IK1 + IKs + IKr + INa + INaK + ICab + Istim), # Current  normalized by capacitance
        D(na_i) ~ -(IfNa + INab + INa + 3 * INaCa + 3 * INaK) * ACAP_F / Vmyo,
        D(k_i) ~ -(IfK + Ito + IK1 + IKs + IKr + Istim - 2 * INaK) * ACAP_F / Vmyo,
    ]

    sys = ODESystem(eqs, t; name)
    for s2 in (barsys, capdesys, camkiisys, icasys, inasys, iksys, sersys, naksys)
        sys = extend(sys, s2; name)
    end
    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end
