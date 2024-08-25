# Electrophysiology of neonatal rat CMC model
# Model of ECC of rat neonatal ventricular myocyte 2009
# Code & model: Topi Korhonen, University of Oulu (topi.korhonen@oulu.fi)
#
# PLEASE MENTION THE FOLLOWING REFERENCE WHEN USING THIS CODE OR PART OF IT:
# Korhonen et al. "Model of excitation-contraction coupling of rat neonatal
# ventricular myocytes" Biophys J. 2009, Feb; 96(3):1189-1209

function get_nak_sys(vm, Nai, Nao, Ko; name=:naksys)
    @parameters begin
        INaKmax = 2.7μA / cm^2
        KmNaiNaK = 18600μM
        nNaK = 3.2
        KmKoNaK = 1500μM
    end

    @variables INaK(t)
    sigma = 1 / 7 * expm1(Nao / 67300μM)
    fNaK = inv(1 + 0.1245 * exp(-0.1vm * iVT) + 0.0365sigma * exp(-vm * iVT))
    fKo = hil(Ko, KmKoNaK)
    fNai = hil(Nai, KmNaiNaK, nNaK)
    eqs = INaK ~ INaKmax * fNaK * fKo * fNai
    return ODESystem(eqs, t; name)
end

function build_neonatal_ecc_sys(;
    simplify=true,
    rSR_true=6μm,
    rSL_true=10.5μm,
    name=:neonataleccsys
)
    @parameters begin
        Ca_o = 1796μM
        Na_o = 154578μM
        K_o = 5366μM
        Mg_i = 1000μM
        ROS = 0μM
        ISO = 0μM
        ATP = 5000μM
        Cm = 1μF / cm^2
        Acap = 4π * rSL_true^2
        VSR = 0.043 * 1.5 * 1.4pL
        VNSR = 0.9 * VSR
        VJSR = VSR - VNSR
        Vmyo = 4 // 3 * π * (rSL_true^3 - rSR_true^3)
        ACAP_VMYO_F = Acap * Cm / Faraday / Vmyo
        Istim = 0
    end

    @variables begin
        Na_i(t) = 13838.37602μM
        K_i(t) = 150952.75035μM
        CaNSR(t) = 619.09843μM
        CaJSR(t) = 613.87556μM
        vm(t) = -68.79268mV
        E_Na(t)
        E_K(t)
        E_Ca(t)
    end

    barsys = get_bar_sys(; ATP, ISO)
    @unpack LCCa_PKAp, LCCb_PKAp, fracPLBp, TnI_PKAp, IKUR_PKAp = barsys
    capdesys = get_ca_pde_sys(; TnI_PKAp, rSR_true, rSL_true)
    @unpack Cai_sub_SL, Cai_sub_SR, Cai_mean, JCa_SL, JCa_SR = capdesys
    camkiisys = get_camkii_sys(; ROS, Ca=Cai_mean)
    ICa_scale = get_ICa_scalep(LCCb_PKAp)
    ncxsys = get_ncx_sys(Na_i, Cai_sub_SL, Na_o, Ca_o, ICa_scale)
    @unpack INaCa = ncxsys
    lccsys = get_lcc_sys(Cai_sub_SL, Ca_o, vm, ICa_scale)
    @unpack ICaL = lccsys
    tccsys = get_tcc_sys(vm, E_Ca)
    @unpack ICaT = tccsys
    ibgsys = get_ibg_sys(vm, E_Na, E_Ca)
    @unpack ICab, INab = ibgsys
    ifsys = get_if_sys(vm, E_Na, E_K)
    @unpack IfNa, IfK, If = ifsys
    inasys = get_ina_sys(vm, E_Na)
    @unpack INa = inasys
    iksys = get_ik_eqs(vm, E_K, K_i, K_o, Na_i, Na_o, IKUR_PKAp)
    @unpack IK1, Ito, IKs, IKr = iksys
    ryrsys = get_ryr_sys(Cai_sub_SR, CaJSR)
    @unpack Jrel = ryrsys
    sercasys = get_serca_sys(Cai_sub_SR, CaNSR, CaJSR, fracPLBp)
    @unpack Jup, Jleak, Jtr, betaSR = sercasys
    naksys = get_nak_sys(vm, Na_i, Na_o, K_o)
    @unpack INaK = naksys

    eqs = [
        E_Na ~ nernst(Na_o, Na_i),
        E_K ~ nernst(K_o, K_i),
        E_Ca ~ nernst(Ca_o, Cai_sub_SL, 2),
        JCa_SL ~ (2 * INaCa - ICaL - ICaT - ICab) * ACAP_VMYO_F * Vmyo,
        JCa_SR ~ Jleak - Jup + Jrel,
        D(CaJSR) * VJSR ~ betaSR * (-Jrel + Jtr),
        D(CaNSR) * VNSR ~ Jup - Jleak - Jtr,
        D(vm) * Cm ~ INab + INaCa + ICaL + ICaT + If + Ito + IK1 + IKs + IKr + INa + INaK + ICab + Istim,
        D(Na_i) ~ -(IfNa + INab + INa + 3 * INaCa + 3 * INaK) * ACAP_VMYO_F,
        D(K_i) ~ -(IfK + Ito + IK1 + IKs + IKr + Istim - 2 * INaK) * ACAP_VMYO_F
    ]
    sys = ODESystem(eqs, t; name)

    for s2 in (barsys, capdesys, camkiisys, ncxsys, lccsys, tccsys, ibgsys, ifsys, inasys, iksys, ryrsys, sercasys, naksys)
        sys = extend(s2, sys; name)
    end

    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end
