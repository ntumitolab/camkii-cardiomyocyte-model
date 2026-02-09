#===
Electrophysiology of neonatal rat CMC model
# Model of ECC of rat neonatal ventricular myocyte 2009
# Code & model: Topi Korhonen, University of Oulu (topi.korhonen@oulu.fi)
# PLEASE MENTION THE FOLLOWING REFERENCE WHEN USING THIS CODE OR PART OF IT: Korhonen et al. "Model of excitation-contraction coupling of rat neonatal ventricular myocytes" Biophys J. 2009, Feb; 96(3):1189-1209
# https://pmc.ncbi.nlm.nih.gov/articles/PMC2716686/
===#
function get_nak_eqs(na_i, na_o, k_o, vm)
    @parameters begin
        INaKmax = 2.7μAμF
        KmNaiNaK = 18.6mM
        nNaK = 3.2
        KmKoNaK = 1.5mM
    end
    @variables INaK(t)
    sigma = 1 / 7 * expm1(na_o / 67.3mM)
    fNaK = inv(1 + 0.1245 * exp(-0.1vm * iVT) + 0.0365 * sigma * exp(-vm * iVT))
    inak = INaKmax * fNaK * hil(k_o, KmKoNaK) * hil(na_i, KmNaiNaK, nNaK)
    eqs_inak = [INaK ~ inak]
    return (; eqs_inak, INaK)
end

function get_nak_sys(na_i, na_o, k_o, vm; name=:naksys)
    @unpack eqs_inak = get_nak_eqs(na_i, na_o, k_o, vm)
    return System(eqs_inak, t; name)
end

function build_neonatal_ecc_sys(;
    rSR_true=6μm,
    rSL_true=10.5μm,
    dx=0.1μm,
    name=:neonataleccsys,
    )
    @parameters begin
        Istim(t) = 0μAμF
        ca_o = 1.796mM
        na_o = 154.578mM
        k_o = 5.366mM
        mg_i = 1mM
        ROS = 0μM
        ISO = 0μM
        ATP = 5mM
        # cell geometry
        Cm = 1μF / cm^2
        Acap = 4π * rSL_true^2
        Vmyo = 4π / 3 * (rSL_true^3 - rSR_true^3) # 3.944 pL
        ACAP_F = Acap * Cm / Faraday
        V_sub_SR = 4π / 3 * ((rSR_true + dx)^3 - (rSR_true)^3) # 0.046 pL
        V_sub_SL = 4π / 3 * (rSL_true^3 - (rSL_true - dx)^3)   # 0.137 pL
    end

    @variables begin
        na_i(t) = 13.83837602mM
        k_i(t) = 150.95275035mM
        vm(t) = -68.79268mV
        JCa_SL(t)
        JCa_SR(t)
    end

    @unpack LCCa_PKAp, LCCb_PKAp, fracPLBp, TnI_PKAp, IKUR_PKAp, eqs_bar = get_bar_eqs_reduced(ISO)
    @unpack Cai_sub_SL, Cai_sub_SR, Cai_mean, eqs_cai = get_ca_pde_eqs(; JCa_SR, JCa_SL, TnI_PKAp, rSR_true, rSL_true, dx)
    @unpack eqs_camkii = get_camkii_simp_eqs(;Ca=Cai_mean, ROS)
    @unpack INaCa, ICaL, ICaT, ICab, eqs_ica = get_ica_eqs(na_i, Cai_sub_SL, na_o, ca_o, vm; LCCb_PKAp)
    @unpack eqs_ina, INa, INab, E_Na = get_ina_eqs(; na_i, na_o, vm)
    @unpack eqs_ik, IK1, Ito, IKs, IKr, IfNa, IfK, If = get_ik_eqs(; na_i, k_i, na_o, k_o, vm, IKUR_PKAp, E_Na)
    @unpack eqs_sr, Jrel, Jup, Jleak, Jtr, JCa_SR = get_ser_eqs(Cai_sub_SR; fracPLB_CKp=0, fracPLBp, RyR_CKp=0.2, V_sub_SR)
    @unpack eqs_inak, INaK = get_nak_eqs(na_i, na_o, k_o, vm)

    eqs = [
        D(vm) ~ -(INab + INaCa + ICaL + ICaT + If + Ito + IK1 + IKs + IKr + INa + INaK + ICab + Istim), ## Currents are normalized by capacitance
        D(na_i) ~ -(IfNa + INab + INa + 3 * INaCa + 3 * INaK) * ACAP_F / Vmyo,
        D(k_i) ~ -(IfK + Ito + IK1 + IKs + IKr + Istim - 2 * INaK) * ACAP_F / Vmyo,
        JCa_SL ~ (2 * INaCa - ICaL - ICaT - ICab) * ACAP_F / 2 / V_sub_SL,
    ]

    eqs = vcat(eqs, eqs_bar, eqs_cai, eqs_camkii, eqs_ica, eqs_ina, eqs_ik, eqs_sr, eqs_inak)
    return System(eqs, t; name)
end
