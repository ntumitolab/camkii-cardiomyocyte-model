# Electrophysiology of neonatal rat CMC model
# Model of ECC of rat neonatal ventricular myocyte 2009
# Code & model: Topi Korhonen, University of Oulu (topi.korhonen@oulu.fi)
#
# PLEASE MENTION THE FOLLOWING REFERENCE WHEN USING THIS CODE OR PART OF IT:
# Korhonen et al. "Model of excitation-contraction coupling of rat neonatal
# ventricular myocytes" Biophys J. 2009, Feb; 96(3):1189-1209
#
# ONLY FOR ACADEMIC USE, DO NOT DISTRIBUTE
#
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using NaNMath

"Calcium buffered by troponin and calmodulin"
function beta_cai(Ca;
    TnI_PKAp=0,
    TrpnTotal = 35μM,
    CmdnTotal = 50μM,
    KmCmdn = 2.38μM,
    KmTrpn = 0.5μM,
    fracTnIp0 = 0.062698 # Baseline effect
)
    fPKA_TnI = 1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIp0) # Max effect +61%
    KmTrpnNew = KmTrpn / fPKA_TnI
    return inv(1 + TrpnTotal * KmTrpnNew / (Ca + KmTrpnNew)^2 + CmdnTotal * KmCmdn / (Ca + KmCmdn)^2)
end

"Calcium diffusion between sarcolemma (SL) and sarcoplasmic reticulum (SR)"
function get_ca_pde_sys(;
    Cai_sub_SR_default=0.2556μM,
    Cai_sub_SL_default=0.25151μM,
    dx=0.1μm,
    rSR_true=6μm,
    rSL_true=10.5μm,
    TnI_PKAp=0,
    name=:capdesys
)
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial indices
    m = length(j)
    @variables Cai(t)[1:m] Cai_mean(t) Cai_sub_SR(t) Cai_sub_SL(t) JCa_SR(t) JCa_SL(t)
    @parameters begin
        Dca = 7μm^2 / ms
        V_sub_SR = 4 / 3 * pi * ((rSR_true + dx)^3 - (rSR_true)^3)
        V_sub_SL = 4 / 3 * pi * (rSL_true^3 - (rSL_true - dx)^3)
        TrpnTotal = 35μM
        CmdnTotal = 50μM
        KmCmdn = 2.38μM
        KmTrpn = 0.5μM
        fracTnIp0 = 0.062698 # Baseline effect
    end
    eqs = [
        Cai_mean ~ sum(collect(Cai)) / m,
        Cai_sub_SR ~ Cai[1],
        Cai_sub_SL ~ Cai[m],
        D(Cai[1]) ~ (Dca / (j[1] * dx^2) * ((1 + j[1]) * Cai[2] - 2 * j[1] * Cai[1] + (j[1] - 1) * Cai[1]) + JCa_SR / V_sub_SR) * beta_cai(Cai[1]; TnI_PKAp, TrpnTotal, CmdnTotal, KmCmdn, KmTrpn, fracTnIp0),
        D(Cai[m]) ~ (Dca / (j[m] * dx^2) * ((1 + j[m]) * Cai[m] - 2 * j[m] * Cai[m] + (j[m] - 1) * Cai[m-1]) + JCa_SL / V_sub_SL) * beta_cai(Cai[m]; TnI_PKAp, TrpnTotal, CmdnTotal, KmCmdn, KmTrpn, fracTnIp0),
    ]

    defaults = [Cai[1] => Cai_sub_SR_default, Cai[m] => Cai_sub_SL_default]

    for i in 2:m-1
        eq = D(Cai[i]) ~ (Dca / (j[i] * dx^2) * ((1 + j[i]) * Cai[i+1] - 2 * j[i] * Cai[i] + (j[i] - 1) * Cai[i-1])) * beta_cai(Cai[i]; TnI_PKAp, TrpnTotal, CmdnTotal, KmCmdn, KmTrpn, fracTnIp0)
        push!(eqs, eq)
        push!(defaults, Cai[i] => (Cai_sub_SR_default + Cai_sub_SL_default)/2)
    end

    return ODESystem(eqs, t; name, defaults)
end

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
function get_ncx_sys(nai, cai, nao, cao, vm, ICa_scale=1; name=:ncxsys)
    @parameters begin
        fNaCa = 1
        kNaCa = 2.2680e-016 * μA / cm^2 / μM^4
        dNaCa = 1e-16 / μM^4
        gamma = 0.5
    end
    @variables INaCa(t)
    return ODESystem([INaCa ~ ICa_scale * kNaCa * ((exp(iVT * gamma * vm) * nai^3 * cao - exp(iVT * (gamma - 1) * vm) * cai * nao^3 * fNaCa) / (1 + dNaCa * (nao^3 * cai * fNaCa + nai^3 * cao)))], t; name)
end

"L-type calcium current (LCC)"
function get_lcc_sys(cai, cao, vm, ICa_scale=1; name=:lccsys)
    @parameters GCaL = 1.3125e-4 * 0.8 * 0.6 / ms
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
    taufca = 10ms
    kfca = 1 - (fcainf > i_fca) * (V > -60)

    eqs = [
        ICaL ~ ICa_scale * i_d * i_f * i_fca * ghkVm(GCaL, vm, cai, 0.341 * cao, 2),
        D(i_d) * taud ~ dinf - i_d,
        D(i_f) * tauf ~ finf - i_f,
        D(i_fca) * taufca ~ kfca * (fcainf - i_fca),
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
        ginf ~ expit(-(V + 66) / 6)
        taug ~ (1 + 40 * expit(-(V + 65) * 0.08)) * ms,
        D(i_b) * taub ~ binf - i_b,
        D(i_g) * taug ~ ginf - i_g,
        ICaT ~ gCaT * i_b * i_g * (vm - E_Ca + 106.5mV),
    ]
    return ODESystem(eqs, t; name)
end

"Background current"
function get_ibg_sys(vm, E_Na, E_Ca; name=:ibgsys)
    @parameters gCab = 0.0008mS / cm^2
    @parameters gNab = 0.0026mS / cm^2
    @variables ICab(t) INab(t)
    eqs = [
        ICab ~ gCab * (vm - E_Ca),
        INab ~ gNab * (vm - E_Na),
    ]
    return ODESystem(eqs, t; name)
end

"Funny current (If)"
function get_if_sys(vm, E_Na, E_K; name=:ifsys)
    @parameters gf = 0.021mS / cm^2 fNa = 0.2
    @variables IfNa(t) IfK(t) If(t) yinf(t) tauy(t) i_y(t) = 0.07192
    fK = 1 - fNa
    V = vm * Volt / mV # Convert voltage to mV
    eqs = [
        yinf ~ expit(-(V + 78.65) / 6.33),
        tauy ~ 1 / (0.11885 * exp((V + 75) / 28.37) + 0.56236 * exp(-(V + 75) / 14.19)),
        IfNa ~ gf * fNa * i_y * (vm - E_Na),
        IfK ~ gf * fK * i_y * (vm - E_K),
        If ~ IfNa + IfK,
        D(i_y) * tauy ~ yinf - i_y
    ]
    return ODESystem(eqs, t; name)
end

"Fast sodium current (INa)"
function get_ina_sys(vm, E_Na; name=:inasys)
    @parameters gNa = 35mS / cm^2
    @variables begin
        INa(t)
        Naminf(t)
        Nahinf(t)
        Najinf(t)
        Nataum(t)
        Natauh(t)
        i_Nam(t) = 0.0250
        i_Nah(t) = 0.22242
        i_Naj(t) = 0.19081
    end

    V = vm * Volt / mV # Convert voltage to mV
    NatauhHI = 0.4537ms * expit((V + 10.66) / 11.1)
    NatauhLOW = 3.49ms / (0.135 * exp((V + 80) / -6.8) + 3.56 * exp(0.079V) + 3.1e5 * exp(0.35V))
    NataujHI = 11.63ms * (1 + exp(-0.1 * (V + 32))) / exp(-2.535e-7V)
    NataujLOW = 3.49ms / ((V + 37.78) / (1 + exp(0.311 * (V + 79.23))) * (-127140 * exp(0.2444V) - 3.474e-5 * exp(-0.04391V)) + 0.1212 * exp(-0.01052V) / (1 + exp(-0.1378 * (V + 40.14))))

    eqs = [
        Naminf ~ expit((V + 45) / 6.5),
        Nahinf ~ expit(-(V + 76.1) / 6.07),
        Najinf ~ Nahinf,
        Nataum ~ 1.36ms / (3.2 * exprel(-0.1 * (V + 47.13)) + 0.08 * exp(-V / 11)),
        Natauh ~ ifelse(V >= -40, NatauhHI, NatauhLOW),
        Natauj ~ ifelse(V >= -40, NataujHI, NataujLOW),
        INa ~ gNa * i_Nam^3 * i_Nah * i_Naj * (vm - E_Na),
        D(i_Nam) * Nataum ~ Naminf - i_Nam,
        D(i_Nah) * Natauh ~ Nahinf - i_Nah,
        D(i_Naj) * Natauj ~ Najinf - i_Naj,
    ]
    return ODESystem(eqs, t; name)
end

"Potassium currents"
function get_ik_eqs(vm, E_K, K_i, K_o, Na_i, Na_o, IKur_PKAp=0; name=:iksys)
    V = vm * Volt / mV # Convert voltage to mV

    # IK1: time-independent
    @parameters gK1 = 0.0515mS / cm^2 * hil(K_o, 210μM)
    @variables IK1(t)
    vk1 = vm - E_K - 6.1373mV

    # Ito: Where does this come from? Perhaps here: https://modeldb.science/262081
    @parameters gt = 0.1mS / cm^2 f_is = 0.706
    @variables begin
        Ito(t)
        i_r(t) = 0.00702
        i_s(t) = 0.9660
        i_sslow(t) = 0.22156
        sinf(t)
        rinf(t)
        slowinf(t)
        taur(t)
        taus(t)
        tausslow(t)
    end

    # IKs (and IKur)
    @parameters GKs = 0.05mS / cm^2
    @variables begin
        IKs(t)
        i_nKs(t) = 0.09243
        nKsinf(t)
        nKstau(t) = 750ms
    end

    # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    fracIKurp0 = 0.437635       # Derived quantity (IKur_PKAp(baseline)/IKurtot)
    fracIKurpISO = 0.718207     # Derived quantity (IKur_PKAp(ISO)/IKurtot)
    a_Kur = (1.20 - 1) / (fracIKurpISO / fracIKurp0 - 1)
    fracIKuravail = (1 - a_Kur) + a_Kur * (IKur_PKAp / fracIKurp0)  # +20# with 0.1 uM ISO
    alphan = 0.00000481333 * (V + 26.5) / (-expm1(-0.128 * (V + 26.5)))
    betan = 0.0000953333 * exp(-0.038 * (V + 26.5))

    # IKr
    @parameters begin
        GKr = 0.06mS / cm^2
        kf = 0.023761 / ms
        kb = 0.036778 / ms
    end

    @variables begin
        IKr(t)
        E_Kr(t)
        CK0(t)
        i_CK1(t) = 0.00188
        i_CK2(t) = 0.00977
        i_OK(t) = 0.26081
        i_IK(t) = 0.07831
    end

    alphaa0 = 0.022348 / ms * exp(0.01176 * V)
    betaa0 = 0.047002 / ms * exp(-0.0631 * V)
    alphaa1 = 0.013733 / ms * exp(0.038198 * V)
    betaa1 = 0.0000689 / ms * exp(-0.04178 * V)
    alphai_mERG = 0.090821 / ms * exp(0.023391 * V)
    betai_mERG = 0.006497 / ms * exp(-0.03268 * V)

    eqs = [
        IK1 ~ gK1 * vk1 / (0.1653 + exp(0.0319 / mV * vk1)),
        sinf ~ expit((vm + 31.97156mV) / -4.64291mV),
        rinf ~ expit((vm - 3.55716mV) / 14.61299mV),
        slowinf ~ sinf,
        taur ~ inv(45.16 * exp(0.03577 / mV * (vm + 50mV)) + 98.9 * exp(-0.1 / mV * (vm + 38mV))),
        taus ~ (0.35 * exp(-(((vm + 70mV) / 15mV)^2)) + 0.035) - 26.9ms,
        tausslow ~ (3.7 * exp(-(((vm + 70mV) / 30mV)^2)) + 0.035) + 37.4ms,
        Ito ~ gt * i_r * (f_is * i_s + (1 - f_is) * i_sslow) * (vm - E_K),
        D(i_r) * taur ~ rinf - i_r,
        D(i_s) * taus ~ sinf - i_s,
        D(i_sslow) * tausslow ~ slowinf - i_sslow,
        IKs ~ GKs * i_nKs^2 * (vm - E_K) * fracIKuravail * 2,
        nKsinf ~ alphan / (alphan + betan),
        D(i_nKs) * nKstau ~ nKsinf - i_nKs,
        E_Kr ~ nerst(0.98 * K_o + 0.02 * Na_o, 0.98 * K_i + 0.02 * Na_i),
        IKr ~ i_OK * GKr * (vm - E_Kr),
        CK0 ~ 1 - (i_CK1 + i_CK2 + i_OK + i_IK),
        D(i_CK1) ~ (alphaa0 * CK0 - betaa0 * i_CK1 + kb * i_CK2 - kf * i_CK1),
        D(i_CK2) ~ (kf * i_CK1 - kb * i_CK2 + betaa1 * i_OK - alphaa1 * i_CK2),
        D(i_OK) ~ (alphaa1 * i_CK2 - betaa1 * i_OK + betai_mERG * i_IK - alphai_mERG * i_OK),
        D(i_IK) ~ (alphai_mERG * i_OK - betai_mERG * i_IK),
    ]
    return ODESystem(eqs, t; name)
end

function get_ryr_sys(Ca, CaJSR; name=:ryrsys)
    @parameters begin
        nRyR = 4
        nu1RyR = 0.01 / ms
        kaposRyR = 1 / ms
        kanegRyR = 0.16 / ms
    end
    @variables begin
        i_PO1(t) = 0.0037
        i_PC1(t)
        Jrel(t)
    end
    KmRyR = (1.35 * 2.6 * expit(-(Ca - 530μM) / 200μM) + 1.5 - 0.9 - 0.3 - 0.05) * μM
    eqs = [
        1 ~ i_PO1 + i_PC1,
        Jrel ~ nu1RyR * i_PO1 * (CaJSR - Ca),
        D(PO1_RyR) ~ kaposRyR * hil(Cai_sub_SR, KmRyR, nRyR) * PC1_RyR - kanegRyR * PO1_RyR,
    ]
    return ODESystem(eqs, t; name)
end

function get_serca_sys(Cai, CaNSR, CaJSR, fracPLB_CKp=0, fracPLBp=0, RyR_CKp=0; name=:sercasys)
    @parameters begin
        VmaxfSR = 0.9996μM / ms
        VmaxrSR = VmaxfSR
        KmfSR = 0.5μM
        KmrSR = 7000 * KmfSR
        HfSR = 2
        HrSR = 1 * Hf
        kSRleak = 5e-6 / ms
        fracPKA_PLBo = 1 - 0.079755
    end

    @variables Jup(t) Jleak(t) Jtr(t)

    fCKII_PLB = (1 - 0.5 * fracPLB_CKp)  # Max effect: fCKII_PLB=0.5
    PLB_PKAn = 1 - fracPLBp
    fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * (100 - 55.31) / 100 + 55.31 / 100
    # Select smaller value (resulting in max reduction of Kmf)
    Kmfp = 2 * Kmf * min(fCKII_PLB, fPKA_PLB)  #fCKII_PLB
    fSR = NaNMath.pow(Cai / Kmfp, HfSR)
    rSR = NaNMath.pow(CaNSR / KmrSR, HrSR)

    kleak = (1 / 2 + 5 * RyR_CKp / 2) * kSRleak
    eqs = [
        Jup ~ (VmaxfSR * fSR - VmaxrSR * rSR) / (1 + fSR + rSR),
        Jleak ~ kleak * (CaNSR - Cai),
        Jtr ~ (CaNSR - CaJSR) / 200ms
    ]
    return ODESystem(eqs, t; name)
end

function build_neonatal_ecc_sys(;name=:neonataleccsys)
    @parameters begin
        Ca_o = 1796μM
        Na_o = 154578μM
        K_o = 5366μM
        Mg_i = 1000μM
        ROS = 0μM
        ISO = 0μM
        ATP = 5000μM
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

    capdesys = get_ca_pde_sys()
    @unpack Cai_sub_SL, Cai_sub_SR, Cai_mean = capdesys
    camkiisys = get_camkii_eqs(Cai_mean, ROS)
    barsys = get_bar_sys(ATP, ISO)
    @unpack LCCa_PKAp, LCCb_PKAp, fracPLBp, TnI_PKAp = barsys
    ICa_scale = get_ICa_scalep

    eqs = [
        E_Na ~ nernst(Na_o, Na_i),
        E_K ~ nernst(K_o, K_i),
        E_Ca ~ nernst(Ca_o, Cai_sub_SL, 2),
    ]
    return ODESystem(eqs, t; name)
end
