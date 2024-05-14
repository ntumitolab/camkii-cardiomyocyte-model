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
using NaNMath

"Calcium buffered by troponin and calmodulin"
function beta_cai(Ca; TnI_PKAp=0)

    @parameters begin
        TrpnTotal=35μM
        CmdnTotal=50μM
        KmCmdn=2.38μM
        KmTrpn=0.5μM
        fracTnIp0=0.062698 # Baseline effect
    end

    fPKA_TnI = 1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIp0) # Max effect +61%
    KmTrpnNew = KmTrpn / fPKA_TnI
    return inv(1 + TrpnTotal * KmTrpnNew / (Ca + KmTrpnNew)^2 + CmdnTotal * KmCmdn / (Ca + KmCmdn)^2)
end

"Calcium diffusion between sarcolemma (SL) and sarcoplasmic reticulum (SR)"
function get_ca_pde_eqs(;
    dx=0.1μm,
    rSR_true=6μm,
    rSL_true=10.5μm,
    V_sub_SR=4 / 3 * pi * ((rSR_true + dx)^3 - (rSR_true)^3),
    V_sub_SL=4 / 3 * pi * (rSL_true^3 - (rSL_true - dx)^3),
    TnI_PKAp=0,
)
    @variables t
    D = Differential(t)
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial indices
    m = length(j)
    @variables Cai(t)[1:m] Cai_mean(t) Cai_sub_SR(t) Cai_sub_SL(t) JCa_SR(t) JCa_SL(t)
    @parameters Dca = 7μm^2 / ms
    eqs = [
        Cai_mean ~ sum(collect(Cai)) / m,
        Cai_sub_SR ~ Cai[1],
        Cai_sub_SL ~ Cai[m],
        D(Cai[1]) ~ (Dca / (j[1] * dx^2) * ((1 + j[1]) * Cai[2] - 2 * j[1] * Cai[1] + (j[1] - 1) * Cai[1]) + JCa_SR / V_sub_SR) * beta_cai(Cai[1]; TnI_PKAp),
        D(Cai[m]) ~ (Dca / (j[m] * dx^2) * ((1 + j[m]) * Cai[m] - 2 * j[m] * Cai[m] + (j[m] - 1) * Cai[m-1]) + JCa_SL / V_sub_SL) * beta_cai(Cai[m]; TnI_PKAp),
    ]

    for i in 2:m-1
        eq = D(Cai[i]) ~ (Dca / (j[i] * dx^2) * ((1 + j[i]) * Cai[i+1] - 2 * j[i] * Cai[i] + (j[i] - 1) * Cai[i-1])) * beta_cai(Cai[i]; TnI_PKAp)
        push!(eqs, eq)
    end

    return eqs
end

function build_neonatal_ecc_eqs(;
    LCCb_PKAp=0,  # Fraction of LCC phosphorylated by PKAPLB_CKp
    PLBT17p=0,
    PLBp=0,
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
        Cai_sub_SR(t)
        E_Na(t)
        E_K(t)
        E_Ca(t)
        vm(t)
    end

    D = Differential(t)

    # NCX
    @parameters begin
        fNaCa = 1
        kNaCa = 2.2680e-016 * μA / cm^2 / μM^4
        dNaCa = 1e-16 / μM^4
        gamma = 0.5
    end
    @variables INaCa(t)

    # LCC
    @parameters GCaL = 1.3125e-4 * 0.8 * 0.6 / ms
    @variables ICaL(t) i_d(t) i_f(t) i_fca(t)

    # L-type calcium current (LCC)
    dinf = expit((vm + 11.1mV) / 7.2mV)
    alphad = 1.4 * expit((vm + 35mV) / 13mV) + 0.25
    betad = 1.4 * expit(-(vm + 5mV) / 5mV)
    gammad = expit((vm - 50mV) / 20mV)
    taud = (alphad * betad + gammad) * ms
    finf = expit(-(vm + 23.3mV) / 5.4mV)
    tauf = (1125 * exp(-(V + 27mV)^2 / 240mV^2) + 165 * expit((vm - 25mV) / 10mV) + 120) * ms
    alphafca = hilr(Cai_sub_SL, 0.000325 * 1.5, 8)
    betafca = 0.1 * expit(-(Cai_sub_SL - 0.0005mM) / 0.0001mM)
    gammafca = 0.2 * expit(-(Cai_sub_SL - 0.00075mM) / 0.0008mM)
    fcainf = (alphafca + betafca + gammafca + 0.23) / 1.46
    taufca = 10ms
    kfca = 1 - (fcainf > i_fca) * (vm > -60mV)

    # T-Type calcium current (TCC)
    @parameters gCaT = 0.2mS / cm^2
    @variables ICaT(t) i_b(t) i_g(t)
    binf = expit((vm + 37.49098mV) / 5.40634mV)
    taub = (0.6 + 5.4 * expit(-(vm + 100mV) * 0.03 / mV)) * ms
    ginf = expit(-(vm + 66mV) / 6mV)
    taug = (1 + 40 * expit(-(vm + 65mV) * 0.08 / mV)) * ms

    # Background calcium and sodium currents
    @parameters gCab = 0.0008mS / cm^2
    @parameters gNab = 0.0026mS / cm^2
    @variables ICab(t) INab(t)

    # If
    @parameters gf = 0.021mS / cm^2 fNa = 0.2
    @variables IfNa(t) IfK(t) If(t) i_y(t)
    fK = 1 - fNa
    yinf = expit(-(vm + 78.65mV) / 6.33mV)
    tauy = expit(-(vm + 75mV) / 28.37mV) / 0.11885 + 0.56236 * exp((vm + 75mV) / 14.19mV)

    # IK1
    @parameters gK1 = 0.0515mS / cm^2 * hil(K_o, 210μM)
    @variables IK1(t)
    vk1 = vm - E_K - 6.1373mV

    # Ito
    # Where does this come from? Perhaps here: https://modeldb.science/262081
    @parameters gt = 0.1mS / cm^2
    @variables Ito(t) i_r(t) i_s(t) i_sslow(t)
    sinf = expit(-(vm + 31.97156mV) / 4.64291mV)
    rinf = expit((vm - 3.55716mV) / 14.61299mV)
    slowinf = sinf
    taur = inv(45.16 * exp(0.03577 / mV * (vm + 50mV)) + 98.9 * exp(-0.1 / mV * (vm + 38mV)))
    taus = (0.35 * exp(-(((vm + 70mV) / 15mV)^2)) + 0.035) - 26.9ms
    tausslow = (3.7 * exp(-(((vm + 70mV) / 30mV)^2)) + 0.035) + 37.4ms

    # INa
    @parameters gNa = 35mS / cm^2
    @variables INa(t) i_Nam(t) i_Nah(t) i_Naj(t)
    Naminf = expit((vm + 45mV) / 6.5mV)
    Nahinf = expit(-(vm + 76.1mV) / 6.07mV)
    Najinf = Nahinf
    Nataum = 0.00136second / (3.2 * exprel(-0.1/mV *(vm + 47.13mV))+0.08*exp(-vm/11mV))
    NatauhHI = 0.0004537second * expit(( vm + 10.66mV)/11.1mV)
    NatauhLOW = 0.00349second /( 0.135 * exp((vm + 80mV) / -6.8mV) + 3.56 * exp(0.079/mV * vm) + 3.1e5 * exp(0.35/mV * vm))
    Natauh = ifelse(vm >= -40mV, NatauhHI, NatauhLOW)
    NataujHI = 0.01163second * (1 + exp(-0.1/mV * (vm + 32mV))) / exp(-2.535e-7/mV * vm)
    NataujLOW = 0.00349second / ((vm + 37.78mV)/mV / (1 + exp(0.311/mV * (vm + 79.23mV))) * (-127140 * exp(0.2444/mV * vm) - 3.474e-5 * exp(-0.04391/mV * vm)) + 0.1212 * exp(-0.01052/mV * vm) / (1 + exp(-0.1378/mV * (vm + 40.14mV))))
    Natauj = ifelse(vm >= -40mV, NataujHI, NataujLOW)

    # RyR
    @parameters begin
        nRyR = 4
        nu1RyR = 0.01/ms
        kaposRyR = 3/ms
        kanegRyR = 0.48/ms
    end
    @variables Jrel(t) Ca_JSR(t) Cai_sub_SR(t) PO1_RyR(t) PC1_RyR(t) KmRyR(t)

    KmRyR = (1.35 * 2.6 * expit(-(Ca_JSR - 530μM) / 200μM) + 1.5 - 0.9 - 0.3 - 0.05) * μM

    # SERCA and PLB
    @parameters begin
        VmaxfSR = 0.9996μM/ms
        VmaxrSR = VmaxfSR
        KmfSR = 0.5μM
        KmrSR = 7000 * KmfSR
        HfSR = 2
        HrSR = 1 * Hf
        kSRleak = 5e-6/ms
        fracPKA_PLBo = 1 - 0.079755
        PLBtot = 106μM
    end

    @variables Jup(t) Cai_sub_SR(t) CaNSR(t)

    PLB_CKp = PLBT17p / PLBtot
    fCKII_PLB = (1 - 0.5 * PLB_CKp)  # Max effect: fCKII_PLB=0.5
    PLB_PKAn = (PLBtot - PLBp) / PLBtot
    fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * (100 - 55.31) / 100 + 55.31 / 100
    # Phosphorylation will decrease KmfSR (increase affinity)
    Kmfp = KmfSR * min(fCKII_PLB, fPKA_PLB)
    fSR = NaNMath.pow(Cai_sub_SR / Kmfp, HfSR)
    rSR = NaNMath.pow(CaNSR / KmrSR, HrSR)


    eqs = [
        E_Na ~ nernst(Na_o, Na_i),
        E_K ~ nernst(K_o, K_i),
        E_Ca ~ nernst(Ca_o, Cai_sub_SL, 2),
        INaCa ~ ICa_scalep * kNaCa * ((exp(iVT * gamma * vm) * Na_i^3 * Ca_o - exp(iVT * (gamma - 1) * vm) * Cai_sub_SL * Na_o^3 * fNaCa) / (1 + dNaCa * (Na_o^3 * Cai_sub_SL * fNaCa + Na_i^3 * Ca_o))),
        ICaL ~ ICa_scalep * i_d * i_f * i_fca * ghkVm(GCaL, vm, Cai_sub_SL, 0.341 * Ca_o, 2),
        D(i_d) * taud ~ dinf - i_d,
        D(i_f) * tauf ~ finf - i_f,
        D(i_fca) * taufca ~ kfca * (fcainf - i_fca),
        ICaT ~ gCaT * i_b * i_g * (vm - E_Ca + 106.5mV),
        D(i_b) * taub ~ binf - i_b,
        D(i_g) * taug ~ ginf - i_g,
        ICab ~ gCab * (V - E_Ca),
        INab ~ gNab * (V - E_Na),
        IfNa ~ gf * fNa * i_y * (vm - E_Na),
        IfK ~ gf * fK * i_y * (vm - E_K),
        If ~ IfNa + IfK,
        D(i_y) * tauy ~ yinf - i_y,
        IK1 ~ gK1 * vk1 / (0.1653 + exp(0.0319 / mV * vk1)),
        Ito ~ gt * i_r * (0.706 * i_s + 0.294 * i_sslow) * (vm - E_K),
        D(i_r) * taur ~ rinf - i_r,
        D(i_s) * taus ~ sinf - i_s,
        D(i_sslow) * tausslow ~ slowinf - i_sslow,
        INa ~ gNa * i_Nam^3 * i_Nah * i_Naj * (vm - E_Na)
        D(i_Nam) ~ (Naminf - i_Nam) / Nataum,
        D(i_Nah) ~ (Nahinf - i_Nah) / Natauh,
        D(i_Naj) ~ (Najinf - i_Naj) / Natauj,
        Jrel ~ nu1RyR * PO1_RyR * (Ca_JSR - Cai_sub_SR),
        1 ~ PO1_RyR + PC1_RyR,
        D(PO1_RyR) ~ kaposRyR * hil(Cai_sub_SR, KmRyR, nRyR) * PC1_RyR - kanegRyR * PO1_RyR,
        Jup ~ (VmaxfSR * fSR - VmaxrSR * rSR) / (1 + fSR + rSR),

    ]

    return eqs
end
