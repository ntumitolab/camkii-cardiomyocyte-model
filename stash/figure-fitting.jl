using Catalyst
using DifferentialEquations
using Plots
using ModelingToolkit
using Sundials
using Statistics
using XLSX
import NaNMath as nm
Plots.default(fmt=:png)

function beta_cai(Cai, TnI_PKAp)
    fracTnIpo = 0.062698  # Derived quantity (TnI_PKAp(baseline)/TnItot)
    fPKA_TnI = (1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIpo)) # Max effect +61#
    TrpnTotal = 35
    KmTrpn = 0.5 / fPKA_TnI
    CmdnTotal = 50
    KmCmdn = 2.38
    return inv(1 + TrpnTotal * KmTrpn / (Cai + KmTrpn)^2 + CmdnTotal * KmCmdn / (Cai + KmCmdn)^2)
end

function get_Morotti_equations()
    @variables t
    D = Differential(t)

    # CaMDyad / Cyt / SL variables
    #@variables CaM_dyad(t) Ca2CaM_dyad(t) Ca4CaM_dyad(t) CaMB_dyad(t) Ca2CaMB_dyad(t) Ca4CaMB_dyad(t) Pb2_dyad(t) Pb_dyad(t) Pt_dyad(t) Pt2_dyad(t)
    #@variables Pa_dyad(t) Ca4CaN_dyad(t) CaMCa4CaN_dyad(t) Ca2CaMCa4CaN_dyad(t) Ca4CaMCa4CaN_dyad(t)
    @variables CaM_sl(t) Ca2CaM_sl(t) Ca4CaM_sl(t) CaMB_sl(t) Ca2CaMB_sl(t) Ca4CaMB_sl(t) Pb2_sl(t) Pb_sl(t) Pt_sl(t) Pt2_sl(t)
    @variables Pa_sl(t) Ca4CaN_sl(t) CaMCa4CaN_sl(t) Ca2CaMCa4CaN_sl(t) Ca4CaMCa4CaN_sl(t)
    @variables CaM_cyt(t) Ca2CaM_cyt(t) Ca4CaM_cyt(t) CaMB_cyt(t) Ca2CaMB_cyt(t) Ca4CaMB_cyt(t) Pb2_cyt(t) Pb_cyt(t) Pt_cyt(t) Pt2_cyt(t)
    @variables Pa_cyt(t) Ca4CaN_cyt(t) CaMCa4CaN_cyt(t) Ca2CaMCa4CaN_cyt(t) Ca4CaMCa4CaN_cyt(t)

    # CaMKII variables
    @variables LCC_PKAp(t) RyR2809p(t) RyR2815p(t) PLBT17p(t) LCC_CKslp(t) #LCC_CKdyadp(t)
    @variables I1p_PP1(t) Pb_sl(t) Pt_sl(t) Pt2_sl(t) Pa_sl(t) #Pb_dyad(t) Pt_dyad(t) Pt2_dyad(t) Pa_dyad(t)

    #BAR variables
    @variables LR(t) LRG(t) RG(t) b1AR_S464(t) b1AR_S301(t) GsaGTPtot(t) GsaGDP(t) Gsby(t) AC_GsaGTP(t) PDEp(t)
    @variables cAMPtot(t) RC_I(t) RCcAMP_I(t) RCcAMPcAMP_I(t) RcAMPcAMP_I(t) PKACI(t) PKACI_PKI(t) RC_II(t) RCcAMP_II(t) RCcAMPcAMP_II(t)
    @variables RcAMPcAMP_II(t) PKACII(t) PKACII_PKI(t) I1p_PP1(t) I1ptot(t) PLBp(t) PLMp(t) LCCap(t) LCCbp(t) RyRp(t)
    @variables TnIp(t) KS79(t) KS80(t) KSp(t) CFTRp(t) KURp(t)

    # ECC variables
    @variables i_CaJSR(t) CaNSR(t) V(t) Nai(t) Ki(t) i_b(t) i_g(t) i_d(t) i_f(t) i_fca(t) i_y(t) i_r(t) i_s(t) i_sslow(t)
    @variables i_Nam(t) i_Nah(t) i_Naj(t) i_PO1(t) i_PO2(t) i_PC2(t) i_nKs(t) i_CK1(t) i_CK2(t) i_OK(t) i_IK(t)
    #@variables Cai_mean(t)

    ## Neonatal Rat
    # Cell geometry
    rSR_true = 6
    rSL_true = 10.5

    # Ca diffusion grid
    dx = 0.1
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial index of Cai diffusion
    # Cai
    m = length(j)
    @variables Cai(t)[1:m]

    Cai_sub_SR = Cai[1]
    Cai_sub_SL = Cai[m]
    Cai_mean = mean(skipmissing(Cai))


    ## CamDyad/ CaMSL/ CamCyt
    ## Parameters
    Mg = 1      # [mM]
    K = 135     # [mM]s
    Btot_dyad = 0
    CaMKIItot_dyad = 120         # [uM]
    CaNtot_dyad = 3e-3 / 8.293e-4  # [uM]
    PP1tot_dyad = 96.5           # [uM]

    ## Parameters for Cyt and SL
    Btot = 24.2     # [uM]
    CaMKIItot = 120 * 8.293e-4  # [uM]
    CaNtot = 3e-3             # [uM]
    PP1tot = 0.57             # [uM]

    ## Parameters
    # Ca/CaM parameters
    if Mg <= 1
        Kd02 = 0.0025 * (1 + K / 0.94 - Mg / 0.012) * (1 + K / 8.1 + Mg / 0.022)  # [uM^2]
        Kd24 = 0.128 * (1 + K / 0.64 + Mg / 0.0014) * (1 + K / 13.0 - Mg / 0.153) # [uM^2]
    else
        Kd02 = 0.0025 * (1 + K / 0.94 - 1 / 0.012 + (Mg - 1) / 0.060) * (1 + K / 8.1 + 1 / 0.022 + (Mg - 1) / 0.068)   # [uM^2]
        Kd24 = 0.128 * (1 + K / 0.64 + 1 / 0.0014 + (Mg - 1) / 0.005) * (1 + K / 13.0 - 1 / 0.153 + (Mg - 1) / 0.150)  # [uM^2]
    end
    k20 = 10               # [s^-1]
    k02 = k20 / Kd02         # [uM^-2 s^-1]
    k42 = 500              # [s^-1]
    k24 = k42 / Kd24         # [uM^-2 s^-1]

    # CaM buffering (B) parameters
    k0Boff = 0.0014        # [s^-1]
    k0Bon = k0Boff / 0.2   # [uM^-1 s^-1] kon = koff/Kd
    k2Boff = k0Boff / 100    # [s^-1]
    k2Bon = k0Bon          # [uM^-1 s^-1]
    k4Boff = k2Boff        # [s^-1]
    k4Bon = k0Bon          # [uM^-1 s^-1]

    # using thermodynamic constraints
    k20B = k20 / 100 # [s^-1] thermo constraint on loop 1
    k02B = k02     # [uM^-2 s^-1]
    k42B = k42     # [s^-1] thermo constraint on loop 2
    k24B = k24     # [uM^-2 s^-1]

    # CaMKII parameters
    # Wi Wa Wt Wp
    kbi = 2.2      # [s^-1] (Ca4CaM dissocation from Wb)
    kib = kbi / 33.5e-3 # [uM^-1 s^-1]
    kib2 = kib
    kb2i = kib2 * 5
    kb24 = k24
    kb42 = k42 * 33.5e-3 / 5
    kpp1 = 1.72     # [s^-1] (PP1-dep dephosphorylation rates)
    Kmpp1 = 11.5    # [uM]
    kta = kbi / 1000  # [s^-1] (Ca4CaM dissociation from Wt)
    kat = kib       # [uM^-1 s^-1] (Ca4CaM reassociation with Wa)
    kt42 = k42 * 33.5e-6 / 5
    kt24 = k24
    kat2 = kib
    kt2a = kib * 5

    # CaN parameters
    kcanCaoff = 1              # [s^-1]
    kcanCaon = kcanCaoff / 0.5               # [uM^-1 s^-1]
    kcanCaM4on = 46            # [uM^-1 s^-1]
    kcanCaM4off = 1.3e-3       # [s^-1]
    kcanCaM2on = kcanCaM4on
    kcanCaM2off = 2508 * kcanCaM4off
    kcanCaM0on = kcanCaM4on
    kcanCaM0off = 165 * kcanCaM2off
    k02can = k02
    k20can = k20 / 165
    k24can = k24
    k42can = k20 / 2508

    #=
    ## Dyad Fluxes

    Ca_Dyad = Ca_j*100
    # CaM Reaction fluxes
    B_dyad = Btot_dyad - CaMB_dyad - Ca2CaMB_dyad - Ca4CaMB_dyad
    rcn02_dyad = k02*Ca_Dyad^2*CaM_dyad - k20*Ca2CaM_dyad
    rcn24_dyad = k24*Ca_Dyad^2*Ca2CaM_dyad - k42*Ca4CaM_dyad
    # CaM buffer fluxes
    rcn02B_dyad = k02B*Ca_Dyad^2*CaMB_dyad - k20B*Ca2CaMB_dyad
    rcn24B_dyad = k24B*Ca_Dyad^2*Ca2CaMB_dyad - k42B*Ca4CaMB_dyad
    rcn0B_dyad = k0Bon*CaM_dyad*B_dyad - k0Boff*CaMB_dyad
    rcn2B_dyad = k2Bon*Ca2CaM_dyad*B_dyad - k2Boff*Ca2CaMB_dyad
    rcn4B_dyad = k4Bon*Ca4CaM_dyad*B_dyad - k4Boff*Ca4CaMB_dyad
    # CaN reaction fluxes
    Ca2CaN_dyad = CaNtot_dyad - Ca4CaN_dyad - CaMCa4CaN_dyad - Ca2CaMCa4CaN_dyad - Ca4CaMCa4CaN_dyad
    rcnCa4CaN_dyad = kcanCaon*Ca_Dyad^2*Ca2CaN_dyad - kcanCaoff*Ca4CaN_dyad
    rcn02CaN_dyad = k02can*Ca_Dyad^2*CaMCa4CaN_dyad - k20can*Ca2CaMCa4CaN_dyad
    rcn24CaN_dyad = k24can*Ca_Dyad^2*Ca2CaMCa4CaN_dyad - k42can*Ca4CaMCa4CaN_dyad
    rcn0CaN_dyad = kcanCaM0on*CaM_dyad*Ca4CaN_dyad - kcanCaM0off*CaMCa4CaN_dyad
    rcn2CaN_dyad = kcanCaM2on*Ca2CaM_dyad*Ca4CaN_dyad - kcanCaM2off*Ca2CaMCa4CaN_dyad
    rcn4CaN_dyad = kcanCaM4on*Ca4CaM_dyad*Ca4CaN_dyad - kcanCaM4off*Ca4CaMCa4CaN_dyad
    # CaMKII reaction fluxes
    Pi_dyad = 1 - Pb2_dyad - Pb_dyad - Pt_dyad - Pt2_dyad - Pa_dyad
    rcnCKib2_dyad = kib2*Ca2CaM_dyad*Pi_dyad - kb2i*Pb2_dyad
    rcnCKb2b_dyad = kb24*Ca_Dyad^2*Pb2_dyad - kb42*Pb_dyad
    rcnCKib_dyad = kib*Ca4CaM_dyad*Pi_dyad - kbi*Pb_dyad
    T_dyad = Pb_dyad + Pt_dyad + Pt2_dyad + Pa_dyad
    kbt_dyad = 0.055*T_dyad + 0.0074*T_dyad^2 + 0.015*T_dyad^3
    rcnCKbt_dyad = kbt_dyad*Pb_dyad - kpp1*PP1tot_dyad*Pt_dyad/(Kmpp1+CaMKIItot_dyad*Pt_dyad)
    rcnCKtt2_dyad = kt42*Pt_dyad - kt24*Ca_Dyad^2*Pt2_dyad
    rcnCKta_dyad = kta*Pt_dyad - kat*Ca4CaM_dyad*Pa_dyad
    rcnCKt2a_dyad = kt2a*Pt2_dyad - kat2*Ca2CaM_dyad*Pa_dyad
    rcnCKt2b2_dyad = kpp1*PP1tot_dyad*Pt2_dyad/(Kmpp1+CaMKIItot_dyad*Pt2_dyad)
    rcnCKai_dyad = kpp1*PP1tot_dyad*Pa_dyad/(Kmpp1+CaMKIItot_dyad*Pa_dyad)
    =#

    ## SL Fluxes

    Ca_SL = Cai_sub_SL
    # CaM Reaction fluxes
    B_sl = Btot - CaMB_sl - Ca2CaMB_sl - Ca4CaMB_sl
    rcn02_sl = k02 * Ca_SL^2 * CaM_sl - k20 * Ca2CaM_sl
    rcn24_sl = k24 * Ca_SL^2 * Ca2CaM_sl - k42 * Ca4CaM_sl
    # CaM buffer fluxes
    rcn02B_sl = k02B * Ca_SL^2 * CaMB_sl - k20B * Ca2CaMB_sl
    rcn24B_sl = k24B * Ca_SL^2 * Ca2CaMB_sl - k42B * Ca4CaMB_sl
    rcn0B_sl = k0Bon * CaM_sl * B_sl - k0Boff * CaMB_sl
    rcn2B_sl = k2Bon * Ca2CaM_sl * B_sl - k2Boff * Ca2CaMB_sl
    rcn4B_sl = k4Bon * Ca4CaM_sl * B_sl - k4Boff * Ca4CaMB_sl
    # CaN reaction fluxes
    Ca2CaN_sl = CaNtot - Ca4CaN_sl - CaMCa4CaN_sl - Ca2CaMCa4CaN_sl - Ca4CaMCa4CaN_sl
    rcnCa4CaN_sl = kcanCaon * Ca_SL^2 * Ca2CaN_sl - kcanCaoff * Ca4CaN_sl
    rcn02CaN_sl = k02can * Ca_SL^2 * CaMCa4CaN_sl - k20can * Ca2CaMCa4CaN_sl
    rcn24CaN_sl = k24can * Ca_SL^2 * Ca2CaMCa4CaN_sl - k42can * Ca4CaMCa4CaN_sl
    rcn0CaN_sl = kcanCaM0on * CaM_sl * Ca4CaN_sl - kcanCaM0off * CaMCa4CaN_sl
    rcn2CaN_sl = kcanCaM2on * Ca2CaM_sl * Ca4CaN_sl - kcanCaM2off * Ca2CaMCa4CaN_sl
    rcn4CaN_sl = kcanCaM4on * Ca4CaM_sl * Ca4CaN_sl - kcanCaM4off * Ca4CaMCa4CaN_sl
    # CaMKII reaction fluxes
    Pi_sl = 1 - Pb2_sl - Pb_sl - Pt_sl - Pt2_sl - Pa_sl
    rcnCKib2_sl = kib2 * Ca2CaM_sl * Pi_sl - kb2i * Pb2_sl
    rcnCKb2b_sl = kb24 * Ca_SL^2 * Pb2_sl - kb42 * Pb_sl
    rcnCKib_sl = kib * Ca4CaM_sl * Pi_sl - kbi * Pb_sl
    T_sl = Pb_sl + Pt_sl + Pt2_sl + Pa_sl
    kbt_sl = 0.055 * T_sl + 0.0074 * T_sl^2 + 0.015 * T_sl^3
    rcnCKbt_sl = kbt_sl * Pb_sl - kpp1 * PP1tot * Pt_sl / (Kmpp1 + CaMKIItot * Pt_sl)
    rcnCKtt2_sl = kt42 * Pt_sl - kt24 * Ca_SL^2 * Pt2_sl
    rcnCKta_sl = kta * Pt_sl - kat * Ca4CaM_sl * Pa_sl
    rcnCKt2a_sl = kt2a * Pt2_sl - kat2 * Ca2CaM_sl * Pa_sl
    rcnCKt2b2_sl = kpp1 * PP1tot * Pt2_sl / (Kmpp1 + CaMKIItot * Pt2_sl)
    rcnCKai_sl = kpp1 * PP1tot * Pa_sl / (Kmpp1 + CaMKIItot * Pa_sl)

    ## Cyt Fluxes

    Ca_Cyt = Cai_mean
    # CaM Reaction fluxes
    B_cyt = Btot - CaMB_cyt - Ca2CaMB_cyt - Ca4CaMB_cyt
    rcn02_cyt = k02 * Ca_Cyt^2 * CaM_cyt - k20 * Ca2CaM_cyt
    rcn24_cyt = k24 * Ca_Cyt^2 * Ca2CaM_cyt - k42 * Ca4CaM_cyt
    # CaM buffer fluxes
    rcn02B_cyt = k02B * Ca_Cyt^2 * CaMB_cyt - k20B * Ca2CaMB_cyt
    rcn24B_cyt = k24B * Ca_Cyt^2 * Ca2CaMB_cyt - k42B * Ca4CaMB_cyt
    rcn0B_cyt = k0Bon * CaM_cyt * B_cyt - k0Boff * CaMB_cyt
    rcn2B_cyt = k2Bon * Ca2CaM_cyt * B_cyt - k2Boff * Ca2CaMB_cyt
    rcn4B_cyt = k4Bon * Ca4CaM_cyt * B_cyt - k4Boff * Ca4CaMB_cyt
    # CaN reaction fluxes
    Ca2CaN_cyt = CaNtot - Ca4CaN_cyt - CaMCa4CaN_cyt - Ca2CaMCa4CaN_cyt - Ca4CaMCa4CaN_cyt
    rcnCa4CaN_cyt = kcanCaon * Ca_Cyt^2 * Ca2CaN_cyt - kcanCaoff * Ca4CaN_cyt
    rcn02CaN_cyt = k02can * Ca_Cyt^2 * CaMCa4CaN_cyt - k20can * Ca2CaMCa4CaN_cyt
    rcn24CaN_cyt = k24can * Ca_Cyt^2 * Ca2CaMCa4CaN_cyt - k42can * Ca4CaMCa4CaN_cyt
    rcn0CaN_cyt = kcanCaM0on * CaM_cyt * Ca4CaN_cyt - kcanCaM0off * CaMCa4CaN_cyt
    rcn2CaN_cyt = kcanCaM2on * Ca2CaM_cyt * Ca4CaN_cyt - kcanCaM2off * Ca2CaMCa4CaN_cyt
    rcn4CaN_cyt = kcanCaM4on * Ca4CaM_cyt * Ca4CaN_cyt - kcanCaM4off * Ca4CaMCa4CaN_cyt
    # CaMKII reaction fluxes
    Pi_cyt = 1 - Pb2_cyt - Pb_cyt - Pt_cyt - Pt2_cyt - Pa_cyt
    rcnCKib2_cyt = kib2 * Ca2CaM_cyt * Pi_cyt - kb2i * Pb2_cyt
    rcnCKb2b_cyt = kb24 * Ca_Cyt^2 * Pb2_cyt - kb42 * Pb_cyt
    rcnCKib_cyt = kib * Ca4CaM_cyt * Pi_cyt - kbi * Pb_cyt
    T_cyt = Pb_cyt + Pt_cyt + Pt2_cyt + Pa_cyt
    kbt_cyt = 0.055 * T_cyt + 0.0074 * T_cyt^2 + 0.015 * T_cyt^3
    #rcnCKbt_cyt = Pb_cyt - kpp1*PP1tot*Pt_cyt/(Kmpp1+CaMKIItot*Pt_cyt)
    rcnCKbt_cyt = kbt_cyt * Pb_cyt - kpp1 * PP1tot * Pt_cyt / (Kmpp1 + CaMKIItot * Pt_cyt)
    rcnCKtt2_cyt = kt42 * Pt_cyt - kt24 * Ca_Cyt^2 * Pt2_cyt
    rcnCKta_cyt = kta * Pt_cyt - kat * Ca4CaM_cyt * Pa_cyt
    rcnCKt2a_cyt = kt2a * Pt2_cyt - kat2 * Ca2CaM_cyt * Pa_cyt
    rcnCKt2b2_cyt = kpp1 * PP1tot * Pt2_cyt / (Kmpp1 + CaMKIItot * Pt2_cyt)
    rcnCKai_cyt = kpp1 * PP1tot * Pa_cyt / (Kmpp1 + CaMKIItot * Pa_cyt)


    ## Ordinary Differential Equations
    Vmyo = 2.1454e-11           # [L]
    Vdyad = 1.7790e-014         # [L]
    VSL = 6.6013e-013           # [L]
    kSLmyo = 8.587e-15          # [L/msec]
    CaMKIItotDyad = 120         # [uM]
    BtotDyad = 1.54 / 8.293e-4    # [uM]
    #CaMtotDyad = CaM_dyad+Ca2CaM_dyad+Ca4CaM_dyad+CaMB_dyad+Ca2CaMB_dyad+Ca4CaMB_dyad+CaMKIItotDyad*(Pb2_dyad+Pb_dyad+Pt_dyad+Pt2_dyad)+CaMCa4CaN_dyad+Ca2CaMCa4CaN_dyad+Ca4CaMCa4CaN_dyad
    #Bdyad = BtotDyad - CaMtotDyad                                                  # [uM dyad]
    #J_cam_dyadSL = 1e-3*(k0Boff*CaM_dyad - k0Bon*Bdyad*CaM_sl)                     # [uM/msec dyad]
    #J_ca2cam_dyadSL = 1e-3*(k2Boff*Ca2CaM_dyad - k2Bon*Bdyad*Ca2CaM_sl)            # [uM/msec dyad]
    #J_ca4cam_dyadSL = 1e-3*(k2Boff*Ca4CaM_dyad - k4Bon*Bdyad*Ca4CaM_sl)            # [uM/msec dyad]
    J_cam_SLmyo = kSLmyo * (CaM_sl - CaM_cyt)                                          # [umol/msec]
    J_ca2cam_SLmyo = kSLmyo * (Ca2CaM_sl - Ca2CaM_cyt)                                 # [umol/msec]
    J_ca4cam_SLmyo = kSLmyo * (Ca4CaM_sl - Ca4CaM_cyt)                                 # [umol/msec]

    #=
    # CaMDyad equations
    CaMdyad_eqs = [
        D(CaM_dyad) ~ (1e-3*(-rcn02_dyad - rcn0B_dyad - rcn0CaN_dyad)-J_cam_dyadSL),                                                                      # du[1]
        D(Ca2CaM_dyad) ~ (1e-3*(rcn02_dyad - rcn24_dyad - rcn2B_dyad - rcn2CaN_dyad + CaMKIItot_dyad.*(-rcnCKib2_dyad + rcnCKt2a_dyad))-J_ca2cam_dyadSL),      # du[2]
        D(Ca4CaM_dyad) ~ (1e-3*(rcn24_dyad - rcn4B_dyad - rcn4CaN_dyad + CaMKIItot_dyad.*(-rcnCKib_dyad+rcnCKta_dyad))-J_ca4cam_dyadSL),                       # du[3]
        D(CaMB_dyad) ~ 1e-3*(rcn0B_dyad-rcn02B_dyad),                       # du[4]
        D(Ca2CaMB_dyad) ~ 1e-3*(rcn02B_dyad + rcn2B_dyad - rcn24B_dyad),    # du[5]
        D(Ca4CaMB_dyad) ~ 1e-3*(rcn24B_dyad + rcn4B_dyad),                  # du[6]
        # CaMKII equations
        D(Pb2_dyad) ~ 1e-3*(rcnCKib2_dyad - rcnCKb2b_dyad + rcnCKt2b2_dyad),     # du[7]
        D(Pb_dyad) ~ 1e-3*(rcnCKib_dyad + rcnCKb2b_dyad - rcnCKbt_dyad),         # du[8]
        D(Pt_dyad) ~ 1e-3*(rcnCKbt_dyad-rcnCKta_dyad-rcnCKtt2_dyad),             # du[9]
        D(Pt2_dyad) ~ 1e-3*(rcnCKtt2_dyad-rcnCKt2a_dyad-rcnCKt2b2_dyad),         # du[10]
        D(Pa_dyad) ~ 1e-3*(rcnCKta_dyad+rcnCKt2a_dyad-rcnCKai_dyad),             # du[11]
        # CaN equations
        D(Ca4CaN_dyad) ~ 1e-3*(rcnCa4CaN_dyad - rcn0CaN_dyad - rcn2CaN_dyad - rcn4CaN_dyad),      # du[12]
        D(CaMCa4CaN_dyad) ~ 1e-3*(rcn0CaN_dyad - rcn02CaN_dyad),                        # du[13]
        D(Ca2CaMCa4CaN_dyad) ~ 1e-3*(rcn2CaN_dyad+rcn02CaN_dyad-rcn24CaN_dyad),              # du[14]
        D(Ca4CaMCa4CaN_dyad) ~ 1e-3*(rcn4CaN_dyad+rcn24CaN_dyad)                       # du[15]
    ]
    =#
    # CaMSL equations
    CaMSL_eqs = [
        D(CaM_sl) ~ (1e-3 * (-rcn02_sl - rcn0B_sl - rcn0CaN_sl) - J_cam_SLmyo / VSL), #+ J_cam_dyadSL*Vdyad/VSL                                                                    # du[1]
        D(Ca2CaM_sl) ~ (1e-3 * (rcn02_sl - rcn24_sl - rcn2B_sl - rcn2CaN_sl + CaMKIItot .* (-rcnCKib2_sl + rcnCKt2a_sl)) - J_ca2cam_SLmyo / VSL), #+ J_ca2cam_dyadSL*Vdyad/VSL       # du[2]
        D(Ca4CaM_sl) ~ (1e-3 * (rcn24_sl - rcn4B_sl - rcn4CaN_sl + CaMKIItot .* (-rcnCKib_sl + rcnCKta_sl)) - J_ca4cam_SLmyo / VSL),  #+ J_ca4cam_dyadSL*Vdyad/VSL                     # du[3]
        D(CaMB_sl) ~ 1e-3 * (rcn0B_sl - rcn02B_sl),                     # du[4]
        D(Ca2CaMB_sl) ~ 1e-3 * (rcn02B_sl + rcn2B_sl - rcn24B_sl),    # du[5]
        D(Ca4CaMB_sl) ~ 1e-3 * (rcn24B_sl + rcn4B_sl),                # du[6]
        # CaMKII equations
        D(Pb2_sl) ~ 1e-3 * (rcnCKib2_sl - rcnCKb2b_sl + rcnCKt2b2_sl),     # du[7]
        D(Pb_sl) ~ 1e-3 * (rcnCKib_sl + rcnCKb2b_sl - rcnCKbt_sl),         # du[8]
        D(Pt_sl) ~ 1e-3 * (rcnCKbt_sl - rcnCKta_sl - rcnCKtt2_sl),             # du[9]
        D(Pt2_sl) ~ 1e-3 * (rcnCKtt2_sl - rcnCKt2a_sl - rcnCKt2b2_sl),         # du[10]
        D(Pa_sl) ~ 1e-3 * (rcnCKta_sl + rcnCKt2a_sl - rcnCKai_sl),             # du[11]
        # CaN equations
        D(Ca4CaN_sl) ~ 1e-3 * (rcnCa4CaN_sl - rcn0CaN_sl - rcn2CaN_sl - rcn4CaN_sl),      # du[12]
        D(CaMCa4CaN_sl) ~ 1e-3 * (rcn0CaN_sl - rcn02CaN_sl),                        # du[13]
        D(Ca2CaMCa4CaN_sl) ~ 1e-3 * (rcn2CaN_sl + rcn02CaN_sl - rcn24CaN_sl),              # du[14]
        D(Ca4CaMCa4CaN_sl) ~ 1e-3 * (rcn4CaN_sl + rcn24CaN_sl)                       # du[15]
    ]
    # CaMCyt equations
    CaMcyt_eqs = [
        D(CaM_cyt) ~ (1e-3 * (-rcn02_cyt - rcn0B_cyt - rcn0CaN_cyt) + J_cam_SLmyo / Vmyo),                                                                    # du[1]
        D(Ca2CaM_cyt) ~ (1e-3 * (rcn02_cyt - rcn24_cyt - rcn2B_cyt - rcn2CaN_cyt + CaMKIItot .* (-rcnCKib2_cyt + rcnCKt2a_cyt)) + J_ca2cam_SLmyo / Vmyo),       # du[2]
        D(Ca4CaM_cyt) ~ (1e-3 * (rcn24_cyt - rcn4B_cyt - rcn4CaN_cyt + CaMKIItot .* (-rcnCKib_cyt + rcnCKta_cyt)) + J_ca4cam_SLmyo / Vmyo),                       # du[3]
        D(CaMB_cyt) ~ 1e-3 * (rcn0B_cyt - rcn02B_cyt),                      # du[4]
        D(Ca2CaMB_cyt) ~ 1e-3 * (rcn02B_cyt + rcn2B_cyt - rcn24B_cyt),    # du[5]
        D(Ca4CaMB_cyt) ~ 1e-3 * (rcn24B_cyt + rcn4B_cyt),                 # du[6]
        # CaMKII equations
        D(Pb2_cyt) ~ 1e-3 * (rcnCKib2_cyt - rcnCKb2b_cyt + rcnCKt2b2_cyt),     # du[7]
        D(Pb_cyt) ~ 1e-3 * (rcnCKib_cyt + rcnCKb2b_cyt - rcnCKbt_cyt),         # du[8]
        D(Pt_cyt) ~ 1e-3 * (rcnCKbt_cyt - rcnCKta_cyt - rcnCKtt2_cyt),             # du[9]
        D(Pt2_cyt) ~ 1e-3 * (rcnCKtt2_cyt - rcnCKt2a_cyt - rcnCKt2b2_cyt),         # du[10]
        D(Pa_cyt) ~ 1e-3 * (rcnCKta_cyt + rcnCKt2a_cyt - rcnCKai_cyt),             # du[11]
        # CaN equations
        D(Ca4CaN_cyt) ~ 1e-3 * (rcnCa4CaN_cyt - rcn0CaN_cyt - rcn2CaN_cyt - rcn4CaN_cyt),      # du[12]
        D(CaMCa4CaN_cyt) ~ 1e-3 * (rcn0CaN_cyt - rcn02CaN_cyt),                        # du[13]
        D(Ca2CaMCa4CaN_cyt) ~ 1e-3 * (rcn2CaN_cyt + rcn02CaN_cyt - rcn24CaN_cyt),              # du[14]
        D(Ca4CaMCa4CaN_cyt) ~ 1e-3 * (rcn4CaN_cyt + rcn24CaN_cyt)                       # du[15]
    ]
    ## For adjusting Ca buffering in EC coupling model
    #JCaDyad = 1e-3*(2*CaMKIItot_dyad*(rcnCKtt2_dyad-rcnCKb2b_dyad) - 2*(rcn02_dyad+rcn24_dyad+rcn02B_dyad+rcn24B_dyad+rcnCa4CaN_dyad+rcn02CaN_dyad+rcn24CaN_dyad))   # [uM/msec]
    #JCaSL = 1e-3*(2*CaMKIItot*(rcnCKtt2_sl-rcnCKb2b_sl) - 2*(rcn02_sl+rcn24_sl+rcn02B_sl+rcn24B_sl+rcnCa4CaN_sl+rcn02CaN_sl+rcn24CaN_sl))   # [uM/msec]
    #JCaCyt = 1e-3*(2*CaMKIItot*(rcnCKtt2_cyt-rcnCKb2b_cyt) - 2*(rcn02_cyt+rcn24_cyt+rcn02B_cyt+rcn24B_cyt+rcnCa4CaN_cyt+rcn02CaN_cyt+rcn24CaN_cyt))   # [uM/msec]


    ## CaMKII
    ## RATE CONSTANTS and KM VALUES
    # L-Type Ca Channel (LTCC) parameters
    k_ckLCC = 0.4                   # [s^-1]
    k_pp1LCC = 0.1103               # [s^-1]
    k_pkaLCC = 13.5                 # [s^-1]
    k_pp2aLCC = 10.1                # [s^-1]
    KmCK_LCC = 12                   # [uM]
    KmPKA_LCC = 21                  # [uM]
    KmPP2A_LCC = 47                 # [uM]
    KmPP1_LCC = 9                   # [uM]

    # Ryanodine Receptor (RyR) parameters
    k_ckRyR = 0.4                   # [s^-1]
    k_pkaRyR = 1.35                 # [s^-1]
    k_pp1RyR = 1.07                 # [s^-1]
    k_pp2aRyR = 0.481               # [s^-1]

    # Basal RyR phosphorylation (numbers based on param estimation)
    kb_2809 = 0.51                  # [uM/s] - PKA site
    kb_2815 = 0.35                  # [uM/s] - CaMKII site

    KmCK_RyR = 12                   # [uM]
    KmPKA_RyR = 21                  # [uM]
    KmPP1_RyR = 9                   # [uM]
    KmPP2A_RyR = 47                 # [uM]

    # Phospholamban (PLB) parameters
    k_ckPLB = 8e-3                  # [s^-1]
    k_pp1PLB = 0.0428               # [s^-1]

    KmCK_PLB = 12
    KmPP1_PLB = 9

    # Okadaic Acid inhibition params (based on Huke/Bers [2008])
    # Want to treat OA as non-competitive inhibitor of PP1 and PP2A
    Ki_OA_PP1 = 0.78                # [uM] - Values from fit
    Ki_OA_PP2A = 0.037              # [uM] - Values from fit

    # Default PKA level
    PKAc = 95.6 * 0.54

    ## Parameters for CaMKII module
    LCCtotDyad = 31.4 * 0.9      # [uM] - Total Dyadic [LCC] - (umol/l dyad)
    LCCtotSL = 0.0846          # [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
    RyRtot = 382.6             # [uM] - Total RyR (in Dyad)
    PP1_dyad = 95.7            # [uM] - Total dyadic [PP1]
    PP1_SL = 0.57              # [uM] - Total Subsarcolemmal [PP1]
    PP2A_dyad = 95.76          # [uM] - Total dyadic PP2A
    OA = 0                     # [uM] - PP1/PP2A inhibitor Okadaic Acid
    plb_val = 106 # MOUSE
    PLBtot = plb_val           # [uM] - Total [PLB] in cytosolic units

    ## OA inhibition term (non-competitive) for PP1 and PP2A
    OA_PP1 = 1 / (1 + (OA / Ki_OA_PP1)^3)
    OA_PP2A = 1 / (1 + (OA / Ki_OA_PP2A)^3)

    CaMKIItotDyad = 120             # [uM]
    CaMKIItotSL = 120 * 8.293e-4      # [uM]
    PP1_PLBtot = 0.89               # [uM] - [umol/L cytosol]

    ## ODE EQUATIONS
    # LTCC states (note: PP2A is acting on PKA site and PP1 on CKII site)

    # Variables related to camdyad_ODEfile
    #CaMKIIact_Dyad = CaMKIItotDyad .* (Pb_dyad + Pt_dyad + Pt2_dyad + Pa_dyad)
    CaMKIIact_SL = CaMKIItotSL .* (Pb_sl + Pt_sl + Pt2_sl + Pa_sl)
    PP1_PLB_avail = 1 - I1p_PP1 / PP1_PLBtot + 0.081698
    # CaMKII phosphorylation of Dyadic LCCs
    #LCC_CKdyadn = LCCtotDyad - LCC_CKdyadp
    #LCCDyad_PHOS = (k_ckLCC*CaMKIIact_Dyad*LCC_CKdyadn)/(KmCK_LCC+LCC_CKdyadn)
    #LCCDyad_DEPHOS = (k_pp1LCC*PP1_dyad*LCC_CKdyadp)/(KmPP1_LCC+LCC_CKdyadp)*OA_PP1
    LCC_CKsln = LCCtotSL - LCC_CKslp
    LCCSL_PHOS = (k_ckLCC * CaMKIIact_SL * LCC_CKsln) / (KmCK_LCC + LCC_CKsln)
    LCCSL_DEPHOS = (k_pp1LCC * PP1_SL * LCC_CKslp) / (KmPP1_LCC + LCC_CKslp) * OA_PP1
    LCC_PKAn = LCCtotDyad - LCC_PKAp
    RyR2815n = RyRtot - RyR2815p
    RyR_BASAL = kb_2815 * RyR2815n
    #RyR_PHOS = (k_ckRyR*CaMKIIact_Dyad*RyR2815n)/(KmCK_RyR+RyR2815n)
    RyR_PP1_DEPHOS = (k_pp1RyR * PP1_dyad * RyR2815p) / (KmPP1_RyR + RyR2815p) * OA_PP1
    RyR_PP2A_DEPHOS = (k_pp2aRyR * PP2A_dyad * RyR2815p) / (KmPP2A_RyR + RyR2815p) * OA_PP2A
    RyR2809n = RyRtot - RyR2809p
    PP1_PLB = PP1_dyad * PP1_PLB_avail  # Inhibitor-1 regulation of PP1_dyad included here
    PLBT17n = PLBtot - PLBT17p
    #PLB_PHOS = (k_ckPLB*PLBT17n*CaMKIIact_Dyad)/(KmCK_PLB+PLBT17n)
    PLB_DEPHOS = (k_pp1PLB * PP1_PLB * PLBT17p) / (KmPP1_PLB + PLBT17p) * OA_PP1

    CaMKII_eqs = [
        #D(LCC_CKdyadp) ~ (LCCDyad_PHOS - LCCDyad_DEPHOS)*1e-3, # du[2]
        # CaMKII phosphorylation of Sub-sarcolemmal LCCs
        D(LCC_CKslp) ~ (LCCSL_PHOS - LCCSL_DEPHOS) * 1e-3, # du[6]
        # PKA phosphorylation (currently unused elsewhere)
        D(LCC_PKAp) ~ ((k_pkaLCC * PKAc * LCC_PKAn) / (KmPKA_LCC + LCC_PKAn) - (k_pp2aLCC * PP2A_dyad * LCC_PKAp) / (KmPP2A_LCC + LCC_PKAp) * OA_PP2A) * 1e-3, # du[1]
        # RyR states
        D(RyR2815p) ~ (RyR_BASAL - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS) * 1e-3, # du[4] + RyR_PHOS
        # PKA phosphorylation of Ser 2809 on RyR (currently unused elsewhere)
        D(RyR2809p) ~ (kb_2809 * RyR2809n + (k_pkaRyR * PKAc * RyR2809n) / (KmPKA_RyR + RyR2809n) - (k_pp1RyR * PP1_dyad * RyR2809p) / (KmPP1_RyR + RyR2809p) * OA_PP1) * 1e-3, # du[3]
        # PLB states
        D(PLBT17p) ~ (-PLB_DEPHOS) * 1e-3 # du[5]PLB_PHOS
    ]


    ## BAR
    # Drug concentrations
    Ligtot = 0.0               # [uM] - SET LIGAND CONCENTRATION (0 or 0.1)
    FSK = 0
    IBMX = 0
    LCCtotBA = 0.025           # [uM] - [umol/L cytosol]
    plb_val = 106                # MOUSE
    PP1_PLBtot = 0.89          # [uM] - [umol/L cytosol]
    PLMtotBA = 48              # [uM] - [umol/L cytosol] MOUSE
    PLBtotBA = plb_val                     # [uM] - [umol/L cytosol]
    ISO = Ligtot

    ## b-AR module
    b1ARtot = 0.00528        # (uM) total b1-AR protein # MOUSE
    #b1ARtot=0.028  # RABBIT
    kf_LR = 1              # (1/[uM ms]) forward rate for ISO binding to b1AR
    kr_LR = 0.285          # (1/ms) reverse rate for ISO binding to b1AR
    kf_LRG = 1              # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
    kr_LRG = 0.062          # (1/ms) reverse rate for ISO:b1AR association with Gs
    kf_RG = 1              # (1/[uM ms]) forward rate for b1AR association with Gs
    kr_RG = 33             # (1/ms) reverse rate for b1AR association with Gs
    Gstot = 3.83           # (uM) total Gs protein
    k_G_act = 16e-3          # (1/ms) rate constant for Gs activation
    k_G_hyd = 0.8e-3         # (1/ms) rate constant for G-protein hydrolysis
    k_G_reassoc = 1.21       # (1/[uM ms]) rate constant for G-protein reassociation
    kf_bARK = 1.1e-6         # (1/[ms]) forward rate for b1AR phosphorylation by b1ARK
    kr_bARK = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
    kf_PKA = 3.6e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
    kr_PKA = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by PKA
    b1ARact = b1ARtot - b1AR_S464 - b1AR_S301
    b1AR = b1ARact - LR - LRG - RG
    Gs = Gstot - LRG - RG - Gsby
    bARK_desens = kf_bARK * (LR + LRG)
    bARK_resens = kr_bARK * b1AR_S464
    PKA_desens = kf_PKA * PKACI * b1ARact
    PKA_resens = kr_PKA * b1AR_S301
    G_act = k_G_act * (RG + LRG)
    G_hyd = k_G_hyd * GsaGTPtot
    G_reassoc = k_G_reassoc * GsaGDP * Gsby

    bar_eqs = [
        D(LR) ~ (kf_LR * ISO * b1AR - kr_LR * LR + kr_LRG * LRG - kf_LRG * LR * Gs),    # du[1]
        D(LRG) ~ (kf_LRG * LR * Gs - kr_LRG * LRG - k_G_act * LRG),                 # du[2]
        D(RG) ~ (kf_RG * b1AR * Gs - kr_RG * RG - k_G_act * RG),                    # du[3]
        D(b1AR_S464) ~ (bARK_desens - bARK_resens),                         # du[4]
        D(b1AR_S301) ~ (PKA_desens - PKA_resens),                           # du[5]
        D(GsaGTPtot) ~ (G_act - G_hyd),                                     # du[6]
        D(GsaGDP) ~ (G_hyd - G_reassoc),                                    # du[7]
        D(Gsby) ~ (G_act - G_reassoc)                                       # du[8]
    ]

    ## cAMP module
    ACtot = 70.57e-3        # (uM) total adenylyl cyclase # MOUSE
    # ACtot=47e-3  # RABBIT
    ATP = 5e3             # (uM) total ATP
    k_AC_basal = 0.2e-3          # (1/ms) basal cAMP generation rate by AC
    Km_AC_basal = 1.03e3          # (uM) basal AC affinity for ATP
    Kd_AC_Gsa = 0.4             # (uM) Kd for AC association with Gsa
    kf_AC_Gsa = 1               # (1/[uM ms]) forward rate for AC association with Gsa
    kr_AC_Gsa = Kd_AC_Gsa       # (1/ms) reverse rate for AC association with Gsa
    k_AC_Gsa = 8.5e-3          # (1/ms) basal cAMP generation rate by AC:Gsa
    Km_AC_Gsa = 315.0           # (uM) AC:Gsa affinity for ATP
    Kd_AC_FSK = 44.0            # (uM) Kd for FSK binding to AC
    k_AC_FSK = 7.3e-3          # (1/ms) basal cAMP generation rate by AC:FSK
    Km_AC_FSK = 860.0           # (uM) AC:FSK affinity for ATP
    PDEtot = 22.85e-3        # (uM) total phosphodiesterase
    k_cAMP_PDE = 5e-3            # (1/ms) cAMP hydrolysis rate by PDE
    k_cAMP_PDEp = 2 * k_cAMP_PDE    # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
    Km_PDE_cAMP = 1.3             # (uM) PDE affinity for cAMP
    Kd_PDE_IBMX = 30.0            # (uM) Kd_R2cAMP_C for IBMX binding to PDE
    k_PKA_PDE = 7.5e-3          # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
    k_PP_PDE = 1.5e-3          # (1/ms) rate constant for PDE dephosphorylation by phosphatases
    cAMP = cAMPtot - (RCcAMP_I + 2 * RCcAMPcAMP_I + 2 * RcAMPcAMP_I) - (RCcAMP_II + 2 * RCcAMPcAMP_II + 2 * RcAMPcAMP_II)
    AC = ACtot - AC_GsaGTP
    GsaGTP = GsaGTPtot - AC_GsaGTP
    AC_FSK = FSK * AC / Kd_AC_FSK
    AC_ACT_BASAL = k_AC_basal * AC * ATP / (Km_AC_basal + ATP)
    AC_ACT_GSA = k_AC_Gsa * AC_GsaGTP * ATP / (Km_AC_Gsa + ATP)
    AC_ACT_FSK = k_AC_FSK * AC_FSK * ATP / (Km_AC_FSK + ATP)
    PDE_IBMX = PDEtot * IBMX / Kd_PDE_IBMX
    PDE = PDEtot - PDE_IBMX - PDEp
    PDE_ACT = k_cAMP_PDE * PDE * cAMP / (Km_PDE_cAMP + cAMP) + k_cAMP_PDEp * PDEp * cAMP / (Km_PDE_cAMP + cAMP)

    cAMP_eqs = [
        D(AC_GsaGTP) ~ (kf_AC_Gsa * GsaGTP * AC - kr_AC_Gsa * AC_GsaGTP),
        D(PDEp) ~ (k_PKA_PDE * PKACII * PDE - k_PP_PDE * PDEp),
        D(cAMPtot) ~ (AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT)
    ]


    ## PKA module
    PKItot = 0.18             # (uM) total PKI
    kf_RC_cAMP = 1                # (1/[uM ms]) Kd for PKA RC binding to cAMP
    kf_RCcAMP_cAMP = 1                # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
    kf_RcAMPcAMP_C = 4.375            # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
    kf_PKA_PKI = 1                # (1/[uM ms]) Ki for PKA inhibition by PKI
    kr_RC_cAMP = 1.64             # (1/ms) Kd for PKA RC binding to cAMP
    kr_RCcAMP_cAMP = 9.14             # (1/ms) Kd for PKA RC:cAMP binding to cAMP
    kr_RcAMPcAMP_C = 1                # (1/ms) Kd for PKA R:cAMPcAMP binding to C
    kr_PKA_PKI = 2e-4             # (1/ms) Ki for PKA inhibition by PKI
    epsilon = 10               # (-) AKAP-mediated scaling factor
    PKI = PKItot - PKACI_PKI - PKACII_PKI

    PKA_eqs = [
        D(RC_I) ~ (-kf_RC_cAMP * RC_I * cAMP + kr_RC_cAMP * RCcAMP_I),
        D(RCcAMP_I) ~ (-kr_RC_cAMP * RCcAMP_I + kf_RC_cAMP * RC_I * cAMP - kf_RCcAMP_cAMP * RCcAMP_I * cAMP + kr_RCcAMP_cAMP * RCcAMPcAMP_I),
        D(RCcAMPcAMP_I) ~ (-kr_RCcAMP_cAMP * RCcAMPcAMP_I + kf_RCcAMP_cAMP * RCcAMP_I * cAMP - kf_RcAMPcAMP_C * RCcAMPcAMP_I + kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI),
        D(RcAMPcAMP_I) ~ (-kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI + kf_RcAMPcAMP_C * RCcAMPcAMP_I),
        D(PKACI) ~ (-kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI + kf_RcAMPcAMP_C * RCcAMPcAMP_I - kf_PKA_PKI * PKACI * PKI + kr_PKA_PKI * PKACI_PKI),
        D(PKACI_PKI) ~ (-kr_PKA_PKI * PKACI_PKI + kf_PKA_PKI * PKACI * PKI),
        D(RC_II) ~ (-kf_RC_cAMP * RC_II * cAMP + kr_RC_cAMP * RCcAMP_II),
        D(RCcAMP_II) ~ (-kr_RC_cAMP * RCcAMP_II + kf_RC_cAMP * RC_II * cAMP - kf_RCcAMP_cAMP * RCcAMP_II * cAMP + kr_RCcAMP_cAMP * RCcAMPcAMP_II),
        D(RCcAMPcAMP_II) ~ (-kr_RCcAMP_cAMP * RCcAMPcAMP_II + kf_RCcAMP_cAMP * RCcAMP_II * cAMP - kf_RcAMPcAMP_C * RCcAMPcAMP_II + kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII),
        D(RcAMPcAMP_II) ~ (-kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII + kf_RcAMPcAMP_C * RCcAMPcAMP_II),
        D(PKACII) ~ (-kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII + kf_RcAMPcAMP_C * RCcAMPcAMP_II - kf_PKA_PKI * PKACII * PKI + kr_PKA_PKI * PKACII_PKI),
        D(PKACII_PKI) ~ (-kr_PKA_PKI * PKACII_PKI + kf_PKA_PKI * PKACII * PKI)
    ]


    ## I-1/PP1 module
    I1tot = 0.3             # (uM) total inhibitor 1
    k_PKA_I1 = 60e-3           # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
    Km_PKA_I1 = 1.0             # (uM) Km for I-1 phosphorylation by type 1 PKA
    Vmax_PP2A_I1 = 14.0e-3         # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
    Km_PP2A_I1 = 1.0             # (uM) Km for I-1 dephosphorylation by PP2A
    Ki_PP1_I1 = 1.0e-3          # (uM) Ki for PP1 inhibition by I-1
    kf_PP1_I1 = 1               # (uM) Ki for PP1 inhibition by I-1
    PP1tot = PP1_PLBtot      # PP1tot = 0.89  # (uM) total phosphatase 1
    kr_PP1_I1 = Ki_PP1_I1       # (uM) Ki for PP1 inhibition by I-1
    I1 = I1tot - I1ptot
    PP1 = PP1tot - I1p_PP1
    I1p = I1ptot - I1p_PP1
    I1_phosph = k_PKA_I1 * PKACI * I1 / (Km_PKA_I1 + I1)
    I1_dephosph = Vmax_PP2A_I1 * I1ptot / (Km_PP2A_I1 + I1ptot)

    PP1_eqs = [
        D(I1p_PP1) ~ (kf_PP1_I1 * PP1 * I1p - kr_PP1_I1 * I1p_PP1),
        D(I1ptot) ~ (I1_phosph - I1_dephosph)
    ]


    ## PLB module
    PLBtot = PLBtotBA   # [uM]
    k_PKA_PLB = 54e-3   # [1/ms]
    Km_PKA_PLB = 21     # [uM]
    k_PP1_PLB = 8.5e-3  # [1/ms]
    Km_PP1_PLB = 7.0    # [uM]

    PLB = PLBtot - PLBp
    PLB_phosph = k_PKA_PLB * PKACI * PLB / (Km_PKA_PLB + PLB)
    PLB_dephosph = k_PP1_PLB * PP1 * PLBp / (Km_PP1_PLB + PLBp)

    PLB_eqs = [
        D(PLBp) ~ (PLB_phosph - PLB_dephosph)
    ]


    ## PLM module (included 09/18/12) MOUSE
    PLMtot = PLMtotBA   # [uM]
    k_PKA_PLM = 54e-3   # [1/ms]
    Km_PKA_PLM = 21     # [uM]
    k_PP1_PLM = 8.5e-3  # [1/ms]
    Km_PP1_PLM = 7.0    # [uM]
    PLM = PLMtot - PLMp
    PLM_phosph = k_PKA_PLM * PKACI * PLM / (Km_PKA_PLM + PLM)
    PLM_dephosph = k_PP1_PLM * PP1 * PLMp / (Km_PP1_PLM + PLMp)

    PLM_eqs = [
        D(PLMp) ~ (PLM_phosph - PLM_dephosph)
    ]

    ## LCC module
    PKAIItot = 0.059        # (uM) total type 2 PKA # MOUSE
    LCCtot = LCCtotBA       # [uM]
    PKACII_LCCtot = 0.025   # [uM]
    PP1_LCC = 0.025         # [uM]
    PP2A_LCC = 0.025        # [uM]
    k_PKA_LCC = 54e-3       # [1/ms]
    Km_PKA_LCC = 21         # [uM]
    k_PP1_LCC = 8.52e-3     # [1/ms] RABBIT, MOUSE
    Km_PP1_LCC = 3          # [uM]
    k_PP2A_LCC = 10.1e-3    # [1/ms]
    Km_PP2A_LCC = 3         # [uM]
    PKACII_LCC = (PKACII_LCCtot / PKAIItot) * PKACII
    LCCa = LCCtot - LCCap
    LCCa_phosph = epsilon * k_PKA_LCC * PKACII_LCC * LCCa / (Km_PKA_LCC + epsilon * LCCa)
    LCCa_dephosph = epsilon * k_PP2A_LCC * PP2A_LCC * LCCap / (Km_PP2A_LCC + epsilon * LCCap)
    LCCb = LCCtot - LCCbp
    LCCb_phosph = epsilon * k_PKA_LCC * PKACII_LCC * LCCb / (Km_PKA_LCC + epsilon * LCCb)
    LCCb_dephosph = epsilon * k_PP1_LCC * PP1_LCC * LCCbp / (Km_PP1_LCC + epsilon * LCCbp)

    LCC_eqs = [
        D(LCCap) ~ (LCCa_phosph - LCCa_dephosph),
        D(LCCbp) ~ (LCCb_phosph - LCCb_dephosph)
    ]

    ## RyR module (not included in Yang-Saucerman)
    PKAIIryrtot = 0.034         # [uM]
    PP1ryr = 0.034              # [uM]
    PP2Aryr = 0.034             # [uM]
    kcat_pka_ryr = 54e-3        # [1/ms]
    Km_pka_ryr = 21             # [uM]
    kcat_pp1_ryr = 8.52e-3      # [1/ms]
    Km_pp1_ryr = 7              # [uM]
    kcat_pp2a_ryr = 10.1e-3     # [1/ms]
    Km_pp2a_ryr = 4.1           # [uM]
    RyRtot = 0.135              # [uM]

    PKACryr = (PKAIIryrtot / PKAIItot) * PKACII
    RyR = RyRtot - RyRp
    RyRPHOSPH = epsilon * kcat_pka_ryr * PKACryr * RyR / (Km_pka_ryr + epsilon * RyR)
    RyRDEPHOSPH1 = epsilon * kcat_pp1_ryr * PP1ryr * RyRp / (Km_pp1_ryr + epsilon * RyRp)
    RyRDEPHOSPH2A = epsilon * kcat_pp2a_ryr * PP2Aryr * RyRp / (Km_pp2a_ryr + epsilon * RyRp)

    RyRp_eqs = [
        D(RyRp) ~ (RyRPHOSPH - RyRDEPHOSPH1 - RyRDEPHOSPH2A)
    ]


    ## TnI module
    PP2A_TnI = 0.67         # [uM]
    k_PKA_TnI = 54e-3       # [1/ms]
    Km_PKA_TnI = 21         # [uM]
    k_PP2A_TnI = 10.1e-3    # [1/ms]
    Km_PP2A_TnI = 4.1       # [uM]
    TnItot = 70             # [uM]
    TnI = TnItot - TnIp
    TnI_phosph = k_PKA_TnI * PKACI * TnI / (Km_PKA_TnI + TnI)
    TnI_dephosph = k_PP2A_TnI * PP2A_TnI * TnIp / (Km_PP2A_TnI + TnIp)

    TnI_eqs = [
        D(TnIp) ~ (TnI_phosph - TnI_dephosph)
    ]


    ## Iks module (not present in mouse)
    IKstot = 0.025

    IKs_eqs = [
        D(KS79) ~ 0,  # ydot(27) not ODE
        D(KS80) ~ 0,  # ydot(28) not ODE
        D(KSp) ~ 0  # ydot(29)
    ]


    ## CFTR module (included 04/30/10)
    ICFTRtot = 0.025

    CFTR_eqs = [
        D(CFTRp) ~ 0  #CFTRphos - CFTRdephos  # ydot(30)
    ]

    ## Ikur module (included 04/10/12) MOUSE
    PKAII_KURtot = 0.025    # [uM]
    PP1_KURtot = 0.025      # [uM]
    k_pka_KUR = 54e-3       # [1/ms]
    Km_pka_KUR = 21         # [uM]
    k_pp1_KUR = 8.52e-3     # [1/ms]
    Km_pp1_KUR = 7          # [uM]
    IKurtot = 0.025         # [uM]
    KURn = IKurtot - KURp   # Non-phos = tot - phos
    PKAC_KUR = (PKAII_KURtot / PKAIItot) * PKACII     # (PKA_KURtot/PKAIItot)*PKAIIact
    KURphos = epsilon * KURn * PKAC_KUR * k_pka_KUR / (Km_pka_KUR + epsilon * KURn)
    KURdephos = PP1_KURtot * k_pp1_KUR * epsilon * KURp / (Km_pp1_KUR + epsilon * KURp)
    Ikur_eqs = [
        D(KURp) ~ (KURphos - KURdephos)
    ]


    ## ECC

    ## Adjusting Variables
    LCCtotDyad = 31.4 * 0.9        # [uM] - Total Dyadic [LCC] - (umol/l dyad)
    RyRtot = 382.6              # [uM] - Total RyR (in Dyad)
    plb_val = 106                 # MOUSE
    LCCtotBA = 0.025            # [uM] - [umol/L cytosol]
    TnItotBA = 70               # [uM] - [umol/L cytosol]
    IKurtotBA = 0.025           # [uM] - [umol/L cytosol] MOUSE
    PLBtot = plb_val                        # [uM] - Total [PLB] in cytosolic units
    PLBtotBA = plb_val                      # [uM] - [umol/L cytosol]


    RyR_CKp = RyR2815p / RyRtot
    PLB_CKp = PLBT17p / PLBtot
    LCCa_PKAp = LCCap / LCCtotBA
    LCCb_PKAp = LCCbp / LCCtotBA
    PLB_PKAn = (PLBtotBA - PLBp) / PLBtotBA
    TnI_PKAp = TnIp / TnItotBA
    IKur_PKAp = KURp / IKurtotBA

    # -------------------------------------------------------------------------
    # Model of ECC of rat neonatal ventricular myocyte 2009
    # Code & model: Topi Korhonen, University of Oulu (topi.korhonen@oulu.fi)
    #
    # PLEASE MENTION THE FOLLOWING REFERENCE WHEN USING THIS CODE OR PART OF IT:
    # Korhonen et al. "Model of excitation-contraction coupling of rat neonatal
    # ventricular myocytes" Biophys J. 2009, Feb; 96(3):1189-1209
    #
    # ONLY FOR ACADEMIC USE, DO NOT DISTRIBUTE
    # -------------------------------------------------------------------------

    # Index numbers, put the Cai equations to the end -> increase i_Cai_sub_SR when
    # adding other odes

    # Physical constants
    F = 96.5
    T = 305
    R = 8.314
    Cm = 1.0

    # Ion concentrations in DMEM
    Cao = 1796
    Nao = 154578
    Ko = 5366

    # Cell geometry
    rSR_true = 6
    rSL_true = 10.5

    # Ca diffusion grid
    dx = 0.1
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial index of Cai diffusion

    # More cell geometry
    V_sub_SR = 4 / 3 * pi * (rSR_true + dx)^3 / 1000 - 4 / 3 * pi * (rSR_true)^3 / 1000 # pl
    V_sub_SL = 4 / 3 * pi * rSL_true^3 / 1000 - 4 / 3 * pi * (rSL_true - dx)^3 / 1000 #pl
    Acap = 4 * pi * rSL_true^2 * 1e-8 # cm^2
    VSR = 0.043 * 1.5 * 1.4
    VNSR = 0.9 * VSR
    VJSR = VSR - VNSR
    Vmyo = 4 / 3 * pi * rSL_true^3 / 1000 - 4 / 3 * pi * rSR_true^3 / 1000


    # Ca buffers
    csqntot = 24750
    Kmcsqn = 800
    betaSR = inv(1 + csqntot * Kmcsqn ./ (i_CaJSR + Kmcsqn) .^ 2)


    # NCX from Pandit rat model
    fNaCa = 1
    kNaCa = 2.2680e-016
    dNaCa = 1e-16
    gamma = 0.5

    # INaK
    INaKmax = 2.7
    KmNai = 18600
    nNaK = 3.2
    KmKo = 1500


    # PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph)
    ICa_scale = 0.95#5.25
    fracLCCbp0 = 0.250657 # Derived quantity - (LCCbp(baseline)/LCCbtot)
    fracLCCbpISO = 0.525870 # Derived quantity - (LCCbp(ISO)/LCCbtot)
    a_favail = (1.56 - 1) / (fracLCCbpISO / fracLCCbp0 - 1) # fracLCCbp ISO (x1.56 o.1 ISO)
    favail = (1 - a_favail) + a_favail * (LCCb_PKAp / fracLCCbp0)  # Test (max x2.52 100# phosph)
    ICa_scalep = ICa_scale * favail


    # --------------------------------------------------------
    # Equations used in odes
    # --------------------------------------------------------

    ENa = R * T / F * nm.log((Nao) / (Nai))
    EK = R * T / F * nm.log((Ko) / (Ki))
    ECa = R * T / 2 / F * nm.log(Cao / Cai_sub_SL)

    # NCX
    INaCa = ICa_scalep * kNaCa * ((exp(0.03743 * gamma .* V) .* Nai .^ 3 .* Cao - exp(0.03743 * (gamma - 1) .* V) .* Nao^3 .* Cai_sub_SL .* fNaCa) / (1 + dNaCa * (Nao^3 .* Cai_sub_SL .* fNaCa + Nai .^ 3 .* Cao)))

    # L-type calcium current

    GCaL = 1.3125e-4 * 0.8 * 0.6
    ICaL = ICa_scalep * GCaL * i_d * i_f * i_fca * 4 * V * F^2 / R / T * (Cai_sub_SL * exp(2 * V * F / R / T) - 0.341 * Cao) / (exp(2 * V * F / R / T) - 1)
    dinf = 1 / (1 + exp((11.1 + V) / -7.2))
    alphad = 1.4 / (1 + exp((-35 - V) / 13)) + 0.25
    betad = 1.4 / (1 + exp((V + 5) / 5))
    gammad = 1 / (1 + exp((50 - V) / 20))
    taud = alphad * betad + gammad
    finf = 1 / (1 + exp((V + 23.3) / 5.4))
    tauf = 1125 * exp(-(V + 27)^2 / 240) + 165 / (1 + exp((25 - V) / 10)) + 120
    #a = finf/tauf * (1-junc_mode2)
    #b = 1/8 * (1-finf) / tauf * junc_mode2
    alphafca = 1 / (1 + (Cai_sub_SL * 1e-3 / (0.000325 * 1.5))^8)
    betafca = 0.1 / (1 + exp((Cai_sub_SL * 1e-3 - 0.0005) / 0.0001))
    gammafca = 0.2 / (1 + exp((Cai_sub_SL * 1e-3 - 0.00075) / 0.0008))
    fcainf = (alphafca + betafca + gammafca + 0.23) / 1.46
    taufca = 10 # modif
    kfca = 1 - (fcainf > i_fca) * (V > -60)

    # T-Type
    gCaT = 0.2
    binf = 1 / (1 + exp(-(V + 37.49098) / 5.40634))
    taub = 0.6 + 5.4 / (1 + exp((V + 100) * 0.03))
    ginf = 1 / (1 + exp((V + 66) / 6))
    taug = 1 + 40 / (1 + exp((V + 65) * 0.08))
    ICaT = gCaT * i_b * i_g * (V - ECa + 106.5)

    # Cab
    gCab = 0.0008
    ICab = gCab * (V - ECa)

    # Nab
    gNab = 0.0026
    INab = gNab * (V - ENa)

    # If
    gf = 0.021
    fNa = 0.2
    fK = 1 - fNa
    yinf = 1 / (1 + exp((V + 78.65) ./ 6.33)) # Fitted
    tauy = 1 / (0.11885 .* exp((V + 75) ./ 28.37) + 0.56236 .* exp((V + 75) ./ -14.19)) .* 1000
    IfNa = gf * i_y * fNa * (V - ENa)
    IfK = gf * i_y * fK * (V - EK)
    If = gf * i_y * (fNa * (V - ENa) + fK * (V - EK))

    # IK1
    IK1 = 0.0515 .* (Ko / (Ko + 210)) .* ((V - EK - 6.1373) ./ (0.1653 + exp(0.0319 * (V - EK - 6.1373))))

    # Ito
    gt = 0.1
    sinf = 1 ./ (1 + exp((V + 31.97156) ./ 4.64291))
    rinf = 1 ./ (1 + exp((V - 3.55716) ./ -14.61299))
    slowinf = sinf
    taur = 1 ./ (45.16 .* exp(0.03577 .* (V + 50)) + 98.9 .* exp(-0.1 .* (V + 38))) * 1000
    taus = (0.35 .* exp(-(((V + 70) / 15) .^ 2)) + 0.035) * 1000 - 26.9
    tausslow = (3.7 .* exp(-(((V + 70) / 30) .^ 2)) + 0.035) * 1000 + 37.4
    Ito = gt * i_r .* (0.706 .* i_s + 0.294 .* i_sslow) * (V - EK)

    # INa
    gNa = 35
    Naminf = 1 / (1 + exp((V + 45) / -6.5))
    Nahinf = 1 / (1 + exp((V + 76.1) / 6.07))
    Najinf = Nahinf
    Nataum = 0.00136 / (0.32 * (V + 47.13) / (1 - exp(-0.1 * (V + 47.13))) + 0.08 * exp(-V / 11))
    Natauh = ifelse(V >= -40, 0.0004537 .* (1 + exp((V + 10.66) ./ -11.1)), 0.00349 ./ (0.135 .* exp((V + 80) ./ -6.8) + 3.56 .* exp(0.079 .* V) + 3.1e5 .* exp(0.35 .* V)))
    Natauj = ifelse(V >= -40, 0.01163 .* (1 + exp(-0.1 .* (V + 32))) ./ exp(-2.535e-7 .* V), 0.00349 ./ ((V + 37.78) ./ (1 + exp(0.311 .* (V + 79.23))) .* (-127140 .* exp(0.2444 .* V) - 3.474e-5 .* exp(-0.04391 .* V)) + 0.1212 .* exp(-0.01052 .* V) ./ (1 + exp(-0.1378 .* (V + 40.14)))))

    INa = gNa * i_Nam .^ 3 * i_Nah * i_Naj * (V - ENa)

    # RyR
    n = 4
    PC1 = 1 - i_PO1
    kapos = 3 / 3
    kaneg = 0.48 / 3
    nu1 = 0.01 * 3 / 3
    KmRyR = 1.35 * 2.6 ./ (1 + exp((i_CaJSR - 530) ./ 200)) + 1.5 - 0.9 - 0.3 - 0.05
    Jrel = nu1 * (i_PO1) * (i_CaJSR - Cai_sub_SR)


    # SR Ca-ATPase
    Vmaxf = 0.9996
    Vmaxr = Vmaxf
    Kmf = 0.5
    Kmr = 7000 * Kmf
    Hf = 2
    Hr = 1 * Hf
    k = 5e-6

    # CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
    fCKII_PLB = (1 - 0.5 * PLB_CKp)  # Max effect: fCKII_PLB=0.5
    fracPKA_PLBo = 1 - 0.079755
    fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * (100 - 55.31) / 100 + 55.31 / 100
    # Select smaller value (resulting in max reduction of Kmf)
    Kmfp = ifelse(fCKII_PLB < fPKA_PLB, Kmf * fCKII_PLB, Kmf * fPKA_PLB) * 2   #fCKII_PLB

    Jup = (Vmaxf .* (Cai_sub_SR ./ Kmfp) .^ Hf - Vmaxr .* (CaNSR ./ Kmr) .^ Hr) ./ (1 + (Cai_sub_SR ./ Kmfp) .^ Hf + (CaNSR ./ Kmr) .^ Hr)

    # Jleak
    kleak = (1 / 2 + 5 * RyR_CKp / 2) * k
    Jleak = kleak * (CaNSR - Cai_sub_SR)

    # IKs
    GKs = 0.05
    # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    fracIKurp0 = 0.437635       # Derived quantity (IKur_PKAp(baseline)/IKurtot)
    fracIKurpISO = 0.718207     # Derived quantity (IKur_PKAp(ISO)/IKurtot)
    a_Kur = (1.20 - 1) / (fracIKurpISO / fracIKurp0 - 1)
    fracIKuravail = (1 - a_Kur) + a_Kur * (IKur_PKAp / fracIKurp0)  # +20# with 0.1 uM ISO
    IKs = GKs * i_nKs^2 * (V - EK) * fracIKuravail * 2
    alphan = 0.00000481333 * (V + 26.5) / (1 - exp(-0.128 * (V + 26.5)))
    betan = 0.0000953333 * exp(-0.038 * (V + 26.5))
    nKsinf = alphan / (alphan + betan)
    nKstau = 750


    # Na/K pump current
    sigma = 1 / 7 * (exp(Nao / 67300) - 1)
    fNaK = 1 / (1 + 0.1245 * exp(-0.1 * V * F / R / T) + 0.0365 * sigma * exp(-V * F / R / T))
    INaK = INaKmax * fNaK * Ko / (1 + (KmNai / Nai)^(nNaK)) / (Ko + KmKo)

    # IKr
    GKr = 0.06
    kf = 0.023761
    kb = 0.036778
    IKr = i_OK * GKr * (V - R * T / F * nm.log((0.98 * Ko + 0.02 * Nao) / (0.98 * Ki + 0.02 * Nai)))
    CK0 = 1 - (i_CK1 + i_CK2 + i_OK + i_IK)
    alphaa0 = 0.022348 * exp(0.01176 * V)
    betaa0 = 0.047002 * exp(-0.0631 * V)
    alphaa1 = 0.013733 * exp(0.038198 * V)
    betaa1 = 0.0000689 * exp(-0.04178 * V)
    alphai_mERG = 0.090821 * exp(0.023391 * V)
    betai_mERG = 0.006497 * exp(-0.03268 * V)

    Jtr = (CaNSR - i_CaJSR) / 200

    # Ca fluxes
    JCa_SL = (2 * INaCa - ICaL - ICaT - ICab) * Acap * Cm / 2 / 1 / F * 1e6
    JCa_SR = Jleak - Jup + Jrel

    # -------------------------------------------------------------------
    # Differential equations
    # -------------------------------------------------------------------
    # Cai
    m = length(j)
    @variables Cai(t)[1:m]
    # Diffusion coefficient
    @parameters Dca = 7 # mum^2/ms set to achive correct diff. speed 0.31 mum/ms

    eqs = [D(Cai[1]) ~ (Dca / (j[1] * dx^2) * ((1 + j[1]) * Cai[2] - 2 * j[1] * Cai[1] + (j[1] - 1) * Cai[1]) + JCa_SR / V_sub_SR) * beta_cai(Cai[1], TnI_PKAp)]


    #eqs = [eq, Cai_sub_SR~Cai[1], Cai_sub_SL~Cai[m], Ca_j~Cai[23]]

    for i in 2:m-1
        eq = D(Cai[i]) ~ (Dca / (j[i] * dx^2) * ((1 + j[i]) * Cai[i+1] - 2 * j[i] * Cai[i] + (j[i] - 1) * Cai[i-1])) * beta_cai(Cai[i], TnI_PKAp)
        push!(eqs, eq)
    end

    eq_end = D(Cai[m]) ~ (Dca / (j[m] * dx^2) * ((1 + j[m]) * Cai[m] - 2 * j[m] * Cai[m] + (j[m] - 1) * Cai[m-1]) + JCa_SL / V_sub_SL) * beta_cai(Cai[m], TnI_PKAp)
    push!(eqs, eq_end)

    # Istim
    freq = 1                     # [Hz] - CHANGE DEPENDING ON FREQUENCY (1, <=1, -40) (2, <=2, -35) (3, <=2, -30)
    cycleLength = 1e3 / freq      # [ms]

    Istim = ifelse(mod(t, cycleLength) <= 1, -40, 0.0)

    # Other odes
    SR_eqs = [
        D(i_CaJSR) ~ betaSR * (-Jrel + Jtr) / VJSR,
        D(CaNSR) ~ (Jup - Jleak - Jtr) / VNSR,
        D(V) ~ -(INab + INaCa + ICaL + ICaT + If + Ito + IK1 + IKs + IKr + INa + INaK + ICab + Istim) / Cm
    ]
    # Nai, Ki
    Na_K_eqs = [
        D(Nai) ~ -(IfNa + INab + INa + 3 * INaCa + 3 * INaK) * Acap * Cm / F * 1e6 / Vmyo,
        D(Ki) ~ -(IfK + Ito + IK1 + IKs + IKr + Istim - 2 * INaK) * Acap * Cm / F * 1e6 / Vmyo
    ]

    # T-Type
    Ttype_eqs = [
        D(i_b) ~ (binf - i_b) / taub,
        D(i_g) ~ (ginf - i_g) / taug
    ]

    # L-type calcium current
    Ltype_eqs = [
        D(i_d) ~ (dinf - i_d) / taud,
        D(i_f) ~ (finf - i_f) / tauf,
        #D(i_f) ~ a-(a+b)*i_f*1.5,
        D(i_fca) ~ kfca * (fcainf - i_fca) / taufca
    ]


    # If
    IF_eqs = [
        D(i_y) ~ (yinf - i_y) / tauy
    ]


    # IKto
    IKto_eqs = [
        D(i_r) ~ (rinf - i_r) / taur,
        D(i_s) ~ (sinf - i_s) / taus,
        D(i_sslow) ~ (slowinf - i_sslow) / tausslow
    ]

    # INa
    INa_eqs = [
        D(i_Nam) ~ (Naminf - i_Nam) / Nataum / 1000,
        D(i_Nah) ~ (Nahinf - i_Nah) / Natauh / 1000,
        D(i_Naj) ~ (Najinf - i_Naj) / Natauj / 1000
    ]


    # RyR
    ryr_eqs = [
        D(i_PO1) ~ kapos * (Cai_sub_SR)^n / ((Cai_sub_SR)^n + KmRyR^n) * PC1 - kaneg * i_PO1,
        D(i_PO2) ~ 0,
        D(i_PC2) ~ 0
    ]


    # IKs
    InKs_eqs = [
        D(i_nKs) ~ (nKsinf - i_nKs) / nKstau
    ]


    # Rapid delayed rectifier K current (mERG)
    IKr_eqs = [
        D(i_CK1) ~ (alphaa0 * CK0 - betaa0 * i_CK1 + kb * i_CK2 - kf * i_CK1),
        D(i_CK2) ~ (kf * i_CK1 - kb * i_CK2 + betaa1 * i_OK - alphaa1 * i_CK2),
        D(i_OK) ~ (alphaa1 * i_CK2 - betaa1 * i_OK + betai_mERG * i_IK - alphai_mERG * i_OK),
        D(i_IK) ~ (alphai_mERG * i_OK - betai_mERG * i_IK)
    ]


    return vcat(eqs, SR_eqs, Na_K_eqs, Ttype_eqs, Ltype_eqs, IF_eqs, IKto_eqs, INa_eqs, ryr_eqs, InKs_eqs, IKr_eqs, CaMSL_eqs, CaMcyt_eqs, #CaMdyad_eqs,
        CaMKII_eqs, bar_eqs, cAMP_eqs, PKA_eqs, PP1_eqs, PLB_eqs, PLM_eqs, LCC_eqs, RyRp_eqs, TnI_eqs, IKs_eqs, CFTR_eqs, Ikur_eqs)
end

eq_morotti = get_Morotti_equations()

@named osys = ODESystem(eq_morotti)

#osys = structural_simplify(osys)

#@variables t Cai_mean(t)
@variables t Cai(t)[1:45]
Cai_mean = mean(skipmissing(Cai))

##Chemical Reaction of CaMKII Activity (Including OX states)
ca_model = @reaction_network begin
    ##  Two Ca2+ ions bind to C or N-lobe.
    (k_1C_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_2C_on / (k_1C_off + k_2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_1C_off * k_2C_off / (k_1C_off + k_2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), CaM0 <--> Ca2CaM_C
    (k_1N_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_2N_on / (k_1N_off + k_2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_1N_off * k_2N_off / (k_1N_off + k_2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), CaM0 <--> Ca2CaM_N
    (k_1C_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_2C_on / (k_1C_off + k_2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_1C_off * k_2C_off / (k_1C_off + k_2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), Ca2CaM_C <--> Ca4CaM
    (k_1N_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_2N_on / (k_1N_off + k_2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_1N_off * k_2N_off / (k_1N_off + k_2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), Ca2CaM_N <--> Ca4CaM
    ##  Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII complex.
    (k_K1C_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_K2C_on / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_K1C_off * k_K2C_off / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), CaM0_CaMK <--> Ca2CaM_C_CaMK
    (k_K1N_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_K2N_on / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_K1N_off * k_K2N_off / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), CaM0_CaMK <--> Ca2CaM_N_CaMK
    (k_K1C_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_K2C_on / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_K1C_off * k_K2C_off / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), Ca2CaM_C_CaMK <--> Ca4CaM_CaMK
    (k_K1N_on * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart) * k_K2N_on / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)), k_K1N_off * k_K2N_off / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), Ca2CaM_N_CaMK <--> Ca4CaM_CaMK
    ##  Binding of Ca to CaM-CaMKIIP.
    (k_K1C_on * k_K2C_on / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)) * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart), k_K1C_off * k_K2C_off / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), CaM0_CaMKP <--> Ca2CaM_C_CaMKP
    (k_K1N_on * k_K2N_on / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)) * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart), k_K1N_off * k_K2N_off / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), CaM0_CaMKP <--> Ca2CaM_N_CaMKP
    (k_K1C_on * k_K2C_on / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)) * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart), k_K1C_off * k_K2C_off / (k_K1C_off + k_K2C_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), Ca2CaM_C_CaMKP <--> Ca4CaM_CaMKP
    (k_K1N_on * k_K2N_on / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart)) * ($Cai_mean)^2 * (t <= tstop) * (t >= tstart), k_K1N_off * k_K2N_off / (k_K1N_off + k_K2N_on * ($Cai_mean) * (t <= tstop) * (t >= tstart))), Ca2CaM_N_CaMKP <--> Ca4CaM_CaMKP
    ##  Binding of CaM to CaMKII or CaMII-P
    (kCaM0_on, kCaM0_off), CaM0 + CaMK <--> CaM0_CaMK
    (kCaM2C_on, kCaM2C_off), Ca2CaM_C + CaMK <--> Ca2CaM_C_CaMK
    (kCaM2N_on, kCaM2N_off), Ca2CaM_N + CaMK <--> Ca2CaM_N_CaMK
    (kCaM4_on, kCaM4_off), Ca4CaM + CaMK <--> Ca4CaM_CaMK
    (kCaM0P_on, kCaM0P_off), CaM0 + CaMKP <--> CaM0_CaMKP
    (kCaM2CP_on, kCaM2CP_off), Ca2CaM_C + CaMKP <--> Ca2CaM_C_CaMKP
    (kCaM2NP_on, kCaM2NP_off), Ca2CaM_N + CaMKP <--> Ca2CaM_N_CaMKP
    (kCaM4P_on, kCaM4P_off), Ca4CaM + CaMKP <--> Ca4CaM_CaMKP
    ##  Phosphorylation CaMXCaMKII -> CaMXCaMKIIP.
    (k_phosCaM * (Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + Ca4CaM_CaMK + Ca4CaM_CaMKP + CaMKP + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, k_PB), Ca2CaM_C_CaMK <--> Ca2CaM_C_CaMKP
    (k_phosCaM * (Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + Ca4CaM_CaMK + Ca4CaM_CaMKP + CaMKP + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, k_OXPOX), Ca2CaM_N_CaMK <--> Ca2CaM_N_CaMKP
    #(kbi*kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)/(kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)+0.01851e3),k_PB), Ca2CaM_C_CaMK <--> Ca2CaM_C_CaMKP
    #(kbi*kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)/(kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)+0.01851e3),k_OXPOX),Ca2CaM_N_CaMK <--> Ca2CaM_N_CaMKP
    ##  Dephosphorylation CaMKP -> CaMK
    k_dephospho, CaMKP --> CaMK
    ## Adjustment for Oxdization

    (k_phosCaM * (Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + Ca4CaM_CaMK + Ca4CaM_CaMKP + CaMKP + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, k_PB), Ca4CaM_CaMK <--> Ca4CaM_CaMKP #
    (k_phosCaM * (Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + Ca4CaM_CaMK + Ca4CaM_CaMKP + CaMKP + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, k_OXPOX), Ca4CaM_CaMKOX <--> Ca4CaM_CaMKPOX #
    #(kbi*kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)/(kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)+0.01851e3),k_PB), Ca4CaM_CaMK <--> Ca4CaM_CaMKP #
    #(kbi*kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)/(kbi/kib/(1/(Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+Ca4CaM_CaMK+Ca4CaM_CaMKP+CaMKP+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)-1)+0.01851e3),k_OXPOX), Ca4CaM_CaMKOX <--> Ca4CaM_CaMKPOX #
    (ROS * k_BOX, k_OXB), Ca4CaM_CaMK <--> Ca4CaM_CaMKOX
    (ROS * k_POXP, k_OXPP), Ca4CaM_CaMKP <--> Ca4CaM_CaMKPOX
end

###########################  Parameters  ###########################
CaMT = 30 #Total calmodulin concentration.
CaMKII_T = 70 #Total CaMKII concentration.

binding_To_PCaMK = 0.1  ## 0.1
decay_CaM = 3000 # seconds
phospho_rate = 1
phosphatase = 1

rn_osys = convert(ODESystem, ca_model)
@named sys = extend(osys, rn_osys)

sys = structural_simplify(sys)

@unpack Cai, i_CaJSR, CaNSR, V, Nai, Ki, i_b, i_g, i_d, i_f, i_fca, i_y, i_r, i_s, i_sslow,
i_Nam, i_Nah, i_Naj, i_PO1, i_PO2, i_PC2, i_nKs, i_CK1, i_CK2, i_OK, i_IK,                                                  # ecc_ODEfile
#CaM_dyad, Ca2CaM_dyad, Ca4CaM_dyad, CaMB_dyad, Ca2CaMB_dyad, Ca4CaMB_dyad, Pb2_dyad, Pb_dyad,
#Pt_dyad, Pt2_dyad, Pa_dyad, Ca4CaN_dyad, CaMCa4CaN_dyad, Ca2CaMCa4CaN_dyad, Ca4CaMCa4CaN_dyad,                              # camdyad_ODEfile
CaM_sl, Ca2CaM_sl, Ca4CaM_sl, CaMB_sl, Ca2CaMB_sl, Ca4CaMB_sl, Pb2_sl, Pb_sl,
Pt_sl, Pt2_sl, Pa_sl, Ca4CaN_sl, CaMCa4CaN_sl, Ca2CaMCa4CaN_sl, Ca4CaMCa4CaN_sl,                                            # camsl_ODEfile
CaM_cyt, Ca2CaM_cyt, Ca4CaM_cyt, CaMB_cyt, Ca2CaMB_cyt, Ca4CaMB_cyt, Pb2_cyt, Pb_cyt,
Pt_cyt, Pt2_cyt, Pa_cyt, Ca4CaN_cyt, CaMCa4CaN_cyt, Ca2CaMCa4CaN_cyt, Ca4CaMCa4CaN_cyt,                                     # camcyt_ODEfile
LCC_PKAp, RyR2809p, RyR2815p, PLBT17p, LCC_CKslp,                                                                           # camkii_ODEfile LCC_CKdyadp,
LR, LRG, RG, b1AR_S464, b1AR_S301, GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, PDEp, cAMPtot, RC_I, RCcAMP_I,
RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII,                         # bar_ODEfile
PKACII_PKI, I1p_PP1, I1ptot, PLBp, PLMp, LCCap, LCCbp, RyRp, TnIp, KS79, KS80, KSp, CFTRp, KURp,
CaMK, CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, Ca4CaM_CaMKP,
CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, CaMKP, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, k_1C_on, k_1C_off, k_2C_on, k_2C_off,
k_1N_on, k_1N_off, k_2N_on, k_2N_off, k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off, k_K1N_on, k_K1N_off,
k_K2N_on, k_K2N_off, kCaM0_on, kCaM2C_on, kCaM2N_on, kCaM4_on, kCaM0_off, kCaM2C_off, kCaM2N_off, kCaM4_off,
kCaM0P_on, kCaM2CP_on, kCaM2NP_on, kCaM4P_on, kCaM0P_off, kCaM2CP_off, kCaM2NP_off, kCaM4P_off, # kbi, kib,
k_PB, k_OXPOX, k_dephospho, ROS, k_phosCaM, CaMKII_T, k_BOX, k_OXB, k_POXP, k_OXPP, tstop, tstart = sys

tspan = (0.0, 300e3)

u0 = [Cai[1] => 0.2556, Cai[2] => 0.25574, Cai[3] => 0.25587, Cai[4] => 0.25599, Cai[5] => 0.25609,
    Cai[6] => 0.25618, Cai[7] => 0.25625, Cai[8] => 0.25631, Cai[9] => 0.25636, Cai[10] => 0.25639,
    Cai[11] => 0.25642, Cai[12] => 0.25643, Cai[13] => 0.25642, Cai[14] => 0.25641, Cai[15] => 0.25639,
    Cai[16] => 0.25635, Cai[17] => 0.25631, Cai[18] => 0.25625, Cai[19] => 0.25619, Cai[20] => 0.25611,
    Cai[21] => 0.25602, Cai[22] => 0.25593, Cai[23] => 0.25583, Cai[24] => 0.25571, Cai[25] => 0.25559,
    Cai[26] => 0.25546, Cai[27] => 0.25532, Cai[28] => 0.25517, Cai[29] => 0.25502, Cai[30] => 0.25485,
    Cai[31] => 0.25468, Cai[32] => 0.2545, Cai[33] => 0.25431, Cai[34] => 0.25412, Cai[35] => 0.25392,
    Cai[36] => 0.25371, Cai[37] => 0.25349, Cai[38] => 0.25326, Cai[39] => 0.25303, Cai[40] => 0.2528,
    Cai[41] => 0.25255, Cai[42] => 0.2523, Cai[43] => 0.25204, Cai[44] => 0.25178, Cai[45] => 0.25151,
    i_CaJSR => 613.87556, CaNSR => 619.09843, V => -68.79268, Nai => 13838.37602, Ki => 150952.75035,
    i_b => 0.00305, i_g => 0.61179, i_d => 0.00033, i_f => 0.99869, i_fca => 0.9911, i_y => 0.07192,
    i_r => 0.00702, i_s => 0.96604, i_sslow => 0.22156, i_Nam => 0.02506, i_Nah => 0.22242, i_Naj => 0.19081,
    i_PO1 => 0.0037, i_PO2 => 0.0, i_PC2 => 0.0, i_nKs => 0.09243, i_CK1 => 0.00188, i_CK2 => 0.00977,
    i_OK => 0.26081, i_IK => 0.07831, CaM_sl => 0.03744, Ca2CaM_sl => 0.00031, Ca4CaM_sl => 0.0,
    CaMB_sl => 4.20703, Ca2CaMB_sl => 10.08438, Ca4CaMB_sl => 0.00111, Pb2_sl => 6.0e-5, Pb_sl => 0.0, Pt_sl => 0.0,
    Pt2_sl => 0.0, Pa_sl => 0.0, Ca4CaN_sl => 0.00037, CaMCa4CaN_sl => 0.0, Ca2CaMCa4CaN_sl => 0.0,
    Ca4CaMCa4CaN_sl => 0.00141, CaM_cyt => 0.03779, Ca2CaM_cyt => 0.00031, Ca4CaM_cyt => 0.0, CaMB_cyt => 2.5048,
    Ca2CaMB_cyt => 2.789, Ca4CaMB_cyt => 0.00032, Pb2_cyt => 6.0e-5, Pb_cyt => 0.0, Pt_cyt => 0.0, Pt2_cyt => 0.0,
    Pa_cyt => 0.0, Ca4CaN_cyt => 0.00069, CaMCa4CaN_cyt => 0.0, Ca2CaMCa4CaN_cyt => 1.0e-5, Ca4CaMCa4CaN_cyt => 1.0e-5,
    LCC_PKAp => 16.45439, RyR2809p => 297.35744, RyR2815p => 130.5212, PLBT17p => 1.87769,
    LCC_CKslp => 1.0e-5, LR => 0.0, LRG => 0.0, RG => 0.00048, b1AR_S464 => 0.0, b1AR_S301 => 0.00065,
    GsaGTPtot => 0.00961, GsaGDP => 0.00063, Gsby => 0.01002, AC_GsaGTP => 0.00142, PDEp => 0.00223,
    cAMPtot => 1.02286, RC_I => 0.80424, RCcAMP_I => 0.14186, RCcAMPcAMP_I => 0.00449, RcAMPcAMP_I => 0.22889,
    PKACI => 0.08583, PKACI_PKI => 0.14356, RC_II => 0.051, RCcAMP_II => 0.009, RCcAMPcAMP_II => 0.00028,
    RcAMPcAMP_II => 0.0577, PKACII => 0.02159, PKACII_PKI => 0.0361, I1p_PP1 => 0.07292, I1ptot => 0.07301,
    PLBp => 8.49358, PLMp => 5.62885, LCCap => 0.0055, LCCbp => 0.00628, RyRp => 0.02763, TnIp => 4.41392,
    KURp => 0.01095, KS79 => 0.00153, KS80 => 0.00153, KSp => 0.00184, CFTRp => 0.00406,
    CaM0 => 1000, Ca2CaM_C => 0.0, Ca2CaM_N => 0.0, Ca4CaM => 0.0,
    CaM0_CaMK => 0.0, Ca2CaM_C_CaMK => 0.0, Ca2CaM_N_CaMK => 0.0, Ca4CaM_CaMK => 0.0,
    CaM0_CaMKP => 0.0, Ca2CaM_C_CaMKP => 0.0, Ca2CaM_N_CaMKP => 0.0, Ca4CaM_CaMKP => 0.0,
    CaMK => 0.0, CaMKP => CaMKII_T, Ca4CaM_CaMKOX => 0.0, Ca4CaM_CaMKPOX => 0.0]


ks = [k_1C_on => 5e-3, k_1C_off => 50e-3, k_2C_on => 10e-3, k_2C_off => 10e-3, k_1N_on => 100e-3,
    k_1N_off => 2000e-3, k_2N_on => 200e-3, k_2N_off => 500e-3, k_K1C_on => 44e-3, k_K1C_off => 33e-3,
    k_K2C_on => 44e-3, k_K2C_off => 0.8e-3, k_K1N_on => 76e-3, k_K1N_off => 300e-3, k_K2N_on => 76e-3,
    k_K2N_off => 20e-3, kCaM0_on => 3.8e-6, kCaM2C_on => 0.92e-3, kCaM2N_on => 0.12e-3, kCaM4_on => 30e-3, #20
    kCaM0_off => 5.5e-3, kCaM2C_off => 6.8e-3, kCaM2N_off => 1.7e-3, kCaM4_off => 1.5e-3,
    kCaM0P_on => 3.8e-6 * binding_To_PCaMK, kCaM2CP_on => 0.92e-3 * binding_To_PCaMK,
    kCaM2NP_on => 0.12e-3 * binding_To_PCaMK, kCaM4P_on => 30e-3 * binding_To_PCaMK,
    kCaM0P_off => 1 / decay_CaM, kCaM2CP_off => 1 / decay_CaM, kCaM2NP_off => 1 / decay_CaM,
    kCaM4P_off => 1 / decay_CaM, k_dephospho => (1 / 6000) * phosphatase, k_phosCaM => 2e-3 * phospho_rate, CaMKII_T => 70,
    k_BOX => 2.91e-4, k_PB => 0.00003, k_OXPOX => 0.00003,
    k_OXB => 2.23e-5, k_POXP => 2.91e-4, k_OXPP => 2.23e-5, tstop => 260e3, tstart => 130e3, ROS => 0]

oprob = ODEProblem(sys, u0, tspan, ks)


## To generate reference figures
psmap = Dict(k => i for (i, k) in enumerate(parameters(sys)))
usmap = Dict(k => i for (i, k) in enumerate(states(sys)))

# Ca4CaM_CaMK+Ca4CaM_CaMKP+Ca4CaM_CaMKOX+Ca4CaM_CaMKPOX+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+CaMKP
idxROS = psmap[ROS]
idxCaM0 = usmap[CaM0]
idxB = usmap[Ca4CaM_CaMK]
idxP = usmap[Ca4CaM_CaMKP]
idxOX = usmap[Ca4CaM_CaMKOX]
idxPOX = usmap[Ca4CaM_CaMKPOX]
idxP0 = usmap[CaM0_CaMKP]
idxP2C = usmap[Ca2CaM_C_CaMKP]
idxP2N = usmap[Ca2CaM_N_CaMKP]
idxKP = usmap[CaMKP]
idxK = usmap[CaMK]
idxK0 = usmap[CaM0_CaMK]
idxK2C = usmap[Ca2CaM_C_CaMK]
idxK2N = usmap[Ca2CaM_N_CaMK]

# Warm up
sol = solve(oprob, TRBDF2(), abstol=1e-9, reltol=1e-9, alg_hints=[:stiff], tstops=0:1000:tspan[end])


actK = [Ca4CaM_CaMK, Ca4CaM_CaMKP, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, CaMKP]

plot(sol, idxs=[sum(actK)], linewidth=1.2, xlabel="Time(ms)", ylabel="Concentration(uM)", label="simulation", denseplot=false, c=:steelblue2, title="CaMKII Activity (ISO=0.0/ROS=0/1Hz)", size=(800, 600))

B = values(sol(250e3)[usmap[Ca4CaM_CaMK]])
P = values(sol(250e3)[usmap[Ca4CaM_CaMKP]])
OX = values(sol(250e3)[usmap[Ca4CaM_CaMKOX]])
POX = values(sol(250e3)[usmap[Ca4CaM_CaMKPOX]])
P0 = values(sol(250e3)[usmap[CaM0_CaMKP]])
P2C = values(sol(250e3)[usmap[Ca2CaM_C_CaMKP]])
P2N = values(sol(250e3)[usmap[Ca2CaM_N_CaMKP]])
KP = values(sol(250e3)[usmap[CaMKP]])
K = values(sol(250e3)[usmap[CaMK]])
K0 = values(sol(250e3)[usmap[CaM0_CaMK]])
K2C = values(sol(250e3)[usmap[Ca2CaM_C_CaMK]])
K2N = values(sol(250e3)[usmap[Ca2CaM_N_CaMK]])

XLSX.openxlsx("figure fitting.xlsx", mode="rw") do f # if write a new file -> w, write an existing file -> rw
    sheet = f[2]
    sheet["A20"] = 1000
    sheet["B20"] = (P + P0 + P2C + P2N + KP + POX) / (B + OX + P + POX + P0 + P2C + P2N + KP + K0 + K2C + K2N)
    sheet["C20"] = (P + P0 + P2C + P2N + KP + POX) / (B + OX + P + POX + P0 + P2C + P2N + KP + K0 + K2C + K2N + K)
end

print((P + P0 + P2C + P2N + KP + POX) / (B + OX + P + POX + P0 + P2C + P2N + KP + K0 + K2C + K2N + K))
print((P + P0 + P2C + P2N + KP + POX) / (B + OX + P + POX + P0 + P2C + P2N + KP + K0 + K2C + K2N))

kinase_act = map(0.0:0.001:0.002) do cam0#ros
    t = 330e3
    p = copy(oprob.p)
    p[idxCaM0] = cam0
    sol = solve(remake(oprob, p=p), TRBDF2(), abstol=1e-9, reltol=1e-9, alg_hints=[:stiff], tstops=0:1000:tspan[end])
    B = sol(t)[idxB]
    P = sol(t)[idxP]
    OX = sol(t)[idxOX]
    POX = sol(t)[idxPOX]
    P0 = sol(t)[idxP0]
    P2C = sol(t)[idxP2C]
    P2N = sol(t)[idxP2N]
    KP = sol(t)[idxKP]
    K = sol(t)[idxK]
    K0 = sol(t)[idxK0]
    K2C = sol(t)[idxK2C]
    K2N = sol(t)[idxK2N]
    (B + P + P0 + P2C + P2N + KP + OX + POX) / (B + OX + P + POX + P0 + P2C + P2N + KP + K0 + K2C + K2N)
end

#print(sol(0)[14])

XLSX.openxlsx("kinase_act.xls", mode="rw") do f # if write a new file -> w, write an existing file -> rw
    sheet = f[1]

    # Used for renaming the file for convience
    #XLSX.rename!(sheet, "new_sheet")

    # This will add output along A column
    #sheet["H1"] = "(B+P+P0+P2C+P2N+KP+OX+POX) / (B+OX+P+POX+P0+P2C+P2N+KP)#0:0.1:1"
    sheet["H2", dim=1] = kinase_act
end

plot(scatter(x=[0, 0.01, 0.1, 1]), kinase_act * 100)


sol = solve(oprob, TRBDF2(), abstol=1e-9, reltol=1e-9, alg_hints=[:stiff], tstops=0:1000:tspan[end])
#FBDF, maxiters=Int(1e8)
## colors: orchid, slateblue1, tan1, darkseagreen4, indianred3, navajowhite4
#steelblue2(blue),goldenrod1(oranyellow),lightcoral(red),lightseagreen,slateblue1(purple),cadetblue4,pink1
actK = [Ca4CaM_CaMK, Ca4CaM_CaMKP, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, CaMKP]
plot(sol, idxs=mean(skipmissing(Cai)), linewidth=1.5, title="Calcium Transient (3Hz)", xlabel="Time(ms)", ylabel="[Ca2+](uM)", label="ISO=0.0", xlim=(3999, 8001))

plot(sol, idxs=mean(skipmissing(Cai)), linewidth=1.5, title="Calcium Transient (1Hz)", xlabel="Time(ms)", ylabel="[Ca2+](mM)", label="ISO=0.1", xlim=(0, 1.05e5), denseplot=false)

plot(sol, idxs=CaMK, linewidth=1.5, title="Action Potential (Control)", xlabel="Time(ms)", ylabel="Voltage (mV)", label="ISO=0", denseplot=false)

plot!(sol, idxs=[Ca4CaM_CaMK + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + CaMKP], linewidth=1.2, xlabel="Time(ms)", ylabel="Concentration(uM)", title="CaMKII Activity (ISO=0.0/Freq=1Hz)", label="ROS=0uM", denseplot=false, ylim=(0, 33)) #ROS=0uM
#
plot!(sol, idxs=[Ca4CaM_CaMK + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + CaMKP], linewidth=1, xlabel="Time(ms)", ylabel="Concentration(uM)", title="CaMKII Activity (ISO=0.1/ROS=0uM)", label="1Hz", denseplot=false, c=:dodgerblue) #ROS=1uM
#
plot!(sol, idxs=[Ca4CaM_CaMK + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + CaMKP], linewidth=1, xlabel="Time(ms)", ylabel="Concentration(uM)", title="CaMKII Activity (ISO=0.1/ROS=0uM)", label="2Hz", denseplot=false, c=:orchid) #ROS=2uM

plot(sol, idxs=[Ca4CaM_CaMK + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + CaMKP], linewidth=1, xlabel="Time(ms)", ylabel="Concentration(uM)", title="CaMKII Activity (ISO=0.1/ROS=0uM)", label="3Hz", denseplot=false, c=:slateblue1, ylim=(-1, 35), size=(800, 600)) #ROS=3uM


#, linewidth=1
plot(sol, idxs=[sum(actK)], linewidth=1.2, xlabel="Time(ms)", ylabel="Concentration(uM)", label="250s", denseplot=false, c=:steelblue2, title="CaMKII Activity (ISO=0.0/ROS=0/2Hz)", size=(800, 600))
plot!(sol, idxs=[sum(actK)], linewidth=1.2, xlabel="Time(ms)", ylabel="Concentration(uM)", label="200s", denseplot=false, c=:lightcoral)
plot!(sol, idxs=[sum(actK)], linewidth=1.2, xlabel="Time(ms)", ylabel="Concentration(uM)", label="150s", denseplot=false, c=:goldenrod1)
plot!(sol, idxs=[sum(actK)], linewidth=1.2, xlabel="Time(ms)", ylabel="Concentration(uM)", label="100s", denseplot=false, c=:lightseagreen)
plot!(sol, idxs=[sum(actK)], linewidth=1.2, xlabel="Time(ms)", ylabel="Concentration(uM)", label="50s", denseplot=false, c=:slateblue1)
