### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 67a68b30-9d7f-11ef-1c30-bf3107341c1b
begin
	using ProgressLogging
	using Catalyst
	using DifferentialEquations
	using Plots
	using ModelingToolkit
	using Sundials
	using Statistics
	using XLSX
	import NaNMath as nm
	using ModelingToolkit: t_nounits as t, D_nounits as D
end

# ╔═╡ 0583b4d2-6c54-4621-8d98-3df7a69e306d
function beta_cai(Cai, TnI_PKAp)
    fracTnIpo = 0.062698  # Derived quantity (TnI_PKAp(baseline)/TnItot)
    fPKA_TnI = (1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIpo)) # Max effect +61#
    TrpnTotal = 35
    KmTrpn = 0.5 / fPKA_TnI
    CmdnTotal = 50
    KmCmdn = 2.38
    return inv(1 + TrpnTotal * KmTrpn / (Cai + KmTrpn)^2 + CmdnTotal * KmCmdn / (Cai + KmCmdn)^2)
end

# ╔═╡ 447af9d0-daa7-4855-be20-c38dac16e8a0
function get_Morotti_equations()
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
    # Cell geometry (μm)
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

# ╔═╡ ee3c7f42-4e1d-4ee5-b0b2-ffbe8b4fd197
eq_morotti = get_Morotti_equations();

# ╔═╡ e9271d60-9abb-49b2-90ea-d5c543834e8e
@named osys = ODESystem(eq_morotti, t);

# ╔═╡ 34816d55-a8a9-411f-9a03-adeceb6518f9
function get_camkii_model()
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
	
end

# ╔═╡ fac108bf-63be-4524-8622-71ce4cddab60
rn_osys = convert(ODESystem, get_camkii_model());

# ╔═╡ a798d3ac-5841-4811-9009-776f6f02b7de
begin
	@named sys = extend(osys, rn_osys);
	sys = structural_simplify(sys)
end

# ╔═╡ d8e639a5-c0cb-476c-b5a4-e3e841bcf4a0
 tspan = (0.0, 300e3)

# ╔═╡ 838f5859-0e71-42f6-97b7-6f00732438a2
oprob = let

CaMT = 30 #Total calmodulin concentration.
CaMKII_T = 70 #Total CaMKII concentration.

binding_To_PCaMK = 0.1  ## 0.1
decay_CaM = 3000 # seconds
phospho_rate = 1
phosphatase = 1

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

ks = [k_1C_on => 5e-3, k_1C_off => 50e-3, k_2C_on => 10e-3, k_2C_off => 10e-3, k_1N_on => 100e-3, k_1N_off => 2000e-3, k_2N_on => 200e-3, k_2N_off => 500e-3, k_K1C_on => 44e-3, k_K1C_off => 33e-3, k_K2C_on => 44e-3, k_K2C_off => 0.8e-3, k_K1N_on => 76e-3, k_K1N_off => 300e-3, k_K2N_on => 76e-3, k_K2N_off => 20e-3, kCaM0_on => 3.8e-6, kCaM2C_on => 0.92e-3, kCaM2N_on => 0.12e-3, kCaM4_on => 30e-3, kCaM0_off => 5.5e-3, kCaM2C_off => 6.8e-3, kCaM2N_off => 1.7e-3, kCaM4_off => 1.5e-3, kCaM0P_on => 3.8e-6 * binding_To_PCaMK, kCaM2CP_on => 0.92e-3 * binding_To_PCaMK, kCaM2NP_on => 0.12e-3 * binding_To_PCaMK, kCaM4P_on => 30e-3 * binding_To_PCaMK, kCaM0P_off => 1 / decay_CaM, kCaM2CP_off => 1 / decay_CaM, kCaM2NP_off => 1 / decay_CaM, kCaM4P_off => 1 / decay_CaM, k_dephospho => (1 / 6000) * phosphatase, k_phosCaM => 2e-3 * phospho_rate, CaMKII_T => 70, k_BOX => 2.91e-4, k_PB => 0.00003, k_OXPOX => 0.00003, k_OXB => 2.23e-5, k_POXP => 2.91e-4, k_OXPP => 2.23e-5, tstop => 260e3, tstart => 130e3, ROS => 0]


oprob = ODEProblem(sys, u0, tspan, ks)
end

# ╔═╡ c31a50c7-b116-485a-aa4c-57bc37ab3812
alg = TRBDF2()

# ╔═╡ 6c80bb29-dda0-4f1a-8865-47bb5f2bab72
@time sol = solve(oprob, alg, abstol=1e-6, reltol=1e-6, tstops=0:1000:tspan[end], progress=true)

# ╔═╡ ef4446ea-5210-4586-a7c2-3ab21a3594f0
length(sol)

# ╔═╡ 94bfa7bc-458c-4fc1-b31d-64261d8fb6a6
plot(sol, idxs=sys.Cai[45]*1000, tspan=(200e3, 201e3), xlabel="Time (ms)", ylabel="Conc. (nM)", label="Ca (SL)", ylims=(250, 600))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Catalyst = "479239e8-5488-4da2-87a7-35f2df7eef83"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"
XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[compat]
Catalyst = "~14.4.1"
DifferentialEquations = "~7.15.0"
ModelingToolkit = "~9.49.0"
NaNMath = "~1.0.2"
Plots = "~1.40.8"
ProgressLogging = "~0.1.4"
Statistics = "~1.11.1"
Sundials = "~4.26.0"
XLSX = "~0.10.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "84005ab61a098e7dcd69a24335edda1a95463c6f"

[[deps.ADTypes]]
git-tree-sha1 = "eea5d80188827b35333801ef97a40c2ed653b081"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.9.0"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "b392ede862e506d451fc1616e79aa6f4c673dab8"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.38"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "492681bc44fac86804706ddb37da10880a2bd528"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.10.4"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "a2c85f53ddcb15b4099da59867868bd40f005579"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.7.5"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "d434647f798823bcae510aee0bc0401927f64391"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "1.1.1"
weakdeps = ["BandedMatrices"]

    [deps.BlockArrays.extensions]
    BlockArraysBandedMatricesExt = "BandedMatrices"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "BoundaryValueDiffEqFIRK", "BoundaryValueDiffEqMIRK", "BoundaryValueDiffEqShooting", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "98da8bd76b89a4ae1b8dde9fc6dcd75dcd6b5282"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.12.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.BoundaryValueDiffEqCore]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "ForwardDiff", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "768423f5fa51adc47ea5bd2540497fc785c2f0cc"
uuid = "56b672f2-a5fe-4263-ab2d-da677488eb3a"
version = "1.0.1"

[[deps.BoundaryValueDiffEqFIRK]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "35e1e7822d1c77d85ecf568606ca64d60fbd39de"
uuid = "85d9eb09-370e-4000-bb32-543851f73618"
version = "1.0.2"

[[deps.BoundaryValueDiffEqMIRK]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "e1fa0dee3d8eca528ab96e765a52760fd7466ffa"
uuid = "1a22d4ce-7765-49ea-b6f2-13c8438986a6"
version = "1.0.1"

[[deps.BoundaryValueDiffEqShooting]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "fac04445ab0fdfa29b62d84e1af6b21334753a94"
uuid = "ed55bfe0-3725-4db6-871e-a1dc9f42a757"
version = "1.0.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+2"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "0157e592151e39fa570645e2b2debcdfb8a0f112"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.4.3"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.Catalyst]]
deps = ["Combinatorics", "DataStructures", "DiffEqBase", "DocStringExtensions", "DynamicPolynomials", "DynamicQuantities", "Graphs", "JumpProcesses", "LaTeXStrings", "Latexify", "LinearAlgebra", "MacroTools", "ModelingToolkit", "Parameters", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SymbolicUtils", "Symbolics", "Unitful"]
git-tree-sha1 = "dafd059c1d80ba02006978cd05ea7264ff1f362e"
uuid = "479239e8-5488-4da2-87a7-35f2df7eef83"
version = "14.4.1"

    [deps.Catalyst.extensions]
    CatalystBifurcationKitExtension = "BifurcationKit"
    CatalystCairoMakieExtension = "CairoMakie"
    CatalystGraphMakieExtension = "GraphMakie"
    CatalystHomotopyContinuationExtension = "HomotopyContinuation"
    CatalystStructuralIdentifiabilityExtension = "StructuralIdentifiability"

    [deps.Catalyst.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
    GraphMakie = "1ecd5474-83a3-4783-bb4f-06765db800d2"
    HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327"
    StructuralIdentifiability = "220ca800-aa68-49bb-acd8-6037fa93a544"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodecInflate64]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "d981a6e8656b1e363a2731716f46851a2257deb7"
uuid = "6309b1aa-fc58-479c-8956-599a07234577"
version = "0.1.3"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "13951eb68769ad1cd460cdb2e64e5e95f1bf123d"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.27.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "c7acce7a7e1078a20a285211dd73cd3941a871d6"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.0"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "3e532ae5ac37b4be60228d63e28e0bc9d5900bcc"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.1"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "PrecompileTools"]
git-tree-sha1 = "3faae67b8899797592335832fccf4b3c80bb04fa"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.15"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "066f60231c1b0ae2905ffd2651e207accd91f627"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.48.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "TruncatedStacktraces"]
git-tree-sha1 = "f8eefbb7e910f59087c4bb09ce670f235758ee4a"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.158.3"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "Functors", "LinearAlgebra", "Markdown", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "b1f970a2873a2cf76ce35fb0ed2b755a11b31052"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.1.0"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ab1e6515ce15f01316a9825b02729fefa51726bd"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.23.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "d55af9d6b51c54f81ae30d1a463206d32cc4c24a"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.15.0"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "ba137efeddd4b6e6a7154f2a92d2922a0057486f"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.6.21"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = "Enzyme"
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = "ForwardDiff"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = "PolyesterForwardDiff"
    DifferentiationInterfaceReverseDiffExt = "ReverseDiff"
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DispatchDoctor]]
deps = ["MacroTools", "Preferences"]
git-tree-sha1 = "453df2ce2aef7de59c69a56d31dcd2ec3384dd77"
uuid = "8d63f2c5-f18a-4cf2-ba9d-b3f60fc568c8"
version = "0.4.17"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.DispatchDoctor.extensions]
    DispatchDoctorChainRulesCoreExt = "ChainRulesCore"
    DispatchDoctorEnzymeCoreExt = "EnzymeCore"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3101c32aab536e7a27b1763c0797dba151b899ad"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.113"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "490392af2c7d63183bfa2c8aaa6ab981c5ba7561"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.14"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "bbf1ace0781d9744cb697fb856bd2c3f6568dadb"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.0"

[[deps.DynamicQuantities]]
deps = ["DispatchDoctor", "TestItems", "Tricks"]
git-tree-sha1 = "9f826f051e3d2c76f924d39c9aa4fcfd18cc9256"
uuid = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
version = "1.2.0"

    [deps.DynamicQuantities.extensions]
    DynamicQuantitiesLinearAlgebraExt = "LinearAlgebra"
    DynamicQuantitiesMeasurementsExt = "Measurements"
    DynamicQuantitiesScientificTypesExt = "ScientificTypes"
    DynamicQuantitiesUnitfulExt = "Unitful"

    [deps.DynamicQuantities.weakdeps]
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    ScientificTypes = "321657f4-b219-11e9-178b-2701a2544e81"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "04c777af6ef65530a96ab68f0a81a4608113aa1d"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.5"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "380053d61bb9064d6aa4a9777413b40429c79901"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.2.0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "3f03d94c71126b6cfe20d3cbcc41c5cd27e1c419"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.4"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "cbf5edddb61a43669710cbc2241bc08b36d9e660"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.4"

[[deps.FastPower]]
git-tree-sha1 = "58c3431137131577a7c379d00fea00be524338fb"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.1"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FindFirstFunctions]]
git-tree-sha1 = "670e1d9ceaa4a3161d32fe2d2fb2177f8d78b330"
uuid = "64ca27bc-2ba2-4a57-88aa-44e436879224"
version = "1.4.1"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "b10bdafd1647f57ace3885143936749d61638c3b"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.26.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a9ce73d3c827adab2d70bf168aaece8cce196898"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.37"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "64d8e93700c7a3f28f717d265382d52fac9fa1c1"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.12"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee28ddcd5517d54e417182fec3886e7412d3926f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f31929b9e67066bee48eec8b03c0df47d31a74b3"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.8+0"

[[deps.GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random", "libblastrampoline_jll"]
git-tree-sha1 = "c4f9c87b74aedf20920034bd4db81d0bffc527d2"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.3.14"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1dc470db8b1131cfc7fb4c115de89fe391b9e780"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "bc3f416a965ae61968c20d0ad867556367f2817d"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.9"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InputBuffers]]
git-tree-sha1 = "d5c278bee2efd4fda62725a4a794d7e5f55e14f1"
uuid = "0c81fc1b-5583-44fc-8770-48be1e1cca08"
version = "1.0.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Glob", "PrecompileTools", "TOML", "Tokenize"]
git-tree-sha1 = "59cf7ad64f1b0708a4fa4369879d33bad3239b56"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "1.0.62"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "c3a2cb6f968404ed3b1d5382bbdd7b7d83966598"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.14.0"
weakdeps = ["FastBroadcast"]

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "4f20a2df85a9e5d55c9e84634bbf808ed038cabd"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.8"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "360f6039babd6e4d6364eff0d4fc9120834a2d9a"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.2.1"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6ce1e19f3aec9b59186bdf06cdf3c4fc5f5f3e6"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.50.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "61dfdba58e585066d8bce214c5a51eaa0539f269"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "591de175461afd8323aa24b7686062574527aa3a"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.36.2"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8084c25a250e00ae427a379a5b607e7aed96a2dd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.171"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "16a726dba99685d9e94c8d0a8f655383121fc608"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.0.1"
weakdeps = ["BandedMatrices"]

    [deps.MatrixFactorizations.extensions]
    MatrixFactorizationsBandedMatricesExt = "BandedMatrices"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "BlockArrays", "Combinatorics", "CommonSolve", "Compat", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "DynamicQuantities", "ExprTools", "Expronicon", "FindFirstFunctions", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "Graphs", "InteractiveUtils", "JuliaFormatter", "JumpProcesses", "Latexify", "Libdl", "LinearAlgebra", "MLStyle", "NaNMath", "NonlinearSolve", "OffsetArrays", "OrderedCollections", "PrecompileTools", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "SciMLStructures", "Serialization", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "URIs", "UnPack", "Unitful"]
git-tree-sha1 = "db309b354ade7fe3deb97d95ca241ed8f1096e91"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "9.49.0"

    [deps.ModelingToolkit.extensions]
    MTKBifurcationKitExt = "BifurcationKit"
    MTKChainRulesCoreExt = "ChainRulesCore"
    MTKDeepDiffsExt = "DeepDiffs"
    MTKHomotopyContinuationExt = "HomotopyContinuation"
    MTKLabelledArraysExt = "LabelledArrays"

    [deps.ModelingToolkit.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DeepDiffs = "ab62b9b5-e342-54a8-a765-a90f495de1a6"
    HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "90077f1e79de8c9c7c8a90644494411111f4e07b"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.5.2"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "4d8944f32db2b07a2bdf8477e878bcb9c9ea2308"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.15.1"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d9b79c4eed437421ac4285148fcadf42e0700e89"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.4"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "cd892f12371c287dc50d6ad3af075b088b6f2d48"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.89.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "8e3c5978d0531a961f70d2f2730d1d16ed3bbd12"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.1.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "b4498d40bf35da0b6d22652ff2e9d8820590b3c6"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.1.2"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "5e8c500a80674850543394ce3c745b73ad51fea0"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.10.0"
weakdeps = ["EnzymeCore"]

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"

[[deps.OrdinaryDiffEqDefault]]
deps = ["DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport"]
git-tree-sha1 = "c8223e487d58bef28a3535b33ddf8ffdb44f46fb"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.1.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "e63ec633b1efa99e3caa2e26a01faaa88ba6cef9"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.1.0"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "TruncatedStacktraces"]
git-tree-sha1 = "4dbce3f9e6974567082ce5176e21aab0224a69e9"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.1.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "f63938b8e9e5d3a05815defb3ebdbdcf61ec0a74"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.1.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "048bcccc8f59c20d5b4ad268eef4d7d21c005a94"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.2.1"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["DiffEqBase", "FastBroadcast", "FastPower", "GenericLinearAlgebra", "GenericSchur", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polynomials", "RecursiveArrayTools", "Reexport", "RootedTrees", "SciMLOperators", "Symbolics"]
git-tree-sha1 = "5735f4c094dff311f5064d1a351da9669e4647e3"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.2.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "a7cc74d3433db98e59dc3d58bc28174c6c290adf"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.1.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "925a91583d1ab84f1f0fea121be1abf1179c5926"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.1.1"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "103e017ff186ac39d731904045781c9bacfca2b0"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.1.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport"]
git-tree-sha1 = "9f8f52aad2399d7714b400ff9d203254b0a89c4a"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.1.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "0f81a77ede3da0dc714ea61e81c76b25db4ab87a"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.1.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "d4bb32e09d6b68ce2eb45fb81001eab46f60717a"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.2.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "590561f3af623d5485d070b4d7044f8854535f5a"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.2.1"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "e4be6539f4aaae8db1f29fcfdf6ef817df1f25cf"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.2.2"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "ef44754f10e0dfb9bb55ded382afed44cd94ab57"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.1.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "StaticArrays"]
git-tree-sha1 = "a8b7f8107c477e07c6a6c00d1d66cac68b801bbc"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.1.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport"]
git-tree-sha1 = "da525d277962a1b76102c79f30cb0c31e13fe5b9"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.1.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "332f9d17d0229218f66a73492162267359ba85e9"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.1.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "41c09d9c20877546490f907d8dffdd52690dd65f"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.1.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "96b47cdd12cb4ce8f70d701b49f855271a462bd4"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.2.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "f6683803a58de600ab7a26d2f49411c9923e9721"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.1.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "7dbe4ac56f930df5e9abd003cedb54e25cbbea86"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.2.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "348fd6def9a88518715425025eadd58517017325"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.1.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "1b0d894c880e25f7d0b022d7257638cf8ce5b311"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.1.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "4e8b8c8b81df3df17e2eb4603115db3b30a88235"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.1.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "96552f7d4619fabab4038a29ed37dd55e9eb513a"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.1.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "81d7841e73e385b9925d5c8e4427f2adcdda55db"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.1.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6d38fea02d983051776a856b7df75b30cf9a3c1f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.16"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "1a9cfb2dc2c2f1bd63f1906d72af39a79b49b736"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.11"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "6c62ce45f268f3f958821a1e5192cf91c75ae89c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.24"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "6f4dca5fd8e97087a76b7ab8384d1c3086ace0b7"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.3"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "6db1a75507051bc18bfa131fbc7c3f169cc4b2f6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.23"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RootedTrees]]
deps = ["LaTeXStrings", "Latexify", "LinearAlgebra", "Preferences", "RecipesBase", "Requires"]
git-tree-sha1 = "c0c464d3063e46e4128d21fd677ca575ace44fdc"
uuid = "47965b36-3f3e-11e9-0dcf-4570dfd42a8c"
version = "2.23.1"
weakdeps = ["Plots"]

    [deps.RootedTrees.extensions]
    PlotsExt = "Plots"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "f102316e5c958b425ef530ee51c7c8a1def55d1f"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.58.1"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "f66048bb969e67bd7d1bdd03cd0b81219642bbd0"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.1"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "6149620767866d4b0f0f7028639b6e661b6a1e44"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.12"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "25514a6f200219cd1073e4ff23a6324e4a7efe64"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "44021f3efc023be3871195d8ad98b865001a2fa1"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.12.3"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "6914df6005bab9940e2a96879a97a43e1fb1ce78"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "0.6.8"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerDataInterpolationsExt = "DataInterpolations"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    DataInterpolations = "82cc6244-b520-54b8-b5a6-8a565e85f1d0"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "b906758c107b049b6b71599b9f928d9b14e5554a"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.23.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DataStructures", "DocStringExtensions", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "f37f046636f8dc353a39279abfefe296db212171"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.8"
weakdeps = ["Colors"]

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsColorsExt = "Colors"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "87d51a3ee9a4b0d2fe054bdd3fc2436258db2603"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.1.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "777657803913ffc7e8cc20f0fd04b634f871af8f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.8"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "920acf6ae36c86f23969fea1d317e040dbfccf53"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.4.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "bf4bad73c80e058b1d53788ff520e10c8bad7c9d"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.70.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "e87efb31e5360cb223a151c2398903dc2faeb32b"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.26.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "91db7ed92c66f81435fe880947171f1212936b14"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.3+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "20cf607cafb31f922bce84d60379203e7a126911"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.34"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "04e9157537ba51dad58336976f8d04b9ab7122f0"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.7.2"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IfElse", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "0caef7687abf7094132fa3112bf5514c36a99226"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.18.3"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestItems]]
git-tree-sha1 = "42fd9023fef18b9b78c8343a4e2f3813ffbcefcb"
uuid = "1c621080-faea-4a02-84b6-bbd5e436b8fe"
version = "1.0.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3a6f063d690135f5c1ba351412c82bae4d1402bf"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.25"

[[deps.Tokenize]]
git-tree-sha1 = "468b4685af4abe0e9fd4d7bf495a6554a6276e75"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.29"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "be986ad9dac14888ba338c2554dcfec6939e1393"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.1"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "4ab62a49f1d8d9548a1c8d1a75e5f55cf196f64e"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.71"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XLSX]]
deps = ["Artifacts", "Dates", "EzXML", "Printf", "Tables", "ZipArchives", "ZipFile"]
git-tree-sha1 = "7fca49e6dbb35b7b7471956c2a9d3d921360a00f"
uuid = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"
version = "0.10.4"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "6a451c6f33a176150f315726eba8b92fbfdb9ae7"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.4+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "15e637a697345f6743674f1322beefbc5dcd5cfc"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.ZipArchives]]
deps = ["ArgCheck", "CodecInflate64", "CodecZlib", "InputBuffers", "PrecompileTools", "TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "f7fd93aa03f519c25b8b328693f2d36ce01220a9"
uuid = "49080126-0e18-4c2a-b176-c102e4b3760c"
version = "2.4.0"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "f492b7fe1698e623024e873244f10d89c95c340a"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.10.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═67a68b30-9d7f-11ef-1c30-bf3107341c1b
# ╠═0583b4d2-6c54-4621-8d98-3df7a69e306d
# ╠═447af9d0-daa7-4855-be20-c38dac16e8a0
# ╠═ee3c7f42-4e1d-4ee5-b0b2-ffbe8b4fd197
# ╠═e9271d60-9abb-49b2-90ea-d5c543834e8e
# ╠═34816d55-a8a9-411f-9a03-adeceb6518f9
# ╠═fac108bf-63be-4524-8622-71ce4cddab60
# ╠═a798d3ac-5841-4811-9009-776f6f02b7de
# ╠═d8e639a5-c0cb-476c-b5a4-e3e841bcf4a0
# ╠═838f5859-0e71-42f6-97b7-6f00732438a2
# ╠═c31a50c7-b116-485a-aa4c-57bc37ab3812
# ╠═6c80bb29-dda0-4f1a-8865-47bb5f2bab72
# ╠═ef4446ea-5210-4586-a7c2-3ab21a3594f0
# ╠═94bfa7bc-458c-4fc1-b31d-64261d8fb6a6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
