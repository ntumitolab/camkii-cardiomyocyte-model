#=
    This module describes the beta-adrenergic signaling pathway in mouse
    ventricular myocyte, and this file was built upon the code developeded
    by Yang and Saucerman.
    Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
    feed-forward regulator of Ca2+ in beta-adrenergic signaling,
    accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular
    Cardiology 52, 1048-1055.
=#

using Catalyst
using ModelingToolkit

function make_b1ar_rn()
    # Input variables
    @variables t ISO(t) ATP(t)

    rn = @reaction_network begin
        # BAR
        (kf_LR * $ISO, kr_LR), b1AR <--> LR
        (kf_LRG, kr_LRG), LR + Gs <--> LRG
        (kf_RG, kr_RG), b1AR + Gs <--> RG
        k_G_act, (RG, LRG) --> b1AR + Gsby + GsaGTP
        k_G_hyd, GsaGTP --> GsaGDP
        k_G_reassoc, Gsby + GsaGDP --> Gs
        kf_bARK, LR --> b1AR_S464
        kf_bARK, LRG --> b1AR_S464 + Gs
        kr_bARK, b1AR_S464 --> b1AR
        # AC & PKA
        kf_PKA * PKACI, (b1AR, LR) --> b1AR_S301
        kf_PKA * PKACI, (RG, LRG) --> b1AR_S301 + Gs
        kr_PKA, b1AR_S301 --> b1AR
        (kf_AC_Gsa, kr_AC_Gsa), AC + GsaGTP <--> AC_GsaGTP
        (k_PKA_PDE * PKACII, k_PP_PDE), PDE <--> PDEp
        (mm($ATP, k_AC_basal * AC, Km_AC_basal) + mm($ATP, k_AC_Gsa * AC_GsaGTP, Km_AC_Gsa)), 0 --> cAMP
        mm(cAMP, k_cAMP_PDE * PDE + k_cAMP_PDEp * PDEp, Km_PDE_cAMP), cAMP --> 0
        # cAMP and PKAC
        (kf_RC_cAMP, kr_RC_cAMP), RC_I + cAMP <--> RCcAMP_I
        (kf_RC_cAMP, kr_RC_cAMP), RC_II + cAMP <--> RCcAMP_II
        (kf_RCcAMP_cAMP, kr_RCcAMP_cAMP), RCcAMP_I + cAMP <--> RCcAMPcAMP_I
        (kf_RCcAMP_cAMP, kr_RCcAMP_cAMP), RCcAMP_II + cAMP <--> RCcAMPcAMP_II
        (kf_RcAMPcAMP_C, kr_RcAMPcAMP_C), RCcAMPcAMP_II <--> RcAMPcAMP_I + PKACI
        (kf_RcAMPcAMP_C, kr_RcAMPcAMP_C), RCcAMPcAMP_II <--> RcAMPcAMP_II + PKACII
        (kf_PKA_PKI, kr_PKA_PKI), PKACI + PKI <--> PKACI_PKI
        (kf_PKA_PKI, kr_PKA_PKI), PKACII + PKI <--> PKACII_PKI
        # PKA downstream
        mm(I1, k_PKA_I1 * PKACI, Km_PKA_I1), I1 => I1p
        mm(I1p, Vmax_PP2A_I1, Km_PP2A_I1), I1 <= I1p
        (kf_PP1_I1, kr_PP1_I1), PP1 + I1p <--> I1p_PP1
        mm(PLB, k_PKA_PLB * PKACI, Km_PKA_PLB), PLB => PLBp
        mm(PLBp, k_PP1_PLB * PP1, Km_PP1_PLB), PLBp => PLB
        mm(PLM, k_PKA_PLM * PKACI, Km_PKA_PLM), PLM => PLMp
        mm(PLMp, k_PP1_PLM * PP1, Km_PP1_PLM), PLMp => PLMp
        mm(ϵ * LCCa, k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII, Km_PKA_LCC), LCCa => LCCap
        mm(ϵ * LCCap, k_PP2A_LCC * PP2A_LCC, Km_PP2A_LCC), LCCap => LCCa
        mm(ϵ * LCCb, k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII, Km_PKA_LCC), LCCb => LCCbp
        mm(ϵ * LCCbp, k_PP1_LCC * PP1_LCC, Km_PP1_LCC), LCCbp => LCCb
        mm(ϵ * RyR, kcat_pka_ryr * (PKAIIryrtot / PKAIItot) * PKACII, Km_pka_ryr), RyR => RyRp
        mm(ϵ * RyRp, kcat_pp1_ryr * PP1ryr, Km_pp1_ryr), RyRp => RyR
        mm(ϵ * RyRp, kcat_pp2a_ryr * PP2Aryr, Km_pp2a_ryr), RyRp => RyR
        mm(TnI, k_PKA_TnI * PKACI, Km_PKA_TnI), TnI => TnIp
        mm(TnIp, k_PP2A_TnI * PP2A_TnI, Km_PP2A_TnI), TnIp => TnI
        mm(ϵ * Kur, (PKAII_KURtot / PKAIItot) * PKACII, Km_pka_KUR), Kur => Kurp
        mm(ϵ * Kurp, PP1_KURtot * k_pp1_KUR, Km_pp1_KUR), Kurp => Kur
    end

    setdefaults!(rn, [
        :kf_LR => 1,            # (1/[uM ms]) forward rate for ISO binding to b1AR
        :kr_LR => 0.285,        # (1/ms) reverse rate for ISO binding to b1AR
        :kf_LRG => 1,           # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
        :kr_LRG => 0.062,       # (1/ms) reverse rate for ISO:b1AR association with Gs
        :kf_RG => 1,            # (1/[uM ms]) forward rate for b1AR association with Gs
        :kr_RG => 33,           # (1/ms) reverse rate for b1AR association with Gs
        :k_G_act => 16e-3,      # (1/ms) rate constant for Gs activation
        :k_G_hyd => 0.8e-3,     # (1/ms) rate constant for G-protein hydrolysis
        :k_G_reassoc => 1.21,   # (1/[uM ms]) rate constant for G-protein reassociation
        :kf_bARK => 1.1e-6,     # (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
        :kr_bARK => 2.2e-6,     # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
        :kf_PKA => 3.6e-6,      # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
        :kr_PKA => 2.2e-6,          # (1/ms) reverse rate for b1AR phosphorylation by PKA
        :k_AC_basal => 0.2e-3,     # (1/ms) basal cAMP generation rate by AC
        :Km_AC_basal => 1030.0,    # (uM) basal AC affinity for ATP
        :kf_AC_Gsa => 1,           # (1/[uM ms]) forward rate for AC association with Gsa
        :kr_AC_Gsa => 0.4,          # (1/ms) reverse rate for AC association with Gsa
        :k_AC_Gsa => 8.5e-3,        # (1/ms) basal cAMP generation rate by AC:Gsa
        :Km_AC_Gsa => 315.0,       # (uM) AC:Gsa affinity for ATP
        :k_cAMP_PDE => 5e-3,       # (1/ms) cAMP hydrolysis rate by PDE
        :k_cAMP_PDEp => 10e-3,     # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
        :Km_PDE_cAMP => 1.3,       # (uM) PDE affinity for cAMP
        :k_PKA_PDE => 7.5e-3,      # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
        :k_PP_PDE => 1.5e-3,       # (1/ms) rate constant for PDE dephosphorylation by phosphatases
        :kf_RC_cAMP => 1,          # (1/[uM ms]) Kd for PKA RC binding to cAMP
        :kf_RCcAMP_cAMP => 1,      # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
        :kf_RcAMPcAMP_C => 4.375,  # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
        :kf_PKA_PKI => 1,          # (1/[uM ms]) Ki for PKA inhibition by PKI
        :kr_RC_cAMP => 1.64,       # (1/ms) Kd for PKA RC binding to cAMP
        :kr_RCcAMP_cAMP => 9.14,   # (1/ms) Kd for PKA RC:cAMP binding to cAMP
        :kr_RcAMPcAMP_C => 1,      # (1/ms) Kd for PKA R:cAMPcAMP binding to C
        :kr_PKA_PKI => 2e-4,       # (1/ms) Ki for PKA inhibition by PKI
        :ϵ => 10,                  # (-) AKAP-mediated scaling factor
        :k_PKA_I1 => 60e-3,        # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
        :Km_PKA_I1 => 1.0,         # (uM) Km for I-1 phosphorylation by type 1 PKA
        :Vmax_PP2A_I1 => 14.0e-3,  # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
        :Km_PP2A_I1 => 1.0,        # (uM) Km for I-1 dephosphorylation by PP2A
        :kf_PP1_I1 => 1,           # (uM) Ki for PP1 inhibition by I-1
        :kr_PP1_I1 => 1.0e-3,      # (uM) Ki for PP1 inhibition by I-1
        :k_PKA_PLB => 54e-3,       # k_pka_plb     [1/ms]
        :Km_PKA_PLB => 21,         # Km_pka_plb    [uM]
        :k_PP1_PLB => 8.5e-3,      # k_pp1_plb     [1/ms]
        :Km_PP1_PLB => 7.0,        # Km_pp1_plb    [uM]
        :k_PKA_PLM => 54e-3,       # k_pka_plb     [1/ms]
        :Km_PKA_PLM => 21,         # Km_pka_plb    [uM]
        :k_PP1_PLM => 8.5e-3,      # k_pp1_plb     [1/ms]
        :Km_PP1_PLM => 7.0,        # Km_pp1_plb    [uM]
        :PKACII_LCCtot => 0.025,   # PKAIIlcctot   [uM]
        :PP1_LCC => 0.025,         # PP1lcctot     [uM]
        :PP2A_LCC => 0.025,        # PP2Alcctot    [uM]
        :k_PKA_LCC => 54e-3,       # k_pka_lcc     [1/ms]
        :Km_PKA_LCC => 21,         # Km_pka_lcc    [uM]
        :k_PP1_LCC => 8.52e-3,     # k_pp1_lcc     [1/ms] RABBIT, MOUSE
        :Km_PP1_LCC => 3,          # Km_pp1_lcc    [uM]
        :k_PP2A_LCC => 10.1e-3,    # k_pp2a_lcc    [1/ms]
        :Km_PP2A_LCC => 3,         # Km_pp2a_lcc   [uM]
        :PKAIIryrtot => 0.034,     # PKAIIryrtot   [uM]
        :PP1ryr => 0.034,          # PP1ryr        [uM]
        :PP2Aryr => 0.034,         # PP2Aryr       [uM]
        :kcat_pka_ryr => 54e-3,    # kcat_pka_ryr  [1/ms]
        :Km_pka_ryr => 21,         # Km_pka_ryr    [uM]
        :kcat_pp1_ryr => 8.52e-3,  # kcat_pp1_ryr  [1/ms]
        :Km_pp1_ryr => 7,          # Km_pp1_ryr    [uM]
        :kcat_pp2a_ryr => 10.1e-3, # kcat_pp2a_ryr [1/ms]
        :Km_pp2a_ryr => 4.1,       # Km_pp2a_ryr   [uM]
    ])

    return rn
end
