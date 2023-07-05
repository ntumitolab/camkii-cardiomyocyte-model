#=
    This module describes the beta-adrenergic signaling pathway in mouse
    ventricular myocyte, and this file was built upon the code developeded
    by Yang and Saucerman.
    Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
    feed-forward regulator of Ca2+ in beta-adrenergic signaling,
    accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular
    Cardiology 52, 1048-1055.
=#

using ModelingToolkit

function get_bar_equations()
    @variables t
    D = Differential(t)

    # Drug concentrations
    @parameters (ISO=0, FSK=0, IBMX=0)

    @parameters b1ARtot = 0.00528 # (uM) total b1-AR protein (MOUSE); RABBIT=0.028
    @parameters Gstot = 3.83    # (uM) total Gs protein

    @parameters kf_LR = 1       # (1/[uM ms]) forward rate for ISO binding to b1AR
    @parameters kr_LR = 0.285   # (1/ms) reverse rate for ISO binding to b1AR
    @parameters kf_LRG = 1       # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
    @parameters kr_LRG = 0.062   # (1/ms) reverse rate for ISO:b1AR association with Gs
    @parameters kf_RG = 1       # (1/[uM ms]) forward rate for b1AR association with Gs
    @parameters kr_RG = 33      # (1/ms) reverse rate for b1AR association with Gs
    @parameters k_G_act = 16e-3          # (1/ms) rate constant for Gs activation
    @parameters k_G_hyd = 0.8e-3         # (1/ms) rate constant for G-protein hydrolysis
    @parameters k_G_reassoc = 1.21           # (1/[uM ms]) rate constant for G-protein reassociation
    @parameters kf_bARK = 1.1e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
    @parameters kr_bARK = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
    @parameters kf_PKA = 3.6e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
    @parameters kr_PKA = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by PKA

    @variables b1AR(t) b1ARact(t) b1AR_S464(t) b1AR_S301(t) LR(t) LRG(t) RG(t) Gs(t) Gsby(t) bARK_desens(t) bARK_resens(t) PKA_desens(t) PKA_resens(t) G_act(t) GsaGTPtot(t) GsaGDP(t)

    @parameters ACtot = 70.57e-3       # (uM) total adenylyl cyclase # MOUSE
    # ACtot=47e-3; # RABBIT
    @parameters ATP = 5e3            # (uM) total ATP
    @parameters k_AC_basal = 0.2e-3         # (1/ms) basal cAMP generation rate by AC
    @parameters Km_AC_basal = 1.03e3         # (uM) basal AC affinity for ATP

    @parameters Kd_AC_Gsa = 0.4            # (uM) Kd for AC association with Gsa
    @parameters kf_AC_Gsa = 1              # (1/[uM ms]) forward rate for AC association with Gsa
    @parameters kr_AC_Gsa = Kd_AC_Gsa      # (1/ms) reverse rate for AC association with Gsa

    @parameters k_AC_Gsa = 8.5e-3         # (1/ms) basal cAMP generation rate by AC:Gsa
    @parameters Km_AC_Gsa = 315.0          # (uM) AC:Gsa affinity for ATP

    @parameters Kd_AC_FSK = 44.0           # (uM) Kd for FSK binding to AC
    @parameters k_AC_FSK = 7.3e-3         # (1/ms) basal cAMP generation rate by AC:FSK
    @parameters Km_AC_FSK = 860.0          # (uM) AC:FSK affinity for ATP
    @parameters PDEtot = 22.85e-3       # (uM) total phosphodiesterase
    @parameters k_cAMP_PDE = 5e-3           # (1/ms) cAMP hydrolysis rate by PDE
    @parameters k_cAMP_PDEp = 2 * k_cAMP_PDE   # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
    @parameters Km_PDE_cAMP = 1.3            # (uM) PDE affinity for cAMP
    @parameters Kd_PDE_IBMX = 30.0           # (uM) Kd_R2cAMP_C for IBMX binding to PDE
    @parameters k_PKA_PDE = 7.5e-3         # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
    @parameters k_PP_PDE = 1.5e-3         # (1/ms) rate constant for PDE dephosphorylation by phosphatases


    @variables cAMP(t) AC(t) GsaGTP(t) AC_FSK(t) AC_ACT_BASAL(t) AC_ACT_GSA(t) AC_ACT_FSK(t) RC_I(t) RCcAMP_I(t) RCcAMPcAMP_I(t) RcAMPcAMP_I(t) PKACI(t) PKACI_PKI(t) RC_II(t) RCcAMP_II(t) RCcAMPcAMP_II(t) RCcAMPcAMP_II(t) RcAMPcAMP_II(t) PKACII(t) PKACII_PKI(t) AC_GsaGTP(t) PDE_IBMX(t) PDE_ACT(t) PDEp(t) PDEtot(t) cAMPtot(t)

    # PKA module
    @parameters PKItot          = 0.18           # (uM) total PKI
    @parameters kf_RC_cAMP      = 1              # (1/[uM ms]) Kd for PKA RC binding to cAMP
    @parameters kf_RCcAMP_cAMP  = 1              # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
    @parameters kf_RcAMPcAMP_C  = 4.375          # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
    @parameters kf_PKA_PKI      = 1              # (1/[uM ms]) Ki for PKA inhibition by PKI
    @parameters kr_RC_cAMP      = 1.64           # (1/ms) Kd for PKA RC binding to cAMP
    @parameters kr_RCcAMP_cAMP  = 9.14           # (1/ms) Kd for PKA RC:cAMP binding to cAMP
    @parameters kr_RcAMPcAMP_C  = 1              # (1/ms) Kd for PKA R:cAMPcAMP binding to C
    @parameters kr_PKA_PKI      = 2e-4           # (1/ms) Ki for PKA inhibition by PKI
    @parameters epsilon         = 10             # (-) AKAP-mediated scaling factor

    @variables PKI(t) PKACI_PKI(t) PKACII_PKI(t)


     # PP1tot = pin(8) # PP1tot = 0.89  # (uM) total phosphatase 1
     @parameters PP1tot          = 0.89           # (uM) total phosphatase 1
     @parameters I1tot           = 0.3            # (uM) total inhibitor 1
     @parameters k_PKA_I1        = 60e-3          # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
     @parameters Km_PKA_I1       = 1.0            # (uM) Km for I-1 phosphorylation by type 1 PKA
     @parameters Vmax_PP2A_I1    = 14.0e-3        # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
     @parameters Km_PP2A_I1      = 1.0            # (uM) Km for I-1 dephosphorylation by PP2A

     @parameters Ki_PP1_I1       = 1.0e-3         # (uM) Ki for PP1 inhibition by I-1
     @parameters kf_PP1_I1       = 1              # (uM) Ki for PP1 inhibition by I-1
     @parameters kr_PP1_I1       = Ki_PP1_I1      # (uM) Ki for PP1 inhibition by I-1

     @variables I1(t) PP1(t) I1p(t) I1_phosph(t) I1_dephosph(t) I1ptot(t) I1p_PP1(t)

    bar_eqs = [
        # beta-adrenergic receptor
        b1ARtot ~ b1ARact + b1AR_S464 + b1AR_S301,
        b1ARact ~ b1AR + LR + LRG + RG,
        Gstot ~ Gs + LRG + RG + Gsby,
        D(LR) ~ kf_LR * ISO * b1AR - kr_LR * LR + kr_LRG * LRG - kf_LRG * LR * Gs,
        D(LRG) ~ kf_LRG * LR * Gs - kr_LRG * LRG - k_G_act * LRG,
        D(RG) ~ kf_RG * b1AR * Gs - kr_RG * RG - k_G_act * RG,
        bARK_desens ~ kf_bARK * (LR + LRG),
        bARK_resens ~ kr_bARK * b1AR_S464,
        PKA_desens ~ kf_PKA * PKACI * b1ARact,
        PKA_resens ~ kr_PKA * b1AR_S301,
        D(b1AR_S464) ~ bARK_desens - bARK_resens,
        D(b1AR_S301) ~ PKA_desens - PKA_resens,
        G_act ~ k_G_act * (RG + LRG),
        G_hyd ~ k_G_hyd * GsaGTPtot,
        G_reasso ~ k_G_reassoc * GsaGDP * Gsby,
        D(GsaGTPtot) ~ G_act - G_hyd,
        D(GsaGDP) ~ G_hyd - G_reassoc,
        D(Gsby) ~ G_act - G_reassoc,
    ]

    camp_eqs = [
        cAMPtot ~ cAMP + RCcAMP_I + 2 * RCcAMPcAMP_I + 2 * RcAMPcAMP_I + RCcAMP_II + 2 * RCcAMPcAMP_II + 2 * RcAMPcAMP_II,
        ACtot ~ AC + AC_GsaGTP,
        GsaGTPtot ~ GsaGTP + AC_GsaGTP,
        D(AC_GsaGTP) ~ kf_AC_Gsa * GsaGTP * AC - kr_AC_Gsa * AC_GsaGTP,
        AC_FSK ~ FSK * AC / Kd_AC_FSK,
        AC_ACT_BASAL ~ k_AC_basal * AC * hil(ATP, Km_AC_basal),
        AC_ACT_GSA ~ k_AC_Gsa * AC_GsaGTP * hil(ATP, Km_AC_Gsa),
        AC_ACT_FSK ~ k_AC_FSK * AC_FSK * hil(ATP, Km_AC_FSK),
        PDE_IBMX ~ PDEtot * IBMX / Kd_PDE_IBMX,
        PDEtot ~ PDE + PDE_IBMX + PDEp,
        D(PDEp) ~ k_PKA_PDE * PKACII * PDE - k_PP_PDE * PDEp,
        PDE_ACT ~ k_cAMP_PDE * PDE * hil(cAMP, Km_PDE_cAMP) + k_cAMP_PDEp * PDEp * hil(cAMP, Km_PDE_cAMP),
        D(cAMPtot) ~ AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT,
        PKI ~ PKItot - PKACI_PKI - PKACII_PKI,
        D(RC_I) ~ - kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I,
        D(RCcAMP_I) ~ - kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I,
        D(RCcAMPcAMP_I) ~ - kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI,
        D(RcAMPcAMP_I) ~ - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I,
        D(PKACI) ~ - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI,
        D(PKACI_PKI) ~ - kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI,
        D(RC_II) ~ - kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II,
        D(RCcAMP_II) ~ - kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II,
        D(RCcAMPcAMP_II) ~ - kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII,
        D(RcAMPcAMP_II) ~ - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II,
        D(PKACII) ~ - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI,
        D(PKACII_PKI) ~ - kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI
    ]

    # PP1
    pp1_eqs = [
        I1 ~ I1tot - I1ptot,
        PP1 ~ PP1tot - I1p_PP1,
        I1p ~ I1ptot - I1p_PP1,
        I1_phosph ~ k_PKA_I1*PKACI*hil(I1,Km_PKA_I1),
        I1_dephosph ~ Vmax_PP2A_I1*hil(I1ptot,Km_PP2A_I1),
        D(I1p_PP1) ~ kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1,
        D(I1ptot) ~ I1_phosph - I1_dephosph
    ]

    return vcat(bar_eqs, camp_eqs, pp1_eqs)
end
