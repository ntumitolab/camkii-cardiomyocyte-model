# Beta adrenergic system
using Catalyst
using ModelingToolkit

function build_bar_sys(;
    b1ARtot = 0.00528,
    Gstot = 3.83,
    ATP = 5e3,
    PDEtot = 22.85e-3,
    PKItot = 0.18,
    I1tot = 0.3,
    PLBtot = 106,
    PLMtot = 48,
    TnItot = 70,
    LCCtot = 0.025,
    PKACII_LCCtot = 0.025,
    PKAIItot = 0.059,
    PKAII_KURtot = 0.025,
    IKurtot = 0.025,
    PP1_KURtot = 0.025,      # [uM]
    PP1tot = 0.89,
    remove_conserved=true
)
    @variables t ISO(t)
    rn = @reaction_network begin
        # Beta adrenergic receptor
        ($ISO * kf_LR, kr_LR), b1AR <--> LR  # Ligand-receptor
        (kf_LRG, kr_LRG), LR + Gs <--> LRG   # G protein association
        (kf_RG, kr_RG), b1AR + Gs <--> RG    # G protein association
        k_G_act, RG --> b1AR + GsaGTP + Gsby
        k_G_act, LRG --> b1AR + GsaGTP + Gsby
        k_G_hyd, GsaGTP --> GsaGDP
        k_G_reassoc, GsaGDP + Gsby --> Gs
        kf_bARK, LR --> b1AR_S464
        kf_bARK, LRG --> b1AR_S464 + Gs
        kr_bARK, b1AR_S464 --> b1AR
        (kf_PKA * PKACI, kr_PKA), b1AR <--> b1AR_S301
        kf_PKA * PKACI, LR --> b1AR_S301
        kf_PKA * PKACI, LRG --> b1AR_S301 + Gs
        kf_PKA * PKACI, RG --> b1AR_S301 + Gs
        # Adenylyl cyclase
        (kf_AC_Gsa, kr_AC_Gsa), AC + GsaGTP <--> AC_GsaGTP
        k_G_hyd, AC_GsaGTP --> AC + GsaGDP
        (k_PKA_PDE * PKACII, k_PP_PDE), PDE <--> PDEp
        mm(ATP, k_AC_Gsa * AC_GsaGTP, Km_AC_Gsa) + mm(ATP, k_AC_basal * AC, Km_AC_basal), 0 --> cAMP
        mm(cAMP, k_cAMP_PDE * PDE + k_cAMP_PDEp * PDEp, Km_PDE_cAMP), cAMP --> 0
        (kf_RC_cAMP, kr_RC_cAMP), RC_I + cAMP <--> RCcAMP_I
        (kf_RCcAMP_cAMP, kr_RCcAMP_cAMP), RCcAMP_I + cAMP <--> RCcAMPcAMP_I
        (kf_RcAMPcAMP_C, kr_RcAMPcAMP_C), RCcAMPcAMP_I <--> RcAMPcAMP_I + PKACI
        (kf_PKA_PKI, kr_PKA_PKI), PKACI + PKI <--> PKACI_PKI
        (kf_RC_cAMP, kr_RC_cAMP), RC_II + cAMP <--> RCcAMP_II
        (kf_RCcAMP_cAMP, kr_RCcAMP_cAMP), RCcAMP_II + cAMP <--> RCcAMPcAMP_II
        (kf_RcAMPcAMP_C, kr_RcAMPcAMP_C), RCcAMPcAMP_II <--> RcAMPcAMP_II + PKACII
        (kf_PKA_PKI, kr_PKA_PKI), PKACII + PKI <--> PKACII_PKI

        # I-1/PP1/PLB/TnI module
        mm(I1, k_PKA_I1 * PKACI, Km_PKA_I1), I1 => I1p
        mm(I1p, Vmax_PP2A_I1, Km_PP2A_I1), I1p => I1
        (kf_PP1_I1, kr_PP1_I1), I1p + PP1 <--> I1p_PP1
        mm(PLB, k_PKA_PLB * PKACI, Km_PKA_PLB), PLB => PLBp
        mm(PLBp, k_PP1_PLB * PP1, Km_PP1_PLB), PLBp => PLB
        mm(PLM, k_PKA_PLM * PKACI, Km_PKA_PLM), PLM => PLMp
        mm(PLMp, k_PP1_PLM * PP1, Km_PP1_PLM), PLMp => PLM
        mm(TnI, k_PKA_TnI * PKACI, Km_PKA_TnI), TnI => TnIp
        mm(TnIp, PP2A_TnI * k_PP2A_TnI, Km_PP2A_TnI), TnIp => TnI

        # LCC modification
        mm(LCCa, k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII, Km_PKA_LCC/epsilon), LCCa => LCCap
        mm(LCCap, k_PP2A_LCC * PP2A_LCC, Km_PP2A_LCC/epsilon), LCCap => LCCa
        mm(LCCb, k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII, Km_PKA_LCC/epsilon), LCCb => LCCbp
        mm(LCCbp, k_PP1_LCC * PP1_LCC, Km_PP1_LCC/ epsilon), LCCbp => LCCb

        # Ikur modification
        mm(KURn, k_pka_KUR * (PKAII_KURtot / PKAIItot) * PKACII, Km_pka_KUR/epsilon), KURn => KURp
        mm(KURp, PP1_KURtot * k_pp1_KUR, Km_pp1_KUR/epsilon), KURp => KURn
    end

    setdefaults!(rn, [
        :kf_LR => 1,                    # (1/[uM ms]) forward rate for ISO binding to b1AR
        :kr_LR => 0.285,                # (1/ms) reverse rate for ISO binding to b1AR
        :kf_LRG => 1,                   # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
        :kr_LRG => 0.062,               # (1/ms) reverse rate for ISO:b1AR association with Gs
        :kf_RG => 1,                    # (1/[uM ms]) forward rate for b1AR association with Gs
        :kr_RG => 33,                   # (1/ms) reverse rate for b1AR association with Gs
        :k_G_act => 16e-3,              # (1/ms) rate constant for Gs activation
        :k_G_hyd => 0.8e-3,             # (1/ms) rate constant for G-protein hydrolysis
        :k_G_reassoc => 1.21,           # (1/[uM ms]) rate constant for G-protein reassociation
        :kf_bARK => 1.1e-6,             # (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
        :kr_bARK => 2.2e-6,             # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
        :kf_PKA => 3.6e-6,              # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
        :kr_PKA => 2.2e-6,              # (1/ms) reverse rate for b1AR phosphorylation by PKA
        :k_AC_basal => 0.2e-3,          # (1/ms) basal cAMP generation rate by AC
        :Km_AC_basal => 1.03e3,         # (uM) basal AC affinity for ATP
        :kr_AC_Gsa => 0.4,              # (uM) AC dissociation with Gsa
        :kf_AC_Gsa => 1,                # (1/[uM ms]) forward rate for AC association with Gsa
        :k_AC_Gsa => 8.5e-3,            # (1/ms) basal cAMP generation rate by AC:Gsa
        :Km_AC_Gsa => 315.0,            # (uM) AC:Gsa affinity for ATP
        :k_cAMP_PDE => 5e-3,            # (1/ms) cAMP hydrolysis rate by PDE
        :k_cAMP_PDEp => 10e-3,          # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
        :Km_PDE_cAMP => 1.3,            # (uM) PDE affinity for cAMP
        :k_PKA_PDE => 7.5e-3,           # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
        :k_PP_PDE => 1.5e-3,            # (1/ms) rate constant for PDE dephosphorylation by phosphatases
        :kf_RC_cAMP => 1,               # (1/[uM ms]) Kd for PKA RC binding to cAMP
        :kf_RCcAMP_cAMP => 1,           # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
        :kf_RcAMPcAMP_C => 4.375,       # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
        :kf_PKA_PKI => 1,               # (1/[uM ms]) Ki for PKA inhibition by PKI
        :kr_RC_cAMP => 1.64,            # (1/ms) Kd for PKA RC binding to cAMP
        :kr_RCcAMP_cAMP => 9.14,        # (1/ms) Kd for PKA RC:cAMP binding to cAMP
        :kr_RcAMPcAMP_C => 1,           # (1/ms) Kd for PKA R:cAMPcAMP binding to C
        :kr_PKA_PKI => 2e-4,            # (1/ms) Ki for PKA inhibition by PKI
        :epsilon => 10,                 # (-) AKAP-mediated scaling factor
        :k_PKA_I1 => 60e-3,              # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
        :Km_PKA_I1 => 1.0,               # (uM) Km for I-1 phosphorylation by type 1 PKA
        :Vmax_PP2A_I1 => 14.0e-3,        # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
        :Km_PP2A_I1 => 1.0,              # (uM) Km for I-1 dephosphorylation by PP2A
        :kr_PP1_I1 => 1.0e-3,            # (uM) Ki for PP1 inhibition by I-1
        :kf_PP1_I1 => 1,                 # (uM) kf for PP1 inhibition by I-1
        :k_PKA_PLB => 54e-3,   # [1/ms]
        :Km_PKA_PLB => 21,     # [uM]
        :k_PP1_PLB => 8.5e-3,  # [1/ms]
        :Km_PP1_PLB => 7.0,    # [uM]
        :k_PKA_PLM => 54e-3,   # [1/ms]
        :Km_PKA_PLM => 21,     # [uM]
        :k_PP1_PLM => 8.5e-3,  # [1/ms]
        :Km_PP1_PLM => 7.0,    # [uM]
        :PP1_LCC => 0.025,         # [uM]
        :PP2A_LCC => 0.025,        # [uM]
        :k_PKA_LCC => 54e-3,       # [1/ms]
        :Km_PKA_LCC => 21,         # [uM]
        :k_PP1_LCC => 8.52e-3,     # [1/ms] RABBIT, MOUSE
        :Km_PP1_LCC => 3,          # [uM]
        :k_PP2A_LCC => 10.1e-3,    # [1/ms]
        :Km_PP2A_LCC => 3,         # [uM]
        :PP2A_TnI => 0.67,         # [uM]
        :k_PKA_TnI => 54e-3,       # [1/ms]
        :Km_PKA_TnI => 21,         # [uM]
        :k_PP2A_TnI => 10.1e-3,    # [1/ms]
        :Km_PP2A_TnI => 4.1,       # [uM]
        :k_pka_KUR => 54e-3,       # [1/ms]
        :Km_pka_KUR => 21,         # [uM]
        :k_pp1_KUR => 8.52e-3,     # [1/ms]
        :Km_pp1_KUR => 7,          # [uM]
    ])

    return convert(ODESystem, rn; remove_conserved)
end
