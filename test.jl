using CaMKIIModel: μM, ms, Hz
using Catalyst

rn = @reaction_network begin

    @equations begin
        LCCa_PKAp ~ LCCap / LCCtot
        LCCb_PKAp ~ LCCbp / LCCtot
        fracPLBp ~ PLBp / PLBtot
        TnI_PKAp ~ TnIp / TnItot
        IKUR_PKAp ~ KURp / IKurtot
    end

    # Beta adrenergic receptor and G protein
    ($ISO * kf_LR, kr_LR), b1AR <--> LR
    (kf_LRG, kr_LRG), LR + Gs <--> LRG
    (kf_RG, kr_RG), b1AR + Gs <--> RG
    k_G_act, (RG, LRG) --> b1AR + GsaGTP + Gsby
    k_G_hyd, GsaGTP --> GsaGDP
    k_G_reassoc, GsaGDP + Gsby --> Gs
    # Ligand-mediated inactivation
    (kf_bARK, kr_bARK), LR <--> b1AR_S464
    kf_bARK, LRG --> b1AR_S464 + Gs
    # PKA-mediated receptor inactivation
    kr_PKA, b1AR <-- b1AR_S301
    kf_PKA * PKACI, (LR, b1AR) --> b1AR_S301
    kf_PKA * PKACI, (LRG, RG) --> b1AR_S301 + Gs
    # Adenylate cyclase and PDE
    (kf_AC_Gsa, kr_AC_Gsa), AC + GsaGTP <--> AC_GsaGTP
    k_G_hyd, AC_GsaGTP --> AC + GsaGDP
    (k_PKA_PDE * PKACII, k_PP_PDE), PDE <--> PDEp
    mm($ATP, k_AC_Gsa * AC_GsaGTP, Km_AC_Gsa) + mm($ATP, k_AC_basal * AC, Km_AC_basal), 0 --> cAMP
    mm(cAMP, k_cAMP_PDE * PDE + k_cAMP_PDEp * PDEp, Km_PDE_cAMP), cAMP => 0
    # cAMP activating PKA
    (kf_RC_cAMP, kr_RC_cAMP), RC_I + cAMP <--> RCcAMP_I
    (kf_RC_cAMP, kr_RC_cAMP), RC_II + cAMP <--> RCcAMP_II
    (kf_RCcAMP_cAMP, kr_RCcAMP_cAMP), RCcAMP_I + cAMP <--> RCcAMPcAMP_I
    (kf_RCcAMP_cAMP, kr_RCcAMP_cAMP), RCcAMP_II + cAMP <--> RCcAMPcAMP_II
    (kf_RcAMPcAMP_C, kr_RcAMPcAMP_C), RCcAMPcAMP_I <--> RcAMPcAMP_I + PKACI
    (kf_RcAMPcAMP_C, kr_RcAMPcAMP_C), RCcAMPcAMP_II <--> RcAMPcAMP_II + PKACII
    (kf_PKA_PKI, kr_PKA_PKI), PKACI + PKI <--> PKACI_PKI
    (kf_PKA_PKI, kr_PKA_PKI), PKACII + PKI <--> PKACII_PKI
    # PKA modifications
    mm(I1, k_PKA_I1 * PKACI, Km_PKA_I1), I1 => I1p
    mm(I1p, Vmax_PP2A_I1, Km_PP2A_I1), I1p => I1
    (kf_PP1_I1, kr_PP1_I1), I1p + PP1 <--> I1p_PP1
    mm(PLB, k_PKA_PLB * PKACI, Km_PKA_PLB), PLB => PLBp
    mm(PLBp, k_PP1_PLB * PP1, Km_PP1_PLB), PLBp => PLB
    mm(PLM, k_PKA_PLM * PKACI, Km_PKA_PLM), PLM => PLMp
    mm(PLMp, k_PP1_PLM * PP1, Km_PP1_PLM), PLMp => PLM
    mm(TnI, k_PKA_TnI * PKACI, Km_PKA_TnI), TnI => TnIp
    mm(TnIp, k_PP2A_TnI * PP2A_TnI, Km_PP2A_TnI), TnIp => TnI
    mm(LCCa, k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII, Km_PKA_LCC / epsilon), LCCa => LCCap
    mm(LCCap, k_PP2A_LCC * PP2A_LCC, Km_PP2A_LCC / epsilon), LCCap => LCCa
    mm(LCCb, k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII, Km_PKA_LCC / epsilon), LCCb => LCCbp
    mm(LCCbp, k_PP1_LCC * PP1_LCC, Km_PP1_LCC / epsilon), LCCbp => LCCb
    mm(KURn, k_pka_KUR*(PKAII_KURtot / PKAIItot) * PKACII, Km_pka_KUR / epsilon), KURn => KURp
    mm(KURp, PP1_KURtot * k_pp1_KUR, Km_pp1_KUR / epsilon), KURp => KURn
end

return rn
