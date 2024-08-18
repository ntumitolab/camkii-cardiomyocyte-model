# Beta adrenergic system activated by isoproterenol
using Catalyst
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

function get_bar_sys(;
    ATP=5000μM,
    ISO=0μM,
    b1ARtot=0.00528μM,
    Gstot=3.83μM,
    PDEtot = 22.85e-3μM,
    PKItot = 0.18μM,
    I1tot = 0.3μM,
    PLBtot = 106μM,
    PLMtot = 48μM,
    TnItot = 70μM,
    LCCtot = 0.025μM,
    IKurtot = 0.025μM,
    PP1tot = 0.89μM,
    RItot = 1.18μM,
    RIItot = 0.118μM,
    ACtot = 70.57e-3μM,
    PKACII_LCCtot = 0.025μM,
    PKAIItot = 0.059μM,
    PKAII_KURtot = 0.025μM,
    PP1_KURtot = 0.025μM,
    remove_conserved=true,
)
    rn = @reaction_network begin
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
        mm(KURn, k_pka_KUR * (PKAII_KURtot / PKAIItot) * PKACII, Km_pka_KUR / epsilon), KURn => KURp
        mm(KURp, PP1_KURtot * k_pp1_KUR, Km_pp1_KUR / epsilon), KURp => KURn
    end

    setdefaults!(rn, [
        :kf_LR => 1 / (μM * ms),                    # forward rate for ISO binding to b1AR
        :kr_LR => 285Hz,                            # reverse rate for ISO binding to b1AR
        :kf_LRG => 1 / (μM * ms),                   # forward rate for ISO:b1AR association with Gs
        :kr_LRG => 62Hz,                            # reverse rate for ISO:b1AR association with Gs
        :kf_RG => 1 / (μM * ms),                    # forward rate for b1AR association with Gs
        :kr_RG => 33000Hz,                          # reverse rate for b1AR association with Gs
        :k_G_act => 16Hz,                           # rate constant for Gs activation
        :k_G_hyd => 0.8Hz,                          # rate constant for G-protein hydrolysis
        :k_G_reassoc => 1.21 / (μM * ms),           # rate constant for G-protein reassociation
        :kf_bARK => 1.1Hz / mM,                     # forward rate for b1AR phosphorylation by b1ARK
        :kr_bARK => 2.2e-3Hz,                       # reverse rate for b1AR phosphorylation by b1ARK
        :kf_PKA => 3.6Hz / mM,                      # forward rate for b1AR phosphorylation by PKA
        :kr_PKA => 2.2e-3Hz,                        # reverse rate for b1AR phosphorylation by PKA
        :k_AC_basal => 0.2Hz,                       # basal cAMP generation rate by AC
        :Km_AC_basal => 1.03mM,                     # basal AC affinity for ATP
        :kr_AC_Gsa => 0.4μM,                        # AC dissociation with Gsa
        :kf_AC_Gsa => 1 / (μM * ms),                # forward rate for AC association with Gsa
        :k_AC_Gsa => 8.5Hz,                         # basal cAMP generation rate by AC:Gsa
        :Km_AC_Gsa => 315.0μM,                      # AC:Gsa affinity for ATP
        :k_cAMP_PDE => 5Hz,                         # cAMP hydrolysis rate by PDE
        :k_cAMP_PDEp => 10Hz,                       # cAMP hydrolysis rate by phosphorylated PDE
        :Km_PDE_cAMP => 1.3μM,                      # PDE affinity for cAMP
        :k_PKA_PDE => 7.5Hz,                        # rate constant for PDE phosphorylation by type 1 PKA
        :k_PP_PDE => 1.5Hz,                         # rate constant for PDE dephosphorylation by phosphatases
        :kf_RC_cAMP => 1 / (μM * ms),               # Kd for PKA RC binding to cAMP
        :kf_RCcAMP_cAMP => 1 / (μM * ms),           # Kd for PKA RC:cAMP binding to cAMP
        :kf_RcAMPcAMP_C => 4.375 / (μM * ms),       # Kd for PKA R:cAMPcAMP binding to C
        :kf_PKA_PKI => 1 / (μM * ms),               # Ki for PKA inhibition by PKI
        :kr_RC_cAMP => 1.64 / ms,                   # rate constant for PKA RC unbinding to cAMP
        :kr_RCcAMP_cAMP => 9.14 / ms,               # rate constant for PKA RC:cAMP unbinding to cAMP
        :kr_RcAMPcAMP_C => 1 / ms,                  # rate constant for PKA R:cAMPcAMP unbinding to C
        :kr_PKA_PKI => 0.2Hz,                       # reverse rate for PKA inhibition by PKI
        :epsilon => 10,                             # AKAP-mediated scaling factor
        :k_PKA_I1 => 60Hz,                          # rate constant for I-1 phosphorylation by type 1 PKA
        :Km_PKA_I1 => 1.0μM,                        # Km for I-1 phosphorylation by type 1 PKA
        :Vmax_PP2A_I1 => 14Hz,                      # Vmax for I-1 dephosphorylation by PP2A
        :Km_PP2A_I1 => 1.0μM,                       # Km for I-1 dephosphorylation by PP2A
        :kr_PP1_I1 => 1.0Hz,                        # Ki for PP1 inhibition by I-1
        :kf_PP1_I1 => 1μM,                          # kf for PP1 inhibition by I-1
        :k_PKA_PLB => 54Hz,
        :Km_PKA_PLB => 21μM,
        :k_PP1_PLB => 8.5Hz,
        :Km_PP1_PLB => 7.0μM,
        :k_PKA_PLM => 54Hz,
        :Km_PKA_PLM => 21μM,
        :k_PP1_PLM => 8.5Hz,
        :Km_PP1_PLM => 7.0μM,
        :PP1_LCC => 0.025μM,
        :PP2A_LCC => 0.025μM,
        :k_PKA_LCC => 54Hz,
        :Km_PKA_LCC => 21μM,
        :k_PP1_LCC => 8.52Hz,
        :Km_PP1_LCC => 3μM,
        :k_PP2A_LCC => 10.1Hz,
        :Km_PP2A_LCC => 3μM,
        :PP2A_TnI => 0.67μM,
        :k_PKA_TnI => 54Hz,
        :Km_PKA_TnI => 21μM,
        :k_PP2A_TnI => 10.1Hz,
        :Km_PP2A_TnI => 4.1μM,
        :k_pka_KUR => 54Hz,
        :Km_pka_KUR => 21μM,
        :k_pp1_KUR => 8.52Hz,
        :Km_pp1_KUR => 7μM,
        :PKACII_LCCtot => 0.025μM,
        :PKAIItot => 0.059μM,
        :PKAII_KURtot => 0.025μM,
        :PP1_KURtot => 0.025μM,
    ])

    @unpack b1AR, LR, LRG, RG, b1AR_S464, b1AR_S301, Gs, Gsby, GsaGTP, GsaGDP, AC, AC_GsaGTP, PDE, PDEp, cAMP, RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKI, PKACI_PKI, PKACII_PKI, PKACII, RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, I1, I1p, I1p_PP1, PP1, PLB, PLBp, PLM, PLMp, LCCa, LCCap, LCCb, LCCbp, TnI, TnIp, KURn, KURp = rn

    setdefaults!(rn, [
        b1AR => b1ARtot - LR - LRG - RG - b1AR_S464 - b1AR_S301,
        LR => 6.0e-5μM,
        LRG => 0.00294μM,
        RG => 7.0e-5μM,
        b1AR_S464 => 0.00047μM,
        b1AR_S301 => 0.0011μM,
        Gs => Gstot - Gsby - LRG - RG,
        GsaGTP => 0.06028μM - 0.00814μM,
        GsaGDP => 0.00066μM,
        Gsby => 0.06071μM,
        AC => ACtot - AC_GsaGTP,
        AC_GsaGTP => 0.00814μM,
        PDE => PDEtot - PDEp,
        PDEp => 0.00589μM,
        cAMP => 1.50399μM,
        RC_I => RItot - RCcAMP_I - RCcAMPcAMP_I - RcAMPcAMP_I,
        RCcAMP_I => 0.28552μM,
        RCcAMPcAMP_I => 0.04698μM,
        RcAMPcAMP_I => 0.53564μM,
        PKACI => 0.38375μM,
        PKI => PKItot - PKACI_PKI - PKACII_PKI,
        PKACI_PKI => 0.15239μM,
        RC_II => RIItot - RCcAMP_II - RCcAMPcAMP_II - RcAMPcAMP_II,
        RCcAMP_II => 0.00934μM,
        RCcAMPcAMP_II => 0.00154μM,
        RcAMPcAMP_II => 0.09691μM,
        PKACII => 0.06938μM,
        PKACII_PKI => 0.02753μM,
        I1 => I1tot - I1p - I1p_PP1,
        PP1 => PP1tot - I1p_PP1,
        I1p => 0.00033μM,
        I1p_PP1 => 0.19135μM,
        PLB => PLBtot - PLBp,
        PLBp => 98.33936μM,
        PLM => PLMtot - PLMp,
        PLMp => 41.19479μM,
        LCCa => LCCtot - LCCap,
        LCCb => LCCtot - LCCbp,
        TnI => TnItot - TnIp,
        KURn => IKurtot - KURp,
        TnIp => 60.75646μM,
        LCCap => 0.01204μM,
        LCCbp => 0.01313μM,
        KURp => 0.01794μM,
    ])

    return rn
end
