# Beta adrenergic system activated by isoproterenol
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

function get_bar_sys(;
    ATP=5000μM,
    ISO=0μM,
    name=:bar_sys,
    simplify=false
)
    @parameters begin
        b1ARtot = 0.00528μM
        Gstot = 3.83μM
        PDEtot = 22.85e-3μM
        PKItot = 0.18μM
        I1tot = 0.3μM
        PP1tot = 0.89μM
        RItot = 1.18μM
        RIItot = 0.118μM
        ACtot = 70.57e-3μM
        PKACII_LCCtot = 0.025μM
        PKAIItot = 0.059μM
        PKAII_KURtot = 0.025μM
        PP1_KURtot = 0.025μM
        LCCtot = 0.025μM
        PLBtot = 106μM
        PLMtot = 48μM
        TnItot = 70μM
        IKurtot = 0.025μM
        kf_LR = 1 / (μM * ms)                    # forward rate for ISO binding to b1AR
        kr_LR = 285Hz                          # reverse rate for ISO binding to b1AR
        kf_LRG = 1 / (μM * ms)                   # forward rate for ISO:b1AR association with Gs
        kr_LRG = 62Hz                            # reverse rate for ISO:b1AR association with Gs
        kf_RG = 1 / (μM * ms)                   # forward rate for b1AR association with Gs
        kr_RG = 33000Hz                          # reverse rate for b1AR association with Gs
        k_G_act = 16Hz                           # rate constant for Gs activation
        k_G_hyd = 0.8Hz                          # rate constant for G-protein hydrolysis
        k_G_reassoc = 1.21 / (μM * ms)           # rate constant for G-protein reassociation
        kf_bARK = 1.1Hz / mM                    # forward rate for b1AR phosphorylation by b1ARK
        kr_bARK = 2.2e-3Hz                       # reverse rate for b1AR phosphorylation by b1ARK
        kf_PKA = 3.6Hz / mM                      # forward rate for b1AR phosphorylation by PKA
        kr_PKA = 2.2e-3Hz                        # reverse rate for b1AR phosphorylation by PKA
        k_AC_basal = 0.2Hz                       # basal cAMP generation rate by AC
        Km_AC_basal = 1.03mM                     # basal AC affinity for ATP
        kr_AC_Gsa = 0.4μM                        # AC dissociation with Gsa
        kf_AC_Gsa = 1 / (μM * ms)                # forward rate for AC association with Gsa
        k_AC_Gsa = 8.5Hz                        # basal cAMP generation rate by AC:Gsa
        Km_AC_Gsa = 315.0μM                      # AC:Gsa affinity for ATP
        k_cAMP_PDE = 5Hz                         # cAMP hydrolysis rate by PDE
        k_cAMP_PDEp = 10Hz                       # cAMP hydrolysis rate by phosphorylated PDE
        Km_PDE_cAMP = 1.3μM                      # PDE affinity for cAMP
        k_PKA_PDE = 7.5Hz                        # rate constant for PDE phosphorylation by type 1 PKA
        k_PP_PDE = 1.5Hz                         # rate constant for PDE dephosphorylation by phosphatases
        kf_RC_cAMP = 1 / (μM * ms)               # Kd for PKA RC binding to cAMP
        kf_RCcAMP_cAMP = 1 / (μM * ms)           # Kd for PKA RC:cAMP binding to cAMP
        kf_RcAMPcAMP_C = 4.375 / (μM * ms)       # Kd for PKA R:cAMPcAMP binding to C
        kf_PKA_PKI = 1 / (μM * ms)               # Ki for PKA inhibition by PKI
        kr_RC_cAMP = 1.64 / ms                   # rate constant for PKA RC unbinding to cAMP
        kr_RCcAMP_cAMP = 9.14 / ms               # rate constant for PKA RC:cAMP unbinding to cAMP
        kr_RcAMPcAMP_C = 1 / ms                  # rate constant for PKA R:cAMPcAMP unbinding to C
        kr_PKA_PKI = 0.2Hz                       # reverse rate for PKA inhibition by PKI
        epsilon = 10                             # AKAP-mediated scaling factor
        k_PKA_I1 = 60Hz                          # rate constant for I-1 phosphorylation by type 1 PKA
        Km_PKA_I1 = 1.0μM                        # Km for I-1 phosphorylation by type 1 PKA
        Vmax_PP2A_I1 = 14Hz                      # Vmax for I-1 dephosphorylation by PP2A
        Km_PP2A_I1 = 1.0μM                       # Km for I-1 dephosphorylation by PP2A
        kr_PP1_I1 = 1.0Hz                        # Ki for PP1 inhibition by I-1
        kf_PP1_I1 = 1μM                          # kf for PP1 inhibition by I-1
        k_PKA_PLB = 54Hz
        Km_PKA_PLB = 21μM
        k_PP1_PLB = 8.5Hz
        Km_PP1_PLB = 7.0μM
        k_PKA_PLM = 54Hz
        Km_PKA_PLM = 21μM
        k_PP1_PLM = 8.5Hz
        Km_PP1_PLM = 7.0μM
        PP1_LCC = 0.025μM
        PP2A_LCC = 0.025μM
        k_PKA_LCC = 54Hz
        Km_PKA_LCC = 21μM
        k_PP1_LCC = 8.52Hz
        Km_PP1_LCC = 3μM
        k_PP2A_LCC = 10.1Hz
        Km_PP2A_LCC = 3μM
        PP2A_TnI = 0.67μM
        k_PKA_TnI = 54Hz
        Km_PKA_TnI = 21μM
        k_PP2A_TnI = 10.1Hz
        Km_PP2A_TnI = 4.1μM
        k_pka_KUR = 54Hz
        Km_pka_KUR = 21μM
        k_pp1_KUR = 8.52Hz
        Km_pp1_KUR = 7μM
    end

    sts = @variables begin
        LR(t) = 6.0e-5μM
        LRG(t) = 0.00294μM
        RG(t) = 7.0e-5μM
        GsaGTP(t) = 0.05214μM
        GsaGDP(t) = 0.00066μM
        b1AR_S464(t) = 0.00047μM
        b1AR_S301(t) = 0.0011μM
        AC_GsaGTP(t) = 0.00814μM
        PDEp(t) = 0.00589μM
        cAMP(t) = 1.50399μM
        RCcAMP_I(t) = 0.28552μM
        RCcAMP_II(t) = 0.00934μM
        RCcAMPcAMP_I(t) = 0.04698μM
        RCcAMPcAMP_II(t) = 0.00154μM
        PKACI(t) = 0.38375μM
        PKACII(t) = 0.06938μM
        PKACI_PKI(t) = 0.15239μM
        PKACII_PKI(t) = 0.02753μM
        I1p(t) = 0.00033μM
        I1p_PP1(t) = 0.19135μM
        PLBp(t) = 98.33936μM
        PLMp(t) = 41.19479μM
        TnIp(t) = 60.75646μM
        LCCap(t) = 0.01204μM
        LCCbp(t) = 0.01313μM
        KURp(t) = 0.01794μM
    end

    conservedvars = @variables begin
        b1AR(t)
        Gs(t)
        Gsby(t)
        AC(t)
        PDE(t)
        RC_I(t)
        RC_II(t)
        RcAMPcAMP_I(t)
        RcAMPcAMP_II(t)
        PKI(t)
        I1(t)
        PP1(t)
        PLB(t)
        PLM(t)
        TnI(t)
        LCCa(t)
        LCCb(t)
        KURn(t)
    end

    conservedeqs = [
        b1ARtot ~ b1AR + LR + LRG + RG + b1AR_S464 + b1AR_S301,
        Gstot ~ LRG + RG + Gs + GsaGDP + GsaGTP + AC_GsaGTP,
        Gstot ~ LRG + RG + Gs + Gsby,
        ACtot ~ AC + AC_GsaGTP,
        PDEtot ~ PDE + PDEp,
        RItot ~ RC_I + RCcAMP_I + RCcAMPcAMP_I + RcAMPcAMP_I,
        RIItot ~ RC_II + RCcAMP_II + RCcAMPcAMP_II + RcAMPcAMP_II,
        RcAMPcAMP_I ~ PKACI + PKACI_PKI,
        RcAMPcAMP_II ~ PKACII + PKACII_PKI,
        PKItot ~ PKI + PKACI_PKI+ PKACII_PKI,
        I1tot ~ I1 + I1p + I1p_PP1,
        PP1tot ~ PP1 + I1p_PP1,
        PLBtot ~ PLB + PLBp,
        PLMtot ~ PLM + PLMp,
        TnItot ~ TnI + TnIp,
        LCCtot ~ LCCa + LCCap,
        LCCtot ~ LCCb + LCCbp,
        IKurtot ~ KURn + KURp,
    ]

    rates = merge(Dict(sts .=> t - t), Dict(conservedvars .=> t - t))

    # Observables
    @variables begin
        LCCa_PKAp(t)
        LCCb_PKAp(t)
        fracPLBp(t)
        fracBLMp(t)
        TnI_PKAp(t)
        IKUR_PKAp(t)
    end

    obseqs = [
        LCCa_PKAp ~ LCCap / LCCtot,
        LCCb_PKAp ~ LCCbp / LCCtot,
        fracPLBp ~ PLBp / PLBtot,
        fracBLMp ~ PLMp / PLMtot,
        TnI_PKAp ~ TnIp / TnItot,
        IKUR_PKAp ~ KURp / IKurtot,
    ]

    # G-protein receptor
    add_rate!(rates, ISO * kf_LR, [b1AR], kr_LR, [LR]) # b1AR <--> LR
    add_rate!(rates, kf_LRG, [LR, Gs], kr_LRG, [LRG]) # LR + Gs <--> LRG
    add_rate!(rates, kf_RG, [b1AR, Gs], kr_RG, [RG]) # b1AR + Gs <--> RG
    add_rate!(rates, k_G_act, [RG], 0, [b1AR, GsaGTP, Gsby]) # (RG, LRG) --> b1AR + GsaGTP + Gsby
    add_rate!(rates, k_G_act, [LRG], 0, [b1AR, GsaGTP, Gsby])
    add_rate!(rates, k_G_hyd, [GsaGTP], 0, [GsaGDP])  # GsaGTP --> GsaGDP
    add_rate!(rates, k_G_reassoc, [GsaGDP, Gsby], 0, [Gs]) # GsaGDP + Gsby --> Gs
    # Ligand-mediated inactivation
    add_rate!(rates, kf_bARK, [LR], kr_bARK, [b1AR_S464]) # LR <--> b1AR_S464
    add_rate!(rates, kf_bARK, [LRG], 0, [b1AR_S464, Gs]) # LRG --> b1AR_S464 + Gs
    # PKA-mediated receptor inactivation
    add_rate!(rates, kf_PKA * PKACI, [b1AR], kr_PKA, [b1AR_S301]) # b1AR <--> b1AR_S301
    add_rate!(rates, kf_PKA * PKACI, [LR], 0, [b1AR_S301]) # LR --> b1AR_S301
    add_rate!(rates, kf_PKA * PKACI, [LRG], 0, [b1AR_S301, Gs]) # LRG --> b1AR_S301 + Gs
    # Adenylate cyclase and PDE
    add_rate!(rates, kf_AC_Gsa, [AC, GsaGTP], kr_AC_Gsa, [AC_GsaGTP]) # AC + GsaGTP <--> AC_GsaGTP
    add_rate!(rates, k_G_hyd, [AC_GsaGTP], 0, [AC, GsaGDP]) # AC_GsaGTP --> AC + GsaGDP
    add_rate!(rates, k_PKA_PDE * PKACII, [PDE], k_PP_PDE, [PDEp]) # PDE <--> PDEp
    v = k_AC_Gsa * AC_GsaGTP * hil(ATP, Km_AC_Gsa) + k_AC_basal * AC * hil(ATP, Km_AC_basal) - (k_cAMP_PDE * PDE + k_cAMP_PDEp * PDEp) * hil(cAMP, Km_PDE_cAMP)
    add_raw_rate!(rates, v, [], [cAMP]) # 0 <--> cAMP
    # cAMP disinhibition of PKA
    add_rate!(rates, kf_RC_cAMP, [RC_I, cAMP], kr_RC_cAMP, [RCcAMP_I]) # RC_I + cAMP <--> RCcAMP_I
    add_rate!(rates, kf_RC_cAMP, [RC_II, cAMP], kr_RC_cAMP, [RCcAMP_II]) # RC_II + cAMP <--> RCcAMP_II
    add_rate!(rates, kf_RCcAMP_cAMP, [RCcAMP_I, cAMP], kr_RCcAMP_cAMP, [RCcAMPcAMP_I]) # CcAMP_I + cAMP <--> RCcAMPcAMP_I
    add_rate!(rates, kf_RCcAMP_cAMP, [RCcAMP_II, cAMP], kr_RCcAMP_cAMP, [RCcAMPcAMP_II]) # RCcAMP_II + cAMP <--> RCcAMPcAMP_II
    add_rate!(rates, kf_RcAMPcAMP_C, [RCcAMPcAMP_I], kr_RcAMPcAMP_C, [RcAMPcAMP_I, PKACI]) # RCcAMPcAMP_I <--> RcAMPcAMP_I + PKACI
    add_rate!(rates, kf_RcAMPcAMP_C, [RCcAMPcAMP_II], kr_RcAMPcAMP_C, [RcAMPcAMP_II, PKACII]) # RCcAMPcAMP_II <--> RcAMPcAMP_II + PKACII
    add_rate!(rates, kf_PKA_PKI, [PKACI, PKI], kr_PKA_PKI, [PKACI_PKI]) # PKACI + PKI <--> PKACI_PKI
    add_rate!(rates, kf_PKA_PKI, [PKACII, PKI], kr_PKA_PKI, [PKACII_PKI]) # PKACII + PKI <--> PKACII_PKI
    # PKA modifications
    v = k_PKA_I1 * PKACI * hil(I1, Km_PKA_I1) - Vmax_PP2A_I1 * hil(I1p, Km_PP2A_I1)
    add_raw_rate!(rates, v, [I1], [I1p]) # I1 <=> I1p
    add_rate!(rates, kf_PP1_I1, [I1p, PP1], kr_PP1_I1, [I1p_PP1]) # I1p + PP1 <--> I1p_PP1
    v = k_PKA_PLB * PKACI * hil(PLB, Km_PKA_PLB) - k_PP1_PLB * PP1 * hil(PLBp, Km_PP1_PLB)
    add_raw_rate!(rates, v, [PLB], [PLBp]) # PLB <=> PLBp
    v = k_PKA_PLM * PKACI * hil(PLM, Km_PKA_PLM) - k_PP1_PLM * PP1 * hil(PLMp, Km_PP1_PLM)
    add_raw_rate!(rates, v, [PLM], [PLMp]) # PLM <=> PLMp
    v = k_PKA_TnI * PKACI * hil(TnI, Km_PKA_TnI) - k_PP2A_TnI * PP2A_TnI * hil(TnIp, Km_PP2A_TnI)
    add_raw_rate!(rates, v, [TnI], [TnIp]) # TnI <=> TnIp
    v = k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII * hil(LCCa * epsilon, Km_PKA_LCC) - k_PP2A_LCC * PP2A_LCC * hil(LCCap * epsilon, Km_PP2A_LCC)
    add_raw_rate!(rates, v, [LCCa], [LCCap]) # LCCa <=> LCCap
    v = k_PKA_LCC * (PKACII_LCCtot / PKAIItot) * PKACII * hil(LCCb * epsilon, Km_PKA_LCC) - k_PP1_LCC * PP1_LCC * hil(LCCbp * epsilon, Km_PP1_LCC)
    add_raw_rate!(rates, v, [LCCb], [LCCbp]) # LCCb <=> LCCbp
    v = k_pka_KUR * (PKAII_KURtot / PKAIItot) * PKACII * hil(KURn * epsilon, Km_pka_KUR) - PP1_KURtot * k_pp1_KUR * hil(KURp * epsilon, Km_pp1_KUR)
    add_raw_rate!(rates, v, [KURn], [KURp])

    rateeqs = [D(s) ~ rates[s] for s in sts]

    sys = ODESystem([rateeqs; conservedeqs; obseqs], t; name)

    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end
