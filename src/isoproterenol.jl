# Beta adrenergic system activated by isoproterenol
using ModelingToolkit

function get_bar_eqs(ATP=5000μM, ISO=0μM,)
    @parameters begin
        b1ARtot = 0.00528μM
        Gstot = 3.83μM
        PDEtot = 22.85e-3μM
        PKItot = 0.18μM
        I1tot = 0.3μM
        PLBtot = 106μM
        PLMtot = 48μM
        TnItot = 70μM
        LCCtot = 0.025μM
        PKACII_LCCtot = 0.025μM
        PKAIItot = 0.059μM
        PKAII_KURtot = 0.025μM
        IKurtot = 0.025μM
        PP1_KURtot = 0.025μM
        PP1tot = 0.89μM
        RItot = 1.18μM
        RIItot = 0.118μM
        ACtot = 70.57e-3μM
        kf_LR = 1 / (μM * ms)                    # (1/[uM ms]) forward rate for ISO binding to b1AR
        kr_LR = 0.285 / ms                # (1/ms) reverse rate for ISO binding to b1AR
        kf_LRG = 1 / (μM * ms)                   # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
        kr_LRG = 0.062 / ms               # (1/ms) reverse rate for ISO:b1AR association with Gs
        kf_RG = 1 / (μM * ms)                    # (1/[uM ms]) forward rate for b1AR association with Gs
        kr_RG = 33 / ms                   # (1/ms) reverse rate for b1AR association with Gs
        k_G_act = 16e-3 / ms              # (1/ms) rate constant for Gs activation
        k_G_hyd = 0.8e-3 / ms             # (1/ms) rate constant for G-protein hydrolysis
        k_G_reassoc = 1.21 / (μM * ms)           # (1/[uM ms]) rate constant for G-protein reassociation
        kf_bARK = 1.1e-6 / (μM * ms)             # (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
        kr_bARK = 2.2e-6 / ms             # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
        kf_PKA = 3.6e-6 / (μM * ms)              # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
        kr_PKA = 2.2e-6 / ms              # (1/ms) reverse rate for b1AR phosphorylation by PKA
        k_AC_basal = 0.2e-3 / ms          # (1/ms) basal cAMP generation rate by AC
        Km_AC_basal = 1.03e3μM         # (uM) basal AC affinity for ATP
        kr_AC_Gsa = 0.4μM              # (uM) AC dissociation with Gsa
        kf_AC_Gsa = 1 / (μM * ms)                # (1/[uM ms]) forward rate for AC association with Gsa
        k_AC_Gsa = 8.5e-3 / ms         # (1/ms) basal cAMP generation rate by AC:Gsa
        Km_AC_Gsa = 315.0μM          # (uM) AC:Gsa affinity for ATP
        k_cAMP_PDE = 5e-3 / ms            # (1/ms) cAMP hydrolysis rate by PDE
        k_cAMP_PDEp = 10e-3 / ms          # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
        Km_PDE_cAMP = 1.3μM            # (uM) PDE affinity for cAMP
        k_PKA_PDE = 7.5e-3 / ms           # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
        k_PP_PDE = 1.5e-3 / ms            # (1/ms) rate constant for PDE dephosphorylation by phosphatases
        kf_RC_cAMP = 1 / (μM * ms)               # (1/[uM ms]) Kd for PKA RC binding to cAMP
        kf_RCcAMP_cAMP = 1 / (μM * ms)           # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
        kf_RcAMPcAMP_C = 4.375 / (μM * ms)       # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
        kf_PKA_PKI = 1 / (μM * ms)               # (1/[uM ms]) Ki for PKA inhibition by PKI
        kr_RC_cAMP = 1.64 / ms          # (1/ms) Kd for PKA RC binding to cAMP
        kr_RCcAMP_cAMP = 9.14 / ms      # (1/ms) Kd for PKA RC:cAMP binding to cAMP
        kr_RcAMPcAMP_C = 1 / ms         # (1/ms) Kd for PKA R:cAMPcAMP binding to C
        kr_PKA_PKI = 2e-4 / ms          # (1/ms) Ki for PKA inhibition by PKI
        epsilon = 10                  # (-) AKAP-mediated scaling factor
        k_PKA_I1 = 60e-3 / ms              # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
        Km_PKA_I1 = 1.0μM               # (uM) Km for I-1 phosphorylation by type 1 PKA
        Vmax_PP2A_I1 = 14.0e-3 / ms        # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
        Km_PP2A_I1 = 1.0μM              # (uM) Km for I-1 dephosphorylation by PP2A
        kr_PP1_I1 = 1.0e-3 / ms            # (uM) Ki for PP1 inhibition by I-1
        kf_PP1_I1 = 1μM                 # (uM) kf for PP1 inhibition by I-1
        k_PKA_PLB = 54e-3 / ms
        Km_PKA_PLB = 21μM
        k_PP1_PLB = 8.5e-3 / ms  # [1/ms]
        Km_PP1_PLB = 7.0μM    # [uM]
        k_PKA_PLM = 54e-3 / ms   # [1/ms]
        Km_PKA_PLM = 21μM     # [uM]
        k_PP1_PLM = 8.5e-3 / ms  # [1/ms]
        Km_PP1_PLM = 7.0μM    # [uM]
        PP1_LCC = 0.025μM         # [uM]
        PP2A_LCC = 0.025μM        # [uM]
        k_PKA_LCC = 54e-3 / ms       # [1/ms]
        Km_PKA_LCC = 21μM         # [uM]
        k_PP1_LCC = 8.52e-3 / ms     # [1/ms] RABBIT, MOUSE
        Km_PP1_LCC = 3μM          # [uM]
        k_PP2A_LCC = 10.1e-3 / ms    # [1/ms]
        Km_PP2A_LCC = 3μM         # [uM]
        PP2A_TnI = 0.67μM         # [uM]
        k_PKA_TnI = 54e-3 / ms       # [1/ms]
        Km_PKA_TnI = 21μM         # [uM]
        k_PP2A_TnI = 10.1e-3 / ms    # [1/ms]
        Km_PP2A_TnI = 4.1μM       # [uM]
        k_pka_KUR = 54e-3 / ms       # [1/ms]
        Km_pka_KUR = 21μM         # [uM]
        k_pp1_KUR = 8.52e-3 / ms     # [1/ms]
        Km_pp1_KUR = 7μM          # [uM]
    end

    sts = @variables begin
        t
        b1AR(t)  ## Conserved
        LR(t) = 6.0e-5μM
        LRG(t) = 0.00294μM
        RG(t) = 7.0e-5μM
        Gs(t)    ## Conserved
        GsaGTP(t) = 0.06028μM - 0.00814μM
        GsaGDP(t) = 0.00066μM
        Gsby(t) = 0.06071μM
        b1AR_S464(t) = 0.00047μM
        b1AR_S301(t) = 0.0011μM
        AC(t)  ## Conserved
        AC_GsaGTP(t) = 0.00814μM
        PDE(t)  ## Conserved
        PDEp(t) = 0.00589μM
        cAMP(t) = 1.50399μM
        RC_I(t) ## conserved = 0.31134μM
        RCcAMP_I(t) = 0.28552μM
        RCcAMPcAMP_I(t) = 0.04698μM
        RcAMPcAMP_I(t) = 0.53564μM
        PKACI(t) = 0.38375μM
        PKI(t)  ## conserved
        PKACI_PKI(t) = 0.15239μM
        RC_II(t) ## conserved = 0.01018μM
        RCcAMP_II(t) = 0.00934μM
        RCcAMPcAMP_II(t) = 0.00154μM
        RcAMPcAMP_II(t) = 0.09691μM
        PKACII(t) = 0.06938μM
        PKACII_PKI(t) = 0.02753μM
        I1(t) ## conserved
        I1p(t) = 0.00033μM
        PP1(t) ## conserved
        I1p_PP1(t) = 0.19135μM
        PLB(t) ## conserved
        PLBp(t) = 98.33936μM
        PLM(t) ## conserved
        PLMp(t) = 41.19479μM
        TnI(t) ## conserved
        TnIp(t) = 60.75646μM
        LCCa(t) ## conserved
        LCCap(t) = 0.01204μM
        LCCb(t) ## conserved
        LCCbp(t) = 0.01313μM
        KURn(t) ## conserved
        KURp(t) = 0.01794μM
        # Phosphorylation proportions
        LCCa_PKAp(t)
        LCCb_PKAp(t)
    end

    D = Differential(t)
    rates = Dict(sts .=> Num(0))  ## Record accumulated rates

    # Beta adrenergic receptor
    add_rate!(rates, ISO * kf_LR * b1AR - kr_LR * LR, [b1AR], [LR])
    add_rate!(rates, kf_LRG * LR * Gs - kr_LRG * LRG, [LR, Gs], [LRG])
    add_rate!(rates, kf_RG * b1AR * Gs - kr_RG * RG, [b1AR, Gs], [RG])
    add_rate!(rates, RG*k_G_act, [RG], [b1AR, GsaGTP, Gsby])
    add_rate!(rates, LRG*k_G_act, [LRG], [b1AR, GsaGTP, Gsby])
    add_rate!(rates, k_G_hyd * GsaGTP, [GsaGTP], [GsaGDP])
    add_rate!(rates, k_G_reassoc * GsaGDP * Gsby, [GsaGDP, Gsby], [Gs])
    # Receptor inactivation
    add_rate!(rates, kf_bARK * LR, [LR], [b1AR_S464])
    add_rate!(rates, kf_bARK * LRG, [LRG], [b1AR_S464])
    add_rate!(rates, kr_bARK * b1AR_S464, [b1AR_S464], [b1AR])
    # PKA-mediated receptor inactivation
    add_rate!(rates, kf_PKA * PKACI * b1AR - kr_PKA * b1AR_S301, [b1AR], [b1AR_S301])
    add_rate!(rates, kf_PKA * PKACI * LR, [LR], [b1AR_S301])
    add_rate!(rates, kf_PKA * PKACI * LRG, [LRG], [b1AR_S301, Gs])
    add_rate!(rates, kf_PKA * PKACI * RG, [RG], [b1AR_S301, Gs])
    # Adenylate cyclase and PDE
    add_rate!(rates, kf_AC_Gsa * AC * GsaGTP - kr_AC_Gsa * AC_GsaGTP, [AC, GsaGTP], [AC_GsaGTP])
    add_rate!(rates, k_G_hyd * AC_GsaGTP, [AC_GsaGTP], [AC, GsaGDP])
    add_rate!(rates, k_PKA_PDE * PKACII * PDE - k_PP_PDE * PDEp, [PDE], [PDEp])
    add_rate!(rates, k_AC_Gsa * AC_GsaGTP * hil(ATP, Km_AC_Gsa), [], [cAMP])
    add_rate!(rates, k_AC_basal * AC * hil(ATP, Km_AC_basal), [], [cAMP])
    add_rate!(rates, (k_cAMP_PDE * PDE + k_cAMP_PDEp * PDEp) * hil(cAMP, Km_PDE_cAMP), [cAMP], [])
    # cAMP activating PKA
    add_rate!(rates, kf_RC_cAMP * RC_I * cAMP - kr_RC_cAMP * RCcAMP_I, [RC_I, cAMP], [RCcAMP_I])
    add_rate!(rates, kf_RCcAMP_cAMP * RCcAMP_I * cAMP - kr_RCcAMP_cAMP * RCcAMPcAMP_I, [RCcAMP_I, cAMP], [RCcAMPcAMP_I])
    add_rate!(rates, kf_RcAMPcAMP_C * RCcAMPcAMP_I - kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI, [RCcAMPcAMP_I], [RcAMPcAMP_I, PKACI])
    add_rate!(rates, kf_PKA_PKI*PKACI* PKI - kr_PKA_PKI * PKACI_PKI, [PKACI, PKI], [PKACI_PKI])
    add_rate!(rates, kf_RC_cAMP * RC_II * cAMP - kr_RC_cAMP * RCcAMP_II, [RC_II, cAMP], [RCcAMP_II])
    add_rate!(rates, kf_RCcAMP_cAMP * RCcAMP_II * cAMP - kr_RCcAMP_cAMP * RCcAMPcAMP_II, [RCcAMP_II, cAMP], [RCcAMPcAMP_II])
    add_rate!(rates, kf_RcAMPcAMP_C * RCcAMPcAMP_II - kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII, [RCcAMPcAMP_II], [RcAMPcAMP_II, PKACII])
    add_rate!(rates, kf_PKA_PKII*PKACII* PKI - kr_PKA_PKI * PKACII_PKI, [PKACII, PKI], [PKACII_PKI])

    # PKA modifications
    add_rate!(rates, k_PKA_I1 * PKACI * hil(I1, Km_PKA_I1) - Vmax_PP2A_I1 * hil(I1p, Km_PP2A_I1), [I1], [I1p])
    add_rate!(rates, kf_PP1_I1 * I1p * PP1 - kr_PP1_I1 * I1p_PP1, [I1p, PP1], [I1p_PP1])
    add_rate!(rates, k_PKA_PLB * PKACI * hil(PLB, Km_PKA_PLB) - k_PP1_PLB * PP1 * hil(PLBp, Km_PP1_PLB), [PLB], [PLBp])
    add_rate!(rates, k_PKA_PLM * PKACI * hil(PLM, Km_PKA_PLM) - k_PP1_PLM * PP1 * hil(PLMp, Km_PKA_PLM), [PLM], [PLMp])
    add_rate!(rates, k_PKA_TnI * PKACI * hil(TnI, Km_PKA_TnI) - k_PP2A_TnI * PP2A_TnI * hil(TnIp, Km_PP2A_TnI), [TnI], [TnIp])

    pkac2 = (PKACII_LCCtot / PKAIItot) * PKACII
    kmlcc = Km_PKA_LCC / epsilon
    add_rate!(rates, k_PKA_LCC * pkac2 * hil(LCCa, kmlcc) - k_PP2A_LCC * PP2A_LCC * hil(LCCap, Km_PP2A_LCC / epsilon), [LCCa] ,[LCCap])
    add_rate!(rates, k_PKA_LCC * pkac2 * hil(LCCb, kmlcc) - k_PP1_LCC * PP1_LCC * hil(LCCbp, Km_PP1_LCC / epsilon), [LCCb], [LCCbp])

    pkac2 = (PKAII_KURtot / PKAIItot) * PKACII
    add_rate!(rates, k_pka_KUR * pkac2 * hil(KURn, Km_pka_KUR / epsilon) - PP1_KURtot * k_pp1_KUR * hil(KURp, Km_pp1_KUR / epsilon), [KURn], [KURp])

    conservedeqs = [
        b1AR ~ b1ARtot - LR - LRG - RG - b1AR_S464 - b1AR_S301,
        Gs ~ Gstot - Gsby - LRG - RG,
        AC ~ ACtot - AC_GsaGTP,
        PDE ~ PDEtot - PDEp,
        RC_I ~ RItot - RCcAMP_I - RCcAMPcAMP_I - RcAMPcAMP_I,
        RC_II ~ RIItot - RCcAMP_I - RCcAMPcAMP_I - RcAMPcAMP_I,
        PKI ~ PKItot - PKACI_PKI - PKACII_PKI,
        PP1 ~ PP1tot - I1p_PP1,
        I1 ~ I1tot - I1p - I1p_PP1,
        PLB ~ PLBtot - PLBp,
        PLM ~ PLMtot - PLMp,
        LCCa ~ LCCtot - LCCap,
        LCCb ~ LCCtot - LCCbp,
        TnI ~ TnItot - TnIp,
        KURn ~ IKurtot - KURp,
    ]

    propeqs = [
        LCCa_PKAp ~ LCCap/LCCtot,
        LCCb_PKAp ~ LCCbp/LCCtot,
    ]

    odeeqs = [D(x) ~ rates[x] for x in (LR, LRG, RG, GsaGTP, GsaGDP, Gsby, b1AR_S464, b1AR_S301, AC_GsaGTP, PDEp, cAMP, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, I1p, I1p_PP1, PLBp, PLMp, TnIp, LCCap, LCCbp, KURp)]

    return [conservedeqs; propeqs; odeeqs]
end
