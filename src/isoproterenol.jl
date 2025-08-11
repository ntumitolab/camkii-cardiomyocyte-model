"Beta-adrenergic system activated by isoproterenol"
function get_bar_sys(ATP=5000μM, ISO=0μM; name=:bar_sys, simplify=false)
    @parameters begin
        b1ARtot = 5.28nM
        Gstot = 3.83μM
        PDEtot = 22.85nM
        PKItot = 180nM
        I1tot = 300nM
        PP1totBA = 890nM
        RItot = 1180nM
        RIItot = 118nM
        ACtot = 70.57nM
        PKACII_LCCtotBA = 25nM
        PKAIItot = 59nM
        PKAII_KURtot = 25nM
        PP1_KURtot = 25nM
        LCCtotBA = 25nM
        PLBtotBA = 106μM
        PLMtotBA = 48μM
        TnItotBA = 70μM
        IKurtotBA = 25nM
        RyRtotBA = 135nM                        # Total RyR content for the BAR module
        PKAII_RyRtot = 34nM
        kf_LR = 1 / μM / ms                     # forward rate for ISO binding to b1AR
        kr_LR = 285Hz                           # reverse rate for ISO binding to b1AR
        kf_LRG = 1 / μM / ms                    # forward rate for ISO:b1AR association with Gs
        kr_LRG = 62Hz                           # reverse rate for ISO:b1AR association with Gs
        kf_RG = 1 / μM / ms                     # forward rate for b1AR association with Gs
        kr_RG = 33/ms                        # reverse rate for b1AR association with Gs
        k_G_act = 16Hz                          # rate constant for Gs activation
        k_G_hyd = 0.8Hz                         # rate constant for G-protein hydrolysis
        k_G_reassoc = 1.21 / μM / ms            # rate constant for G-protein reassociation
        kf_bARK = 1.1e-3Hz                      # forward rate for b1AR phosphorylation by b1ARK
        kr_bARK = 2.2e-3Hz                      # reverse rate for b1AR phosphorylation by b1ARK
        kf_PKA = 3.6Hz / mM                     # forward rate for b1AR phosphorylation by PKA
        kr_PKA = 2.2e-3Hz                       # reverse rate for b1AR phosphorylation by PKA
        k_AC_basal = 0.2Hz                      # basal cAMP generation rate by AC
        Km_AC_basal = 1.03mM                    # basal AC affinity for ATP
        kr_AC_Gsa = 0.4 / ms                    # AC dissociation with Gsa
        kf_AC_Gsa = 1 / μM / ms                 # forward rate for AC association with Gsa
        k_AC_Gsa = 8.5Hz                        # basal cAMP generation rate by AC:Gsa
        Km_AC_Gsa = 315.0μM                     # AC:Gsa affinity for ATP
        k_cAMP_PDE = 5Hz                        # cAMP hydrolysis rate by PDE
        k_cAMP_PDEp = 10Hz                      # cAMP hydrolysis rate by phosphorylated PDE
        Km_PDE_cAMP = 1.3μM                     # PDE affinity for cAMP
        k_PKA_PDE = 7.5Hz                       # rate constant for PDE phosphorylation by type 1 PKA
        k_PP_PDE = 1.5Hz                        # rate constant for PDE dephosphorylation by phosphatases
        kf_RC_cAMP = 1 / μM / ms                # Kd for PKA RC binding to cAMP
        kf_RCcAMP_cAMP = 1 / μM / ms            # Kd for PKA RC:cAMP binding to cAMP
        kf_RcAMPcAMP_C = 4.375 / ms             # Kd for PKA R:cAMPcAMP binding to C
        kf_PKA_PKI = 1 / μM / ms                # Ki for PKA inhibition by PKI
        kr_RC_cAMP = 1.64 / ms                  # rate constant for PKA RC unbinding to cAMP
        kr_RCcAMP_cAMP = 9.14 / ms              # rate constant for PKA RC:cAMP unbinding to cAMP
        kr_RcAMPcAMP_C = 1 / μM / ms            # rate constant for PKA R:cAMPcAMP binding to C
        kr_PKA_PKI = 0.2Hz                      # reverse rate for PKA inhibition by PKI
        epsilon = 10                            # AKAP-mediated scaling factor
        k_PKA_I1 = 60Hz / μM                    # rate constant for I-1 phosphorylation by type 1 PKA
        Km_PKA_I1 = 1.0μM                       # Km for I-1 phosphorylation by type 1 PKA
        Vmax_PP2A_I1 = 14Hz                     # Vmax for I-1 dephosphorylation by PP2A
        Km_PP2A_I1 = 1μM                        # Km for I-1 dephosphorylation by PP2A
        kr_PP1_I1 = 1 / ms                      # Ki for PP1 inhibition by I-1
        kf_PP1_I1 = 1 / μM / ms                 # kf for PP1 inhibition by I-1
        k_PKA_PLB = 54Hz / μM
        Km_PKA_PLB = 21μM
        k_PP1_PLB = 8.5Hz / μM
        Km_PP1_PLB = 7.0μM
        k_PKA_PLM = 54Hz / μM
        Km_PKA_PLM = 21μM
        k_PP1_PLM = 8.5Hz / μM
        Km_PP1_PLM = 7.0μM
        PP1_LCC = 0.025μM
        PP2A_LCC = 0.025μM
        k_PKA_LCC = 54Hz / μM
        Km_PKA_LCC = 21μM
        k_PP1_LCC = 8.52Hz / μM
        Km_PP1_LCC = 3μM
        k_PP2A_LCC = 10.1Hz / μM
        Km_PP2A_LCC = 3μM
        PP2A_TnI = 0.67μM
        k_PKA_TnI = 54Hz / μM
        Km_PKA_TnI = 21μM
        k_PP2A_TnI = 10.1Hz / μM
        Km_PP2A_TnI = 4.1μM
        k_pka_KUR = 54Hz / μM
        Km_pka_KUR = 21μM
        k_pp1_KUR = 8.52Hz / μM
        Km_pp1_KUR = 7μM
        PP1_RyR = 34nM
        PP2A_RyR = 34nM
        kcat_pka_RyR = 54Hz / μM
        Km_pka_RyR = 21μM
        kcat_pp1_RyR = 8.52Hz / μM
        Km_pp1_RyR = 7μM
        kcat_pp2a_RyR = 10.1Hz / μM
        Km_pp2a_RyR = 4.1μM
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
        RyRp(t) = 0μM
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
        RyRn(t)
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
        PKItot ~ PKI + PKACI_PKI + PKACII_PKI,
        I1tot ~ I1 + I1p + I1p_PP1,
        PP1totBA ~ PP1 + I1p_PP1,
        PLBtotBA ~ PLB + PLBp,
        PLMtotBA ~ PLM + PLMp,
        TnItotBA ~ TnI + TnIp,
        LCCtotBA ~ LCCa + LCCap,
        LCCtotBA ~ LCCb + LCCbp,
        IKurtotBA ~ KURn + KURp,
        RyRtotBA ~ RyRp + RyRn
    ]

    rates = merge(Dict(sts .=> Num(0)), Dict(conservedvars .=> Num(0)))

    # Observables
    @variables begin
        LCCa_PKAp(t)
        LCCb_PKAp(t)
        fracPLBp(t)
        fracPLMp(t)
        TnI_PKAp(t)
        IKUR_PKAp(t)
        RyR_PKAp(t)
    end

    obseqs = [
        LCCa_PKAp ~ LCCap / LCCtotBA,
        LCCb_PKAp ~ LCCbp / LCCtotBA,
        fracPLBp ~ PLBp / PLBtotBA,
        fracPLMp ~ PLMp / PLMtotBA,
        TnI_PKAp ~ TnIp / TnItotBA,
        IKUR_PKAp ~ KURp / IKurtotBA,
        RyR_PKAp ~ RyRp / RyRtotBA
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
    vf = k_AC_Gsa * AC_GsaGTP * hil(ATP, Km_AC_Gsa) + k_AC_basal * AC * hil(ATP, Km_AC_basal)
    vr = (k_cAMP_PDE * PDE + k_cAMP_PDEp * PDEp) * hil(cAMP, Km_PDE_cAMP)
    add_raw_rate!(rates, vf - vr, [], [cAMP]) # 0 <--> cAMP
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
    vf = k_PKA_I1 * PKACI * hil(I1, Km_PKA_I1)
    vr = Vmax_PP2A_I1 * hil(I1p, Km_PP2A_I1)
    add_raw_rate!(rates, vf - vr, [I1], [I1p]) # I1 <=> I1p
    add_rate!(rates, kf_PP1_I1, [I1p, PP1], kr_PP1_I1, [I1p_PP1]) # I1p + PP1 <--> I1p_PP1
    vf = k_PKA_PLB * PKACI * hil(PLB, Km_PKA_PLB)
    vr = k_PP1_PLB * PP1 * hil(PLBp, Km_PP1_PLB)
    add_raw_rate!(rates, vf - vr, [PLB], [PLBp]) # PLB <=> PLBp
    vf = k_PKA_PLM * PKACI * hil(PLM, Km_PKA_PLM)
    vr = k_PP1_PLM * PP1 * hil(PLMp, Km_PP1_PLM)
    add_raw_rate!(rates, vf - vr, [PLM], [PLMp]) # PLM <=> PLMp
    vf = k_PKA_TnI * PKACI * hil(TnI, Km_PKA_TnI)
    vr = k_PP2A_TnI * PP2A_TnI * hil(TnIp, Km_PP2A_TnI)
    add_raw_rate!(rates, vf - vr, [TnI], [TnIp]) # TnI <=> TnIp
    vf = k_PKA_LCC * (PKACII_LCCtotBA / PKAIItot) * PKACII * hil(LCCa * epsilon, Km_PKA_LCC)
    vr = k_PP2A_LCC * PP2A_LCC * hil(LCCap * epsilon, Km_PP2A_LCC)
    add_raw_rate!(rates, vf - vr, [LCCa], [LCCap]) # LCCa <=> LCCap
    vf = k_PKA_LCC * (PKACII_LCCtotBA / PKAIItot) * PKACII * hil(LCCb * epsilon, Km_PKA_LCC)
    vr = k_PP1_LCC * PP1_LCC * hil(LCCbp * epsilon, Km_PP1_LCC)
    add_raw_rate!(rates, vf - vr, [LCCb], [LCCbp]) # LCCb <=> LCCbp
    vf = k_pka_KUR * (PKAII_KURtot / PKAIItot) * PKACII * hil(KURn * epsilon, Km_pka_KUR)
    vr = k_pp1_KUR * PP1_KURtot * hil(KURp * epsilon, Km_pp1_KUR)
    add_raw_rate!(rates, vf - vr, [KURn], [KURp]) # KURn <=> KURp
    vf = kcat_pka_RyR * (PKAII_RyRtot / PKAIItot) * PKACII * hil(RyRn * epsilon, Km_pka_RyR)
    vr = kcat_pp1_RyR * PP1_RyR * hil(RyRp * epsilon, Km_pp1_RyR) + kcat_pp2a_RyR * PP2A_RyR * hil(RyRp * epsilon, Km_pp2a_RyR)
    add_raw_rate!(rates, vf - vr, [RyRn], [RyRp]) # RyRn <=> RyRp

    rateeqs = [D(s) ~ rates[s] for s in sts]

    sys = System([rateeqs; conservedeqs; obseqs], t; name)
    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end

function get_bar_eqs_reduced(ISO=0μM)
    @parameters begin
        PKACI_basal = 0.0734  ## basal activity
        PKACI_activated = 0.1994
        PKACI_KM = 0.0139μM
        PKACII_basal = 0.1840  ## basal activity
        PKACII_activated = 0.3444
        PKACII_KM = 0.01025μM
        PP1_basal = 0.8927
        PP1_activated = 0.0492
        PP1_KI = 0.00636μM
        PLBp_basal = 0.0830
        PLBp_activated = 0.7923
        PLBp_KM = 0.00594μM
        PLBp_nHill = 1.839
        PLMp_basal = 0.1177
        PLMp_activated = 0.6617
        PLMp_KM = 0.00818μM
        PLMp_nHill = 1.3710
        TnIp_basal = 0.0675
        TnIp_activated = 0.7482
        TnIp_KM = 0.007856μM
        TnIp_nHill = 1.6973
        LCCap_basal = 0.2208
        LCCap_activated = 0.2334
        LCCap_KM = 0.007263μM
        LCCbp_basal = 0.2520
        LCCbp_activated = 0.2456
        LCCbp_KM = 0.00696μM
        KURp_basal = 0.4394
        KURp_activated = 0.2557
        KURp_KM = 0.00558μM
        RyRp_basal = 0.2054
        RyRp_activated = 0.2399
        RyRp_KM = 0.00751μM
    end

    vs = @variables begin
        LCCa_PKAp(t)
        LCCb_PKAp(t)
        fracPKACI(t)
        fracPKACII(t)
        fracPP1(t)
        fracPLBp(t)
        fracPLMp(t)
        TnI_PKAp(t)
        IKUR_PKAp(t)
        RyR_PKAp(t)
    end

    ## Fitted activities
    eqs_bar = [
        fracPKACI ~ PKACI_basal + PKACI_activated * hil(ISO, PKACI_KM),
        fracPKACII ~ PKACII_basal + PKACII_activated * hil(ISO, PKACII_KM),
        fracPP1 ~ PP1_basal + PP1_activated * hilr(ISO, PP1_KI),
        fracPLBp ~ PLBp_basal + PLBp_activated * hil(ISO, PLBp_KM, PLBp_nHill),
        fracPLMp ~ PLMp_basal + PLMp_activated * hil(ISO, PLMp_KM, PLMp_nHill),
        TnI_PKAp ~ TnIp_basal + TnIp_activated * hil(ISO, TnIp_KM, TnIp_nHill),
        LCCa_PKAp ~ LCCap_basal + LCCap_activated * hil(ISO, LCCap_KM),
        LCCb_PKAp ~ LCCbp_basal + LCCbp_activated * hil(ISO, LCCbp_KM),
        IKUR_PKAp ~ KURp_basal + KURp_activated * hil(ISO, KURp_KM),
        RyR_PKAp ~ RyRp_basal + RyRp_activated * hil(ISO, RyRp_KM),
    ]

    return (; eqs_bar, fracPKACI, fracPKACII, fracPP1, fracPLBp, fracPLMp, TnI_PKAp, LCCa_PKAp, LCCb_PKAp, IKUR_PKAp, RyR_PKAp)
end

"Algebraic fitted beta-adrenergic system"
function get_bar_sys_reduced(ISO=0μM; name=:bar_sys)
    @unpack eqs_bar = get_bar_eqs_reduced(ISO)
    return ODESystem(eqs_bar, t; name)
end
