module Module_BAR
    using DifferentialEquations
    using ModelingToolkit
    #= 
    This module describes the beta-adrenergic signaling pathway in mouse 
    ventricular myocyte, and this file was built upon the code developeded
    by Yang and Saucerman.
    Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
    feed-forward regulator of Ca2+ in beta-adrenergic signaling,
    accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular 
    Cardiology 52, 1048-1055.
    =#

    function morotti_barODE!(du,u,p,t)
    ## State variables
        (LR, LRG, RG, b1AR_S464, b1AR_S301, GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, PDEp,
         cAMPtot, RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, RC_II, 
         RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, I1p_PP1, I1ptot,
         PLBp, PLMp, LCCap, LCCbp, RyRp, TnIp, KS79, KS80, KSp, CFTRp, KURp) = u

    ## Parameters
        (ISO, LCCtot, RyRtot, PLBtot, TnItot, IKstot, ICFTRtot, PP1tot, IKurtot, PLMtot) = p

    ## Drug Concentrations
        # ISO = pin(1)    # (uM) isoproterenol concentration - Ltot
        FSK = 0         # (uM) forskolin concentration
        IBMX = 0        # (uM) IBMX concentration
    
    ## b-AR module

        b1ARtot         = 0.00528        # (uM) total b1-AR protein # MOUSE
        #b1ARtot=0.028; # RABBIT
        kf_LR           = 1              # (1/[uM ms]) forward rate for ISO binding to b1AR
        kr_LR           = 0.285          # (1/ms) reverse rate for ISO binding to b1AR
        kf_LRG          = 1              # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
        kr_LRG          = 0.062          # (1/ms) reverse rate for ISO:b1AR association with Gs
        kf_RG           = 1              # (1/[uM ms]) forward rate for b1AR association with Gs
        kr_RG           = 33             # (1/ms) reverse rate for b1AR association with Gs
        Gstot           = 3.83           # (uM) total Gs protein
        k_G_act         = 16e-3          # (1/ms) rate constant for Gs activation
        k_G_hyd         = 0.8e-3         # (1/ms) rate constant for G-protein hydrolysis
        k_G_reassoc     = 1.21           # (1/[uM ms]) rate constant for G-protein reassociation
        kf_bARK         = 1.1e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
        kr_bARK         = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
        kf_PKA          = 3.6e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
        kr_PKA          = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by PKA
        
        @variables b1ARact(t) b1ARact ~ b1ARtot - b1AR_S464 - b1AR_S301
        @variables b1AR(t) b1AR ~ b1ARact - LR - LRG - RG
        @variables Gs(t) Gs ~ Gstot - LRG - RG - Gsby
        
        du[1] = kf_LR*ISO*b1AR - kr_LR*LR + kr_LRG*LRG - kf_LRG*LR*Gs   # dLR
        du[2] = kf_LRG*LR*Gs - kr_LRG*LRG - k_G_act*LRG                 # dLRG
        du[3] = kf_RG*b1AR*Gs - kr_RG*RG - k_G_act*RG                   # dRG

        
        @variables bARK_desens(t) bARK_desens ~ kf_bARK*(LR+LRG)
        @variables bARK_resens(t) bARK_resens ~ kr_bARK*b1AR_S464
        @variables PKA_desens(t) PKA_desens ~ kf_PKA*PKACI*b1ARact
        @variables PKA_resens(t) PKA_resens ~ kr_PKA*b1AR_S301

        du[4] = bARK_desens - bARK_resens   # db1AR_S464
        du[5] = PKA_desens - PKA_resens     # db1AR_S301

        @variables G_act(t) G_act ~ k_G_act*(RG+LRG)
        @variables G_hyd(t)  G_hyd ~ k_G_hyd*GsaGTPtot
        @variables G_reasso(t) G_reasso ~ k_G_reassoc*GsaGDP*Gsby
        du[6] = G_act - G_hyd
        du[7] = G_hyd - G_reassoc
        du[8] = G_act - G_reassoc

    ## cAMP module

        ACtot           = 70.57e-3       # (uM) total adenylyl cyclase # MOUSE
        # ACtot=47e-3; # RABBIT
        ATP             = 5e3            # (uM) total ATP
        k_AC_basal      = 0.2e-3         # (1/ms) basal cAMP generation rate by AC
        Km_AC_basal     = 1.03e3         # (uM) basal AC affinity for ATP
        
        Kd_AC_Gsa       = 0.4            # (uM) Kd for AC association with Gsa
        kf_AC_Gsa       = 1              # (1/[uM ms]) forward rate for AC association with Gsa
        kr_AC_Gsa       = Kd_AC_Gsa      # (1/ms) reverse rate for AC association with Gsa
        
        k_AC_Gsa        = 8.5e-3         # (1/ms) basal cAMP generation rate by AC:Gsa
        Km_AC_Gsa       = 315.0          # (uM) AC:Gsa affinity for ATP
        
        Kd_AC_FSK       = 44.0           # (uM) Kd for FSK binding to AC
        k_AC_FSK        = 7.3e-3         # (1/ms) basal cAMP generation rate by AC:FSK
        Km_AC_FSK       = 860.0          # (uM) AC:FSK affinity for ATP
        
        PDEtot          = 22.85e-3       # (uM) total phosphodiesterase
        k_cAMP_PDE      = 5e-3           # (1/ms) cAMP hydrolysis rate by PDE
        k_cAMP_PDEp     = 2*k_cAMP_PDE   # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
        Km_PDE_cAMP     = 1.3            # (uM) PDE affinity for cAMP
        
        Kd_PDE_IBMX     = 30.0           # (uM) Kd_R2cAMP_C for IBMX binding to PDE
        k_PKA_PDE       = 7.5e-3         # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
        k_PP_PDE        = 1.5e-3         # (1/ms) rate constant for PDE dephosphorylation by phosphatases
        
        @variables cAMP(t) cAMP ~ cAMPtot - (RCcAMP_I+2*RCcAMPcAMP_I+2*RcAMPcAMP_I) - (RCcAMP_II+2*RCcAMPcAMP_II+2*RcAMPcAMP_II)
        @variables AC(t) AC ~ ACtot-AC_GsaGTP
        @variables GsaGTP(t) GsaGTP ~ GsaGTPtot - AC_GsaGTP
        du[9] = kf_AC_Gsa*GsaGTP*AC - kr_AC_Gsa*AC_GsaGTP       # dAC_GsaGTP 
        
        @variables AC_FSK(t) AC_FSK ~ FSK*AC/Kd_AC_FSK
        @variables AC_ACT_BASAL(t) AC_ACT_BASAL ~ k_AC_basal*AC*ATP/(Km_AC_basal+ATP)
        @variables AC_ACT_GSA(t) AC_ACT_GSA ~ k_AC_Gsa*AC_GsaGTP*ATP/(Km_AC_Gsa+ATP)
        @variables AC_ACT_FSK(t) AC_ACT_FSK ~ k_AC_FSK*AC_FSK*ATP/(Km_AC_FSK+ATP)
        
        PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX
        @variables PDE(t) PDE ~ PDEtot - PDE_IBMX - PDEp
        du[10] = k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp           # dPDEp
        @variables PDE_ACT(t) PDE_ACT ~ k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP)        
        du[11] = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT # dcAMPtot 

    ## PKA module

        PKItot          = 0.18           # (uM) total PKI
        kf_RC_cAMP      = 1              # (1/[uM ms]) Kd for PKA RC binding to cAMP
        kf_RCcAMP_cAMP  = 1              # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
        kf_RcAMPcAMP_C  = 4.375          # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
        kf_PKA_PKI      = 1              # (1/[uM ms]) Ki for PKA inhibition by PKI
        kr_RC_cAMP      = 1.64           # (1/ms) Kd for PKA RC binding to cAMP
        kr_RCcAMP_cAMP  = 9.14           # (1/ms) Kd for PKA RC:cAMP binding to cAMP
        kr_RcAMPcAMP_C  = 1              # (1/ms) Kd for PKA R:cAMPcAMP binding to C
        kr_PKA_PKI      = 2e-4           # (1/ms) Ki for PKA inhibition by PKI
        epsilon         = 10             # (-) AKAP-mediated scaling factor
        
        @variables PKI(t) PKI ~ PKItot - PKACI_PKI - PKACII_PKI

        # dRC_I
        du[12] = - kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I
        # dRCcAMP_I
        du[13] = - kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I
        # dRCcAMPcAMP_I
        du[14] = - kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI
        # dRcAMPcAMP_I
        du[15] = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I
        # dPKACI
        du[16] = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI
        # dPKACI_PKI
        du[17] = - kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI
        # dRC_II
        du[18] = - kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II
        # dRCcAMP_II
        du[19] = - kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II
        # dRCcAMPcAMP_II
        du[20] = - kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII
        # dRcAMPcAMP_II
        du[21] = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II
        # dPKACII
        du[22] = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI
        # dPKACII_PKI
        du[23] = - kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI

    ## I-1/PP1 module

        # PP1tot = pin(8) # PP1tot = 0.89  # (uM) total phosphatase 1
        I1tot           = 0.3            # (uM) total inhibitor 1
        k_PKA_I1        = 60e-3          # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
        Km_PKA_I1       = 1.0            # (uM) Km for I-1 phosphorylation by type 1 PKA
        Vmax_PP2A_I1    = 14.0e-3        # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
        Km_PP2A_I1      = 1.0            # (uM) Km for I-1 dephosphorylation by PP2A
        
        Ki_PP1_I1       = 1.0e-3         # (uM) Ki for PP1 inhibition by I-1
        kf_PP1_I1       = 1              # (uM) Ki for PP1 inhibition by I-1
        kr_PP1_I1       = Ki_PP1_I1      # (uM) Ki for PP1 inhibition by I-1
        
        @variables I1(t) I1 ~ I1tot - I1ptot
        @variables PP1(t) PP1 ~ PP1tot - I1p_PP1
        @variables I1p(t) I1p ~ I1ptot - I1p_PP1
        @variables I1_phosph(t) I1_phosph ~ k_PKA_I1*PKACI*I1/(Km_PKA_I1+I1)
        @variables I1_dephosph(t) I1_dephosph ~ Vmax_PP2A_I1*I1ptot/(Km_PP2A_I1+I1ptot)
        
        du[24] = kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1      # dI1p_PP1
        du[25] = I1_phosph - I1_dephosph                    # dI1ptot

    ## PLB module

        # PLBtot = pin(4)   #p(41) = PLBtot # PLBtot        [uM]
        k_PKA_PLB = 54e-3   #p(44) = 54     # k_pka_plb     [1/ms]
        Km_PKA_PLB = 21     #p(45) = 21     # Km_pka_plb    [uM]
        k_PP1_PLB = 8.5e-3  #p(46) = 8.5    # k_pp1_plb     [1/ms]
        Km_PP1_PLB = 7.0    #p(47) = 7.0    # Km_pp1_plb    [uM]
        
        @variables PLB(t) PLB ~ PLBtot - PLBp
        @variables PLB_phosph(t) PLB_phosph ~ k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB)
        @variables PLB_dephosph(t) PLB_dephosph ~ k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp)
        du[26] = PLB_phosph - PLB_dephosph  # dPLBp -> output

    ## PLM module (included 09/18/12) MOUSE

        # PLMtot = pin(10)  # p(102) = PLMtot # PLMtot        [uM]
        k_PKA_PLM = 54e-3   # p(103) = 54     # k_pka_plb     [1/ms]
        Km_PKA_PLM = 21     # p(104) = 21     # Km_pka_plb    [uM]
        k_PP1_PLM = 8.5e-3  # p(105) = 8.5    # k_pp1_plb     [1/ms]
        Km_PP1_PLM = 7.0    # p(106) = 7.0    # Km_pp1_plb    [uM]
        
        @variables PLM(t) PLM ~ PLMtot - PLMp
        @variables PLM_phosph(t) PLM_phosph ~ k_PKA_PLM*PKACI*PLM/(Km_PKA_PLM+PLM)
        @variables PLM_dephosph(t) PLM_dephosph ~ k_PP1_PLM*PP1*PLMp/(Km_PP1_PLM+PLMp)
        du[27] = PLM_phosph - PLM_dephosph  # dPLMp -> output

    ## LCC module

        PKAIItot        = 0.059 # (uM) total type 2 PKA # MOUSE
        # PKAIItot=0.084 # RABBIT
        
        # LCCtot = pin(2)       # p(53) = LCCtot    # LCCtot        [uM]
        PKACII_LCCtot = 0.025   #p(54) = 0.025      # PKAIIlcctot   [uM]
        PP1_LCC = 0.025         #p(55) = 0.025      # PP1lcctot     [uM]
        PP2A_LCC = 0.025        #p(56) = 0.025      # PP2Alcctot    [uM]
        k_PKA_LCC = 54e-3       #p(57) = 54         # k_pka_lcc     [1/ms]
        Km_PKA_LCC = 21         #p(58) = 21#*1.6    # Km_pka_lcc    [uM]
        k_PP1_LCC = 8.52e-3     #p(59) = 8.52       # k_pp1_lcc     [1/ms] RABBIT, MOUSE
                                #p(59) = 8.5        # k_pp1_lcc     [1/sec] RAT
        Km_PP1_LCC = 3          #p(60) = 3          # Km_pp1_lcc    [uM]
        k_PP2A_LCC = 10.1e-3    #p(61) = 10.1       # k_pp2a_lcc    [1/ms]
        Km_PP2A_LCC = 3         #p(62) = 3          # Km_pp2a_lcc   [uM]
        
        @variables PKACII_LCC(t) PKACII_LCC ~ (PKACII_LCCtot/PKAIItot)*PKACII
        @variables LCCa(t) LCCa ~ LCCtot - LCCap
        @variables LCCa_phosph(t) LCCa_phosph ~ epsilon*k_PKA_LCC*PKACII_LCC*LCCa/(Km_PKA_LCC+epsilon*LCCa)
        @variables LCCa_dephosph(t) LCCa_dephosp ~ epsilon*k_PP2A_LCC*PP2A_LCC*LCCap/(Km_PP2A_LCC+epsilon*LCCap)
        du[28] = LCCa_phosph - LCCa_dephosph # dLCCap -> output
        @variables LCCb(t) LCCb ~ LCCtot - LCCbp
        @variables LCCb_phosph(t) LCCb_phosph ~ epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb)
        @variables LCCb_dephosph(t) LCCb_dephosph ~ epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp)
        du[29] = LCCb_phosph - LCCb_dephosph # dLCCbp -> output

    ## RyR module (not included in Yang-Saucerman)

        # RyRtot = pin(3)       #p(63) = RyRtot # RyRtot        [uM]
        PKAIIryrtot = 0.034     #p(64) = 0.034  # PKAIIryrtot   [uM]
        PP1ryr = 0.034          #p(65) = 0.034  # PP1ryr        [uM]
        PP2Aryr = 0.034         #p(66) = 0.034  # PP2Aryr       [uM]
        kcat_pka_ryr = 54e-3    #p(67) = 54     # kcat_pka_ryr  [1/ms]
        Km_pka_ryr = 21         #p(68) = 21     # Km_pka_ryr    [uM]
        kcat_pp1_ryr = 8.52e-3  #p(69) = 8.52   # kcat_pp1_ryr  [1/ms]
        Km_pp1_ryr = 7          #p(70) = 7      # Km_pp1_ryr    [uM]
        kcat_pp2a_ryr = 10.1e-3 #p(71) = 10.1   # kcat_pp2a_ryr [1/ms]
        Km_pp2a_ryr = 4.1       #p(72) = 4.1    # Km_pp2a_ryr   [uM]
        
        PKACryr = (PKAIIryrtot/PKAIItot)*PKACII
        RyR = RyRtot-RyRp
        RyRPHOSPH = epsilon*kcat_pka_ryr*PKACryr*RyR/(Km_pka_ryr+epsilon*RyR)
        RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*RyRp/(Km_pp1_ryr+epsilon*RyRp)
        RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*RyRp/(Km_pp2a_ryr+epsilon*RyRp)
        du[30] = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A # dRyRp -> output

    ## TnI module

        # TnItot = pin(5)   #p(73) = TnItot # TnItot        [uM]
        PP2A_TnI = 0.67                     # PP2Atni       [uM]
        k_PKA_TnI = 54e-3                   # kcat_pka_tni  [1/ms]
        Km_PKA_TnI = 21                     # Km_pka_tni    [uM]
        k_PP2A_TnI = 10.1e-3                # kcat_pp2a_tni [1/ms]
        Km_PP2A_TnI = 4.1                   # Km_pp2a_tni   [uM]
        
        @variables TnI(t) TnI ~ TnItot - TnIp
        @variables TnI_phosph(t) TnI_phosph ~ k_PKA_TnI*PKACI*TnI/(Km_PKA_TnI+TnI)
        @variables TnI_dephosph(t) TnI_dephosph ~ k_PP2A_TnI*PP2A_TnI*TnIp/(Km_PP2A_TnI+TnIp)
        du[31] = TnI_phosph - TnI_dephosph # dTnIp -> output

    ## Iks module (not present in mouse)

        # IKstot = pin(6)
        # p(79) = IKstot # Iks_tot       [uM]
        # p(80) = 0.025  # Yotiao_tot    [uM]
        # p(81) = 0.1e-3 # K_yotiao      [uM] ** apply G589D mutation here **
        # p(82) = 0.025  # PKAII_ikstot  [uM]
        # p(83) = 0.025  # PP1_ikstot    [uM]
        # p(84) = 54     # k_pka_iks     [1/sec]
        # p(85) = 21     # Km_pka_iks    [uM]
        # p(86) = 8.52   # k_pp1_iks     [1/sec]
        # p(87) = 7      # Km_pp1_iks    [uM]
        # 
        # IksYot = y(27)*y(28)/p(81)           # [uM]
        # ydot(27) = p(79) - IksYot - y(27)    # [uM]
        # ydot(28) = p(80) - IksYot - y(28)    # [uM]
        # PKACiks = (IksYot/p(79))*(p(82)/p(34))*y(18)
        # PP1iks = (IksYot/p(79))*p(83)
        # Iks = p(79)-y(29)
        # IKS_PHOSPH = p(40)*p(84)*PKACiks*Iks/(p(85)+p(40)*Iks)
        # IKS_DEPHOSPH = p(40)*p(86)*PP1iks*y(29)/(p(87)+p(40)*y(29))
        # ydot(29) = IKS_PHOSPH-IKS_DEPHOSPH

        du[32] = 0 # dKS79  not ODE
        du[33] = 0 # dKS80  not ODE
        du[34] = 0 # dKSp -> output -> 0

    ## CFTR module (included 04/30/10)

        # ICFTRtot = pin(7)         #p(88) = ICFTRtot   # ICFTR_tot      [uM]
        # PKAII_CFTRtot = 0.025     #p(89) = 0.025      # PKAII_CFTRtot [uM]
        # PP1_CFTRtot = 0.025       #p(90) = 0.025      # PP1_CFTRtot   [uM]
        # k_pka_CFTR = 54e-3        #p(91) = 54         # k_pka_CFTR    [1/ms]
        # Km_pka_CFTR = 8.5         #p(92) = 8.5        # Km_pka_CFTR   [uM]
        # k_pp1_CFTR = 8.52e-3      #p(93) = 8.52       # k_pp1_CFTR    [1/ms]
        # Km_pp1_CFTR = 7           #p(94) = 7          # Km_pp1_CFTR   [uM]
        # 
        # CFTRn = ICFTRtot - CFTRp                      # Non-phos = tot - phos
        # PKAC_CFTR = (PKAII_CFTRtot/PKAIItot)*PKACII   # (PKACFTRtot/PKAIItot)*PKAIIact
        # CFTRphos = epsilon*CFTRn*PKAC_CFTR*k_pka_CFTR/(Km_pka_CFTR+epsilon*CFTRn)
        # CFTRdephos = PP1_CFTRtot*k_pp1_CFTR*epsilon*CFTRp/(Km_pp1_CFTR + epsilon*CFTRp)
        du[35] = 0 # dCFTRp -> output -> 0

    ## Ikur module (included 04/10/12) MOUSE

        # IKurtot = pin(9)      # p(95) = IKurtot   # Ikur_tot     [uM]
        PKAII_KURtot = 0.025    # p (96) = 0.025    # PKAII_KURtot [uM]
        PP1_KURtot = 0.025      # p(97) = 0.025     # PP1_KURtot   [uM]
        k_pka_KUR = 54e-3       # p(98) = 54        # k_pka_KUR    [1/ms]
        Km_pka_KUR = 21         # p(99) = 21        # Km_pka_KUR   [uM]
        k_pp1_KUR = 8.52e-3     # p(100) = 8.52     # k_pp1_KUR    [1/ms]
        Km_pp1_KUR = 7          # p(101) = 7        # Km_pp1_KUR   [uM]
        
        @variables KURn(t) KURn ~ IKurtot - KURp  # Non-phos = tot - phos
        @variables PKAC_KUR(t) PKAC_KUR ~ (PKAII_KURtot/PKAIItot)*PKACII    # (PKA_KURtot/PKAIItot)*PKAIIact
        @variables KURphos(t) KURphos ~ epsilon*KURn*PKAC_KUR*k_pka_KUR/(Km_pka_KUR+epsilon*KURn)
        @variables KURdephos(t) KURdephos ~ PP1_KURtot*k_pp1_KUR*epsilon*KURp/(Km_pp1_KUR+epsilon*KURp)
        du[36] = KURphos - KURdephos # dKURp -> output

        return du
    end

    export morotti_barODE

end
