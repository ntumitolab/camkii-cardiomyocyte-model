using ModelingToolkit

"""
CaMKII system with ROS activation (full model)
"""
function get_camkii_eqs(;
    Ca=0μM,
    ROS=0μM,
    binding_To_PCaMK=0,   ## 0.1 for T287D mutation
)

    @independent_variables t
    D = Differential(t)

    @parameters begin
        CAM_T = 30μM            ## Total calmodulin Concentration
        CAMKII_T = 70μM         ## Total CaMKII Concentration
        k_1C_on = 5Hz / μM      ## 1.2-9.6uM-1Hz
        k_1C_off = 50Hz         ## 10-70 Hz
        k_2C_on = 10Hz / μM     ## 5-25uM-1Hz.
        k_2C_off = 10Hz         ## 8.5-10Hz.
        ## N-lobe
        k_1N_on = 100Hz / μM    ## 25-260uM-1Hz
        k_1N_off = 2000Hz       ## 1000-4000 Hz
        k_2N_on = 200Hz / μM    ## 50-300uM-1Hz.
        k_2N_off = 500Hz        ## 500-1000.Hz

        ## Ca2+ binding to CaM-CAMKII (KCaM)
        ## C-lobe
        k_K1C_on = 44Hz / μM
        k_K1C_off = 33Hz
        k_K2C_on = 44Hz / μM
        k_K2C_off = 0.8Hz ## 0.49-4.9Hz
        ## N-lobe
        k_K1N_on = 76Hz / μM
        k_K1N_off = 300Hz
        k_K2N_on = 76Hz / μM
        k_K2N_off = 20Hz ## 6-60Hz

        ## CaM binding to CaMKII
        kCaM0_on = 3.8e-3Hz / μM ## Changed to Pepke's value (Chang: 3.8)
        kCaM0_off = 5.5Hz
        kCaM2C_on = 0.5Hz / μM  # 0.92 μM-1Hz
        kCaM2C_off = 6.8Hz
        kCaM2N_on = 0.12Hz / μM
        kCaM2N_off = 1.7Hz
        kCaM4_on = 15Hz / μM  # 14-60 uM-1Hz
        kCaM4_off = 1.5Hz  # 1.1 - 2.3 Hz
        kCaM0P_off = inv(3second)
        kCaM2CP_off = inv(3second)
        kCaM2NP_off = inv(3second)
        kCaM4P_off = inv(3second)
        k_phosCaM = 5Hz # 30Hz
        k_dephospho = inv(6second)
        k_P1_P2 = inv(60second)
        k_P2_P1 = inv(15second)

        ## Oxidation / reduction of Met
        k_B_OX = 291Hz / mM
        k_P_OXP = 291Hz / mM
        k_OX_B = inv(45second)
        k_OXP_P = inv(45second)
    end

    ## State variables
    sts = @variables begin
        Ca2CaM_C(t) = 0
        Ca2CaM_N(t) = 0
        Ca4CaM(t) = 0
        CaM0_CaMK(t) = 0
        Ca2CaM_C_CaMK(t) = 0
        Ca2CaM_N_CaMK(t) = 0
        Ca4CaM_CaMK(t) = 0
        CaM0_CaMKP(t) = 0
        Ca2CaM_C_CaMKP(t) = 0
        Ca2CaM_N_CaMKP(t) = 0
        Ca4CaM_CaMKP(t) = 0
        Ca4CaM_CaMKOX(t) = 0
        Ca4CaM_CaMKPOX(t) = 0
        CaMKP(t) = 0
        CaMKP2(t) = 0
        CaMKPOX(t) = 0
        CaMKOX(t) = 0
    end

    ## Dependent variables
    @variables begin
        CaMKAct(t)          ## Active CaMKII fraction
        CaM0(t)             ## Apo CaM
        CaMK(t)             ## Apo (inactive) CaMKII
    end

    ## Ca binding/unbinding reaction rates
    function _ca_cam(ca, k1on, k1off, k2on, k2off)
        den = ca * k2on + k1off
        return (ca * ca * k1on * k2on / den, k1off * k2off / den)
    end

    ## Rates counter
    rates = Dict()
    ## Two Ca2+ ions bind to C (high affinity) or N (low affinity)-lobe of CaM
    kon, koff = _ca_cam(Ca, k_1C_on, k_1C_off, k_2C_on, k_2C_off)
    ## CaM0 <--> Ca2CaM_C
    add_rate!(rates, kon, CaM0, koff, Ca2CaM_C)
    ## Ca2CaM_N <--> Ca4CaM
    add_rate!(rates, kon, Ca2CaM_N, koff, Ca4CaM)
    kon, koff = _ca_cam(Ca, k_1N_on, k_1N_off, k_2N_on, k_2N_off)
    ## CaM0 <--> Ca2CaM_N
    add_rate!(rates, kon, CaM0, koff, Ca2CaM_N)
    ## Ca2CaM_C <--> Ca4CaM
    add_rate!(rates, kon, Ca2CaM_C, koff, Ca4CaM)

    ## Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII(P) complexes
    kon, koff = _ca_cam(Ca, k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off)
    for (a, b) in zip((CaM0_CaMK, Ca2CaM_N_CaMK, CaM0_CaMKP, Ca2CaM_N_CaMKP), (Ca2CaM_C_CaMK, Ca4CaM_CaMK, Ca2CaM_C_CaMKP, Ca4CaM_CaMKP))
        add_rate!(rates, kon, a, koff, b)
    end
    kon, koff = _ca_cam(Ca, k_K1N_on, k_K1N_off, k_K2N_on, k_K2N_off)
    for (a, b) in zip((CaM0_CaMK, Ca2CaM_C_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP), (Ca2CaM_N_CaMK, Ca4CaM_CaMK, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP))
        add_rate!(rates, kon, a, koff, b)
    end

    ## CaM binding to CaMKII / CaMkII-P / CaMkII-POX / CaMkII-OX
    add_rate!(rates, kCaM0_on, [CaM0, CaMK], kCaM0_off, [CaM0_CaMK]) # CaM0 + CaMK <--> CaM0_CaMK
    add_rate!(rates, kCaM2C_on, [Ca2CaM_C, CaMK], kCaM2C_off, [Ca2CaM_C_CaMK]) # Ca2CaM_C + CaMK <--> Ca2CaM_C_CaMK
    add_rate!(rates, kCaM2N_on, [Ca2CaM_N, CaMK], kCaM2N_off, [Ca2CaM_N_CaMK]) # Ca2CaM_N + CaMK <--> Ca2CaM_N_CaMK
    add_rate!(rates, kCaM4_on, [Ca4CaM, CaMK], kCaM4_off, [Ca4CaM_CaMK]) # Ca4CaM + CaMK <--> Ca4CaM_CaMK
    add_rate!(rates, kCaM0_on * binding_To_PCaMK, [CaM0, CaMKP], kCaM0P_off, [CaM0_CaMKP])  # CaM0 + CaMKP <--> CaM0_CaMKP
    add_rate!(rates, kCaM2C_on * binding_To_PCaMK, [Ca2CaM_C, CaMKP], kCaM2CP_off, [Ca2CaM_C_CaMKP]) # Ca2CaM_C + CaMKP <--> Ca2CaM_C_CaMKP
    add_rate!(rates, kCaM2N_on * binding_To_PCaMK, [Ca2CaM_N, CaMKP], kCaM2NP_off, [Ca2CaM_N_CaMKP]) # Ca2CaM_N + CaMKP <--> Ca2CaM_N_CaMKP
    add_rate!(rates, kCaM4_on * binding_To_PCaMK, [Ca4CaM, CaMKP], kCaM4P_off, [Ca4CaM_CaMKP]) # Ca4CaM + CaMKP <--> Ca4CaM_CaMKP
    add_rate!(rates, kCaM4_on, [Ca4CaM, CaMKOX], kCaM4_off, [Ca4CaM_CaMKOX]) # Ca4CaM + CaMKOX <--> Ca4CaM_CaMKOX
    add_rate!(rates, kCaM4_on * binding_To_PCaMK, [Ca4CaM, CaMKPOX], kCaM4P_off, [Ca4CaM_CaMKPOX]) # Ca4CaM + CaMKPOX <--> Ca4CaM_CaMKPOX

    ## Auto-phosphorylation of CaMKII
    # (Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, Ca4CaM_CaMKOX) <--> (Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, Ca4CaM_CaMKPOX)
    kphos = k_phosCaM * CaMKAct
    add_rate!(rates, kphos, Ca2CaM_C_CaMK, k_dephospho, Ca2CaM_C_CaMKP)
    add_rate!(rates, kphos, Ca2CaM_N_CaMK, k_dephospho, Ca2CaM_N_CaMKP)
    add_rate!(rates, kphos, Ca4CaM_CaMK, k_dephospho, Ca4CaM_CaMKP)
    add_rate!(rates, kphos, Ca4CaM_CaMKOX, k_dephospho, Ca4CaM_CaMKPOX)
    ## Second phosphorylation of CaMKII-P
    add_rate!(rates, k_P1_P2, CaMKP, k_P2_P1, CaMKP2) # CaMKP <--> CaMKP2
    ## Dephosphorylation of CaMKII-P
    # (CaMKP, CaMKPOX) --> (CaMK, CaMKOX)
    add_rate!(rates, k_dephospho, CaMKP, 0, CaMK)
    add_rate!(rates, k_dephospho, CaMKPOX, 0, CaMKOX)
    ## Redox reactions by ROS and reductases
    add_rate!(rates, k_B_OX * ROS, Ca4CaM_CaMK, k_OX_B, Ca4CaM_CaMKOX) # Ca4CaM_CaMK <--> Ca4CaM_CaMKOX
    add_rate!(rates, k_P_OXP * ROS, Ca4CaM_CaMKP, k_OXP_P, Ca4CaM_CaMKPOX) # Ca4CaM_CaMKP <--> Ca4CaM_CaMKPOX
    add_rate!(rates, k_OX_B, CaMKOX, 0, CaMK)
    add_rate!(rates, k_OXP_P, CaMKPOX, 0, CaMKP)

    rateeqs = [D(s) ~ rates[s] for s in sts]
    eqs = [
        CAMKII_T ~ CaMK + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaMKP + CaMKP2 + CaMKPOX + CaMKOX,
        CAM_T ~ CaM0 + Ca2CaM_C + Ca2CaM_N + Ca4CaM + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX,
        CaMKAct ~ (1 - CaMK / CAMKII_T)
    ]
    eqs_camkii = [rateeqs; eqs]
    return (; eqs_camkii, CaMKAct)
end

function get_camkii_sys(;
    Ca=0μM,
    ROS=0μM,
    binding_To_PCaMK=0,   ## 0.1 for T287D mutation
    name=:camkii_sys)
    @independent_variables t
    @unpack eqs_camkii = get_camkii_eqs(; Ca, ROS, binding_To_PCaMK)
    return System(eqs_camkii, t; name)
end

"""
Simplified CaMKII system with one-step activation of CaMK
"""
function get_camkii_simp_eqs(;
    Ca=0μM,
    ROS=0μM,
    binding_To_PCaMK=0, ## 0.1 for T287D mutation
    binding_To_OCaMK=0
)

    @independent_variables t
    D = Differential(t)

    @parameters begin
        r_CaMK = 79.5360Hz              ## Inverse of time scale of CaMK <--> CaMKB reaction (adjustable) # 3Hz
        kb_CaMKP = 3.0667Hz             ## CaMCa dissociation rate of CaMKP --> CaMKA (adjustable)       # 0.3Hz
        kfa2_CaMK = 0.2651              ## Maximal binding ratio by CaM-Ca2 (adjustable)
        kfa4_CaMK = 0.1636              ## Maximal binding ratio by CaM-Ca4 (adjustable)
        kfb_CaMK = 0.001                ## Basal binding by CaM (adjustable)
        kmCa2_CaMK = 0.7385μM           ## Half-saturation calcium concentration for CaM-Ca2 binding (adjustable)
        kmCa4_CaMK = 1.2515μM           ## Half-saturation calcium concentration for CaM-Ca4 binding (adjustable)
        kphos_CaMK = 1.8438Hz           ## Autophosphorylation rate ## 12.5Hz
        kdeph_CaMK = inv(12.0138second) ## Dephosphorylation rate ## inv(6 second)
        k_P1_P2 = inv(112.7635second)   ## Second autophosphorylation rate ## inv(60second)
        k_P2_P1 = inv(28.1950second)    ## Second dephosphorylation rate
        kox_CaMK = inv(45second) / 50μM ## 291Hz / mM   ## Oxidation by H2O2 (adjustable)
        krd_CaMK = inv(45second)        ## Reduction rate
        KActScale = 1.0                 ## CaMKII activation scaling factor (adjustable)
    end

    sts = @variables begin
        CaMKB(t) = 0.008728     ## Bound to CaMCa
        CaMKBOX(t) = 0          ## Bound to CaMCa, oxidized
        CaMKP(t) = 0.003916     ## Bound to CaMCa, autophosphorylated
        CaMKPOX(t) = 0          ## Bound to CaMCa, autophosphorylated, oxidized
        CaMKA(t) = 0.007833     ## Unbound, autophosphorylated
        CaMKA2(t) = 0           ## Unbound, phosphorylated, slow component
        CaMKAOX(t) = 0          ## Unbound, autophosphorylated, oxidized
        CaMKOX(t) = 0           ## Unbound, oxidized
    end

    @variables begin
        CaMK(t)             ## Inactive CaMKII
        CaMKAct(t)          ## Active CaMKII fraction
        CaMKBInf(t)         ## Steady-state CaMKB fraction
        fracCaMKPhos(t)     ## Total phosphorylated CaMKII fraction
        fracCaMKOx(t)       ## Total oxidized CaMKII fraction
    end

    rates = Dict()

    ## CaMK(OX) <--CaM,Ca--> CaMKB(OX)
    kf = r_CaMK * CaMKBInf
    kb = r_CaMK * (1 - CaMKBInf)
    add_rate!(rates, kf, CaMK, kb, CaMKB)
    add_rate!(rates, kf * binding_To_OCaMK, CaMKOX, kb, CaMKBOX)
    ## Auto-phosphorylation of CaMKII: B(OX) <--> P(OX)
    kphos = kphos_CaMK * CaMKAct
    add_rate!(rates, kphos, CaMKB, kdeph_CaMK, CaMKP)
    add_rate!(rates, kphos, CaMKBOX, kdeph_CaMK, CaMKPOX)
    ## CaMKA(OX) <-CaM,Ca-> CaMKP(OX)
    add_rate!(rates, kf * binding_To_PCaMK, CaMKA, kb_CaMKP, CaMKP)
    add_rate!(rates, kf * binding_To_PCaMK, CaMKAOX, kb_CaMKP, CaMKPOX)
    ## Second phosphorylation of Apo CaMKII-P
    # CaMKA <--> CaMKA2
    add_rate!(rates, k_P1_P2, CaMKA, k_P2_P1, CaMKA2)

    ## Dephosphorylation of Apo CaMK-P
    ## CaMKA(OX) --> CaMK(OX)
    add_rate!(rates, kdeph_CaMK, CaMKA, 0, CaMK)
    add_rate!(rates, kdeph_CaMK, CaMKAOX, 0, CaMKOX)

    ## Redox reactions by ROS and reductases
    ## CaMKX <--> CaMKOX
    add_rate!(rates, kox_CaMK * ROS, CaMKB, krd_CaMK, CaMKBOX)
    add_rate!(rates, kox_CaMK * ROS, CaMKP, krd_CaMK, CaMKPOX)
    add_rate!(rates, krd_CaMK, CaMKOX, 0, CaMK)
    add_rate!(rates, krd_CaMK, CaMKAOX, 0, CaMKA)

    rateeqs = [D(s) ~ rates[s] for s in sts]
    eqs = [
        CaMKBInf ~ kfb_CaMK + kfa2_CaMK * hil(Ca, kmCa2_CaMK, 2) + kfa4_CaMK * hil(Ca, kmCa4_CaMK, 4),
        fracCaMKPhos ~ CaMKP + CaMKPOX + CaMKA + CaMKA2 + CaMKAOX,
        fracCaMKOx ~ CaMKBOX + CaMKPOX + CaMKAOX + CaMKOX,
        CaMKAct ~ KActScale * (CaMKB + CaMKBOX + CaMKP + CaMKPOX + CaMKA + CaMKA2 + CaMKAOX + CaMKOX),
        1 ~ CaMK + CaMKB + CaMKBOX + CaMKP + CaMKPOX + CaMKA + CaMKA2 + CaMKAOX + CaMKOX,
    ]
    eqs_camkii = [rateeqs; eqs]
    return (; eqs_camkii, CaMKAct)
end

function get_camkii_simp_sys(;
    Ca=0μM,
    ROS=0μM,
    binding_To_PCaMK=0,
    binding_To_OCaMK=0,
    name=:camkii_sys)

    @independent_variables t
    @unpack eqs_camkii = get_camkii_simp_eqs(; Ca, ROS, binding_To_PCaMK, binding_To_OCaMK)
    return System(eqs_camkii, t; name)
end

"""
Diamond-shaped CaM binding to calcium with rapid equilibrium assumption
"""
function get_camkii_dia_eqs(;
    Ca=0μM,
    ROS=0μM,
    binding_To_PCaMK=0, ## 0.1 for T287D mutation
    binding_To_OCaMK=0,
)

    @independent_variables t
    D = Differential(t)

    ## Parameters from Chang et. al 2019
    @parameters begin
        CAM_T = 30μM            ## Total calmodulin Concentration
        CAMKII_T = 70μM         ## Total CaMKII Concentration
        ## Ca2+ binding to CaM
        ## C-lobe
        k_1C_on = 5Hz / μM      ## 1.2-9.6uM-1Hz
        k_1C_off = 50Hz         ## 10-70 Hz
        k_2C_on = 10Hz / μM     ## 5-25uM-1Hz
        k_2C_off = 10Hz         ## 8.5-10Hz
        Keq_CamC = k_1C_on * k_2C_on / (k_1C_off * k_2C_off) ## 0.1/uM^2
        ## N-lobe
        k_1N_on = 100Hz / μM    ## 25-260uM-1Hz
        k_1N_off = 2000Hz       ## 1000-4000 Hz
        k_2N_on = 200Hz / μM    ## 50-300uM-1Hz
        k_2N_off = 500Hz        ## 500-1000Hz
        Keq_CamN = k_1N_on * k_2N_on / (k_1N_off * k_2N_off) ## 0.02/uM^2

        ## Ca2+ binding to CaM-CAMKII (KCaM)
        ## C-lobe
        k_K1C_on = 44Hz / μM
        k_K1C_off = 33Hz
        k_K2C_on = 44Hz / μM
        k_K2C_off = 0.8Hz ## 0.49-4.9Hz
        Keq_KCamC = k_K1C_on * k_K2C_on / (k_K1C_off * k_K2C_off) ##  73/uM^2
        ## N-lobe
        k_K1N_on = 76Hz / μM
        k_K1N_off = 300Hz
        k_K2N_on = 76Hz / μM
        k_K2N_off = 20Hz ## 6-60Hz
        Keq_KCamN = k_K1N_on * k_K2N_on / (k_K1N_off * k_K2N_off) ## 0.96/uM^2
        ## CaM binding to CaMKII
        kCaM0_on = 3.8e-3Hz / μM ## Changed to Pepke's value from Chang's 3.8
        kCaM0_off = 5.5Hz
        kCaM2C_on = 0.5Hz / μM  # 0.92 μM-1Hz
        kCaM2C_off = 6.8Hz
        kCaM2N_on = 0.12Hz / μM
        kCaM2N_off = 1.7Hz
        kCaM4_on = 15Hz / μM  # 14-60 uM-1Hz
        kCaM4_off = 1.5Hz  # 1.1 - 2.3 Hz
        kb_CaMKP = 3Hz     ## CaMCa dissociation rate of CaMKP --> CaMKA (adjustable)       # 0.3Hz
        kCaM0P_off = kb_CaMKP
        kCaM2CP_off = kb_CaMKP
        kCaM2NP_off = kb_CaMKP
        kCaM4P_off = kb_CaMKP
        kphos_CaMK = 2Hz           ## Autophosphorylation rate ## 12.5Hz
        kdeph_CaMK = inv(12second) ## Dephosphorylation rate ## inv(6 second)
        k_P1_P2 = inv(60second)
        k_P2_P1 = inv(15second)

        ## Oxidation / reduction of Met
        krd_CaMK = inv(45second)        ## Reduction rate
        kox_CaMK = krd_CaMK / 50μM      ## 291Hz / mM   ## Oxidation by H2O2 (adjustable)
        KActScale = 1.0                 ## CaMKII activation scaling
    end

    sts = @variables begin
        CaMKB(t) = 0μM
        CaMKBOX(t) = 0μM
        CaMKP(t) = 0μM
        CaMKPOX(t) = 0μM
        CaMKA(t) = 0μM
        CaMKA2(t) = 0μM
        CaMKAOX(t) = 0μM
        CaMKOX(t) = 0μM
    end

    ## Dependent variables
    @variables begin
        CaMKAct(t)
        CaMK(t)
        ## calmodulin species
        CaM(t)
        CaM0(t)
        CaM2C(t)
        CaM2N(t)
        CaM4(t)
        ## bound CaMKII species
        CaMKB0(t)
        CaMKB2C(t)
        CaMKB2N(t)
        CaMKB4(t)
        CaMKBOX0(t)
        CaMKBOX2C(t)
        CaMKBOX2N(t)
        CaMKBOX4(t)
        ## Bound and phosphorylated CaMKII species
        CaMKP0(t)
        CaMKP2C(t)
        CaMKP2N(t)
        CaMKP4(t)
        CaMKPOX0(t)
        CaMKPOX2C(t)
        CaMKPOX2N(t)
        CaMKPOX4(t)
        ## CaM fractions
        fKCaM0(t)
        fKCaM2C(t)
        fKCaM2N(t)
        fKCaM4(t)
        fCaM0(t)
        fCaM2C(t)
        fCaM2N(t)
        fCaM4(t)
    end

    ## CaM fractions under rapid ca binding
    function _cam_fractions(ca, keq_c, keq_n)
        w0 = 1
        w2c = ca * ca * keq_c
        w2n = ca * ca * keq_n
        w4 = w2c * w2n
        wsum = w0 + w2c + w2n + w4
        f0 = w0 / wsum
        f2C = w2c / wsum
        f2N = w2n / wsum
        f4 = w4 / wsum
        return (f0, f2C, f2N, f4)
    end

    ## CaM fractions of CaM0, CaM2C, CaM2N, and CaM4
    f0, f2C, f2N, f4 = _cam_fractions(Ca, Keq_CamC, Keq_CamN)
    fK0, fK2C, fK2N, fK4 = _cam_fractions(Ca, Keq_KCamC, Keq_KCamN)

    rates = Dict()

    ## CaMK (OX) <--> CaMKB (OX)
    k0b = kCaM0_on * fCaM0 + kCaM2C_on * fCaM2C + kCaM2N_on * fCaM2N + kCaM4_on * fCaM4
    kb0 = kCaM0_off * fKCaM0 + kCaM2C_off * fKCaM2C + kCaM2N_off * fKCaM2N + kCaM4_off * fKCaM4
    add_rate!(rates, k0b, [CaM, CaMK], kb0, CaMKB)
    add_rate!(rates, binding_To_OCaMK * k0b, [CaM, CaMKOX], kb0, CaMKBOX)
    ## CaMKB (OX) <--> CaMKP (OX)
    kphos = kphos_CaMK * CaMKAct * KActScale
    add_rate!(rates, kphos, CaMKB, kdeph_CaMK, CaMKP)
    add_rate!(rates, kphos, CaMKBOX, kdeph_CaMK, CaMKPOX)
    ## CaMKP (OX) <--> CaMKA (OX)
    kfa = kCaM0P_off * fKCaM0 + kCaM2CP_off * fKCaM2C + kCaM2NP_off * fKCaM2N + kCaM4P_off * fKCaM4
    kaf = binding_To_PCaMK * k0b
    add_rate!(rates, kfa, CaMKP, kaf, [CaMKA, CaM])
    add_rate!(rates, kfa, CaMKPOX, kaf, [CaMKAOX, CaM])
    ## CaMKA <--> CaMKA2
    add_rate!(rates, k_P1_P2, CaMKA, k_P2_P1, CaMKA2)
    ## CaMKA (OX) --> CaMK (OX)
    add_rate!(rates, kdeph_CaMK, CaMKA, 0, CaMK)
    add_rate!(rates, kdeph_CaMK, CaMKAOX, 0, CaMKOX)
    ## CaMKB <--> CaMKBOX
    kox = kox_CaMK * ROS
    add_rate!(rates, kox, CaMKB, krd_CaMK, CaMKBOX)
    ## CaMKP <--> CaMKPOX
    add_rate!(rates, kox, CaMKP, krd_CaMK, CaMKPOX)
    ## CaMKAOX <--> CaMKA
    add_rate!(rates, 0, CaMKA, krd_CaMK, CaMKAOX)
    ## CaMK -->  CaMKOX
    add_rate!(rates, 0, CaMK, krd_CaMK, CaMKOX)

    rateeqs = [D(s) ~ rates[s] for s in sts]
    eqs = [
        CAMKII_T ~ CaMK + CaMKB + CaMKBOX + CaMKP + CaMKPOX + CaMKA + CaMKA2 + CaMKAOX + CaMKOX,
        CaMKAct ~ (CAMKII_T - CaMK - CaMKB0) / CAMKII_T,
        CAM_T ~ CaM + CaMKB + CaMKBOX + CaMKP + CaMKPOX,
        fCaM0 ~ f0,
        fCaM2C ~ f2C,
        fCaM2N ~ f2N,
        fCaM4 ~ f4,
        fKCaM0 ~ fK0,
        fKCaM2C ~ fK2C,
        fKCaM2N ~ fK2N,
        fKCaM4 ~ fK4,
        CaM0 ~ fCaM0 * CaM,
        CaM2C ~ fCaM2C * CaM,
        CaM2N ~ fCaM2N * CaM,
        CaM4 ~ fCaM4 * CaM,
        CaMKB0 ~ fKCaM0 * CaMKB,
        CaMKB2C ~ fKCaM2C * CaMKB,
        CaMKB2N ~ fKCaM2N * CaMKB,
        CaMKB4 ~ fKCaM4 * CaMKB,
        CaMKBOX0 ~ fKCaM0 * CaMKBOX,
        CaMKBOX2C ~ fKCaM2C * CaMKBOX,
        CaMKBOX2N ~ fKCaM2N * CaMKBOX,
        CaMKBOX4 ~ fKCaM4 * CaMKBOX,
        CaMKP0 ~ fKCaM0 * CaMKP,
        CaMKP2C ~ fKCaM2C * CaMKP,
        CaMKP2N ~ fKCaM2N * CaMKP,
        CaMKP4 ~ fKCaM4 * CaMKP,
        CaMKPOX0 ~ fKCaM0 * CaMKPOX,
        CaMKPOX2C ~ fKCaM2C * CaMKPOX,
        CaMKPOX2N ~ fKCaM2N * CaMKPOX,
        CaMKPOX4 ~ fKCaM4 * CaMKPOX,
    ]
    eqs_camkii = [eqs; rateeqs]
    return (; eqs_camkii, CaMKAct)
end

function get_camkii_dia_sys(;
    Ca=0μM,
    ROS=0μM,
    binding_To_PCaMK=0,
    binding_To_OCaMK=0,
    name=:camkii_sys)

    @independent_variables t
    @unpack eqs_camkii = get_camkii_dia_eqs(; Ca, ROS, binding_To_PCaMK, binding_To_OCaMK)
    return System(eqs_camkii, t; name)
end
