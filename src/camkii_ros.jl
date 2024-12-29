"""
CaMKII system with ROS activation

CaMKII model: "Mechanisms of Ca2+/calmodulin-dependent kinase II activation in single dendritic spines" (Chang et al., 2019, Nature Communications); https://doi.org/10.1038/s41467-019-10694-z

ROS activation model: Oxidized Calmodulin Kinase II Regulates Conduction Following Myocardial Infarction: A Computational Analysis (Christensen et al. 2009); https://doi.org/10.1371/journal.pcbi.1000583
"""
function get_camkii_sys(Ca=0μM;
    ROS=0μM,
    binding_To_PCaMK=0,   ## 0.1 for T287D mutation
    name=:camkii_sys,
    simplify=false
)
    @parameters begin
        CAM_T = 30μM            ## Total calmodulin Concentration
        CAMKII_T = 70μM         ## Total CaMKII Concentration
        k_1C_on = 5Hz / μM      ## 1.2-9.6uM-1s-1
        k_1C_off = 50Hz         ## 10-70 s-1
        k_2C_on = 10Hz / μM     ## 5-25uM-1s-1.
        k_2C_off = 10Hz         ## 8.5-10s-1.
        ## N-lobe
        k_1N_on = 100Hz / μM    ## 25-260uM-1s-1
        k_1N_off = 2000Hz       ## 1000-4000 s-1
        k_2N_on = 200Hz / μM    ## 50-300uM-1s-1.
        k_2N_off = 500Hz        ## 500-1000.s-1

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
        k_K2N_off = 20Hz ## 6-60s-1

        ## CaM binding to CaMKII
        kCaM0_on = 3.8Hz / mM
        kCaM0_off = 5.5Hz
        kCaM2C_on = 0.5Hz / μM  # 0.92 μM-1s-1
        kCaM2C_off = 6.8Hz
        kCaM2N_on = 0.12Hz / μM
        kCaM2N_off = 1.7Hz
        kCaM4_on = 15Hz / μM  # 14-60 uM-1s-1
        kCaM4_off = 1.5Hz  # 1.1 - 2.3 s-1
        kCaM0P_on = kCaM0_on * binding_To_PCaMK
        kCaM2CP_on = kCaM2C_on * binding_To_PCaMK
        kCaM2NP_on = kCaM2N_on * binding_To_PCaMK
        kCaM4P_on = kCaM4_on * binding_To_PCaMK
        kCaM0P_off = inv(3second)
        kCaM2CP_off = inv(3second)
        kCaM2NP_off = inv(3second)
        kCaM4P_off = inv(3second)
        k_phosCaM = 5Hz # 30Hz
        k_dephospho = inv(6second)
        k_P1_P2 = inv(60second)
        k_P2_P1 = inv(15second)

        ## Oxidation / reduction of Met
        k_BOX = 291Hz / mM
        k_POXP = 291Hz / mM
        k_OXB = inv(45second)
        k_OXPP = inv(45second)
    end

    sts = @variables begin
        Ca2CaM_C(t) = 84.33nM
        Ca2CaM_N(t) = 8.578nM
        Ca4CaM(t) = 0.02nM
        CaM0_CaMK(t) = 1μM
        Ca2CaM_C_CaMK(t) = 398.68nM
        Ca2CaM_N_CaMK(t) = 16.6nM
        Ca4CaM_CaMK(t) = 6.66nM
        CaM0_CaMKP(t) = 623.72nM
        Ca2CaM_C_CaMKP(t) = 1.01μM
        Ca2CaM_N_CaMKP(t) = 10.75nM
        Ca4CaM_CaMKP(t) = 15.81nM
        Ca4CaM_CaMKOX(t) = 0mM
        Ca4CaM_CaMKPOX(t) = 0mM
        CaMKP(t) = 3.3μM
        CaMKP2(t) = 831.43nM
        CaMKPOX(t) = 0mM
        CaMKOX(t) = 0mM
    end

    conservedvars = @variables begin
        CaM0(t)
        CaMK(t)
    end

    conservedeqs = [
        CAMKII_T ~ CaMK + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaMKP + CaMKP2 + CaMKPOX + CaMKOX,
        CAM_T ~ CaM0 + Ca2CaM_C + Ca2CaM_N + Ca4CaM + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX,
    ]

    rates = merge(Dict(sts .=> Num(0)), Dict(conservedvars .=> Num(0)))

    # Observables
    @variables CaMKAct(t)
    obseqs = [CaMKAct ~ (1 - CaMK / CAMKII_T)]

    "Ca binding/unbinding reaction rates"
    function _ca_cam(ca, k1on, k1off, k2on, k2off)
        fcaon = hil(ca * k2on, k1off)
        fcaoff = 1 - fcaon
        return (ca * k1on * fcaon, k2off * fcaoff)
    end

    ## Two Ca2+ ions bind to C (high affinity) or N (low affinity)-lobe of CaM
    kon, koff = _ca_cam(Ca, k_1C_on, k_1C_off, k_2C_on, k_2C_off)
    add_rate!(rates, kon, [CaM0], koff, [Ca2CaM_C]) # CaM0 <--> Ca2CaM_C
    add_rate!(rates, kon, [Ca2CaM_N], koff, [Ca4CaM]) # Ca2CaM_C <--> Ca4CaM
    kon, koff = _ca_cam(Ca, k_1N_on, k_1N_off, k_2N_on, k_2N_off)
    add_rate!(rates, kon, [CaM0], koff, [Ca2CaM_N]) # CaM0 <--> Ca2CaM_N
    add_rate!(rates, kon, [Ca2CaM_C], koff, [Ca4CaM]) # Ca2CaM_N <--> Ca4CaM
    ## Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII(P) complex
    kon, koff = _ca_cam(Ca, k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off)
    add_rate!(rates, kon, [CaM0_CaMK], koff, [Ca2CaM_C_CaMK]) # CaM0_CaMK <--> Ca2CaM_C_CaMK
    add_rate!(rates, kon, [Ca2CaM_N_CaMK], koff, [Ca4CaM_CaMK]) # Ca2CaM_C_CaMK <--> Ca4CaM_CaMK
    add_rate!(rates, kon, [CaM0_CaMKP], koff, [Ca2CaM_C_CaMKP]) # CaM0_CaMKP <--> Ca2CaM_C_CaMKP
    add_rate!(rates, kon, [Ca2CaM_N_CaMKP], koff, [Ca4CaM_CaMKP]) # Ca2CaM_C_CaMKP <--> Ca4CaM_CaMKP
    kon, koff = _ca_cam(Ca, k_K1N_on, k_K1N_off, k_K2N_on, k_K2N_off)
    add_rate!(rates, kon, [CaM0_CaMK], koff, [Ca2CaM_N_CaMK]) # CaM0_CaMK <--> Ca2CaM_N_CaMK
    add_rate!(rates, kon, [Ca2CaM_C_CaMK], koff, [Ca4CaM_CaMK]) # Ca2CaM_N_CaMK <--> Ca4CaM_CaMK
    add_rate!(rates, kon, [CaM0_CaMKP], koff, [Ca2CaM_N_CaMKP]) # CaM0_CaMKP <--> Ca2CaM_N_CaMKP
    add_rate!(rates, kon, [Ca2CaM_C_CaMKP], koff, [Ca4CaM_CaMKP]) # Ca2CaM_N_CaMKP <--> Ca4CaM_CaMKP

    ## CaM binding to CaMKII / CaMkII-P / CaMkII-POX / CaMkII-OX
    add_rate!(rates, kCaM0_on, [CaM0, CaMK], kCaM0_off, [CaM0_CaMK]) # CaM0 + CaMK <--> CaM0_CaMK
    add_rate!(rates, kCaM2C_on, [Ca2CaM_C, CaMK], kCaM2C_off, [Ca2CaM_C_CaMK]) # Ca2CaM_C + CaMK <--> Ca2CaM_C_CaMK
    add_rate!(rates, kCaM2N_on, [Ca2CaM_N, CaMK], kCaM2N_off, [Ca2CaM_N_CaMK]) # Ca2CaM_N + CaMK <--> Ca2CaM_N_CaMK
    add_rate!(rates, kCaM4_on, [Ca4CaM, CaMK], kCaM4_off, [Ca4CaM_CaMK]) # Ca4CaM + CaMK <--> Ca4CaM_CaMK
    add_rate!(rates, kCaM0P_on, [CaM0, CaMKP], kCaM0P_off, [CaM0_CaMKP])  # CaM0 + CaMKP <--> CaM0_CaMKP
    add_rate!(rates, kCaM2CP_on, [Ca2CaM_C, CaMKP], kCaM2CP_off, [Ca2CaM_C_CaMKP]) # Ca2CaM_C + CaMKP <--> Ca2CaM_C_CaMKP
    add_rate!(rates, kCaM2NP_on, [Ca2CaM_N, CaMKP], kCaM2NP_off, [Ca2CaM_N_CaMKP]) # Ca2CaM_N + CaMKP <--> Ca2CaM_N_CaMKP
    add_rate!(rates, kCaM4P_on, [Ca4CaM, CaMKP], kCaM4P_off, [Ca4CaM_CaMKP]) # Ca4CaM + CaMKP <--> Ca4CaM_CaMKP
    add_rate!(rates, kCaM4_on, [Ca4CaM, CaMKOX], kCaM4_off, [Ca4CaM_CaMKOX]) # Ca4CaM + CaMKOX <--> Ca4CaM_CaMKOX
    add_rate!(rates, kCaM4P_on, [Ca4CaM, CaMKPOX], kCaM4P_off, [Ca4CaM_CaMKPOX]) # Ca4CaM + CaMKPOX <--> Ca4CaM_CaMKPOX

    ## Auto-phosphorylation of CaMKII
    # (Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, Ca4CaM_CaMKOX) <--> (Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, Ca4CaM_CaMKPOX)
    kphos = k_phosCaM * CaMKAct
    add_rate!(rates, kphos, [Ca2CaM_C_CaMK], k_dephospho, [Ca2CaM_C_CaMKP])
    add_rate!(rates, kphos, [Ca2CaM_N_CaMK], k_dephospho, [Ca2CaM_N_CaMKP])
    add_rate!(rates, kphos, [Ca4CaM_CaMK], k_dephospho, [Ca4CaM_CaMKP])
    add_rate!(rates, kphos, [Ca4CaM_CaMKOX], k_dephospho, [Ca4CaM_CaMKPOX])
    ## Second phosphorylation of CaMKII-P
    add_rate!(rates, k_P1_P2, [CaMKP], k_P2_P1, [CaMKP2]) # CaMKP <--> CaMKP2
    ## Dephosphorylation of CaMKII-P
    # (CaMKP, CaMKPOX) --> (CaMK, CaMKOX)
    add_rate!(rates, k_dephospho, [CaMKP], 0, [CaMK])
    add_rate!(rates, k_dephospho, [CaMKPOX], 0, [CaMKOX])
    ## Redox reactions by ROS and reductases
    add_rate!(rates, k_BOX * ROS, [Ca4CaM_CaMK], k_OXB, [Ca4CaM_CaMKOX]) # Ca4CaM_CaMK <--> Ca4CaM_CaMKOX
    add_rate!(rates, k_POXP * ROS, [Ca4CaM_CaMKP], k_OXPP, [Ca4CaM_CaMKPOX]) # Ca4CaM_CaMKP <--> Ca4CaM_CaMKPOX
    add_rate!(rates, k_OXB, [CaMKOX], 0, [CaMK])
    add_rate!(rates, k_OXPP, [CaMKPOX], 0, [CaMKP])

    rateeqs = [D(s) ~ rates[s] for s in sts]
    sys = ODESystem([rateeqs; conservedeqs; obseqs], t; name)
    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end

"Simplified CaMKII system"
function get_camkii_simp_sys(Ca=0μM;
    ROS=0μM,
    binding_To_PCaMK=0,
    binding_To_OCaMK=0,
    name=:camkii_sys,
    simplify=false
)

    @parameters begin
        r_CaMK = 3Hz                ## Inverse of time scale of CaMK <--> CaMKB reaction (adjustable)
        kb_CaMKP = inv(3second)     ## Dissociation rate of CaMKP --> CaMKA (adjustable)
        kfa_CaMK = 0.4393           ## Activated forward rate CaMK --> CaMKB (adjustable)
        kfb_CaMK = 0.0056           ## Basal forward rate of CaMK --> CaMKB (adjustable)
        kmCa_CaMK = 0.9716μM        ## Half-activation concentration of calcium
        nCa_CaMK = 2.293            ## Hill coefficient of calcium
        kphos_CaMK = 5Hz           ## Phosphorylation rate (originally 30Hz)
        kdeph_CaMK = inv(6second)   ## Dephosphorylation rate
        k_P1_P2 = inv(60second)     ## Second autophosphorylation rate
        k_P2_P1 = inv(15second)     ## Second dephosphorylation rate
        kox_CaMK = 291Hz / mM       ## Oxidation rate
        krd_CaMK = inv(45second)    ## Reduction rate
    end

    #==
    States: bound/unbound, autophosphorylated/dephosphorylated, oxidized/reduced
    ==#
    sts = @variables begin
        CaMKB(t) = 0.008728     ## Bound to CaMCa
        CaMKBOX(t) = 0          ## Bound to CaMCa, oxidized
        CaMKP(t) = 0.003916     ## Bound to CaMCa, autophosphorylated
        CaMKPOX(t) = 0          ## Bound to CaMCa, autophosphorylated, oxidized
        CaMKA(t) = 0.007833     ## Unbound, autophosphorylated
        CaMKA2(t) = 0.001958    ## Unbound, doublely phosphorylated
        CaMKAOX(t) = 0          ## Unbound, autophosphorylated, oxidized
        CaMKOX(t) = 0           ## Unbound, oxidized
    end

    conservedvars = @variables begin
        CaMK(t)             ## Inactive CaMKII
    end

    rates = merge(Dict(sts .=> Num(0)), Dict(conservedvars .=> Num(0)))

    # Observables
    @variables CaMKAct(t)
    eqs = [
        CaMKAct ~ 1 - CaMK,
        1 ~ CaMK + CaMKB + CaMKBOX + CaMKP + CaMKPOX + CaMKA + CaMKA2 + CaMKAOX + CaMKOX,
    ]

    ## CaMK(OX) <--CaM,Ca--> CaMKB(OX)
    ## Calc the steady-state bound fraction (binf) for forward / backward reaction rates
    binf = kfa_CaMK * hil(Ca, kmCa_CaMK, nCa_CaMK) + kfb_CaMK
    kf = r_CaMK * binf
    kb = r_CaMK * (1 - binf)
    add_rate!(rates, kf, [CaMK], kb, [CaMKB])
    add_rate!(rates, kf * binding_To_OCaMK, [CaMKOX], kb, [CaMKBOX])

    ## CaMKA(OX) <-CaM,Ca-> CaMKP(OX)
    add_rate!(rates, kf * binding_To_PCaMK, [CaMKA], kb_CaMKP, [CaMKP])
    add_rate!(rates, kf * binding_To_PCaMK, [CaMKAOX], kb_CaMKP, [CaMKPOX])

    ## Auto-phosphorylation of CaMKII
    ## B(OX) <--> P(OX)
    kphos = kphos_CaMK * CaMKAct
    add_rate!(rates, kphos, [CaMKB], kdeph_CaMK, [CaMKP])
    add_rate!(rates, kphos, [CaMKBOX], kdeph_CaMK, [CaMKPOX])

    ## Second phosphorylation of Apo CaMKII-P
    # CaMKA <--> CaMKA2
    add_rate!(rates, k_P1_P2, [CaMKA], k_P2_P1, [CaMKA2])

    ## Dephosphorylation of Apo CaMK-P
    ## CaMKA(OX) --> CaMK(OX)
    add_rate!(rates, kdeph_CaMK, [CaMKA], 0, [CaMK])
    add_rate!(rates, kdeph_CaMK, [CaMKAOX], 0, [CaMKOX])

    ## Redox reactions by ROS and reductases
    add_rate!(rates, kox_CaMK * ROS, [CaMKB], krd_CaMK, [CaMKBOX])
    add_rate!(rates, kox_CaMK * ROS, [CaMKP], krd_CaMK, [CaMKPOX])
    add_rate!(rates, krd_CaMK, [CaMKOX], 0, [CaMK])
    add_rate!(rates, krd_CaMK, [CaMKAOX], 0, [CaMKA])

    rateeqs = [D(s) ~ rates[s] for s in sts]
    sys = ODESystem([rateeqs; eqs], t; name)
    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end
