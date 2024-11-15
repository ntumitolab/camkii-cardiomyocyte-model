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
        k_phosCaM = 12.6Hz # 30Hz
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
    # (Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, Ca4CaM_CaMKOX) --> (Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, Ca4CaM_CaMKPOX)
    k = k_phosCaM * CaMKAct
    add_rate!(rates, k, [Ca2CaM_C_CaMK], 0, [Ca2CaM_C_CaMKP])
    add_rate!(rates, k, [Ca2CaM_N_CaMK], 0, [Ca2CaM_N_CaMKP])
    add_rate!(rates, k, [Ca4CaM_CaMK], 0, [Ca4CaM_CaMKP])
    add_rate!(rates, k, [Ca4CaM_CaMKOX], 0, [Ca4CaM_CaMKPOX])
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

"CaMKII system with ROS activation and fast calcium binding"
function get_camkii_fast_ca_binding_sys(Ca=0μM;
    ROS=0μM,
    binding_To_PCaMK=0,
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
        k_2N_off = 500Hz        ## 500->1000.s-1

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
        kCaM2C_on = 0.5Hz / μM # 0.92 μM-1s-1
        kCaM2C_off = 6.8Hz
        kCaM2N_on = 0.12Hz / μM
        kCaM2N_off = 1.7Hz
        kCaM4_on = 15Hz / μM # 14-60 uM-1s-1
        kCaM4_off = 1.5Hz # 1.1 - 2.3 s-1
        kCaM0P_on = kCaM0_on * binding_To_PCaMK
        kCaM2CP_on = kCaM2C_on * binding_To_PCaMK
        kCaM2NP_on = kCaM2N_on * binding_To_PCaMK
        kCaM4P_on = kCaM4_on * binding_To_PCaMK
        kCaM0P_off = inv(3second)
        kCaM2CP_off = inv(3second)
        kCaM2NP_off = inv(3second)
        kCaM4P_off = inv(3second)
        k_phosCaM = 12.6Hz # 30Hz
        k_dephospho = inv(6second)
        k_P1_P2 = inv(60second)
        k_P2_P1 = inv(15second)

        ## Oxidation / reduction
        k_BOX = 291Hz / mM
        k_POXP = 291Hz / mM
        k_OXB = inv(45second)
        k_OXPP = inv(45second)
    end

    sts = @variables begin
        KCaM(t) = 0
        PCaM(t) = 0
        OCaM(t) = 0
        OPCaM(t) = 0
        CaMKP(t) = 0
        CaMKP2(t) = 0
        CaMKPOX(t) = 0
        CaMKOX(t) = 0
    end

    conservedvars = @variables begin
        CaM(t)
        CaMK(t)
    end

    rates = merge(Dict(sts .=> Num(0)), Dict(conservedvars .=> Num(0)))

    # Observables
    @variables begin
        CaMKAct(t)
        CaM0(t)
        CaM2C(t)
        CaM2N(t)
        CaM4(t)
        CaM0_CaMK(t)
        CaM2C_CaMK(t)
        CaM2N_CaMK(t)
        CaM4_CaMK(t)
        CaM0_CaMKP(t)
        CaM2C_CaMKP(t)
        CaM2N_CaMKP(t)
        CaM4_CaMKP(t)
        CaM0_CaMKOX(t)
        CaM2C_CaMKOX(t)
        CaM2N_CaMKOX(t)
        CaM4_CaMKOX(t)
        CaM0_CaMKPOX(t)
        CaM2C_CaMKPOX(t)
        CaM2N_CaMKPOX(t)
        CaM4_CaMKPOX(t)
    end
    eqs = [
        CaMKAct ~ (1 - CaMK / CAMKII_T),
        CAMKII_T ~ CaMK + KCaM + PCaM + OCaM + OPCaM + CaMKP + CaMKP2 + CaMKPOX + CaMKOX,
        CAM_T ~ CaM + KCaM + PCaM + OCaM + OPCaM,
    ]

    "Ca binding/unbinding reaction equilibrium"
    _ca_eq(ca, k1on, k1off, k2on, k2off) = ca^2 * (k1on * k2on) / (k1off * k2off)

    ## Two Ca2+ ions bind to C (high affinity) or N (low affinity)-lobe of CaM
    keqc = _ca_eq(Ca, k_1C_on, k_1C_off, k_2C_on, k_2C_off)
    keqn = _ca_eq(Ca, k_1N_on, k_1N_off, k_2N_on, k_2N_off)
    den = (1 + keqc) * (1 + keqn)
    push!(eqs, CaM0 ~ CaM / den, CaM2C ~ CaM * keqc / den, CaM2N ~ CaM * keqn / den, CaM4 ~ CaM * keqc * keqn / den)

    ## Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII(P) complex
    keqc = _ca_eq(Ca, k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off)
    keqn = _ca_eq(Ca, k_K1N_on, k_K1N_off, k_K2N_on, k_K2N_off)
    den = (1 + keqc) * (1 + keqn)
    push!(eqs, CaM0_CaMK ~ KCaM / den, CaM2C_CaMK ~ KCaM * keqc / den, CaM2N_CaMK ~ KCaM * keqn / den, CaM4_CaMK ~ KCaM * keqc * keqn / den)
    push!(eqs, CaM0_CaMKP ~ PCaM / den, CaM2C_CaMKP ~ PCaM * keqc / den, CaM2N_CaMKP ~ PCaM * keqn / den, CaM4_CaMKP ~ PCaM * keqc * keqn / den)
    push!(eqs, CaM0_CaMKOX ~ OCaM / den, CaM2C_CaMKOX ~ OCaM * keqc / den, CaM2N_CaMKOX ~ OCaM * keqn / den, CaM4_CaMKOX ~ OCaM * keqc * keqn / den)
    push!(eqs, CaM0_CaMKPOX ~ OPCaM / den, CaM2C_CaMKPOX ~ OPCaM * keqc / den, CaM2N_CaMKPOX ~ OPCaM * keqn / den, CaM4_CaMKPOX ~ OPCaM * keqc * keqn / den)

    ## CaM + CaMK <--> KCaM
    vf1 = CaMK * (kCaM0_on * CaM0 + kCaM2C_on * CaM2C + kCaM2N_on * CaM2N + kCaM4_on * CaM4)
    vr1 = kCaM0_off * CaM0_CaMK + kCaM2C_off * CaM2C_CaMK + kCaM2N_off * CaM2N_CaMK + kCaM4_off * CaM4_CaMK
    add_raw_rate!(rates, vf1 - vr1, [CaMK, CaM], [KCaM])

    ## CaM + CaMKP <--> PCaM
    vf2 = CaMKP * (kCaM0P_on * CaM0 + kCaM2CP_on * CaM2C + kCaM2NP_on * CaM2N + kCaM4P_on * CaM4)
    vr2 = kCaM0P_off * CaM0_CaMKP + kCaM2CP_off * CaM2C_CaMKP + kCaM2NP_off * CaM2N_CaMKP + kCaM4P_off * CaM4_CaMKP
    add_raw_rate!(rates, vf2 - vr2, [CaMKP, CaM], [PCaM])

    ## CaM + CaMKOX <--> OCaM
    vf3 = CaMKOX * (kCaM0_on * CaM0 + kCaM2C_on * CaM2C + kCaM2N_on * CaM2N + kCaM4_on * CaM4)
    vr3 = kCaM0_off * CaM0_CaMKOX + kCaM2C_off * CaM2C_CaMKOX + kCaM2N_off * CaM2N_CaMKOX + kCaM4_off * CaM4_CaMKOX
    add_raw_rate!(rates, vf3 - vr3, [CaMKOX, CaM], [OCaM])

    ## CaM + CaMKPOX <--> OPCaM
    vf4 = CaMKPOX * (kCaM0P_on * CaM0 + kCaM2CP_on * CaM2C + kCaM2NP_on * CaM2N + kCaM4P_on * CaM4)
    vr4 = kCaM0P_off * CaM0_CaMKPOX + kCaM2CP_off * CaM2C_CaMKPOX + kCaM2NP_off * CaM2N_CaMKPOX + kCaM4P_off * CaM4_CaMKPOX
    add_raw_rate!(rates, vf4 - vr4, [CaMKPOX, CaM], [OPCaM])

    ## Auto-phosphorylation of CaMKII
    k = k_phosCaM * CaMKAct
    add_rate!(rates, k, [KCaM], 0, [PCaM])
    add_rate!(rates, k, [OCaM], 0, [OPCaM])

    ## Second phosphorylation of CaMKII-P
    add_rate!(rates, k_P1_P2, [CaMKP], k_P2_P1, [CaMKP2]) # CaMKP <--> CaMKP2

    ## Dephosphorylation of CaMKII-P: (CaMKP, CaMKPOX) --> (CaMK, CaMKOX)
    add_rate!(rates, k_dephospho, [CaMKP], 0, [CaMK])
    add_rate!(rates, k_dephospho, [CaMKPOX], 0, [CaMKOX])

    ## Redox reactions by ROS and reductases
    add_rate!(rates, k_BOX * ROS, [KCaM], k_OXB, [OCaM])
    add_rate!(rates, k_POXP * ROS, [PCaM], k_OXPP, [OPCaM])
    add_rate!(rates, k_OXB, [CaMKOX], 0, [CaMK])
    add_rate!(rates, k_OXPP, [CaMKPOX], 0, [CaMKP])

    rateeqs = [D(s) ~ rates[s] for s in sts]
    sys = ODESystem([rateeqs; eqs], t; name)
    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end
