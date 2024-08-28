# CaMKII system with ROS activation

function get_camkii_sys(;
    Ca=0μM, ROS=0μM,
    binding_To_PCaMK=0.1,
    decay_CaM=3,
    phospho_rate=1Hz,
    phosphatase=1Hz,
    name=:camkii_sys
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

        ## Ca2+ binding to CaM-CAMKII(K)
        ## C-lobe
        k_K1C_on = 44Hz / μM
        k_K1C_off = 33Hz
        k_K2C_on = 44Hz / μM
        k_K2C_off = 0.8Hz ## 0.49-4.9Hz
        ## N-lobe
        k_K1N_on = 76Hz / μM
        k_K1N_off = 300Hz
        k_K2N_on = 76Hz / μM
        k_K2N_off = 20Hz ## 6-60-1

        ## CaM binding to CaMKII
        kCaM0_on = 3.8e-3Hz / μM
        kCaM2C_on = 0.92Hz / μM
        kCaM2N_on = 0.12Hz / μM
        kCaM4_on = 30Hz / μM
        kCaM0_off = 5.5Hz
        kCaM2C_off = 6.8Hz
        kCaM2N_off = 1.7Hz
        kCaM4_off = 1.5Hz
        kCaM0P_on = 3.8e-3Hz / μM * binding_To_PCaMK
        kCaM2CP_on = 0.92Hz / μM * binding_To_PCaMK
        kCaM2NP_on = 0.12Hz / μM * binding_To_PCaMK
        kCaM4P_on = 30Hz / μM * binding_To_PCaMK
        kCaM0P_off = 1Hz / decay_CaM
        kCaM2CP_off = 1Hz / decay_CaM
        kCaM2NP_off = 1Hz / decay_CaM
        kCaM4P_off = 1Hz / decay_CaM
        k_phosCaM = 30 * phospho_rate
        k_dephospho = (1 / 6) * phosphatase
        k_P1_P2 = 1 / 60Hz
        k_P2_P1 = (1 / 6) * 0.25Hz

        ## Oxidation / reduction
        k_OXPOX = 0.03Hz / μM
        k_POXP = 0.291Hz / μM
        k_OXB = 2.23e-2Hz
        k_OXPP = 2.23e-2Hz
    end

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
        CaM0(t) ## Conserved
        CaMK(t) ## Conserved
        CaMKAct(t) ## Observable
    end

    rates = Dict(sts .=> t - t)

    "Ca binding/unbinding reaction rates"
    function _ca_cam(ca, k1on, k1off, k2on, k2off)
        fcaon = hil(ca * k2on, k1off)
        fcaoff = 1 - fcaon
        return (ca * k1on * fcaon, k2off * fcaoff)
    end

    ## Two Ca2+ ions bind to C (high affinity) or N (low affinity)-lobe of CaM
    kon, koff = _ca_cam(Ca, k_1C_on, k_1C_off, k_2C_on, k_2C_off)
    add_rate!(rates, kon, [CaM0], koff, [Ca2CaM_C]) # CaM0 <--> Ca2CaM_C
    add_rate!(rates, kon, [Ca2CaM_C], koff, [Ca4CaM]) # Ca2CaM_C <--> Ca4CaM
    kon, koff = _ca_cam(Ca, k_1N_on, k_1N_off, k_2N_on, k_2N_off)
    add_rate!(rates, kon, [CaM0], koff, [Ca2CaM_N]) # CaM0 <--> Ca2CaM_N
    add_rate!(rates, kon, [Ca2CaM_N], koff, [Ca4CaM]) # Ca2CaM_N <--> Ca4CaM
    ## Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII(P) complex
    kon, koff = _ca_cam(Ca, k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off)
    add_rate!(rates, kon, [CaM0_CaMK], koff, [Ca2CaM_C_CaMK]) # CaM0_CaMK <--> Ca2CaM_C_CaMK
    add_rate!(rates, kon, [Ca2CaM_C_CaMK], koff, [Ca4CaM_CaMK]) # Ca2CaM_C_CaMK <--> Ca4CaM_CaMK
    add_rate!(rates, kon, [CaM0_CaMKP], koff, [Ca2CaM_C_CaMKP]) # CaM0_CaMKP <--> Ca2CaM_C_CaMKP
    add_rate!(rates, kon, [Ca2CaM_C_CaMKP], koff, [Ca4CaM_CaMKP]) # Ca2CaM_C_CaMKP <--> Ca4CaM_CaMKP
    kon, koff = _ca_cam(Ca, k_K1N_on, k_K1N_off, k_K2N_on, k_K2N_off)
    add_rate!(rates, kon, [CaM0_CaMK], koff, [Ca2CaM_N_CaMK]) # CaM0_CaMK <--> Ca2CaM_N_CaMK
    add_rate!(rates, kon, [Ca2CaM_N_CaMK], koff, [Ca4CaM_CaMK]) # Ca2CaM_N_CaMK <--> Ca4CaM_CaMK
    add_rate!(rates, kon, [CaM0_CaMKP], koff, [Ca2CaM_N_CaMKP]) # CaM0_CaMKP <--> Ca2CaM_N_CaMKP
    add_rate!(rates, kon, [Ca2CaM_N_CaMKP], koff, [Ca4CaM_CaMKP]) # Ca2CaM_N_CaMKP <--> Ca4CaM_CaMKP

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
    add_rate!(rates, k_OXPOX * ROS, [Ca4CaM_CaMK], k_OXB, [Ca4CaM_CaMKOX]) # Ca4CaM_CaMK <--> Ca4CaM_CaMKOX
    add_rate!(rates, k_POXP * ROS, [Ca4CaM_CaMKP], k_OXPP, [Ca4CaM_CaMKPOX]) # Ca4CaM_CaMKP <--> Ca4CaM_CaMKPOX
    add_rate!(rates, k_OXB, [CaMKOX], 0, [CaMK])
    add_rate!(rates, k_OXPP, [CaMKPOX], 0, [CaMKP])

    eqs = [
        CaMKAct ~ (1 - CaMK / CAMKII_T),
        CAMKII_T ~ CaMK + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaMKP + CaMKP2 + CaMKPOX + CaMKOX,
        CAM_T ~ CaM0 + Ca2CaM_C + Ca2CaM_N + Ca4CaM + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX,
    ]

    rateeqs = [D(s) ~ rates[s] for s in sts if s ∉ (CaMKAct, CAMKII_T, CAM_T)]

    return ODESystem([eqs; rateeqs], t; name)
end
