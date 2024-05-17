# CaMKII system with ROS activation
using Catalyst
using ModelingToolkit

function get_camkii_eqs(
    Ca, ROS=0.0μM;
    cam_total=30μM, ## Total calmodulin Concentration
    camkii_total=70μM, ## Total CaMKII Concentration
    binding_To_PCaMK=0.1,
    decay_CaM=3,
    phospho_rate=1Hz,
    phosphatase=1Hz)

    @parameters begin
        CAM_T = cam_total
        CAMKII_T = camkii_total

        ## Ca2+ binding to CaM-CAMKII
        ## C lobe
        k_1C_on = 5Hz / μM ## 1.2-9.6uM-1s-1
        k_1C_off = 50Hz ## 10-70 s-1
        k_2C_on = 10Hz / μM ## 5-25uM-1s-1.
        k_2C_off = 10Hz ## 8.5-10s-1.

        ## N-lobe
        k_1N_on = 100Hz / μM ## 25-260uM-1s-1
        k_1N_off = 2000Hz ## 1000-4000 s-1
        k_2N_on = 200Hz / μM ## 50-300uM-1s-1.
        k_2N_off = 500Hz ## 500->1000.s-1

        ## Ca2+ binding to CaM-CAMKII(K)
        ## C-lobe
        k_K1C_on = 44Hz / μM
        k_K1C_off = 33Hz
        k_K2C_on = 44Hz / μM
        k_K2C_off = 0.8Hz ## 0.49-4.9s-1

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
        k_POXP = 0.291Hz / μM #uM-1s-1  0.291
        k_OXB = 2.23e-2Hz
        k_OXPP = 2.23e-2Hz
    end

    @variables begin
        t
        CaM0(t)
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
        CaMK(t)
        CaMKP(t) = 0
        CaMKP2(t) = 0
        CaMKPOX(t) = 0
        CaMKOX(t) = 0
        CaMKII_act(t)
    end

    D = Differential(t)
    _konoff(k1on, k1off, k2on, k2off) = (Ca^2 * k1on * k2on / (Ca * k2on + k1off), k1off * k2off / (Ca * k2on + k1off))

    ## Two Ca2+ ions bind to C (high affinity) or N (low affinity)-lobe of CaM
    kcon, kcoff = _konoff(k_1C_on, k_1C_off, k_2C_on, k_2C_off)
    knon, knoff = _konoff(k_1N_on, k_1N_off, k_2CNon, k_2N_off)

    # Record rates
    #===
    rates = Dict(Ca=>Num(0))
    rates[Ca] += 1
    rates[Ca] -= dCa
    ===#

    v_CaM0_to_Ca2CaM_C = kcon * CaM0 - kcoff * Ca2CaM_C
    v_Ca2CaM_C_to_Ca4CaM = kcon * Ca2CaM_C - kcoff * Ca4CaM
    v_CaM0_to_Ca2CaM_N = knon * CaM0 - knoff * Ca2CaM_N
    v_Ca2CaM_N_to_Ca4CaM = knon * Ca2CaM_N - knoff * Ca4CaM

    ## Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII(P) complex
    kkcon, kkcoff = _konoff(k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off)
    kknon, kknoff = _konoff(k_K1N_on, k_K1N_off, k_K2CNon, k_K2N_off)

    v_CaM0_CaMK_to_Ca2CaM_C_CaMK = kkcon * CaM0_CaMK - kkcoff * Ca2CaM_C_CaMK
    v_Ca2CaM_C_CaMK_to_Ca4CaM_CaMK = kkcon * Ca2CaM_C_CaMK - kkcoff * Ca4CaM_CaMK
    v_CaM0_CaMKP_to_Ca2CaM_C_CaMKP = kkcon * CaM0_CaMKP - kkcoff * Ca2CaM_C_CaMKP
    v_Ca2CaM_C_CaMKP_to_Ca4CaM_CaMKP = kkcon * Ca2CaM_C_CaMKP - kkcoff * Ca4CaM_CaMKP
    v_CaM0_CaMK_to_Ca2CaM_N_CaMK = kknon * CaM0_CaMK - kknoff * Ca2CaM_N_CaMK
    v_Ca2CaM_N_CaMK_to_Ca4CaM_CaMK = kknon * Ca2CaM_N_CaMK - kknoff * Ca4CaM_CaMK
    v_CaM0_CaMKP_to_Ca2CaM_N_CaMKP = kknon * CaM0_CaMKP - kknoff * Ca2CaM_N_CaMKP
    v_Ca2CaM_N_CaMKP_to_Ca4CaM_CaMKP = kknon * Ca2CaM_N_CaMKP - kknoff * Ca4CaM_CaMKP

    ## CaM binding to CaMKII / CaMkII-P / CaMkII-POX / CaMkII-OX
    v_CaM0_CaMK = kCaM0_on * CaM0 * CaMK - kCaM0_off * CaM0_CaMK
    v_Ca2CaM_C_CaMK = kCaM2C_on * Ca2CaM_C * CaMK - kCaM2C_off * Ca2CaM_C_CaMK
    v_Ca2CaM_N_CaMK = kCaM2N_on * Ca2CaM_N * CaMK - kCaM2N_off * Ca2CaM_N_CaMK
    v_Ca4CaM_CaMK = kCaM4_on * Ca4CaM * CaMK - kCaM4_off * Ca4CaM_CaMK
    v_CaM0_CaMKP = kCaM0P_on * CaM0 * CaMKP - kCaM0P_off * CaM0_CaMKP
    v_Ca2CaM_C_CaMKP = kCaM2CP_on * Ca2CaM_C * CaMKP - kCaM2CP_off * Ca2CaM_C_CaMKP
    v_Ca2CaM_N_CaMKP = kCaM2NP_on * Ca2CaM_N * CaMKP - kCaM2NP_off * Ca2CaM_N_CaMKP
    v_Ca4CaM_CaMKP = kCaM4P_on * Ca4CaM * CaMKP - kCaM4P_off * Ca4CaM_CaMKP
    v_Ca4CaM_CaMKOX = kCaM4_on * Ca4CaM * CaMKOX - kCaM4_off * Ca4CaM_CaMKOX
    v_Ca4CaM_CaMKPOX = kCaM4P_on * Ca4CaM * CaMKPOX - kCaM4P_off * Ca4CaM_CaMKPOX

    # Phosphorylation of CaMKII
    kphos = k_phosCaM * CaMKII_act
    vphos_Ca2CaM_C_CaMK = kphos * Ca2CaM_C_CaMK
    vphos_Ca2CaM_N_CaMK = kphos * Ca2CaM_N_CaMK

    eqs = [
        CAM_T ~ CaM0 + Ca2CaM_C + Ca2CaM_N + Ca4CaM + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX,
        CAMKII_T ~ CaMK + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaMKP + CaMKP2 + CaMKPOX + CaMKOX,
        CaMKII_act ~ (CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX + CaMKP + CaMKP2 + CaMKPOX + CaMKOX) / CAMKII_T,
        # D(CaM0)  ## Conserved
        # D(CaMK)  ## Conserved
        D(Ca2CaM_C) ~ v_CaM0_to_Ca2CaM_C - v_Ca2CaM_C_to_Ca4CaM - v_Ca2CaM_C_CaMK - v_Ca2CaM_C_CaMKP,
        D(Ca2CaM_N) ~ v_CaM0_to_Ca2CaM_N - v_Ca2CaM_N_to_Ca4CaM - v_Ca2CaM_N_CaMK - v_Ca2CaM_N_CaMKP,
        D(Ca4CaM) ~ v_Ca2CaM_C_to_Ca4CaM + v_Ca2CaM_N_to_Ca4CaM - v_Ca4CaM_CaMK - v_Ca4CaM_CaMKP - v_Ca4CaM_CaMKOX - v_Ca4CaM_CaMKPOX,
        D(CaM0_CaMK) ~ -v_CaM0_CaMK_to_Ca2CaM_C_CaMK - v_CaM0_CaMK_to_Ca2CaM_N_CaMK + v_CaM0_CaMK,
        D(Ca2CaM_C_CaMK) ~ v_CaM0_CaMK_to_Ca2CaM_C_CaMK - v_Ca2CaM_C_CaMK_to_Ca4CaM_CaMK + v_Ca2CaM_C_CaMK - vphos_Ca2CaM_C_CaMK,
        D(Ca2CaM_N_CaMK) ~ v_CaM0_CaMK_to_Ca2CaM_N_CaMK - v_Ca2CaM_N_CaMK_to_Ca4CaM_CaMK + v_Ca2CaM_N_CaMK - vphos_Ca2CaM_N_CaMK,
        D(Ca4CaM_CaMK) ~ v_Ca2CaM_C_CaMK_to_Ca4CaM_CaMK + v_Ca2CaM_N_CaMK_to_Ca4CaM_CaMK + v_Ca4CaM_CaMK,
        D(Ca2CaM_C_CaMKP) ~ v_CaM0_CaMKP_to_Ca2CaM_C_CaMKP - v_Ca2CaM_C_CaMKP_to_Ca4CaM_CaMKP + v_Ca2CaM_C_CaMKP + vphos_Ca2CaM_C_CaMK,
        D(Ca2CaM_N_CaMKP) ~ v_CaM0_CaMKP_to_Ca2CaM_N_CaMKP - v_Ca2CaM_N_CaMKP_to_Ca4CaM_CaMKP + v_Ca2CaM_N_CaMKP + vphos_Ca2CaM_N_CaMK,
        D(Ca4CaM_CaMKP) ~ v_Ca2CaM_C_CaMKP_to_Ca4CaM_CaMKP + v_Ca2CaM_N_CaMKP_to_Ca4CaM_CaMKP + v_Ca4CaM_CaMKP,
        D(CaMKP) ~ -v_CaM0_CaMKP - v_Ca2CaM_C_CaMKP - v_Ca2CaM_N_CaMKP - v_Ca4CaM_CaMKP,
        D(CaM0_CaMKP) ~ -v_CaM0_CaMKP_to_Ca2CaM_C_CaMKP - v_CaM0_CaMKP_to_Ca2CaM_N_CaMKP + v_CaM0_CaMKP,
        D(CaMKOX) ~ -v_Ca4CaM_CaMKOX,
        D(CaMKPOX) ~ -v_Ca4CaM_CaMKPOX,
        D(Ca4CaM_CaMKOX) ~ v_Ca4CaM_CaMKOX,
        D(Ca4CaM_CaMKPOX) ~ v_Ca4CaM_CaMKPOX,
    ]

    return eqs
end

function get_camkii_rn(Ca, ROS=0.0μM;
    cam_total=30μM, ## Total calmodulin Concentration
    camkii_total=70μM, ## Total CaMKII Concentration
    binding_To_PCaMK=0.1,
    decay_CaM=3,
    phospho_rate=1Hz,
    phosphatase=1Hz,
)

    rn = @reaction_network begin
        ## Two Ca2+ ions bind to C (high affinity) or N (low affinity)-lobe of CaM
        (mm($Ca * k_2C_on, $Ca * k_1C_on, k_1C_off), mmr($Ca * k_2C_on, k_2C_off, k_1C_off)), (CaM0, Ca2CaM_C) <--> (Ca2CaM_C, Ca4CaM)
        (mm($Ca * k_2N_on, $Ca * k_1N_on, k_1N_off), mmr($Ca * k_2N_on, k_2N_off, k_1N_off)), (CaM0, Ca2CaM_N) <--> (Ca2CaM_N, Ca4CaM)

        ## Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII(P) complex
        (mm($Ca * k_K2C_on, $Ca * k_K1C_on, k_K1C_off), mmr($Ca * k_K2C_on, k_K2C_off, k_K1C_off)), (CaM0_CaMK, Ca2CaM_C_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP) <--> (Ca2CaM_C_CaMK, Ca4CaM_CaMK, Ca2CaM_C_CaMKP, Ca4CaM_CaMKP)
        (mm($Ca * k_K2N_on, $Ca * k_K1N_on, k_K1N_off), mmr($Ca * k_K2N_on, k_K2N_off, k_K1N_off)), (CaM0_CaMK, Ca2CaM_N_CaMK, CaM0_CaMKP, Ca2CaM_N_CaMKP) <--> (Ca2CaM_N_CaMK, Ca4CaM_CaMK, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP)

        ## CaM binding to CaMKII / CaMkII-P / CaMkII-POX / CaMkII-OX
        (kCaM0_on, kCaM0_off), CaM0 + CaMK <--> CaM0_CaMK
        (kCaM2C_on, kCaM2C_off), Ca2CaM_C + CaMK <--> Ca2CaM_C_CaMK
        (kCaM2N_on, kCaM2N_off), Ca2CaM_N + CaMK <--> Ca2CaM_N_CaMK
        (kCaM4_on, kCaM4_off), Ca4CaM + CaMK <--> Ca4CaM_CaMK
        (kCaM0P_on, kCaM0P_off), CaM0 + CaMKP <--> CaM0_CaMKP
        (kCaM2CP_on, kCaM2CP_off), Ca2CaM_C + CaMKP <--> Ca2CaM_C_CaMKP
        (kCaM2NP_on, kCaM2NP_off), Ca2CaM_N + CaMKP <--> Ca2CaM_N_CaMKP
        (kCaM4P_on, kCaM4P_off), Ca4CaM + CaMKP <--> Ca4CaM_CaMKP
        (kCaM4P_on, kCaM4P_off), Ca4CaM + CaMKPOX <--> Ca4CaM_CaMKPOX
        (kCaM4_on, kCaM4_off), Ca4CaM + CaMKOX <--> Ca4CaM_CaMKOX

        ## Phosphorylation of CaMKII
        k_phosCaM * (CaMKP + CaMKP2 + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP + Ca4CaM_CaMKOX + Ca4CaM_CaMKPOX) / $camkii_total, (Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, Ca4CaM_CaMKOX) --> (Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, Ca4CaM_CaMKPOX)

        ## Second phosphorylation of CaMKII-P
        (k_P1_P2, k_P2_P1), CaMKP <--> CaMKP2

        ## Dephosphorylation of CaMKII-P
        k_dephospho, CaMKP --> CaMK

        ## Oxidation by ROS
        $ROS * k_OXPOX, Ca4CaM_CaMK --> Ca4CaM_CaMKOX
        $ROS * k_POXP, Ca4CaM_CaMKP --> Ca4CaM_CaMKPOX

        ## Reduction
        k_OXPP, (CaMKPOX, Ca4CaM_CaMKP) --> (CaMKP, Ca4CaM_CaMKPOX)
        k_OXB, (CaMKOX, Ca4CaM_CaMK) --> (CaMK, Ca4CaM_CaMKOX)
    end

    setdefaults!(rn, [
        ## ICS
        :CaM0 => cam_total,
        :Ca2CaM_C => 0μM,
        :Ca2CaM_N => 0μM,
        :Ca4CaM => 0μM,
        :CaM0_CaMK => 0μM,
        :Ca2CaM_C_CaMK => 0μM,
        :Ca2CaM_N_CaMK => 0μM,
        :Ca4CaM_CaMK => 0μM,
        :CaM0_CaMKP => 0μM,
        :Ca2CaM_C_CaMKP => 0μM,
        :Ca2CaM_N_CaMKP => 0μM,
        :Ca4CaM_CaMKP => 0μM,
        :CaMK => camkii_total,
        :CaMKP => 0μM,
        :CaMKP2 => 0μM,
        :Ca4CaM_CaMKOX => 0μM,
        :Ca4CaM_CaMKPOX => 0μM,
        :CaMKPOX => 0μM,
        :CaMKOX => 0μM,

        ## Ca2+ binding to CaM-CAMKII
        ## C lobe
        :k_1C_on => 5Hz / μM, ## 1.2-9.6uM-1s-1
        :k_1C_off => 50Hz, ## 10-70 s-1
        :k_2C_on => 10Hz / μM, ## 5-25uM-1s-1.
        :k_2C_off => 10Hz, ## 8.5-10s-1.

        ## N-lobe
        :k_1N_on => 100Hz / μM, ## 25-260uM-1s-1
        :k_1N_off => 2000Hz, ## 1000-4000 s-1
        :k_2N_on => 200Hz / μM, ## 50-300uM-1s-1.
        :k_2N_off => 500Hz, ## 500->1000.s-1

        ## Ca2+ binding to CaM-CAMKII(K)
        ## C-lobe
        :k_K1C_on => 44Hz / μM,
        :k_K1C_off => 33Hz,
        :k_K2C_on => 44Hz / μM,
        :k_K2C_off => 0.8Hz, ## 0.49-4.9s-1

        ## N-lobe
        :k_K1N_on => 76Hz / μM,
        :k_K1N_off => 300Hz,
        :k_K2N_on => 76Hz / μM,
        :k_K2N_off => 20Hz, ## 6-60-1

        ## CaM binding to CaMKII
        :kCaM0_on => 3.8e-3Hz / μM,
        :kCaM2C_on => 0.92Hz / μM,
        :kCaM2N_on => 0.12Hz / μM,
        :kCaM4_on => 30Hz / μM,
        :kCaM0_off => 5.5Hz,
        :kCaM2C_off => 6.8Hz,
        :kCaM2N_off => 1.7Hz,
        :kCaM4_off => 1.5Hz,
        :kCaM0P_on => 3.8e-3Hz / μM * binding_To_PCaMK,
        :kCaM2CP_on => 0.92Hz / μM * binding_To_PCaMK,
        :kCaM2NP_on => 0.12Hz / μM * binding_To_PCaMK,
        :kCaM4P_on => 30Hz / μM * binding_To_PCaMK,
        :kCaM0P_off => 1Hz / decay_CaM,
        :kCaM2CP_off => 1Hz / decay_CaM,
        :kCaM2NP_off => 1Hz / decay_CaM,
        :kCaM4P_off => 1Hz / decay_CaM,
        :k_phosCaM => 30 * phospho_rate,
        :k_dephospho => (1 / 6) * phosphatase,
        :k_P1_P2 => 1 / 60Hz,
        :k_P2_P1 => (1 / 6) * 0.25Hz,

        ## Oxidation / reduction
        :k_OXPOX => 0.03Hz / μM,
        :k_POXP => 0.291Hz / μM, #uM-1s-1  0.291
        :k_OXB => 2.23e-2Hz,
        :k_OXPP => 2.23e-2Hz,
    ])

    return rn
end
