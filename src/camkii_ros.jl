# CaMKII system with ROS activation
using Catalyst

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
        :k_1C_on => 5Hz/μM, ## 1.2-9.6uM-1s-1
        :k_1C_off => 50Hz, ## 10-70 s-1
        :k_2C_on => 10Hz/μM, ## 5-25uM-1s-1.
        :k_2C_off => 10Hz, ## 8.5-10s-1.

        ## N-lobe
        :k_1N_on => 100Hz/μM, ## 25-260uM-1s-1
        :k_1N_off => 2000Hz, ## 1000-4000 s-1
        :k_2N_on => 200Hz/μM, ## 50-300uM-1s-1.
        :k_2N_off => 500Hz, ## 500->1000.s-1

        ## Ca2+ binding to CaM-CAMKII(K)
        ## C-lobe
        :k_K1C_on => 44Hz/μM,
        :k_K1C_off => 33Hz,
        :k_K2C_on => 44Hz/μM,
        :k_K2C_off => 0.8Hz, ## 0.49-4.9s-1

        ## N-lobe
        :k_K1N_on => 76Hz/μM,
        :k_K1N_off => 300Hz,
        :k_K2N_on => 76Hz/μM,
        :k_K2N_off => 20Hz, ## 6-60-1

        ## CaM binding to CaMKII
        :kCaM0_on => 3.8e-3Hz/μM,
        :kCaM2C_on => 0.92Hz/μM,
        :kCaM2N_on => 0.12Hz/μM,
        :kCaM4_on => 30Hz/μM,
        :kCaM0_off => 5.5Hz,
        :kCaM2C_off => 6.8Hz,
        :kCaM2N_off => 1.7Hz,
        :kCaM4_off => 1.5Hz,
        :kCaM0P_on => 3.8e-3Hz/μM * binding_To_PCaMK,
        :kCaM2CP_on => 0.92Hz/μM * binding_To_PCaMK,
        :kCaM2NP_on => 0.12Hz/μM * binding_To_PCaMK,
        :kCaM4P_on => 30Hz/μM * binding_To_PCaMK,
        :kCaM0P_off => 1Hz / decay_CaM,
        :kCaM2CP_off => 1Hz / decay_CaM,
        :kCaM2NP_off => 1Hz / decay_CaM,
        :kCaM4P_off => 1Hz / decay_CaM,
        :k_phosCaM => 30 * phospho_rate,
        :k_dephospho => (1 / 6)* phosphatase,
        :k_P1_P2 => 1 / 60Hz,
        :k_P2_P1 => (1 / 6) * 0.25Hz,

        ## Oxidation / reduction
        :k_OXPOX => 0.03Hz/μM,
        :k_POXP => 0.291Hz/μM, #uM-1s-1  0.291
        :k_OXB => 2.23e-2Hz,
        :k_OXPP => 2.23e-2Hz,
    ])

    return rn
end
