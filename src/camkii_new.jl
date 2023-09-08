using Catalyst
using ModelingToolkit

# :CaMT => 30e-6 # Total calmodulin concentration.

function make_camkii_rn(;
    decay_CaM = 3000,
    binding_To_PCaMK = 0.1,
    phospho_rate = 1,
    phosphatase = 1
)
    # Chemical Reactions
    rn = @reaction_network begin
        @variables Ca(t)

        # (d*50e-9, d), 0 <--> Ca

        # Two Ca2+ ions bind to C or N-lobe.
        (mm(k_2C_on * Ca, k_1C_on * Ca, k_1C_off), mmr(k_2C_on * Ca, k_2C_off, k_1C_off)), CaM0 <--> Ca2CaM_C
        (mm(k_2N_on * Ca, k_1N_on * Ca, k_1N_off), mmr(k_2N_on * Ca, k_2N_off, k_1N_off)), CaM0 <--> Ca2CaM_N
        (mm(k_2C_on * Ca, k_1C_on * Ca, k_1C_off), mmr(k_2C_on * Ca, k_2C_off, k_1C_off)), Ca2CaM_C <--> Ca4CaM
        (mm(k_2N_on * Ca, k_1N_on * Ca, k_1N_off), mmr(k_2N_on * Ca, k_2N_off, k_1N_off)), Ca2CaM_N <--> Ca4CaM

        # Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII complex.
        (mm(k_K2C_on * Ca, k_K1C_on * Ca, k_K1C_off), mmr(k_K2C_on * Ca, k_K2C_off, k_K1C_off)), (CaM0_CaMK, CaM0_CaMKP) <--> (Ca2CaM_C_CaMK, Ca2CaM_C_CaMKP)
        (mm(k_K2N_on * Ca, k_K1N_on * Ca, k_K1N_off), mmr(k_K2N_on * Ca, k_K2N_off, k_K1N_off)), (CaM0_CaMK, CaM0_CaMKP) <--> (Ca2CaM_N_CaMK, Ca2CaM_N_CaMKP)
        (mm(k_K2C_on * Ca, k_K1C_on * Ca, k_K1C_off), mmr(k_K2C_on * Ca, k_K2C_off, k_K1C_off)), (Ca2CaM_C_CaMK, Ca2CaM_C_CaMKP) <--> (Ca4CaM_CaMK, Ca4CaM_CaMKP)
        (mm(k_K2N_on * Ca, k_K1N_on * Ca, k_K1N_off), mmr(k_K2N_on * Ca, k_K2N_off, k_K1N_off)), (Ca2CaM_N_CaMK, Ca2CaM_N_CaMKP) <--> (Ca4CaM_CaMK, Ca4CaM_CaMKP)

        ##  Binding of CaM to CaMKII or CaMII-P
        (kCaM0_on, kCaM0_off), CaM0 + CaMK <--> CaM0_CaMK
        (kCaM2C_on, kCaM2C_off), Ca2CaM_C + CaMK <--> Ca2CaM_C_CaMK
        (kCaM2N_on, kCaM2N_off), Ca2CaM_N + CaMK <--> Ca2CaM_N_CaMK
        (kCaM4_on, kCaM4_off), Ca4CaM + CaMK <--> Ca4CaM_CaMK
        (kCaM0P_on, kCaM0P_off), CaM0 + CaMKP <--> CaM0_CaMKP
        (kCaM2CP_on, kCaM2CP_off), Ca2CaM_C + CaMKP <--> Ca2CaM_C_CaMKP
        (kCaM2NP_on, kCaM2NP_off), Ca2CaM_N + CaMKP <--> Ca2CaM_N_CaMKP
        (kCaM4P_on, kCaM4P_off), Ca4CaM + CaMKP <--> Ca4CaM_CaMKP

        ##  Phosphorylation CaMXCaMKII -> CaMXCaMKIIP.
        k_phosCaM * (CaMKP + CaMKP2 + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, Ca2CaM_C_CaMK --> Ca2CaM_C_CaMKP
        k_phosCaM * (CaMKP + CaMKP2 + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, Ca2CaM_N_CaMK --> Ca2CaM_N_CaMKP
        k_phosCaM * (CaMKP + CaMKP2 + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, Ca4CaM_CaMK --> Ca4CaM_CaMKP

        ##  Dephosphorylation CaMKP -> CaMK
        k_dephospho, CaMKP --> CaMK
        ##  Second phosphorylation state (P2) CaMKP <-> CaMKP2
        (k_P1_P2, k_P2_P1), CaMKP <--> CaMKP2
    end

    setdefaults!(rn, [
        :CaMKII_T => 70e-6 # Total CaMKII concentration.
        :k_phosCaM => 30e-3 * phospho_rate
        :k_dephospho => 1e-3/6 * phosphatase
        :k_1C_on => 5000      # 1/ms/uM
        :k_1C_off => 50e-3    # 1/ms
        :k_2C_on => 10000     # 1/ms/uM
        :k_2C_off => 10e-3    # 1/ms
        :k_1N_on => 100e3     # 1/ms/uM
        :k_1N_off => 2000e-3  # 1/ms
        :k_2N_on => 200e3     # 1/ms/uM
        :k_2N_off => 500e-3   # 1/ms
        :k_K1C_on => 44e3     # 1/ms/uM
        :k_K1C_off => 33e-3   # 1/ms
        :k_K2C_on => 44e3    # 1/ms/uM
        :k_K2C_off => 0.8e-3  # 1/ms
        :k_K1N_on => 76e3     # 1/ms/uM
        :k_K1N_off => 300e-3  # 1/ms
        :k_K2N_on => 76e3     # 1/ms/uM
        :k_K2N_off => 20e-3   # 1/ms
        :kCaM0_on => 3.8
        :kCaM0_off => 5.5e-3
        :kCaM2C_on => 0.92e3
        :kCaM2C_off => 6.8e-3
        :kCaM2N_on => 0.12e3
        :kCaM2N_off => 1.7e-3
        :kCaM4_on => 30e3
        :kCaM4_off => 1.5e-3
        :kCaM0P_on => 3.8 * binding_To_PCaMK
        :kCaM0P_off => inv(decay_CaM)
        :kCaM2CP_on => 0.92e3 * binding_To_PCaMK
        :kCaM2CP_off => inv(decay_CaM)
        :kCaM2NP_on => 0.12e3 * binding_To_PCaMK
        :kCaM2NP_off => inv(decay_CaM)
        :kCaM4P_on => 30e3 * binding_To_PCaMK
        :kCaM4P_off => inv(decay_CaM)
        :k_P1_P2 => 1e-3 / 60
        :k_P2_P1 => 1e-3 / 6 * 0.25
    ])

    return rn
end
