using Catalyst
using DifferentialEquations
using ModelingToolkit
using Plots

function make_camodel()
    ##  Two Ca2+ ions bind to C or N-lobe.
    caMon(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off) = k_1C_on*Ca^2*k_2C_on/(k_1C_off+k_2C_on*Ca)
    caMoff(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off) = k_1C_off*k_2C_off/(k_1C_off+k_2C_on*Ca)

    ca_model = @reaction_network begin
        (dCa*CaResting, Ca), 0 <--> Ca

        ##  Two Ca2+ ions bind to C or N-lobe of CaM
        caMon(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), CaM0 --> Ca2CaM_C
        caMoff(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), CaM0 <-- Ca2CaM_C
        caMon(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), CaM0 --> Ca2CaM_N
        caMoff(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), CaM0 <-- Ca2CaM_N
        caMon(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), Ca2CaM_C --> Ca4CaM
        caMoff(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), Ca2CaM_C <-- Ca4CaM
        caMon(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), Ca2CaM_N --> Ca4CaM
        caMoff(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), Ca2CaM_N <-- Ca4CaM

        ##  Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII complex
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMK --> Ca2CaM_C_CaMK
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMK <-- Ca2CaM_C_CaMK
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMK --> Ca2CaM_N_CaMK
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMK <-- Ca2CaM_N_CaMK
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMK --> Ca4CaM_CaMK
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMK <-- Ca4CaM_CaMK
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMK --> Ca4CaM_CaMK
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMK <-- Ca4CaM_CaMK

        ##  Binding of Ca to CaM-CaMKIIP.
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMKP --> Ca2CaM_C_CaMKP
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMKP <-- Ca2CaM_C_CaMKP
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMKP --> Ca2CaM_N_CaMKP
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMKP <-- Ca2CaM_N_CaMKP
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMKP --> Ca4CaM_CaMKP
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMKP <-- Ca4CaM_CaMKP
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMKP --> Ca4CaM_CaMKP
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMKP <-- Ca4CaM_CaMKP
    end

    return ca_model
end

make_camodel()
