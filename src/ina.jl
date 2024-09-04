# Sodium current
"Fast sodium current (INa)"
function get_ina_sys(vm, E_Na; name=:inasys)
    @parameters gNa = 35mS / cm^2
    @variables begin
        INa(t)
        Naminf(t)
        Nahinf(t)
        Najinf(t)
        Nataum(t)
        Natauh(t)
        Natauj(t)
        i_Nam(t) = 0.0250
        i_Nah(t) = 0.22242
        i_Naj(t) = 0.19081
    end

    V = vm * Volt / mV # Convert voltage to mV
    NatauhHI = 0.4537ms * expit((V + 10.66) / 11.1)
    NatauhLOW = 3.49ms / (0.135 * exp((V + 80) / -6.8) + 3.56 * exp(0.079V) + 3.1e5 * exp(0.35V))
    NataujHI = 11.63ms * (1 + exp(-0.1 * (V + 32))) / exp(-2.535e-7V)
    NataujLOW = 3.49ms / ((V + 37.78) / (1 + exp(0.311 * (V + 79.23))) * (-127140 * exp(0.2444V) - 3.474e-5 * exp(-0.04391V)) + 0.1212 * exp(-0.01052V) / (1 + exp(-0.1378 * (V + 40.14))))

    eqs = [
        Naminf ~ expit((V + 45) / 6.5),
        Nahinf ~ expit(-(V + 76.1) / 6.07),
        Najinf ~ Nahinf,
        Nataum ~ 1.36ms / (3.2 * exprel(-0.1 * (V + 47.13)) + 0.08 * exp(-V / 11)),
        Natauh ~ ifelse(V >= -40, NatauhHI, NatauhLOW),
        Natauj ~ ifelse(V >= -40, NataujHI, NataujLOW),
        INa ~ gNa * i_Nam^3 * i_Nah * i_Naj * (vm - E_Na),
        D(i_Nam) ~ (Naminf - i_Nam)/Nataum,
        D(i_Nah) ~ (Nahinf - i_Nah)/Natauh,
        D(i_Naj) ~ (Najinf - i_Naj)/Natauj,
    ]
    return ODESystem(eqs, t; name)
end
