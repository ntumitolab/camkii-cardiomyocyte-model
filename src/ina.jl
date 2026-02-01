function get_ina_eqs(nai, nao, vm)
    @parameters begin
        gNa = 12.8mSμF
        gNab = 0.0026mSμF
    end
    @variables begin
        INa(t)
        i_Nam(t) = 0.0327    # Fast Na gating (activation)
        i_Nah(t) = 0.9909    # Fast Na gating (inactivation)
        i_Naj(t) = 0.9941    # Fast Na gating (slow inactivation)
        INab(t)
        E_Na(t)
    end

    v = vm
    mα = 0.32 / 0.1 * exprel(-0.1 * (v + 47.13))
    mβ = 0.08 * exp(-v / 11)
    hα = 0.135 * exp(-(v + 80) / 6.8)
    hβ = 7.6923 / (1 + exp(-(v + 10.66) / 11.1))
    jα = max((-127140 * exp(0.2444v) - 3.474e-5 * exp(-0.04391v)) * (v + 37.78) / (1 + exp(0.311 * (v + 79.23))), 0)
    jβ = 0.3 * exp(-2.535e-7v) / (1 + exp(-0.1 * (v + 32)))

    eqs_ina = [
        INa ~ gNa * i_Nam^3 * i_Nah * i_Naj * (vm - E_Na),
        D(i_Nam) ~ inv(ms) * (mα - i_Nam * (mα + mβ)),
        D(i_Nah) ~ inv(ms) * (hα - i_Nah * (hα + hβ)),
        D(i_Naj) ~ inv(ms) * (jα - i_Naj * (jα + jβ)),
        E_Na ~ nernst(nao, nai),
        INab ~ gNab * (vm - E_Na),
    ]

    return (; eqs_ina, INa, INab)
end

"Fast sodium current (INa) and background sodium current"
function get_ina_sys(nai, nao, vm; name=:inasys)
    @unpack eqs_ina = get_ina_eqs(nai, nao, vm)
    return System(eqs_ina, t; name)
end
