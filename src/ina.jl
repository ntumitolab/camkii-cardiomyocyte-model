"Fast sodium current (INa) and background sodium current"
function get_ina_sys(nai, nao, vm; name=:inasys)
    @parameters begin
        gNa = 12.8mSμF⁻¹
        gNab = 0.0026mSμF⁻¹
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
    mα = 0.32/ms / 0.1 * exprel((v + 47.13mV) * inv(-10mV))
    mβ = 0.08/ms * exp(v * inv(-11mV))
    hα = 0.135/ms * exp((v + 80mV) * inv(-6.8mV))
    hβ = 7.6923/ms / (1 + exp((v + 10.66mV) * inv(-11.1mV)))
    jα = max((-127140 * exp(0.2444/mV * v) - 3.474e-5 * exp(-0.04391/mV * v)) * (v + 37.78mV) / (1 + exp(0.311/mV * (v + 79.23mV))) / (mV * ms), 0)
    jβ = 0.3/ms * exp(-2.535e-7/mV * v) / (1 + exp((v + 32mV) * inv(-10mV)))

    eqs = [
        INa ~ gNa * i_Nam^3 * i_Nah * i_Naj * (vm - E_Na),
        D(i_Nam) ~ mα - i_Nam * (mα + mβ),
        D(i_Nah) ~ hα - i_Nah * (hα + hβ),
        D(i_Naj) ~ inv(ms) * (jα - i_Naj * (jα + jβ)),
        E_Na ~ nernst(nao, nai),
        INab ~ gNab * (vm - E_Na),
    ]

    return ODESystem(eqs, t; name)
end
