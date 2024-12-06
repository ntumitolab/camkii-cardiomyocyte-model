"Fast sodium current (INa) and background sodium current"
function get_ina_sys(nai, nao, vm; name=:inasys)
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

    ishigh = (v >= -40)

    hαhi = 0
    hβhi = 7.6923 * expit((v + 10.66) / 11.1)
    hαlo = 0.135 * exp(-(v + 80) / 6.8)
    hβlo = 3.56 * exp(0.079v) + 3.1E5 * exp(0.35v)
    hα = ishigh * hαhi + (1 - ishigh) * hαlo
    hβ = ishigh * hβhi + (1 - ishigh) * hβlo

    jαhi = 0
    jβhi = 0.3 * exp(-2.535e-7v) * expit(0.1 * (v + 32))
    jαlo = (-127140 * exp(0.2444v) - 3.474e-5 * exp(-0.04391v)) * (v + 37.78) * expit(-0.311 * (v + 79.23))
    jβlo = 0.212 * exp(-0.01052v) * expit(0.1378 * (v + 40.14))
    jα = ishigh * jαhi + (1 - ishigh) * jαlo
    jβ = ishigh * jβhi + (1 - ishigh) * jβlo

    eqs = [
        INa ~ gNa * i_Nam^3 * i_Nah * i_Naj * (vm - E_Na),
        D(i_Nam) ~ inv(ms) * (mα - i_Nam * (mα + mβ)),
        D(i_Nah) ~ inv(ms) * (hα - i_Nah * (hα + hβ)),
        D(i_Naj) ~ inv(ms) * (jα - i_Naj * (jα + jβ)),
        E_Na ~ nernst(nao, nai),
        INab ~ gNab * (vm - E_Na),
    ]

    return ODESystem(eqs, t; name)
end
