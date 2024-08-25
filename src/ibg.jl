# Background currents
"Background current"
function get_ibg_sys(vm, E_Na, E_Ca; name=:ibgsys)
    @parameters gCab = 0.0008mS / cm^2
    @parameters gNab = 0.0026mS / cm^2
    @variables ICab(t) INab(t)
    eqs = [
        ICab ~ gCab * (vm - E_Ca),
        INab ~ gNab * (vm - E_Na),
    ]
    return ODESystem(eqs, t; name)
end

"Funny current (If)"
function get_if_sys(vm, E_Na, E_K; name=:ifsys)
    @parameters gf = 0.021mS / cm^2 fNa = 0.2
    @variables IfNa(t) IfK(t) If(t) yinf(t) tauy(t) i_y(t) = 0.07192
    fK = 1 - fNa
    V = vm * Volt / mV # Convert voltage to mV
    eqs = [
        yinf ~ expit(-(V + 78.65) / 6.33),
        tauy ~ 1 / (0.11885 * exp((V + 75) / 28.37) + 0.56236 * exp(-(V + 75) / 14.19)),
        IfNa ~ gf * fNa * i_y * (vm - E_Na),
        IfK ~ gf * fK * i_y * (vm - E_K),
        If ~ IfNa + IfK,
        D(i_y) * tauy ~ yinf - i_y
    ]
    return ODESystem(eqs, t; name)
end
