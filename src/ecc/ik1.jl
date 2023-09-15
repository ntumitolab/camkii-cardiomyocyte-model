# I_k1: Time-Independent K Current (I_ki)
using ModelingToolkit

function get_iki_eqs(vm, ek, Ko, Kcoeff=1, CKIIflag=0)
    aki = 1.02/(1+exp(0.2385*(vm-ek-59.215)))
    bki =(0.49124*exp(0.08032*(vm+5.476-ek))+exp(0.06175*(vm-ek-594.31)))/(1 + exp(-0.5143*(vm-ek+4.753)))
    kiss = aki/(aki+bki)

    @parameters GK1 = 0.3*sqrt(Ko/5.4) * (1 - 0.5 * CKIIflag)

    @variables t I_ki(t)

    return [
        I_ki ~ Kcoeff * GK1 * kiss * (vm - ek)
    ]
end
