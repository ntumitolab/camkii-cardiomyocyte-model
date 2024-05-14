# I_kr: Rapidly Activating K Current

using ModelingToolkit

function get_ikr_eqs(vm, ek, Ko, Kcoeff=1)
    xrss = expit((vm+50)/7.5)
    tauxr = 1/(1.38e-3*(vm+7)/(1-exp(-0.123*(vm+7)))+6.1e-4*(vm+10)/(exp(0.145*(vm+10))-1))
    rkr = expit(-(vm+33)/22.4)

    @variables t iKrx(t) I_kr(t)
    @variables Gkr = 0.03*sqrt(Ko/5.4)

    D = Differential(t)

    return [
        D(iKrx) ~ (xrss- iKrx )/tauxr
        I_kr ~ Kcoeff * Gkr * iKrx * rkr * (vm - ek)
    ]
end
