# I_kp: Plateau K current

using ModelingToolkit

function get_ikp_eqs(vm, ek, Fjunc, Fsl, Kcoeff=1)
    kp_kp = expit((vm - 7.488)/5.98)
    @parameters Gkp = 0.001
    @variables t I_kp_junc(t) I_kp_sl(t) I_kp(t)
    return [
        I_kp_junc ~ Kcoeff*Fjunc*Gkp*kp_kp*(vm-ek)
        I_kp_sl ~ Kcoeff*Fsl*Gkp*kp_kp*(vm-ek)
        I_kp ~ I_kp_junc+I_kp_sl
    ]
end
