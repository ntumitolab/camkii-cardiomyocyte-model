# I_cabk: Ca Background Current
using ModelingToolkit

function get_icabk_eqs(Fjunc, Fsl, vm, eca_junc, eca_sl, CaffeineFlag=0)
    @parameters GCaB = 3 * 2.513e-4 * (1-CaffeineFlag) # [uA/uF] changed from rabbit (2.513e-4)
    @variables t I_cabk_junc(t) I_cabk_sl(t) I_cabk(t)

    return [
        I_cabk_junc ~ Fjunc * GCaB * (vm - eca_junc)
        I_cabk_sl ~ Fsl * GCaB * (vm - eca_sl)
        I_cabk ~ I_cabk_junc + I_cabk_sl
    ]
end
