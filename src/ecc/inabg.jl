# Background Na current
using ModelingToolkit

function get_inabg_eqs(vm, ena_junc, ena_sl, Fjunc, Fsl, loop=0)
    RyRp_WT_mean = 0.2101
    RyRp_OE_mean = 0.7387   # Derived (1 Hz, no loop)
    RyRp_OEloop_min = 0.7033 # Derived (1 Hz, OE loop)
    delta_loop = loop * ((3 / (RyRp_OE_mean - RyRp_WT_mean)) * RyR_CKp - (3 / (RyRp_OE_mean - RyRp_WT_mean)) * RyRp_WT_mean)

    @parameters GNaB = (4.5) * 0.297e-3 * (1 + delta_loop) # [mS/uF] changed from rabbit
    @variables t I_nabk_junc(t) I_nabk_sl(t) I_nabk(t)

    return [
        I_nabk_junc ~ Fjunc * GNaB * (vm - ena_junc)
        I_nabk_sl ~ Fsl * GNaB * (vm - ena_sl)
        I_nabk ~ I_nabk_junc + I_nabk_sl
    ]
end
