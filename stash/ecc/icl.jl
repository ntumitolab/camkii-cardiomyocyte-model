# I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
using ModelingToolkit

function get_icl_eqs(vm, ecl, Cajunc, Casl, Fjunc, Fsl)

    @parameters GClCa = 0.109625    # [mS/uF]
    @parameters GClB = 9e-3         # [mS/uF]
    @parameters KdClCa = 100e-3     # [mM]

    @variables I_ClCa_junc(t) I_ClCa_sl(t) I_ClCa(t) I_Clbk(t)

    return [
        I_ClCa_junc ~ Fjunc*GClCa * hil(Cajunc, KdClCa) * (vm-ecl)
        I_ClCa_sl ~ Fsl*GClCa * hil(Casl, KdClCa) * (vm-ecl)
        I_ClCa ~ I_ClCa_junc + I_ClCa_sl
        I_Clbk ~ GClB * (vm-ecl)
    ]
end
