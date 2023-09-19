# I_pca: Sarcolemmal Ca Pump Current
using ModelingToolkit

_ipca(ca, KmPCa) = hil(ca, KmPCa, 1.6)

function get_ipca_eqs(Cajunc, Casl, Fjunc, Fsl,)
    @parameters IbarSLCaP = 0.0673 # [uA/uF]
    @parameters Q10SLCaP = 2.35
    @parameters KmPCa = 0.5e-3 # [mM]
    @variables t I_pca_junc(t) I_pca_sl(t) I_pca(t)

    imax = Q10SLCaP^Qpow * IbarSLCaP

    return [
        I_pca_junc ~ Fjunc * imax * _ipca(Cajunc, KmPCa)
        I_pca_sl ~ Fsl * imax * _ipca(Casl, KmPCa)
        I_pca ~ I_pca_junc + I_pca_sl
    ]
end
