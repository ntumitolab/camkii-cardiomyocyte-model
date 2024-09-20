# I_pca: Sarcolemmal Ca Pump Current


function get_ipca_eqs(Cajunc, Casl, Fjunc, Fsl,)
    @parameters IbarSLCaP = 0.0673 # [uA/uF]
    @parameters Q10SLCaP = 2.35
    @parameters KmPCa = 0.5e-3 # [mM]
    @variables t I_pca_junc(t) I_pca_sl(t) I_pca(t)

    imax = Q10SLCaP^Qpow * IbarSLCaP

    _ipca(ca) = hil(ca, KmPCa, 1.6)

    return [
        I_pca_junc ~ Fjunc * imax * _ipca(Cajunc)
        I_pca_sl ~ Fsl * imax * _ipca(Casl)
        I_pca ~ I_pca_junc + I_pca_sl
    ]
end
