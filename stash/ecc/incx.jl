# I_ncx: Na/Ca Exchanger flux
using ModelingToolkit

function _incx(vm, Nai, Nao, Cai, Cao, KmNai, KmNao, KmCai, KmCao, Kdact, nu, ksat)
    ka = hill(Cai, Kdact, 3)
    s1 = exp(nu * vm * FoRT) * Nai^3 * Cao
    s2 = exp((nu-1) * vm * FoRT) * Nao^3 * Cai
    s3 = KmCai * Nao^3 * (1 + (Nai/KmNai)^3) + KmNao^3 * Cai + KmNai^3 * Cao * (1 + Cai/KmCai) + KmCao * Nai^3 + Nai^3 * Cao + Nao^3 * Cai
    s4 = s3 * (1 + ksat * exp((nu-1) * vm * FoRT))
    return (s1 - s2) / s4
end

function get_incx_eqs(Cajunc, Casl, Cao, Najunc, Nasl, Nao, Fjunc_ncx, Fsl_ncx, CKIIflag=0)

    @parameters IbarNCX = 1 + 0.5 * CKIIflag # [uA/uF] changed from rabbit (9)
    @parameters KmCai = 3.59e-3    #[mM]
    @parameters KmCao = 1.3        #[mM]
    @parameters KmNai = 12.29      #[mM]
    @parameters KmNao = 87.5       #[mM]
    @parameters ksat = 0.27        #[none]
    @parameters nu = 0.35          #[none]
    @parameters Kdact = 1/2*0.256e-3; #[mM] changed from rabbit
    @parameters Q10NCX = 1.57      #[none]

    @variables t I_ncx_junc(t) I_ncx_sl(t) I_ncx(t)

    return [
        I_ncx_junc ~ Fjunc_ncx * IbarNCX * Q10NCX^Qpow * _incx(vm, Najunc, Nao, Cajunc, Cao, KmNai, KmNao, KmCai, KmCao, Kdact, nu, ksat)
        I_ncx_sl ~ Fsl_ncx * IbarNCX * Q10NCX^Qpow * _incx(vm, Nasl, Nao, Casl, Cao, KmNai, KmNao, KmCai, KmCao, Kdact, nu, ksat)
        I_ncx ~ I_ncx_junc + I_ncx_sl
    ]
end
