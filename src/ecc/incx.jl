# I_ncx: Na/Ca Exchanger flux
using ModelingToolkit

function _incx(vm, Nai, Nao, Cai, Cao, KmNai, KmNao, KmCai, KmCao, Kdact, nu, ksat)
    ka = hill(Cai, Kdact, 3)
    s1 = exp(nu * vm * FoRT) * Nai^3 * Cao
    s2 = exp((nu-1) * vm * FoRT) * Nao^3 * Cai
    s3 =
end

function get_incx_eqs(Cajunc, Casl, Cao, Najunc, Nasl, Nao)

end
