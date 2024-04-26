# Cell geometry for calcium diffusion in neonatal mouse CMC model
using ModelingToolkit

# Calcium nuffer term
function beta_cai(Cai, TnI_PKAp)
    fracTnIpo = 0.062698  # Derived quantity (TnI_PKAp(baseline)/TnItot)
    fPKA_TnI = (1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIpo)) # Max effect +61#
    TrpnTotal = 35
    KmTrpn = 0.5 / fPKA_TnI
    CmdnTotal = 50
    KmCmdn = 2.38
    return inv(1 + TrpnTotal * KmTrpn / (Cai + KmTrpn)^2 + CmdnTotal * KmCmdn / (Cai + KmCmdn)^2)
end

function get_capde_eqs(;dx = 0.1, rSR_true = 6, rSL_true = 10.5)
    # Ca diffusion grid
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial indices
    m = length(j)
    @variables t Cai(t)[1:m] Cai_mean(t) Cai_sub_SR(t) Cai_sub_SL(t)
    @parameters Dca = 7 # mum^2/ms set to achive correct diff. speed 0.31 mum/ms

    eqs = [
        Cai_mean ~ sum(collect(Cai))/m,
        Cai_sub_SR ~ Cai[1],
        Cai_sub_SL ~ Cai[m],

    ]


    return eqs
end
