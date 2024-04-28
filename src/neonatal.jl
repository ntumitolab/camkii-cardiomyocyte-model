# Neonatal mouse CMC model
using ModelingToolkit

# Calcium buffer term
function beta_cai(Cai, TnI_PKAp=0)
    fracTnIpo = 0.062698  # Derived quantity (TnI_PKAp(baseline)/TnItot)
    fPKA_TnI = (1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIpo)) # Max effect +61#
    TrpnTotal = 35
    KmTrpn = 0.5 / fPKA_TnI
    CmdnTotal = 50
    KmCmdn = 2.38
    return inv(1 + TrpnTotal * KmTrpn / (Cai + KmTrpn)^2 + CmdnTotal * KmCmdn / (Cai + KmCmdn)^2)
end

function build_neonatal_model(; dx = 0.1, rSR_true = 6, rSL_true = 10.5)

    @constants begin
        Vmyo = 2.1454e-11Liter           # [L]
        Vdyad = 1.7790e-014Liter         # [L]
        VSL = 6.6013e-013Liter           # [L]
        kSLmyo = 8.587e-15Liter/ms       # [L/msec]
        CaMKIItotDyad = 120μM            # [uM]
        BtotDyad = 1.54 / 8.293e-4μM     # [uM]
    end

    @parameters begin
        Cao = 1796μM   # [uM]
        Nao = 154578μM # [uM]
        Ko = 5366μM    # [uM]
        Mg = 1mM       # [mM]
        K = 135mM      # [mM]
    end

    # Ca diffusion
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial indices
    m = length(j)
    V_sub_SR = 4 / 3 * pi * (rSR_true + dx)^3 / 1000 - 4 / 3 * pi * (rSR_true)^3 / 1000 # pl
    V_sub_SL = 4 / 3 * pi * rSL_true^3 / 1000 - 4 / 3 * pi * (rSL_true - dx)^3 / 1000 #pl
    @variables t Cai(t)[1:m] Cai_mean(t) Cai_sub_SR(t) Cai_sub_SL(t) JCa_SR(t) JCa_SL(t)
    @parameters Dca = 7 # mum^2/ms set to achive correct diff. speed 0.31 mum/ms
    D = Differential(t)
    eqs = [
        Cai_mean ~ sum(collect(Cai))/m,
        Cai_sub_SR ~ Cai[1],
        Cai_sub_SL ~ Cai[m],
        D(Cai[1]) ~ (Dca / (j[1] * dx^2) * ((1 + j[1]) * Cai[2] - 2 * j[1] * Cai[1] + (j[1] - 1) * Cai[1]) + JCa_SR / V_sub_SR) * beta_cai(Cai[1], TnI_PKAp),
        D(Cai[m]) ~ (Dca / (j[m] * dx^2) * ((1 + j[m]) * Cai[m] - 2 * j[m] * Cai[m] + (j[m] - 1) * Cai[m-1]) + JCa_SL / V_sub_SL) * beta_cai(Cai[m], TnI_PKAp),
    ]

    for i in 2:m-1
        eq = D(Cai[i]) ~ (Dca / (j[i] * dx^2) * ((1 + j[i]) * Cai[i+1] - 2 * j[i] * Cai[i] + (j[i] - 1) * Cai[i-1])) * beta_cai(Cai[i], TnI_PKAp)
        push!(eqs, eq)
    end


    return eqs
end
