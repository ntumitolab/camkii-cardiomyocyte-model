# Ca diffusion model
using ModelingToolkit

function get_ca_pde_eqs(;
    dx=0.1μm,
    rSR_true=6μm,
    rSL_true=10.5μm,
    V_sub_SR = 4 / 3 * pi * ((rSR_true + dx)^3 - (rSR_true)^3),
    V_sub_SL = 4 / 3 * pi * (rSL_true^3 - (rSL_true - dx)^3),
    TnI_PKAp = 0,
)
    @variables t
    D = Differential(t)
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial indices
    m = length(j)
    @variables Cai(t)[1:m] Cai_mean(t) Cai_sub_SR(t) Cai_sub_SL(t) JCa_SR(t) JCa_SL(t)
    @parameters Dca = 7μm^2/ms
    eqs = [
        Cai_mean ~ sum(collect(Cai)) / m,
        Cai_sub_SR ~ Cai[1],
        Cai_sub_SL ~ Cai[m],
        D(Cai[1]) ~ (Dca / (j[1] * dx^2) * ((1 + j[1]) * Cai[2] - 2 * j[1] * Cai[1] + (j[1] - 1) * Cai[1]) + JCa_SR / V_sub_SR) * beta_cai(Cai[1]; TnI_PKAp),
        D(Cai[m]) ~ (Dca / (j[m] * dx^2) * ((1 + j[m]) * Cai[m] - 2 * j[m] * Cai[m] + (j[m] - 1) * Cai[m-1]) + JCa_SL / V_sub_SL) * beta_cai(Cai[m]; TnI_PKAp),
    ]

    for i in 2:m-1
        eq = D(Cai[i]) ~ (Dca / (j[i] * dx^2) * ((1 + j[i]) * Cai[i+1] - 2 * j[i] * Cai[i] + (j[i] - 1) * Cai[i-1])) * beta_cai(Cai[i]; TnI_PKAp)
        push!(eqs, eq)
    end

    return eqs
end
