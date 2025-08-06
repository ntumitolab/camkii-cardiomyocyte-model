"Calcium diffusion between sarcolemma (SL) and sarcoplasmic reticulum (SR)"
function get_ca_pde_eqs(;
    Cai_default=121.13nM,
    dx=0.1μm,
    rSR_true=6μm,
    rSL_true=10.5μm,
    TnI_PKAp=0,
    JCa_SR=0,
    JCa_SL=0,
)
    rSR = rSR_true + 0.5 * dx
    rSL = rSL_true - 0.5 * dx
    j = round(rSR / dx):1:round(rSL / dx) # Spatial indices
    m = length(j)
    @variables Cai(t)[1:m] Cai_mean(t) Cai_sub_SR(t) Cai_sub_SL(t)
    @parameters begin
        Dca = 7μm^2 / ms
        TrpnTotal = 35μM  # Half of adult rat CMC (70 μM)
        CmdnTotal = 30μM  # Adjusted from 50μM
        KmCmdn = 2.38μM
        KmTrpn = 0.5μM
        fracTnIp0 = 0.062698 # Baseline effect
    end
    # Effect 0.96 ~ 1.61
    fPKA_TnI = 1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIp0)
    KmTrpnNew = KmTrpn / fPKA_TnI

    "Calcium buffering power by troponin and calmodulin"
    beta_cai(Ca) = inv(1 + TrpnTotal * KmTrpnNew / (Ca + KmTrpnNew)^2 + CmdnTotal * KmCmdn / (Ca + KmCmdn)^2)

    eqs_cai = [
        Cai_mean ~ sum(collect(Cai)) / m,
        Cai_sub_SR ~ Cai[1],
        Cai_sub_SL ~ Cai[m],
        D(Cai[1]) ~ (Dca / (j[1] * dx^2) * ((1 + j[1]) * Cai[2] - 2 * j[1] * Cai[1] + (j[1] - 1) * Cai[1]) + JCa_SR) * beta_cai(Cai[1]),
        D(Cai[m]) ~ (Dca / (j[m] * dx^2) * ((1 + j[m]) * Cai[m] - 2 * j[m] * Cai[m] + (j[m] - 1) * Cai[m-1]) + JCa_SL) * beta_cai(Cai[m]),
    ]

    defaults_cai = [Cai[1] => Cai_default, Cai[m] => Cai_default]

    for i in 2:m-1
        push!(eqs_cai, D(Cai[i]) ~ (Dca / (j[i] * dx^2) * ((1 + j[i]) * Cai[i+1] - 2 * j[i] * Cai[i] + (j[i] - 1) * Cai[i-1])) * beta_cai(Cai[i]))
        push!(defaults_cai, Cai[i] => Cai_default)
    end
    return (; eqs_cai, defaults_cai)
end

function get_ca_pde_sys(;
    Cai_default=121.13nM,
    dx=0.1μm,
    rSR_true=6μm,
    rSL_true=10.5μm,
    TnI_PKAp=0,
    JCa_SR=0,
    JCa_SL=0,
    name=:capdesys
)
    @unpack eqs_cai, defaults_cai = get_ca_pde_eqs(; Cai_default, dx, rSR_true, rSL_true, TnI_PKAp, JCa_SR, JCa_SL)
    return System(eqs_cai, t; name, defaults_cai)
end
