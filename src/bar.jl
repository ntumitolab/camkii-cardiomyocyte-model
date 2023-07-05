#=
    This module describes the beta-adrenergic signaling pathway in mouse
    ventricular myocyte, and this file was built upon the code developeded
    by Yang and Saucerman.
    Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
    feed-forward regulator of Ca2+ in beta-adrenergic signaling,
    accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular
    Cardiology 52, 1048-1055.
=#

using Catalyst

function get_bar_equations()
    @variables t
    D = Differential(t)

    # Drug concentrations
    @parameters (ISO=0, FSK=0, IBMX=0)

    @parameters b1ARtot = 0.00528  # (uM) total b1-AR protein (MOUSE); RABBIT=0.028
    @parameters Gstot = 3.83       # (uM) total Gs protein
    @variables b1AR(t) b1ARact(t) b1AR_S464(t) b1AR_S301(t) LR(t) LRG(t) RG(t) Gs(t) Gsby(t)

    eqs = [
        b1ARtot ~ b1ARact + b1AR_S464 + b1AR_S301,
        b1ARact ~ b1AR + LR + LRG + RG,
        Gstot ~ Gs + LRG + RG + Gsby,
        D(LR) ~ kf_LR*ISO*b1AR - kr_LR*LR + kr_LRG*LRG - kf_LRG*LR*Gs,
        D(LRG) ~ kf_LRG*LR*Gs - kr_LRG*LRG - k_G_act*LRG,
        D(RG) ~ kf_RG*b1AR*Gs - kr_RG*RG - k_G_act*RG
    ]

    return eqs
end
