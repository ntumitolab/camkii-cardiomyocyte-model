# I_kur - IK,slow

using ModelingToolkit

function get_ikur_eqs(vm, ek, Kcoeff=1)
    # I_kur
    xurss = expit((vm+15)/14)
    yurss = expit(-(vm+48)/6.2)
    tauxur = 0.95+0.05*exp(-0.08*vm)
    tauyur1 = 400+900*exp(-((vm+55)/16)^2)-250/(1+exp(-(vm+60)/8))
    tauyur2 = 400+900*exp(-((vm+55)/16)^2)+550/(1+exp(-(vm+60)/8))

    # I_ss
    xssss = xurss
    tauxss = 70*exp(-((vm+43)/30)^2)+14

    # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    fracIKurp0 = 0.437635  # Derived quantity (IKur_PKAp(baseline)/IKurtot)
    fracIKurpISO = 0.718207 # Derived quantity (IKur_PKAp(ISO)/IKurtot)
    a_Kur = (1.20-1)/(fracIKurpISO/fracIKurp0-1);
    fracIKuravail = (1-a_Kur)+a_Kur*(IKur_PKAp/fracIKurp0) # +20% with 0.1 uM ISO

    @variables t
    @variables iKslowx(t) # y(84)
    @variables iKslowy1(t) # y(85)
    @variables iKslowy2(t) # y(87)
    @variables iKss(t)
    @variables I_kur1(t) I_kur2(t) I_kur(t) I_ss(t)
    @parameters Gkur1 = 1.1*0.16 # fast
    @parameters Gkur2 = 0.14 # slow
    @parameters Gss = 0.15 # [mS/uF] only in MOUSE

    D = Differential(t)

    return [
        D(iKslowx) ~ (xurss-iKslowx)/tauxur
        D(iKslowy1) ~ (yurss-iKslowy1)/tauyur1
        D(iKslowy2) ~ (yurss-iKslowy2)/tauyur2
        I_kur1 ~ Kcoeff*fracIKuravail*Gkur1*iKslowx* iKslowy1*(vm-ek)
        I_kur2 ~ Kcoeff*Gkur2*iKslowx * iKslowy2*(vm-ek)
        I_kur ~ I_kur1 + I_kur2
        D(iKss) ~ (xssss-iKss)/tauxss
        I_ss ~ Kcoeff*Gss*iKss*(vm-ek)
    ]

end
