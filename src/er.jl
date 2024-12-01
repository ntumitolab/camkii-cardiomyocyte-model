"Sarcoplasmic reticulum system"
function get_ser_sys(Cai_sub_SR; fracPLB_CKp=0, fracPLBp=0, RyR_CKp=0.2, V_sub_SR=0.046pL, name=:sersys)
    @parameters begin
        VSR = 0.0903pL
        VNSR = 0.9 * VSR
        VJSR = VSR - VNSR
        # RyR
        kRyR = 10Hz
        kaposRyR = 1000Hz
        kanegRyR = 160Hz
        RyRsensitivity = 1.0
        # SERCA
        VmaxSR = 0.9996mM * Hz
        KmfSR = 0.5μM
        KmrSR = 3.5mM
        kSRleak = 0.005Hz
        fracPKA_PLBo = 1 - 0.079755
        ktrCaSR = 5Hz
        csqntot = 24.750mM
        Kmcsqn = 0.8mM
    end

    @variables begin
        CaNSR(t) = 830μM
        CaJSR(t) = 830μM
        PO1RyR(t) = 0.0026
        PC1RyR(t)
        Jrel(t)
        Jup(t)
        Jleak(t)
        Jtr(t)
        betaSR(t)
        JCa_SR(t)
        KmRyR(t)
    end

    fCKII_PLB = (1 - 0.5 * fracPLB_CKp)  # Max effect: fCKII_PLB=0.5
    fPKA_PLB = ((1 - fracPLBp) / fracPKA_PLBo) * (1 - 0.5531) + 0.5531
    # Select the smaller value (resulting in max reduction of Kmf)
    Kmfp = KmfSR * min(fCKII_PLB, fPKA_PLB)  # fCKII_PLB
    fSR = (Cai_sub_SR / Kmfp)^2
    rSR = (CaNSR / KmrSR)^2
    kleak = (1 + 5 * RyR_CKp) * kSRleak / 2

    eqs = [
        1 ~ PO1RyR + PC1RyR,
        Jrel ~ kRyR * PO1RyR * (CaJSR - Cai_sub_SR),
        KmRyR ~ (3.51 / (1 + exp((CaJSR - 530μM) / 200μM)) + 0.25) * 1μM,
        D(PO1RyR) ~ kaposRyR * hil(Cai_sub_SR, KmRyR / RyRsensitivity, 4) * PC1RyR - kanegRyR * PO1RyR,
        Jup ~ VmaxSR * (fSR - rSR) / (1 + fSR + rSR),
        Jleak ~ kleak * (CaNSR - Cai_sub_SR),
        Jtr ~ ktrCaSR * (CaNSR - CaJSR),
        betaSR ~ inv(1 + csqntot * Kmcsqn / (CaJSR + Kmcsqn)^2),
        JCa_SR ~ (Jleak - Jup + Jrel) * (1pL / V_sub_SR),
        D(CaJSR) ~ betaSR * (-Jrel + Jtr) * (1pL / VJSR),
        D(CaNSR) ~ (Jup - Jleak - Jtr) * (1pL / VNSR),
    ]
    return ODESystem(eqs, t; name)
end
