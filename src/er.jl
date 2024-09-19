"Sarcoplasmic reticulum system"
function get_ser_sys(Cai_sub_SR; fracPLB_CKp=0, fracPLBp=0, RyR_CKp=0.2, V_sub_SR=0.046pL, name=:sersys)
    @parameters begin
        VSR = 0.093pL
        VNSR = 0.9 * VSR
        VJSR = VSR - VNSR
        # RyR
        kRyR = 10Hz
        kaposRyR = 1000Hz
        kanegRyR = 160Hz
        # SERCA
        VmaxSR = 0.9996 * μM / ms
        KmfSR = 0.5μM
        KmrSR = 7000 * KmfSR
        kSRleak = 5e-6 / ms
        fracPKA_PLBo = 1 - 0.079755
        ktrCaSR = inv(200ms)
        csqntot = 24750μM
        Kmcsqn = 800μM
    end

    @variables begin
        CaNSR(t) = 619.09843μM
        CaJSR(t) = 613.87556μM
        PO1RyR(t) = 0.0037
        PC1RyR(t)
        Jrel(t)
        Jup(t)
        Jleak(t)
        Jtr(t)
        betaSR(t)
        JCa_SR(t)
        KmRyR(t)
        dPO1RyR(t)
    end

    fCKII_PLB = (1 - 0.5 * fracPLB_CKp)  # Max effect: fCKII_PLB=0.5
    fPKA_PLB = ((1 - fracPLBp) / fracPKA_PLBo) * (1 - 0.5531) + 0.5531
    # Select the smaller value (resulting in max reduction of Kmf)
    Kmfp = 2 * KmfSR * min(fCKII_PLB, fPKA_PLB)  #fCKII_PLB
    fSR = (Cai_sub_SR / Kmfp)^2
    rSR = (CaNSR / KmrSR)^2
    kleak = (1 + 5 * RyR_CKp) * kSRleak / 2

    eqs = [
        1 ~ PO1RyR + PC1RyR,
        Jrel ~ kRyR * PO1RyR * (CaJSR - Cai_sub_SR),
        KmRyR ~ (3.51 / (1 + exp((CaJSR - 530μM) / 200μM)) + 0.25) * μM,
        dPO1RyR ~ kaposRyR * hil(Cai_sub_SR, KmRyR, 4) * PC1RyR - kanegRyR * PO1RyR,
        D(PO1RyR) ~ dPO1RyR,
        Jup ~ VmaxSR * (fSR - rSR) / (1 + fSR + rSR),
        Jleak ~ kleak * (CaNSR - Cai_sub_SR),
        Jtr ~ ktrCaSR * (CaNSR - CaJSR),
        betaSR ~ inv(1 + csqntot * Kmcsqn / (CaJSR + Kmcsqn)^2),
        JCa_SR ~ Jleak - Jup + Jrel,
        D(CaJSR) ~ betaSR * (-Jrel * V_sub_SR + Jtr * VNSR)/VJSR,
        D(CaNSR) ~ (Jup - Jleak) * V_sub_SR / VNSR - Jtr,
    ]
    return ODESystem(eqs, t; name)
end
