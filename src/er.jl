"Sarcoplasmic reticulum system"
function get_ser_sys(Cai; fracPLB_CKp=0, fracPLBp=0, RyR_CKp=0, name=:sersys)
    @parameters begin
        VSR = 0.043 * 1.5 * 1.4pL
        VNSR = 0.9 * VSR
        VJSR = VSR - VNSR
        # RyR
        nRyR = 4
        nu1RyR = 0.01 // ms
        kaposRyR = 1 // ms
        kanegRyR = 0.16 // ms
        # SERCA
        VmaxfSR = 0.9996 * (μM // ms)
        VmaxrSR = VmaxfSR
        KmfSR = 0.5μM
        KmrSR = 7000 * KmfSR
        HfSR = 2
        HrSR = 1 * HfSR
        kSRleak = 5e-6 // ms
        fracPKA_PLBo = 1 - 0.079755
        ktrCaSR = 1 // 200ms
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
    end

    KmRyR = (1.35 * 2.6 * expit(-(Cai - 530μM) / 200μM) + 1.5 - 0.9 - 0.3 - 0.05) * μM
    fCKII_PLB = (1 - 0.5 * fracPLB_CKp)  # Max effect: fCKII_PLB=0.5
    fPKA_PLB = ((1 - fracPLB) / fracPKA_PLBo) * (1 - 0.5531) + 0.5531
    # Select smaller value (resulting in max reduction of Kmf)
    Kmfp = 2 * KmfSR * min(fCKII_PLB, fPKA_PLB)  #fCKII_PLB
    fSR = NaNMath.pow(Cai / Kmfp, HfSR)
    rSR = NaNMath.pow(CaNSR / KmrSR, HrSR)
    kleak = (1 / 2 + 5 * RyR_CKp / 2) * kSRleak

    return ODESystem([
        1 ~ PO1RyR + PC1RyR,
        Jrel ~ nu1RyR * PO1RyR * (CaJSR - Ca),
        D(PO1RyR) ~ kaposRyR * hil(Ca, KmRyR, nRyR) * PC1RyR - kanegRyR * PO1RyR,
        Jup ~ (VmaxfSR * fSR - VmaxrSR * rSR) / (1 + fSR + rSR),
        Jleak ~ kleak * (CaNSR - Cai),
        Jtr ~ ktrCaSR * (CaNSR - CaJSR),
        betaSR ~ inv(1 + csqntot * Kmcsqn / (CaJSR + Kmcsqn)^2),
        JCa_SR ~ Jleak - Jup + Jrel,
        D(CaJSR) ~ betaSR * (-Jrel + Jtr) / VJSR,
        D(CaNSR) ~ (Jup - Jleak - Jtr) / VNSR,
    ], t; name)
end
