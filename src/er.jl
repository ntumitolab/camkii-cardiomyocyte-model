"Sarcoplasmic reticulum system"
function get_ser_sys(Cai_sub_SR; fracPLB_CKp=0, fracPLBp=0, RyR_CKp=0, name=:sersys)
    @parameters begin
        VSR = 0.043 * 1.5 * 1.4pL
        VNSR = 0.9 * VSR
        VJSR = 0.1 * VSR
        # RyR
        nu1RyR = 10Hz
        kaposRyR = 1000Hz
        kanegRyR = 160Hz
        # SERCA
        VmaxfSR = 0.9996 * μM / ms
        VmaxrSR = VmaxfSR
        KmfSR = 0.5μM
        KmrSR = 7000 * KmfSR
        HfSR = 2
        HrSR = 1 * HfSR
        kSRleak = 5e-6 / ms
        fracPKA_PLBo = 1 - 0.079755
        ktrCaSR = 5Hz
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
    # Select smaller value (resulting in max reduction of Kmf)
    Kmfp = 2 * KmfSR * min(fCKII_PLB, fPKA_PLB)  #fCKII_PLB
    fSR = NaNMath.pow(Cai_sub_SR / Kmfp, HfSR)
    rSR = NaNMath.pow(CaNSR / KmrSR, HrSR)
    kleak = (1 + 5 * RyR_CKp) * kSRleak / 2

    return ODESystem([
        1 ~ PO1RyR + PC1RyR,
        Jrel ~ nu1RyR * PO1RyR * (CaJSR - Cai_sub_SR),
        KmRyR ~ (3.51 * inv(1 + exp((CaJSR - 530μM) / 200μM)) + 0.25) * μM,
        dPO1RyR ~ kaposRyR * hil(Cai_sub_SR, KmRyR, 4) * PC1RyR - kanegRyR * PO1RyR,
        D(PO1RyR) ~ dPO1RyR,
        Jup ~ (VmaxfSR * fSR - VmaxrSR * rSR) / (1 + fSR + rSR),
        Jleak ~ kleak * (CaNSR - Cai_sub_SR),
        Jtr ~ ktrCaSR * (CaNSR - CaJSR),
        betaSR ~ inv(1 + csqntot * Kmcsqn / (CaJSR + Kmcsqn)^2),
        JCa_SR ~ Jleak - Jup + Jrel,
        D(CaJSR) ~ betaSR * (-Jrel + Jtr) / VJSR,
        D(CaNSR) ~ (Jup - Jleak - Jtr) / VNSR,
    ], t; name)
end
