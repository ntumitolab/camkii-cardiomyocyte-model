"Potassium currents"
function get_ik_sys(k_i, k_o, na_i, na_o, vm; IKUR_PKAp=0, name=:iksys)
    @parameters begin
        # IK1: time-independent
        GK1 = 0.0515mSμF
        # Ito: Where does this come from? Perhaps here: https://modeldb.science/262081
        Gt = 0.1mSμF
        f_is = 0.706
        # IKs (and IKur)
        GKs = 0.05mSμF
        nKstau = 750ms
        # IKr
        GKr = 0.06mSμF
        kf = 0.023761 / ms
        kb = 0.036778 / ms
        # Hyperpolarization activated current (Funny current, If)
        Gf = 0.021mSμF
        fNa = 0.2
    end

    @variables begin
        E_Na(t)
        E_K(t)
        IK1(t)
        Ito(t)
        i_r(t) = 0.00702
        i_s(t) = 0.9660
        i_sslow(t) = 0.22156
        sinf(t)
        rinf(t)
        slowinf(t)
        taur(t)
        taus(t)
        tausslow(t)
        IKs(t)
        i_nKs(t) = 0.09243
        nksinf(t)
        IKr(t)
        E_Kr(t)
        i_CK0(t)
        i_CK1(t) = 0.00188
        i_CK2(t) = 0.00977
        i_OK(t) = 0.26081
        i_IK(t) = 0.07831
        IfNa(t)
        IfK(t)
        If(t)
        yinf(t)
        tauy(t)
        i_y(t) = 0.07192
    end

    V = vm  # mV is the base unit now
    vk1 = V - E_K - 6.1373

    # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    fracIKurp0 = 0.437635       # Derived quantity (IKur_PKAp(baseline)/IKurtot)
    fracIKurpISO = 0.718207     # Derived quantity (IKur_PKAp(ISO)/IKurtot)
    a_Kur = (1.20 - 1) / (fracIKurpISO / fracIKurp0 - 1)
    fracIKuravail = (1 - a_Kur) + a_Kur * (IKUR_PKAp / fracIKurp0)  # +20# with 0.1 uM ISO
    alphan = 0.00000481333 / 0.128 * exprel(-0.128 * (V + 26.5))
    betan = 0.0000953333 * exp(-0.038 * (V + 26.5))

    # IKr
    alphaa0 = 0.022348/ms * exp(0.01176 * V)
    betaa0 = 0.047002/ms * exp(-0.0631 * V)
    alphaa1 = 0.013733/ms * exp(0.038198 * V)
    betaa1 = 0.0000689/ms * exp(-0.04178 * V)
    alphai_mERG = 0.090821/ms * exp(0.023391 * V)
    betai_mERG = 0.006497/ms * exp(-0.03268 * V)

    rc0c1 = alphaa0 * i_CK0 - betaa0 * i_CK1
    rc1c2 = kf * i_CK1 - kb * i_CK2
    rc2o = alphaa1 * i_CK2 - betaa1 * i_OK
    roi = alphai_mERG * i_OK - betai_mERG * i_IK

    fK = 1 - fNa

    return System([
            E_Na ~ nernst(na_o, na_i),
            E_K ~ nernst(k_o, k_i),
            IK1 ~ GK1 * hil(k_o, 210μM) * expit(-0.0319 * vk1, vk1, 0.1653),
            sinf ~ expit((V + 31.97156) / -4.64291),
            rinf ~ expit((V - 3.55716) / 14.61299),
            slowinf ~ sinf,
            taur ~ 1000ms / (45.16 * exp(0.03577 * (V + 50)) + 98.9 * exp(-0.1 * (V + 38))),
            taus ~ (350 * exp(-(((V + 70) / 15)^2)) + 35 - 26.9) * ms,
            tausslow ~ (3700 * exp(-(((V + 70) / 30)^2)) + 35 + 37.4) * ms,
            Ito ~ Gt * i_r * (f_is * i_s + (1 - f_is) * i_sslow) * (vm - E_K),
            D(i_r) ~ (rinf - i_r) / taur,
            D(i_s) ~ (sinf - i_s) / taus,
            D(i_sslow) ~ (slowinf - i_sslow) / tausslow,
            IKs ~ GKs * i_nKs^2 * (vm - E_K) * fracIKuravail * 2,
            nksinf ~ alphan / (alphan + betan),
            D(i_nKs) ~ (nksinf - i_nKs) / nKstau,
            E_Kr ~ nernst(0.98 * k_o + 0.02 * na_o, 0.98 * k_i + 0.02 * na_i),
            IKr ~ GKr * i_OK * (vm - E_Kr),
            1 ~ i_CK0 + i_CK1 + i_CK2 + i_OK + i_IK,
            D(i_CK1) ~ rc0c1 - rc1c2,
            D(i_CK2) ~ rc1c2 - rc2o,
            D(i_OK) ~ rc2o - roi,
            D(i_IK) ~ roi,
            yinf ~ expit(-(V + 78.65) / 6.33),
            tauy ~ 1000ms / (0.11885 * exp((V + 75) / 28.37) + 0.56236 * exp(-(V + 75) / 14.19)),
            IfNa ~ Gf * fNa * i_y * (vm - E_Na),
            IfK ~ Gf * fK * i_y * (vm - E_K),
            If ~ IfNa + IfK,
            D(i_y) ~ (yinf - i_y) / tauy
        ], t; name)
end
