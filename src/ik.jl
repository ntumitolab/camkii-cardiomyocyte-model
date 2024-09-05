# Potassium currents
"Potassium currents"
function get_ik_eqs(vm, E_K, K_i, K_o, Na_i, Na_o, IKur_PKAp=0; name=:iksys)
    V = vm * Volt / mV # Convert voltage to mV

    # IK1: time-independent
    @parameters gK1 = 0.0515mS / cm^2 * hil(K_o, 210μM)
    @variables IK1(t)
    vk1 = vm - E_K - 6.1373mV

    # Ito: Where does this come from? Perhaps here: https://modeldb.science/262081
    @parameters gt = 0.1mS / cm^2 f_is = 0.706
    @variables begin
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
    end

    # IKs (and IKur)
    @parameters GKs = 0.05mS / cm^2 nKstau = 750ms
    @variables begin
        IKs(t)
        i_nKs(t) = 0.09243
    end

    # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    fracIKurp0 = 0.437635       # Derived quantity (IKur_PKAp(baseline)/IKurtot)
    fracIKurpISO = 0.718207     # Derived quantity (IKur_PKAp(ISO)/IKurtot)
    a_Kur = (1.20 - 1) / (fracIKurpISO / fracIKurp0 - 1)
    fracIKuravail = (1 - a_Kur) + a_Kur * (IKur_PKAp / fracIKurp0)  # +20# with 0.1 uM ISO
    alphan = 0.00000481333 * (V + 26.5) / (-expm1(-0.128 * (V + 26.5)))
    betan = 0.0000953333 * exp(-0.038 * (V + 26.5))

    # IKr
    @parameters begin
        GKr = 0.06mS / cm^2
        kf = 0.023761 / ms
        kb = 0.036778 / ms
    end

    @variables begin
        IKr(t)
        E_Kr(t)
        CK0(t)
        i_CK1(t) = 0.00188
        i_CK2(t) = 0.00977
        i_OK(t) = 0.26081
        i_IK(t) = 0.07831
    end

    alphaa0 = 0.022348 / ms * exp(0.01176 * V)
    betaa0 = 0.047002 / ms * exp(-0.0631 * V)
    alphaa1 = 0.013733 / ms * exp(0.038198 * V)
    betaa1 = 0.0000689 / ms * exp(-0.04178 * V)
    alphai_mERG = 0.090821 / ms * exp(0.023391 * V)
    betai_mERG = 0.006497 / ms * exp(-0.03268 * V)

    rc0c1 = alphaa0 * CK0 - betaa0 * i_CK1
    rc1c2 = kf * i_CK1 - kb * i_CK2
    rc2o = alphaa1 * i_CK2 - betaa1 * i_OK
    roi = alphai_mERG * i_OK - betai_mERG * i_IK

    eqs = [
        IK1 ~ gK1 * vk1 / (0.1653 + exp(0.0319 / mV * vk1)),
        sinf ~ expit((vm + 31.97156mV) / -4.64291mV),
        rinf ~ expit((vm - 3.55716mV) / 14.61299mV),
        slowinf ~ sinf,
        taur ~ inv(45.16 * exp(0.03577 / mV * (vm + 50mV)) + 98.9 * exp(-0.1 / mV * (vm + 38mV))),
        taus ~ (0.35 * exp(-(((vm + 70mV) / 15mV)^2)) + 0.035) - 26.9ms,
        tausslow ~ (3.7 * exp(-(((vm + 70mV) / 30mV)^2)) + 0.035) + 37.4ms,
        Ito ~ gt * i_r * (f_is * i_s + (1 - f_is) * i_sslow) * (vm - E_K),
        D(i_r) ~ (rinf - i_r)/taur,
        D(i_s) ~ (sinf - i_s)/taus,
        D(i_sslow) ~ (slowinf - i_sslow)/tausslow,
        IKs ~ GKs * i_nKs^2 * (vm - E_K) * fracIKuravail * 2,
        D(i_nKs) ~ (alphan / (alphan + betan) - i_nKs)/nKstau,
        E_Kr ~ nernst(0.98 * K_o + 0.02 * Na_o, 0.98 * K_i + 0.02 * Na_i),
        IKr ~ GKr * i_OK * (vm - E_Kr),
        1 ~ CK0 + i_CK1 + i_CK2 + i_OK + i_IK,
        D(i_CK1) ~ rc0c1 - rc1c2,
        D(i_CK2) ~ rc1c2 - rc2o,
        D(i_OK) ~ rc2o - roi,
        D(i_IK) ~ roi,
    ]
    return ODESystem(eqs, t; name)
end
