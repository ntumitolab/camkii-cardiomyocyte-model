# LCC Markov model based on Mahajan et al. (2008)

using Catalyst
using ModelingToolkit

function lcc_markov_sys(mode=1)

    @variables t v(t) (k(t))[1:6] (kp(t))[1:6] (s(t))[1:2] (sp(t))[1:2] α(t) β(t) Ca(t)

    # LCC Current Fixed Parameters
    @parameters taupo = 1          # [ms] - Time constant of activation
    @parameters TBa = 450          # [ms] - Time constant of slow inactivation
    @parameters s1o = 0.0221
    @parameters k1o = ifelse(mode == 1, 0.03, 0.3)
    @parameters kop = 2.5e-3       # [mM] - Calcium inactivation concentration
    @parameters cpbar = 8e-3       # [mM]
    @parameters tca = 78.0312
    @parameters ICa_scale = 5.25
    @parameters recoveryReduc = 3
    @parameters SSAshift = 0
    @parameters SSIshift = 0

    # Tranisition Rates (20 rates)
    @parameters r1 = 0.3                            # [1/ms] - Opening rate
    @parameters r2 = ifelse(mode == 1, 3, 3 / 8)        # [1/ms] - closing rate
    @parameters s1p = 0.00195                       # [ms] - Inactivation rate
    @parameters k1p = 0.00413                       # [ms] - Inactivation rate
    @parameters k2 = 0.0001                         # [ms] - Inactivation rate
    @parameters k2p = 0.00224                       # [ms] - Inactivation rate

    poss = expit((v + SSAshift) / 8)
    fca = hil(Ca, kop)
    rv = 10 + 4954 * exp(v / 15.6)
    prLCC = 1 - expit((v + 40) / 4)
    psLCC = expit((v + 40 + SSIshift) / 11.32)
    TCa = 0.1 + tca * hil(cpbar, Ca, 2)
    tauCa = (rv - TCa) * prLCC + TCa
    tauBa = (rv - TBa) * prLCC + TBa

    rateEqs = [
        α ~ poss / taupo,
        β ~ (1 - poss) / taupo,
        s[1] ~ s1o * fca,
        sp[1] ~ s1p,
        s[2] ~ s[1] * (k[2] / k[1]) * (r1 / r2),
        sp[2] ~ sp[1] * (kp[2] / kp[1]) * (r1 / r2),
        k[1] ~ k1o * fca,
        kp[1] ~ k1p,
        k[2] ~ k2,
        kp[2] ~ k2p,
        k[3] ~ hil(exp(-(v + 40) / 3)) / 3,
        kp[3] ~ k[3],
        k[4] ~ k[3] * (α / β) * (k[1] / k[2]) * (k[5] / k[6]),
        kp[4] ~ kp[3] * (α / β) * (kp[1] / kp[2]) * (kp[5] / kp[6]),
        k[5] ~ (1 - prLCC) / (tauCa * recoveryReduc),
        kp[5] ~ (1 - prLCC) / (tauBa * recoveryReduc),
        k[6] ~ fca * psLCC / tauCa,
        kp[6] ~ psLCC / tauBa,
    ]

    rn = @reaction_network begin
        ($α, $β), C2Ca <--> C1Ca
        ($r1, $r2), C1Ca <--> OCa
        ($(k[4]), $(k[3])), I2Ca <--> I1Ca
        ($(s[2]), $(s[1])), I1Ca <--> OCa
        ($(kp[4]), $(kp[3])), I2Ba <--> I1Ba
        ($(sp[2]), $(sp[1])), I1Ba <--> OCa
        ($(kp[6]), $(kp[5])), I2Ba <--> C2Ca
        ($(k[6]), $(k[5])), I2Ca <--> C2Ca
        ($(kp[2]), $(kp[1])), I1Ba <--> C1Ca
        ($(k[2]), $(k[1])), I1Ca <--> C1Ca
    end

    osys = convert(ODESystem, rn, remove_conserved=true)
    lccsys = extend(osys, ODESystem(rateEqs, name=:rates))
    return lccsys
end
