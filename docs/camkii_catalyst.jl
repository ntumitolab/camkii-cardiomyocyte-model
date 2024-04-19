# # Isolated CaMKII response
using Catalyst
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DiffEqCallbacks

# Model
function make_camodel(;
    carest=50e-9,
    cam_total=30e-6, ## Total calmodulin concentration
    camkii_total=70e-6, ## Total CaMKII concentration
    binding_To_PCaMK=0.1,
    decay_CaM=3,
    phospho_rate=1,
    phosphatase=1,
    decay_calcium = 10.0,
)
    ##  Two Ca2+ ions bind to C or N-lobe.
    caMon(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off) = Ca^2 * k_1C_on * k_2C_on / (k_1C_off + k_2C_on * Ca)
    caMoff(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off) = k_1C_off * k_2C_off / (k_1C_off + k_2C_on * Ca)

    rn = @reaction_network begin
        (fCa * CaResting, dCa * dCaRev), 0 <--> Ca

        ##  Two Ca2+ ions bind to C or N-lobe of CaM
        caMon(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), CaM0 --> Ca2CaM_C
        caMoff(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), CaM0 <-- Ca2CaM_C
        caMon(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), CaM0 --> Ca2CaM_N
        caMoff(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), CaM0 <-- Ca2CaM_N
        caMon(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), Ca2CaM_C --> Ca4CaM
        caMoff(Ca, k_1C_on, k_2C_on, k_1C_off, k_2C_off), Ca2CaM_C <-- Ca4CaM
        caMon(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), Ca2CaM_N --> Ca4CaM
        caMoff(Ca, k_1N_on, k_2N_on, k_1N_off, k_2N_off), Ca2CaM_N <-- Ca4CaM

        ##  Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII complex
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMK --> Ca2CaM_C_CaMK
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMK <-- Ca2CaM_C_CaMK
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMK --> Ca2CaM_N_CaMK
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMK <-- Ca2CaM_N_CaMK
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMK --> Ca4CaM_CaMK
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMK <-- Ca4CaM_CaMK
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMK --> Ca4CaM_CaMK
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMK <-- Ca4CaM_CaMK

        ##  Binding of Ca to CaM-CaMKIIP.
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMKP --> Ca2CaM_C_CaMKP
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), CaM0_CaMKP <-- Ca2CaM_C_CaMKP
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMKP --> Ca2CaM_N_CaMKP
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), CaM0_CaMKP <-- Ca2CaM_N_CaMKP
        caMon(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMKP --> Ca4CaM_CaMKP
        caMoff(Ca, k_K1C_on, k_K2C_on, k_K1C_off, k_K2C_off), Ca2CaM_C_CaMKP <-- Ca4CaM_CaMKP
        caMon(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMKP --> Ca4CaM_CaMKP
        caMoff(Ca, k_K1N_on, k_K2N_on, k_K1N_off, k_K2N_off), Ca2CaM_N_CaMKP <-- Ca4CaM_CaMKP

        ##  Binding of CaM to CaMKII or CaMII-P
        (kCaM0_on, kCaM0_off), CaM0 + CaMK <--> CaM0_CaMK
        (kCaM2C_on, kCaM2C_off), Ca2CaM_C + CaMK <--> Ca2CaM_C_CaMK
        (kCaM2N_on, kCaM2N_off), Ca2CaM_N + CaMK <--> Ca2CaM_N_CaMK
        (kCaM4_on, kCaM4_off), Ca4CaM + CaMK <--> Ca4CaM_CaMK
        (kCaM0P_on, kCaM0P_off), CaM0 + CaMKP <--> CaM0_CaMKP
        (kCaM2CP_on, kCaM2CP_off), Ca2CaM_C + CaMKP <--> Ca2CaM_C_CaMKP
        (kCaM2NP_on, kCaM2NP_off), Ca2CaM_N + CaMKP <--> Ca2CaM_N_CaMKP
        (kCaM4P_on, kCaM4P_off), Ca4CaM + CaMKP <--> Ca4CaM_CaMKP

        ## Phosphorylation of CaMKII
        k_phosCaM * (CaMKP + CaMKP2 + CaM0_CaMK + Ca2CaM_C_CaMK + Ca2CaM_N_CaMK + Ca4CaM_CaMK + CaM0_CaMKP + Ca2CaM_C_CaMKP + Ca2CaM_N_CaMKP + Ca4CaM_CaMKP) / CaMKII_T, (Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK) --> (Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP)

        ## Dephosphorylation
        k_dephospho, CaMKP --> CaMK
        ## Second phosphorylation
        (k_P1_P2, k_P2_P1), CaMKP <--> CaMKP2
    end

    setdefaults!(rn, [
        :Ca => carest,
        :CaM0 => cam_total,
        :Ca2CaM_C => 0.0,
        :Ca2CaM_N => 0.0,
        :Ca4CaM => 0.0,
        :CaM0_CaMK => 0.0,
        :Ca2CaM_C_CaMK => 0.0,
        :Ca2CaM_N_CaMK => 0.0,
        :Ca4CaM_CaMK => 0.0,
        :CaM0_CaMKP => 0.0,
        :Ca2CaM_C_CaMKP => 0.0,
        :Ca2CaM_N_CaMKP => 0.0,
        :Ca4CaM_CaMKP => 0.0,
        :CaMK => camkii_total,
        :CaMKP => 0.0,
        :CaMKP2 => 0.0,
        :dCa => decay_calcium,
        :dCaRev => 1.0,
        :fCa => decay_calcium,
        :CaResting => carest,
        :k_1C_on => 5e6, ## 1.2-9.6uM-1s-1
        :k_1C_off => 50, ## 10-70 s-1
        :k_2C_on => 10e6, ## 5-25uM-1s-1.
        :k_2C_off => 10, ## 8.5-10s-1.

        ## N-lobe
        :k_1N_on => 100e6, ## 25-260uM-1s-1
        :k_1N_off => 2000, ## 1000-4000 s-1
        :k_2N_on => 200e6, ## 50-300uM-1s-1.
        :k_2N_off => 500, ## 500->1000.s-1

        ## Ca2+ binding to CaM-CAMKII(K)
        ## C-lobe
        :k_K1C_on => 44e6,
        :k_K1C_off => 33,
        :k_K2C_on => 44e6,
        :k_K2C_off => 0.8, ## 0.49-4.9s-1

        ## N-lobe
        :k_K1N_on => 76e6,
        :k_K1N_off => 300,
        :k_K2N_on => 76e6,
        :k_K2N_off => 20, ## 6-60-1
        :kCaM0_on => 3.8e3,
        :kCaM2C_on => 0.92e6,
        :kCaM2N_on => 0.12e6,
        :kCaM4_on => 30e6,
        :kCaM0_off => 5.5,
        :kCaM2C_off => 6.8,
        :kCaM2N_off => 1.7,
        :kCaM4_off => 1.5,
        :kCaM0P_on => 3.8e3 * binding_To_PCaMK,
        :kCaM2CP_on => 0.92e6 * binding_To_PCaMK,
        :kCaM2NP_on => 0.12e6 * binding_To_PCaMK,
        :kCaM4P_on => 30e6 * binding_To_PCaMK,
        :kCaM0P_off => 1 / decay_CaM,
        :kCaM2CP_off => 1 / decay_CaM,
        :kCaM2NP_off => 1 / decay_CaM,
        :kCaM4P_off => 1 / decay_CaM,
        :k_phosCaM => 30 * phospho_rate,
        :k_dephospho => (1 / 6) * phosphatase,
        :k_P1_P2 => 1 / 60,
        :k_P2_P1 => (1 / 6) * 0.25,
        :CaMKII_T => camkii_total,
    ])
    return rn
end

# ## First sumulation
rn = make_camodel()
@unpack Ca = rn

# Adding calcium
# Difference from the original model:
# Calcium pulses started from 100 sec to 250 sec.
function make_ca_events(;
    starttime=100,
    step=1/3,
    endtime=250,
    peakCa = 4e-6
    )
    affect! = (integrator) -> begin
        t = integrator.t - starttime
        i = t / step
        strength = peakCa * 0.5 * (1 + exp(-i / 5))
        integrator[Ca] += strength
        set_proposed_dt!(integrator, 0.01)
    end
    return PresetTimeCallback(starttime:step:endtime, affect!)
end

# Solve the problem
tspan = (0.0, 400.0)
oprob = ODEProblem(rn, [], tspan)
alg = TRBDF2()
sol = solve(oprob, alg, callback=make_ca_events())

# Calcium
plot(sol, idxs=Ca, tspan=(100, 150))

# Components
@unpack CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2 = rn
plot(
    sol,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600)
)

# Active CaMKII
plot(sol, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="3Hz")
plot!(sol, idxs=Ca*10, plotdensity=100_000, label="Ca * 10")

# Change frequency to 2Hz
sol2 = solve(oprob, alg, callback=make_ca_events(step=1/2), progress=true)

# Calcium
plot(sol2, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol2,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600),
)

# Active CaMKII
plot(sol2, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="2Hz")
plot!(sol2, idxs=Ca*10, plotdensity=100_000, label="Ca * 10")

# Change frequency to 1Hz
sol3 = solve(oprob, alg, callback=make_ca_events(step=1), progress=true)

# Calcium
plot(sol3, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol3,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600),
)

# Active CaMKII
plot(sol3, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="1Hz")
plot!(sol3, idxs=Ca*10, plotdensity=100_000, label="Ca * 10")

# Frequency-dependent response
plot(sol, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), title="Act. CaMKII", label="3Hz")
plot!(sol2, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="2Hz")
plot!(sol3, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="1Hz")

# ## More realistic calcium pulses
# Exponential growth at the rising phase
# and exponential decay at the falling phase

function make_ca_waves(rn;
    starttime=100,
    period=1/3,
    endtime=250,
    peakCa = 4e-6,
    strength = 2,
    )

    @unpack Ca, dCaRev = rn
    caidx = findfirst(isequal(Ca), states(rn))

    rise! = (integrator) -> begin
        integrator.ps[dCaRev] *= -strength
        set_proposed_dt!(integrator, 0.01)
    end

    fallcond = (u, t, integrator) -> begin
        u[caidx] - peakCa
    end

    fallaffect! = (integrator) -> begin
        integrator.ps[dCaRev] = 1.0
        set_proposed_dt!(integrator, 0.01)
    end

    riseevents = PresetTimeCallback(starttime:period:endtime, rise!)
    fallevents = ContinuousCallback(fallcond, fallaffect!)
    return CallbackSet(riseevents, fallevents)
end

# Solve the problem
cawave3hz = make_ca_waves(rn)
sol = solve(oprob, alg, callback=cawave3hz)
plot(sol, idxs=Ca, tspan=(200, 201))

# Components
@unpack CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2 = rn
plot(
    sol,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600)
)

# Active CaMKII
plot(sol, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="3Hz")
plot!(sol, idxs=Ca*10, plotdensity=100_000, label="Ca * 10")

# Change frequency to 2Hz
sol2 = solve(oprob, alg, callback=make_ca_waves(rn; period=1/2))

# Calcium
plot(sol2, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol2,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600),
)

# Active CaMKII
plot(sol2, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="2Hz")
plot!(sol2, idxs=Ca*10, plotdensity=100_000, label="Ca * 10")

# Change frequency to 1Hz
sol3 = solve(oprob, alg, callback=make_ca_waves(rn; period=1))

# Calcium
plot(sol3, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol3,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600),
)

# Active CaMKII
plot(sol3, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="1Hz")
plot!(sol3, idxs=Ca*10, plotdensity=100_000, label="Ca * 10")

# Frequency-dependent response
plot(sol, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), title="Act. CaMKII", label="3Hz")
plot!(sol2, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="2Hz")
plot!(sol3, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="1Hz")
