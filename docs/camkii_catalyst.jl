# # Isolated CaMKII response
using Catalyst
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DiffEqCallbacks
using CaMKIIModel: nM, μM, Hz, get_camkii_rn

# Exponential decay calcium model
@parameters begin
    CaResting = 50nM
    dCa = 10.0Hz
    dCaRev = 1
end
@variables t Ca(t)
D = Differential(t)
eq = [D(Ca) ~ -dCa * dCaRev * (Ca - CaResting * dCaRev)]
@named osys = ODESystem(eq, t)

# ## First simulation
@parameters ROS=0μM
rn = get_camkii_rn(Ca, ROS)

## TODO: Why Ca(t) is eleiminated?
rnsys = convert(ODESystem, rn, remove_conserved=true)



@named connsys = ODESystem([osys.Ca ~ rnsys.Ca], t)

@named sys = compose(ODESystem([osys.Ca ~ rnsys.Ca], t, name=:composed), [])

@named odesys = extend(casys, rnsys)

observed(odesys)

# Calcium pulses
function make_ca_events(;
    starttime=200,
    step=1 / 3,
    endtime=300,
    peakCa=4μM
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
oprob = ODEProblem(odesys, [], tspan)
alg = TRBDF2()
sol = solve(oprob, alg, callback=make_ca_events())

# Calcium
plot(sol, idxs=Ca, tspan=(100, 150))

# Components
@unpack CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, CaMKOX, CaMKII_T, CaMKPOX = rn

camkiistates= [CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, CaMKPOX, CaMKOX]

plot(
    sol,
    idxs=camkiistates,
    size=(800, 600)
)

camk2act = sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, CaMKPOX, CaMKOX]) / CaMKII_T

# Active CaMKII
plot(sol, idxs=camk2act, label="Act. CaMKII", title="3Hz", ylims=(0.0, 1.0))

# Change frequency to 2Hz
sol2 = solve(oprob, alg, callback=make_ca_events(step=1 / 2))

# Calcium
plot(sol2, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol2,
    idxs=camkiistates,
    size=(800, 600),
)

# Active CaMKII
plot(sol2, idxs=camk2act, label="Act. CaMKII", title="2Hz", ylims=(0.0, 1.0))

# Change frequency to 1Hz
sol3 = solve(oprob, alg, callback=make_ca_events(step=1))

# Calcium
plot(sol3, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol3,
    idxs=camkiistates,
    size=(800, 600),
)

# Active CaMKII
plot(sol3, idxs=camk2act, label="Act. CaMKII", title="1Hz", ylims=(0.0, 1.0))

# Frequency-dependent response
plot(sol, idxs=camk2act, title="Act. CaMKII", label="3Hz")
plot!(sol2, idxs=camk2act, label="2Hz")
plot!(sol3, idxs=camk2act, label="1Hz", ylims=(0.0, 1.0))

# ## More realistic calcium pulses
# Exponential growth at the rising phase
# and exponential decay at the falling phase
function make_ca_waves(rn;
    starttime=100,
    period=1 / 3,
    endtime=250,
    peakCa=4e-6,
    strength=5,
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

plot(
    sol,
    idxs=camkiistates,
    size=(800, 600)
)

# Active CaMKII
plot(sol, idxs=camk2act, label="Act. CaMKII", title="3Hz", ylims=(0, 1))

# Change frequency to 2Hz
cawave2hz = make_ca_waves(rn; period=1 / 2)
sol2 = solve(oprob, alg, callback=cawave2hz)

# Calcium
plot(sol2, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol2,
    idxs=camkiistates,
    size=(800, 600),
)

# Active CaMKII
plot(sol2, idxs=camk2act, label="Act. CaMKII", title="2Hz", ylims=(0, 1))

# Change frequency to 1Hz
cawave1hz = make_ca_waves(rn; period=1)
sol3 = solve(oprob, alg, callback=cawave1hz)

# Calcium
plot(sol3, idxs=Ca, tspan=(100, 150))

# Components
plot(
    sol3,
    idxs=camkiistates,
    size=(800, 600),
)

# Active CaMKII
plot(sol3, idxs=camk2act, label="Act. CaMKII", title="1Hz", ylims=(0, 1))

# Frequency-dependent response
plot(sol, idxs=camk2act, title="Act. CaMKII", label="3Hz")
plot!(sol2, idxs=camk2act, label="2Hz")
plot!(sol3, idxs=camk2act, label="1Hz", ylims=(0, 1))

# ## Effect of ROS
tspan = (0.0, 400.0)
osys = convert(ODESystem, rn, remove_conserved=true)
oprob_ROS0 = ODEProblem(osys, [ROS=>0.0], tspan)
oprob_ROS1 = ODEProblem(osys, [ROS=>1.0], tspan)
oprob_ROS2 = ODEProblem(osys, [ROS=>2.0], tspan)
alg = TRBDF2()

sol_ROS0 = solve(oprob_ROS0, alg, callback=cawave1hz)
sol_ROS1 = solve(oprob_ROS1, alg, callback=cawave1hz)
sol_ROS2 = solve(oprob_ROS2, alg, callback=cawave1hz)

plot(
    sol_ROS0,
    idxs=camkiistates,
    size=(800, 600)
)

plot(
    sol_ROS1,
    idxs=camkiistates,
    size=(800, 600)
)

#---
plot(sol_ROS0, idxs=camk2act, label="ROS = 0.0", title="Act. CaMKII (1Hz)")
plot!(sol_ROS1, idxs=camk2act, label="ROS = 1.0")
plot!(sol_ROS2, idxs=camk2act, label="ROS = 2.0", ylim=(0.0, 1.0), size=(400, 400))

savefig("ROS.png")
