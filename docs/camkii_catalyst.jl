# # Isolated CaMKII response
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel: nM, μM, Hz, get_camkii_eqs

# ## Setup model
# Exponential decay calcium model
function ca_decay_eqs(;
    carest=50nM,
    decay_calcium=10.0Hz,
)
    @parameters begin
        CaResting = carest
        dCa = decay_calcium
        dCaRev = 1
    end
    @variables t Ca(t)
    D = Differential(t)
    eqs = [D(Ca) ~ -dCa * dCaRev * (Ca - CaResting * dCaRev)]
    return eqs
end

# CaMKII model
@parameters ROS=0μM
@variables t Ca(t)
eqs = get_camkii_eqs(Ca, ROS)
caeqs = ca_decay_eqs()
sys = ODESystem([eqs; caeqs], t, name=:sys) |> structural_simplify

tspan = (0.0, 400.0)
prob = ODEProblem(sys, [Ca => sys.CaResting], tspan)

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

# ## First simulation

# Solve the problem
tspan = (0.0, 400.0)
alg = TRBDF2()
sol = solve(prob, alg, callback=make_ca_events())

# Calcium
plot(sol, idxs=Ca, tspan=(200, 300))

# Components
@unpack CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, CaMKOX, CaMKII_act, CaMKPOX = sys

camkiistates= [CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2, Ca4CaM_CaMKOX, Ca4CaM_CaMKPOX, CaMKPOX, CaMKOX]

plot(sol, idxs=camkiistates, size=(800, 600))

# Active CaMKII
plot(sol, idxs=CaMKII_act, label="Act. CaMKII", title="3Hz", ylims=(0.0, 1.0))

# Change frequency to 2Hz
sol2 = solve(prob, alg, callback=make_ca_events(step=1 / 2))

# Calcium
plot(sol2, idxs=Ca, tspan=(200, 300))

# Components
plot(sol2,idxs=camkiistates, size=(800, 600), legends=:topleft)

# Active CaMKII
plot(sol2, idxs=CaMKII_act, label="Act. CaMKII", title="2Hz", ylims=(0.0, 1.0))

# Change frequency to 1Hz
@time sol3 = solve(prob, TRBDF2(), callback=make_ca_events(step=1))

# Calcium
plot(sol3, idxs=Ca, tspan=(200, 300))

# Components
plot(sol3, idxs=camkiistates, size=(800, 600), legend=:topleft)

# Active CaMKII
plot(sol3, idxs=CaMKII_act, label="Act. CaMKII", title="1Hz", ylims=(0.0, 1.0))

# Frequency-dependent response
plot(sol, idxs=CaMKII_act, title="Act. CaMKII", label="3Hz")
plot!(sol2, idxs=CaMKII_act, label="2Hz")
plot!(sol3, idxs=CaMKII_act, label="1Hz", ylims=(0.0, 1.0))

# ## More realistic calcium pulses
# Exponential growth at the rising phase
# and exponential decay at the falling phase
function make_ca_waves(sys;
    starttime=200,
    period=1 / 3,
    endtime=300,
    peakCa=4μM,
    strength=5,
)

    @unpack Ca, dCaRev = sys
    caidx = findfirst(isequal(Ca), unknowns(sys))

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
cawave3hz = make_ca_waves(sys)
sol = solve(prob, alg, callback=cawave3hz)
plot(sol, idxs=Ca, tspan=(200, 201))

# Components
plot(sol, idxs=camkiistates, size=(800, 600), legend=:topleft)

# Active CaMKII
plot(sol, idxs=CaMKII_act, label="Act. CaMKII", title="3Hz", ylims=(0, 1))

# Change frequency to 2Hz
cawave2hz = make_ca_waves(sys; period=1 / 2)
sol2 = solve(prob, alg, callback=cawave2hz)

# Calcium
plot(sol2, idxs=Ca, tspan=(200, 201))

# Components
plot(sol2, idxs=camkiistates, size=(800, 600), legend=:topleft)

# Active CaMKII
plot(sol2, idxs=CaMKII_act, label="Act. CaMKII", title="2Hz", ylims=(0, 1))

# Change frequency to 1Hz
cawave1hz = make_ca_waves(sys; period=1)
sol3 = solve(prob, alg, callback=cawave1hz)

# Calcium
plot(sol3, idxs=Ca, tspan=(200, 202))

# Components
plot(sol3,idxs=camkiistates, size=(800, 600),legend=:topleft)

# Active CaMKII
plot(sol3, idxs=CaMKII_act, label="Act. CaMKII", title="1Hz", ylims=(0, 1))

# Frequency-dependent response
plot(sol, idxs=CaMKII_act, title="Act. CaMKII", label="3Hz")
plot!(sol2, idxs=CaMKII_act, label="2Hz")
plot!(sol3, idxs=CaMKII_act, label="1Hz", ylims=(0, 1))

# ## Effect of ROS
oprob_ROS0 = remake(prob, p=[ROS=>0.0])
oprob_ROS1 = remake(prob, p=[ROS=>1.0μM])
oprob_ROS2 = remake(prob, p=[ROS=>2.0μM])

sol_ROS0 = solve(oprob_ROS0, alg, callback=cawave1hz)
sol_ROS1 = solve(oprob_ROS1, alg, callback=cawave1hz)
sol_ROS2 = solve(oprob_ROS2, alg, callback=cawave1hz)

plot(sol_ROS0, idxs=camkiistates, size=(800, 600),legend=:topleft)

plot(sol_ROS1, idxs=camkiistates, size=(800, 600),legend=:topleft)

#---
plot(sol_ROS0, idxs=CaMKII_act, label="ROS = 0.0μM", title="Act. CaMKII (1Hz)")
plot!(sol_ROS1, idxs=CaMKII_act, label="ROS = 1.0μM")
plot!(sol_ROS2, idxs=CaMKII_act, label="ROS = 2.0μM", ylim=(0.0, 1.0), size=(600, 600))
