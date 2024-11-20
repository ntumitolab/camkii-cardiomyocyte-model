# # Caffeine Effects
# Caffeine increase RyR opening sensitivity to luminal and subspace calcium
# In this model, we decrease the mid saturation sub-SR calcium concentration for the opeing rate
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel
Plots.default(lw=1.5)

#---
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 500.0
prob = ODEProblem(sys, [], tend)
stimstart = 100.0
stimend = 300.0
alg = FBDF()
function add_coffee_affect!(integrator)
    integrator.ps[sys.RyRsensitivity] = 10
end

@unpack Istim = sys
callback = CallbackSet(build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart), PresetTimeCallback(200.0, add_coffee_affect!))

#---
sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol, idxs=sys.vm*1000, tspan=(199, 210), title="Action potential", lab=false)

#---
plot(sol, idxs=sys.Cai_sub_SR*1E6, tspan=(195, 210), title="Calcium transient", lab=false, ylabel="Subspace calcium (nM)")

#---
plot(sol, idxs=sys.Cai_sub_SR*1E6, tspan=(295, 300), title="Calcium transient", lab=false, ylabel="Subspace calcium (nM)")

#---
plot(sol, idxs=[sys.CaNSR*1E3, sys.CaJSR*1E3], tspan=(195, 210), title="SR Calcium", lab=["NSR" "JSR"], ylabel="SR calcium (μM)")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

#---
plot(sol, idxs=sys.CaMKAct, tspan=(195, 210), title="Active CaMKII")
