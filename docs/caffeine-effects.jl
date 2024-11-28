# # Caffeine Effects
# Caffeine increase RyR opening sensitivity to luminal and subspace calcium
# In this model, we decrease the mid saturation sub-SR calcium concentration for the opening rate
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
alg = TRBDF2()
function add_coffee_affect!(integrator)
    integrator.ps[sys.RyRsensitivity] = 10
end

@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
# Add caffeine at t = 200
callback_caf = CallbackSet(build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart), PresetTimeCallback(200.0, add_coffee_affect!))

#---
sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)
sol_caf = solve(prob, alg; callback=callback_caf, abstol=1e-6, reltol=1e-6)

#---
plot(sol, idxs=sys.vm * 1000, title="Action potential", lab="Ctl")
plot!(sol_caf, idxs=sys.vm * 1000, lab="Caf", tspan=(195, 210), ylabel="Voltage (mV)")

#---
plot(sol, idxs=sys.Cai_sub_SR * 1E6, title="Calcium transient (During caffeine addition)", lab="Ctl")
plot!(sol_caf, idxs=sys.Cai_sub_SR * 1E6, tspan=(195, 210), lab="Caf", ylabel="Subspace calcium (nM)")

#---
plot(sol, idxs=sys.PO1RyR, title="RyR open (During caffeine addition)", lab="Ctl")
plot!(sol_caf, idxs=sys.PO1RyR, tspan=(195, 210), lab="Caf", ylabel="Open probability", ylims=(0, 1))

#---
plot(sol, idxs=sys.Cai_sub_SR * 1E6, title="Calcium transient (After caffeine addition)", lab="Ctl", ylabel="Subspace calcium (nM)")
plot!(sol_caf, idxs=sys.Cai_sub_SR * 1E6, tspan=(295, 300), lab="Caf")

#---
plot(sol, idxs=sys.CaJSR * 1E3, title="SR Calcium (During caffeine addition)", lab="Ctl", ylabel="SR calcium (μM)")
plot!(sol_caf, idxs=sys.CaJSR * 1E3, tspan=(195, 210), lab="Caf", ylims=(0, 850))

#---
plot(sol, idxs=sys.Jrel, title="Ca flux", lab="Ctl  (Jrel)")
plot!(sol_caf, idxs=sys.Jrel, lab="Caf (Jrel)", ylabel="μM/ms")
plot!(sol, idxs=sys.Jup, lab="Ctl (Jup)")
plot!(sol_caf, idxs=sys.Jup, lab="Caf (Jup)", tspan=(195, 210), ylabel="μM/ms")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII", lab="Ctl")
plot!(sol_caf, idxs=sys.CaMKAct, lab="Caf")
