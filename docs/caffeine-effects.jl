# # Caffeine Effects
# Caffeine increase RyR opening sensitivity to luminal and subspace calcium
# In this model, we decrease the mid saturation sub-SR calcium concentration for the opening rate
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel
using CaMKIIModel: second
Plots.default(lw=1.5)

#---
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 500.0second
prob = ODEProblem(sys, [], tend)
stimstart = 100.0second
stimend = 300.0second
alg = TRBDF2()
function add_coffee_affect!(integrator)
    integrator.ps[sys.RyRsensitivity] = 10
end

@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
# Add caffeine at t = 200 econd
callback_caf = CallbackSet(build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart), PresetTimeCallback(200.0second, add_coffee_affect!))

#---
sol = solve(prob, alg; callback)
sol_caf = solve(prob, alg; callback=callback_caf)

#---
plot(sol, idxs=sys.vm, title="Action potential", lab="Ctl")
plot!(sol_caf, idxs=sys.vm, lab="Caf", tspan=(198second, 205second), ylabel="Voltage (mV)")

#---
plot(sol, idxs=sys.Cai_sub_SR * 1000, title="Calcium transient (During caffeine addition)", lab="Ctl")
plot!(sol_caf, idxs=sys.Cai_sub_SR * 1000, tspan=(198second, 205second), lab="Caf", ylabel="Subspace calcium (nM)")

#---
plot(sol, idxs=sys.PO1RyR, title="RyR open (During caffeine addition)", lab="Ctl")
plot!(sol_caf, idxs=sys.PO1RyR, tspan=(198second, 205second), lab="Caf", ylabel="Open probability", ylims=(0, 1), xlabel="Time (ms)")

#---
plot(sol, idxs=sys.Cai_sub_SR * 1000, title="Calcium transient (After caffeine addition)", lab="Ctl", ylabel="Subspace calcium (nM)")
plot!(sol_caf, idxs=sys.Cai_sub_SR * 1000, tspan=(198second, 205second), lab="Caf", xlabel="Time (ms)")

#---
plot(sol, idxs=sys.CaJSR, title="SR Calcium (During caffeine addition)", lab="Ctl", ylabel="SR calcium (μM)")
plot!(sol_caf, idxs=sys.CaJSR, tspan=(198second, 205second), lab="Caf", ylims=(0, 850), xlabel="Time (ms)")

#---
plot(sol, idxs=sys.Jrel, title="Ca flux", lab="Ctl  (Jrel)")
plot!(sol_caf, idxs=sys.Jrel, lab="Caf (Jrel)", tspan=(198second, 205second), ylabel="μM/ms", xlabel="Time (ms)")

#---
plot(sol, idxs=sys.CaMKAct*100, title="Active CaMKII", lab="Ctl")
plot!(sol_caf, idxs=sys.CaMKAct*100, lab="Caf", ylabel="CaMKII activity (%)", xlabel="Time (ms)")
