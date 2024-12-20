# # Caffeine Effects
# Caffeine increase RyR opening sensitivity to luminal and subspace calcium
# In this model, we decrease the mid saturation sub-SR calcium concentration for the opening rate
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel
using CaMKIIModel: second, metre, Farad
Plots.default(lw=1.5)

#---
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)
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
@time sol = solve(prob, alg; callback)
@time sol_caf = solve(prob, alg; callback=callback_caf)

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

# Add caffeine in the beginning of the simulation
# Add caffeine and nifedipine in the beginning of the simulation (nifedipine blocks 90% of L-type calcium channel)
prob = ODEProblem(sys, [], tend)
prob_caf = ODEProblem(sys, [sys.RyRsensitivity => 10], tend)
prob_nif_caf = ODEProblem(sys, [sys.RyRsensitivity => 10, sys.GCaL => 6.3e-6 * (metre^3 / second / Farad)], tend)
@time sol = solve(prob, alg; callback)
@time sol_caf = solve(prob_caf, alg; callback)
@time sol_nif_caf = solve(prob_nif_caf, alg; callback)

#---
plot(sol, idxs=sys.vm, title="Action potential", lab="Ctl")
plot!(sol_caf, idxs=sys.vm, lab="Caf")
plot!(sol_nif_caf, idxs=sys.vm, lab="Caf + Nif", tspan=(198second, 205second), ylabel="Voltage (mV)", xlabel="Time (ms)")

#---
plot(sol, idxs=sys.Cai_sub_SR * 1000, title="Calcium transient", lab="Ctl")
plot!(sol_caf, idxs=sys.Cai_sub_SR * 1000, lab="Caf", ylabel="Subspace calcium (nM)")
plot!(sol_nif_caf, idxs=sys.Cai_sub_SR * 1000, tspan=(198second, 205second), lab="Caf + Nif", ylabel="Subspace calcium (nM)", xlabel="Time (ms)")


#---
plot(sol, idxs=sys.PO1RyR, title="RyR opening", lab="Ctl")
plot!(sol_caf, idxs=sys.PO1RyR, lab="Caf")
plot!(sol_nif_caf, idxs=sys.PO1RyR, tspan=(198second, 205second), lab="Caf + Nif", ylabel="Open probability", ylims=(0, 1), xlabel="Time (ms)")

#---
plot(sol, idxs=sys.CaJSR, title="SR Calcium", lab="Ctl", ylabel="SR calcium (μM)")
plot!(sol_caf, idxs=sys.CaJSR, lab="Caf")
plot!(sol_nif_caf, idxs=sys.CaJSR, tspan=(198second, 205second), lab="Caf + Nif", ylims=(0, 850), xlabel="Time (ms)")

#---
plot(sol, idxs=sys.CaMKAct*100, title="CaMKII", lab="Ctl")
plot!(sol_caf, idxs=sys.CaMKAct*100, lab="Caf")
plot!(sol_nif_caf, idxs=sys.CaMKAct*100, lab="Caf + Nif", ylabel="Active fraction (%)", xlabel="Time (ms)")
