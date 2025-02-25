# # Caffeine Effects
# Caffeine increase RyR opening sensitivity to luminal and subspace calcium
# In this model, we decrease the mid saturation sub-SR calcium concentration for the opening rate
using ModelingToolkit
using DifferentialEquations
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
alg = FBDF()
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
i = (sys.t/1000, sys.vm)
plot(sol, idxs=i, title="Action potential", lab="Ctl", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, lab="Caf", tspan=(198second, 205second), ylabel="Voltage (mV)", xlabel="Time (s)")

#---
i = (sys.t/1000, sys.Cai_sub_SR * 1000)
plot(sol, idxs=i, title="Calcium transient (During caffeine addition)", lab="Ctl", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, tspan=(198second, 205second), lab="Caf", ylabel="Subspace calcium (nM)", xlabel="Time (s)")

#---
i = (sys.t/1000, sys.PO1RyR)
plot(sol, idxs=i, title="RyR open (During caffeine addition)", lab="Ctl", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, tspan=(198second, 205second), lab="Caf", ylabel="Open probability", ylims=(0, 1), xlabel="Time (s)")

#---
i = (sys.t/1000, sys.Cai_sub_SR * 1000)
plot(sol, idxs=i, title="Calcium transient (After caffeine addition)", lab="Ctl", ylabel="Subspace calcium (nM)", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, lab="Caf", xlabel="Time (s)", tspan=(198second, 205second))

#---
i = (sys.t/1000, sys.CaJSR)
plot(sol, idxs=i, title="SR Calcium (During caffeine addition)", lab="Ctl", ylabel="SR calcium (μM)", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, tspan=(198second, 205second), lab="Caf", ylims=(0, 850), xlabel="Time (s)")

#---
i = (sys.t/1000, sys.Jrel)
plot(sol, idxs=sys.Jrel, title="Ca flux", lab="Ctl  (Jrel)", tspan=(198second, 205second))
plot!(sol_caf, idxs=sys.Jrel, lab="Caf (Jrel)", tspan=(198second, 205second), ylabel="μM/ms", xlabel="Time (s)")

#---
i = (sys.t/1000, sys.CaMKAct*100)
plot(sol, idxs=i, title="Active CaMKII", lab="Ctl")
plot!(sol_caf, idxs=i, lab="Caf", ylabel="CaMKII activity (%)", xlabel="Time (s)")

# Add caffeine in the beginning of the simulation
# Add caffeine and nifedipine in the beginning of the simulation (nifedipine blocks 90% of L-type calcium channel)
prob = ODEProblem(sys, [], tend)
prob_caf = ODEProblem(sys, [sys.RyRsensitivity => 10], tend)
prob_nif_caf = ODEProblem(sys, [sys.RyRsensitivity => 10, sys.GCaL => 6.3e-6 * (metre^3 / second / Farad)], tend)
@time sol = solve(prob, alg; callback)
@time sol_caf = solve(prob_caf, alg; callback)
@time sol_nif_caf = solve(prob_nif_caf, alg; callback)

#---
i = (sys.t/1000, sys.vm)
tspan = (198second, 205second)
plot(sol, idxs=i, title="Action potential", lab="Ctl"; tspan)
plot!(sol_caf, idxs=i, lab="Caf"; tspan)
plot!(sol_nif_caf, idxs=i, lab="Caf + Nif", tspan=tspan, ylabel="Voltage (mV)", xlabel="Time (s)")

#---
i = (sys.t/1000, sys.Cai_sub_SR * 1000)
tspan=(198second, 205second)
plot(sol, idxs=i, title="Calcium transient", lab="Ctl";tspan)
plot!(sol_caf, idxs=i, lab="Caf", ylabel="Subspace calcium (nM)";tspan)
plot!(sol_nif_caf, idxs=i, lab="Caf + Nif", ylabel="Subspace calcium (nM)", xlabel="Time (s)";tspan)

#---
i = (sys.t/1000, sys.CaJSR)
tspan=(198second, 205second)
plot(sol, idxs=i, title="SR Calcium", lab="Ctl", ylabel="SR calcium (μM)"; tspan)
plot!(sol_caf, idxs=i, lab="Caf"; tspan)
plot!(sol_nif_caf, idxs=i, lab="Caf + Nif", ylims=(0, 850), xlabel="Time (s)"; tspan)

#---
i = (sys.t/1000, sys.CaMKAct*100)
plot(sol, idxs=i, title="CaMKII", lab="Ctl")
plot!(sol_caf, idxs=i, lab="Caf")
plot!(sol_nif_caf, idxs=i, lab="Caf + Nif", ylabel="Active fraction (%)", xlabel="Time (s)")
