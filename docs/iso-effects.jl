# # Effects of isoproterenol
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using CaMKIIModel
using CaMKIIModel: second, μM
Plots.default(lw=1.5)

#---
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)
tend = 500.0second
prob = ODEProblem(sys, [], tend)
stimstart = 100.0second
stimend = 300.0second
alg = FBDF()

# ## Without isoproterenol
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, alg; callback)

#---
plot(sol, idxs=sys.vm, tspan=(295second, 300second), title="Action potential", xlabel="Time (ms)")

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299second, 300second), title="Calcium transcient", xlabel="Time (ms)", ylabel="Conc. (μM)")

#---
plot(sol, idxs=sys.CaMKAct*100, title="Active CaMKII", label=false, ylabel="Active fraction (%)" , xlabel="Time (ms)")

# ## 1uM isoproterenol
prob2 = remake(prob, p=[sys.ISO => 1μM])
sol2 = solve(prob2, alg; callback)

#---
plot(sol2, idxs=sys.vm, tspan=(295second, 300second), title="Action potential", xlabel="Time (ms)")

#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299second, 300second), title="Calcium transcient", xlabel="Time (ms)", ylabel="Conc. (μM)")

#---
plot(sol2, idxs=sys.CaMKAct*100, title="Active CaMKII", label=false, ylabel="Active fraction (%)" , xlabel="Time (ms)")

# ## Comparison
plot(sol, idxs=sys.Cai_mean, title="Calcium transcient", lab="ISO (0uM)")
plot!(sol2, idxs=sys.Cai_mean, tspan=(299second, 300second), lab="ISO (1uM)", xlabel="Time (ms)", ylabel="Conc. (μM)")

#---
plot(sol, idxs=sys.CaMKAct*100, title="Active CaMKII", lab="ISO (0uM)")
plot!(sol2, idxs=sys.CaMKAct*100, lab="ISO (1uM)", ylabel="Active fraction (%)" , xlabel="Time (ms)")

#---
plot(sol, idxs=sys.vm, tspan=(297second, 300second), title="Action potential", lab="ISO (0uM)")
plot!(sol2, idxs=sys.vm, tspan=(297second, 300second), lab="ISO (1uM)", xlabel="Time (ms)", ylabel="Voltage (mV)")
