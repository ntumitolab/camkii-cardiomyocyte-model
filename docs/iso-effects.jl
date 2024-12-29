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
alg = TRBDF2()

# ## Without isoproterenol
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, alg; callback)

#---
i = (sys.t/1000, sys.vm)
plot(sol, idxs=i, tspan=(295second, 300second), title="Action potential", xlabel="Time (s)")

#---
plot(sol, idxs=(sys.t/1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(299second, 300second), title="Calcium transcient", xlabel="Time (s)", ylabel="Conc. (μM)", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"])

#---
plot(sol, idxs=(sys.t/1000, sys.CaMKAct*100), title="Active CaMKII", label=false, ylabel="Active fraction (%)" , xlabel="Time (s)")

# ## 0.1uM isoproterenol
prob2 = remake(prob, p=[sys.ISO => 0.1μM])
sol2 = solve(prob2, alg; callback)

#---
plot(sol2, idxs=(sys.t/1000, sys.vm), tspan=(295second, 300second), title="Action potential", xlabel="Time (s)")

#---
plot(sol2, idxs=(sys.t/1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(299second, 300second), title="Calcium transcient", xlabel="Time (s)", ylabel="Conc. (μM)", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"])

#---
plot(sol2, idxs=(sys.t/1000, sys.CaMKAct*100), title="Active CaMKII", label=false, ylabel="Active fraction (%)" , xlabel="Time (s)")

# ## Comparison
i = (sys.t/1000, sys.Cai_mean * 1000)
tspan = (299second, 300second)
plot(sol, idxs=i, title="Calcium transcient", lab="ISO (0uM)"; tspan)
plot!(sol2, idxs=i, lab="ISO (0.1uM)", xlabel="Time (s)", ylabel="Conc. (μM)"; tspan)

#---
i = (sys.t/1000, sys.CaMKAct*100)
plot(sol, idxs=i, title="Active CaMKII", lab="ISO (0uM)")
plot!(sol2, idxs=i, lab="ISO (0.1uM)", ylabel="Active fraction (%)" , xlabel="Time (s)")

#---
i = (sys.t/1000, sys.vm)
plot(sol, idxs=i, tspan=(297second, 300second), title="Action potential", lab="ISO (0uM)")
plot!(sol2, idxs=i, tspan=(297second, 300second), lab="ISO (1uM)", xlabel="Time (ms)", ylabel="Voltage (mV)")
