# # Effects of isoproterenol
using ModelingToolkit
using OrdinaryDiffEq
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

# ## Without isoproterenol
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
@time sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol, idxs=sys.vm, tspan=(295, 300), title="Action potential")

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transcient")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

# ## 1uM isoproterenol
prob2 = remake(prob, p=[sys.ISO => 1E-3])
@time sol2 = solve(prob2, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol2, idxs=sys.vm, tspan=(295, 300), title="Action potential")

#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transcient")

#---
plot(sol2, idxs=sys.CaMKAct, title="Active CaMKII")


# ## Comparison
plot(sol, idxs=sys.Cai_mean, tspan=(299, 300), title="Calcium transcient", lab="ISO (0uM)")
plot!(sol2, idxs=sys.Cai_mean, tspan=(299, 300), lab="ISO (1uM)")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII", lab="ISO (0uM)")
plot!(sol2, idxs=sys.CaMKAct, lab="ISO (1uM)")

#---
plot(sol, idxs=sys.vm, tspan=(297, 300), title="Action potential", lab="ISO (0uM)")
plot!(sol2, idxs=sys.vm, tspan=(297, 300), lab="ISO (1uM)")
