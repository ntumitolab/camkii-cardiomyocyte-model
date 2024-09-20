#===
# Period vs active CaMKII
===#
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel

sys = build_neonatal_ecc_sys(simplify=true)
tend = 300.0
prob = ODEProblem(sys, [], tend)

# 1Hz
@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend; period=1)
alg = FBDF()
@time sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol, idxs=sys.vm, tspan=(295, 300), title="Action potential")
#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transcient")
#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

# 2Hz
callback = build_stim_callbacks(Istim, tend; period=1/2)
@time sol2 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol2, idxs=sys.vm, tspan=(298, 300), title="Action potential")
#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transcient")
#---
plot(sol2, idxs=sys.CaMKAct, title="Active CaMKII")

# 3Hz
callback = build_stim_callbacks(Istim, tend; period=1/3)
@time sol3 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol3, idxs=sys.vm, tspan=(298, 300), title="Action potential")
# Na channel not fully recovered
plot(sol3, idxs=[sys.i_Nam, sys.i_Nah, sys.i_Naj], tspan=(299, 300), title="Sodium gating")
#---
plot(sol3, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transcient")
#---
plot(sol3, idxs=sys.CaMKAct, title="Active CaMKII")

# Comparing 1-3 Hz
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII", lab="1Hz")
plot!(sol2, idxs=sys.CaMKAct, lab="2Hz")
plot!(sol3, idxs=sys.CaMKAct, lab="3Hz", ylim=(0.0, 0.8))
