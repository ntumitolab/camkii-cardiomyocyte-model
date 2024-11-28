# # ROS effects
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using CaMKIIModel
Plots.default(lw=1.5)

# ## Setup system
# Electrical stimulation starts at `t`=100 seconds and ends at `t`=300 seconds.
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 500.0
prob = ODEProblem(sys, [], tend)
stimstart = 100.0
stimend = 300.0
@unpack Istim = sys
alg = TRBDF2()

# ## No ROS
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol, idxs=sys.vm, tspan=(298, 300), title="Action potential")

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transcient")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

# ## ROS 1uM
prob2 = remake(prob, p=[sys.ROS => 1e-3])
sol2 = solve(prob2, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol2, idxs=sys.vm, tspan=(298, 300), title="Action potential")

#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(298, 300), title="Calcium transcient")

#---
plot(sol2, idxs=sys.CaMKAct, title="Active CaMKII")

# ## ROS 5uM
prob3 = remake(prob, p=[sys.ROS => 5e-3])
sol3 = solve(prob3, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol3, idxs=sys.vm, tspan=(298, 300), title="Action potential")

#---
plot(sol3, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(298, 300), title="Calcium transcient")

#---
plot(sol3, idxs=sys.CaMKAct, title="Active CaMKII")

# ## Comparisons
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII", lab="ROS=0", ylims=(0, 1))
plot!(sol2, idxs=sys.CaMKAct, lab="ROS=1uM")
plot!(sol3, idxs=sys.CaMKAct, lab="ROS=5uM")
