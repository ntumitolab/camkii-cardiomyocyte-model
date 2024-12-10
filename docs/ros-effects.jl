# # ROS effects
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using CaMKIIModel
using CaMKIIModel: second, μM
Plots.default(lw=1.5)

# ## Setup system
# Electrical stimulation starts at `t`=100 seconds and ends at `t`=300 seconds.
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)
tend = 500.0second
prob = ODEProblem(sys, [], tend)
stimstart = 100.0second
stimend = 300.0second
@unpack Istim = sys
alg = FBDF()

# ## No ROS
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, alg; callback)

#---
plot(sol, idxs=sys.vm, tspan=(298second, 300second), title="Action potential")

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299second, 300second), title="Calcium transcient")

#---
plot(sol, idxs=sys.CaMKAct*100, title="CaMKII", xlabel="Time (ms)", ylabel="Active fraction (%)", label=false)

# ## ROS 0.1uM
prob2 = remake(prob, p=[sys.ROS => 0.1μM])
@time sol2 = solve(prob2, alg; callback)

#---
plot(sol2, idxs=sys.vm, tspan=(298second, 300second), title="Action potential")

#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(298second, 300second), title="Calcium transcient")

#---
plot(sol2, idxs=sys.CaMKAct*100, title="CaMKII", xlabel="Time (ms)", ylabel="Active fraction (%)",label=false)

# ## ROS 1uM
prob3 = remake(prob, p=[sys.ROS => 1μM])
@time sol3 = solve(prob3, alg; callback)

#---
plot(sol3, idxs=sys.vm, tspan=(298second, 300second), title="Action potential")

#---
plot(sol3, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(298second, 300second), title="Calcium transcient")

#---
plot(sol3, idxs=sys.CaMKAct*100, title="Active CaMKII", label=false)

# ## Comparisons
plot(sol, idxs=sys.CaMKAct*100, title="Active CaMKII", lab="ROS=0")
plot!(sol2, idxs=sys.CaMKAct*100, lab="ROS=0.1uM")
plot!(sol3, idxs=sys.CaMKAct*100, lab="ROS=1uM", xlabel="Time (ms)", ylabel="Active fraction (%)")
