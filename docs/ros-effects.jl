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
i = (sys.t / 1000, sys.vm)
tspan=(298second, 300second)
plot(sol, idxs=i, title="Action potential"; tspan)

#---
i = (sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean])
tspan=(298second, 300second)
plot(sol, idxs=i, title="Calcium transcient", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"]; tspan)

#---
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol, idxs=i, title="CaMKII", xlabel="Time (s)", ylabel="Active fraction (%)", label=false)

# ## ROS 0.1uM
prob2 = remake(prob, p=[sys.ROS => 0.1μM])
@time sol2 = solve(prob2, alg; callback)

#---
i = (sys.t / 1000, sys.vm)
tspan=(298second, 300second)
plot(sol2, idxs=i, title="Action potential"; tspan)

#---
i = (sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean])
tspan=(298second, 300second)
plot(sol2, idxs=i, title="Calcium transcient", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"]; tspan)

#---
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol2, idxs=i, title="CaMKII", xlabel="Time (s)", ylabel="Active fraction (%)", label=false)

# ## ROS 1uM
prob3 = remake(prob, p=[sys.ROS => 1μM])
@time sol3 = solve(prob3, alg; callback)

#---
i = (sys.t / 1000, sys.vm)
tspan=(298second, 300second)
plot(sol3, idxs=i, title="Action potential"; tspan)

#---
i = (sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean])
tspan=(298second, 300second)
plot(sol3, idxs=i, title="Calcium transcient", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"]; tspan)

#---
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol3, idxs=i, title="CaMKII", xlabel="Time (s)", ylabel="Active fraction (%)", label=false)

# ## Comparisons
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol, idxs=i, title="Active CaMKII", lab="ROS 0uM")
plot!(sol2, idxs=i, lab="ROS 0.1uM")
plot!(sol3, idxs=i, lab="ROS 1uM", xlabel="Time (s)", ylabel="Active fraction (%)")
