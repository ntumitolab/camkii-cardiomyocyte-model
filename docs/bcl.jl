# # Pacing frequency
using ModelingToolkit
using DifferentialEquations
using Plots
using CaMKIIModel
Plots.default(lw=1.5)

# Setup system
# Electrical stimulation starts at `t`=100 seconds and ends at `t`=300 seconds
sys = build_neonatal_ecc_sys(simplify=true)
tend = 500.0
prob = ODEProblem(sys, [], tend)
stimstart = 100.0
stimend = 300.0
@unpack Istim = sys
alg = FBDF()

# ## 1Hz
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
@time sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

# ## 2Hz
callback = build_stim_callbacks(Istim, stimend; period=1/2, starttime=stimstart)
@time sol2 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol2, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol2, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

#---
plot(sol2, idxs=sys.CaMKAct, title="Active CaMKII")

# ## 3Hz
callback = build_stim_callbacks(Istim, stimend; period=1/3, starttime=stimstart)
@time sol3 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol3, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol3, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol3, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

#---
plot(sol3, idxs=sys.CaMKAct, title="Active CaMKII")

# ## Comparing 1-3 Hz
# Action potential
plot(sol, idxs=sys.vm*1000, title="Action potential", lab="1Hz")
plot!(sol2, idxs=sys.vm*1000, lab="2Hz")
plot!(sol3, idxs=sys.vm*1000, lab="3Hz", tspan=(299, 300), xlabel="Time (sec.)", ylabel="Voltage (mV)")

# CaMKII activity
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII", lab="1Hz")
plot!(sol2, idxs=sys.CaMKAct, lab="2Hz")
plot!(sol3, idxs=sys.CaMKAct, lab="3Hz", ylim=(0.0, 0.8))
