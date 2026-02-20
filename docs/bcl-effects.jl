# # Pacing simulations
# Responses of calcium transients and electrophysiology to pacing frequencies.

using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CSV
using StatsBase: mean
using DataFrames
using CaMKIIModel
using CaMKIIModel: ms
Plots.default(lw=1.5)

# ## Setup the ODE system
# Electrical stimulation starts at `t`=100 sec and ends at `t`=300 sec.
@time @mtkcompile sys = build_neonatal_ecc_sys()
tend = 500.0*1000ms
prob = ODEProblem(sys, [], tend)
stimstart = 100.0*1000ms
stimend = 300.0*1000ms
@unpack Istim = sys
alg = KenCarp47()

# ## Single pulse
callback = build_stim_callbacks(Istim, stimstart + 1*1000ms; period=10*1000ms, starttime=stimstart)

@time sol = solve(prob, alg; callback)

plot(sol, idxs=(sys.t / 1000 - 100, sys.vm), title="Action potential (single pulse)", ylabel="mV", xlabel="Time (s)", label=false, tspan=(100*1000ms, 103*1000ms))

#---
savefig("single-pulse.pdf")

#---
plot(sol, idxs=(sys.t / 1000 - 100, sys.Cai_mean), tspan=(100*1000ms, 103*1000ms), title="Calcium transient", ylabel="Conc. (μM)", xlabel="Time (s)", label="Avg Ca")

#---
savefig("single-cat.pdf")

# ## 1Hz pacing
callback = build_stim_callbacks(Istim, stimend; period=1*1000ms, starttime=stimstart)
@time sol = solve(prob, alg; callback);

plot(sol, idxs=(sys.t / 1000, sys.vm), title="Action potential", ylabel="mV", xlabel="Time (s)", label=false)

#---
plot(sol, idxs=(sys.t / 1000 - 299, sys.vm), title="Action potential", tspan=(299*1000ms, 300*1000ms), ylabel="mV", xlabel="Time (s)", label=false)

#---
plot(sol, idxs=(sys.t / 1000 - 299, [sys.IK1, sys.Ito, sys.IKs, sys.IKr, sys.If]), tspan=(299*1000ms, 300*1000ms), ylabel="μA/μF", xlabel="Time (s)", label=["IK1" "Ito" "IKs" "IKr" "If"])

#---
plot(sol, idxs=(sys.t / 1000 - 299, [sys.ICaL, sys.INaCa, sys.ICaT, sys.ICab]), tspan=(299*1000ms, 300*1000ms), ylabel="μA/μF", xlabel="Time (s)", label=["ICaL" "INaCa" "ICaT" "ICab"])

#---
plot(sol, idxs=(sys.t / 1000 - 299, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(299*1000ms, 300*1000ms), title="Calcium transient", ylabel="μM", xlabel="Time (s)", label=["CaSR" "CaSL" "CaAvg"])

#---
plot(sol, idxs=(sys.t / 1000, sys.CaMKAct * 100), title="Active CaMKII", ylabel="Active CaMKII (%)", xlabel="Time (s)", label=false)

# ### 3D surface plot
xx = 1:44
yy = range(299*1000ms, 300*1000ms, length=100)
zz = [sol(t, idxs=sys.Cai[u]) for t in yy, u in xx]

surface(xx, (yy .- 299*1000ms) ./ 1000, zz, colorbar=:none, yguide="(sec.)", zguide="Ca Conc. (μM)", xticks=false, size=(600, 600))
annotate!(3, 0, 0.65, "Ca (SL)")
annotate!(41, 0.25, 0.58, "Ca (SR)")

#---
savefig("3d-surface.pdf")

# ### Mean calcium
xx = 1:44
yy = [sol(299.22*1000ms, idxs=sys.Cai[u]) for u in xx]
avg = mean(yy)

plot(xx, yy, label="Cai", ylabel="Ca Conc. (μM)", xlabel="Compartment (SL to SR)", legend=:topright)
hline!([avg], linestyle=:dash, label="Avg Cai")

#---
savefig("cai-spatial.pdf")

# ## 2Hz pacing
callback = build_stim_callbacks(Istim, stimend; period=0.5*1000ms, starttime=stimstart)
@time sol2 = solve(prob, alg; callback)

plot(sol2, idxs=(sys.t / 1000, sys.vm), title="Action potential", ylabel="mV", xlabel="Time (s)", label=false)

#---
plot(sol2, idxs=(sys.t / 1000 - 299, sys.vm), title="Action potential", tspan=(299*1000ms, 300*1000ms), ylabel="mV", xlabel="Time (s)", label=false)

#---
plot(sol2, idxs=(sys.t / 1000 - 299, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(299*1000ms, 300*1000ms), title="Calcium transient", ylabel="Concentration (μM)", xlabel="Time (s)", label=["CaSSR" "CaSL" "CaAvg"])

#---
plot(sol2, idxs=(sys.t / 1000, sys.CaMKAct * 100), title="Active CaMKII", ylabel="Active CaMKII (%)", xlabel="Time (s)", label=false)

# ## Comparing 1 and 2 Hz pacing
idxs = (sys.t / 1000 - 299, sys.vm)
plot(sol, idxs=idxs, title="Action potential", lab="1Hz", tspan=(299*1000ms, 300*1000ms))
plot!(sol2, idxs=idxs, lab="2Hz", tspan=(299*1000ms, 300*1000ms), xlabel="Time (s)", ylabel="Voltage (mV)")

#---
savefig("bcl-ap.pdf")

#---
idxs = (sys.t / 1000 - 299, sys.Cai_mean)
plot(sol, idxs=idxs, title="Calcium transient", lab="1Hz", tspan=(299*1000ms, 300*1000ms))
plot!(sol2, idxs=idxs, lab="2Hz", tspan=(299*1000ms, 300*1000ms), xlabel="Time (s)", ylabel="Concentration (μM)")

#---
savefig("bcl-cat.pdf")

#---
idxs = (sys.t / 1000, sys.CaMKAct * 100)
plot(sol, idxs=idxs, title="CaMKII", lab="1Hz")
plot!(sol2, idxs=idxs, lab="2Hz", xlabel="Time (s)", ylabel="Active fraction (%)")

#---
savefig("bcl-camkact.pdf")
