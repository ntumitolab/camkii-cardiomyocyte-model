# # Pacing frequency
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CSV
using DataFrames
using CurveFit
using CaMKIIModel
using CaMKIIModel: second
Plots.default(lw=1.5)

# ## Experiments
freqdf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-freq.csv"), DataFrame)
ts = 0:5:205
onehz = freqdf[!, "1Hz (Mean)"]
onehz_error = freqdf[!, "1Hz (SD)"] ./ sqrt.(freqdf[!, "1Hz (N)"])
twohz = freqdf[!, "2Hz (Mean)"]
twohz_error = freqdf[!, "2Hz (SD)"] ./ sqrt.(freqdf[!, "2Hz (N)"])

plot(ts, onehz, yerr=onehz_error, lab="1 Hz", color=:blue, markerstrokecolor=:blue)
plot!(ts, twohz, yerr=twohz_error, lab="2 Hz", color=:red, markerstrokecolor=:red)
plot!(title="Experiment", xlabel="Time (s)", ylabel="CaMKII activity (AU)")

#---
savefig("pacing-frequency-exp.pdf")

# ## Simulations
@time "Build system" @mtkcompile sys = build_neonatal_ecc_sys()
tend = 205.0second
@time "Build problem" prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
stimstart = 30.0second
stimend = 120.0second
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)

alg = KenCarp47()
@time "Solve problem" sol1 = solve(prob, alg; callback)

callback2 = build_stim_callbacks(Istim, stimend; period=0.5second, starttime=stimstart)
@time "Solve problem" sol2 = solve(prob, alg; callback=callback2)
idxs = (sys.t / 1000, sys.CaMKAct * 100)

plot(sol1, idxs=idxs, lab="1 Hz", color=:blue)
plot!(sol2, idxs=idxs, lab="2 Hz", color=:red)
plot!(title="Simulation", xlabel="Time (s)", ylabel="CaMKII activity (%)")

#---
savefig("pacing-frequency-sim.pdf")

# ## Phosphorylated fraction
idxs = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2) * 100)
plot(sol1, idxs=idxs, lab="1 Hz", color=:blue)
plot!(sol2, idxs=idxs, lab="2 Hz", color=:red)
plot!(title="Simulation", xlabel="Time (s)", ylabel="Phosphorylated CaMKII (%)")

#---
savefig("pacing-frequency-phos.pdf")
