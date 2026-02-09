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

# ## Decay rates
# Fit against an exponential decay model.
decay_model(p, x) = @. p[1] * exp(-x / p[2]) + p[3]

# Data from experiments
# Record 50 seconds after pacing ends
ts = collect(range(0.0, stop=50.0, step=5.0))
ydata_1hz = onehz[24:34]
ydata_2hz = twohz[24:34]

# Simulation points
ysim_1hz = sol1(stimend:5second:stimend+50second ; idxs=sys.CaMKAct * 100).u
ysim_2hz = sol2(stimend:5second:stimend+50second ; idxs=sys.CaMKAct * 100).u

# Fit data to an exponential decay model
fit_1hz = solve(CurveFitProblem(ts, ydata_1hz), ExpSumFitAlgorithm(n=1, withconst=true))
fit_2hz = solve(CurveFitProblem(ts, ydata_2hz), ExpSumFitAlgorithm(n=1, withconst=true))
fit_1hz_sim = solve(CurveFitProblem(ts, ysim_1hz), ExpSumFitAlgorithm(n=1, withconst=true))
fit_2hz_sim = solve(CurveFitProblem(ts, ysim_2hz), ExpSumFitAlgorithm(n=1, withconst=true))

# Fitting results (experiments)
p1 = plot(ts, ydata_1hz, label="Exp 1 Hz")
plot!(p1, ts, predict(fit_1hz), label="Fit", linestyle=:dash)
p2 = plot(ts, ydata_2hz, label="Exp 2 Hz")
plot!(p2, ts, predict(fit_2hz), label="Fit", linestyle=:dash)
p3 = plot(ts, ysim_1hz, label="Sim 1 Hz")
plot!(p3, ts, predict(fit_1hz_sim), label="Fit", linestyle=:dash)
p4 = plot(ts, ysim_2hz, label="Sim 2 Hz")
plot!(p4, ts, predict(fit_2hz_sim), label="Fit", linestyle=:dash)
plot(p1, p2, p3, p4, layout=(2,2), xlabel="Time (s)", ylabel="CaMKII activity (AU)")


# Calculate time scales (tau) from fit parameters
tau_exp_1hz = inv(-fit_1hz.u.λ[])
tau_exp_2hz = inv(-fit_2hz.u.λ[])
tau_sim_1hz = inv(-fit_1hz_sim.u.λ[])
tau_sim_2hz = inv(-fit_2hz_sim.u.λ[])

println("The time scales for experiments: ")
for (tau, freq) in zip((tau_exp_1hz, tau_exp_2hz), (1, 2))
    println("$freq Hz pacing is $(round(tau; digits=2)) seconds.")
end
println("The time scales for simulations: ")
for (tau, freq) in zip((tau_sim_1hz, tau_sim_2hz), (1, 2))
    println("$freq Hz pacing is $(round(tau; digits=2)) seconds.")
end
