# # Fig4: Pacing frequency
# CaMKII activation with 1Hz vs 2Hz pacing

using Model
using Model: second
using CSV
using CurveFit
using DataFrames
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEqSDIRK
using Plots
using SteadyStateDiffEq
Plots.default(lw=1.5)

# ## Experiment
freqdf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-freq-normalized.csv"), DataFrame)
ts = 0:5:205
onehz = freqdf[!, "1Hz (Mean)"]
onehz_error = freqdf[!, "1Hz (SD)"] ./ sqrt.(freqdf[!, "1Hz (N)"])
twohz = freqdf[!, "2Hz (Mean)"]
twohz_error = freqdf[!, "2Hz (SD)"] ./ sqrt.(freqdf[!, "2Hz (N)"])

fig4a = plot(ts, onehz, yerr=onehz_error, lab="1 Hz", color=:blue, markerstrokecolor=:blue)
plot!(fig4a, ts, twohz, yerr=twohz_error, lab="2 Hz", color=:red, markerstrokecolor=:red)
plot!(fig4a, title="A", xlabel="Time (s)", ylabel="CaMKAR (R/R0)", titlelocation = :left)

# ## Simulation
@time "Build system" sys = Model.DEFAULT_SYS
tend = 205.0second
@time "Build problem" prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
stimstart = 30.0second
stimend = 120.0second
alg = FBDF()

callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time "Solve problem" sol1 = solve(prob, alg; callback)

callback2 = build_stim_callbacks(Istim, stimend; period=0.5second, starttime=stimstart)
@time "Solve problem" sol2 = solve(prob, alg; callback=callback2)

idxs = (sys.t / 1000, sys.CaMKAct)
idxs_phos = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2))
fig4b = plot(sol1, idxs=idxs, lab="1 Hz", color=:blue)
plot!(fig4b, sol1, idxs=idxs_phos, lab="1 Hz (phos.)", color=:blue, linestyle=:dash)
plot!(fig4b, sol2, idxs=idxs, lab="2 Hz", color=:red)
plot!(fig4b, sol2, idxs=idxs_phos, lab="2 Hz (phos.)", color=:red, linestyle=:dash)
plot!(fig4b, title="B", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left)

# Action potential
idxs = (sys.t / 1000, sys.vm)
tspan = (100.0second, 101.0second)
fig4c = plot(sol1, idxs=idxs, lab="1 Hz", color=:blue, tspan=tspan)
plot!(fig4c, sol2, idxs=idxs, lab="2 Hz", color=:red, tspan=tspan)
plot!(fig4c, title="C", xlabel="Time (s)", ylabel="Membrane potential (mV)", titlelocation = :left)

# Calcium transient
idxs = (sys.t / 1000, sys.Cai_mean * 1000)
fig4d = plot(sol1, idxs=idxs, lab="1 Hz", color=:blue, tspan=tspan)
plot!(fig4d, sol2, idxs=idxs, lab="2 Hz", color=:red, tspan=tspan)
plot!(fig4d, title="D", xlabel="Time (s)", ylabel="Mean intracellular Ca (nM)", titlelocation = :left)

# Combine all panels
fig4 = plot(fig4a, fig4b, fig4c, fig4d, layout=(2, 2), size=(800, 600))

#---
savefig(fig4, "fig4.png")
savefig(fig4, "fig4.pdf")

# ## Decay rates
# Fit against an exponential decay model.
# Data from experiments: record 50 seconds after pacing ends

ts = collect(range(0.0, stop=50.0, step=5.0))
ydata_1hz = onehz[24:34]
ydata_2hz = twohz[24:34]

ysim_1hz = sol1(stimend:5second:stimend+50second; idxs=sys.CaMKAct).u
ysim_2hz = sol2(stimend:5second:stimend+50second; idxs=sys.CaMKAct).u

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

# Decay time scales (tau) from fit parameters
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
