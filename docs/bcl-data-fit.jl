# # Pacing data
# Experiments vs simulations.
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CSV
using DataFrames
using LsqFit
using CaMKIIModel
using CaMKIIModel: second
Plots.default(lw=1.5)

#---
@time "Build system" sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)
tend = 500second
@time "Build problem" prob = ODEProblem(sys, [], tend)
stimstart = 100second
stimend = 300second
@unpack Istim = sys
alg = KenCarp47()

#===
## Pacing duration and CaMKII activity

### Experiments

30 seconds resting + N seconds 1Hz pacing + resting.
===#
durationdf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-duration.csv"), DataFrame)
ts = durationdf[!, "Time(sec)"]
fifteen = durationdf[!, "1Hz 15sec (Mean)"]
fifteen_error = durationdf[!, "1Hz 15sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 15sec (N)"])
thirty = durationdf[!, "1Hz 30sec (Mean)"] .+ 0.25
thirty_error = durationdf[!, "1Hz 30sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 30sec (N)"])
sixty = durationdf[!, "1Hz 60sec (Mean)"]
sixty_error = durationdf[!, "1Hz 60sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 60sec (N)"])
ninety = durationdf[!, "1Hz 90sec (Mean)"] .- 0.25
ninety_error = durationdf[!, "1Hz 90sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 90sec (N)"])

## 30 sec timeseries +0.25 and 90 sec timeseries -0.25 for a consistent baseline before pacing.
plot(ts, fifteen, yerr=fifteen_error, lab="15 sec", color=:blue, markerstrokecolor=:blue)
plot!(ts, thirty, yerr=thirty_error, lab="30 sec", color=:red, markerstrokecolor=:red)
plot!(ts, sixty, yerr=sixty_error, lab="60 sec", color=:orange, markerstrokecolor=:orange)
plot!(ts, ninety, yerr=ninety_error, lab="90 sec", color=:green, markerstrokecolor=:green)
plot!(title="Experiment", xlabel="Time (s)", ylabel="CaMKII activity (AU)")

#---
savefig("pacing-duration-exp.pdf")
savefig("pacing-duration-exp.png")

# ### Simulation
stimstart = 30second
callback15 = build_stim_callbacks(Istim, stimstart + 15second; period=1second, starttime=stimstart)
sol15 = solve(prob, alg; callback=callback15)
callback30 = build_stim_callbacks(Istim, stimstart + 30second; period=1second, starttime=stimstart)
sol30 = solve(prob, alg; callback=callback30)
callback60 = build_stim_callbacks(Istim, stimstart + 60second; period=1second, starttime=stimstart)
sol60 = solve(prob, alg; callback=callback60)
callback90 = build_stim_callbacks(Istim, stimstart + 90second; period=1second, starttime=stimstart)
sol90 = solve(prob, alg; callback=callback90)
idxs = (sys.t / 1000, sys.CaMKAct * 100)

plot(sol15, idxs=idxs, tspan=(0second, 205second), lab="15 sec", color=:blue)
plot!(sol30, idxs=idxs, tspan=(0second, 205second), lab="30 sec", color=:red)
plot!(sol60, idxs=idxs, tspan=(0second, 205second), lab="60 sec", color=:orange)
plot!(sol90, idxs=idxs, tspan=(0second, 205second), lab="90 sec", color=:green)
plot!(title="Simulation", xlabel="Time (s)", ylabel="CaMKII activity (%)")

#===
### Decay rates

Fit against an exponential decay model.
===#
@. decay_model(x, p) = p[1] * exp(-x / p[2]) + p[3]
rmse(fit) = sqrt(sum(abs2, fit.resid) / length(fit.resid))

xdata_15 = ts[10:20] .- ts[10]
ydata_15 = fifteen[10:20]
xdata_30 = ts[13:23] .- ts[13]
ydata_30 = thirty[13:23]
xdata_60 = ts[19:29] .- ts[19]
ydata_60 = sixty[19:29]
xdata_90 = ts[25:35] .- ts[25]
ydata_90 = ninety[25:35]

p0 = [1.0, 2.0, 13.0]

@time fit_15 = curve_fit(decay_model, xdata_15, ydata_15, p0, autodiff=:forwarddiff)
@time fit_30 = curve_fit(decay_model, xdata_30, ydata_30, p0, autodiff=:forwarddiff)
@time fit_60 = curve_fit(decay_model, xdata_60, ydata_60, p0, autodiff=:forwarddiff)
@time fit_90 = curve_fit(decay_model, xdata_90, ydata_90, p0, autodiff=:forwarddiff)

#---
@show rmse(fit_15)
@show rmse(fit_30)
@show rmse(fit_60)
@show rmse(fit_90);

# Simulation decay fits

#---
println("The time scale for: ")
for (fit, dur) in zip((fit_15, fit_30, fit_60, fit_90), (15, 30, 60, 90))
    println("  $dur sec pacing is $(round(coef(fit)[2]; digits=2)) seconds.")
end

#---
plot([15, 30, 60, 90], [coef(fit_15)[2], coef(fit_30)[2], coef(fit_60)[2], coef(fit_90)[2]], xlabel="Pacing Duration (s)", ylabel="Decay Time Scale (s)", title="Decay Time Scale vs Pacing Duration", marker=:circle, label=false, color=:black, xlims=(0, 100))

#---
savefig("pacing-decay-exp.pdf")
savefig("pacing-decay-exp.png")



#---
savefig("pacing-duration-sim.pdf")

# ### Phosphorylated fraction
idxs = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2) * 100)
plot(sol15, idxs=idxs, tspan=(0second, 205second), lab="15 sec", color=:blue)
plot!(sol30, idxs=idxs, tspan=(0second, 205second), lab="30 sec", color=:red)
plot!(sol60, idxs=idxs, tspan=(0second, 205second), lab="60 sec", color=:orange)
plot!(sol90, idxs=idxs, tspan=(0second, 205second), lab="90 sec", color=:green)
plot!(title="Simulation", xlabel="Time (s)", ylabel="Phosphorylated CaMKII (%)")

#---
savefig("pacing-duration-phos.pdf")

#===
## Pacing frequency and CaMKII activity

### Experiments
===#
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

# ### Simulations
tend = 205.0second
prob = ODEProblem(sys, [], tend)
stimstart = 30.0second
stimend = 120.0second
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
sol1 = solve(prob, alg; callback)

callback2 = build_stim_callbacks(Istim, stimend; period=0.5second, starttime=stimstart)
sol2 = solve(prob, alg; callback=callback2)
idxs = (sys.t / 1000, sys.CaMKAct * 100)

plot(sol1, idxs=idxs, lab="1 Hz", color=:blue)
plot!(sol2, idxs=idxs, lab="2 Hz", color=:red)
plot!(title="Simulation", xlabel="Time (s)", ylabel="CaMKII activity (%)")

#---
savefig("pacing-frequency-sim.pdf")

#---
idxs = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2) * 100)
plot(sol1, idxs=idxs, lab="1 Hz", color=:blue)
plot!(sol2, idxs=idxs, lab="2 Hz", color=:red)
plot!(title="Simulation", xlabel="Time (s)", ylabel="Phosphorylated CaMKII (%)")

#---
savefig("pacing-frequency-phos.pdf")
