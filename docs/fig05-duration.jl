# # Fig5: Pacing durations
# Experiments vs simulations in a range of pacing durations

using Model
using Model: second, Hz
using CSV
using CurveFit
using DataFrames
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using Plots
using SteadyStateDiffEq
Plots.default(lw=1.5)

# ## Experiment
# 30 sec resting + N sec 1Hz pacing + resting.
durationdf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-duration-normalized.csv"), DataFrame)
ts = durationdf[!, "Time(sec)"]
fifteen = durationdf[!, "1Hz 15sec (Mean)"]
fifteen_error = durationdf[!, "1Hz 15sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 15sec (N)"])
thirty = durationdf[!, "1Hz 30sec (Mean)"]
thirty_error = durationdf[!, "1Hz 30sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 30sec (N)"])
sixty = durationdf[!, "1Hz 60sec (Mean)"]
sixty_error = durationdf[!, "1Hz 60sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 60sec (N)"])
ninety = durationdf[!, "1Hz 90sec (Mean)"]
ninety_error = durationdf[!, "1Hz 90sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 90sec (N)"])

# Note that 30 sec timeseries have a lower baseline and 90 sec timeseries have a higher baseline before pacing.
fig5a = plot(ts, fifteen, yerr=fifteen_error, lab="15 sec", color=:blue, markerstrokecolor=:blue)
plot!(fig5a, ts, thirty, yerr=thirty_error, lab="30 sec", color=:red, markerstrokecolor=:red)
plot!(fig5a, ts, sixty, yerr=sixty_error, lab="60 sec", color=:orange, markerstrokecolor=:orange)
plot!(fig5a, ts, ninety, yerr=ninety_error, lab="90 sec", color=:green, markerstrokecolor=:green)
plot!(fig5a, title="A", xlabel="Time (s)", ylabel="CaMKAR (R/R0)", titlelocation = :left)

# ## Simulation
# Following the experimental protocol, we simulate 30 sec resting + N sec 1Hz pacing + resting.
@time "Build system" sys = Model.DEFAULT_SYS
tend = 210second
@time "Build problem" prob = ODEProblem(sys, [], (0second, tend))
alg = KenCarp47()

stimstart = 30second
@unpack Istim, CaMKAct = sys
idxs = (sys.t / 1000, CaMKAct)

@time "Solve problems" sols = map((15, 30, 60, 90)) do dur
    cb = build_stim_callbacks(Istim, stimstart + dur * second; period=1second, starttime=stimstart)
    sol = solve(prob, alg; callback=cb)
end

fig5b = plot()
for (sol, dur, color) in zip(sols, (15, 30, 60, 90), (:blue, :red, :orange, :green))
    plot!(fig5b, sol, idxs=idxs, tspan=(0second, 205second), lab="$dur sec", color=color)
end
plot!(fig5b, title="B", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left)

# Phosphorylated CaMKII fraction
idxs_phos = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2))
fig5c = plot()
for (sol, dur, color) in zip(sols, (15, 30, 60, 90), (:blue, :red, :orange, :green))
    plot!(fig5c, sol, idxs=idxs_phos, tspan=(0second, 205second), lab="$dur sec", color=color)
end
plot!(fig5c, title="C", xlabel="Time (s)", ylabel="Phosphorylated CaMKII fraction", titlelocation = :left)

# ## Decay rates
# Fit data from experiments and simulations against an exponential decay model.
# Record for 50 seconds after pacing ends.
ts = collect(range(0.0, stop=50.0, step=5.0)) ## in seconds

ydata_15 = fifteen[10:20]
ydata_30 = thirty[13:23]
ydata_60 = sixty[19:29]
ydata_90 = ninety[25:35]

# Simulation points
ysim_15 = sols[1](stimstart+15second:5second:stimstart+15second+50second ; idxs=sys.CaMKAct * 100).u
ysim_30 = sols[2](stimstart+30second:5second:stimstart+30second+50second ; idxs=sys.CaMKAct * 100).u
ysim_60 = sols[3](stimstart+60second:5second:stimstart+60second+50second ; idxs=sys.CaMKAct * 100).u
ysim_90 = sols[4](stimstart+90second:5second:stimstart+90second+50second ; idxs=sys.CaMKAct * 100).u

# Fit
fit_15 = solve(CurveFitProblem(ts, ydata_15), ExpSumFitAlgorithm(n=1, withconst=true))
fit_30 = solve(CurveFitProblem(ts, ydata_30), ExpSumFitAlgorithm(n=1, withconst=true))
fit_60 = solve(CurveFitProblem(ts, ydata_60), ExpSumFitAlgorithm(n=1, withconst=true))
fit_90 = solve(CurveFitProblem(ts, ydata_90), ExpSumFitAlgorithm(n=1, withconst=true))

fit_sim_15 = solve(CurveFitProblem(ts, ysim_15), ExpSumFitAlgorithm(n=1, withconst=true))
fit_sim_30 = solve(CurveFitProblem(ts, ysim_30), ExpSumFitAlgorithm(n=1, withconst=true))
fit_sim_60 = solve(CurveFitProblem(ts, ysim_60), ExpSumFitAlgorithm(n=1, withconst=true))
fit_sim_90 = solve(CurveFitProblem(ts, ysim_90), ExpSumFitAlgorithm(n=1, withconst=true))

# ### Fitting results (experiments)
p1 = plot(ts, ydata_15, label="Exp 15 sec")
plot!(p1, ts, predict(fit_15), label="Fit", linestyle=:dash)
p2 = plot(ts, ydata_30, label="Exp 30 sec")
plot!(p2, ts, predict(fit_30), label="Fit", linestyle=:dash)
p3 = plot(ts, ydata_60, label="Exp 60 sec")
plot!(p3, ts, predict(fit_60), label="Fit", linestyle=:dash)
p4 = plot(ts, ydata_90, label="Exp 90 sec")
plot!(p4, ts, predict(fit_90), label="Fit", linestyle=:dash)
plot(p1, p2, p3, p4, layout=(2,2), xlabel="Time (s)", ylabel="CaMKII activity (AU)")

# ### Fitting results (simulations)
p1s = plot(ts, ysim_15, label="Sim 15 sec")
plot!(p1s, ts, predict(fit_sim_15), label="Fit", linestyle=:dash)
p2s = plot(ts, ysim_30, label="Sim 30 sec")
plot!(p2s, ts, predict(fit_sim_30), label="Fit", linestyle=:dash)
p3s = plot(ts, ysim_60, label="Sim 60 sec")
plot!(p3s, ts, predict(fit_sim_60), label="Fit", linestyle=:dash)
p4s = plot(ts, ysim_90, label="Sim 90 sec")
plot!(p4s, ts, predict(fit_sim_90), label="Fit", linestyle=:dash)
plot(p1s, p2s, p3s, p4s, layout=(2,2), xlabel="Time (s)", ylabel="CaMKII activity (%)")

# ### Decay time scales (tau)
tau_exp_15 = inv(-fit_15.u.λ[])
tau_exp_30 = inv(-fit_30.u.λ[])
tau_exp_60 = inv(-fit_60.u.λ[])
tau_exp_90 = inv(-fit_90.u.λ[])
tau_sim_15 = inv(-fit_sim_15.u.λ[])
tau_sim_30 = inv(-fit_sim_30.u.λ[])
tau_sim_60 = inv(-fit_sim_60.u.λ[])
tau_sim_90 = inv(-fit_sim_90.u.λ[])

println("The time scales for experiments: ")
for (tau, dur) in zip((tau_exp_15, tau_exp_30, tau_exp_60, tau_exp_90), (15, 30, 60, 90))
    println("$dur sec pacing is $(round(tau; digits=2)) seconds.")
end

println("The time scales for simulations: ")
for (tau, dur) in zip((tau_sim_15, tau_sim_30, tau_sim_60, tau_sim_90), (15, 30, 60, 90))
    println("$dur sec pacing is $(round(tau; digits=2)) seconds.")
end

# Plot pacing time vs decay time scale
pacing_durations = [15.0, 30.0, 60.0, 90.0]
tau_experiments = [tau_exp_15, tau_exp_30, tau_exp_60, tau_exp_90]
tau_simulations = [tau_sim_15, tau_sim_30, tau_sim_60, tau_sim_90]
fig5d = plot(pacing_durations, tau_experiments, label="Experiments", marker=:circle, color=:blue)
plot!(fig5d, pacing_durations, tau_simulations, label="Simulations", marker=:square, color=:red)
plot!(fig5d, title="D", xlabel="Pacing Duration (s)", ylabel="Decay Time Scale (s)", titlelocation = :left)

# ## Combine all panels
fig5 = plot(fig5a, fig5b, fig5c, fig5d, layout=(2, 2), size=(800, 600))

#---
savefig(fig5, "fig5.png")
savefig(fig5, "fig5.pdf")
