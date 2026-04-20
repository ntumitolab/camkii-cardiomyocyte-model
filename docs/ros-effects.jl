# # ROS effects
using Model
using Model: second, μM
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using StatsPlots
using CSV
using DataFrames
import Dates
using CurveFit
Plots.default(lw=1.5)

# ## Setup model
@time "Build system" sys = Model.DEFAULT_SYS
@unpack Istim = sys
tend = 205second
stimstart = 30second
stimend = 120second
newkox_CaMK = inv(45second) / 50μM
@time "Build problem" prob = ODEProblem(sys, [sys.kox_CaMK => inv(45second) / 50μM], tend)
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
alg = KenCarp47()

# ## Comparisons
# ROS (H2O2) 0uM vs 50uM vs 200uM
@time "Solve problem" sol = solve(prob, alg; callback)

prob2 = remake(prob, p=[sys.ROS => 50μM])
@time "Solve problem" sol2 = solve(prob2, alg; callback)

prob3 = remake(prob, p=[sys.ROS => 200μM])
@time "Solve problem" sol3 = solve(prob3, alg; callback);

# ## Comparisons
i = (sys.t / 1000, sys.CaMKAct)
plot(sol, idxs=i, lab="ROS (-)", color=:blue)
plot!(sol2, idxs=i, lab="ROS 50uM", color=:red)
plot!(sol3, idxs=i, lab="ROS 200uM", color=:green)
plot!(xlabel="Time (s)", ylabel="Active CaMKII fraction", title="Simulation")

#---
savefig("ros-camkii.pdf")
savefig("ros-camkii.png")
# Oxidized and autophosphorylated fraction
iox = (sys.t / 1000, 100 * (sys.CaMKBOX + sys.CaMKPOX + sys.CaMKAOX + sys.CaMKOX))
plot(sol, idxs=iox, lab="ROS (-), OX", color=:blue)
plot!(sol2, idxs=iox, lab="ROS 50uM, OX", color=:red)
plot!(sol3, idxs=iox, lab="ROS 200uM, OX", color=:green)
iphos = (sys.t / 1000, 100 * (sys.CaMKP + sys.CaMKA + sys.CaMKA2))
plot!(sol, idxs=iphos, lab="ROS (-), phos", color=:blue, linestyle=:dash)
plot!(sol2, idxs=iphos, lab="ROS 50uM, phos", color=:red, linestyle=:dash)
plot!(sol3, idxs=iphos, lab="ROS 200uM, phos", color=:green, linestyle=:dash)
plot!(xlabel="Time (s)", ylabel="CaMKII fraction", title="Simulation")

#---
savefig("ros-camkiiox-pox.pdf")
savefig("ros-camkiiox-pox.png")

# Proportions for oxidized and autophosphorylated fractions in active CaMKII
# Using StatsPlots for grouped bar plot at the end of the stimulation period (120 seconds)
timepoint = 120second
pidx = (sys.CaMKP + sys.CaMKA + sys.CaMKA2) / sys.CaMKAct
oidx = (sys.CaMKBOX + sys.CaMKPOX + sys.CaMKAOX + sys.CaMKOX) / sys.CaMKAct
pvals = [s(timepoint, idxs=pidx) for s in (sol, sol2, sol3)]
ovals = [s(timepoint, idxs=oidx) for s in (sol, sol2, sol3)]
groupedbar(["0uM", "50uM", "200uM"], [pvals ovals], bar_position=:stack, label=["Phosphorylated" "Oxidized"], color=[:blue :red :green], title="Simulations", ylabel="Fraction of active CaMKII", xlabel="ROS concentration", legend=:topleft)
#---
savefig("ros-camkii-groupbar.pdf")
savefig("ros-camkii-groupbar.png")

# ## Experimental data
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical.csv"), DataFrame)
ts = Dates.value.(chemicaldf[!, "Time"]) ./ 10^9
ctl = chemicaldf[!, "Ctrl Mean"]
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"])

#---
ros50 = chemicaldf[!, "H2O2 50uM Mean"]
ros50_error = chemicaldf[!, "H2O2 50uM SD"] ./ sqrt.(chemicaldf[!, "H2O2 50uM N"])
ros200 = chemicaldf[!, "H2O2 200uM Mean"]
ros200_error = chemicaldf[!, "H2O2 200uM SD"] ./ sqrt.(chemicaldf[!, "H2O2 200uM N"])

#---
plot(ts, ctl, yerr=ctl_error, lab="Control", color=:blue, markerstrokecolor=:blue)
plot!(ts, ros50, yerr=ros50_error, lab="H2O2 50uM", color=:red, markerstrokecolor=:red)
plot!(ts, ros200, yerr=ros200_error, lab="H2O2 200uM", color=:green, markerstrokecolor=:green)
plot!(xlabel="Time (s)", ylabel="CaMKII activity (A.U.)", title="Experiment")

#---
savefig("ros-exp.pdf")
savefig("ros-exp.png")

# ## Decay kinectics
# Fit against an exponential decay model.
# `decay_model(p, x) = @. p[1] * exp(-x / p[2]) + p[3]`

# Data from experiments (starts from 120 seconds till the end of the experiment).
xdata = collect(range(0.0, step=5.0, length=18))
ydata_ctl = ctl[25:end]
ydata_50 = ros50[25:end]
ydata_200 = ros200[25:end]
ydata_ctl_sim = sol(stimend:5second:tend, idxs=sys.CaMKAct * 100).u
ydata_01_sim = sol2(stimend:5second:tend, idxs=sys.CaMKAct * 100).u
ydata_05_sim = sol3(stimend:5second:tend, idxs=sys.CaMKAct * 100).u

# Fit data to an exponential decay model
fit_ctl = solve(CurveFitProblem(xdata, ydata_ctl), ExpSumFitAlgorithm(n=1, withconst=true))
fit_50 = solve(CurveFitProblem(xdata, ydata_50), ExpSumFitAlgorithm(n=1, withconst=true))
fit_200 = solve(CurveFitProblem(xdata, ydata_200), ExpSumFitAlgorithm(n=1, withconst=true))
fit_ctl_sim = solve(CurveFitProblem(xdata, ydata_ctl_sim), ExpSumFitAlgorithm(n=1, withconst=true))
fit_01_sim = solve(CurveFitProblem(xdata, ydata_01_sim), ExpSumFitAlgorithm(n=1, withconst=true))
fit_05_sim = solve(CurveFitProblem(xdata, ydata_05_sim), ExpSumFitAlgorithm(n=1, withconst=true))

# Calculate time scales (tau) from fit parameters
tau_exp_ctl = inv(-fit_ctl.u.λ[])
tau_exp_50 = inv(-fit_50.u.λ[])
tau_exp_200 = inv(-fit_200.u.λ[])
tau_sim_ctl = inv(-fit_ctl_sim.u.λ[])
tau_sim_01 = inv(-fit_01_sim.u.λ[])
tau_sim_05 = inv(-fit_05_sim.u.λ[])

println("The time scales for experiments: ")
for (tau, freq) in zip((tau_exp_ctl, tau_exp_50, tau_exp_200), (0, 50, 200))
    println("$freq uM ROS is $(round(tau; digits=2)) seconds.")
end

println("\nThe time scales for simulations: ")
for (tau, freq) in zip((tau_sim_ctl, tau_sim_01, tau_sim_05), (0, 50, 200))
    println("$freq uM ROS is $(round(tau; digits=2)) seconds.")
end
