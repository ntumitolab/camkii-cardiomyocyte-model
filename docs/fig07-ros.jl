# # Fig7: ROS effects
using CSV
using CurveFit
using DataFrames
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEqSDIRK
using Plots
using StatsPlots
using Model
using Model: second, μM
Plots.default(lw=1.5)

# ## Experimental data
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical.csv"), DataFrame)
ts = 0:5:205
ctl = chemicaldf[!, "Ctrl Mean"]
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"])

#---
ros50 = chemicaldf[!, "H2O2 50uM Mean"]
ros50_error = chemicaldf[!, "H2O2 50uM SD"] ./ sqrt.(chemicaldf[!, "H2O2 50uM N"])
ros200 = chemicaldf[!, "H2O2 200uM Mean"]
ros200_error = chemicaldf[!, "H2O2 200uM SD"] ./ sqrt.(chemicaldf[!, "H2O2 200uM N"])

#---
fig7a = plot(ts, ctl, yerr=ctl_error, lab="Control", color=:blue, markerstrokecolor=:blue)
plot!(fig7a, ts, ros50, yerr=ros50_error, lab="50μM H2O2", color=:red, markerstrokecolor=:red)
plot!(fig7a, ts, ros200, yerr=ros200_error, lab="200μM H2O2", color=:green, markerstrokecolor=:green)
plot!(fig7a, xlabel="Time (s)", ylabel="CaMKAR measurement (AU)", title="A", titlelocation=:left)

# ## Simulation
@time "Build system" sys = Model.DEFAULT_SYS
@unpack Istim = sys
tend = 205second
stimstart = 30second
stimend = 120second
@time "Build problem" prob = ODEProblem(sys, [], tend)
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
alg = FBDF()

prob2 = remake(prob, p=[sys.ROS => 50μM])
prob3 = remake(prob, p=[sys.ROS => 200μM])
@time "Solve problem" sol = solve(prob, alg; callback)
@time "Solve problem" sol2 = solve(prob2, alg; callback)
@time "Solve problem" sol3 = solve(prob3, alg; callback)

idxs = (sys.t / 1000, sys.CaMKAct)
fig7b = plot(sol, idxs=idxs, lab="Control", color=:blue)
plot!(fig7b, sol2, idxs=idxs, lab="50μM H2O2", color=:red)
plot!(fig7b, sol3, idxs=idxs, lab="200μM H2O2", color=:green)
plot!(fig7b, title="B", titlelocation=:left, xlabel="Time (s)", ylabel="Active CaMKII fraction")

# ### Oxidized and autophosphorylated fraction
iox = (sys.t / 1000, (sys.CaMKBOX + sys.CaMKPOX + sys.CaMKAOX + sys.CaMKOX))
fig7c = plot(sol, idxs=iox, lab="Control, OX", color=:blue)
plot!(fig7c, sol2, idxs=iox, lab="50μM H2O2, OX", color=:red)
plot!(fig7c, sol3, idxs=iox, lab="200μM H2O2, OX", color=:green)
iphos = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2))
plot!(fig7c, sol, idxs=iphos, lab="Control, phos", color=:blue, linestyle=:dash)
plot!(fig7c, sol2, idxs=iphos, lab="50μM H2O2, phos", color=:red, linestyle=:dash)
plot!(fig7c, sol3, idxs=iphos, lab="200μM H2O2, phos", color=:green, linestyle=:dash)
plot!(fig7c, xlabel="Time (s)", ylabel="Active CaMKII fraction", title="C", titlelocation=:left, legend=:topright)

# Proportions for oxidized and autophosphorylated fractions in active CaMKII
# Using StatsPlots for grouped bar plot at the end of the stimulation period (120 seconds)
timepoint = 120second
pidx = (sys.CaMKP + sys.CaMKA + sys.CaMKA2) / sys.CaMKAct
oidx = (sys.CaMKBOX + sys.CaMKPOX + sys.CaMKAOX + sys.CaMKOX) / sys.CaMKAct
pvals = [s(timepoint, idxs=pidx) for s in (sol, sol2, sol3)]
ovals = [s(timepoint, idxs=oidx) for s in (sol, sol2, sol3)]
fig7d = groupedbar(["0uM", "50uM", "200uM"], [pvals ovals], bar_position=:stack, label=["Phosphorylated" "Oxidized"], color=[:blue :red :green], title="D", ylabel="Fraction of active CaMKII", xlabel="H2O2 concentration", legend=:topleft, titlelocation=:left)

# ## Save figure
plot(fig7a, fig7b, fig7c, fig7d, layout=(2, 2), size=(1000, 800))
savefig("fig07.png")
savefig("fig07.pdf")

# ## Decay kinectics
# Fit against an exponential decay model.

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
    println("$freq uM H2O2 is $(round(tau; digits=2)) seconds.")
end

println("\nThe time scales for simulations: ")
for (tau, freq) in zip((tau_sim_ctl, tau_sim_01, tau_sim_05), (0, 50, 200))
    println("$freq uM H2O2 is $(round(tau; digits=2)) seconds.")
end
