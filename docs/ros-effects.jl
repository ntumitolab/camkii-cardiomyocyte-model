# # ROS effects
using ModelingToolkit
using OrdinaryDiffEq, DiffEqCallbacks
using Plots
using CSV
using DataFrames
import Dates
using CurveFit
using CaMKIIModel
using CaMKIIModel: second, μM
Plots.default(lw=1.5)

# ## Setup model
@time @mtkcompile sys = build_neonatal_ecc_sys()
tend = 205second
@time prob = ODEProblem(sys, [], tend)
stimstart = 30second
stimend = 120second
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
alg = FBDF()

# ## Comparisons
@time sol = solve(prob, alg; callback)

## ROS (H2O2) 0.1uM
prob2 = remake(prob, p=[sys.ROS => 0.1μM])
@time sol2 = solve(prob2, alg; callback)

## ROS (H2O2) 0.5uM
prob3 = remake(prob, p=[sys.ROS => 0.5μM])
@time sol3 = solve(prob3, alg; callback);

# ## Comparisons
i = (sys.t / 1000, sys.CaMKAct * 100)
plot(sol, idxs=i, lab="ROS (-)")
plot!(sol2, idxs=i, lab="ROS 0.1uM")
plot!(sol3, idxs=i, lab="ROS 0.5uM")
plot!(xlabel="Time (s)", ylabel="Active fraction (%)", title="Simulation")

#---
savefig("ros-camkii.pdf")

# Oxidized fraction
i = (sys.t / 1000, 100 * (sys.CaMKBOX + sys.CaMKPOX + sys.CaMKAOX + sys.CaMKOX ))
plot(sol, idxs=i, lab="ROS (-)")
plot!(sol2, idxs=i, lab="ROS 0.1uM")
plot!(sol3, idxs=i, lab="ROS 0.5uM")
plot!(xlabel="Time (s)", ylabel="Oxidized fraction (%)", title="Simulation")

#---
savefig("ros-camkiiox.pdf")

# Autophosphorylated fraction
i = (sys.t / 1000, 100 * (sys.CaMKP + sys.CaMKA + sys.CaMKA2))
plot(sol, idxs=i, lab="ROS (-)")
plot!(sol2, idxs=i, lab="ROS 0.1uM")
plot!(sol3, idxs=i, lab="ROS 0.5uM")
plot!(xlabel="Time (s)", ylabel="Phosphorylated fraction (%)", title="Simulation")

#---
savefig("ros-camkiip.pdf")

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

# ## Decay rates
# Fit against an exponential decay model.
decay_model(p, x) = @. p[1] * exp(-x / p[2]) + p[3]

# Data from experiments (starts from 120 seconds till the end of the experiment).
xdata = collect(range(0.0, step=5.0, length=18))
ydata_ctl = ctl[25:end]
ydata_50 = ros50[25:end]
ydata_200 = ros200[25:end]

# Fit data to an exponential decay model
fit_ctl = solve(CurveFitProblem(xdata, ydata_ctl), ExpSumFitAlgorithm(n=1, withconst=true))
fit_50 = solve(CurveFitProblem(xdata, ydata_50), ExpSumFitAlgorithm(n=1, withconst=true))
fit_200 = solve(CurveFitProblem(xdata, ydata_200), ExpSumFitAlgorithm(n=1, withconst=true))

# Calculate time scales (tau) from fit parameters
tau_exp_ctl = inv(-fit_ctl.u.λ[])
tau_exp_50 = inv(-fit_50.u.λ[])
tau_exp_200 = inv(-fit_200.u.λ[])

println("The time scales for experiments: ")
for (tau, freq) in zip((tau_exp_ctl, tau_exp_50, tau_exp_200), (0, 50, 200))
    println("$freq uM ROS is $(round(tau; digits=2)) seconds.")
end
