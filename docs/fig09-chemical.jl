# # Chemical Effects
# Validating and predicting the model
# AS 100397: CaMKII competitive inhibitor at the catalytic site
# CalA: PP1 inhibitor, prolonging phosphorylated CaMKII activity
using Model
using Model: second, Hz
using CSV
using DataFrames
using DiffEqCallbacks
using ModelingToolkit
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using Plots
using Plots.Measures
Plots.default(lw=1.5)

# ## Experimental data
# Data on CaMKII activity with chemical treatment
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical.csv"), DataFrame)
ts = 0:5:205
# Normalize the mean and SD relative to the initial value
ctl_initial = chemicaldf[!, "Ctrl Mean"][1]
as10093_initial = chemicaldf[!, "AS 10093 Mean"][1]
cala_initial = chemicaldf[!, "CalA Mean"][1]

ctl = chemicaldf[!, "Ctrl Mean"] ./ ctl_initial
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"]) ./ ctl_initial
as10093 = chemicaldf[!, "AS 10093 Mean"] ./ as10093_initial
as10093_error = chemicaldf[!, "AS 10093 SD"] ./ sqrt.(chemicaldf[!, "AS 10093 N"]) ./ as10093_initial
cala = chemicaldf[!, "CalA Mean"] ./ cala_initial
cala_error = chemicaldf[!, "CalA SD"] ./ sqrt.(chemicaldf[!, "CalA N"]) ./ cala_initial

fig9a = plot(ts, ctl, yerr=ctl_error, lab="Control", color=:blue, markerstrokecolor=:blue)
plot!(fig9a, ts, as10093, yerr=as10093_error, lab="AS 100397", color=:green, markerstrokecolor=:green)
plot!(fig9a, ts, cala, yerr=cala_error, lab="CalA", color=:red, markerstrokecolor=:red)
plot!(fig9a, xlabel="Time (s)", ylabel="CaMKAR (R/R0)", title= "A", titlelocation=:left)

# ## Simulation results
# AS 100397: suppressing CaMKII activation by 90%
# AS 100397 (0.5x): suppressing CaMKII activation by 50%
# CalA: prolonging CaMKII dephosphorylation time by 2-fold
# CalA (2x): prolonging CaMKII dephosphorylation time by 4-fold
@time "Build system" sys = Model.DEFAULT_SYS
tend = 205second
stimstart = 30second
stimend = 120second
@time "Build ctl problem" prob_ctl = ODEProblem(sys, [], tend)
@time "Build AS problem" prob_as = remake(prob_ctl, p=[sys.KActScale => 0.1])
@time "Build AS (0.5x) problem" prob_as05 = remake(prob_ctl, p=[sys.KActScale => 0.5])
@time "Build calA problem" prob_cala = remake(prob_ctl, p=[sys.kdeph_CaMK => inv(24second)])
@time "Build calA (++) problem" prob_cala2 = remake(prob_ctl, p=[sys.kdeph_CaMK => inv(100second)])

alg = KenCarp47()
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time "Solve problem" sol_ctl = solve(prob_ctl, alg; callback)
@time "Solve problem" sol_as = solve(prob_as, alg; callback)
@time "Solve problem" sol_as05 = solve(prob_as05, alg; callback)
@time "Solve problem" sol_cala = solve(prob_cala, alg; callback)
@time "Solve problem" sol_cala2 = solve(prob_cala2, alg; callback)

idxs = (sys.t / 1000, sys.CaMKAct)
fig9b = plot(sol_ctl, idxs=idxs, lab="Control", color=:blue)
plot!(fig9b, sol_as, idxs=idxs, lab="KAct 0.1x", color=:green)
plot!(fig9b, sol_cala, idxs=idxs, lab="CalA", color=:red)
plot!(fig9b, sol_as05, idxs=idxs, lab="KAct 0.5x", color=:cyan)
plot!(fig9b, sol_cala2, idxs=idxs, lab="CalA (++)", color=:orange)
plot!(fig9b, title="B", titlelocation=:left, xlabel="Time (s)", ylabel="Active CaMKII fraction")

# ## Save figure
plot(fig9a, fig9b, layout=(1, 2), size=(900, 400), bottom_margin = 5mm, left_margin = 8mm)
savefig("fig9.png")
savefig("fig9.pdf")
