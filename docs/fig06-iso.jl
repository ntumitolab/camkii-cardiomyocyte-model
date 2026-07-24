# # Fig6: Effects of isoproterenol
using Model
using Model: second, μM
using CSV
using CurveFit
using DataFrames
using DiffEqCallbacks
using ModelingToolkit
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using Plots
using SteadyStateDiffEq
Plots.default(lw=1.5)

# ## Experimental data
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical-normalized.csv"), DataFrame)
ts = 0:5:205
ctl = chemicaldf[!, "Ctrl Mean"]
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"])
iso = chemicaldf[!, "isoproterenol 100nM Mean"]
iso_error = chemicaldf[!, "isoproterenol 100nM SD"] ./ sqrt.(chemicaldf[!, "isoproterenol 100nM N"])

fig6a = plot(ts, ctl, yerr=ctl_error, lab="Control", color=:blue, markerstrokecolor=:blue)
plot!(fig6a, ts, iso, yerr=iso_error, lab="ISO 100nM", color=:red, markerstrokecolor=:red)
plot!(fig6a, xlabel="Time (s)", ylabel="CaMKAR (R/R0)", title="A", titlelocation=:left)

# ## Simulation
@time "Build system" sys = Model.DEFAULT_SYS
tend = 205second
@time "Build problem" prob = ODEProblem(sys, [], tend)
## Bump up the ICaL scaling factor to simulate ISO effect (previously 1.56x, now 2.5x)
prob2 = remake(prob, p=[sys.ISO => 0.1μM, sys.ICa_scale_ISO => 2.5])
stimstart = 30second
stimend = 120second
alg = KenCarp47()

@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, alg; callback)
@time sol2 = solve(prob2, alg; callback)

# ### CaMKII activity
idxs = (sys.t / 1000, sys.CaMKAct)
fig6b = plot(sol, idxs=idxs, lab="Control", color=:blue, title="B", titlelocation=:left)
plot!(fig6b, sol2, idxs=idxs, lab="ISO 100 nM", color=:red, xlabel="Time (s)", ylabel="Active CaMKII fraction")

# ### Action potential
idxs = (sys.t / 1000, sys.vm)
tspan = (100second, 101second)
fig6c = plot(sol, idxs=idxs, lab="Control", color=:blue, title="C", titlelocation=:left, tspan=tspan)
plot!(fig6c, sol2, idxs=idxs, lab="ISO 100 nM", color=:red, xlabel="Time (s)", ylabel="Membrane potential (mV)", tspan=tspan)

# ### Calcium transient
idxs = (sys.t / 1000, sys.Cai_mean * 1000)
fig6d = plot(sol, idxs=idxs, lab="Control", color=:blue, title="D", titlelocation=:left, tspan=tspan)
plot!(fig6d, sol2, idxs=idxs, lab="ISO 100 nM", color=:red, xlabel="Time (s)", ylabel="Mean intracellular Ca (nM)", tspan=tspan)

#---
ca_ctl = sol(tspan[1]:1:tspan[2], idxs=sys.Cai_mean)
println(extrema(ca_ctl))
ca_iso = sol2(tspan[1]:1:tspan[2], idxs=sys.Cai_mean)
println(extrema(ca_iso))

# ## Combined figure
fig6 = plot(fig6a, fig6b, fig6c, fig6d, layout=(2, 2), size=(800, 600))

#---
savefig(fig6, "fig6.png")
savefig(fig6, "fig6.pdf")
