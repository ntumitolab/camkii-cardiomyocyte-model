# # Effects of isoproterenol
using Model
using Model: second, μM
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CSV
using DataFrames
import Dates
Plots.default(lw=1.5)

# ## Setup model
@time "Build system" sys = Model.DEFAULT_SYS
tend = 205second
@time "Build problem" prob = ODEProblem(sys, [], tend)
stimstart = 30second
stimend = 120second
alg = KenCarp47()

# ## Without isoproterenol
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, alg; callback)

#---
i = (sys.t / 1000, sys.vm)
plot(sol, idxs=i, tspan=(100second, 101second), title="Action potential", xlabel="Time (s)")

#---
plot(sol, idxs=(sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(100second, 101second), title="Calcium transient", xlabel="Time (s)", ylabel="Conc. (μM)", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"])

#---
plot(sol, idxs=(sys.t / 1000, sys.CaMKAct), title="Active CaMKII", label=false, ylabel="Active CaMKII fraction", xlabel="Time (s)")

# ## 0.1uM isoproterenol
prob2 = remake(prob, p=[sys.ISO => 0.1μM])
@time sol2 = solve(prob2, alg; callback)

#---
plot(sol2, idxs=(sys.t / 1000, sys.vm), tspan=(100second, 101second), title="Action potential", xlabel="Time (s)")

#---
plot(sol2, idxs=(sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(100second, 101second), title="Calcium transcient", xlabel="Time (s)", ylabel="Conc. (μM)", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"])

# ## Comparison
i = (sys.t / 1000, sys.Cai_mean)
tspan = (100second, 101second)
plot(sol, idxs=i, title="Calcium transcient", lab="ISO (-)"; tspan)
plot!(sol2, idxs=i, lab="ISO (0.1uM)", xlabel="Time (s)", ylabel="Concentration (μM)"; tspan)

#---
savefig("iso-caT.pdf")
savefig("iso-caT.png")
# Maximal and minimal calcium concentrations in calcium transients.
ca_ctl = sol(tspan[1]:1:tspan[2], idxs=sys.Cai_mean)
println(extrema(ca_ctl))

#---
ca_iso = sol2(tspan[1]:1:tspan[2], idxs=sys.Cai_mean)
println(extrema(ca_iso))

#---
i = (sys.t / 1000, sys.CaMKAct)
plot(sol, idxs=i, title="", lab="ISO (-)")
plot!(sol2, idxs=i, lab="ISO (0.1uM)", ylabel="Active CaMKII fraction", xlabel="Time (s)")

#---
savefig("iso-camkact.pdf")
savefig("iso-camkact.png")

#---
i = (sys.t / 1000, sys.vm)
tspan = (100second, 101second)
plot(sol, idxs=i, title="Action potential", lab="ISO (-)"; tspan)
plot!(sol2, idxs=i, lab="ISO (0.1uM)", xlabel="Time (s)", ylabel="Voltage (mV)"; tspan)
#---
savefig("iso-ap.pdf")
savefig("iso-ap.png")

# ## Experimental data
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical.csv"), DataFrame)
ts = Dates.value.(chemicaldf[!, "Time"]) ./ 10^9
ctl = chemicaldf[!, "Ctrl Mean"]
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"])

#---
iso = chemicaldf[!, "isoproterenol 100nM Mean"]
iso_error = chemicaldf[!, "isoproterenol 100nM SD"] ./ sqrt.(chemicaldf[!, "isoproterenol 100nM N"])

#---
plot(ts, ctl, yerr=ctl_error, lab="Control", color=:blue, markerstrokecolor=:blue)
plot!(ts, iso, yerr=iso_error, lab="ISO 100nM", color=:red, markerstrokecolor=:red)
plot!(xlabel="Time (s)", ylabel="CaMKII activity (A.U.)")

#---
savefig("iso-exp.pdf")
savefig("iso-exp.png")
