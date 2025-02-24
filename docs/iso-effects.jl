# # Effects of isoproterenol
using ModelingToolkit
using DifferentialEquations
using Plots
using CSV
using DataFrames
using Dates
using CaMKIIModel
using CaMKIIModel: second, μM
Plots.default(lw=1.5)

# ## Setup model
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)
tend = 205second
prob = ODEProblem(sys, [], tend)
stimstart = 30second
stimend = 120second
alg = TRBDF2()

# ## Without isoproterenol
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, alg; callback)

#---
i = (sys.t/1000, sys.vm)
plot(sol, idxs=i, tspan=(100second, 101second), title="Action potential", xlabel="Time (s)")

#---
plot(sol, idxs=(sys.t/1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(100second, 101second), title="Calcium transient", xlabel="Time (s)", ylabel="Conc. (μM)", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"])

#---
plot(sol, idxs=(sys.t/1000, sys.CaMKAct*100), title="Active CaMKII", label=false, ylabel="Active fraction (%)" , xlabel="Time (s)")

# ## 0.1uM isoproterenol
prob2 = remake(prob, p=[sys.ISO => 0.1μM])
sol2 = solve(prob2, alg; callback)

#---
plot(sol2, idxs=(sys.t/1000, sys.vm), tspan=(100second, 101second), title="Action potential", xlabel="Time (s)")

#---
plot(sol2, idxs=(sys.t/1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean]), tspan=(100second, 101second), title="Calcium transcient", xlabel="Time (s)", ylabel="Conc. (μM)", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"])

#---
plot(sol2, idxs=(sys.t/1000, sys.CaMKAct*100), title="Active CaMKII", label=false, ylabel="Active fraction (%)" , xlabel="Time (s)")

# ## Comparison
i = (sys.t/1000, sys.Cai_mean)
tspan = (100second, 101second)
plot(sol, idxs=i, title="Calcium transcient", lab="ISO (-)"; tspan)
plot!(sol2, idxs=i, lab="ISO (0.1uM)", xlabel="Time (s)", ylabel="Concentration (μM)"; tspan)

#---
savefig("iso-caT.pdf")

#---
i = (sys.t/1000, sys.CaMKAct*100)
plot(sol, idxs=i, title="Active CaMKII", lab="ISO (-)")
plot!(sol2, idxs=i, lab="ISO (0.1uM)", ylabel="Active fraction (%)" , xlabel="Time (s)")

#---
savefig("iso-camkact.pdf")

#---
i = (sys.t/1000, sys.vm)
tspan = (100second, 101second)
plot(sol, idxs=i, title="Action potential", lab="ISO (-)"; tspan)
plot!(sol2, idxs=i, lab="ISO (0.1uM)", xlabel="Time (ms)", ylabel="Voltage (mV)"; tspan)

# ## Experimental data
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical.csv"), DataFrame)
ts = Dates.value.(chemicaldf[!, "Time"]) ./ 10^9
ctl = chemicaldf[!, "Ctrl Mean"]
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"])

iso = chemicaldf[!, "isoproterenol 100nM Mean"]
iso_error = chemicaldf[!, "isoproterenol 100nM SD"] ./ sqrt.(chemicaldf[!, "isoproterenol 100nM N"])

plot(ts, ctl, yerr=ctl_error, lab="Control", color=:blue, markerstrokecolor=:blue)
plot!(ts, iso, yerr=iso_error, lab="ISO 100nM", color=:red, markerstrokecolor=:red)
plot!(xlabel="Time (sec.)", ylabel="CaMKII activity (A.U.)")

#---
savefig("iso-exp.pdf")
