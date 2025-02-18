# # ROS effects
using ModelingToolkit
using OrdinaryDiffEq
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
@unpack Istim = sys
alg = FBDF()

# ## No ROS
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, alg; callback)

#---
i = (sys.t / 1000, sys.vm)
tspan=(100second, 101second)
plot(sol, idxs=i, title="Action potential"; tspan)

#---
i = (sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean])
tspan=(100second, 101second)
plot(sol, idxs=i, title="Calcium transcient", label=["Ca (Sub SR)" "Ca (Sub SL)" "Ca (avg)"]; tspan)

#---
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol, idxs=i, title="CaMKII", xlabel="Time (s)", ylabel="Active fraction (%)", label=false)

# ## ROS 0.1uM
prob2 = remake(prob, p=[sys.ROS => 0.1μM])
@time sol2 = solve(prob2, alg; callback)

#---
i = (sys.t / 1000, sys.vm)
tspan=(100second, 101second)
plot(sol2, idxs=i, title="Action potential"; tspan)

#---
i = (sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean])
tspan=(100second, 101second)
plot(sol2, idxs=i, title="Calcium transcient", label=["Ca (Sub SR)" "Ca (Sub SL)" "Ca (avg)"]; tspan)

#---
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol2, idxs=i, title="CaMKII", xlabel="Time (s)", ylabel="Active fraction (%)", label=false)

# ## ROS 1uM
prob3 = remake(prob, p=[sys.ROS => 1μM])
@time sol3 = solve(prob3, alg; callback)

#---
i = (sys.t / 1000, sys.vm)
tspan=(100second, 101second)
plot(sol3, idxs=i, title="Action potential"; tspan)

#---
i = (sys.t / 1000, [sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean])
tspan=(100second, 101second)
plot(sol3, idxs=i, title="Calcium transcient", label=["Ca (SR)" "Ca (SL)" "Ca (avg)"]; tspan)

#---
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol3, idxs=i, title="CaMKII", xlabel="Time (s)", ylabel="Active fraction (%)", label=false)

# ## Comparisons
i = (sys.t / 1000,sys.CaMKAct*100)
plot(sol, idxs=i, title="Active CaMKII", lab="ROS (-)")
plot!(sol2, idxs=i, lab="ROS 0.1uM")
plot!(sol3, idxs=i, lab="ROS 1uM", xlabel="Time (s)", ylabel="Active fraction (%)")

#---
savefig("ros-camkii.pdf")

# ## Experimental data
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical.csv"), DataFrame)
ts = Dates.value.(chemicaldf[!, "Time"]) ./ 10^9
ctl = chemicaldf[!, "Ctrl Mean"]
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"])

ros50 = chemicaldf[!, "H2O2 50uM Mean"]
ros50_error = chemicaldf[!, "H2O2 50uM SD"] ./ sqrt.(chemicaldf[!, "H2O2 50uM N"])
ros200 = chemicaldf[!, "H2O2 200uM Mean"]
ros200_error = chemicaldf[!, "H2O2 200uM SD"] ./ sqrt.(chemicaldf[!, "H2O2 200uM N"])

plot(ts, ctl, yerr=ctl_error, lab="Control", color=:blue, markerstrokecolor=:blue)
plot!(ts, ros50, yerr=ros50_error, lab="H2O2 50uM", color=:red, markerstrokecolor=:red)
plot!(ts, ros200, yerr=ros200_error, lab="H2O2 200uM", color=:green, markerstrokecolor=:green)
plot!(xlabel="Time (s)", ylabel="CaMKII activity (A.U.)")

#---
savefig("ros-exp.pdf")
