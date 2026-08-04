#===
# Mutant CaMKII

Increased sensitivty to CaM-Ca binding by increasing the forward binding rate 100%, mimicking R275H mutation.

Fitting results:

- Basal activity: 4.2598803909827387e-19
- Maximal activity by CaM-Ca2 binding: 0.31337501524655215
- Half saturation Ca concentration for CaM-Ca2 binding: 0.4732164328053309 μM
- Maximal activity by CaM-Ca4 binding: 0.11553337613420915
- Half saturation Ca concentration for CaM-Ca4 binding: 0.37215724954233154 μM
- RMSE: 0.0013323008899588052
===#
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

#---
@time "Build system" sys = Model.DEFAULT_SYS
tend = 205.0second
@time "Build problem" prob = ODEProblem(sys, [], tend)
@time "Build mutant problem" prob_mut = remake(prob, p=[sys.kfa2_CaMK => 0.3134, sys.kmCa2_CaMK => 0.4732μM, sys.kfa4_CaMK => 0.1155, sys.kmCa4_CaMK => 0.3722μM, sys.kfb_CaMK => 0.0])

@unpack Istim = sys
stimstart = 30.0second
stimend = 120.0second
alg = KenCarp47()

callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
callback2 = build_stim_callbacks(Istim, stimend; period=0.5second, starttime=stimstart)
@time "Solve problem" sol1 = solve(prob, alg; callback)
@time "Solve problem" sol2 = solve(prob, alg; callback=callback2)
@time "Solve mutant problem" sol_mut = solve(prob_mut, alg; callback)
@time "Solve mutant problem" sol2_mut = solve(prob_mut, alg; callback=callback2)

iact = (sys.t / 1000, sys.CaMKAct)
iphos = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2))
fig10a = plot(sol1, idxs=iact, lab="WT (1 Hz)", color=:blue)
plot!(fig10a, sol1, idxs=iphos, lab=false, color=:blue, linestyle=:dot)
plot!(fig10a, sol2, idxs=iact, lab="WT (2 Hz)", color=:red)
plot!(fig10a, sol2, idxs=iphos, lab=false, color=:red, linestyle=:dot)
plot!(fig10a, sol_mut, idxs=iact, lab="Mutant (1 Hz)", color=:green)
plot!(fig10a, sol_mut, idxs=iphos, lab=false, color=:green, linestyle=:dot)
plot!(fig10a, sol2_mut, idxs=iact, lab="Mutant (2 Hz)", color=:orange)
plot!(fig10a, sol2_mut, idxs=iphos, lab=false, color=:orange, linestyle=:dot)
plot!(fig10a, title="A", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left, ylims=(0, 1))

#---
icat = (sys.t / 1000 - 100, sys.Cai_mean * 1000)
tspan = (100.0second, 101.0second)
fig10b = plot(sol1, idxs=icat, lab="WT (1 Hz)", color=:blue, tspan=tspan)
plot!(fig10b, sol2, idxs=icat, lab="WT (2 Hz)", color=:red, tspan=tspan)
plot!(fig10b, sol_mut, idxs=icat, lab="Mutant (1 Hz)", color=:green, tspan=tspan)
plot!(fig10b, sol2_mut, idxs=icat, lab="Mutant (2 Hz)", color=:orange, tspan=tspan)
plot!(fig10b, title="B", xlabel="Time (s)", ylabel="Cytosolic Ca concentration (nM)", titlelocation = :left)

#---
prob_iso = remake(prob, p=[sys.ISO => 0.1μM, sys.ICa_scale_ISO => 2.5])
prob_iso_mut = remake(prob_mut, p=[sys.ISO => 0.1μM, sys.ICa_scale_ISO => 2.5])
@time "Solve problem" sol_iso = solve(prob_iso, alg; callback)
@time "Solve mutant problem" sol_iso_mut = solve(prob_iso_mut, alg; callback)

fig10c = plot(sol1, idxs=iact, lab="WT (w/o ISO)", color=:blue)
plot!(fig10c, sol1, idxs=iphos, lab=false, color=:blue, linestyle=:dot)
plot!(fig10c, sol_iso, idxs=iact, lab="WT (ISO)", color=:red)
plot!(fig10c, sol_iso, idxs=iphos, lab=false, color=:red, linestyle=:dot)
plot!(fig10c, sol_mut, idxs=iact, lab="Mutant (w/o ISO)", color=:green)
plot!(fig10c, sol_mut, idxs=iphos, lab=false, color=:green, linestyle=:dot)
plot!(fig10c, sol_iso_mut, idxs=iact, lab="Mutant (ISO)", color=:orange)
plot!(fig10c, sol_iso_mut, idxs=iphos, lab=false, color=:orange, linestyle=:dot)
plot!(fig10c, title="C", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left, ylims=(0, 1))

#---
icat = (sys.t / 1000 - 100, sys.Cai_mean * 1000)
tspan = (100.0second, 101.0second)
fig10d = plot(sol1, idxs=icat, lab="WT (w/o ISO)", color=:blue, tspan=tspan)
plot!(fig10d, sol_iso, idxs=icat, lab="WT (ISO)", color=:red, tspan=tspan)
plot!(fig10d, sol_mut, idxs=icat, lab="Mutant (w/o ISO)", color=:green, tspan=tspan)
plot!(fig10d, sol_iso_mut, idxs=icat, lab="Mutant (ISO)", color=:orange, tspan=tspan)
plot!(fig10d, title="D", xlabel="Time (s)", ylabel="Cytosolic Ca concentration (nM)", titlelocation = :left)

#---
prob_ros = remake(prob, p=[sys.ROS => 50μM])
prob_ros_mut = remake(prob_mut, p=[sys.ROS => 50μM])
@time "Solve problem" sol_ros = solve(prob_ros, alg; callback)
@time "Solve mutant problem" sol_ros_mut = solve(prob_ros_mut, alg; callback)
fig10e = plot(sol1, idxs=iact, lab="WT (w/o ROS)", color=:blue)
plot!(fig10e, sol_ros, idxs=iact, lab="WT (ROS)", color=:red)
plot!(fig10e, sol_mut, idxs=iact, lab="Mutant (w/o ROS)", color=:green)
plot!(fig10e, sol_ros_mut, idxs=iact, lab="Mutant (ROS)", color=:orange)
plot!(fig10e, title="E", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left, ylims=(0, 1))

#---
fig10f = plot(sol1, idxs=icat, lab="WT (w/o ROS)", color=:blue, tspan=tspan)
plot!(fig10f, sol_ros, idxs=icat, lab="WT (ROS)", color=:red, tspan=tspan)
plot!(fig10f, sol_mut, idxs=icat, lab="Mutant (w/o ROS)", color=:green, tspan=tspan)
plot!(fig10f, sol_ros_mut, idxs=icat, lab="Mutant (ROS)", color=:orange, tspan=tspan)
plot!(fig10f, title="F", xlabel="Time (s)", ylabel="Cytosolic Ca concentration (nM)", titlelocation = :left)

#---
fig10 = plot(fig10a, fig10b, fig10c, fig10d, fig10e, fig10f, layout=(3, 2), size=(1000, 1000))

savefig(fig10, "fig10.png")
savefig(fig10, "fig10.pdf")

# ## Decay time constants
# Fit against an exponential decay model
# Data from experiments: record 50 seconds after pacing ends
ts = collect(range(0.0, stop=50.0, step=5.0))
ysim_1hz = sol1(stimend:5second:stimend+50second; idxs=sys.CaMKAct).u
ysim_1hz_mut = sol_mut(stimend:5second:stimend+50second; idxs=sys.CaMKAct).u

fit_1hz_sim = solve(CurveFitProblem(ts, ysim_1hz), ExpSumFitAlgorithm(n=1, withconst=true))
fit_1hz_mut_sim = solve(CurveFitProblem(ts, ysim_1hz_mut), ExpSumFitAlgorithm(n=1, withconst=true))

# Decay time scales (tau) from fit parameters
tau_sim_1hz = inv(-fit_1hz_sim.u.λ[])
tau_sim_1hz_mut = inv(-fit_1hz_mut_sim.u.λ[])

println("Decay time scale (tau) for WT (1 Hz): ", tau_sim_1hz, " s")
println("Decay time scale (tau) for Mutant (1 Hz): ", tau_sim_1hz_mut, " s")
