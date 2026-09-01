#===
# Mutant CaMKII

Increased sensitivty to CaM-Ca binding by increasing the forward binding rate 100%, mimicking R275H mutation.

Fits:

- Basal activity: 4.2598803909827387e-19
- Maximal activity by CaM-Ca2 binding: 0.31337501524655215
- Half saturation Ca concentration for CaM-Ca2 binding: 0.4732164328053309 μM
- Maximal activity by CaM-Ca4 binding: 0.11553337613420915
- Half saturation Ca concentration for CaM-Ca4 binding: 0.37215724954233154 μM
- RMSE: 0.0013323008899588052
===#
using Model
using Model: second, μM, Hz
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

stim_freq = [0.5Hz, 1.0Hz, 1.5Hz, 2.0Hz]
stim_periods = inv.(stim_freq)
stim_durations = [15.0second, 30.0second, 60.0second, 90.0second, 120.0second]

function solve_prob(prob; stim_period=1.0second, stim_duration=90.0second, stimstart=30.0second, alg=KenCarp47(), add_iso=false, add_ros=false)
    stimend = stimstart + stim_duration
    @unpack Istim = sys
    callback = build_stim_callbacks(Istim, stimend; period=stim_period, starttime=stimstart)
    sol = solve(remake(prob; p=[sys.ISO => (add_iso ? 0.1μM : 0.0μM), sys.ROS => (add_ros ? 100μM : 0.0μM)]), alg; callback)
    ts = range(0.0, stop=50, step=5) .* second .+ stimend
    return sol(ts, idxs=sys.CaMKAct)
end

# ## Pacing frequency sweep
@time begin
    sol_15s_wt = [solve_prob(prob; stim_period=stim_period, stim_duration=15.0second).u[1] for stim_period in stim_periods]
    sol_30s_wt = [solve_prob(prob; stim_period=stim_period, stim_duration=30.0second).u[1] for stim_period in stim_periods]
    sol_60s_wt = [solve_prob(prob; stim_period=stim_period, stim_duration=60.0second).u[1] for stim_period in stim_periods]
    sol_90s_wt = [solve_prob(prob; stim_period=stim_period, stim_duration=90.0second).u[1] for stim_period in stim_periods]
    sol_120s_wt = [solve_prob(prob; stim_period=stim_period, stim_duration=120.0second).u[1] for stim_period in stim_periods]
    sol_15s_mut = [solve_prob(prob_mut; stim_period=stim_period, stim_duration=15.0second).u[1] for stim_period in stim_periods]
    sol_30s_mut = [solve_prob(prob_mut; stim_period=stim_period, stim_duration=30.0second).u[1] for stim_period in stim_periods]
    sol_60s_mut = [solve_prob(prob_mut; stim_period=stim_period, stim_duration=60.0second).u[1] for stim_period in stim_periods]
    sol_90s_mut = [solve_prob(prob_mut; stim_period=stim_period, stim_duration=90.0second).u[1] for stim_period in stim_periods]
    sol_120s_mut = [solve_prob(prob_mut; stim_period=stim_period, stim_duration=120.0second).u[1] for stim_period in stim_periods]
end

fig10a = plot(stim_freq .* 1000 , sol_15s_wt, lab="WT (15 s)", marker=:circle)
plot!(fig10a, stim_freq .* 1000 , sol_30s_wt, lab="WT (30 s)", marker=:circle)
plot!(fig10a, stim_freq .* 1000 , sol_60s_wt, lab="WT (60 s)", marker=:circle)
plot!(fig10a, stim_freq .* 1000 , sol_90s_wt, lab="WT (90 s)", marker=:circle)
plot!(fig10a, stim_freq .* 1000 , sol_120s_wt, lab="WT (120 s)", marker=:circle, )
plot!(fig10a, stim_freq .* 1000 , sol_15s_mut, lab="Mutant (15 s)", marker=:square)
plot!(fig10a, stim_freq .* 1000 , sol_30s_mut, lab="Mutant (30 s)", marker=:square)
plot!(fig10a, stim_freq .* 1000 , sol_60s_mut, lab="Mutant (60 s)", marker=:square)
plot!(fig10a, stim_freq .* 1000 , sol_90s_mut, lab="Mutant (90 s)", marker=:square)
plot!(fig10a, stim_freq .* 1000 , sol_120s_mut, lab="Mutant (120 s)", marker=:square)
plot!(fig10a, xlims=(0.4, 2.2), xlabel="Stimulus frequency (Hz)", ylabel="Max active CaMKII", title="A", titlelocation=:left, legend=:outerright)

# ## Pacing frequency sweep with ISO and ROS
@time begin
    sol_ctl_wt = [solve_prob(prob; stim_period=stim_period).u[1] for stim_period in stim_periods]
    sol_iso_wt = [solve_prob(prob; stim_period=stim_period, add_iso=true).u[1] for stim_period in stim_periods]
    sol_ros_wt = [solve_prob(prob; stim_period=stim_period, add_ros=true).u[1] for stim_period in stim_periods]
    sol_both_wt = [solve_prob(prob; stim_period=stim_period, add_iso=true, add_ros=true).u[1] for stim_period in stim_periods]
    sol_ctl_mut = [solve_prob(prob_mut; stim_period=stim_period).u[1] for stim_period in stim_periods]
    sol_iso_mut = [solve_prob(prob_mut; stim_period=stim_period, add_iso=true).u[1] for stim_period in stim_periods]
    sol_ros_mut = [solve_prob(prob_mut; stim_period=stim_period, add_ros=true).u[1] for stim_period in stim_periods]
    sol_both_mut = [solve_prob(prob_mut; stim_period=stim_period, add_iso=true, add_ros=true).u[1] for stim_period in stim_periods]
end

fig10b = plot(stim_freq .* 1000 , sol_ctl_wt, lab="WT (Ctl)", marker=:circle)
plot!(fig10b, stim_freq .* 1000 , sol_iso_wt, lab="WT (ISO)", marker=:circle)
plot!(fig10b, stim_freq .* 1000 , sol_ros_wt, lab="WT (ROS)", marker=:circle)
plot!(fig10b, stim_freq .* 1000 , sol_both_wt, lab="WT (Both)", marker=:circle)
plot!(fig10b, stim_freq .* 1000 , sol_ctl_mut, lab="Mutant (Ctl)", marker=:square)
plot!(fig10b, stim_freq .* 1000 , sol_iso_mut, lab="Mutant (ISO)", marker=:square)
plot!(fig10b, stim_freq .* 1000 , sol_ros_mut, lab="Mutant (ROS)", marker=:square)
plot!(fig10b, stim_freq .* 1000 , sol_both_mut, lab="Mutant (Both)", marker=:square)
plot!(fig10b, xlims=(0.4, 2.2), xlabel="Stimulus frequency (Hz)", ylabel="Max active CaMKII", title="B", titlelocation=:left, legend=:outerright)

#---
fig10 = plot(fig10a, fig10b, layout=(2, 1), size=(800, 800))
savefig(fig10, "fig10.png")
savefig(fig10, "fig10.pdf")

# ## Old code
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
callback2 = build_stim_callbacks(Istim, stimend; period=0.5second, starttime=stimstart)
@time "Solve problem" sol1 = solve(prob, alg; callback)
@time "Solve problem" sol2 = solve(prob, alg; callback=callback2)
@time "Solve mutant problem" sol_mut = solve(prob_mut, alg; callback)
@time "Solve mutant problem" sol2_mut = solve(prob_mut, alg; callback=callback2)

iact = (sys.t / 1000, sys.CaMKAct)
iphos = (sys.t / 1000, (sys.CaMKP + sys.CaMKA + sys.CaMKA2))
figS3a = plot(sol1, idxs=iact, lab="WT (1 Hz)", color=:blue)
plot!(figS3a, sol1, idxs=iphos, lab=false, color=:blue, linestyle=:dot)
plot!(figS3a, sol2, idxs=iact, lab="WT (2 Hz)", color=:red)
plot!(figS3a, sol2, idxs=iphos, lab=false, color=:red, linestyle=:dot)
plot!(figS3a, sol_mut, idxs=iact, lab="Mutant (1 Hz)", color=:green)
plot!(figS3a, sol_mut, idxs=iphos, lab=false, color=:green, linestyle=:dot)
plot!(figS3a, sol2_mut, idxs=iact, lab="Mutant (2 Hz)", color=:orange)
plot!(figS3a, sol2_mut, idxs=iphos, lab=false, color=:orange, linestyle=:dot)
plot!(figS3a, title="A", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left, ylims=(0, 1))

#---
icat = (sys.t / 1000 - 100, sys.Cai_mean * 1000)
tspan = (100.0second, 101.0second)
figS3b = plot(sol1, idxs=icat, lab="WT (1 Hz)", color=:blue, tspan=tspan)
plot!(figS3b, sol2, idxs=icat, lab="WT (2 Hz)", color=:red, tspan=tspan)
plot!(figS3b, sol_mut, idxs=icat, lab="Mutant (1 Hz)", color=:green, tspan=tspan)
plot!(figS3b, sol2_mut, idxs=icat, lab="Mutant (2 Hz)", color=:orange, tspan=tspan)
plot!(figS3b, title="B", xlabel="Time (s)", ylabel="Cytosolic Ca concentration (nM)", titlelocation = :left)

#---
prob_iso = remake(prob, p=[sys.ISO => 0.1μM, sys.ICa_scale_ISO => 2.5])
prob_iso_mut = remake(prob_mut, p=[sys.ISO => 0.1μM, sys.ICa_scale_ISO => 2.5])
@time "Solve problem" sol_iso = solve(prob_iso, alg; callback)
@time "Solve mutant problem" sol_iso_mut = solve(prob_iso_mut, alg; callback)

figS3c = plot(sol1, idxs=iact, lab="WT (w/o ISO)", color=:blue)
plot!(figS3c, sol1, idxs=iphos, lab=false, color=:blue, linestyle=:dot)
plot!(figS3c, sol_iso, idxs=iact, lab="WT (ISO)", color=:red)
plot!(figS3c, sol_iso, idxs=iphos, lab=false, color=:red, linestyle=:dot)
plot!(figS3c, sol_mut, idxs=iact, lab="Mutant (w/o ISO)", color=:green)
plot!(figS3c, sol_mut, idxs=iphos, lab=false, color=:green, linestyle=:dot)
plot!(figS3c, sol_iso_mut, idxs=iact, lab="Mutant (ISO)", color=:orange)
plot!(figS3c, sol_iso_mut, idxs=iphos, lab=false, color=:orange, linestyle=:dot)
plot!(figS3c, title="C", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left, ylims=(0, 1))

#---
icat = (sys.t / 1000 - 100, sys.Cai_mean * 1000)
tspan = (100.0second, 101.0second)
figS3d = plot(sol1, idxs=icat, lab="WT (w/o ISO)", color=:blue, tspan=tspan)
plot!(figS3d, sol_iso, idxs=icat, lab="WT (ISO)", color=:red, tspan=tspan)
plot!(figS3d, sol_mut, idxs=icat, lab="Mutant (w/o ISO)", color=:green, tspan=tspan)
plot!(figS3d, sol_iso_mut, idxs=icat, lab="Mutant (ISO)", color=:orange, tspan=tspan)
plot!(figS3d, title="D", xlabel="Time (s)", ylabel="Cytosolic Ca concentration (nM)", titlelocation = :left)

#---
prob_ros = remake(prob, p=[sys.ROS => 50μM])
prob_ros_mut = remake(prob_mut, p=[sys.ROS => 50μM])
@time "Solve problem" sol_ros = solve(prob_ros, alg; callback)
@time "Solve mutant problem" sol_ros_mut = solve(prob_ros_mut, alg; callback)
figS3e = plot(sol1, idxs=iact, lab="WT (w/o ROS)", color=:blue)
plot!(figS3e, sol_ros, idxs=iact, lab="WT (ROS)", color=:red)
plot!(figS3e, sol_mut, idxs=iact, lab="Mutant (w/o ROS)", color=:green)
plot!(figS3e, sol_ros_mut, idxs=iact, lab="Mutant (ROS)", color=:orange)
plot!(figS3e, title="E", xlabel="Time (s)", ylabel="Active CaMKII fraction", titlelocation = :left, ylims=(0, 1))

#---
figS3f = plot(sol1, idxs=icat, lab="WT (w/o ROS)", color=:blue, tspan=tspan)
plot!(figS3f, sol_ros, idxs=icat, lab="WT (ROS)", color=:red, tspan=tspan)
plot!(figS3f, sol_mut, idxs=icat, lab="Mutant (w/o ROS)", color=:green, tspan=tspan)
plot!(figS3f, sol_ros_mut, idxs=icat, lab="Mutant (ROS)", color=:orange, tspan=tspan)
plot!(figS3f, title="F", xlabel="Time (s)", ylabel="Cytosolic Ca concentration (nM)", titlelocation = :left)

#---
figS3 = plot(figS3a, figS3b, figS3c, figS3d, figS3e, figS3f, layout=(3, 2), size=(1000, 1000))

savefig(figS3, "figS3.png")
savefig(figS3, "figS3.pdf")

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
