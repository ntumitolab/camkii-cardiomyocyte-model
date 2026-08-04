#===
# Mutant CaMKII

Increased sensitivty to CaM-Ca binding by increasing the forward binding rate 100%.

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
@time "Solve problem" sol1 = solve(prob, alg; callback)
@time "Solve mutant problem" sol_mut = solve(prob_mut, alg; callback)

idxs = (sys.t / 1000, sys.CaMKAct)
fig10a = plot(sol1, idxs=idxs, lab="WT", color=:black)
plot!(fig10a, sol_mut, idxs=idxs, lab="Mutant", color=:red)
