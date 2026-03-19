# # Calcium transient fitting
using ModelingToolkit
using CurveFit
using OrdinaryDiffEq
using Plots
using CaMKIIModel
using CaMKIIModel: second, μM

# ## Setup the ODE system
stimstart = 0.0second
stimend = 100.0second
tend = stimend
@time "Building ODE system" @mtkcompile sys = build_neonatal_ecc_sys()
@time "Building ODE problem" prob = ODEProblem(sys, [], tend)
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time "Solving ODE problem" sol = solve(prob, KenCarp47(); callback)

#---
cai = sol(range(99second, stop=100second, length=101), idxs=sys.Cai_mean).u
ts = collect(0:0.01:1)
plot(ts, cai)

# Fit
prob = CurveFitProblem(ts, cai)
@time sol = solve(prob, RationalPolynomialFitAlgorithm(4, 4))


plot(ts,[cai fitted(sol)], label=["Data" "Fit"], xlabel="Time (s)", ylabel="Ca (μM)", title="Calcium transient fitting")
