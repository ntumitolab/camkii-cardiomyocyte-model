using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true)
tend = 1000.0
prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend; period=0.5)
alg = FBDF()
@time sol = solve(prob, alg; callback, progress=true, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

plot(sol, idxs=sys.vm, tspan=(998, 1000))
plot(sol, idxs=[sys.CaJSR, sys.CaNSR], tspan=(998, 1000))
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(998, 1000))
plot(sol, idxs=sys.Jtr)
