using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using CaMKIIModel
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true)
tend = 30.0
prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend)
@time sol = solve(prob, Rodas5P(); callback, progress=true, abstol=1e-7, reltol=1e-7)

plot(sol, idxs=[sys.vm*1000])

plot(sol, idxs=[sys.CaJSR, sys.CaNSR])
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL])
# RyR instability?
plot(sol, idxs=[sys.betaSR])
