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
@time sol = solve(prob, TRBDF2(); callback, progress=true, abstol=1e-7, reltol=1e-7)

plot(sol, idxs=[sys.vm*1000])

plot(sol, idxs=[sys.CaJSR, sys.CaNSR, sys.Cai_sub_SR])

plot(sol, idxs=[sys.PO1RyR])
@unpack Cai_sub_SR = sys
