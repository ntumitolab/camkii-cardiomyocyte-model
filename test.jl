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
alg = FBDF()
@time sol = solve(prob, FBDF(); callback, progress=true, abstol=1e-8, reltol=1e-8)

plot(sol, idxs=[sys.CaJSR, sys.CaNSR])
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL])
# RyR instability?
plot(sol, idxs=[sys.PO1RyR])
plot(sol, idxs=[sys.dPO1RyR])
plot(sol, idxs=[sys.Jup], tspan=(0.08, 0.10))
