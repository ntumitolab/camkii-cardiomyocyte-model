using ProgressLogging
using ModelingToolkit
using DifferentialEquations
using Plots
using CaMKIIModel
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true)
tend = 1000.0
prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend)
alg = FBDF()
@time sol = solve(prob, alg; callback, progress=true, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))


plot(sol, idxs=sys.vm, tspan=(90, 100))
plot(sol, idxs=[sys.CaJSR, sys.CaNSR], tspan=(90, 100))
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(99, 100))
plot(sol, idxs=sys.Jtr)

u0 = sol[end]
@time sol = solve(remake(prob, u0=u0), alg; progress=true, abstol=1e-9, reltol=1e-9)
