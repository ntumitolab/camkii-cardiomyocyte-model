using ProgressLogging
using ModelingToolkit
using DifferentialEquations
using Plots
using CaMKIIModel
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(dx=4e-7, simplify=true)
tend = 30.0
prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend)
alg = FBDF()
@time sol = solve(prob, alg; callback, progress=true, abstol=1e-7, reltol=1e-7, maxiters=1e8)

plot(sol, idxs=[sys.CaJSR, sys.CaNSR])
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL])
# Sub SR Ca instability even there's no external stimuli
# Jup unstable? PDE unstable?
@unpack JCa_SR, V_sub_SR= sys
plot(sol, idxs=[sys.JCa_SR/V_sub_SR])
plot(sol, idxs=[sys.Jup, sys.Jrel])

u0 = sol[end]
@time sol = solve(remake(prob, u0=u0), alg; progress=true, abstol=1e-9, reltol=1e-9)
