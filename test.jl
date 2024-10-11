using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel
using CaMKIIModel: μM
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true)
sys_r = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 500.0
prob = ODEProblem(sys, [], tend)
prob_r = ODEProblem(sys_r, [], tend)
stimstart = 100.0
stimend = 300.0
alg = FBDF()

@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)

@time sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

@time sol_r = solve(prob_r, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

plot(sol, idxs=sys.PKACII/sys.RIItot, title="PKA activation", lab="PKACII")
plot!(sol, idxs=sys.PKACI/sys.RItot, lab="PKACI")

plot(sol_r, idxs=sys_r.PKACII/sys_r.PKACIItot, title="PKA activation", lab="PKACII")
plot!(sol_r, idxs=sys_r.PKACI/sys_r.PKACItot, lab="PKACI")

plot(sol, idxs=sys.Cai_mean, tspan=(299, 300), title="Calcium transcient", lab="Full model")
plot!(sol_r, idxs=sys_r.Cai_mean, tspan=(299, 300), lab="Reduced model")

plot(sol_r, idxs=[sys_r.Cai_sub_SR, sys_r.Cai_sub_SL, sys_r.Cai_mean], tspan=(299, 300), title="Calcium transcient")


prob2 = remake(prob, p=[sys.ISO => 1E-4])
sol2 = solve(prob2, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

prob2_r = remake(prob_r, p=[sys_r.ISO => 1E-4])
sol2_r = solve(prob2_r, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

plot(sol2, idxs=sys.Cai_mean, tspan=(299, 300), title="Calcium transcient", lab="Full model")
plot!(sol2_r, idxs=sys_r.Cai_mean, tspan=(299, 300), lab="Reduced model")
