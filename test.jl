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
@time sol = solve(prob, Rodas5P(); callback, progress=true, maxiters=1e8)

length(sys.Cai)

[(x, y) for x in 1:3, y in 1:2]

plot(25:0.01:30, 1:length(sys.Cai), (t, x) -> sol(t, idxs=sys.Cai[x]), linetype=:wireframe, size=(800, 800))


plot(sol, idxs=sys.vm*1000)
plot(sol, idxs=[sys.TnI_PKAp])
plot(sol, idxs=[INaCa, ICaL, ICaT, ICab], size=(800, 800), tspan=(20, 30))
