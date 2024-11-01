# # Initial conditions
# No stimulation for 50000 seconds.
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
using CaMKIIModel
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 50000.0
prob = SteadyStateProblem(sys, [])
alg = DynamicSS(Rodas5P())

sol = solve(prob, alg; abstol=1e-8, reltol=1e-8)

for (k, v) in zip(unknowns(sys), sol.u)
    println(k, " => ", v , ",")
end
