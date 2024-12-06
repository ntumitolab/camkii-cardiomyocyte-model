# # Initial conditions
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
using CaMKIIModel
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
prob = SteadyStateProblem(sys, [])
alg = DynamicSS(Rodas5P())

sol = solve(prob, alg; abstol=1e-10, reltol=1e-10)

for (k, v) in zip(unknowns(sys), sol.u)
    println(k, " => ", v , ",")
end
