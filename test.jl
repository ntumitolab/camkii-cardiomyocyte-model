using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 50000.0
prob = ODEProblem(sys, [], tend)
alg = FBDF()
callback = TerminateSteadyState()

sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

unknowns(sys)
sol.u[end]

for (k, v) in zip(unknowns(sys), sol.u[end])
    println(k, " => ", v , ",")
end
