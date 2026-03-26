# # Initial conditions
# Load packages
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq
using Plots
using CaMKIIModel
Plots.default(lw=2)

#---
@time "Build system" sys = CaMKIIModel.DEFAULT_SYS
@time "Build problem" prob = SteadyStateProblem(sys, [])
@time "Solve problem" sol = solve(prob, DynamicSS(KenCarp47()); abstol=1e-10, reltol=1e-10)
for (k, v) in zip(unknowns(sys), sol.u)
    println(k, " => ", v, ",")
end
