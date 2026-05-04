# # Initial conditions
using Model
using ModelingToolkit
using DifferentialEquations
using SteadyStateDiffEq
using Plots
Plots.default(lw=2)

#---
@time "Build system" sys = Model.DEFAULT_SYS
@time "Build problem" sprob = SteadyStateProblem(sys, [])
@time "Solve problem" sol = solve(sprob, DynamicSS(KenCarp4()); abstol=1e-8, reltol=1e-8)

for (k, v) in zip(unknowns(sys), sol.u)
    println(k, " => ", v, ",")
end
