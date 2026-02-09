# # Initial conditions
# Load packages
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CaMKIIModel
Plots.default(lw=2)

#---
@time "Build system" @mtkcompile sys = build_neonatal_ecc_sys()
@time "Build problem" prob = SteadyStateProblem(sys, [])
@time "Solve problem" sol = solve(prob, DynamicSS(KenCarp47()); abstol=1e-10, reltol=1e-10)
for (k, v) in zip(unknowns(sys), sol.u)
    println(k, " => ", v, ",")
end

# ## Disable CaMKA2
# Setting the rate towards (`k_P1P2`) to zero
@time "Build problem" prob2 = remake(prob, p=[sys.k_P1_P2 => 0])
@time "Solve problem" sol2 = solve(prob2, DynamicSS(KenCarp47()); abstol=1e-10, reltol=1e-10)
for (k, v) in zip(unknowns(sys), sol2.u)
    println(k, " => ", v, ",")
end
