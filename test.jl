using CaMKIIModel
using ModelingToolkit
using OrdinaryDiffEq
using BenchmarkTools
using Plots

sys = build_neonatal_ecc_sys(simplify=true)

prob = ODEProblem(sys, [], float(100))

sol = solve(prob, TRBDF2())

# TODO: fix electrophysiology
plot(sol, idxs=sys.vm*1000)

sol[i_CK2]

for s in unknowns(sys)
    println(s, " => ", sol(sol.t[end], idxs=s))
end

for eq in equations(sys)
    println(eq)
end
