using CaMKIIModel
using ModelingToolkit
using OrdinaryDiffEq
using BenchmarkTools
using Plots

sys = build_neonatal_ecc_sys(simplify=true)

unknowns(sys)
observed(sys)

tend = float(1)
prob = ODEProblem(sys, [], tend)

# Unstable
sol = solve(prob, TRBDF2())

sol(sol.t[end], idxs=sys.CK0)

@unpack CK0, i_CK1, i_CK2, i_OK, i_IK, vm = sys

plot(sol, idxs=[CK0, i_CK1, i_CK2, i_OK, i_IK])
plot(sol, idxs=vm*1000)

sol[i_CK2]

for s in unknowns(sys)
    println(s, " => ", sol(sol.t[end], idxs=s))
end

for eq in equations(sys)
    println(eq)
end
