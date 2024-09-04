using CaMKIIModel
using ModelingToolkit
using OrdinaryDiffEq
using BenchmarkTools
using Plots

sys = build_neonatal_ecc_sys()

tend = float(1000)
prob = ODEProblem(sys, [], tend)

# Unstable
sol = solve(prob, TRBDF2(), maxiters=1E8, abstol=1e-7, reltol=1e-7)

@unpack CK0, i_CK1, i_CK2, i_OK, i_IK = sys
plot(sol, idxs=[])

sol[i_CK2]

for s in unknowns(sys)
    println(s)
end

for eq in equations(sys)
    println(eq)
end
