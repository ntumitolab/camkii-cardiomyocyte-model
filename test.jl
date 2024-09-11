using CaMKIIModel
using ModelingToolkit
using OrdinaryDiffEq
using BenchmarkTools
using Plots

sys = build_neonatal_ecc_sys(simplify=true)

prob = ODEProblem(sys, [], float(10))

# TODO: fix unstable (ICaL?)
sol = solve(prob, Rodas5P())

sol[sys.ICaL]

plot(sol, idxs=sys.PO1RyR)

for (k, v) in zip(unknowns(sys), sol[end])
    println(k, " => ", v)
end
