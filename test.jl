using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using CaMKIIModel

sys = build_neonatal_ecc_sys(simplify=true)
tend = 10.0
prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend)
@time sol = solve(prob, Rodas5P(); callback, progress=true)

for (k, v) in zip(unknowns(sys), sol[end])
    println(k, " => ", v)
end
