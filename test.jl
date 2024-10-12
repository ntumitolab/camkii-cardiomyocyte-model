using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using CaMKIIModel
using CaMKIIModel: μM
Plots.default(lw=2)

sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 10000.0
prob = ODEProblem(sys, [], tend)
alg = FBDF()
callback = TerminateSteadyState()
sol = solve(prob, alg; callback) ## warmup

sol[sys.PLBp]
println(sol[sys.CaMKP2][end])
