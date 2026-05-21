using OrdinaryDiffEqSDIRK
using SteadyStateDiffEq
using Plots

@compile_workload begin
    sys = build_neonatal_ecc_sys() |> mtkcompile
    prob = ODEProblem(sys, [], (0.0, 200.0second))
    sprob = SteadyStateProblem(sys, [])
    sol = solve(prob, KenCarp4())
    plot(sol, idxs=sys.CaMKAct)
end
