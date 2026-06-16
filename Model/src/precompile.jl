using CSV
using DataFrames
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using NaNMath
using OrdinaryDiffEqSDIRK
using Plots
using SteadyStateDiffEq

@compile_workload begin
    sys = build_neonatal_ecc_sys() |> mtkcompile
    prob = ODEProblem(sys, [], (0.0, 200.0second))
    sprob = SteadyStateProblem(sys, [])
    sol = solve(prob, KenCarp4())
end
