using PrecompileTools: @compile_workload, @recompile_invalidations

@recompile_invalidations begin
    using ADTypes
    using CSV
    using CurveFit
    using DataFrames
    using DiffEqCallbacks
    using DifferentialEquations
    using ModelingToolkit
    using NaNMath
    using Optimization
    using OptimizationOptimJL
    using OrdinaryDiffEq
    using OrdinaryDiffEqSDIRK
    using Plots
    using StatsBase
    using StatsPlots
end

@compile_workload begin
    DEFAULT_SYS = build_neonatal_ecc_sys() |> mtkcompile
    DEFAULT_PROB = ODEProblem(DEFAULT_SYS, [], (0.0, 100.0second))
end
