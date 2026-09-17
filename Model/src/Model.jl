module Model

using ADTypes
using CSV
using CurveFit
using DataFrames
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using NaNMath
using Optimization
using OptimizationLBFGSB
using OptimizationOptimJL
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using Plots
using Statistics

export build_neonatal_ecc_sys, build_stim_callbacks

include("utils.jl")
include("camkii_ros.jl")
include("isoproterenol.jl")
include("istim.jl")
include("capde.jl")
include("ica.jl")
include("ik.jl")
include("ina.jl")
include("er.jl")
include("ecc_neonatal.jl")
include("precompile.jl")

const DEFAULT_SYS = build_neonatal_ecc_sys() |> mtkcompile
const DEFAULT_PROB = ODEProblem(DEFAULT_SYS, [], (0.0, 100.0second))

end # module Model
