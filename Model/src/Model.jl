module Model

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using NaNMath
using OrdinaryDiffEq
using DiffEqCallbacks
using PrecompileTools: @setup_workload, @compile_workload

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
const DEFAULT_PROB = ODEProblem(DEFAULT_SYS, [], (0.0, 200.0second))

end # module Model
