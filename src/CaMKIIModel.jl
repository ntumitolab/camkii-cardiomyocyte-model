module CaMKIIModel

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using NaNMath
using SciMLBase  # CallbackSet
using DiffEqCallbacks

export get_camkii_sys, get_bar_sys, build_neonatal_ecc_sys, build_stim_callbacks, get_bar_sys_reduced

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

end # module
