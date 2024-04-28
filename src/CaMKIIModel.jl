module CaMKIIModel
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using NaNMath

export build_neonatal_model

include("utils.jl")
include("neonatal.jl")
end # module
