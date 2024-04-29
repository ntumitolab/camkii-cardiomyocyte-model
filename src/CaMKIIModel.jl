module CaMKIIModel
using ModelingToolkit
using NaNMath

export build_neonatal_model

include("utils.jl")
include("neonatal.jl")
end # module
