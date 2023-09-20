module CaMKIIModel

export make_b1ar_rn

include("utils.jl")
include("bar.jl")
include("cam.jl")
include("ecc/ica_markov.jl")
include("ecc/ina_hh.jl")
include("ecc/ina_markov.jl")

end # module CaMKIIModel
