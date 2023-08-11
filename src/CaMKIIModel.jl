module CaMKIIModel

export make_b1ar_rn

include("utils.jl")
include("bar.jl")
include("cam.jl")
include("ions/ica_markov.jl")
include("ions/ina_hh.jl")
include("ions/ina_markov.jl")

end # module CaMKIIModel
