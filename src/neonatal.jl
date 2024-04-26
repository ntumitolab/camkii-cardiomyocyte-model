# Neonatal mouse CMC model
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

include("neonatal_ca_pde.jl")

function build_neonatal_model()

    @constants begin
        Vmyo = 2.1454e-11           # [L]
        Vdyad = 1.7790e-014         # [L]
        VSL = 6.6013e-013           # [L]
        kSLmyo = 8.587e-15          # [L/msec]
        CaMKIItotDyad = 120         # [uM]
        BtotDyad = 1.54 / 8.293e-4  # [uM]
    end

end
