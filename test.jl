using Catalyst
using ModelingToolkit

@parameters begin
    dCa=1
end
@variables t Ca(t) = 1.0
D = Differential(t)
eq = [D(Ca) ~ -dCa * (Ca)]
@named osys = ODESystem(eq, t)

rn = @reaction_network begin
    $Ca, A + B --> C
end

setdefaults!(rn, [:A=>1.0, :B=> 1.0, :C=>0.0])

rnosys = convert(ODESystem, rn)

sys = extend(osys, rnosys)

simpsys = structural_simplify(sys)
