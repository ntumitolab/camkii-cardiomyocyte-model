using Catalyst

var = 3

rn = @reaction_network begin
    @parameters k1 = ($var<2) ? 1 : 2
    k1, A --> 0
end
