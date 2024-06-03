using Catalyst, OrdinaryDiffEq, Plots, ModelingToolkit

rn = @reaction_network begin
    1, P --> Q
end

@unpack P, Q = rn
defaults = [P=>1.0, Q=>0.0]
setdefaults!(rn, defaults)

@variables t total(t)
obseq = [total ~ 2 * P + Q]
rnsys = convert(ODESystem, rn; remove_conserved = true, observed = obseq)

simp = structural_simplify(rnsys)

prob = ODEProblem(simp, [])
