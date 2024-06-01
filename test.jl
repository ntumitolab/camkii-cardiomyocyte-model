using Catalyst, OrdinaryDiffEq, Plots, ModelingToolkit

# set the growth rate to 1.0
@parameters λ = -1.0

# set the initial volume to 1.0
@variables t V(t) = 1.0

# build the ODESystem for dV/dt
D = Differential(t)
eq = [D(V) ~ λ * V]
@named osys = ODESystem(eq, t)

# build the ReactionSystem with no protein initially
rn = @reaction_network begin
    @species P(t) = 1.0 Q(t) = 0.0
    $V, P --> Q
end

conservationlaw_constants(rn)

rnsys = convert(ODESystem, rn, remove_conserved = true)

sys = extend(osys, rnsys)

simp = structural_simplify(sys)

prob = ODEProblem(simp, [])

ts = (0, 10)

sol = solve(prob, Tsit5(); tspan=ts)

plot(sol)
