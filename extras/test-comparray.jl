import ComponentArrays as CA
import ComponentArrays.ComponentArray as CArray
import ComponentArrays.ComponentVector as CVector
using SimpleUnPack
using Optimization
using OptimizationLBFGSB
using ForwardDiff
using ADTypes
using OrdinaryDiffEq

function lotka_volterra!(du, u, p, t)
    @unpack x, y = u
    @unpack α, β, δ, γ = p
    du.x = α * x - β * x * y
    du.y = -δ * y + γ * x * y
    return nothing
end

u0 = CVector(x=1.0, y=1.0)
tsteps = 0.0:0.1:10.0
p = CVector(α=1.5, β=1.0, δ=3.0, γ=1.0)
prob = ODEProblem(lotka_volterra!, u0, (tsteps[begin], tsteps[end]), p)
sol = solve(prob, Tsit5(), saveat=tsteps) ## Reference solution

function loss(x, p)
    ## Cannot use `CArray(prob.p; α=x.α, β=x.β)` directly due to different types (Dual vs Float64)
    _prob = remake(prob, p=CArray(α=x.α, β=x.β, δ=p.δ, γ=p.γ))
    _sol = solve(_prob, Tsit5(), saveat=tsteps, maxiters=100000, verbose=false)
    SciMLBase.successful_retcode(_sol) || return Inf
    return sum(abs2, Array(_sol) .- Array(sol))
end

optf = OptimizationFunction(loss, ADTypes.AutoForwardDiff())
x0 = CArray(α=1.0, β=0.5)
optprob = OptimizationProblem(optf, x0, p)
optsol = solve(optprob, LBFGSB())
