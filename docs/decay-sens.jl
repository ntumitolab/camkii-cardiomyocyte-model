# # Sensitiviy analysis of the CaMKII system
# Influence of parameters on CaMKII activity decay specific time at 1Hz pacing for 300 sec.
# Load packages
using ModelingToolkit
using OrdinaryDiffEq
using CurveFit
using Plots
Plots.default(lw=2)
using Model
using Model: Hz, second

# ## Setup the ODE system
# Electrical stimulation starts at `t`=100 sec and ends at `t`=300 sec.
@time @mtkcompile sys = build_neonatal_ecc_sys()
tend = 500.0second
@time prob = ODEProblem(sys, [], tend)
stimstart = 100.0second
stimend = 300.0second
@unpack Istim = sys

# ## Simulations
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, KenCarp47(); callback, reltol = 1e-8, abstol = 1e-8)

# Default decay time
decay_model(p, x) = @. p[1] * exp(-x / p[2]) + p[3]

function get_decay_time(sol; stimend=300.0second)
    xs = collect(0:1second:50.0second)
    ys = sol(xs .+ stimend, idxs=sys.CaMKAct).u
    fit = solve(CurveFitProblem(xs, ys), ExpSumFitAlgorithm(n=1, withconst=true))
    inv(-fit.u.λ[]) ./ 1000
end

@time tau0 = get_decay_time(sol; stimend=300.0second)
println("The default decay time of CaMKII activity at 1Hz pacing for 300 sec is: ", tau0, " seconds")

# Sensitivity analysis of the solution at `t`=300 sec against parameters.
@unpack r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1 = sys
sensitivities = Dict()

for k in (r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1)
    println("Calculating sensitivity for parameter: ", k)
    original_value = prob.ps[k]
    _prob = remake(prob, p=[k => original_value * 1.01]) ## Increase 1% of the parameter value
    sol = solve(_prob, KenCarp47(); callback, reltol = 1e-8, abstol = 1e-8)
    sensitivities[k] = (get_decay_time(sol; stimend) / tau0 - 1) * 100
end

sensitivities

#---
println("Relative sensitivity of CaMKII activity at 1Hz pacing for 300 sec to parameters:")
for (k, v) in sensitivities
     println("$k : ", v)
end
