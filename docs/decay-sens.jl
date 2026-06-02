# # Sensitiviy analysis of the CaMKII system
# Influence of parameters on CaMKII activity decay specific time at 1Hz pacing for 300 sec.
using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEqSDIRK
using CurveFit
using Plots
Plots.default(lw=2)
using Model
using Model: Hz, second

# ## Setup the ODE system
# Electrical stimulation starts at `t`=100 sec and ends at `t`=300 sec.
@time sys = Model.DEFAULT_SYS
tend = 400.0second
@time prob = ODEProblem(sys, [], tend)
stimstart = 100.0second
stimend = 300.0second
@unpack Istim = sys

# ## Simulations
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, KenCarp4(); callback, reltol = 1e-8, abstol = 1e-8)

# Get the decay time of CaMKII activity from `t`=300 sec to `t`=350 sec by fitting the solution to an exponential function. The decay time is the inverse of the exponent coefficient.
function get_decay_time(sol; stimend=300.0second)
    xs = collect(0:1second:50.0second)
    ys = sol(xs .+ stimend, idxs=sys.CaMKAct).u
    fit = solve(CurveFitProblem(xs, ys), ExpSumFitAlgorithm(n=1, withconst=true))
    inv(-fit.u.λ[]) ./ 1000
end

@time tau0 = get_decay_time(sol; stimend=300.0second)
println("The default decay time of CaMKII activity at 1Hz pacing for 200 sec is: ", tau0, " seconds")

# Sensitivity analysis of the solution at `t`=300 sec against parameters.
@unpack r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1 = sys
params = [r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1]
sensitivities = Dict()

@time for k in params
    println("Calculating sensitivity for parameter: ", k)
    original_value = prob.ps[k]
    _prob = remake(prob, p=[k => original_value * 1.01]) ## Increase 1% of the parameter value
    sol = solve(_prob, KenCarp4(); callback, reltol = 1e-8, abstol = 1e-8)
    sensitivities[k] = (get_decay_time(sol; stimend) / tau0 - 1) * 100
end

sensitivities

#---
println("Relative sensitivity of CaMKII activity at 1Hz pacing for 200 seconds to parameters:")
for (k, v) in sensitivities
     println("$k : ", v)
end
