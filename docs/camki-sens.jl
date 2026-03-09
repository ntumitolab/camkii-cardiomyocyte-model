# # Sensitiviy analysis of the CaMKII system
# Influence of parameters on CaMKII activity decay specific time at 1Hz pacing for 300 sec.
# Load packages
using ModelingToolkit
using OrdinaryDiffEq
using CurveFit
import SciMLSensitivity as SMS
import ForwardDiff as FD
using Plots
Plots.default(lw=2)
using CaMKIIModel
using CaMKIIModel: Hz, second

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
function sens(x; t=300.0second)
    _prob = remake(prob, p=[sys.r_CaMK => x[1], sys.kb_CaMKP => x[2], sys.kphos_CaMK => x[3], sys.kdeph_CaMK => x[4], sys.k_P1_P2 => x[5], sys.k_P2_P1 => x[6]])
    sol = solve(_prob, KenCarp47(); callback, reltol = 1e-8, abstol = 1e-8)
    get_decay_time(sol; stimend)
end

#---
x = [3Hz, inv(3second), 5Hz, inv(6second), inv(60second), inv(15second)];
@time dx = FD.gradient(sens, x)

#---
println("Relative sensitivity of CaMKII activity at 1Hz pacing for 300 sec to parameters:")
for (i, k) in enumerate([r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1])
     println("$k : ", dx[i] * x[i] / tau0)
end
