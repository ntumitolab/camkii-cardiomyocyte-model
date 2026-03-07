# # Sensitiviy analysis of the CaMKII system
# Influence of CaMKII parameters on activity decay rate at 1Hz pacing for 300 sec.
# Load packages
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq
using CaMKIIModel
using CaMKIIModel: ms, Hz, second
import SciMLSensitivity as SMS
import ForwardDiff as FD
using Plots
Plots.default(lw=2)

# ## Setup the ODE system
# Electrical stimulation starts at `t`=100 sec and ends at `t`=300 sec.
@time @mtkcompile sys = build_neonatal_ecc_sys()
tend = 500.0second
@time prob = ODEProblem(sys, [], tend)
stimstart = 100.0second
stimend = 300.0second
@unpack Istim = sys

# ## 1Hz pacing
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time sol = solve(prob, KenCarp47(); callback, reltol = 1e-8, abstol = 1e-8);
# Sensitivity analysis of the solution at `t`=300 sec against all parameters.
@unpack r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1 = sys

function sens(x; t=300.0second)
    _prob = remake(prob, p=[sys.r_CaMK => x[1], sys.kb_CaMKP => x[2], sys.kphos_CaMK => x[3], sys.kdeph_CaMK => x[4], sys.k_P1_P2 => x[5], sys.k_P2_P1 => x[6]])
    sol = solve(_prob, KenCarp47(); callback, reltol = 1e-8, abstol = 1e-8, saveat = [t])
    sol(t, idxs = sys.CaMKAct)
end

#---
x = [3Hz, inv(3second), 5Hz, inv(6second), inv(60second), inv(15second)];
@time dx = FD.gradient(sens, x)

#---
println("Relative sensitivity of CaMKII activity at 1Hz pacing for 300 sec to parameters:")
for (i, k) in enumerate([r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1])
     println("$k : ", dx[i] * x[i] / sol(300.0second, idxs=sys.CaMKAct))
end
