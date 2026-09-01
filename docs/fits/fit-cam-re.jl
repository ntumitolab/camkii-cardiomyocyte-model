# # CaM rapid binding
# Simplify the CaM binding to CaMKII by assuming that the binding of CaM to Ca is rapid and at equilibrium. This allows us to reduce the number of states in the model and focus on the steady-state behavior of CaMKII activation as a function of calcium concentration.
using Model
using Model: μM, hil, second, Hz
using CurveFit
using DiffEqCallbacks
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using ModelingToolkit
using Plots
using SteadyStateDiffEq
Plots.default(lw=1.5)

# ## Full model
# with phospoprlation disabled
@parameters Ca = 0μM ROS = 0μM
@time "Build system" sys = Model.get_camkii_sys(; Ca=Ca, ROS=ROS) |> mtkcompile
@time "Build problem" camprob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])

# Physiological cytosolic calcium levels ranges from 30nM to 10μM.
ca = logrange(0.03μM, 10μM, 1001)
@time "Solve problem" sim = map(ca) do c
    newprob = remake(camprob, p=[Ca => c])
    solve(newprob, DynamicSS(KenCarp47()); abstol=1e-8, reltol=1e-8)
end;

"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k], sim)

# CaMKII system composition across physiological calcium levels.
xopts = (xlabel="Ca (μM)", xscale=:log10, minorgrid=true, xlims=(ca[1], ca[end]))
figs1a = let
    plot(ca, extract(sim, sys.Ca2CaM_C), lab="Ca2CaM_C", ylabel="Conc. (μM)"; xopts...)
    plot!(ca, extract(sim, sys.Ca2CaM_N), lab="Ca2CaM_N")
    plot!(ca, extract(sim, sys.Ca4CaM), lab="Ca4CaM")
    plot!(ca, extract(sim, sys.CaMK), lab="CaMK")
    plot!(ca, extract(sim, sys.CaM0_CaMK), lab="CaM0_CaMK")
    plot!(ca, extract(sim, sys.Ca2CaM_C_CaMK), lab="Ca2CaM_C_CaMK")
    plot!(ca, extract(sim, sys.Ca2CaM_N_CaMK), lab="Ca2CaM_N_CaMK")
    plot!(ca, extract(sim, sys.Ca4CaM_CaMK), lab="Ca4CaM_CaMK", legend=:left)
    plot!(title="A", titlelocation=:left)
end

# ## Rapid CaM binding to Ca
@time "Build system" sys_re = Model.get_camkii_dia_sys(; Ca=Ca, ROS=ROS) |> mtkcompile

observed(sys_re)

@time "Build problem" camprob_re = SteadyStateProblem(sys_re, [sys_re.kphos_CaMK => 0])

@time "Solve problem" sim_re = map(ca) do c
    newprob = remake(camprob_re, p=[Ca => c])
    solve(newprob, DynamicSS(KenCarp47()); abstol=1e-8, reltol=1e-8)
end;

figs1b = let
    plot(ca, extract(sim_re, sys_re.CaM2C), lab="CaM2C", ylabel="Conc. (μM)"; xopts...)
    plot!(ca, extract(sim_re, sys_re.CaM2N), lab="CaM2N")
    plot!(ca, extract(sim_re, sys_re.CaM4), lab="CaM4")
    plot!(ca, extract(sim_re, sys_re.CaMK), lab="CaMK")
    plot!(ca, extract(sim_re, sys_re.CaMKB0), lab="CaMKB0")
    plot!(ca, extract(sim_re, sys_re.CaMKB2C), lab="CaMKB2C")
    plot!(ca, extract(sim_re, sys_re.CaMKB2N), lab="CaMKB2N")
    plot!(ca, extract(sim_re, sys_re.CaMKB4), lab="CaMKB4", legend=:left)
    plot!(title="A", titlelocation=:left)
end
