# # CaM rapid binding
# Simplify the CaM binding to CaMKII by assuming that the binding of CaM to Ca is rapid and at equilibrium. This allows us to reduce the number of states in the model and focus on the steady-state behavior of CaMKII activation as a function of calcium concentration.
using Model
using Model: μM, hil, second, Hz
using ADTypes
using CurveFit
using DiffEqCallbacks
using ForwardDiff
using Zygote
using LinearAlgebra
using ModelingToolkit
using Optimization
using OptimizationOptimJL
using OptimizationLBFGSB
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using Plots
using SteadyStateDiffEq
Plots.default(lw=1.5)

# ## Full model
# with phosphorylation disabled
tend = 1000second
alg = KenCarp47()
@parameters Ca = 0μM ROS = 0μM
@time "Build system" sys = Model.get_camkii_sys(; Ca=Ca, ROS=ROS) |> mtkcompile
@time "Build problem" camprob = ODEProblem(sys, [sys.k_phosCaM => 0], tend)

"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k][end], sim)

"""Solve the problem for a range of calcium concentrations. Assuming smooth transitions."""
function solve_range(prob, ca)
    u0 = prob.u0
    map(ca) do c
        _prob = remake(prob; u0=u0, p=[Ca => c])
        sol = solve(_prob, alg; abstol=1e-10, reltol=1e-10, save_everystep=false, save_start=false)
        u0 = sol.u[end]
        sol
    end
end

# Physiological cytosolic calcium levels ranges from 30nM to 10μM.
ca = logrange(0.03μM, 10μM, 101)
@time "Solve problem" sim = solve_range(camprob, ca)

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
    plot!(title="A", titlelocation=:left, ylims = (0, 70))
end

# ## Rapid CaM binding to Ca
@time "Build system" sys_re = Model.get_camkii_dia_sys(; Ca=Ca, ROS=ROS) |> mtkcompile
@time "Build problem" camprob_re = ODEProblem(sys_re, [sys_re.kphos_CaMK => 0], tend)
@time "Solve problem" sim_re = solve_range(camprob_re, ca)

figs1b = let
    plot(ca, extract(sim_re, sys_re.CaM2C), lab="CaM2C", ylabel="Conc. (μM)"; xopts...)
    plot!(ca, extract(sim_re, sys_re.CaM2N), lab="CaM2N")
    plot!(ca, extract(sim_re, sys_re.CaM4), lab="CaM4")
    plot!(ca, extract(sim_re, sys_re.CaMK), lab="CaMK")
    plot!(ca, extract(sim_re, sys_re.CaMKB0), lab="CaMKB0")
    plot!(ca, extract(sim_re, sys_re.CaMKB2C), lab="CaMKB2C")
    plot!(ca, extract(sim_re, sys_re.CaMKB2N), lab="CaMKB2N")
    plot!(ca, extract(sim_re, sys_re.CaMKB4), lab="CaMKB4")
    plot!(title="B", titlelocation=:left, legend=:left, ylims = (0, 70))
end

#---
figs1c = let
    plot(ca, extract(sim_re, sys_re.fCaM0), lab="fCaM0", ylabel="Conc. (μM)")
    plot!(ca, extract(sim_re, sys_re.fCaM2C), lab="fCaM2C")
    plot!(ca, extract(sim_re, sys_re.fCaM2N), lab="fCaM2N")
    plot!(ca, extract(sim_re, sys_re.fCaM4), lab="fCaM4")
    plot!(ca, extract(sim_re, sys_re.fKCaM0), lab="fKCaM0", linestyle=:dot, ylabel="Conc. (μM)")
    plot!(ca, extract(sim_re, sys_re.fKCaM2C), lab="fKCaM2C", linestyle=:dot)
    plot!(ca, extract(sim_re, sys_re.fKCaM2N), lab="fKCaM2N", linestyle=:dot)
    plot!(ca, extract(sim_re, sys_re.fKCaM4), lab="fKCaM4", linestyle=:dot)
    plot!(title="C", titlelocation=:left, legend=:left; xopts...)
end

#---
figs1d = let
    plot(ca, extract(sim, sys.CaMKAct), lab="Full model", ylabel="Active CaMKII fraction")
    plot!(ca, extract(sim_re, sys_re.CaMKAct), lab="Rapid CaM binding", linestyle=:dash)
    plot!(title="D", titlelocation=:left, legend=:left, ylims = (0, 0.5) ;xopts...)
end

# ## Fitting the rapid CaM binding model
data = (
    CaM0 = extract(sim, sys.CaM0_CaMK),
    CaM2C = extract(sim, sys.Ca2CaM_C),
    CaM2N = extract(sim, sys.Ca2CaM_N),
    CaM4 = extract(sim, sys.Ca4CaM),
    CaMK = extract(sim, sys.CaMK),
    CaMKB0 = extract(sim, sys.CaM0_CaMK),
    CaMKB2C = extract(sim, sys.Ca2CaM_C_CaMK),
    CaMKB2N = extract(sim, sys.Ca2CaM_N_CaMK),
    CaMKB4 = extract(sim, sys.Ca4CaM_CaMK),
)

@unpack KEQ_CAMC, KEQ_CAMN, KEQ_KCAMC, KEQ_KCAMN = camprob_re.f.sys

function loss(theta, data)
    keq_camc = exp10(theta[1])
    keq_camn = exp10(theta[2])
    keq_kcamc = exp10(theta[3])
    keq_kcamn = exp10(theta[4])

    _prob = remake(camprob_re; p=[
            KEQ_CAMC => keq_camc,
            KEQ_CAMN => keq_camn,
            KEQ_KCAMC => keq_kcamc,
            KEQ_KCAMN => keq_kcamn,
        ]
    )

    sim_re = solve_range(_prob, ca)

    loss = sum(abs2, extract(sim_re, sys_re.CaM0) .- data.CaM0) +
           sum(abs2, extract(sim_re, sys_re.CaM2C) .- data.CaM2C) +
           sum(abs2, extract(sim_re, sys_re.CaM2N) .- data.CaM2N) +
           sum(abs2, extract(sim_re, sys_re.CaM4) .- data.CaM4) +
           sum(abs2, extract(sim_re, sys_re.CaMK) .- data.CaMK) +
           sum(abs2, extract(sim_re, sys_re.CaMKB0) .- data.CaMKB0) +
           sum(abs2, extract(sim_re, sys_re.CaMKB2C) .- data.CaMKB2C) +
           sum(abs2, extract(sim_re, sys_re.CaMKB2N) .- data.CaMKB2N) +
           sum(abs2, extract(sim_re, sys_re.CaMKB4) .- data.CaMKB4)
end

# Test the loss function
theta0 = log10.([camprob_re.ps[KEQ_CAMC], camprob_re.ps[KEQ_CAMN], camprob_re.ps[KEQ_KCAMC], camprob_re.ps[KEQ_KCAMN]])

@time loss(theta0, data)
g = ForwardDiff.gradient((theta) -> loss(theta, data), theta0)
# ### Optimization
optf = OptimizationFunction(loss, ADTypes.AutoForwardDiff())
optprob = OptimizationProblem(optf, theta0, data, lb=[-1, -1, -1, -1] + theta0, ub=[1, 1, 1, 1] + theta0)
@time sol = solve(optprob, LBFGSB())

prob_fit = remake(camprob_re; p=[
        KEQ_CAMC => exp10(sol[1]),
        KEQ_CAMN => exp10(sol[2]),
        KEQ_KCAMC => exp10(sol[3]),
        KEQ_KCAMN => exp10(sol[4]),
    ]
)

sim_fit = solve_range(prob_fit, ca)

figs1e = let
    plot(ca, extract(sim, sys.CaMKAct), lab="Full model", ylabel="Active CaMKII fraction")
    plot!(ca, extract(sim_re, sys_re.CaMKAct), lab="Rapid CaM binding", linestyle=:dash)
    plot!(ca, extract(sim_fit, sys_re.CaMKAct), lab="Rapid CaM binding (fitted)", linestyle=:dash)
    plot!(title="D", titlelocation=:left, legend=:left, ylims = (0, 0.5) ;xopts...)
end
