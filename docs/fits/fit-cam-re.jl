# # CaM rapid binding
# Simplify the CaM binding to CaMKII by assuming that the binding of CaM to Ca is rapid and at equilibrium. This allows us to reduce the number of states in the model and focus on the steady-state behavior of CaMKII activation as a function of calcium concentration.
using Model
using Model: μM, hil, second, Hz
using ADTypes
using CurveFit
using DiffEqCallbacks
using ForwardDiff
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
# with phospoprlation disabled
@parameters Ca = 0μM ROS = 0μM
@time "Build system" sys = Model.get_camkii_sys(; Ca=Ca, ROS=ROS) |> mtkcompile
@time "Build problem" camprob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])

# Physiological cytosolic calcium levels ranges from 30nM to 10μM.
# FIXME: use a loop to make it faster
ca = logrange(0.03μM, 10μM, 101)
prob_func(prob, ctx) = remake(prob, p=[Ca => ca[ctx.sim_id]])
ensemble_prob = EnsembleProblem(camprob; prob_func)
@time "Solve problem" sim = solve(ensemble_prob, DynamicSS(KenCarp47()), EnsembleThreads(); trajectories=length(ca), abstol=1e-10, reltol=1e-10)

"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k], sim.u)

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
@time "Build problem" camprob_re = SteadyStateProblem(sys_re, [sys_re.kphos_CaMK => 0])
@time "Solve a single problem" sol_re = solve(camprob_re, DynamicSS(KenCarp47()); abstol=1e-10, reltol=1e-10)

sol_re[sys_re.CaM0]

ensemble_prob_re = EnsembleProblem(camprob_re; prob_func)
@time "Solve problem" sim_re = solve(ensemble_prob_re, DynamicSS(KenCarp47()), EnsembleThreads(); trajectories=length(ca), abstol=1e-10, reltol=1e-10)

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

    ## Parallel ensemble simulation
    function prob_func(prob, ctx)
        i = ctx.sim_id
        remake(prob, p=[
            KEQ_CAMC => keq_camc,
            KEQ_CAMN => keq_camn,
            KEQ_KCAMC => keq_kcamc,
            KEQ_KCAMN => keq_kcamn,
            Ca => ca[i]
            ],
        )
    end

    ## Calculate loss in the output function
    function output_func(sol, ctx)
        SciMLBase.successful_retcode(sol) || return (Inf, false)
        i = ctx.sim_id
        loss = (sol(sys_re.CaM0) - data.CaM0[i])^2 +
                (sol(sys_re.CaM2C) - data.CaM2C[i])^2 +
                (sol(sys_re.CaM2N) - data.CaM2N[i])^2 +
                (sol(sys_re.CaM4) - data.CaM4[i])^2 +
                (sol(sys_re.CaMK) - data.CaMK[i])^2 +
                (sol(sys_re.CaMKB0) - data.CaMKB0[i])^2 +
                (sol(sys_re.CaMKB2C) - data.CaMKB2C[i])^2 +
                (sol(sys_re.CaMKB2N) - data.CaMKB2N[i])^2 +
                (sol(sys_re.CaMKB4) - data.CaMKB4[i])^2
        return (loss, false)
    end

    ensemble_prob_re = EnsembleProblem(camprob_re; prob_func, output_func)
    sim_re = solve(ensemble_prob_re, DynamicSS(KenCarp47()), EnsembleThreads(); trajectories=length(ca), abstol=1e-10, reltol=1e-10)
    loss = sum(sim_re)
end

# Test the loss function
theta0 = log10.([camprob_re.ps[KEQ_CAMC], camprob_re.ps[KEQ_CAMN], camprob_re.ps[KEQ_KCAMC], camprob_re.ps[KEQ_KCAMN]])

@time loss(theta0, data)
g = ForwardDiff.gradient((theta) -> loss(theta, data), theta0)
# ### Optimization
optf = OptimizationFunction(loss, ADTypes.AutoForwardDiff())
optprob = OptimizationProblem(optf, theta0, data, lb=[-1, -1, -1, -1] + theta0, ub=[1, 1, 1, 1] + theta0)
sol = solve(optprob, LBFGSB())
