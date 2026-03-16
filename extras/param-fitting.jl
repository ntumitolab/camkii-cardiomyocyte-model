# # Parameter fitting
# Fitting CaMKII parameters to experimental decay rates
using ModelingToolkit
using SymbolicIndexingInterface
using OrdinaryDiffEq
using Optimization
using OptimizationBBO
using CaMKIIModel
using CaMKIIModel: second

# Setup simulation
tend = 210second
@time "Build system" @mtkcompile sys = build_neonatal_ecc_sys()
@unpack r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1, Istim = sys
@time "Build problem" prob = ODEProblem(sys, [k_P1_P2 => 0], tend)

# Hyperparameters
ps = (;
    stimstart = 30second,
    experiment_durations = [16.48, 16.73, 17.65, 18.08],
    callbacks = [build_stim_callbacks(Istim, stimstart + duration * second; period=1second, starttime=stimstart) for duration in (15.0, 30.0, 60.0, 90.0)],
    kphos_dephos = prob.ps[kphos_CaMK] / prob.ps[kdeph_CaMK]
)

# Decay times for 15sec, 30 sec, 60 sec, and 90 sec pacing
# Fitted in `Pacing durations`
# Loss function
function loss(x, ps; alg=KenCarp47(), stimstart=30second)
    _p = remake(prob, p=[kb_CaMKP => x[1], kdeph_CaMK => x[2], kphos_CaMK => kphos_dephos * x[2]])
    simulation_durations = map(zip((15.0, 30.0, 60.0, 90.0), callbacks)) do (duration, callback)
        ts = collect(range(0.0, stop=50.0, step=5.0))
        saveat = (ts .+ duration) .* second .+ stimstart
        sol = solve(_p, alg; callback=callback, saveat=saveat, maxiters=10000, verbose=false)
        SciMLBase.successful_retcode(sol) || return Inf
        ysim = sol(saveat; idxs=sys.CaMKAct).u
        fit_sim = solve(CurveFitProblem(ts, ysim), ExpSumFitAlgorithm(n=1, withconst=true))
        tau_sim = inv(-fit_sim.u.λ[])
    end
    return sum(abs2, experiment_durations .- simulation_durations)
end

x0 = [prob.ps[kb_CaMKP], prob.ps[kdeph_CaMK]]
@time loss(x0, ps)

optprob = OptimizationProblem(loss, x0, ps, lb=0.1 .* x0, ub=10.0 .* x0)
@time sol = Optimization.solve(optprob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters = 100000)
