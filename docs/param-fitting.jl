# # Parameter fitting
# Fitting CaMKII parameters to experimental decay rates
using ModelingToolkit
using SymbolicIndexingInterface
using OrdinaryDiffEq
using CurveFit
using Optimization
using OptimizationBBO
using CaMKIIModel
using CaMKIIModel: second

# Setup simulation
tend = 210second
@time "Build system" @mtkcompile sys = build_neonatal_ecc_sys()
@unpack r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1, Istim = sys
@time "Build problem" prob = ODEProblem(sys, [k_P1_P2 => 0], tend)

# Decay times for 15sec, 30 sec, 60 sec, and 90 sec pacing
# Fitted in `Pacing durations`
# Loss function
function loss(x; alg = KenCarp47(), stimstart = 30second)
    experimental_durations = [16.48, 16.73, 17.65, 18.08]
    _p = remake(prob, p=[r_CaMK => x[1], kb_CaMKP => x[2], kphos_CaMK => x[3], kdeph_CaMK => x[4]])
    simulation_durations = map((15.0, 30.0, 60.0, 90.0)) do duration
        callback = build_stim_callbacks(Istim, stimstart + duration * second; period=1second, starttime=stimstart)
        ts = collect(range(0.0, stop=50.0, step=5.0))
        saveat = (ts .+ duration) .* second .+ stimstart
        sol = solve(_p, alg; callback=callback, saveat=saveat, maxiters = 10000, verbose = false)
        SciMLBase.successful_retcode(sol) || return Inf
        ysim = sol(saveat; idxs=sys.CaMKAct * 100).u
        fit_sim = solve(CurveFitProblem(ts, ysim), ExpSumFitAlgorithm(n=1, withconst=true))
        tau_sim = inv(-fit_sim.u.λ[])
    end
    return sum(abs2, experimental_durations .- simulation_durations)
end

p0 = [prob.ps[r_CaMK], prob.ps[kb_CaMKP], prob.ps[kphos_CaMK], prob.ps[kdeph_CaMK]]
@time loss(p0)

optf = OptimizationFunction((x, p) -> loss(x))
optprob = OptimizationProblem(optf, p0)
@time sol = Optimization.solve(optprob, Optim.SimulatedAnnealing())
