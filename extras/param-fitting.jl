# # Parameter fitting
# Fitting CaMKII parameters to experimental decay rates
using ModelingToolkit
using OrdinaryDiffEq
using Optimization
using OptimizationLBFGSB
using ADTypes
using ForwardDiff
using CaMKIIModel
using CaMKIIModel: second

# Setup simulation without CaMKII A2 state
tend = 210second
stimstart = 30second
@time "Build system" @mtkcompile sys = build_neonatal_ecc_sys()
@unpack r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1, Istim = sys
@time "Build problem" prob = ODEProblem(sys, [k_P1_P2 => 0], tend; saveat=[stimstart])
@time "Solve problem" sol = solve(prob, FBDF()) ## Baseline simulation

# Hyperparameters
function get_hyp(sol;
    stimstart=30second,
    durations=[15.0, 30.0, 60.0, 90.0],
    experiment_taus=[16.48, 16.73, 17.65, 18.08],
    )
    prob = sol.prob
    sys = prob.f.sys
    @unpack kb_CaMKP, kdeph_CaMK, kphos_CaMK, CaMKAct = sys

    saveats = map(durations) do duration
        stimend = stimstart + duration * second
        range(stimend; step=1second, length=51) |> collect
    end

    return (;
        durations,
        experiment_taus,
        stimstart,
        saveats = saveats,
        kphos_dephos_ratio=prob.ps[kphos_CaMK] / prob.ps[kdeph_CaMK],
        a0=sol[CaMKAct] |> only,
        prob = prob,
        callbacks=[build_stim_callbacks(Istim, stimstart + duration * second; period=1second, starttime=stimstart) for duration in durations]
    )
end
ps = get_hyp(sol)

# Loss function
function loss(x, ps)
    @unpack prob, kphos_dephos_ratio, stimstart, durations, experiment_taus, callbacks, a0, saveats = ps
    @unpack kb_CaMKP, kdeph_CaMK, kphos_CaMK, CaMKAct = ps.prob.f.sys
    errors = zero(eltype(x))
    ## Parallel ensemble simulation
    function prob_func(prob, i, repeat)
        remake(prob, p=[kb_CaMKP => exp10(x[1]), kdeph_CaMK => exp10(x[2]), kphos_CaMK => kphos_dephos_ratio * exp10(x[2])], callback=callbacks[i], saveat=saveats[i])
    end

    ensemble_prob = EnsembleProblem(prob; prob_func)
    sim = solve(ensemble_prob, FBDF(); trajectories=length(durations), maxiters=100000, verbose=false)

    for i in 1:4
        sol = sim[i]
        stimend = stimstart + durations[i] * second
        peakAct = sol(stimend, idxs=CaMKAct)
        SciMLBase.successful_retcode(sol) || return Inf
        ysim = sol[CaMKAct]
        ts = range(0, length=length(ysim), step=1)
        yexpected = a0 .+ (peakAct - a0) .* exp.(-ts ./ experiment_taus[i])
        errors += sum(abs2, ysim .- yexpected) / length(ysim)
    end
    return errors
end

x0 = [log10(prob.ps[kb_CaMKP]), log10(prob.ps[kdeph_CaMK])]
@time loss(x0, ps)

#---
optf = OptimizationFunction(loss, ADTypes.AutoForwardDiff())
optprob = OptimizationProblem(optf, x0, ps, lb=[-0.5, -2] + x0, ub=[2, 2] + x0)
@time sol = solve(optprob, OptimizationLBFGSB.LBFGSB())
sol.objective
sol.original

#---
fit = Dict(kb_CaMKP => exp10(sol.u[1]), kdeph_CaMK => exp10(sol.u[2]), kphos_CaMK => ps.kphos_dephos_ratio * exp10(sol.u[2]))

println("Fitted parameters:")
println("Dissociation time of CaMKP --> CaMKA:" , 1e-3 / fit[kb_CaMKP], " seconds.")
println("Dephosphorylation time of CaMKA:" , 1e-3 / fit[kdeph_CaMK], " seconds.")
println("Phosphorylation rate of CaMKB:" , fit[kphos_CaMK] * 1000, " Hz")
